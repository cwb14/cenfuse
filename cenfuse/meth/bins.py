"""Per sample x context x chromosome: 50-kb class bins (covered cytosines n, summed N_valid v, summed methylated
calls m), streamed from tabix-indexed bedMethyl in 10-Mb regions, in parallel, cached per chromosome. Then the depth
filter and the exclusions drop a sample's unreliable bins."""
from __future__ import annotations

import hashlib
import logging
import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from . import params
from .classes import ClassMap, classify
from .io import read_bedmethyl

log = logging.getLogger("cenmeth")
CACHE_VERSION = 1


@dataclass
class ChromBins:
    """Sparse class bins of one sample x context x chromosome (only non-empty bin x class cells)."""
    bin: np.ndarray   # int32
    cls: np.ndarray   # int16
    n: np.ndarray     # int64 covered cytosines
    v: np.ndarray     # int64 summed N_valid
    m: np.ndarray     # int64 summed methylated calls

    def dense(self, n_classes: int, n_bins: int):
        """(n, v, m) as float arrays of shape (classes, bins)."""
        out = []
        for x in (self.n, self.v, self.m):
            a = np.zeros((n_classes, n_bins))
            np.add.at(a, (self.cls, self.bin), x)
            out.append(a)
        return tuple(out)

    def keep_bins(self, keep: np.ndarray) -> ChromBins:
        k = keep[self.bin]
        return ChromBins(self.bin[k], self.cls[k], self.n[k], self.v[k], self.m[k])


def n_bins(length: int) -> int:
    return -(-length // params.BIN)


def bin_chromosome(path, chrom: str, length: int, segs, other: int, n_classes: int) -> tuple[ChromBins, int]:
    """Stream one chromosome of one bedMethyl file -> (ChromBins, malformed line count)."""
    nb = n_bins(length)
    acc = np.zeros((3, nb * n_classes), dtype=np.int64)
    bad = 0
    pysam.set_verbosity(0)                 # htslib notes (e.g. index age) are reported once by io.read_samples
    with pysam.TabixFile(str(path)) as tb:
        for s in range(0, length, params.CHUNK):
            pos, v, m, b = read_bedmethyl(tb, chrom, s, min(s + params.CHUNK, length), params.MIN_CALLS_SITE)
            bad += b
            if len(pos) == 0:
                continue
            key = (pos // params.BIN) * n_classes + classify(segs, pos, other)
            for i, w in enumerate((None, v, m)):
                acc[i] += np.bincount(key, weights=w, minlength=nb * n_classes).astype(np.int64)
    nz = np.flatnonzero(acc[0])
    return ChromBins((nz // n_classes).astype(np.int32), (nz % n_classes).astype(np.int16),
                     acc[0, nz], acc[1, nz], acc[2, nz]), bad


def fingerprint(paths, extra: str = "") -> str:
    """Fast input fingerprint: name, size, mtime, first and last MiB of each file (not a full checksum)."""
    h = hashlib.md5(f"{CACHE_VERSION}:{params.BIN}:{params.MIN_CALLS_SITE}:{extra}".encode())
    n = params.CACHE_FP_BYTES
    for p in paths:
        p = Path(p)
        st = p.stat()
        h.update(f"{p.name}:{st.st_size}:{st.st_mtime_ns}".encode())
        with open(p, "rb") as f:
            h.update(f.read(n))
            if st.st_size > 2 * n:
                f.seek(-n, 2)
                h.update(f.read())
    return h.hexdigest()[:16]


def _job(a):
    path, chrom, length, segs, other, k, cache_file, key = a
    if cache_file is not None and cache_file.exists():
        try:
            z = np.load(cache_file)
            if str(z["key"]) == key:
                return ChromBins(z["bin"], z["cls"], z["n"], z["v"], z["m"]), 0, True
        except Exception:
            pass
    cb, bad = bin_chromosome(path, chrom, length, segs, other, k)
    if cache_file is not None:
        tmp = cache_file.with_name(cache_file.name + f".{os.getpid()}.tmp.npz")
        np.savez(tmp, key=np.array(key), bin=cb.bin, cls=cb.cls, n=cb.n, v=cb.v, m=cb.m)
        os.replace(tmp, cache_file)
    return cb, bad, False


def collect(samples, chroms: dict[str, int], cmap: ClassMap, ann_key: str, cache_dir, threads: int):
    """All samples x contexts x chromosomes -> {(sample, ctx): {chrom: ChromBins}}. ann_key = fingerprint of the
    annotation inputs (a change re-bins everything)."""
    cache_dir = Path(cache_dir) if cache_dir is not None else None
    if cache_dir is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)
    jobs, ids = [], []
    k = len(cmap.names)
    for smp in samples:
        for ctx, path in smp.files.items():
            fkey = fingerprint([path], f"{ann_key}:{k}")
            for c, L in chroms.items():
                cf = cache_dir / f"{smp.name}.{ctx}.{c}.npz" if cache_dir is not None else None
                jobs.append((path, c, L, cmap.segs.get(c), cmap.other, k, cf, fkey))
                ids.append((smp.name, ctx, c))
    out, bad, cached = {}, 0, 0
    if threads > 1:
        with ProcessPoolExecutor(threads) as ex:
            res = list(ex.map(_job, jobs, chunksize=1))
    else:
        res = [_job(a) for a in jobs]
    for (s, ctx, c), (cb, b, hit) in zip(ids, res):
        out.setdefault((s, ctx), {})[c] = cb
        bad += b
        cached += hit
    if bad:
        log.warning("bedMethyl: %d malformed lines skipped", bad)
    log.info("class bins: %d sample x context x chromosome jobs (%d from cache)", len(jobs), cached)
    return out


@dataclass
class DepthStats:
    sample: str
    median_depth: float
    bins_with_data: int
    dropped_depth: int
    dropped_excluded: int


def depth_filter(bins, samples, chroms: dict[str, int], excl: pd.DataFrame | None):
    """Drop, per sample, CHG-defined bins whose mean depth per covered cytosine is outside DEPTH_LO-DEPTH_HI x the
    sample's genome median, plus its excluded regions; applied to every context. Returns (bins, [DepthStats])."""
    stats, out = [], {}
    for smp in samples:
        chg = bins[(smp.name, "CHG")]
        dep = {}
        for c, L in chroms.items():
            cb = chg[c]
            nn = np.bincount(cb.bin, weights=cb.n, minlength=n_bins(L))
            vv = np.bincount(cb.bin, weights=cb.v, minlength=n_bins(L))
            dep[c] = np.divide(vv, nn, out=np.full(len(nn), np.nan), where=nn > 0)
        alld = np.concatenate(list(dep.values()))
        med = float(np.nanmedian(alld)) if np.isfinite(alld).any() else float("nan")
        keep, n_dep, n_exc = {}, 0, 0
        for c, L in chroms.items():
            d = dep[c]
            bad = np.isfinite(d) & ((d < params.DEPTH_LO * med) | (d > params.DEPTH_HI * med))
            ex = np.zeros(len(d), bool)
            if excl is not None and len(excl):
                e = excl[(excl.chrom == c) & ((excl["sample"] == "") | (excl["sample"] == smp.name))]
                for s0, e0 in zip(e.start.to_numpy(), e.end.to_numpy()):
                    ex[s0 // params.BIN:-(-e0 // params.BIN)] = True
            n_dep += int((bad & ~ex & np.isfinite(d)).sum())
            n_exc += int((ex & np.isfinite(d)).sum())
            keep[c] = ~(bad | ex)
        for ctx in smp.files:
            out[(smp.name, ctx)] = {c: bins[(smp.name, ctx)][c].keep_bins(keep[c]) for c in chroms}
        stats.append(DepthStats(smp.name, med, int(np.isfinite(alld).sum()), n_dep, n_exc))
    return out, stats
