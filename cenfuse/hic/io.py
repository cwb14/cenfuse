"""Input readers: .fai, Hi-C pairs (4DN or simple 6-column), .cool/.mcool via h5py, BED/GFF3 intervals."""
from __future__ import annotations

import gzip
import os
import re
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from pandas.errors import ParserWarning

from . import params

_SKIP_RE = re.compile(r"Skipping line \d+")  # pandas on_bad_lines="warn" message, one entry per skipped line


class InputError(Exception):
    """Fatal problem with an input file."""


def check_readable(paths) -> None:
    """Fail fast, before any long step: every input path must be an existing, readable file."""
    for p in paths:
        if not os.path.isfile(p):
            raise InputError(f"{p}: no such file")
        if not os.access(p, os.R_OK):
            raise InputError(f"{p}: not readable")


def is_cool(path) -> bool:
    return str(path).endswith((".cool", ".mcool"))


def read_fai(path) -> dict[str, int]:
    """Sequence name -> length, in file order, from a samtools .fai (name<TAB>length<TAB>...). Any line that is not
    of that form (e.g. a FASTA given by mistake) is an InputError naming the file and line."""
    out = {}
    try:
        with open(path) as f:
            for i, line in enumerate(f, 1):
                if not line.strip():
                    continue
                x = line.rstrip("\n").split("\t")
                try:
                    L = int(x[1])
                except (IndexError, ValueError):
                    L = -1
                if L < 0:
                    raise InputError(f"{path}: line {i}: not a .fai line (want name<TAB>length<TAB>...): "
                                     f"{line.strip()[:60]!r}")
                out[x[0]] = L
    except UnicodeDecodeError as e:
        raise InputError(f"{path}: not a text .fai file ({e.reason})") from None
    if not out:
        raise InputError(f"{path}: empty .fai")
    return out


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


@dataclass
class PairsFormat:
    kind: str                     # "4dn" (1-based) or "simple6" (0-based)
    columns: list | None
    n_header: int


def pairs_format(path) -> PairsFormat:
    n, cols = 0, None
    with _open(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            n += 1
            if line.startswith("#columns:"):
                cols = line.split(":", 1)[1].split()
    if n == 0:
        return PairsFormat("simple6", None, 0)
    if cols is None or not {"chrom1", "pos1", "chrom2", "pos2"} <= set(cols):
        raise InputError(f"{path}: header lines but no '#columns:' line naming chrom1 pos1 chrom2 pos2")
    return PairsFormat("4dn", cols, n)


@dataclass
class PairStats:
    read: int = 0
    kept: int = 0
    cis: int = 0
    dropped_mapq: int = 0
    dropped_chrom: int = 0
    malformed: int = 0
    mapq_checked: bool = False
    unknown_chroms: list = field(default_factory=list)


def _next_chunk(it):
    """Next chunk from a pandas reader and the number of lines it skipped. pandas' bad-line warnings are
    captured here so they never reach the caller; other warnings are re-emitted unchanged."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ch = next(it, None)
    skipped = 0
    for w in caught:
        if issubclass(w.category, ParserWarning):
            skipped += len(_SKIP_RE.findall(str(w.message)))
        else:
            warnings.warn(w.message, w.category)
    return ch, skipped


def iter_pairs(path, chroms: dict, stats: PairStats, chunksize: int = params.PAIRS_CHUNK):
    """Yield (ci1, p1, ci2, p2) int64 arrays per chunk; ci = index in `chroms` order; positions 0-based.
    Drops pairs with an end off `chroms`, with MAPQ < params.MAPQ (when mapq1/mapq2 exist) or unparsable positions."""
    fmt = pairs_format(path)
    idx = {c: i for i, c in enumerate(chroms)}
    if fmt.kind == "simple6":
        names, use, off = ["chrom1", "pos1", "s1", "chrom2", "pos2", "s2"], ["chrom1", "pos1", "chrom2", "pos2"], 0
    else:
        names, off = fmt.columns, 1
        use = ["chrom1", "pos1", "chrom2", "pos2"] + [c for c in ("mapq1", "mapq2") if c in names]
        stats.mapq_checked = "mapq1" in names and "mapq2" in names
    # No usecols: with usecols pandas silently accepts over-wide rows; without it they are skipped and reported.
    reader = pd.read_csv(path, sep="\t", header=None, names=names, skiprows=fmt.n_header,
                         chunksize=chunksize, dtype={"chrom1": str, "chrom2": str}, on_bad_lines="warn")
    it = iter(reader)
    while True:
        ch, skipped = _next_chunk(it)
        if ch is None:
            break
        ch = ch[use]
        stats.read += len(ch) + skipped
        stats.malformed += skipped
        p1, p2 = pd.to_numeric(ch.pos1, errors="coerce"), pd.to_numeric(ch.pos2, errors="coerce")
        bad = (p1.isna() | p2.isna()).to_numpy()
        stats.malformed += int(bad.sum())
        ok = ~bad
        if stats.mapq_checked:
            q = ((pd.to_numeric(ch.mapq1, errors="coerce") >= params.MAPQ)
                 & (pd.to_numeric(ch.mapq2, errors="coerce") >= params.MAPQ)).to_numpy()
            stats.dropped_mapq += int((ok & ~q).sum())
            ok &= q
        i1, i2 = ch.chrom1.map(idx), ch.chrom2.map(idx)
        inn = (i1.notna() & i2.notna()).to_numpy()
        if len(stats.unknown_chroms) < 5:
            unk = set(ch.chrom1[i1.isna()].unique()) | set(ch.chrom2[i2.isna()].unique())
            stats.unknown_chroms += sorted(unk - set(stats.unknown_chroms))[: 5 - len(stats.unknown_chroms)]
        stats.dropped_chrom += int((ok & ~inn).sum())
        ok &= inn
        if not ok.any():
            continue
        ci1, ci2 = i1[ok].to_numpy(np.int64), i2[ok].to_numpy(np.int64)
        q1, q2 = p1[ok].to_numpy(np.int64) - off, p2[ok].to_numpy(np.int64) - off
        stats.kept += len(ci1)
        stats.cis += int((ci1 == ci2).sum())
        yield ci1, q1, ci2, q2


def _is_gff(path) -> bool:
    p = str(path).lower()
    if p.endswith((".gff", ".gff3", ".gff.gz", ".gff3.gz")):
        return True
    with _open(path) as f:
        return f.readline().startswith("##gff-version")


def _attrs(s: str) -> dict:
    out = {}
    for kv in s.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def read_intervals(path, attrs: bool = True) -> pd.DataFrame:
    """BED (0-based half-open; optional name) or GFF3 (1-based closed). Columns: chrom, start, end, type, attrs.
    attrs=False leaves the attrs column empty (None): a large repeat annotation then costs far less memory.
    A non-integer coordinate is an InputError naming the file and line; lines with too few columns are skipped."""
    rows = []
    try:
        gff = _is_gff(path)
        kind = "GFF3" if gff else "BED"
        with _open(path) as f:
            for i, line in enumerate(f, 1):
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                x = line.rstrip("\n").split("\t")
                try:
                    if gff and len(x) >= 9:
                        rows.append((x[0], int(x[3]) - 1, int(x[4]), x[2], _attrs(x[8]) if attrs else None))
                    elif not gff and len(x) >= 3:
                        a = ({"Name": x[3]} if len(x) > 3 else {}) if attrs else None
                        rows.append((x[0], int(x[1]), int(x[2]), "region", a))
                except ValueError:
                    raise InputError(f"{path}: line {i}: non-integer coordinate in a {kind} line: "
                                     f"{line.strip()[:80]!r}") from None
    except UnicodeDecodeError as e:
        raise InputError(f"{path}: not a text BED/GFF3 file ({e.reason})") from None
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "type", "attrs"])


def _open_cool(path):
    import h5py
    try:
        return h5py.File(path, "r")
    except OSError as e:
        raise InputError(f"{path}: cannot open as HDF5 (.cool/.mcool): {e}") from e


_COOL_GROUPS = ("chroms", "bins", "pixels")


def _require_cool(g, path, where: str) -> None:
    missing = [k for k in _COOL_GROUPS if k not in g]
    if missing:
        raise InputError(f"{path}: {where} lacks {', '.join(missing)} group(s); not a .cool/.mcool")


def _bin_size(f, path) -> int:
    try:
        return int(f.attrs["bin-size"])
    except KeyError as e:
        raise InputError(f"{path}: .cool has no 'bin-size' attribute") from e


def _mcool_resolutions(f) -> list:
    return sorted(int(r) for r in f["resolutions"].keys())


def cool_resolutions(path) -> list:
    with _open_cool(path) as f:
        if "resolutions" in f:
            return _mcool_resolutions(f)
        _require_cool(f, path, "file root")
        return [_bin_size(f, path)]


def read_cool(path, resolution: int, chroms: dict):
    """Counts of a .cool / .mcool resolution on Hi-C engine bins of the same size over `chroms` (in that order).
    Pixels involving sequences absent from `chroms`, or bins beyond a chromosome's .fai length, are dropped."""
    import scipy.sparse as sp

    from .bins import Bins, symmetrize
    with _open_cool(path) as f:
        if "resolutions" in f:
            avail = _mcool_resolutions(f)
            if resolution not in avail:
                raise InputError(f"{path}: resolution {resolution} not in .mcool; available: {avail}")
            where = f"resolutions/{resolution}"
            g = f[where]
        else:
            bs = _bin_size(f, path)
            if bs != resolution:
                raise InputError(f"{path}: .cool has bin-size {bs}, requested {resolution}")
            where, g = "file root", f
        _require_cool(g, path, where)
        try:
            names = [n.decode() if isinstance(n, bytes) else str(n) for n in g["chroms/name"][:]]
            bchrom, bstart = g["bins/chrom"][:], g["bins/start"][:]
            b1, b2, cnt = g["pixels/bin1_id"][:], g["pixels/bin2_id"][:], g["pixels/count"][:]
        except KeyError as e:
            raise InputError(f"{path}: {where} is missing a required table: {e}") from e
    bins = Bins(chroms, resolution)
    idx = {c: i for i, c in enumerate(chroms)}
    ci = np.array([idx.get(n, -1) for n in names])[bchrom]
    loc = bstart // resolution
    valid = (ci >= 0) & (loc < bins.nb[np.maximum(ci, 0)])
    glob = np.where(valid, bins.off[np.maximum(ci, 0)] + loc, -1)
    a, b = glob[b1], glob[b2]
    keep = (a >= 0) & (b >= 0)
    lo, hi = np.minimum(a[keep], b[keep]), np.maximum(a[keep], b[keep])
    M = sp.coo_matrix((cnt[keep].astype(np.float64), (lo, hi)), shape=(bins.N, bins.N)).tocsr()
    return symmetrize(M)
