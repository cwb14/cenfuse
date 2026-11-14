"""The repeat-family line: centromeres built from interspersed repeats, found from the FASTA alone (0.16).

Where centromeres are clusters of particular retrotransposon families (the Tekay-type centromeres of Solanum, the CRM-rich
centromeres of grasses) and no satellite family passes its one-per-chromosome test, the FASTA still holds what the LTR line
reads from an annotation: families concentrated at one locus on most chromosomes. K-mers stand in for the annotation. A
FracMinHash sample of canonical K-mers keeps every copy of each kept K-mer, so its copy number is exact. A K-mer is
concentrated on a chromosome where at least MIN_COPIES of its copies lie there and at least half of them inside one window;
a family K-mer is concentrated on at least half of the chromosomes taking part. Each family K-mer votes once per chromosome
where it is concentrated, at the center of its densest window, and the votes per window are read by the AT line's
one-per-chromosome test (atline.scan_profiles) at windows of 0.1-3 Mb.

Left out: copies inside the satellite engine's arrays (the satellite line's; a K-mer with at least half of its copies
inside them is a satellite K-mer and is dropped), K-mers with a period of at most K // 2 (tandem: microsatellites, telomeric
repeats), K-mers with fewer than 1.5 N copies (N chromosomes: too few to be concentrated on half of them), and windows whose
midpoint lies in a terminal eighth (as in the repeat-density line).

What a pass means: every family K-mer votes at most once per chromosome, so one family that clusters once per chromosome
puts its votes one per chromosome wherever its clusters sit. A pass says that such a family exists, not where the
centromeres are: the test is not calibrated. With each genome's own copies and every repeat unit moved as a whole to a
random place on its chromosome, it still passed at FAM_ALPHA in 13-52% of draws in each of 11 monocentric genomes (0.16.1);
holocentric genomes do not pass it (C. elegans P 0.44, B. mori P 0.71). No family-level null could be built from
K-mers: grouping them into families splits diverged families or merges neighbors (docs/future.md #76). Its peaks are
leads.

Calls: the strongest vote window at the chosen size, widened to its half-height run. The line places only `low` calls, on
chromosomes nothing else places (fuse: fam_only), and flags a `low` sequence call whose family window lies elsewhere
(fam_elsewhere). It never pairs with Hi-C (an uncalibrated test cannot be one of two tested lines) and never moves, opposes
or raises a sequence call."""
from __future__ import annotations

import math

import numpy as np
import pysam

from . import atline, repline
from .hicfuse import END_ZONE
from .util import log

K = 31                    # K-mer length: odd (no K-mer is its own reverse complement), the longest a uint64 holds; a
                          # family's copies at 90% identity still share 4% of their 31-mers, unrelated sequence none
SAMPLE = 32               # keep 1 K-mer in SAMPLE by hash (FracMinHash), at every genome size: the rate the line was
                          # developed and tested at; a kept K-mer keeps every copy
SCALES = (100_000, 300_000, 1_000_000, 3_000_000)   # a x3 series over regional centromeres, fungi to large-genome plants
MIN_COPIES = 3            # the fewest copies on a chromosome that can show a cluster beyond a pair
FAM_ALPHA = 0.01          # conventional level, after the correction for the window sizes tried (as the repeat-density line)
MIN_COPY_FACTOR = 1.5     # a K-mer needs ceil(1.5 N) copies: the least that can be concentrated (3 copies) on half of N
PIECE = 1_000_000         # bp read at a time (a multiple of atline.BIN): streaming, never a whole chromosome string


def periodic(km: np.ndarray, k: int = K, p_max: int = K // 2) -> np.ndarray:
    """True for K-mers with a period p <= p_max (base i equals base i + p throughout): with p_max <= k // 2 they hold their own
    unit at least twice, so their copies count tandem units at one locus. Periodicity survives reverse complementing, so the
    canonical code is tested directly."""
    out = np.zeros(len(km), bool)
    for p in range(1, p_max + 1):
        out |= (km >> np.uint64(2 * p)) == (km & np.uint64((1 << (2 * (k - p))) - 1))
    return out


def sample(fasta, fai: dict, chroms: dict, arrays: dict | None = None) -> dict:
    """The kept copies of the family-candidate K-mers on the analyzed chromosomes: dict(kid, chrom, pos, nk, called) with kid
    the K-mer's index, chrom the chromosome's index in `chroms`, pos 0-based; called = {chrom: non-N bases per atline.BIN
    bin}. fai: {sequence: length} of the whole FASTA (every sequence counts copies); arrays: {chrom: [(start, end)]} 0-based
    half-open (the satellite engine's arrays)."""
    names = list(chroms)
    idx = {c: i for i, c in enumerate(names)}
    arrays = {c: sorted(v) for c, v in (arrays or {}).items()}
    called = {c: np.zeros(L // atline.BIN + 1, np.int64) for c, L in chroms.items()}
    S = np.uint64(SAMPLE)
    KM, CH, PO, AR = [], [], [], []
    with pysam.FastaFile(str(fasta)) as fa:
        for name, L in fai.items():
            ci = idx.get(name, -1)
            for s in range(0, L, PIECE):
                e = min(L, s + PIECE)
                raw = np.frombuffer(fa.fetch(name, s, min(L, e + K - 1)).encode(), dtype=np.uint8)
                km, ok = repline.canonical_kmers(repline.LUT[raw], K)
                km, ok = km[:e - s], ok[:e - s]                        # K-mers starting in [s, e)
                pos = np.flatnonzero(ok)
                km = km[pos]
                sel = (repline.fmix64(km) % S) == 0
                km, pos = km[sel], pos[sel]
                if ci >= 0:
                    up = raw[:e - s] & np.uint8(0xDF)
                    lo, hi = s // atline.BIN, (e - 1) // atline.BIN + 1
                    called[name][lo:hi] += np.bincount((np.arange(e - s) + s) // atline.BIN - lo, up != 78,
                                                       minlength=hi - lo).astype(np.int64)
                    # a K-mer that touches an array is the array's: one starting up to K - 1 bases before it holds its unit too
                    touch = repline._inside([(max(0, a - K + 1), b) for a, b in arrays.get(name) or ()], s, e)
                    AR.append(touch[pos])
                    PO.append((pos + s).astype(np.int64))
                else:
                    AR.append(np.zeros(len(km), bool))
                    PO.append(np.zeros(len(km), np.int64))
                KM.append(km)
                CH.append(np.full(len(km), ci, np.int16))
    km, ch, po, ar = (np.concatenate(x) for x in (KM, CH, PO, AR))
    del KM, CH, PO, AR
    order = np.argsort(km, kind="stable")
    km, ch, po, ar = km[order], ch[order], po[order], ar[order]
    del order
    n = len(km)
    if n == 0:
        return dict(kid=np.zeros(0, np.int64), chrom=np.zeros(0, np.int16), pos=np.zeros(0, np.int64), nk=0, called=called)
    first = np.r_[True, km[1:] != km[:-1]]
    starts = np.flatnonzero(first)
    counts = np.diff(np.r_[starts, n])
    gid = np.repeat(np.arange(len(starts)), counts)
    in_arr = np.bincount(gid, weights=ar, minlength=len(starts))
    c_lo = math.ceil(MIN_COPY_FACTOR * len(names))
    keep = (counts >= c_lo) & (2 * in_arr < counts) & ~periodic(km[starts])
    sel = keep[gid] & (ch >= 0) & ~ar
    new_id = np.cumsum(keep) - 1
    log.debug(f"repfam: {n:,} sampled K-mer copies, {len(starts):,} distinct; {int(keep.sum()):,} kept (>= {c_lo} copies, "
              f"not satellite or tandem), {int(sel.sum()):,} of their copies read")
    return dict(kid=new_id[gid[sel]].astype(np.int64), chrom=ch[sel], pos=po[sel], nk=int(keep.sum()), called=called)


def _densest(key: np.ndarray, W: int):
    """Per group of a sorted key (group << 32 | position): the group, its copies n, the densest W-window's count d and the
    centre of the first densest window."""
    i = np.arange(len(key))
    cnt = np.searchsorted(key, key + W, side="left") - i                  # copies in [pos, pos + W) of the same group
    g = key >> 32
    start = np.flatnonzero(np.r_[True, g[1:] != g[:-1]])
    n = np.diff(np.r_[start, len(key)])
    grp = np.repeat(np.arange(len(start)), n)
    d = np.maximum.reduceat(cnt, start)
    hit = np.flatnonzero(cnt == d[grp])
    best = hit[np.r_[True, grp[hit][1:] != grp[hit][:-1]]]
    pos = key & 0xFFFFFFFF
    centre = (pos[best] + pos[best + cnt[best] - 1]) // 2
    return g[start], n, d, centre


def profiles_of(fam: dict, chroms: dict):
    """profiles(use, k) for atline.scan_profiles: minus the family votes per window of k bins (step 1 bin), +inf where under
    atline.MIN_CALLED of the bases are called or the midpoint lies in a terminal eighth (END_ZONE)."""
    names = list(chroms)
    nc = len(names)
    idx = {c: i for i, c in enumerate(names)}
    g = fam["kid"] * nc + fam["chrom"].astype(np.int64)
    key = np.sort((g << 32) | fam["pos"])

    def profiles(use, k):
        W = k * atline.BIN
        um = np.zeros(nc, bool)
        um[[idx[c] for c in use]] = True
        kk = key[um[(key >> 32) % nc]]
        out = {}
        centre, ch, v = np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0, bool)
        if len(kk):
            grp, n, d, centre = _densest(kk, W)
            kid, ch = grp // nc, grp % nc
            conc = (n >= MIN_COPIES) & (2 * d >= n)
            m = np.bincount(kid[conc], minlength=fam["nk"])
            v = conc & (2 * m >= len(use))[kid]
        for c in use:
            nb = len(fam["called"][c])
            votes = np.bincount(centre[v & (ch == idx[c])] // atline.BIN, minlength=nb)[:nb]
            if nb < k:
                out[c] = np.zeros(0)
                continue
            vs = repline._win(votes.astype(np.float64), k)
            ok = repline._win(fam["called"][c], k) >= atline.MIN_CALLED * k * atline.BIN
            prof = np.where(ok, -vs, np.inf)
            j = np.arange(len(prof))
            L = chroms[c]
            mid = (j * atline.BIN + np.minimum(L, (j + k) * atline.BIN)) / 2
            prof[np.minimum(mid, L - mid) < END_ZONE * L] = np.inf
            out[c] = prof
        return out
    return profiles


def scan(fam: dict, chroms: dict) -> dict:
    """atline.scan_profiles on the family votes: {'scales', 'best', 'p', 'status'}; each call's `gc` field holds its window's
    family votes (positive). Positions are packed into 32 bits: a chromosome of 2^32 bp or more is not read (untested)."""
    if max(chroms.values(), default=0) >= 1 << 32:
        log.warning("repeat families: a chromosome of 4.29 Gb or more; the line is not run (untested)")
        return dict(scales={}, best=None, p=1.0, status="untested")
    res = atline.scan_profiles(profiles_of(fam, chroms), chroms, SCALES, FAM_ALPHA)
    for s in res["scales"].values():
        for x in s["calls"].values():
            x["gc"] = -x["gc"]
    if res["best"] is not None:
        s = res["scales"][res["best"]]
        log.info(f"repeat families: window {res['best']} bp: the strongest {s['n']} vote windows cover {s['covered']} "
                 f"chromosomes (null {s['null_median']:.0f}), the strongest {s['m']} fall one per chromosome at P "
                 f"{res['p']:.2g} after {len(res['scales'])} sizes: {res['status']}")
    return res
