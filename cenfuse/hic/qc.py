"""QC: flags on calls, distance to the satellite centromere, and the read-class split of centromere contacts."""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

from . import params
from .io import PairStats, iter_pairs


def merge_intervals(df) -> dict:
    out = {}
    for c, d in df.groupby("chrom"):
        d = d.sort_values("start")
        ms, me = [], []
        for s, e in zip(d.start.to_numpy(), d.end.to_numpy()):
            if ms and s <= me[-1]:
                me[-1] = max(me[-1], e)
            else:
                ms.append(s)
                me.append(e)
        out[c] = (np.array(ms, dtype=np.int64), np.array(me, dtype=np.int64))
    return out


def inside(iv, chrom, pos):
    pos = np.asarray(pos)
    if chrom not in iv or len(iv[chrom][0]) == 0:
        return np.zeros(len(pos), bool)
    s, e = iv[chrom]
    j = np.searchsorted(s, pos, side="right") - 1
    return (j >= 0) & (pos < e[np.clip(j, 0, len(e) - 1)])


def hub_bins(M, bins, good):
    """Bins whose total coverage and inter-chromosomal contacts are both >= HUB_FOLD x the median (good bins)."""
    cov = np.asarray(M.sum(axis=1)).ravel()
    K = len(bins.names)
    C = sp.csr_matrix((np.ones(bins.N), (np.arange(bins.N), bins.chrom_of)), shape=(bins.N, K))
    per = np.asarray((M @ C).todense())
    trans = cov - per[np.arange(bins.N), bins.chrom_of]
    g = good & (cov > 0)
    return g & (cov >= params.HUB_FOLD * np.median(cov[g])) & (trans >= params.HUB_FOLD * np.median(trans[g]))


def satellite_centres(df, bins) -> dict:
    """Midpoint (bp) of each chromosome's satellite centromere: GFF 'centromere' features if any, else every row
    (a BED of satellite centromeres); several per chromosome -> the longest."""
    d = df[df.type == "centromere"] if (df.type == "centromere").any() else df
    idx = {c: i for i, c in enumerate(bins.names)}
    out = {}
    for c, x in d.groupby("chrom"):
        if c in idx:
            r = x.loc[(x.end - x.start).idxmax()]
            out[idx[c]] = (r.start + r.end) / 2
    return out


def rdna_intervals(df) -> dict:
    keep = df["attrs"].map(lambda a: a.get("rdna", ".") not in ("", ".", None))
    return merge_intervals(df[keep]) if keep.any() else {}


def call_flags(k, pc, bins, hub, rdna_iv, rabl_passed, p_null, converged, cycle_unstable=False) -> list:
    """Flags of one call (flags only: they never change the confidence level)."""
    f = []
    g = int(bins.off[k] + pc["peak"])
    if hub[g]:
        f.append("hub")
    s = int(bins.off[k] + pc["core_lo"])
    e = int(bins.off[k] + pc["core_hi"])
    chrom = bins.names[k]
    if chrom in rdna_iv:
        rs, re_ = rdna_iv[chrom]
        if ((rs < bins.end[e]) & (re_ > bins.start[s])).any():
            f.append("rDNA")
    if pc["peak"] < params.END_WIN or pc["peak"] >= bins.nb[k] - params.END_WIN:
        f.append("telomeric")
    if not rabl_passed:
        f.append("no_rabl")
    if not np.isfinite(p_null):
        f.append("null_untestable")
    if not converged:
        f.append("not_converged")
    if cycle_unstable:
        f.append("cycle_unstable")
    return f


def read_class_split(path, chroms, bins, members, pcs, est, repeat_iv) -> pd.DataFrame:
    """For each call: trans contacts of read ends in its core (and on the rest of the chromosome, 'background',
    outside core +/- HALF_WIN) split by whether the end lies in a repeat; 'cen' = the other end lies in another
    member chromosome's partner window (estimate +/- HALF_WIN)."""
    inwin = np.zeros(bins.N, bool)
    for y, c in est.items():
        lo, hi = max(c - params.HALF_WIN, 0), min(c + params.HALF_WIN, int(bins.nb[y]) - 1)
        inwin[bins.off[y] + lo: bins.off[y] + hi + 1] = True
    zone = np.full(bins.N, -1)                      # 1 core, 0 background, -1 ignore
    for k, pc in pcs.items():
        sl = bins.sl(k)
        z = np.zeros(int(bins.nb[k]), int)
        z[max(pc["core_lo"] - params.HALF_WIN, 0): pc["core_hi"] + params.HALF_WIN + 1] = -1
        z[pc["core_lo"]: pc["core_hi"] + 1] = 1
        zone[sl] = z
    counts = {}
    for ci1, p1, ci2, p2 in iter_pairs(path, chroms, PairStats()):
        for ca, pa, cb, pb in ((ci1, p1, ci2, p2), (ci2, p2, ci1, p1)):
            ok = ((ca != cb) & members[ca] & members[cb]
                  & (pa >= 0) & (pa < bins.lengths[ca]) & (pb >= 0) & (pb < bins.lengths[cb]))
            ga, gb = bins.glob(ca[ok], pa[ok]), bins.glob(cb[ok], pb[ok])
            za = zone[ga]
            keep = za >= 0
            for k in np.unique(ca[ok][keep]):
                m = (ca[ok][keep] == k)
                pos = pa[ok][keep][m]
                rep = inside(repeat_iv, bins.names[int(k)], pos)
                cen = inwin[gb[keep][m]]
                zz = za[keep][m]
                for zv, zn in ((1, "core"), (0, "background")):
                    for rv, rn in ((False, "non_repeat"), (True, "repeat")):
                        sel = (zz == zv) & (rep == rv)
                        key = (bins.names[int(k)], zn, rn)
                        t, c = counts.get(key, (0, 0))
                        counts[key] = (t + int(sel.sum()), c + int(cen[sel].sum()))
    rows = [dict(chrom=k[0], zone=k[1], cls=k[2], trans=v[0], cen=v[1]) for k, v in counts.items()]
    t = pd.DataFrame(rows, columns=["chrom", "zone", "cls", "trans", "cen"])
    t["share"] = t.cen / t.trans.where(t.trans > 0)
    return t
