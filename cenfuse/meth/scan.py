"""One pass over the chromosomes at one geometry. Per sample: the CHG adjusted track. Calibration pass (no direction):
class sums pooled over samples at each anchor and at every window where the class has data away from the anchor, and
each class's genome-wide calls.
Candidate pass (a direction): each sample's picks with their class sums in every context, the plant-specific score,
and the data at each anchor."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import params
from .bins import n_bins
from .candidates import half_width, pick, robust_z, sign
from .track import Geometry, adjusted, class_sums, slide


@dataclass
class Pick:
    sample: str
    chrom: str
    rank: int
    bin: int
    value: float                 # CHG adjusted track at the pick
    left: int                    # half-depth span, first bin
    right: int                   # half-depth span, last bin
    z: float                     # robust z among the chromosome's windows
    num: dict = field(default_factory=dict)   # context -> class sums of deviations at the pick, float[K]
    vok: dict = field(default_factory=dict)   # context -> class calls with a reference, float[K]
    ps: float = float("nan")     # plant-specific score: this sample minus the median of the others
    ps_p: float = float("nan")   # share of the chromosome's windows at least as plant-specific
    ps_n: int = 0                # other samples with data at the pick


@dataclass
class AnchorData:
    sample: str
    chrom: str
    bin: int
    sites: float                 # covered CHG cytosines in the anchor window (all classes)
    calls: float                 # CHG calls with a class reference in the window (all classes)
    num: dict = field(default_factory=dict)
    vok: dict = field(default_factory=dict)


@dataclass
class Pooled:
    """One context x anchored chromosome, summed over samples. Row K (after the K classes) = all classes."""
    anchor_num: np.ndarray       # float[K + 1] at the anchor window
    anchor_vok: np.ndarray
    elsewhere: dict              # class row -> (num, vok) at every window with >= MIN_CALLS calls that lies >= sep
                                 # from the anchor and from both ends


@dataclass
class Scan:
    geo: Geometry
    adj: dict          # sample -> chrom -> float[bins], CHG adjusted track
    pooled: dict       # context -> chrom -> Pooled (calibration pass only)
    totals: dict       # context -> float[K]: windowed calls per class over all chromosomes and samples (class shares)
    picks: list
    at_anchor: list


def median_of_others(stack: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Column median over finite rows and the number of finite rows; NaN where none is finite (no warnings)."""
    k = np.isfinite(stack).sum(axis=0)
    srt = np.sort(stack, axis=0)                     # NaN sorts last
    lo, hi = np.clip((k - 1) // 2, 0, None), np.clip(k // 2, 0, None)
    med = (np.take_along_axis(srt, lo[None], 0)[0] + np.take_along_axis(srt, hi[None], 0)[0]) / 2
    return np.where(k >= 1, med, np.nan), k


def pooled_at(num: np.ndarray, vok: np.ndarray, a: int, sep: int) -> Pooled:
    """num, vok: classes x bins summed over samples; a: anchor bin."""
    num = np.vstack([num, num.sum(axis=0)])
    vok = np.vstack([vok, vok.sum(axis=0)])
    i = np.arange(num.shape[1])
    away = (i >= sep) & (i < len(i) - sep) & (np.abs(i - a) >= sep)
    ok = (vok >= params.MIN_CALLS) & away
    return Pooled(num[:, a].copy(), vok[:, a].copy(),
                  {int(k): (num[k, ok[k]], vok[k, ok[k]]) for k in np.flatnonzero(ok.any(axis=1))})


def scan(bins, samples, chroms: dict[str, int], n_classes: int, geo: Geometry, anchor_bin: dict[str, int],
         direction: str | None) -> Scan:
    adj = {s.name: {} for s in samples}
    pooled, totals, picks, at_anchor = {}, {}, [], []
    for c, L in chroms.items():
        nb = n_bins(L)
        mine, pool = {}, {}
        for smp in samples:
            for ctx in smp.files:                      # CHG first (io.CONTEXTS order)
                n_, v, m = bins[(smp.name, ctx)][c].dense(n_classes, nb)
                num, vok = class_sums(v, m, geo)
                if ctx == "CHG":
                    y = adjusted(num, vok)
                    adj[smp.name][c] = y
                    if direction is not None:
                        mine[smp.name] = []
                        for r, j in enumerate(pick(y, geo, direction), 1):
                            lft, rgt = half_width(y, j, direction)
                            mine[smp.name].append(Pick(smp.name, c, r, j, float(y[j]), lft, rgt, robust_z(y, j)))
                        if c in anchor_bin:
                            j = anchor_bin[c]
                            at_anchor.append(AnchorData(smp.name, c, j, float(slide(n_.sum(axis=0), geo.w)[j]),
                                                        float(vok.sum(axis=0)[j])))
                if direction is not None:
                    for p in mine.get(smp.name, []):
                        p.num[ctx], p.vok[ctx] = num[:, p.bin].copy(), vok[:, p.bin].copy()
                    if c in anchor_bin:
                        a = next(x for x in at_anchor if x.sample == smp.name and x.chrom == c)
                        a.num[ctx], a.vok[ctx] = num[:, anchor_bin[c]].copy(), vok[:, anchor_bin[c]].copy()
                else:
                    totals[ctx] = totals.get(ctx, 0.0) + vok.sum(axis=1)
                    if c in anchor_bin:
                        acc = pool.setdefault(ctx, [np.zeros_like(num), np.zeros_like(vok)])
                        acc[0] += num
                        acc[1] += vok
        for ctx, (num, vok) in pool.items():
            pooled.setdefault(ctx, {})[c] = pooled_at(num, vok, anchor_bin[c], geo.sep)
        if direction is not None:
            s = sign(direction)
            names = [x.name for x in samples]
            for smp_name, plist in mine.items():
                if len(names) < 3 or not plist:
                    picks += plist
                    continue
                med, k = median_of_others(np.stack([adj[o][c] for o in names if o != smp_name]))
                ps = adj[smp_name][c] - med
                fin = ps[np.isfinite(ps)] * s
                for p in plist:
                    p.ps, p.ps_n = float(ps[p.bin]), int(k[p.bin])
                    if np.isfinite(p.ps) and fin.size:
                        p.ps_p = float((fin >= p.ps * s).sum() / fin.size)
                picks += plist
    return Scan(geo, adj, pooled, totals, picks, at_anchor)
