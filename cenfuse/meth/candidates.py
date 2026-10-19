"""Candidates: per sample x chromosome up to MAX_CAND distinct extremes of the adjusted track in the learned
direction, >= 4 D apart; their evidence (half-depth width, robust z, class deviations, signature match, nearest
anchor, plant-specific score) and the anchored report. Labels are 'centromere_like' (signature score >= 0.5:
closer to the anchors' signature than to the background) or 'dip'; never 'centromere'."""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import params
from .bins import n_bins


def sign(direction: str) -> float:
    """+1 when centromeres gain methylation ('up'), -1 when they lose it ('down')."""
    return 1.0 if direction == "up" else -1.0


def pick(y: np.ndarray, geo, direction: str) -> list[int]:
    """Bins of up to MAX_CAND extremes of y in the direction, each >= geo.sep bins from the previous ones; only
    windows that deviate in the direction (beyond 0) qualify."""
    z = sign(direction) * np.asarray(y, float)
    out = []
    while len(out) < params.MAX_CAND and np.isfinite(z).any():
        j = int(np.nanargmax(z))
        if z[j] <= 0:
            break
        out.append(j)
        z[max(0, j - geo.sep + 1):j + geo.sep] = np.nan
    return out


def half_width(y: np.ndarray, j: int, direction: str) -> tuple[int, int]:
    """First and last bin of the contiguous run around j that reaches at least half of y[j] (in the direction);
    a NaN window ends the run."""
    z = sign(direction) * np.asarray(y, float)
    half = z[j] / 2
    lft = j
    while lft > 0 and np.isfinite(z[lft - 1]) and z[lft - 1] >= half:
        lft -= 1
    rgt = j
    while rgt < len(z) - 1 and np.isfinite(z[rgt + 1]) and z[rgt + 1] >= half:
        rgt += 1
    return lft, rgt


def robust_z(y: np.ndarray, j: int) -> float:
    """(y[j] - median) / (1.4826 x MAD) over the chromosome's finite windows."""
    f = y[np.isfinite(y)]
    if f.size < 3:
        return float("nan")
    med = np.median(f)
    mad = 1.4826 * np.median(np.abs(f - med))
    return float((y[j] - med) / mad) if mad > 0 else float("nan")


def signature_match(num: dict, vok: dict, sig: pd.DataFrame, names: list[str]) -> tuple[int, int, float]:
    """(signature entries with data, entries matched one by one, signature score). The score is the weighted
    least-squares slope of the candidate's class deviations on the anchors' mean deviations (weights = calls):
    1 = deviates like the anchors, 0 = like the background. An entry matches one by one when it deviates in its
    learned direction by at least SIG_FRACTION of the anchors' mean. Entries with < MIN_CALLS calls are left out."""
    idx = {n: i for i, n in enumerate(names)}
    with_data = matched = 0
    sxy = sxx = 0.0
    for r in sig.itertuples():
        if r.context not in vok or r.cls not in idx:
            continue
        k = idx[r.cls]
        w = vok[r.context][k]
        if w < params.MIN_CALLS:
            continue
        with_data += 1
        d = num[r.context][k] / w
        matched += bool(np.sign(d) == np.sign(r.mean_dev) and abs(d) >= params.SIG_FRACTION * abs(r.mean_dev))
        sxy += w * d * r.mean_dev
        sxx += w * r.mean_dev ** 2
    return with_data, matched, (sxy / sxx if sxx > 0 else float("nan"))


def table(picks, sig: pd.DataFrame, names: list[str], anchors: dict, geo, D: float,
          lengths: dict[str, int]) -> pd.DataFrame:
    """One row per pick, with every evidence column. anchors: chrom -> dict(name, bin, ...); lengths: chrom -> bp.
    Coordinates are clipped to the chromosome (the last bin is usually shorter than BIN); peak = middle of the
    pick's bin. near_end: the pick lies within geo.sep bins of a chromosome end, so its window may be truncated and
    its class reference ring is one-sided (a flag only: acrocentric centromeres sit there)."""
    b = params.BIN
    rows = []
    for p in picks:
        wd, mt, score = signature_match(p.num, p.vok, sig, names)
        label = "centromere_like" if wd > 0 and score >= params.SIG_FRACTION else "dip"
        a = anchors.get(p.chrom)
        near = a is not None and abs(p.bin - a["bin"]) < geo.sep
        L = lengths[p.chrom]
        w0, w1 = max(0, p.bin - (geo.w - 1) // 2), p.bin + geo.w // 2 + 1
        d0, d1 = p.left * b, min((p.right + 1) * b, L)
        row = dict(sample=p.sample, chrom=p.chrom, rank=p.rank, label=label,
                   peak=(p.bin * b + min((p.bin + 1) * b, L)) // 2,
                   window_start=w0 * b, window_end=min(w1 * b, L), dip_start=d0, dip_end=d1,
                   width=d1 - d0, value=p.value, robust_z=p.z,
                   near_end=bool(min(p.bin, n_bins(L) - 1 - p.bin) < geo.sep),
                   signature_with_data=wd, signature_matched=mt, signature_score=score,
                   plant_specific=p.ps, plant_specific_p=p.ps_p, plant_specific_n=p.ps_n,
                   sample_specific=bool(np.isfinite(p.ps_p) and p.ps_p <= params.SPECIFIC_P),
                   anchor=a["name"] if near else "",
                   anchor_distance=(p.bin - a["bin"]) * b if a is not None else np.nan, D=D)
        for ctx in p.vok:
            dev = np.divide(p.num[ctx], p.vok[ctx], out=np.full(len(names), np.nan),
                            where=p.vok[ctx] >= params.MIN_CALLS)
            for k, n in enumerate(names):
                row[f"{ctx}:{n}"] = dev[k]
        rows.append(row)
    t = pd.DataFrame(rows)
    if t.empty:
        return t
    dev_cols = [c for c in t.columns if ":" in c]
    keep = [c for c in dev_cols if t[c].notna().any()]
    return t.drop(columns=[c for c in dev_cols if c not in keep])


def anchored(scan, anchors: dict, cand: pd.DataFrame, direction: str, sig: pd.DataFrame,
             names: list[str]) -> pd.DataFrame:
    """Per anchor x sample: CHG adjusted value at the anchor, its rank and rank share among the chromosome's
    windows (1 = most extreme in the direction), the nearest candidate within 4 D, and the data there."""
    s = sign(direction)
    b = params.BIN
    best_sites = {}
    for x in scan.at_anchor:
        best_sites[x.chrom] = max(best_sites.get(x.chrom, 0.0), x.sites)
    rows = []
    for x in scan.at_anchor:
        a = anchors[x.chrom]
        y = scan.adj[x.sample][x.chrom]
        fin = np.isfinite(y)
        val = float(y[x.bin])
        rank = int((s * y[fin] > s * val).sum()) + 1 if np.isfinite(val) else np.nan
        wd, mt, score = signature_match(x.num, x.vok, sig, names)
        row = dict(anchor=a["name"], chrom=x.chrom, sample=x.sample, anchor_mid=a["mid"], value=val,
                   window_rank=rank, window_rank_share=rank / fin.sum() if np.isfinite(val) else np.nan,
                   n_windows=int(fin.sum()),
                   signature_with_data=wd, signature_matched=mt, signature_score=score, sites=x.sites, calls=x.calls,
                   sites_vs_best=x.sites / best_sites[x.chrom] if best_sites[x.chrom] > 0 else np.nan,
                   candidate_rank=np.nan, candidate_label="", candidate_distance=np.nan)
        if not cand.empty:
            c = cand[(cand["sample"] == x.sample) & (cand.chrom == x.chrom)]
            c = c[(c.peak // b - x.bin).abs() < scan.geo.sep]
            if len(c):
                c = c.iloc[int(np.argmin((c.peak // b - x.bin).abs().to_numpy()))]
                row.update(candidate_rank=int(c["rank"]), candidate_label=c.label,
                           candidate_distance=int(c.peak - (x.bin * b + b // 2)))
        rows.append(row)
    return pd.DataFrame(rows)
