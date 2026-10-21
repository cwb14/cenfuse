"""Calibration with anchors: which sequence classes mark the known centromeres by their methylation (the species'
signature) and in which direction. A class enters the signature when it is (1) methylated differently at the
anchors: tested against N_PERM random anchor placements on windows where that class has data, Benjamini-Hochberg
across class x context; (2) concentrated at the anchors: its genome-wide share of calls is at most SIG_FRACTION
of its share at the anchors. A class found everywhere carries every dip, true or false, so it cannot mark one; and
(3) detectable: its mean deviation at the anchors is at least MIN_EFFECT, so one candidate window can tell it from
zero."""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from . import params

log = logging.getLogger("cenfuse-meth")
ALL = "all"   # pseudo-class: the adjusted track itself (all classes together)


def anchor_table(anchors: pd.DataFrame, chroms: dict[str, int]) -> dict:
    """One anchor per analysed chromosome: chrom -> dict(name, start, end, mid, bin). Several intervals on one
    chromosome are merged into their span (warned); anchors on sequences that are not analysed are dropped (warned)."""
    out = {}
    skipped = sorted(set(anchors.chrom) - set(chroms))
    if skipped:
        log.warning("anchors on %d sequences not analysed (too short or unplaced) ignored: %s", len(skipped),
                    ", ".join(skipped[:5]))
    for c, g in anchors[anchors.chrom.isin(set(chroms))].groupby("chrom", sort=False):
        if len(g) > 1:
            s, e = int(g.start.min()), int(g.end.max())
            mid, name = (s + e) // 2, f"{c}_anchor"
            log.warning("%d anchor intervals on %s merged into one span %d-%d", len(g), c, s, e)
        else:
            r = g.iloc[0]
            s, e, mid, name = int(r.start), int(r.end), int(r.mid), str(r["name"])
        mid = min(max(mid, 0), chroms[c] - 1)
        out[c] = dict(name=name, start=s, end=e, mid=mid, bin=mid // params.BIN)
    return out


def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q values (NaN p stays NaN)."""
    q = np.full(len(p), np.nan)
    f = np.flatnonzero(np.isfinite(p))
    if f.size == 0:
        return q
    o = f[np.argsort(p[f])]
    m = f.size
    adj = p[o] * m / np.arange(1, m + 1)
    q[o] = np.minimum(1.0, np.minimum.accumulate(adj[::-1])[::-1])
    return q


def perm_floor_ok(m: int) -> bool:
    """Can a lone signal at the smallest possible P, 1 / (N_PERM + 1), pass Benjamini-Hochberg over m tests?"""
    return m / (params.N_PERM + 1) <= params.FDR


def calibrate(pooled: dict, totals: dict, names: list[str], seed: int = params.SEED) -> pd.DataFrame:
    """pooled: context -> chrom -> scan.Pooled; totals: context -> calls per class (scan.Scan.totals). One row per
    context x class (plus 'all'), with the deviation test and the enrichment at the anchors.
    A class is tested on the anchored chromosomes where it has data at the anchor and somewhere else; each of the
    N_PERM placements puts every such chromosome's anchor on a random window where the class has data (>= 4 D from
    the anchor and the ends). P = (1 + placements whose mean deviation is at least as extreme) / (1 + N_PERM)."""
    rng = np.random.default_rng(seed)
    rows = []
    for ctx in sorted(pooled, key=lambda x: (x != "CHG", x)):
        per = pooled[ctx]
        tot = np.append(totals[ctx], totals[ctx].sum())
        shares = [per[c].anchor_vok / per[c].anchor_vok[-1] for c in sorted(per) if per[c].anchor_vok[-1] > 0]
        anchor_share = np.mean(shares, axis=0) if shares else np.full(len(tot), np.nan)
        for k, cls in enumerate(names + [ALL]):
            at, tested, null = [], [], []
            for c in sorted(per):
                pc = per[c]
                if pc.anchor_vok[k] < params.MIN_CALLS:
                    continue
                d = pc.anchor_num[k] / pc.anchor_vok[k]
                at.append(d)
                if k in pc.elsewhere:
                    en, ev = pc.elsewhere[k]
                    i = rng.integers(0, len(en), params.N_PERM)
                    tested.append(d)
                    null.append(en[i] / ev[i])
            p = np.nan
            if len(tested) >= params.MIN_ANCHORS:
                obs, nul = np.mean(tested), np.mean(null, axis=0)
                p = (1 + int((np.abs(nul) >= abs(obs)).sum())) / (1 + params.N_PERM)
            gs = tot[k] / tot[-1] if tot[-1] > 0 else np.nan
            rows.append(dict(context=ctx, cls=cls, anchors_with_data=len(at),
                             mean_dev=float(np.mean(at)) if at else np.nan, chromosomes_tested=len(tested), p=p,
                             anchor_share=anchor_share[k], genome_share=gs))
    t = pd.DataFrame(rows, columns=["context", "cls", "anchors_with_data", "mean_dev", "chromosomes_tested", "p",
                                    "anchor_share", "genome_share"])
    t["q"] = bh(t.p.to_numpy(float))
    m = int(t.p.notna().sum())
    if m and not perm_floor_ok(m):
        log.warning("too many class x context tests (%d) for %d permutations: a lone signal cannot reach the FDR",
                    m, params.N_PERM)
    t["direction"] = np.where(t.mean_dev < 0, "down", np.where(t.mean_dev > 0, "up", ""))
    t["enriched"] = t.genome_share <= params.SIG_FRACTION * t.anchor_share
    t["detectable"] = t.mean_dev.abs() >= params.MIN_EFFECT
    t["in_signature"] = (t.q <= params.FDR) & t.enriched & t.detectable & (t.cls != ALL)
    return t.sort_values(["in_signature", "p"], ascending=[False, True], kind="stable").reset_index(drop=True)


def signature(cal: pd.DataFrame | None) -> pd.DataFrame:
    if cal is None:
        return pd.DataFrame(columns=["context", "cls", "mean_dev", "direction"])
    return cal[cal.in_signature][["context", "cls", "mean_dev", "direction"]].reset_index(drop=True)


def learned_direction(cal: pd.DataFrame | None) -> tuple[str, str]:
    """(direction, source): the CHG adjusted track's direction at the anchors when it passes the FDR, else 'down'
    (the default; plant CENH3 chromatin usually loses CHG)."""
    if cal is not None:
        r = cal[(cal.context == "CHG") & (cal.cls == ALL) & (cal.q <= params.FDR)]
        if len(r):
            return str(r.direction.iloc[0]), "anchors"
    return "down", "default"
