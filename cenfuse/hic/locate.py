"""Prior-free location: telomere-pole score, start, iterated partner windows, peak + half-max core, Rabl test."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats

from . import params
from .score import cluster_score, end_bins, partner_oe, window_bins


def telomere_score(W, bins, good, mean, members):
    """Mean O/E of each bin with the other chromosomes' end windows (the telomere pole of the Rabl nucleus)."""
    return cluster_score(partner_oe(W, bins, good, mean, end_bins(bins, good, members)), bins, good, members)


def start_positions(T, bins, members):
    est = {}
    k = params.START_SMOOTH
    for c in np.where(members)[0]:
        v = pd.Series(T[bins.sl(int(c))]).rolling(k, center=True, min_periods=(k + 1) // 2).mean().to_numpy()
        if np.isfinite(v).any():
            est[int(c)] = int(np.nanargmin(v))
    return est


@dataclass
class Located:
    est: dict
    S: np.ndarray
    oe: dict
    rounds: int
    converged: bool
    cycle: int = 0          # length of the cycle the iteration fell into (0 = settled on a fixed point)
    spread: dict = field(default_factory=dict)   # per chromosome: max - min of its estimate over the cycle states
    # (0 at a fixed point; empty when not converged). `est` = the partner-window centres used to compute S (and oe).
    # Call peaks are always read from S (peak_core), never from est: on a chromosome that oscillates in a cycle,
    # est can differ from argmax(S).


def _total(S, bins, est):
    """Sum over chromosomes of the peak height (max of S) -- how well-defined the cluster is in this state."""
    return float(sum(np.nanmax(v) for v in (S[bins.sl(k)] for k in est) if np.isfinite(v).any()))


def locate(W, bins, good, mean, members, start):
    """Iterate: partner windows = estimates +/- HALF_WIN; each estimate = top bin of the cluster score. Stops at a
    fixed point, or when a state repeats (a cycle of states): then the cycle state with the highest summed peak
    height (max of S per chromosome) is kept (converged, cycle = its length). In a cycle most chromosomes usually
    alternate between neighbouring bins of a flat-topped peak, but a chromosome without a clear cluster can jump
    between distant bins; `spread` records, per chromosome, how far its estimate moves within the cycle."""
    est = dict(start)
    seen, states = {}, []
    for r in range(1, params.MAX_ROUNDS + 1):
        oe = partner_oe(W, bins, good, mean, window_bins(bins, good, est))
        S = cluster_score(oe, bins, good, members)
        key = tuple(sorted(est.items()))
        seen[key] = len(states)
        states.append((dict(est), S, oe))
        new = {}
        for k in est:
            v = S[bins.sl(k)]
            new[k] = int(np.nanargmax(v)) if np.isfinite(v).any() else est[k]
        if new == est:
            return Located(est, S, oe, r, True, spread={k: 0 for k in est})
        nkey = tuple(sorted(new.items()))
        if nkey in seen:
            cyc = states[seen[nkey]:]
            best = max(cyc, key=lambda s: _total(s[1], bins, s[0]))
            spread = {k: max(s[0][k] for s in cyc) - min(s[0][k] for s in cyc) for k in est}
            return Located(best[0], best[1], best[2], r, True, cycle=len(cyc), spread=spread)
        est = new
    last_est, last_S, last_oe = states[-1]      # the state whose S/oe were computed, not the advanced `new`
    return Located(last_est, last_S, last_oe, params.MAX_ROUNDS, False)


def peak_core(S, bins, k):
    """Peak bin and the unbroken run of bins >= median + (peak - median)/2 around it (NaN breaks the run)."""
    v = S[bins.sl(k)]
    ok = np.isfinite(v)
    if not ok.any():
        return None
    j = int(np.nanargmax(v))
    med = float(np.nanmedian(v))
    half = med + (v[j] - med) / 2
    lo = hi = j
    while lo - 1 >= 0 and ok[lo - 1] and v[lo - 1] >= half:
        lo -= 1
    while hi + 1 < len(v) and ok[hi + 1] and v[hi + 1] >= half:
        hi += 1
    return dict(peak=j, core_lo=lo, core_hi=hi, peak_score=float(v[j]), median=med,
                ratio=float(v[j] / med) if med > 0 else float("nan"))


def rabl_test(T, bins, members):
    """Per chromosome: Spearman rho of the telomere-pole score vs distance from the nearest end (Rabl: rho < 0).
    One-sided sign test across chromosomes; p_floor = the smallest P this many chromosomes can give."""
    rows = []
    for k in np.where(members)[0]:
        v = T[bins.sl(int(k))]
        i = np.arange(len(v))
        dist = np.minimum(i, len(v) - 1 - i)
        ok = np.isfinite(v)
        # a constant score has no rank correlation: dropped up front (scipy would warn and return NaN)
        if ok.sum() >= params.RABL_MIN_BINS and np.ptp(v[ok]) > 0:
            rho = stats.spearmanr(dist[ok], v[ok]).correlation
            if np.isfinite(rho):
                rows.append(dict(chrom=bins.names[int(k)], rho=rho))
    t = pd.DataFrame(rows, columns=["chrom", "rho"])
    n = len(t)
    neg = int((t.rho < 0).sum())
    p = stats.binomtest(neg, n, 0.5, alternative="greater").pvalue if n else 1.0
    return dict(n=n, negative=neg, p=float(p), passed=bool(p <= params.RABL_P), p_floor=0.5 ** n if n else 1.0,
                per_chrom=t)
