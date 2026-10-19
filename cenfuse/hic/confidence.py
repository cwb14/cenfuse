"""How sure is each call: Poisson bootstrap, leave one partner out, position-aware null, confidence levels."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor

import numpy as np
import scipy.sparse as sp

from . import params
from .balance import balanced, ice
from .bins import symmetrize
from .score import block_means, cluster_score, partner_oe, window_bins


def poisson_resample(M, rng):
    U = sp.triu(M).tocsr()
    R = sp.csr_matrix((rng.poisson(U.data).astype(np.float64), U.indices, U.indptr), shape=U.shape)
    return symmetrize(R)


def _score_fixed(M, bins, members, est):
    bias, good, _, _ = ice(M)
    W = balanced(M, bias, good)
    mean = block_means(W, bins, good, members)
    return cluster_score(partner_oe(W, bins, good, mean, window_bins(bins, good, est)), bins, good, members)


def bootstrap(M, bins, members, est, n=params.N_BOOT, seed=params.SEED, threads=1):
    """Peak bin of every chromosome in n Poisson resamples (re-balanced, re-scored; partner windows fixed at est)."""
    seeds = np.random.SeedSequence(seed).spawn(n)

    def one(ss):
        S = _score_fixed(poisson_resample(M, np.random.default_rng(ss)), bins, members, est)
        return {k: int(np.nanargmax(S[bins.sl(k)])) if np.isfinite(S[bins.sl(k)]).any() else -1 for k in est}

    with ThreadPoolExecutor(max(threads, 1)) as ex:
        reps = list(ex.map(one, seeds))
    return {k: np.array([r[k] for r in reps]) for k in est}


def leave_one_out(oe, bins, good, members, est):
    """Peak bin of each chromosome with each partner left out in turn."""
    out = {}
    for x in est:
        peaks = []
        for y in oe:
            if y == x:
                continue
            v = cluster_score(oe, bins, good, members, exclude={x: {y}})[bins.sl(x)]
            peaks.append(int(np.nanargmax(v)) if np.isfinite(v).any() else -1)
        out[x] = np.array(peaks)
    return out


def null_test(W, bins, good, mean, members, est, peaks, ratios, n=params.N_NULL, seed=params.SEED):
    """P = (1 + #replicates where the max within NULL_REACH bins of the call, over the chromosome median, reaches the
    observed ratio) / (1 + n). Partner windows go to random bins >= NULL_EXCL from that partner's estimate and from
    its ends; partners without such bins sit out. Fewer than 2 such partners, a non-finite observed
    ratio, or no finite score in any replicate: P = NaN (untestable)."""
    rng = np.random.default_rng(seed)
    elig = {}
    for y, c in est.items():
        loc = np.arange(params.NULL_EXCL, int(bins.nb[y]) - params.NULL_EXCL)
        loc = loc[np.abs(loc - c) >= params.NULL_EXCL]
        if len(loc):
            elig[y] = loc
    if len(elig) < 2:
        return {x: float("nan") for x in est}
    hits = {x: 0 for x in est}
    seen = {x: 0 for x in est}      # replicates in which x had a finite score
    for _ in range(n):
        cent = {y: int(rng.choice(loc)) for y, loc in elig.items()}
        S = cluster_score(partner_oe(W, bins, good, mean, window_bins(bins, good, cent)), bins, good, members)
        for x in est:
            v = S[bins.sl(x)]
            if not np.isfinite(v).any():
                continue
            seen[x] += 1
            med = np.nanmedian(v)
            lo, hi = max(peaks[x] - params.NULL_REACH, 0), min(peaks[x] + params.NULL_REACH + 1, len(v))
            win = v[lo:hi]
            if med > 0 and np.isfinite(win).any() and np.nanmax(win) / med >= ratios[x]:
                hits[x] += 1
    # a non-finite observed ratio or a never-scored chromosome cannot be tested: NaN, not the best possible P
    return {x: (1 + hits[x]) / (1 + n) if np.isfinite(ratios[x]) and seen[x] else float("nan") for x in est}


def level(boot_in_core: float, lopo_stable: bool, p: float) -> str:
    """high = all three; medium = two; low = one (or none while P <= 0.05); none = P > 0.05 (no call).
    Untestable null (NaN): judged on the other two, capped at medium (flag null_untestable)."""
    a, b = boot_in_core >= params.BOOT_IN_CORE, bool(lopo_stable)
    if not np.isfinite(p):
        return "medium" if a and b else "low"
    if p > params.P_NONE:
        return "none"
    k = int(a) + int(b) + int(p <= params.P_HIGH)
    return {3: "high", 2: "medium"}.get(k, "low")
