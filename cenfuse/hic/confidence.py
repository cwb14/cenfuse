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


class _Resampler:
    """poisson_resample(M, rng) for one M, many times: the upper triangle and the symmetric layout are made once, and each
    replicate's draws (the same, in the same order) are written straight into that layout, zeros dropped as the sparse sum
    drops them, instead of a transpose and two sparse sums per replicate (0.14, future.md #62). The layout is checked against
    the sparse sum's once; if they ever differ, poisson_resample itself runs."""

    def __init__(self, M):
        self.M, U = M, sp.triu(M).tocsr()
        n = U.shape[0]
        rows = np.repeat(np.arange(n, dtype=np.int64), np.diff(U.indptr))
        cols = U.indices.astype(np.int64)
        off = rows != cols
        r_all, c_all = np.concatenate([rows, cols[off]]), np.concatenate([cols, rows[off]])
        order = np.lexsort((c_all, r_all))
        ref = symmetrize(sp.csr_matrix((np.ones_like(U.data), U.indices, U.indptr), shape=U.shape))
        self.U, self.shape = U, U.shape
        self.src = np.concatenate([np.arange(len(rows)), np.flatnonzero(off)])[order]
        self.indices = c_all[order].astype(ref.indices.dtype)
        self.indptr = np.concatenate([[0], np.cumsum(np.bincount(r_all, minlength=n))]).astype(ref.indptr.dtype)
        self.exact = bool(np.array_equal(self.indices, ref.indices) and np.array_equal(self.indptr, ref.indptr))

    def __call__(self, rng):
        if not self.exact:
            return poisson_resample(self.M, rng)
        R = sp.csr_matrix((rng.poisson(self.U.data).astype(np.float64)[self.src], self.indices.copy(), self.indptr.copy()),
                          shape=self.shape)
        R.eliminate_zeros()
        return R


def _peaks(S, bins, keys) -> dict:
    return {k: int(np.nanargmax(S[bins.sl(k)])) if np.isfinite(S[bins.sl(k)]).any() else -1 for k in keys}


def bootstrap(M, bins, members, est, n=params.N_BOOT, seed=params.SEED, threads=1, extra=None):
    """Peak bin of every chromosome in n Poisson resamples (re-balanced, re-scored; partner windows fixed at est).
    extra: {name: {chromosome: partner window bins}}, other fixed windows scored on the same resamples (cenfuse's anchored
    test, whose resamples, from the same seed and matrix, are these): then returns (peaks, {name: peak bins of every callable
    chromosome}). A resample, its balancing and its block means are about 97% of a replicate's cost, so the windows share
    them; the resamples are built by _Resampler (0.14, future.md #62; every number unchanged)."""
    seeds = np.random.SeedSequence(seed).spawn(n)
    resample = _Resampler(M)
    ks = [int(k) for k in np.where(members)[0]]

    def one(ss):
        R = resample(np.random.default_rng(ss))
        bias, good, _, _ = ice(R)
        W = balanced(R, bias, good)
        mean = block_means(W, bins, good, members)
        S = cluster_score(partner_oe(W, bins, good, mean, window_bins(bins, good, est)), bins, good, members)
        ext = {name: _peaks(cluster_score(partner_oe(W, bins, good, mean, {y: w[good[w]] for y, w in win.items()}), bins, good,
                                          members), bins, ks) for name, win in (extra or {}).items()}
        return _peaks(S, bins, est), ext

    with ThreadPoolExecutor(max(threads, 1)) as ex:
        reps = list(ex.map(one, seeds))
    out = {k: np.array([r[0][k] for r in reps]) for k in est}
    if extra is None:
        return out
    return out, {name: {k: np.array([r[1][name][k] for r in reps]) for k in ks} for name in extra}


def leave_one_out(oe, bins, good, members, est):
    """Peak bin of each chromosome with each partner left out in turn: the arithmetic of cluster_score(exclude={x: {y}}) on x's
    own bins only. Scoring every chromosome for every pair, as before 0.14, is cubic in the partners: 80 chromosomes took 81 s
    against 1 s now, with the same peaks (future.md #62; cenfuse's anchored test uses this function too, hicfuse.lopo_all)."""
    out = {}
    for x in est:
        sl = bins.sl(x)
        g = good[sl]
        ys = [y for y in oe if y != x]
        peaks = []
        for y in ys:
            rest = [z for z in ys if z != y]
            if not members[x] or not rest:              # cluster_score scores members only, and needs a partner
                peaks.append(-1)
                continue
            A = np.vstack([oe[z][sl] for z in rest])
            ok = ~np.isnan(A)
            n = ok.sum(axis=0)
            v = np.full(A.shape[1], np.nan)
            np.divide(np.where(ok, A, 0.0).sum(axis=0), n, out=v, where=n > 0)
            v[~g] = np.nan
            peaks.append(int(np.nanargmax(v)) if np.isfinite(v).any() else -1)
        out[x] = np.array(peaks)
    return out


DEGENERATE = -1.0        # P of a null that ran but could not score a chromosome (a non-finite observed ratio, or no finite
                         # score in any replicate), told apart from NaN (fewer than 2 partners: it could not run at all;
                         # future.md #69); not a probability, so no test reads it as one


def null_test(W, bins, good, mean, members, est, peaks, ratios, n=params.N_NULL, seed=params.SEED):
    """P = (1 + #replicates where the max within NULL_REACH bins of the call, over the chromosome median, reaches the
    observed ratio) / (1 + n). Partner windows go to random bins >= NULL_EXCL from that partner's estimate and from
    its ends; partners without such bins sit out. Fewer than 2 such partners: P = NaN (untestable); a non-finite observed
    ratio, or no finite score in any replicate: DEGENERATE (the null ran but says nothing here)."""
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
    # a non-finite observed ratio or a never-scored chromosome cannot be tested: DEGENERATE, not the best possible P
    return {x: (1 + hits[x]) / (1 + n) if np.isfinite(ratios[x]) and seen[x] else DEGENERATE for x in est}


def level(boot_in_core: float, lopo_stable: bool, p: float) -> str:
    """high = all three; medium = two; low = one (or none while P <= 0.05); none = P > 0.05 (no call).
    Untestable null (NaN; DEGENERATE): judged on the other two, capped at medium (flags null_untestable; null_degenerate)."""
    a, b = boot_in_core >= params.BOOT_IN_CORE, bool(lopo_stable)
    if not np.isfinite(p) or p < 0:
        return "medium" if a and b else "low"
    if p > params.P_NONE:
        return "none"
    k = int(a) + int(b) + int(p <= params.P_HIGH)
    return {3: "high", 2: "medium"}.get(k, "low")
