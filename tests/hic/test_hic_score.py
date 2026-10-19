from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from cenfuse.hic.balance import balanced, ice
from cenfuse.hic.bins import Bins
from cenfuse.hic.score import block_means, cluster_score, partner_oe, window_bins


def dense_reference(Md, bins, centres, half=4):
    """Dense reference implementation (ICE, per-block O/E, window means, mean over partners)."""
    cov = Md.sum(1)
    good = cov >= 0.2 * np.median(cov[cov > 0])
    W = Md.astype(float).copy()
    W[~good] = 0
    W[:, ~good] = 0
    for _ in range(200):
        s = W.sum(1)
        s = s / s[good].mean()
        s[~good] = 1
        W /= np.outer(s, s)
        if np.abs(s[good] - 1).max() < 1e-3:
            break
    W[~good] = np.nan
    W[:, ~good] = np.nan
    K = len(bins.names)
    OE = np.full_like(W, np.nan)
    for x in range(K):
        for y in range(K):
            if x != y:
                blk = W[bins.sl(x), bins.sl(y)]
                OE[bins.sl(x), bins.sl(y)] = blk / np.nanmean(blk)
    win = {}
    for y, c in centres.items():
        idx = bins.off[y] + np.arange(max(c - half, 0), min(c + half, bins.nb[y] - 1) + 1)
        win[y] = idx[good[idx]]
    S = np.full(bins.N, np.nan)
    for i in range(bins.N):
        x = bins.chrom_of[i]
        if good[i]:
            S[i] = np.nanmean([np.nanmean(OE[i, win[y]]) for y in range(K) if y != x])
    return S, good


def random_matrix(bins, seed=0):
    rng = np.random.default_rng(seed)
    A = rng.poisson(5, (bins.N, bins.N)).astype(float)
    A = np.triu(A) + np.triu(A, 1).T
    A[3, :] = 0
    A[:, 3] = 0                        # one empty bin -> masked
    return A


def test_ice_rows_equal():
    b = Bins({"a": 2000, "b": 1500, "c": 1800}, 100)
    A = random_matrix(b)
    bias, good, it, conv = ice(sp.csr_matrix(A))
    W = balanced(sp.csr_matrix(A), bias, good).toarray()
    rs = W.sum(1)[good]
    assert conv and not good[3] and np.allclose(rs, rs.mean(), rtol=2e-3)


def test_cluster_score_matches_dense_reference():
    b = Bins({"a": 2000, "b": 1500, "c": 1800}, 100)
    A = random_matrix(b, seed=3)
    centres = {0: 8, 1: 5, 2: 12}
    S_ref, good_ref = dense_reference(A, b, centres)
    M = sp.csr_matrix(A)
    bias, good, _, _ = ice(M)
    W = balanced(M, bias, good)
    members = np.ones(3, bool)
    mean = block_means(W, b, good, members)
    S = cluster_score(partner_oe(W, b, good, mean, window_bins(b, good, centres)), b, good, members)
    assert (good == good_ref).all()
    ok = np.isfinite(S_ref)
    assert np.isfinite(S[ok]).all() and np.allclose(S[ok], S_ref[ok], rtol=1e-6)


def test_exclude_partner_changes_only_that_chromosome():
    b = Bins({"a": 2000, "b": 1500, "c": 1800}, 100)
    M = sp.csr_matrix(random_matrix(b, seed=5))
    bias, good, _, _ = ice(M)
    W = balanced(M, bias, good)
    members = np.ones(3, bool)
    oe = partner_oe(W, b, good, block_means(W, b, good, members), window_bins(b, good, {0: 8, 1: 5, 2: 12}))
    S = cluster_score(oe, b, good, members)
    Sx = cluster_score(oe, b, good, members, exclude={0: {1}})
    assert np.allclose(S[b.sl(1)], Sx[b.sl(1)], equal_nan=True)
    assert np.allclose(Sx[b.sl(0)], oe[2][b.sl(0)], equal_nan=True)


def test_partner_oe_never_inf_and_masks_bad_bins():
    """A zero or NaN block mean gives NaN (never inf), masked bins and the partner's own bins give NaN, and the
    values elsewhere are window mean / block mean; no warnings, W is not converted."""
    import warnings
    b = Bins({"a": 2000, "b": 1500, "c": 1800}, 100)
    M = sp.csr_matrix(random_matrix(b, seed=7))
    bias, good, _, _ = ice(M)
    W = balanced(M, bias, good)
    mean = block_means(W, b, good, np.ones(3, bool))
    mean[0, 1] = 0.0                                   # block a-b: mean 0 (s > 0 there, so s / 0 would be inf)
    mean[2, 1] = np.nan                                # block c-b: undefined
    win = window_bins(b, good, {1: 5})
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        v = partner_oe(W, b, good, mean, win)[1]
    assert not np.isinf(v).any() and np.isnan(v).all()
    mean[0, 1] = 2.0
    v = partner_oe(W, b, good, mean, win)[1]
    Wd = W.toarray()
    ref = Wd[:, win[1]].mean(axis=1) / 2.0
    a = b.sl(0)
    assert np.isnan(v[3]) and not good[3]               # bin 3 is empty -> masked
    ok = good[a]
    assert np.allclose(v[a][ok], ref[a][ok], rtol=1e-12) and np.isnan(v[b.sl(1)]).all()
