"""ICE (iterative correction) on a symmetric sparse count matrix."""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from . import params


def ice(M, min_frac=params.ICE_MIN_FRAC, tol=params.ICE_TOL, max_iter=params.ICE_MAX_ITER):
    """Return (bias, good, n_iter, converged). Balanced contact W_ij = M_ij / (bias_i * bias_j) for good bins."""
    M = M.tocsr().astype(np.float64, copy=False)
    cov = np.asarray(M.sum(axis=1)).ravel()
    pos = cov[cov > 0]
    if len(pos) == 0:
        raise ValueError("empty contact matrix")
    good = cov >= min_frac * np.median(pos)
    g = sp.diags(good.astype(float))
    Mg = (g @ M @ g).tocsr()
    good &= np.asarray(Mg.sum(axis=1)).ravel() > 0
    bias = np.ones(M.shape[0])
    for it in range(1, max_iter + 1):
        inv = np.where(good, 1.0 / bias, 0.0)
        s = (Mg @ inv) * inv
        s = s / s[good].mean()
        s[~good] = 1.0
        bias *= s
        if np.abs(s[good] - 1).max() < tol:
            return bias, good, it, True
    return bias, good, max_iter, False


def balanced(M, bias, good):
    inv = sp.diags(np.where(good, 1.0 / bias, 0.0))
    return (inv @ M.tocsr().astype(np.float64, copy=False) @ inv).tocsr()
