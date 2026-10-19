"""Trans observed/expected per chromosome pair and the centromere-cluster score (sparse matrices)."""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from . import params


def block_means(W, bins, good, members):
    """Mean balanced contact of each chromosome-pair block over good x good bins (zeros included)."""
    K = len(bins.names)
    C = sp.csr_matrix((good.astype(float), (np.arange(bins.N), bins.chrom_of)), shape=(bins.N, K))
    S = (C.T @ W @ C).toarray()
    ng = np.asarray(C.sum(axis=0)).ravel()
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = S / np.outer(ng, ng)
    mean[~members, :] = np.nan
    mean[:, ~members] = np.nan
    np.fill_diagonal(mean, np.nan)
    return mean


def window_bins(bins, good, centres, half=params.HALF_WIN):
    out = {}
    for k, c in centres.items():
        lo, hi = max(int(c) - half, 0), min(int(c) + half, int(bins.nb[k]) - 1)
        idx = bins.off[k] + np.arange(lo, hi + 1)
        out[int(k)] = idx[good[idx]]
    return out


def end_bins(bins, good, members, n=params.END_WIN):
    out = {}
    for k in np.where(members)[0]:
        nb = int(bins.nb[k])
        loc = np.unique(np.concatenate([np.arange(min(n, nb)), np.arange(max(nb - n, 0), nb)]))
        idx = bins.off[k] + loc
        out[int(k)] = idx[good[idx]]
    return out


def partner_oe(W, bins, good, mean, windows):
    """For each partner y: O/E of every bin with y's window = mean balanced contact with the window / block mean.
    W is symmetric, so a bin's contacts with the window are the column sums of the window's CSR rows (no copy or
    conversion of W). Masked bins (~good), y's own bins and blocks with no mean > 0 get NaN."""
    W = W.tocsr()                       # no-op (no copy) for the CSR matrices built here
    out = {}
    for y, w in windows.items():
        if len(w) == 0:
            continue
        s = np.asarray(W[w].sum(axis=0)).ravel() / len(w)
        m = mean[bins.chrom_of, y]
        v = np.full(bins.N, np.nan)
        np.divide(s, m, out=v, where=m > 0)        # NaN block mean compares False: stays NaN, no warning
        v[bins.chrom_of == y] = np.nan
        v[~good] = np.nan
        out[int(y)] = v
    return out


def cluster_score(oe, bins, good, members, exclude=None):
    """Mean over partners y != x of the O/E with y's window; exclude = {x: {y, ...}} leaves partners out for x.
    NaN O/E values are skipped (as np.nanmean); a bin with no finite partner value gets NaN. Computed as sum / count
    without the warnings module, which is not thread-safe (the bootstrap calls this from several threads)."""
    S = np.full(bins.N, np.nan)
    for x in np.where(members)[0]:
        ys = [y for y in oe if y != x and not (exclude and y in exclude.get(int(x), ()))]
        if not ys:
            continue
        sl = bins.sl(int(x))
        A = np.vstack([oe[y][sl] for y in ys])
        ok = ~np.isnan(A)
        n = ok.sum(axis=0)
        np.divide(np.where(ok, A, 0.0).sum(axis=0), n, out=S[sl], where=n > 0)
    S[~good] = np.nan
    return S
