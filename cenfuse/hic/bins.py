"""Genome bins, sparse contact matrices, coarsening and the data-derived bin-size rule."""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

from . import params
from .io import iter_pairs


class Bins:
    """Fixed-size bins per chromosome, starting at 0 on each chromosome (last bin may be short)."""

    def __init__(self, chroms: dict, size: int):
        self.names = list(chroms)
        self.lengths = np.array([chroms[c] for c in self.names], dtype=np.int64)
        self.size = int(size)
        self.nb = -(-self.lengths // self.size)
        self.off = np.concatenate([[0], np.cumsum(self.nb)[:-1]]).astype(np.int64)
        self.N = int(self.nb.sum())
        self.chrom_of = np.repeat(np.arange(len(self.names)), self.nb)
        self.local = np.arange(self.N) - self.off[self.chrom_of]
        self.start = self.local * self.size
        self.end = np.minimum(self.start + self.size, self.lengths[self.chrom_of])

    def sl(self, k: int) -> slice:
        return slice(int(self.off[k]), int(self.off[k] + self.nb[k]))

    def glob(self, ci, pos):
        return self.off[np.asarray(ci)] + np.asarray(pos) // self.size


def symmetrize(M):
    M = M.tocsr()
    return (M + M.T - sp.diags(M.diagonal())).tocsr()


def _chunk_matrix(bins: Bins, ci1, p1, ci2, p2):
    """Upper-triangular count matrix of one chunk of pairs, and the number of pairs off the chromosomes."""
    L = bins.lengths
    ok = (p1 >= 0) & (p2 >= 0) & (p1 < L[ci1]) & (p2 < L[ci2])
    a, b = bins.glob(ci1[ok], p1[ok]), bins.glob(ci2[ok], p2[ok])
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    m = sp.coo_matrix((np.ones(len(lo)), (lo, hi)), shape=(bins.N, bins.N)).tocsr()
    return m, int((~ok).sum())


def add_pairs(acc, bins: Bins, ci1, p1, ci2, p2):
    m, n_out = _chunk_matrix(bins, ci1, p1, ci2, p2)
    return (m if acc is None else acc + m), n_out


def _tree_sum(mats):
    """Pairwise (tree) sum of a non-empty list of equal-shape sparse matrices."""
    mats = list(mats)
    while len(mats) > 1:
        nxt = [mats[i] + mats[i + 1] for i in range(0, len(mats) - 1, 2)]
        if len(mats) % 2:
            nxt.append(mats[-1])
        mats = nxt
    return mats[0]


def matrix_from_pairs(path, chroms: dict, bins: Bins, stats, chunksize: int = params.PAIRS_CHUNK):
    """Symmetric float64 count matrix of all kept pairs, and the number of pairs dropped as off-chromosome.
    Per-chunk matrices are collected and tree-summed (collapsing every params.FANIN matrices into one), so
    the cost is not quadratic in the number of chunks and memory stays bounded."""
    parts, n_out = [], 0
    for ci1, p1, ci2, p2 in iter_pairs(path, chroms, stats, chunksize=chunksize):
        m, k = _chunk_matrix(bins, ci1, p1, ci2, p2)
        parts.append(m)
        n_out += k
        if len(parts) >= params.FANIN:
            parts = [_tree_sum(parts)]
    if not parts:
        return sp.csr_matrix((bins.N, bins.N), dtype=np.float64), n_out
    return symmetrize(_tree_sum(parts)).astype(np.float64), n_out


def coarsen(M, fine: Bins, coarse: Bins):
    """Sum a symmetric fine-bin matrix into coarse bins. Only the upper triangle is summed, so each fine pair is
    counted once (also when both ends fall in one coarse bin); the result equals binning the pairs directly."""
    target = coarse.off[fine.chrom_of] + fine.start // coarse.size
    P = sp.csr_matrix((np.ones(fine.N), (np.arange(fine.N), target)), shape=(fine.N, coarse.N))
    U = sp.triu(M.tocsr(), k=0).tocsr()
    return symmetrize(P.T @ U @ P)


def expected_cen_contacts(M, bins: Bins, members: np.ndarray) -> float:
    """Median over member bins of: trans contacts with other member chromosomes x (summed partner-window length /
    summed length of the other member chromosomes); windows = 2 * HALF_WIN + 1 bins."""
    mem = np.where(members)[0]
    if len(mem) < 2:
        return 0.0
    col = np.full(len(bins.names), -1)
    col[mem] = np.arange(len(mem))
    rows = np.where(members[bins.chrom_of])[0]
    C = sp.csr_matrix((np.ones(len(rows)), (rows, col[bins.chrom_of[rows]])), shape=(bins.N, len(mem)))
    R = np.asarray((M @ C).todense())                       # N x (member chromosomes): contacts per chromosome
    win_bp = (2 * params.HALF_WIN + 1) * bins.size
    total_len = bins.lengths[mem].sum()
    vals = []
    for a, x in enumerate(mem):
        Rx = R[bins.sl(int(x))]
        trans = Rx.sum(axis=1) - Rx[:, a]                   # drop the bin's own chromosome
        share = (len(mem) - 1) * win_bp / (total_len - bins.lengths[x])
        vals.append(trans * share)
    return float(np.median(np.concatenate(vals)))


def choose_bin_size(get_matrix, sizes):
    rows = []
    for size in sizes:
        M, b = get_matrix(size)
        members = b.nb >= params.MIN_BINS
        n_call = int(members.sum())
        e = expected_cen_contacts(M, b, members) if n_call >= params.MIN_CALLABLE else 0.0
        rows.append(dict(bin_size=size, n_callable=n_call, expected_cen_contacts=e,
                         meets_rule=bool(n_call >= params.MIN_CALLABLE and e >= params.MIN_EXPECTED_CEN)))
    t = pd.DataFrame(rows)
    ok = t[t.meets_rule]
    usable = t[t.n_callable >= params.MIN_CALLABLE]
    if len(ok):
        size = int(ok.bin_size.min())
    elif len(usable):
        size = int(usable.bin_size.max())
    else:
        size = None
    t["chosen"] = t.bin_size == size
    return size, t
