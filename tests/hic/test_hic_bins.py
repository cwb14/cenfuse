from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from cenfuse.hic import params
from cenfuse.hic.bins import (Bins, add_pairs, choose_bin_size, coarsen, expected_cen_contacts, matrix_from_pairs,
                         symmetrize)
from cenfuse.hic.io import PairStats


def test_bins_layout():
    b = Bins({"a": 1_000, "b": 450}, 200)
    assert list(b.nb) == [5, 3] and b.N == 8 and list(b.off) == [0, 5]
    assert list(b.chrom_of) == [0] * 5 + [1] * 3 and b.end[-1] == 450 and b.sl(1) == slice(5, 8)
    assert list(b.glob(np.array([1]), np.array([399]))) == [6]


def test_add_pairs_upper_triangle_and_bounds():
    b = Bins({"a": 1_000, "b": 450}, 200)
    acc, n_out = add_pairs(None, b, np.array([1, 0, 0]), np.array([10, 10, 999]),
                           np.array([0, 0, 1]), np.array([210, 10, 451]))
    M = symmetrize(acc).toarray()
    assert n_out == 1                     # b:451 is beyond chrom b (450 bp): dropped, not binned into another chrom
    assert M[1, 5] == 1 and M[5, 1] == 1  # b:10 (global bin 5) x a:210 (bin 1)
    assert M[0, 0] == 1                   # cis pair inside bin 0 counted once on the diagonal


def test_coarsen_sums_counts():
    fine = Bins({"a": 800}, 200)
    coarse = Bins({"a": 800}, 400)
    M = sp.csr_matrix(np.arange(16, dtype=float).reshape(4, 4) + np.arange(16, dtype=float).reshape(4, 4).T)
    C = coarsen(M, fine, coarse).toarray()
    assert C.shape == (2, 2) and C[0, 0] == sp.triu(M[:2, :2]).sum() and C[0, 1] == M[:2, 2:].sum()
    assert C[1, 0] == C[0, 1]             # coarse diagonal counts each within-bin pair once (upper triangle)


def test_coarsen_equals_direct_binning():
    rng = np.random.default_rng(0)
    chroms = {"a": 1_000, "b": 750}
    L = np.array([1_000, 750])
    n = 500
    ci1, ci2 = rng.integers(0, 2, n), rng.integers(0, 2, n)
    p1 = (rng.random(n) * L[ci1]).astype(np.int64)
    p2 = (rng.random(n) * L[ci2]).astype(np.int64)
    fine, coarse = Bins(chroms, 200), Bins(chroms, 400)
    via_fine = coarsen(symmetrize(add_pairs(None, fine, ci1, p1, ci2, p2)[0]), fine, coarse)
    direct = symmetrize(add_pairs(None, coarse, ci1, p1, ci2, p2)[0])
    assert np.array_equal(via_fine.toarray(), direct.toarray())
    assert via_fine.diagonal().sum() > 0  # pairs within one coarse bin are present on the diagonal


def test_matrix_from_pairs_crosses_chunks(tmp_path):
    # simple 6-column pairs: read, chrom1, pos1, strand1, chrom2, pos2, strand2 (0-based positions)
    path = tmp_path / "pairs.txt"
    path.write_text(
        "r1\ta\t10\t+\ta\t210\t+\n"      # bins (0, 1)
        "r2\ta\t10\t+\ta\t10\t+\n"       # diagonal bin 0
        "r3\ta\t5000\t+\ta\t10\t+\n"     # a:5000 beyond 1000 bp: dropped
        "r4\tb\t10\t+\ta\t999\t+\n"      # global bins (5, 4)
        "r5\tb\t451\t+\ta\t10\t+\n"      # b:451 beyond 450 bp: dropped
    )
    chroms = {"a": 1_000, "b": 450}
    stats = PairStats()
    M, n_out = matrix_from_pairs(path, chroms, Bins(chroms, 200), stats, chunksize=2)
    dense = M.toarray()
    assert stats.read == 5 and n_out == 2
    assert np.array_equal(dense, dense.T)
    assert dense[0, 0] == 1 and dense[0, 1] == 1 and dense[1, 0] == 1
    assert dense[4, 5] == 1 and dense[5, 4] == 1
    assert dense.sum() == 5 and M.dtype == np.float64   # (0,0) once; (0,1),(1,0); (4,5),(5,4)


def test_matrix_from_pairs_no_match_is_zero_matrix(tmp_path):
    path = tmp_path / "pairs.txt"
    path.write_text("r1\tz\t10\t+\tz\t20\t+\n")   # chromosome not in the genome
    chroms = {"a": 1_000, "b": 450}
    M, n_out = matrix_from_pairs(path, chroms, Bins(chroms, 200), PairStats(), chunksize=2)
    assert M.shape == (8, 8) and M.nnz == 0 and n_out == 0


def test_expected_contacts_formula():
    b = Bins({"a": 20 * 100, "b": 20 * 100, "c": 20 * 100}, 100)
    M = sp.lil_matrix((b.N, b.N))
    for i in range(b.N):                  # every bin: 30 trans contacts spread on the other chromosomes
        for j in range(b.N):
            if b.chrom_of[i] != b.chrom_of[j]:
                M[i, j] = 30 / 40
    e = expected_cen_contacts(sp.csr_matrix(M), b, np.ones(3, bool))
    win = 2 * params.HALF_WIN + 1
    assert np.isclose(e, 30 * (2 * win * 100) / (2 * 2000))   # trans x (2 partners x window bp) / other length


def test_expected_contacts_excludes_non_members():
    # chromosome d is not a member: its 60 contacts with the members must not count, and it is not a centre
    b = Bins({"a": 2000, "b": 2000, "c": 2000, "d": 2000}, 100)
    M = sp.lil_matrix((b.N, b.N))
    for i in range(b.N):
        for j in range(b.N):
            if b.chrom_of[i] != b.chrom_of[j]:
                M[i, j] = 30 / 40         # per bin: 30 contacts with a,b,c and 15 with d
    members = np.array([True, True, True, False])
    e = expected_cen_contacts(sp.csr_matrix(M), b, members)
    win = 2 * params.HALF_WIN + 1
    assert np.isclose(e, 30 * (2 * win * 100) / (2 * 2000))


def _uniform_get(bins_per_chrom, contact):
    """get(size) for three chromosomes of bins_per_chrom bins each, with `contact` between every cross-chromosome bin pair."""
    def get(size):
        b = Bins({"a": bins_per_chrom * size, "b": bins_per_chrom * size, "c": bins_per_chrom * size}, size)
        M = sp.csr_matrix(np.where(b.chrom_of[:, None] != b.chrom_of[None, :], contact, 0.0))
        return M, b
    return get


def test_choose_bin_size_shallow_library_falls_back_to_largest():
    def get(size):
        b = Bins({"a": 40 * size, "b": 40 * size, "c": 40 * size}, size)
        M = sp.csr_matrix((b.N, b.N))                      # no contacts at all
        return M, b
    size, table = choose_bin_size(get, (100, 200))
    assert size == 200 and not table.meets_rule.any() and table.chosen.sum() == 1


def test_choose_bin_size_smallest_meeting_rule_is_chosen():
    # 40 bins per chromosome at 100 bp and 20 at 200 bp, both with expected centromere contacts of 36 >= 25
    size, table = choose_bin_size(_uniform_get(40, 2.0), (100, 200))
    assert size == 100 and table.meets_rule.all()
    assert table.chosen.tolist() == [True, False]


def test_choose_bin_size_no_callable_chromosomes_gives_none():
    # 10 bins per chromosome is below MIN_BINS, so nothing is callable at any size
    size, table = choose_bin_size(_uniform_get(10, 2.0), (100, 200))
    assert size is None and (table.n_callable == 0).all()
    assert not table.chosen.any() and not table.meets_rule.any()
