"""Tests for confidence: bootstrap, leave one partner out, null, levels."""
from __future__ import annotations

import numpy as np
import pytest

from cenfuse.hic.balance import balanced, ice
from cenfuse.hic.bins import Bins, matrix_from_pairs
from cenfuse.hic.confidence import bootstrap, leave_one_out, level, null_test, poisson_resample
from cenfuse.hic.io import PairStats, read_fai
from cenfuse.hic.locate import locate, peak_core, start_positions, telomere_score
from cenfuse.hic.score import block_means
from hic_synth import make_genome


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    truth = make_genome(tmp_path_factory.mktemp("conf"))
    chroms = read_fai(truth["fai"])
    b = Bins(chroms, truth["size"])
    M, _ = matrix_from_pairs(truth["pairs"], chroms, b, PairStats())
    bias, good, _, _ = ice(M)
    W = balanced(M, bias, good)
    members = b.nb >= 20
    mean = block_means(W, b, good, members)
    loc = locate(W, b, good, mean, members, start_positions(telomere_score(W, b, good, mean, members), b, members))
    pcs = {k: peak_core(loc.S, b, k) for k in loc.est}
    boot = bootstrap(M, b, members, loc.est, n=40, threads=4)
    lopo = leave_one_out(loc.oe, b, good, members, loc.est)
    p = null_test(W, b, good, mean, members, loc.est, {k: pc["peak"] for k, pc in pcs.items()},
                  {k: pc["ratio"] for k, pc in pcs.items()}, n=200)
    return truth, b, pcs, boot, lopo, p


def levels(run):
    truth, b, pcs, boot, lopo, p = run
    out = {}
    for k, pc in pcs.items():
        inside = np.mean((boot[k] >= pc["core_lo"]) & (boot[k] <= pc["core_hi"]))
        stable = bool(((lopo[k] >= pc["core_lo"]) & (lopo[k] <= pc["core_hi"])).all())
        out[truth["names"][k]] = (level(inside, stable, p[k]), p[k], inside, stable)
    return out


def test_planted_centromeres_are_high(run):
    lv = levels(run)
    planted = [nm for nm in run[0]["names"] if nm not in run[0]["unattached"]]
    assert sum(lv[nm][0] == "high" for nm in planted) >= 9
    assert all(lv[nm][0] in ("high", "medium") for nm in planted)


def test_unattached_chromosome_is_not_confident(run):
    lv = levels(run)
    lev, p = lv[run[0]["unattached"][0]][:2]
    assert lev not in ("high", "medium") and p > 0.01


def test_poisson_resample_keeps_symmetry_and_scale():
    import scipy.sparse as sp
    M = sp.csr_matrix(np.array([[10.0, 4.0], [4.0, 6.0]]))
    R = poisson_resample(M, np.random.default_rng(0)).toarray()
    assert np.allclose(R, R.T)


def test_level_rules():
    assert level(0.99, True, 0.005) == "high"
    assert level(0.99, False, 0.005) == "medium"
    assert level(0.5, False, 0.03) == "low"
    assert level(0.99, True, 0.2) == "none"
    assert level(0.99, True, float("nan")) == "medium"     # null untestable: capped at medium


def test_null_untestable_on_short_chromosomes():
    b = Bins({"a": 30 * 100, "b": 30 * 100, "c": 30 * 100}, 100)   # 30 bins: no bin >= 20 from both ends
    import scipy.sparse as sp
    rng = np.random.default_rng(1)
    A = rng.poisson(5, (b.N, b.N)).astype(float)
    A = np.triu(A) + np.triu(A, 1).T
    M = sp.csr_matrix(A)
    bias, good, _, _ = ice(M)
    W = balanced(M, bias, good)
    members = np.ones(3, bool)
    p = null_test(W, b, good, block_means(W, b, good, members), members, {0: 15, 1: 15, 2: 15},
                  {0: 15, 1: 15, 2: 15}, {0: 2.0, 1: 2.0, 2: 2.0}, n=10)
    assert all(np.isnan(v) for v in p.values())


def _random_genome(nbins=100, nchrom=3):
    import scipy.sparse as sp
    b = Bins({c: nbins * 100 for c in "abc"[:nchrom]}, 100)
    rng = np.random.default_rng(2)
    A = rng.poisson(5, (b.N, b.N)).astype(float)
    A = np.triu(A) + np.triu(A, 1).T
    M = sp.csr_matrix(A)
    bias, good, _, _ = ice(M)
    W = balanced(M, bias, good)
    members = np.ones(nchrom, bool)
    return b, W, good, block_means(W, b, good, members), members


def test_null_nan_when_observed_ratio_not_finite():
    b, W, good, mean, members = _random_genome()
    est = {0: 50, 1: 50, 2: 50}
    p = null_test(W, b, good, mean, members, est, est, {0: float("nan"), 1: 2.0, 2: 2.0}, n=10)
    assert np.isnan(p[0]) and np.isfinite(p[1]) and np.isfinite(p[2])


def test_null_nan_when_chromosome_never_scored():
    b, W, good, mean, members = _random_genome()
    members[2] = False            # chromosome 2 gets no score in any replicate
    est = {0: 50, 1: 50, 2: 50}
    p = null_test(W, b, good, mean, members, est, est, {0: 2.0, 1: 2.0, 2: 2.0}, n=10)
    assert np.isnan(p[2]) and np.isfinite(p[0])


def test_threaded_bootstrap_is_warning_free_and_leaves_filters_alone():
    """Empty (masked) bins give all-NaN score columns in every replicate. With 4 threads the bootstrap must neither
    warn nor change the process-wide warnings filters (catch_warnings is not thread-safe)."""
    import warnings

    import scipy.sparse as sp
    b = Bins({c: 40 * 100 for c in "abcd"}, 100)
    rng = np.random.default_rng(4)
    A = rng.poisson(5, (b.N, b.N)).astype(float)
    A = np.triu(A) + np.triu(A, 1).T
    for i in (3, 50, 97, 130):                     # one empty bin per chromosome -> masked
        A[i, :] = 0
        A[:, i] = 0
    M = sp.csr_matrix(A)
    members = np.ones(4, bool)
    est = {0: 20, 1: 20, 2: 20, 3: 20}
    before = list(warnings.filters)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        inside = list(warnings.filters)
        for _ in range(5):
            boot = bootstrap(M, b, members, est, n=40, threads=4)
            assert list(warnings.filters) == inside
    assert list(warnings.filters) == before
    assert all(len(v) == 40 and (v >= 0).all() for v in boot.values())
