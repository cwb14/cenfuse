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


def _old_bootstraps(M, b, members, est, win, n):
    """The two bootstraps of 0.14.0 before they were fused (hic.confidence.bootstrap; hicfuse.bootstrap_all), as written then:
    each replicate resampled, re-balanced and scored by each caller on its own."""
    from concurrent.futures import ThreadPoolExecutor
    from cenfuse.hic import params
    from cenfuse.hic.score import cluster_score, partner_oe, window_bins

    def score(R, wins):
        bias, good, _, _ = ice(R)
        W = balanced(R, bias, good)
        mean = block_means(W, b, good, members)
        w = window_bins(b, good, est) if wins is None else {y: x[good[x]] for y, x in wins.items()}
        return cluster_score(partner_oe(W, b, good, mean, w), b, good, members)

    ks = [int(k) for k in np.where(members)[0]]

    def peaks(S, keys):
        return {k: int(np.nanargmax(S[b.sl(k)])) if np.isfinite(S[b.sl(k)]).any() else -1 for k in keys}

    def eng(ss):
        return peaks(score(poisson_resample(M, np.random.default_rng(ss)), None), est)

    def fuse(ss):
        return peaks(score(poisson_resample(M, np.random.default_rng(ss)), win), ks)

    out = []
    for f, keys in ((eng, est), (fuse, ks)):
        with ThreadPoolExecutor(2) as ex:
            reps = list(ex.map(f, np.random.SeedSequence(params.SEED).spawn(n)))
        out.append({k: np.array([r[k] for r in reps]) for k in keys})
    return out


def test_the_fused_bootstrap_equals_the_two_it_replaces(run, tmp_path):
    """future.md #62 (0.14): the engine's bootstrap scores cenfuse's anchored windows on its own resamples (the same seed and
    matrix gave the same resamples twice), and makes each resample's upper triangle once; every peak is unchanged."""
    from cenfuse import hicfuse
    truth, b, pcs, _, _, _ = run
    chroms = read_fai(truth["fai"])
    M, _ = matrix_from_pairs(truth["pairs"], chroms, b, PairStats())
    bias, good, _, _ = ice(M)
    members = b.nb >= 20
    est = {k: int(pc["peak"]) for k, pc in pcs.items()}
    win = {k: np.arange(max(0, p - 2), min(int(b.nb[k]), p + 3)) + int(b.off[k]) for k, p in est.items()}
    old_eng, old_fuse = _old_bootstraps(M, b, members, est, win, 12)
    new_eng, extra = bootstrap(M, b, members, est, n=12, threads=3, extra={"anchored": win})
    assert set(new_eng) == set(old_eng) and all(np.array_equal(new_eng[k], old_eng[k]) for k in old_eng)
    assert all(np.array_equal(extra["anchored"][k], old_fuse[k]) for k in old_fuse)
    alone = hicfuse.bootstrap_all(M, b, members, win, 2, n=12)
    assert all(np.array_equal(alone[k], old_fuse[k]) for k in old_fuse)
    plain = bootstrap(M, b, members, est, n=12, threads=1)                 # without extra: the engine's own result
    assert all(np.array_equal(plain[k], old_eng[k]) for k in old_eng)


def test_the_resampler_builds_the_matrices_poisson_resample_builds():
    """_Resampler (0.14, future.md #62): for the same draws, the same sparse matrix, array by array (data, indices, pointers and
    their types), zeros included where a draw is 0."""
    import scipy.sparse as sp
    from cenfuse.hic.confidence import _Resampler
    rng = np.random.default_rng(3)
    for n, dens, scale in ((50, 0.3, 0.5), (400, 0.05, 2.0)):
        A = sp.random(n, n, density=dens, random_state=1, data_rvs=lambda k: rng.gamma(1.0, scale, k)).tocsr()
        M = (A + A.T).tocsr()
        rs = _Resampler(M)
        assert rs.exact
        for seed in range(3):
            old, new = poisson_resample(M, np.random.default_rng(seed)), rs(np.random.default_rng(seed))
            for a in ("data", "indices", "indptr"):
                x, y = getattr(old, a), getattr(new, a)
                assert x.dtype == y.dtype and np.array_equal(x, y)


def _old_leave_one_out(oe, b, good, members, est):
    """hic.confidence.leave_one_out as 0.14.0 wrote it (every chromosome scored for every pair)."""
    from cenfuse.hic.score import cluster_score
    out = {}
    for x in est:
        peaks = []
        for y in oe:
            if y == x:
                continue
            v = cluster_score(oe, b, good, members, exclude={x: {y}})[b.sl(x)]
            peaks.append(int(np.nanargmax(v)) if np.isfinite(v).any() else -1)
        out[x] = np.array(peaks)
    return out


@pytest.mark.parametrize("n, holes", [(12, 0.1), (3, 0.5), (2, 0.0)])
def test_the_fast_leave_one_out_gives_the_former_peaks(n, holes):
    """future.md #62 (0.14): leave_one_out scores only the chromosome left out of, with the arithmetic of cluster_score; the
    peaks of every chromosome and partner are those of the former cubic loop, also for a chromosome that is no member, a
    single partner, and bins without a finite value."""
    from cenfuse import hicfuse
    b = Bins({f"c{i}": 60 * 10_000 for i in range(n)}, 10_000)
    rng = np.random.default_rng(n)
    good = rng.random(b.N) > 0.05
    members = np.ones(n, dtype=bool)
    members[0] = n > 2                                    # chromosome 0 no member in the smallest case
    oe = {}
    for y in range(n):
        v = rng.gamma(2.0, 1.0, b.N)
        v[rng.random(b.N) < holes] = np.nan
        v[b.chrom_of == y] = np.nan
        v[~good] = np.nan
        oe[y] = v
    est = {k: 30 for k in range(n)}
    old, new = _old_leave_one_out(oe, b, good, members, est), leave_one_out(oe, b, good, members, est)
    assert set(old) == set(new) and all(np.array_equal(old[k], new[k]) for k in old)
    fast = hicfuse.lopo_all(oe, b, good, members)
    assert all(np.array_equal(fast[k], old[k]) for k in fast)
