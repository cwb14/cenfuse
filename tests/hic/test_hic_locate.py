from __future__ import annotations

import numpy as np
import pytest

from cenfuse.hic import params
from cenfuse.hic.balance import balanced, ice
from cenfuse.hic.bins import Bins, matrix_from_pairs
from cenfuse.hic.io import PairStats, read_fai
from cenfuse.hic.locate import locate, peak_core, rabl_test, start_positions, telomere_score
from cenfuse.hic.score import block_means
from hic_synth import make_genome


def prepare(truth):
    chroms = read_fai(truth["fai"])
    b = Bins(chroms, truth["size"])
    M, _ = matrix_from_pairs(truth["pairs"], chroms, b, PairStats())
    bias, good, _, _ = ice(M)
    W = balanced(M, bias, good)
    members = b.nb >= 20
    return b, M, W, good, members, block_means(W, b, good, members)


@pytest.fixture(scope="module")
def rabl(tmp_path_factory):
    truth = make_genome(tmp_path_factory.mktemp("rabl"))
    return truth, prepare(truth)


def test_rabl_polarity_detected(rabl):
    truth, (b, M, W, good, members, mean) = rabl
    T = telomere_score(W, b, good, mean, members)
    r = rabl_test(T, b, members)
    assert r["n"] == 11 and r["passed"] and r["negative"] >= 10


@pytest.mark.parametrize("seed", [1, 2, 3])      # seed 1 ends in a 4-state cycle; seeds 2 and 3 at a fixed point
def test_prior_free_finds_planted_centromeres(tmp_path, seed):
    truth = make_genome(tmp_path, seed=seed)
    b, M, W, good, members, mean = prepare(truth)
    T = telomere_score(W, b, good, mean, members)
    loc = locate(W, b, good, mean, members, start_positions(T, b, members))
    assert loc.converged
    for k, nm in enumerate(truth["names"]):
        if nm in truth["unattached"]:
            continue
        pc = peak_core(loc.S, b, k)
        assert abs(pc["peak"] - truth["cens"][k]) <= 1, (nm, pc["peak"] - truth["cens"][k])
        assert pc["core_lo"] <= truth["cens"][k] <= pc["core_hi"], nm
        assert pc["ratio"] > 2.5, nm


def test_peak_core_stops_at_masked_bins():
    b = Bins({"a": 1000}, 100)
    S = np.array([1, 1, 2, np.nan, 5, 4, 1, 1, 1, 1.0])
    pc = peak_core(S, b, 0)
    assert pc["peak"] == 4 and pc["core_lo"] == 4 and pc["core_hi"] == 5


def test_no_rabl_genome_fails_the_test(tmp_path):
    truth = make_genome(tmp_path, rabl=False, hub=None, n_trans=150_000, n_cis=150_000)
    b, M, W, good, members, mean = prepare(truth)
    r = rabl_test(telomere_score(W, b, good, mean, members), b, members)
    assert not r["passed"]


def test_few_chromosomes_report_p_floor(tmp_path):
    truth = make_genome(tmp_path, lengths=(100, 120, 90, 110, 96), cens=(40, 30, 60, 55, 20), unattached=(),
                        hub=None, sat_shift=None, rdna_on=None, n_trans=150_000, n_cis=150_000)
    b, M, W, good, members, mean = prepare(truth)
    r = rabl_test(telomere_score(W, b, good, mean, members), b, members)
    assert r["n"] == 5 and r["p_floor"] == pytest.approx(0.5 ** 5) and r["p_floor"] > 0.01 and not r["passed"]


@pytest.fixture(scope="module")
def seed1(tmp_path_factory):
    """Seed 1 falls into a 4-state cycle: planted centromeres alternate between neighbouring bins, while the
    unattached chr11 (no cluster) jumps between distant bins (13 and 97)."""
    truth = make_genome(tmp_path_factory.mktemp("seed1"), seed=1)
    b, M, W, good, members, mean = prepare(truth)
    T = telomere_score(W, b, good, mean, members)
    return truth, b, locate(W, b, good, mean, members, start_positions(T, b, members))


def test_cycle_state_estimates_agree_with_peaks(seed1):
    truth, b, loc = seed1
    assert loc.cycle > 0 and loc.converged
    for k, nm in enumerate(truth["names"]):
        if nm in truth["unattached"]:
            continue
        assert abs(loc.est[k] - peak_core(loc.S, b, k)["peak"]) <= 1, nm


def test_cycle_spread_marks_only_the_unattached_chromosome(seed1):
    truth, b, loc = seed1
    assert set(loc.spread) == set(loc.est)
    for k, nm in enumerate(truth["names"]):
        if nm in truth["unattached"]:
            assert loc.spread[k] > params.CYCLE_SPREAD, nm
        else:
            assert loc.spread[k] <= params.CYCLE_SPREAD, (nm, loc.spread[k])


def test_fixed_point_has_zero_spread(tmp_path):
    truth = make_genome(tmp_path, seed=2)                 # seed 2 settles on a fixed point
    b, M, W, good, members, mean = prepare(truth)
    loc = locate(W, b, good, mean, members, start_positions(telomere_score(W, b, good, mean, members), b, members))
    assert loc.cycle == 0 and loc.converged and set(loc.spread) == set(loc.est)
    assert all(v == 0 for v in loc.spread.values())


def test_rabl_test_drops_constant_chromosome():
    b = Bins({"a": 1000, "b": 1000, "c": 1000}, 100)
    x = np.arange(10, dtype=float)
    T = np.concatenate([np.abs(x - 4.5), np.abs(x - 4.5), np.ones(10)])     # a, b: V-shaped (rho defined); c: constant
    r = rabl_test(T, b, np.array([True, True, True]))
    assert r["n"] == 2 and list(r["per_chrom"].chrom) == ["a", "b"] and r["p_floor"] == pytest.approx(0.25)
