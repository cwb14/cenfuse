import numpy as np
import pytest

from cenfuse.meth import params, track


def test_d_formula_reference_points():
    assert track.d_formula(2.18e9, 10) == pytest.approx(2.22e6)                  # maize Mo17 anchor
    assert track.d_formula(1.806e9, 14) == pytest.approx(1.32e6, rel=0.01)      # P. annua


def test_geometry_matches_validated_poa_settings():
    g = track.geometry(1.29e6)
    assert (g.w, g.lo, g.hi, g.sep) == (10, 39, 103, 103)


def test_chromosome_set_excludes_short_scaffolds():
    fai = {**{f"chr{i}": 60_000_000 for i in range(1, 7)}, "scaf1": 500_000}
    names, G, D = track.chromosome_set(fai)
    assert names == [f"chr{i}" for i in range(1, 7)]
    assert D == pytest.approx(track.d_formula(G, 6))
    assert track.chromosome_count([60e6] * 6 + [5e5, 4e5]) == 6


def test_chromosome_set_homogeneous_genome_keeps_all_chromosomes():
    fai = {f"c{i}": 129_000_000 for i in range(14)}                         # 14 x 129 Mb, no scaffolds
    names, G, D = track.chromosome_set(fai)
    assert names == [f"c{i}" for i in range(14)]
    assert D == pytest.approx(track.d_formula(14 * 129e6, 14))
    assert D == pytest.approx(1.32e6, rel=0.01)


def _poa_like():
    """14 chromosomes of 73-320 Mb and 355 scaffolds of 0.2-0.74 Mb."""
    chroms = [320e6, 98e6, 266e6, 146e6, 134e6, 86e6, 101e6, 73e6, 106e6, 113e6, 102e6, 73e6, 89e6, 74e6]
    return chroms, np.linspace(741_000, 200_000, 355).tolist()


def test_chromosome_count_poa_like_ignores_scaffold_tail():
    chroms, scafs = _poa_like()
    lengths = chroms + scafs + [2_000, 10]                                  # tail: a tenfold drop to ignore
    assert track.chromosome_count(lengths) == 14
    fai = {f"chr{i + 1}": int(L) for i, L in enumerate(chroms)}
    fai.update({f"scaf{i}": int(L) for i, L in enumerate(scafs)})
    fai.update({"tail1": 2_000, "tail2": 10})
    names, G, D = track.chromosome_set(fai)
    assert names == [f"chr{i + 1}" for i in range(14)]


def test_chromosome_count_contig_assembly_and_single_sequence():
    contigs = {f"ctg{i}": 1_000_000 for i in range(2000)}                   # 2-Gb genome in 1-Mb contigs
    assert track.chromosome_count(contigs.values()) == 2000
    assert track.chromosome_set(contigs)[2] < 2 * params.BIN                # refused by the CLI in a later task
    assert track.chromosome_count([5e6]) == 1


def test_tenfold_count_is_the_08_rule():
    """cenfuse 0.8's count, kept for --classify-from targets and old output folders: the same four cases."""
    chroms, scafs = _poa_like()
    contigs = {f"ctg{i}": 1_000_000 for i in range(2000)}
    assert track.tenfold_count([60e6] * 6 + [5e5, 4e5]) == 6
    assert track.tenfold_count(chroms + scafs + [2_000, 10]) == 14
    assert track.tenfold_count(contigs.values()) == 2000
    assert track.tenfold_count([5e6]) == 1


def test_gradual_scaffold_tail_is_dropped():
    """24 chromosomes and a 7.7-Mb scaffold, then scaffolds falling gradually to 1 kb (Cissus-like): no tenfold gap
    between neighbours, so 0.8 counted all 669 sequences. The tail starts at the first under a tenth of the median."""
    L = np.concatenate([np.linspace(60e6, 15e6, 24), [7.68e6], np.geomspace(1.6e6, 1e3, 644)])
    assert track.tenfold_count(L) == 669
    assert track.chromosome_count(L) == 25


def test_deep_tail_drop_is_ignored():
    """A tenfold drop far down the tail (1e4 -> 100 bp) made 0.8 keep the whole tail above it; the tail spans 300x."""
    L = np.concatenate([np.linspace(55e6, 25e6, 24), [3.05e6], np.geomspace(5e5, 1e4, 27), [100]])
    assert track.tenfold_count(L) == 52
    assert track.chromosome_count(L) == 24


def test_dot_chromosomes_are_kept():
    dm6 = [32.08e6, 28.11e6, 25.29e6, 23.54e6, 23.51e6, 3.67e6, 1.35e6, 19.5e3]    # 3R 3L 2R X 2L Y 4 M
    assert track.chromosome_count(dm6) == 7                                         # chr4 kept, as in 0.8
    chicken = np.concatenate([np.geomspace(197e6, 2.5e6, 40), [1.5e6, 1.0e6, 0.6e6, 0.4e6, 0.3e6], [16.8e3]])
    assert track.chromosome_count(chicken) == 45                                    # five dot chromosomes, 0.3-1.5 Mb


def test_documented_limit_dot_chromosome_alone_after_a_tenfold_gap():
    """Without the Y and M, chr4 is over tenfold below every other chromosome with nothing between: dropped, as in
    0.8."""
    assert track.chromosome_count([32.08e6, 28.11e6, 25.29e6, 23.54e6, 23.51e6, 1.35e6]) == 5


def test_the_scan_starts_at_the_half_index():
    """A sequence at 44% of the assembly, then a gradual tail: the sequence that carries the total past half is already
    the tail's first (a tenth below the first sequence), so the set is one sequence long."""
    L = np.concatenate([[44.0], np.geomspace(4.0, 0.1, 41)])
    assert track.tenfold_count(L) == 42
    assert track.chromosome_count(L) == 1


def _by_definition(lengths):
    """(0.8 count, median-cut count, 0.9 count) by their plain definitions: one explicit loop each, with np.median
    recomputed for every k (lengths sorted descending)."""
    L = np.sort(np.asarray(lengths, float))[::-1]
    n, half = len(L), int(np.searchsorted(np.cumsum(L), L.sum() / 2))
    n0 = next((i + 1 for i in range(half, n - 1) if L[i] >= 10 * L[i + 1]), n)
    n2 = next((k for k in range(max(half, 1), n) if 10 * L[k] < np.median(L[:k])), n)
    if n2 >= n0:
        return n0, n2, n0
    small = L[n2:n0]
    return n0, n2, n0 if small.max() < 10 * small.min() else n2


def _random_lengths(rng, i):
    """Length vectors of varied shape: karyotypes (chromosomes, dot chromosomes, a scaffold tail), and vectors of tied
    values with exact tenfold ratios (10 L[k] == median: the strict < matters). The first four have 1-4 sequences."""
    if i < 4:
        return rng.uniform(1, 100, i + 1)
    if i % 3 == 0:
        n = int(rng.integers(5, 60))
        return rng.integers(1, 4, n) * 10.0 ** rng.integers(0, 3, n)
    chroms = rng.uniform(20e6, 100e6, rng.integers(2, 25))
    med = np.median(chroms)
    dots = med * rng.uniform(0.04, 0.15, rng.integers(0, 6)) if i % 3 == 1 else []
    tail = med * np.exp(rng.uniform(np.log(1e-4), np.log(0.09), rng.integers(0, 80)))
    return rng.permutation(np.concatenate([chroms, dots, tail]))


def test_counts_equal_their_definitions_on_random_lengths():
    """The vectorised median cut (median(L[:k]) read off two positions of the sorted lengths) equals the loop that
    recomputes np.median for every k, and chromosome_count equals the rule written out in full: on karyotype-like
    vectors, on ties, and on 1-4 sequences. Every branch of the rule must occur, or the test says nothing."""
    rng = np.random.default_rng(11)
    branches = set()
    for i in range(600):
        L = _random_lengths(rng, i)
        n0, n2, n9 = _by_definition(L)
        assert track.tenfold_count(L) == n0, (i, sorted(L, reverse=True))
        assert track._median_cut(np.sort(L)[::-1]) == n2, (i, sorted(L, reverse=True))
        assert track.chromosome_count(L) == n9, (i, sorted(L, reverse=True))
        branches.add("no tail" if n2 >= n0 else "kept: small members within tenfold" if n9 == n0 else "cut at the tail")
    assert branches == {"no tail", "kept: small members within tenfold", "cut at the tail"}


def test_a_million_contigs_are_counted_in_one_pass():
    """No sequence is under a tenth of the median, so the scan visits every one of 10^6: a median per sequence
    would not finish."""
    L = np.random.default_rng(5).lognormal(np.log(3_000), 0.3, 1_000_000)
    assert track.chromosome_count(L) == 1_000_000


def test_chromosome_set_n_overrides_the_count_for_D():
    """cenfuse passes the size of its own chromosome set (a --classify-from target's 0.8 set, a listed set)."""
    fai = {**{f"chr{i}": 60_000_000 for i in range(1, 7)}, "mid": 9_000_000}
    names, G, D = track.chromosome_set(fai, n=3)
    assert D == track.d_formula(G, 3) and "mid" not in names                # 10 D = 13 Mb
    names, G, D = track.chromosome_set(fai)
    assert D == track.d_formula(G, 7) and "mid" in names                    # 7 counted: 10 D = 5.7 Mb


def test_ring_and_slide_match_brute_force():
    rng = np.random.default_rng(0)
    x = rng.random((3, 40))
    lo, hi, w = 3, 7, 4
    ring = track.ring_sum(x, lo, hi)
    sl = track.slide(x, w)
    for i in range(40):
        js = [j for j in range(40) if lo <= abs(i - j) <= hi]
        assert np.allclose(ring[:, i], x[:, js].sum(axis=1))
        a, b = track.window(i, w, 40)
        assert np.allclose(sl[:, i], x[:, a:b].sum(axis=1))


def _two_classes(nb=400, centre=slice(199, 202)):
    v = np.full((2, nb), 1000.0)
    m = np.vstack([np.full(nb, 700.0), np.full(nb, 450.0)])   # class 0 'satellite' 0.70, class 1 'CRM' 0.45
    return v, m, centre


def test_adjusted_ignores_composition_but_sees_class_loss():
    geo = track.geometry(2.5 * params.BIN)          # w = 1 bin, ring 4-10 bins (beyond the 3-bin centre)
    assert (geo.w, geo.lo, geo.hi) == (1, 4, 10)
    v, m, c = _two_classes()
    v[0, c] *= 5                                    # composition change only: more satellite, same levels
    m[0, c] *= 5
    num, vok = track.class_sums(v, m, geo)
    y = track.adjusted(num, vok)
    plain = m.sum(0) / v.sum(0)
    assert plain[200] - plain[100] > 0.05          # plain methylation rises with the satellite share
    assert abs(y[200]) < 1e-9                      # the adjusted track does not move
    v, m, c = _two_classes()
    m[1, c] = 200.0                                 # CRM loses 0.25 in the centre
    num, vok = track.class_sums(v, m, geo)
    y = track.adjusted(num, vok)
    assert y[200] == pytest.approx(-0.125)          # half of the calls lose 0.25
    assert abs(y[100]) < 1e-9


def test_ratio_needs_min_calls():
    r = track.ratio(np.array([10.0, 10.0]), np.array([params.MIN_CALLS - 1, params.MIN_CALLS]))
    assert np.isnan(r[0]) and r[1] == pytest.approx(10 / params.MIN_CALLS)
