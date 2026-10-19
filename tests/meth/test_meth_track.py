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


def test_chromosome_count_poa_like_ignores_scaffold_tail():
    chroms = [320e6, 98e6, 266e6, 146e6, 134e6, 86e6, 101e6, 73e6, 106e6, 113e6, 102e6, 73e6, 89e6, 74e6]
    scafs = np.linspace(741_000, 200_000, 355).tolist()
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
