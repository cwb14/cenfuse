"""The one-per-chromosome test (cenfuse.sat.karyotype) and the cap it sets on satellite calls."""
import numpy as np
import pandas as pd

from cenfuse.sat import centromere, karyotype as K
from cenfuse.sat.params import Params

LEN = [50_000_000] * 12
CHROMS = {f"c{i}": L for i, L in enumerate(LEN)}


def arrays(per_chrom):
    """{chrom: [(start, bp)]} to a table of arrays of family F."""
    rows = [(c, s, s + bp, "F") for c, x in per_chrom.items() for s, bp in x]
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "family"])


def test_null_cover_is_below_complete_cover():
    null = K.null_covered(LEN, 12, draws=2000)
    assert 6 <= np.median(null) <= 9 and null.max() <= 12
    assert K.null_covered([10.0, 90.0], 1, draws=10).tolist() == [1] * 10


def test_one_large_cluster_per_chromosome_passes():
    t = K.satellite_test(arrays({c: [(10_000_000, 500_000)] for c in CHROMS}), {"F"}, CHROMS, Params())
    assert t["status"] == "pass" and t["covered"] == 12 and t["p"] < 0.01


def test_several_clusters_on_few_chromosomes_fail():
    x = {c: [(s, 500_000) for s in range(5_000_000, 45_000_000, 5_000_000)] for c in list(CHROMS)[:3]}   # 8 each on 3
    t = K.satellite_test(arrays(x), {"F"}, CHROMS, Params())
    assert t["status"] == "fail" and t["covered"] == 2


def test_a_small_cluster_on_every_chromosome_does_not_save_large_ones_elsewhere():
    x = {c: [(1_000_000, 30_000)] for c in CHROMS}
    for c in list(CHROMS)[:3]:
        x[c] += [(s, 500_000) for s in range(5_000_000, 45_000_000, 5_000_000)]
    assert K.satellite_test(arrays(x), {"F"}, CHROMS, Params())["status"] == "fail"


def test_too_few_chromosomes_are_untested():
    ch = {f"c{i}": 50_000_000 for i in range(3)}
    t = K.satellite_test(arrays({c: [(10_000_000, 500_000)] for c in ch}), {"F"}, ch, Params())
    assert t["status"] == "untested"


def test_no_arrays_of_the_families():
    t = K.satellite_test(arrays({"c0": [(1, 100)]}), {"other"}, CHROMS, Params())
    assert t["status"] == "none" and K.sentence(t) == "."


def test_fewer_clusters_than_chromosomes_are_compared_with_as_many_draws():
    ch = {f"c{i}": 50_000_000 for i in range(20)}
    t = K.satellite_test(arrays({c: [(10_000_000, 500_000)] for c in list(ch)[:12]}), {"F"}, ch, Params())
    assert t["status"] == "pass" and t["features"] == 12 and t["covered"] == 12 and t["p"] < 0.05       # 12 draws, all apart: unlikely
    assert "all 12 clusters" in K.sentence(t)
    t = K.satellite_test(arrays({c: [(10_000_000, 500_000)] for c in list(ch)[:4]}), {"F"}, ch, Params())
    assert t["status"] == "fail" and t["features"] == 4                                                 # 4 apart says nothing, and proves nothing


def test_too_many_sequences_are_contigs_and_cost_nothing():
    n = K.MAX_SEQUENCES + 1
    t = K.covered_test(["a", "b"], [1000] * n)
    assert t["status"] == "untested" and K.sentence(t) == f"not tested on {n} sequences"


def test_sentence_states_the_numbers():
    t = K.satellite_test(arrays({c: [(10_000_000, 500_000)] for c in CHROMS}), {"F"}, CHROMS, Params())
    assert K.sentence(t).startswith("12/12 chromosomes hold the 12 largest clusters") and K.sentence(t).endswith("pass")


def _final(not_one, te_flag=None, evidence="satellite"):
    import pysam
    return dict(call=dict(cen_start=1000, cen_end=2000, evidence=evidence, main_family="F", main_family_bp=1000, sat_bp=500_000,
                          second_cluster_bp=0, cen_families="F"), not_one=not_one, te_flag=te_flag)


class _Fa:
    def fetch(self, c, s, e):
        return "A" * (e - s)


def _conf(**k):
    x = _final(**k)
    return centromere.finalize("c", 1_000_000, x["call"], _Fa(), None, None, Params(subtel_frac_len=0.0001), te_flag=x["te_flag"],
                               not_one=x["not_one"])


def test_a_failed_test_caps_a_satellite_call_an_ltr_core_does_not_back():
    r = _conf(not_one=True)
    assert r["confidence"] == "low" and "not_one_per_chromosome" in r["flags"]
    r = _conf(not_one=True, te_flag="te_concordant")
    assert r["confidence"] == "high" and "not_one_per_chromosome" not in r["flags"]


def test_only_satellite_calls_are_capped():
    r = _conf(not_one=True, evidence="retrotransposon")
    assert r["confidence"] == "medium" and "not_one_per_chromosome" not in r["flags"]
    assert _conf(not_one=False)["confidence"] == "high"


def test_a_thin_core_takes_only_what_nothing_else_holds(tmp_path):
    from cenfuse.sat.annotate import ARR_COLS
    from types import SimpleNamespace
    import pysam
    from cenfuse.sat import longunit
    L = 10_000_000
    fa = tmp_path / "g.fa"
    fa.write_text(">c1\n" + "A" * 1000 + "\n>c2\n" + "A" * 1000 + "\n>c3\n" + "A" * 1000 + "\n")
    pysam.faidx(str(fa))
    arrays = pd.DataFrame([(0, "c1", 5_000_000, 5_200_000, "+", "A", 1100, 0.9, 180)], columns=ARR_COLS)
    g = SimpleNamespace(name="g", chroms={"c1": L, "c2": L, "c3": L}, fasta=fa)
    thin = lambda s: dict(start=s, end=s + 500_000, k=2, n=9, share=0.2, p=0.5, significant=False, mids=np.array([s + 1, s + 2]), fold=9.0)
    strong = dict(start=7_000_000, end=7_500_000, k=30, n=40, share=0.75, p=1e-9, significant=True, mids=np.arange(7_000_000, 7_000_030),
                  strong=True, fold=30.0)
    old = longunit.detect
    longunit.detect = lambda *a, **k: pd.DataFrame(columns=longunit.LU_COLS)
    try:
        run = lambda **k: centromere.call_regions(g, arrays, {"A"}, None, pd.DataFrame(columns=["chrom", "start", "end"]), tmp_path / "t",
                                                  1, Params(), te=("CRM", {"c3": strong}, True), weak_cores={"c1": thin(1_000_000), "c2": thin(2_000_000)}, **k)[0]
        r = run(not_one=True)
        assert r["c1"]["evidence"] == "retrotransposon" and r["c1"]["confidence"] == "low" and "weak_core" in r["c1"]["flags"]
        assert r["c1"]["overruled_start"] == 5_000_000
        assert r["c2"]["evidence"] == "retrotransposon" and r["c2"]["confidence"] == "low"         # no call: the core takes it
        r = run(not_one=False)
        assert r["c1"]["evidence"] == "satellite" and "weak_core" not in r["c1"]["flags"]          # a satellite that passed stands
        assert r["c2"]["evidence"] == "retrotransposon" and r["c2"]["confidence"] == "low"
    finally:
        longunit.detect = old
