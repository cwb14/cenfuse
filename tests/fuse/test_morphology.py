"""Chromosome morphology and assembly completeness (cenfuse.morphology): Levan classes, the arms measured from the
kinetochore or the call, telomere ends, assembly gaps and the summary fields."""
import pandas as pd
import pysam
import pytest

from cenfuse import morphology as mo
from cenfuse.sat.params import Params


def test_gap_rule_is_the_satellite_engines():
    assert mo.GAP_MIN_N == Params().min_gap_n


@pytest.mark.parametrize("r,cls", [(1.0, "metacentric"), (1.7, "metacentric"), (1.71, "submetacentric"),
                                   (3.0, "submetacentric"), (3.01, "subtelocentric"), (7.0, "subtelocentric"),
                                   (7.01, "acrocentric"), (float("inf"), "acrocentric")])
def test_levan_classes(r, cls):
    assert mo.levan(r) == cls


def test_classes_over_an_interval():
    L = 100
    assert mo.classes(49, 51, L) == ["metacentric"]
    assert mo.classes(2, 4, L) == ["acrocentric"]
    # short arm 20 -> 30: arm ratios 4.0 (subtelocentric) to 2.33 (submetacentric)
    assert mo.classes(20, 30, L) == ["submetacentric", "subtelocentric"]
    # a span over the middle reaches CI 0.5 whatever its edges
    assert mo.classes(40, 60, L) == ["metacentric"]
    assert mo.classes(30, 70, L) == ["metacentric", "submetacentric"]
    assert mo.classes(10, 60, L) == ["metacentric", "submetacentric", "subtelocentric", "acrocentric"]


def _t(rows):
    cols = ["chrom", "chrom_len", "start", "end", "confidence", "cpg_kt_start", "cpg_kt_end"]
    return pd.DataFrame(rows, columns=cols)


def test_columns():
    t = _t([["m", 1000, 451, 550, "high", ".", "."],          # call over the middle
            ["a", 1000, 1, 40, "high", ".", "."],             # call at the start: acrocentric, short arm at the start
            ["e", 1000, 701, 800, "low", 741, 760],           # kinetochore 740-760: point 750, short arm at the end
            ["n", 1000, ".", ".", "none", ".", "."]])        # no call
    chroms = {"m": 1000, "a": 1000, "e": 1000, "n": 1000}
    telo = {"ends": [("m", "start", 0, 50, 50, 0), ("m", "end", 950, 1000, 50, 0), ("a", "end", 980, 1000, 20, 0)]}
    gaps = {"m": [], "a": [(10, 30), (500, 520)], "e": [(600, 610)], "n": []}
    out = mo.columns(t, chroms, telo, gaps).set_index("chrom")
    m, a, e, n = (out.loc[c] for c in "maen")
    assert (m.cen_point, m.short_arm, m.ci, m.arm_ratio, m.morphology) == (500, ".", 0.5, 1.0, "metacentric")
    assert (m.telomeres, m.gaps, m.call_gaps, m.t2t) == ("both", 0, 0, "yes")
    assert a.short_arm == "start" and a.morphology == "acrocentric" and a.ci == 0.02
    assert (a.telomeres, a.gaps, a.call_gaps, a.t2t) == ("end", 2, 1, "no")
    assert e.cen_point == 750 and e.short_arm == "end" and e.ci == 0.25 and e.arm_ratio == 3.0
    assert e.morphology == "submetacentric/subtelocentric"            # the kinetochore spans CI 0.24-0.26
    assert (e.gaps, e.call_gaps) == (1, 0)
    assert (n.cen_point, n.ci, n.morphology, n.call_gaps) == (".", ".", ".", ".")
    assert (n.telomeres, n.gaps, n.t2t) == ("none", 0, "no")


def test_summary():
    t = _t([["m", 1000, 451, 550, "high", ".", "."], ["a", 1000, 1, 40, "high", ".", "."],
            ["e", 1000, 801, 900, "low", ".", "."], ["n", 1000, ".", ".", "none", ".", "."]])
    chroms = {c: 1000 for c in "maen"}
    out = mo.columns(t, chroms, {"ends": []}, {"m": [], "a": [], "e": [(820, 830)], "n": []})
    s = mo.summary(out)
    assert s["karyotype"] == "m:1,t:1,low:1,no_call:1"            # the low call's class is a guess: counted apart
    assert s["t2t_chromosomes"] == "0/4" and s["calls_with_gaps"] == 1 and s["short_arm_at_end"] == 0  # e is low
    out.loc[out.chrom == "e", "confidence"] = "medium"
    assert mo.summary(out)["short_arm_at_end"] == 1


def test_scan_gaps(tmp_path):
    seq = "A" * 95 + "N" * 10 + "C" * 40 + "N" * 5 + "G" * 50 + "n" * 30
    fa = tmp_path / "g.fa"
    fa.write_text(">c1\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    pysam.faidx(str(fa))
    want = [(95, 105), (200, 230)]
    assert mo.scan_gaps(fa, {"c1": len(seq)}) == {"c1": want}
    assert mo.scan_gaps(fa, {"c1": len(seq)}, chunk=100) == {"c1": want}     # a run across a chunk edge is one gap
    assert mo.scan_gaps(fa, {"c1": len(seq)}, chunk=7) == {"c1": want}
