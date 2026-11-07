"""Insertion age as the second line that lets validated Hi-C move a high call (ageline; docs/future.md #31)."""
from types import SimpleNamespace

import numpy as np
import pandas as pd

from cenfuse import ageline, fuse

L = 50_000_000
BIN = 250_000


def _gff(path, elements):
    """elements: (chrom, start, end, clade, family, K2P_T or None)"""
    with open(path, "w") as f:
        f.write("##gff-version 3\n")
        for i, (c, s, e, cl, fa, t) in enumerate(elements):
            age = f";K2P_T={t}" if t is not None else ""
            f.write(f"{c}\tLTRquest\tLTR_retrotransposon\t{s + 1}\t{e}\t.\t+\t.\tID=te{i};clade={cl};family={fa}{age}\n")
    return path


def _sig(cores=None, level="family", seed="famA", lineage="CRM"):
    return SimpleNamespace(cores=cores or {}, level=level, seed=seed, lineage=lineage)


def _elements(young_at=20_000_000, old_at=5_000_000, n=8, fam="famA", chrom="c1"):
    rng = np.random.default_rng(3)
    els = [(chrom, young_at + 50_000 * i, young_at + 50_000 * i + 8_000, "CRM", fam, int(rng.integers(0, 20_000))) for i in range(n)]
    els += [(chrom, old_at + 50_000 * i, old_at + 50_000 * i + 8_000, "CRM", fam, int(rng.integers(300_000, 600_000))) for i in range(n)]
    els += [(chrom, 35_000_000 + 50_000 * i, 35_000_000 + 50_000 * i + 8_000, "Tekay", "famB", 1) for i in range(n)]
    return els


def test_younger_copies_at_the_peak_give_a_small_p(tmp_path):
    t = ageline.load(_gff(tmp_path / "a.gff3", _elements()), {"c1": L})
    m = ageline.family_mask(t, _sig({"c1": dict(families=["famA"])}), "c1")
    r = ageline.compare(t, m, (20_000_000, 20_400_000), (5_000_000, 5_400_000), 4 * BIN)
    assert r["n_peak"] == 8 and r["n_call"] == 8 and r["med_peak"] < 20_000 < 300_000 < r["med_call"] and r["p"] < 0.001
    r = ageline.compare(t, m, (5_000_000, 5_400_000), (20_000_000, 20_400_000), 4 * BIN)       # older at the peak: no
    assert r["p"] > 0.5


def test_too_few_dated_copies_give_no_test(tmp_path):
    t = ageline.load(_gff(tmp_path / "a.gff3", _elements(n=2)), {"c1": L})
    m = ageline.family_mask(t, _sig({"c1": dict(families=["famA"])}), "c1")
    assert np.isnan(ageline.compare(t, m, (20_000_000, 20_100_000), (5_000_000, 5_100_000), 4 * BIN)["p"])


def test_only_the_chromosomes_signature_families_count(tmp_path):
    t = ageline.load(_gff(tmp_path / "a.gff3", _elements()), {"c1": L})
    assert ageline.family_mask(t, _sig({"c1": dict(families=["famB"])}), "c1").sum() == 8
    assert ageline.family_mask(t, _sig({}, level="lineage", seed="CRM", lineage="CRM"), "c1").sum() == 16   # lineage seed
    assert ageline.family_mask(t, _sig({}), "c1").sum() == 16                                               # the family seed
    assert ageline.family_mask(t, _sig({}), "c9").sum() == 0


def test_no_ages_no_table(tmp_path):
    els = [(c, s, e, cl, fa, None) for c, s, e, cl, fa, _ in _elements()]
    assert ageline.load(_gff(tmp_path / "a.gff3", els), {"c1": L}) is None


def test_for_genome_compares_the_anchored_core_with_the_call(tmp_path):
    gff = _gff(tmp_path / "a.gff3", _elements())
    row = dict(anch_core_start=20_000_000, anch_core_end=20_400_000, bin_size=BIN, gate_pass=True, anch_strong=True)
    h = pd.DataFrame([dict(row, chrom="c1"), dict(row, chrom="c2", anch_core_start=-1, anch_core_end=-1),
                      dict(row, chrom="c3", gate_pass=False), dict(row, chrom="c4", anch_strong=False), dict(row, chrom="c5"),
                      dict(row, chrom="c6")])
    calls = pd.DataFrame([dict(chrom="c1", start=5_000_000, end=5_400_000, confidence="high"),
                          dict(chrom="c2", start=1_000, end=2_000, confidence="high")] +
                         [dict(chrom=c, start=5_000_000, end=5_400_000, confidence="high") for c in ("c3", "c4")] +
                         [dict(chrom="c5", start=5_000_000, end=5_400_000, confidence="medium"),       # moves without a second line
                          dict(chrom="c6", start=20_100_000, end=20_300_000, confidence="high")])      # the peak agrees: no tie
    chroms = {f"c{i}": L for i in range(1, 7)}
    r = ageline.for_genome(gff, chroms, _sig({"c1": dict(families=["famA"])}), h, calls)
    assert set(r) == {"c1"} and r["c1"]["p"] < 0.001
    assert ageline.for_genome(gff, chroms, None, h, calls) == {}


def test_the_premise_needs_young_copies_in_the_other_chromosomes_cores(tmp_path):
    """Young copies concentrated in the cores of the other chromosomes (CRM, Bryco): P small; old ones there (Athila-like
    accumulation near centromeres): no premise."""
    rng = np.random.default_rng(1)
    els = []
    for c in ("c1", "c2", "c3", "c4"):
        els += [(c, 10_000_000 + 30_000 * i, 10_000_000 + 30_000 * i + 8_000, "CRM", "famA", int(rng.integers(0, 50_000))) for i in range(12)]
        els += [(c, 30_000_000 + 30_000 * i, 30_000_000 + 30_000 * i + 8_000, "CRM", "famA", int(rng.integers(200_000, 900_000))) for i in range(12)]
    t = ageline.load(_gff(tmp_path / "a.gff3", els), {c: L for c in ("c1", "c2", "c3", "c4")})
    cores = {c: dict(start=9_900_000, end=10_500_000, significant=True, families=["famA"]) for c in ("c1", "c2", "c3", "c4")}
    assert ageline.premise(t, _sig(cores), "c1") < 1e-6
    old = {c: dict(start=29_900_000, end=30_500_000, significant=True, families=["famA"]) for c in ("c1", "c2", "c3", "c4")}
    assert ageline.premise(t, _sig(old), "c1") > 0.5
    lone = {"c1": dict(start=9_900_000, end=10_500_000, significant=True, families=["famA"])}
    assert np.isnan(ageline.premise(t, _sig(lone), "c1"))                    # nothing outside the chromosome itself


def _ev(age):
    s = dict(start=5_000_000, end=5_400_000, evidence="satellite", confidence="high", flags=".", main_family="satA")
    h = dict(bin_size=BIN, gate_pass=True, anch_strong=True, anch_level="high", anch_peak_start=20_000_000, anch_peak_end=20_250_000,
             anch_core_start=20_000_000, anch_core_end=20_400_000, masked=())
    return fuse.decide(fuse.Evidence("c1", L, s, None, h, None, lambda a, b: 300_000, age=age))


def test_younger_copies_at_a_strong_validated_peak_move_a_high_call():
    d = _ev(dict(n_peak=8, n_call=8, med_peak=5_000.0, med_call=400_000.0, p=1e-4, premise=1e-6))
    assert d["call"] == (20_000_000, 20_400_000) and d["relation"] == "offset" and d["confidence"] == "high"
    assert "moved_with_age" in d["flags"] and "age" in d["support"] and d["second"] == ["age"]
    row = fuse.row(fuse.Evidence("c1", L, dict(start=5_000_000, end=5_400_000, evidence="satellite", confidence="high", flags=".",
                                               main_family="satA"), None, None, age=dict(n_peak=8, n_call=8, med_peak=5_000.0,
                                                                                         med_call=400_000.0, p=1e-4, premise=1e-6)), d)
    assert row["age_peak"] == "5000 (8)" and row["age_call"] == "400000 (8)" and row["age_p"] == 1e-4


def test_without_a_significant_age_difference_the_high_call_is_contested():
    for age in (None, dict(n_peak=8, n_call=8, med_peak=5_000.0, med_call=6_000.0, p=0.4, premise=1e-6), dict(n_peak=2, n_call=8, p=np.nan),
                dict(n_peak=8, n_call=8, med_peak=5_000.0, med_call=400_000.0, p=1e-4, premise=0.3),        # premise not shown
                dict(n_peak=8, n_call=8, med_peak=5_000.0, med_call=400_000.0, p=1e-4)):
        d = _ev(age)
        assert d["call"] == (5_000_000, 5_400_000) and d["relation"] == "contested" and "hic_offset" in d["flags"]
