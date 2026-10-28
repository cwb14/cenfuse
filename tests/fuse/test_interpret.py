"""basis and caveat sentences (cenfuse.interpret): one test per relation and per flag the fusion rules or the
satellite engine set."""
from types import SimpleNamespace

import pytest

from cenfuse import interpret as it

BASE = dict(chrom="c1", chrom_len=100_000_000, start=40_000_001, end=42_000_000, length=2_000_000,
            call_source="sequence", relation="sequence_only", confidence="high", support="satellite", against=".",
            flags=".", offset_bp=".", seq_start=40_000_001, seq_end=42_000_000, seq_evidence="satellite",
            seq_confidence="high", seq_flags=".", seq_alt_start=".", seq_alt_end=".", main_family="g_sat01_171",
            te_clade=".", te_start=".", te_end=".", hic_gate=".", hic_anch_peak_start=".", hic_anch_peak_end=".",
            hic_denovo_peak_start=".", hic_denovo_peak_end=".", fold_core_start=".", fold_core_end=".",
            cpg_kt_start=".", cpg_kt_end=".", call_gaps=0)


def row(**kw):
    return SimpleNamespace(**{**BASE, **kw})


def ctx(**kw):
    c = dict(unit={"g_sat01_171": 171}, sat_bp=lambda chrom, s, e: 1_500_000, marker=None, decisive=False, hic=True,
             tes=True)
    return {**c, **kw}


def test_satellite_call_with_every_supporting_line():
    r = row(relation="concordant", support="satellite,te,hic,meth,fold", te_clade="CRM", cpg_kt_start=41_000_001,
            cpg_kt_end=41_060_000)
    b = it.basis(r, ctx())
    assert b.startswith("satellite g_sat01_171 (171-bp unit), 1.50 Mb of centromeric satellite in the call")
    for x in ("the CRM core agrees", "validated Hi-C agrees", "CHG methylation is centromere-like",
              "the arm fold agrees", "CpG kinetochore 41.00-41.06 Mb (60 kb)"):
        assert x in b
    assert it.caveat(r, ctx()) == "."


def test_sequence_evidence_kinds():
    assert it.basis(row(seq_evidence="weak_satellite"), ctx()).startswith("weak satellite g_sat01_171")
    assert it.basis(row(seq_evidence="long_unit_tandem", main_family="long_unit"), ctx()).startswith(
        "long-unit tandem array")
    assert "not a satellite family" in it.caveat(row(seq_evidence="long_unit_tandem"), ctx())
    r = row(seq_evidence="retrotransposon", te_clade="CRM", support="te", main_family=".")
    assert it.basis(r, ctx(decisive=True)) == "CRM retrotransposon core (no satellite call)"
    assert it.caveat(r, ctx(decisive=True)) == "."
    r = row(seq_evidence="retrotransposon", te_clade="Bryco", support=".", confidence="low", main_family=".")
    assert "only CRM cores may decide" in it.caveat(r, ctx(decisive=False))


def test_satellite_engine_flags():
    r = row(seq_flags="competing_cluster", confidence="low", seq_alt_start=70_000_001, seq_alt_end=72_000_000)
    c = it.caveat(r, ctx())
    assert "another satellite cluster (70.00-72.00 Mb) competes" in c
    assert c.endswith("-> CENP-A/CENH3 data would decide")
    assert "or Hi-C when mappable DNA separates them" in it.caveat(r, ctx(hic=False))
    assert "or a CRM core (tes/; plants)" in it.caveat(r, ctx(tes=False))
    r = row(seq_flags="competing_cluster,te_discordant", confidence="low", seq_alt_start=70_000_001,
            seq_alt_end=72_000_000, te_clade="Bryco", te_start=70_500_001, te_end=70_700_000)
    c = it.caveat(r, ctx())                             # the marker's core in the rival cluster: said once, there
    assert c.startswith("another satellite cluster (70.00-72.00 Mb) competes, where the Bryco core lies:")
    assert "lies elsewhere" not in c
    c = it.caveat(row(seq_flags="no_dominant_satellite", confidence="low"), ctx(hic=False))
    assert "no satellite family dominates" in c and c.endswith("-> Hi-C would confirm")
    c = it.caveat(row(seq_flags="no_dominant_satellite", confidence="low", hic_gate="fail"), ctx())
    assert c == "no satellite family dominates the genome; Hi-C failed its test here"
    for flag, text in (("low_satellite", "little centromeric satellite"), ("diffuse", "a broad block"),
                       ("near_end", "near a chromosome end"), ("gene_rich", "gene-rich")):
        assert text in it.caveat(row(seq_flags=flag), ctx()), flag
    r = row(seq_flags="te_discordant", te_clade="CRM", te_start=60_000_001, te_end=61_000_000, against="te")
    c = it.caveat(r, ctx())
    assert c.count("CRM core lies elsewhere (60.00-61.00 Mb)") == 1
    assert "acrocentric karyotype" in it.basis(row(seq_flags="acrocentric"), ctx())
    assert "replaced a weak satellite array" in it.basis(row(seq_flags="weak_satellite_beside_core",
                                                             seq_evidence="retrotransposon", te_clade="CRM"), ctx())
    assert "replaced a satellite array elsewhere" in it.basis(row(seq_flags="satellite_away_from_core",
                                                                  seq_evidence="retrotransposon", te_clade="CRM"), ctx())


def test_moves_and_placements():
    r = row(call_source="hic", relation="offset", flags="moved_with_te", offset_bp=10_000_000, te_clade="CRM",
            support="te,hic", start=60_000_001, end=60_500_000)
    b = it.basis(r, ctx())
    assert b.startswith("moved 10.00 Mb from the sequence call (40.00-42.00 Mb) by a strong validated Hi-C peak")
    assert "second line: the CRM core at the peak" in b
    r = row(call_source="hic", relation="offset", flags=".", offset_bp=10_000_000, te_clade="CRM", support="te,hic",
            seq_confidence="medium")                    # not high: no moved_with flag, the core is still named
    assert "the CRM core at the peak agrees" in it.basis(r, ctx())
    r = row(call_source="hic", relation="offset", flags="moved_with_te+meth", offset_bp=10_000_000, te_clade="CRM")
    b = it.basis(r, ctx())
    assert "the CRM core at the peak" in b and "CHG methylation at the peak" in b
    r = row(call_source="hic", relation="offset", flags="fold_pair,moved_with_fold", offset_bp=500_000_000,
            support="hic,fold")
    b = it.basis(r, ctx())
    assert "the arm fold and Hi-C's prior-free peak" in b and "second line: the arm fold" in b
    assert "the arm fold agrees" not in b
    r = row(call_source="hic", relation="functional_only", seq_start=".", seq_end=".", seq_evidence="none",
            support="hic")
    assert it.basis(r, ctx()).startswith("placed by a strong validated Hi-C peak (no sequence call)")
    r = row(call_source="hic", relation="functional_only", flags="hic_medium_peak", seq_start=".", seq_end=".",
            confidence="low", support="hic")
    assert "a clean medium validated Hi-C peak" in it.basis(r, ctx())
    r = row(relation="concordant", flags="hic_chose_alt", seq_flags="competing_cluster", support="satellite,hic",
            seq_alt_start=70_000_001, seq_alt_end=72_000_000)
    assert "the second satellite cluster, chosen by validated Hi-C" in it.basis(r, ctx())
    assert "satellite cluster" not in it.caveat(r, ctx())


def test_oppositions():
    r = row(relation="contested", flags="hic_offset", against="hic", confidence="medium",
            hic_anch_peak_start=80_000_001, hic_anch_peak_end=80_250_000)
    c = it.caveat(r, ctx())
    assert c.startswith("validated Hi-C points elsewhere (80.00-80.25 Mb), with no second line there")
    assert c.count("Hi-C points elsewhere") == 1 and c.endswith("or CENP-A/CENH3 data, would decide")
    r = row(relation="contested", flags="fold_offset,hic_denovo_offset", against="fold,hic",
            fold_core_start=10_000_001, fold_core_end=15_000_000, hic_denovo_peak_start=12_000_001,
            hic_denovo_peak_end=12_500_000)
    c = it.caveat(r, ctx())
    assert "the arm fold points elsewhere (10.00-15.00 Mb)" in c
    assert "prior-free peak points elsewhere (12.00-12.50 Mb)" in c and c.count("points elsewhere") == 2
    r = row(flags="hic_disagrees_weak", hic_anch_peak_start=80_000_001, hic_anch_peak_end=80_250_000)
    assert "not strong enough to oppose" in it.caveat(r, ctx())
    assert it.caveat(row(flags="hic_near"), ctx()) == "."


def test_gaps_and_no_call():
    assert it.caveat(row(call_gaps=2), ctx()) == "2 assembly gaps inside the call"
    r = row(relation="none", start=".", end=".", confidence="none", seq_start=".", seq_end=".")
    assert it.basis(r, ctx()) == "."
    assert it.caveat(r, ctx()) == "no line found a centromere (satellite-free?)"
    assert it.caveat(r, ctx(hic=False, tes=False)).endswith("-> Hi-C or a CRM core (tes/; plants) can place one")
    assert it.caveat(r, ctx(hic=True)) == "no line found a centromere (satellite-free?)"
    assert it.caveat(SimpleNamespace(**{**vars(r), "hic_gate": "fail"}), ctx()).endswith("Hi-C failed its test here")


def test_3d_placed_calls_say_so():
    r = row(call_source="hic", relation="functional_only", seq_start=".", seq_end=".", seq_evidence="none",
            support="hic", confidence="medium")
    assert it.caveat(r, ctx()) == "rests on the 3D lines alone: no sequence or LTR evidence at this place"
    assert it.caveat(SimpleNamespace(**{**vars(r), "confidence": "high"}), ctx()) == "."


@pytest.mark.parametrize("bp,txt", [(999, "1 kb"), (60_000, "60 kb"), (999_499, "999 kb"), (1_000_000, "1.00 Mb"),
                                    (5_204_904, "5.20 Mb")])
def test_size(bp, txt):
    assert it.size(bp) == txt
