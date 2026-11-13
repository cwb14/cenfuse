"""basis and caveat sentences (cenfuse.interpret): one test per relation and per flag the fusion rules or the
satellite engine set."""
from types import SimpleNamespace

import pytest

from cenfuse import interpret as it

BASE = dict(chrom="c1", chrom_len=100_000_000, start=40_000_001, end=42_000_000, length=2_000_000,
            call_source="sequence", relation="sequence_only", confidence="high", support="satellite", against=".",
            flags=".", offset_bp=".", seq_start=40_000_001, seq_end=42_000_000, seq_evidence="satellite",
            seq_confidence="high", seq_flags=".", seq_alt_start=".", seq_alt_end=".", seq_overruled_start=".",
            seq_overruled_end=".", main_family="g_sat01_171",
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


def moved(**kw):
    """A call the 3D lines moved, `medium`: the row the caveat of a moved call reads."""
    return row(**{**dict(call_source="hic", relation="offset", confidence="medium", offset_bp=5_000_000, support="hic",
                         seq_start=30_000_001, seq_end=30_100_000), **kw})


def test_caveat_of_a_moved_call_says_no_evidence_only_where_there_is_none():
    r = moved()
    assert "no sequence or LTR evidence at this place" in it.caveat(r, ctx())
    assert "rests on the 3D lines alone" in it.caveat(r, ctx())


def test_caveat_of_a_moved_call_reads_the_satellite_its_basis_reports():
    """Moss's moved calls: basis says '21 kb of centromeric satellite inside', the caveat said there was none."""
    r = moved(support="hic,satellite")
    c = ctx(sat_bp=lambda chrom, s, e: 21_000)
    assert "21 kb of centromeric satellite inside" in it.basis(r, c)
    cav = it.caveat(r, c)
    assert "no sequence or LTR evidence" not in cav and "alone" not in cav
    assert "21 kb of centromeric satellite" in cav


def test_caveat_of_a_moved_call_reads_a_core_or_methylation_at_the_peak():
    cav = it.caveat(moved(support="hic,te", te_clade="CRM"), ctx())
    assert "no sequence or LTR evidence" not in cav and "the CRM core" in cav
    cav = it.caveat(moved(support="hic,meth"), ctx())
    assert "no sequence or LTR evidence" not in cav and "CHG methylation" in cav


def test_no_call_where_only_an_unsupported_long_unit_array_stands_says_so():
    r = row(relation="none", start=".", end=".", confidence="none", call_source="none", flags="long_unit_unsupported",
            seq_evidence="long_unit_tandem", seq_confidence="low", support=".", main_family="long_unit")
    c = it.caveat(r, ctx(hic=False, tes=False))
    assert "long-unit tandem array (40.00-42.00 Mb)" in c and "no line" in c
    assert "satellite-free" not in c and c.endswith("-> Hi-C or an LTR core (tes/; plants) can place one")
    assert it.basis(r, ctx()) == "."


def test_a_chain_family_is_not_called_a_satellite():
    r = row(seq_flags="partial_chain", confidence="low", main_family="g_sat01_1790", support="satellite")
    c = ctx(unit={"g_sat01_1790": 1790})
    b = it.basis(r, c)
    assert not b.startswith("satellite") and "g_sat01_1790" in b and "partial copies" in b
    assert "not a tandem satellite" in it.caveat(r, c)


def test_a_chain_caveat_says_what_would_settle_it():
    r = row(seq_flags="partial_chain", confidence="low", support="satellite")
    assert it.caveat(r, ctx(tes=False)).endswith("-> an LTR core (tes/; plants) or CENP-A/CENH3 data would decide")
    assert it.caveat(r, ctx(tes=True)).endswith("-> CENP-A/CENH3 data would decide")


def test_a_family_rejected_by_the_crm_line_is_named():
    r = row(seq_evidence="retrotransposon", te_clade="CRM", support="te", main_family=".", confidence="medium",
            seq_rejected_family="g_sat02_170", seq_rejected_start=10_000_001, seq_rejected_end=11_000_000)
    cav = it.caveat(r, ctx(decisive=True))
    assert "the CRM line rejected g_sat02_170" in cav and "10.00-11.00 Mb" in cav
    assert it.caveat(row(), ctx()) == "."                     # older tables: no such column


def test_copy_spacing_is_named_when_it_differs_from_the_unit():
    c = ctx()
    assert "171-bp unit, copies spaced 155 bp here" in it.basis(row(seq_spacing=155), c)
    assert "copies spaced" not in it.basis(row(seq_spacing=170), c)                      # within LEN_TOL of the unit
    assert "copies spaced" not in it.basis(row(seq_spacing="."), c)
    assert "copies spaced" not in it.basis(row(), c)                                     # older tables: no such column


def test_sequence_evidence_kinds():
    assert it.basis(row(seq_evidence="weak_satellite"), ctx()).startswith("weak satellite g_sat01_171")
    assert it.basis(row(seq_evidence="long_unit_tandem", main_family="long_unit"), ctx()).startswith(
        "long-unit tandem array")
    assert "not a satellite family" in it.caveat(row(seq_evidence="long_unit_tandem"), ctx())
    r = row(seq_evidence="retrotransposon", te_clade="CRM", support="te", main_family=".")
    assert it.basis(r, ctx(decisive=True)) == "CRM retrotransposon core (no satellite call)"
    assert it.caveat(r, ctx(decisive=True)) == "."
    r = row(seq_evidence="retrotransposon", te_clade="Bryco", support=".", confidence="low", main_family=".")
    assert "this lineage does not decide" in it.caveat(r, ctx(decisive=False))


def test_satellite_engine_flags():
    r = row(seq_flags="competing_cluster", confidence="low", seq_alt_start=70_000_001, seq_alt_end=72_000_000)
    c = it.caveat(r, ctx())
    assert "another satellite cluster (70.00-72.00 Mb) competes" in c
    assert c.endswith("-> CENP-A/CENH3 data would decide")
    assert "or Hi-C when mappable DNA separates them" in it.caveat(r, ctx(hic=False))
    assert "or an LTR core (tes/; plants)" in it.caveat(r, ctx(tes=False))
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
    b = it.basis(row(seq_flags="satellite_beside_core", seq_evidence="retrotransposon", te_clade="CRM"), ctx())
    assert b == ("CRM retrotransposon core; the CRM core replaced a satellite array beside it that holds "
                 "significantly fewer CRM elements")
    for flag in ("weak_satellite_beside_core", "satellite_away_from_core", "satellite_beside_core"):
        assert "no satellite call" not in it.basis(row(seq_flags=flag, seq_evidence="retrotransposon", te_clade="CRM"),
                                                   ctx()), flag


def test_overruled_array_is_named():
    """The three sentences that say a CRM core replaced a satellite array end with where that array was, in the position
    style of this module (1-based table columns -> 'a-b Mb'). Without an interval (a table made before the columns
    existed, or no array to name) the sentence is as it was."""
    sentences = {"weak_satellite_beside_core": "the CRM core replaced a weak satellite array beside it",
                 "satellite_away_from_core": "the CRM core replaced a satellite array elsewhere",
                 "satellite_beside_core": "the CRM core replaced a satellite array beside it that holds "
                                          "significantly fewer CRM elements"}
    for flag, text in sentences.items():
        kw = dict(seq_flags=flag, seq_evidence="retrotransposon", te_clade="CRM")
        r = row(seq_overruled_start=21_123_181, seq_overruled_end=21_170_900, **kw)
        assert it.basis(r, ctx()) == f"CRM retrotransposon core; {text} (21.12-21.17 Mb)", flag
        assert it.basis(row(**kw), ctx()) == f"CRM retrotransposon core; {text}", flag          # "."
        older = SimpleNamespace(**{k: v for k, v in vars(row(**kw)).items() if not k.startswith("seq_overruled")})
        assert it.basis(older, ctx()) == f"CRM retrotransposon core; {text}", flag              # no such column


@pytest.mark.parametrize("flag, clade, decisive", [
    ("satellite_away_from_core", "CRM", True), ("satellite_beside_core", "CRM", True),
    ("weak_satellite_beside_core", "CRM", True),
    (".", "CRM", True),                 # a long-unit call that disagrees with a decisive core: replaced, no flag
    (".", "Athila", False),             # ... or with a flags-only marker's core
    (".", ".", False)])                 # a TE line not given in this run (--satellite-dir of a run with tes/)
def test_caveat_names_the_overruled_array(flag, clade, decisive):
    """The array a TE core replaced is known against the call, so the caveat opens with it and its interval, in words
    true for every case that fills seq_overruled_*: a satellite, weak-satellite or long-unit call that a decisive CRM
    core replaced, and a long-unit call that a core of any clade replaced. The rest of the caveat is unchanged. Without
    the columns ('.', or a table made before they existed) the caveat is as it was; the basis keeps its interval."""
    kw = dict(seq_flags=flag, seq_evidence="retrotransposon", te_clade=clade, support="te", main_family=".",
              confidence="medium", seq_confidence="medium", te_start=40_000_001, te_end=42_000_000)
    c = ctx(decisive=decisive, marker=clade)
    before = it.caveat(row(**kw), c)                                                    # seq_overruled_* = "."
    older = SimpleNamespace(**{k: v for k, v in vars(row(**kw)).items() if not k.startswith("seq_overruled")})
    assert it.caveat(older, c) == before and "overruled" not in before
    r = row(seq_overruled_start=21_123_181, seq_overruled_end=21_170_900, **kw)
    clause = f"the tandem array at 21.12-21.17 Mb was overruled by the {clade if clade != '.' else 'TE'} core"
    assert it.caveat(r, c) == (clause if before == "." else f"{clause}; {before}")
    assert "(21.12-21.17 Mb)" in it.basis(r, c) or flag == "."


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
    c = it.caveat(row(flags="hic_near", confidence="low", hic_anch_peak_start=80_000_001, hic_anch_peak_end=80_250_000), ctx())
    assert "validated Hi-C's peak (80.00-80.25 Mb) lies beside the call" in c          # why Hi-C did not raise a low call


def test_gaps_and_no_call():
    assert it.caveat(row(call_gaps=2), ctx()) == "2 assembly gaps inside the call"
    r = row(relation="none", start=".", end=".", confidence="none", seq_start=".", seq_end=".")
    assert it.basis(r, ctx()) == "."
    assert it.caveat(r, ctx()) == "no line found a centromere (satellite-free?)"
    assert it.caveat(r, ctx(hic=False, tes=False)).endswith("-> Hi-C or an LTR core (tes/; plants) can place one")
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


def test_support_tier_says_support_only_and_no_failure():
    r = row(relation="concordant", confidence="medium", support="satellite,hic", flags="hic_support",
            hic_gate="support")
    b = it.basis(r, ctx())
    assert "Hi-C agrees (support only: its peaks find the other chromosomes' centromeres beyond chance" in b
    assert "validated Hi-C agrees" not in b
    assert "Hi-C failed its test here" not in it.caveat(r, ctx())


def test_hic_alone_says_untested_and_what_would_confirm():
    r = row(call_source="hic", relation="functional_only", confidence="low", support="hic", flags="hic_denovo_only",
            seq_start=".", seq_end=".", seq_evidence="none", seq_confidence="none", hic_gate="untested")
    b = it.basis(r, ctx())
    assert "placed by Hi-C's prior-free peak alone (too few reference calls in this genome to test Hi-C against)" in b
    c = it.caveat(r, ctx())
    assert "rests on the 3D lines alone" in c and c.endswith("-> CENP-A/CENH3 data would confirm")


def test_untested_prior_free_peak_elsewhere_is_named():
    r = row(confidence="low", flags="hic_denovo_elsewhere", hic_gate="untested", hic_denovo_peak_start=30_000_001,
            hic_denovo_peak_end=30_010_000)
    c = it.caveat(r, ctx())
    assert "Hi-C's prior-free peak, untested here, lies elsewhere (30.00-30.01 Mb)" in c
    assert c.endswith("-> CENP-A/CENH3 data would decide")


def test_a_call_placed_by_an_at_trough_says_so():
    r = row(call_source="composition", seq_evidence="none", at_gate="pass", at_gc=0.12, support="at", main_family=".",
            confidence="low", relation="sequence_only")
    assert "deepest AT trough" in it.basis(r, ctx()) and "12%" in it.basis(r, ctx())
    assert "where the sequence is AT-rich, not where the kinetochore sits" in it.caveat(r, ctx())
    assert it.caveat(r, ctx(hic=False)).endswith("-> Hi-C would test it")


def test_a_trough_agreeing_with_hic_is_named():
    r = row(call_source="hic", relation="functional_only", flags="hic_denovo_only,at_pair", at_gc=0.08, support="hic,at",
            confidence="medium", seq_evidence="none", fold_gate=".", align_gate=".")
    assert "an AT trough at the peak agrees" in it.basis(r, ctx())
    assert "an AT trough" in it.caveat(r, ctx())


def test_a_trough_elsewhere_is_said():
    r = row(call_source="hic", relation="functional_only", flags="hic_denovo_only,at_elsewhere", against="at", at_start=500_001,
            at_end=502_000, support="hic", confidence="low", seq_evidence="none", fold_gate=".", align_gate=".")
    assert "the AT trough lies elsewhere (0.50-0.50 Mb)" in it.caveat(r, ctx())


def test_a_pair_of_core_and_hic_is_named():
    r = row(call_source="hic", relation="offset", flags="te_hic_pair", te_clade="CRM", support="te,hic", offset_bp=2_000_000,
            confidence="medium", seq_evidence="satellite")
    assert "the CRM core and Hi-C's peak" in it.basis(r, ctx())
    k = row(relation="concordant", flags="te_hic_pair", te_clade="CRM", support="satellite,te,hic")
    assert "Hi-C's peak lies at the CRM core" in it.basis(k, ctx())


def test_a_weak_core_call_says_why_it_is_low():
    r = row(seq_evidence="retrotransposon", te_clade="Tekay", support="te", main_family=".", confidence="low",
            te_seed="fam7", te_core_elements=6, te_chrom_elements=23)
    assert "seeded by fam7" in it.basis(r, ctx(decisive=True)) and "6 of its 23 elements" in it.basis(r, ctx(decisive=True))
    assert "under half of the signature's elements" in it.caveat(r, ctx(decisive=True))


# ---------------------------------------------------------------- 0.16: the genome-level verdict when nothing is confident
def test_the_verdict_names_a_genome_where_no_line_passed_a_test():
    from cenfuse import interpret
    f = dict(called=6, high=0, medium=0, low=6)
    holo = {"sat_one_per_chromosome": "5/6 chromosomes hold the 6 largest clusters (null 4, P 0.24): fail", "te_decisive": ".",
            "at_gate": ".", "hic_gate_pass": "0", "hic_gate_support": "0", "te_hic": ".",
            "cis_line": "trans peaks at the cis peak on 0/6 chromosomes (P 1): fail, tested by trans contacts; acts on 0", "hic_axis": "."}
    v = interpret.verdict(f, holo, dict(hic=True, tes=True))
    assert "no line of evidence passed a test" in v and "holocentric" in v and "could test" not in v
    v = interpret.verdict(f, holo, dict(hic=False, tes=False))
    assert "Hi-C or an LTR annotation (tes/; plants) could test them" in v
    assert interpret.verdict(f, dict(holo, at_gate="pass"), dict(hic=False, tes=False)) is None       # P. falciparum's troughs
    cis = "trans peaks at the cis peak on 12/14 chromosomes (P 1.4e-13): pass, tested by trans contacts; acts on 14"
    assert interpret.verdict(f, dict(holo, cis_line=cis), dict(hic=True, tes=False)) is None          # C. neoformans in G1
    assert interpret.verdict(dict(f, medium=1, low=5), holo, dict(hic=True, tes=True)) is None
    none = interpret.verdict(dict(called=0, high=0, medium=0, low=0), holo, dict(hic=False, tes=True))
    assert none.startswith("- **No centromere placed**") and "Hi-C could test them" in none


def test_passed_lines_reads_every_test_of_the_summary():
    from cenfuse import interpret
    s = {"sat_one_per_chromosome": "22/24 chromosomes hold the 24 largest clusters (null 15, P 0.0001): pass", "te_decisive": "True",
         "at_gate": "pass", "hic_gate_pass": "3", "te_hic": "8/12 chromosomes (P 1e-4): support", "cis_line": "x: pass",
         "hic_axis": "fold", "cpg_gate_pass": "20"}
    assert len(interpret.passed_lines(s)) == 8
    assert interpret.passed_lines({}) == []
