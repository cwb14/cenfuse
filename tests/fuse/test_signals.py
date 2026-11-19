"""Signals, their roles in each call, the verdict sentence and the genome-test table (cenfuse.signals). Rows are built
from minimal dicts holding only the columns a case uses (every other column '.'); the tomato and Poa annua rows and
the tomato summary values are copied from real runs."""
import re
from pathlib import Path

import pandas as pd
import pytest

from cenfuse import signals as sg

DOCS = Path(__file__).parents[2] / "docs"
TOMATO = dict(lines="satellite,te,hic", sat_one_per_chromosome="7/12 chromosomes hold the 12 largest clusters (null 8, "
              "P 0.87): fail", te_marker="Tekay", te_decisive="True", te_specificity="32.7",
              te_verdict="Tekay seed tomato_LTRs_fam00007, 32.7-fold in its cores, significant signature cores on 12/12 "
                         "chromosomes; decisive",
              at_gate=".", at_p=".", at_hic=".", rep_gate="fail", rep_p="0.537", fam_gate="pass", fam_p="0.0004",
              hic_bin="250000", hic_gate_pass="0", hic_gate_support="0", hic_gate_untested="0", hic_chromosomes="12",
              te_hic="1/12 chromosomes (P 0.46): fail", cis_line="peak at the core on 8/12 chromosomes (P 2.4e-05): pass",
              fold_gate="fail", fold_on_denovo="0/12", fold_p="1", align_gate="fail", align_on_denovo="0/12", align_p="1",
              inputs_reference=".", inputs_other=".", inputs_unknown="hic:hic")
POA = dict(lines="satellite,te,hic,meth", sat_one_per_chromosome="12/14 chromosomes hold the 14 largest clusters "
           "(null 8, P 0.0037): pass", te_marker="CRM", te_decisive="True", hic_gate_pass="14", hic_chromosomes="14",
           inputs_other="meth:C11,meth:I3,meth:J10,meth:P3,meth:R1", inputs_unknown="hic:hic")


def calls(*recs):
    """A centromeres table of minimal rows: chrom c1 of 100 Mb, every other column '.' unless given."""
    cols = set().union(*recs) | {"chrom", "chrom_len"}
    return pd.DataFrame([{**{c: "." for c in cols}, "chrom": "c1", "chrom_len": "100000000", **r} for r in recs])


def one(rec, summary=None, hic_bin=None):
    return sg.rows(calls(rec), summary, hic_bin)[0]


T067 = dict(chrom="CM175067.1", chrom_len="93972077", start="39995569", end="40160512", call_source="sequence",
            relation="concordant", confidence="high", support="te,hic", flags="cis_agrees", rule="cis_agrees",
            seq_start="39995569", seq_end="40160512", seq_evidence="long_unit_tandem",
            seq_flags="low_satellite,te_concordant", te_clade="Tekay", te_start="39656185", te_end="40971792",
            hic_gate="fail", cis_peak_start="40500001", cis_peak_end="40750000", hic_anch_core_start="65500001",
            hic_anch_core_end="65750000", hic_anch_level="none", hic_anch_flags="hub", rep_start="37207001",
            rep_end="38310000", fam_start="39328001", fam_end="41327000")
T068 = dict(chrom="CM175068.1", chrom_len="55854373", start="17250001", end="17500000", call_source="hic",
            relation="offset", confidence="low", support="hic", flags="cis_only", rule="cis_only",
            seq_start="19772317", seq_end="20442564", seq_evidence="retrotransposon", seq_flags="minority_core",
            seq_overruled_start="11042231", seq_overruled_end="11222177", te_clade="Tekay", te_start="19772317",
            te_end="20442564", hic_gate="fail", cis_peak_start="17250001", cis_peak_end="17500000",
            hic_anch_core_start="16500001", hic_anch_core_end="17750000", hic_anch_level="high")
T073 = dict(chrom="CM175073.1", chrom_len="69519740", start="27000001", end="27250000", call_source="hic",
            relation="offset", confidence="low", support="hic", flags="cis_only", rule="cis_only",
            seq_start="23226630", seq_end="23331600", seq_evidence="weak_satellite",
            seq_flags="low_satellite,te_concordant,te_age_core", te_clade="Tekay", te_start="22663358",
            te_end="24331825", hic_gate="fail", cis_peak_start="27000001", cis_peak_end="27250000",
            hic_anch_core_start="26250001", hic_anch_core_end="27750000", hic_anch_level="medium")
T074 = dict(chrom="CM175074.1", chrom_len="68864862", start="25500001", end="25750000", call_source="hic",
            relation="offset", confidence="low", support="hic", against="te", flags="cis_trans_pair",
            rule="cis_trans_pair", seq_start="18919178", seq_end="19036466", seq_evidence="weak_satellite",
            seq_flags="low_satellite,te_concordant", te_clade="Tekay", te_start="18180229", te_end="19695247",
            hic_gate="fail", cis_peak_start="25500001", cis_peak_end="25750000", fam_start="17894001",
            fam_end="19911000", hic_anch_core_start="24750001", hic_anch_core_end="25000000", hic_anch_level="high")
PA1A = dict(chrom="Pa1A", chrom_len="319656038", start="163024561", end="169417680", call_source="sequence",
            relation="concordant", confidence="high", support="te,hic", rule="hic_agrees", seq_start="163024561",
            seq_end="169417680", seq_evidence="retrotransposon", seq_flags="satellite_away_from_core,te_age_core",
            seq_overruled_start="146120146", seq_overruled_end="146762975", te_clade="CRM", te_start="163024561",
            te_end="169417680", hic_gate="pass", hic_anch_core_start="167500001", hic_anch_core_end="170000000",
            hic_anch_level="high")
PA1B = dict(chrom="Pa1B", chrom_len="97799933", start="38666466", end="39867844", call_source="sequence",
            relation="concordant", confidence="high", support="satellite,te,hic", rule="hic_agrees",
            seq_start="38666466", seq_end="39867844", seq_evidence="satellite", seq_flags="te_concordant",
            te_clade="CRM", te_start="38337405", te_end="40488993", hic_gate="pass", hic_anch_core_start="39500001",
            hic_anch_core_end="40000000", hic_anch_level="medium")
FASTA_ONLY = dict(start="40000001", end="42000000", call_source="sequence", relation="sequence_only",
                  confidence="low", support="satellite", rule="sequence", seq_start="40000001", seq_end="42000000",
                  seq_evidence="satellite")
NO_CALL = dict(call_source="none", relation="none", confidence="none", rule="none")


def test_places_are_half_open():
    r = one(dict(FASTA_ONLY, seq_start="4000001", seq_end="5000000"))
    assert r.places["sat"] == (4_000_000, 5_000_000)


def test_hic_token_names_the_readout():
    cis = one(dict(call_source="hic", start="1000001", end="1250000", support="hic", rule="cis_only", flags="cis_only"),
              TOMATO, 250_000)
    assert cis.roles["flank"] == "used"
    agrees = one(dict(FASTA_ONLY, support="satellite,hic", rule="hic_agrees", hic_gate="pass"), POA, 500_000)
    assert agrees.roles["hic"] == "used"
    hub = one(dict(call_source="hic", start="1000001", end="1250000", support="hic", rule="hic_denovo_only",
                   flags="hic_denovo_only"), {"hic_chromosomes": "3", "hic_gate_untested": "3"}, 50_000)
    assert hub.roles["hub"] == "used"


def test_roles_tomato_073():
    r = one(T073, TOMATO, 250_000)
    assert (r.roles["flank"], r.roles["te"], r.roles["sat"]) == ("used", "elsewhere", "not")


def test_methylation_from_other_plants_is_not_used():
    summary = {"inputs_other": "meth:C11"}
    r = one(FASTA_ONLY, summary)
    assert r.roles["meth"] == "not"
    assert sg.role_label("meth", "not", summary) == "not used (other plants)"


@pytest.mark.parametrize("rec,summary,hic_bin,text", [
    (T067, TOMATO, 250_000, "Tekay and Hi-C agree, on a long tandem repeat"),
    (T068, TOMATO, 250_000, "Placed by Hi-C alone; Tekay cluster 2.3 Mb away disagrees"),
    (T074, TOMATO, 250_000, "Placed by two Hi-C measures; Tekay cluster 5.8 Mb away disagrees"),
    (PA1A, POA, 500_000, "CRM cluster, backed by Hi-C; satellite array 16.3 Mb away overruled"),
    (PA1B, POA, 500_000, "Satellite, CRM and Hi-C agree"),
    (FASTA_ONLY, {"lines": "satellite"}, None, "Satellite array alone"),
    (NO_CALL, {"lines": "satellite"}, None, "No centromere called"),
], ids=["tomato067", "tomato068", "tomato074", "Pa1A", "Pa1B", "fasta_only", "no_call"])
def test_verdict_examples(rec, summary, hic_bin, text):
    assert one(rec, summary, hic_bin).verdict == text


def _rule_names() -> list:
    """The `rule` values docs/outputs.md lists: the backticked names of its **`rule`.** paragraph that start an item
    (followed by a comma or a parenthesis), not the words in backticks inside an item's explanation."""
    text = (DOCS / "outputs.md").read_text()
    para = text.split("**`rule`.**", 1)[1].split("\n\n", 1)[0]
    return re.findall(r"`([a-z_]+)`(?=,|\s\()", para)


# each rule's minimal row: where the call comes from, the lines that agree and that point elsewhere, its flags
RULE_ROWS = {
    "sequence": ("sequence", "satellite", ".", "."),
    "hic_agrees": ("sequence", "satellite,hic", ".", "."),
    "hic_support": ("sequence", "satellite,hic", ".", "hic_support"),
    "hic_end": ("sequence", "satellite,hic", ".", "hic_end"),
    "hic_chose_alt": ("sequence", "satellite,hic", ".", "hic_chose_alt"),
    "hic_moved": ("hic", "hic", ".", "."),
    "hic_offset": ("sequence", "satellite", "hic", "hic_offset"),
    "hic_strong_peak": ("hic", "hic", ".", "."),
    "hic_medium_peak": ("hic", "hic", ".", "hic_medium_peak"),
    "hic_denovo_only": ("hic", "hic", ".", "hic_denovo_only"),
    "fold_pair": ("hic", "hic,fold", ".", "fold_pair"),
    "align_pair": ("hic", "hic,align", ".", "align_pair"),
    "fold_offset": ("sequence", "satellite", "fold", "fold_offset"),
    "align_offset": ("sequence", "satellite", "align", "align_offset"),
    "hic_denovo_offset": ("sequence", "satellite", "hic", "hic_denovo_offset"),
    "te_hic_pair": ("hic", "te,hic", ".", "te_hic_pair"),
    "cis_agrees": ("sequence", "satellite,hic", ".", "cis_agrees"),
    "cis_only": ("hic", "hic", ".", "cis_only"),
    "cis_elsewhere": ("sequence", "satellite", "hic", "cis_elsewhere"),
    "cis_trans_pair": ("hic", "hic", ".", "cis_trans_pair"),
    "cis_end_pair": ("hic", "hic", ".", "cis_end_pair"),
    "at_only": ("composition", "at", ".", "at_only"),
    "at_denovo_pair": ("composition", "at,hic", ".", "at_denovo_pair"),
    "at_pair": ("hic", "at,hic", ".", "hic_denovo_only,at_pair"),
    "at_chose_alt": ("sequence", "at", ".", "at_chose_alt"),
    "rep_only": ("composition", "repeats", ".", "rep_only"),
    "rep_pair": ("hic", "repeats,hic", ".", "cis_only,rep_pair"),
    "fam_only": ("composition", "families", ".", "fam_only"),
    "fam_over_low": ("composition", "families", ".", "fam_over_low"),
    "cpg_elsewhere": ("sequence", "satellite", "cpg", "cpg_elsewhere"),
    "meth_support": ("sequence", "satellite,meth", ".", "."),
}
FAR = ("5000001", "5500000")      # every signal's own place, 3.5 Mb from the call
PLACES = ("te", "hic_anch_core", "hic_denovo_peak", "cis_peak", "fold_core", "align_core", "cpg_alt", "at", "rep", "fam")


def minimal(rule: str) -> dict:
    if rule == "none":
        return NO_CALL
    src, sup, ag, flags = RULE_ROWS[rule]
    rec = dict(start="1000001", end="1500000", call_source=src, relation="concordant", confidence="low", support=sup,
               against=ag, flags=flags, rule=rule, te_clade="CRM")
    if src == "sequence":
        rec.update(seq_start="1000001", seq_end="1500000", seq_evidence="satellite")
    if rule == "sequence":                    # a second cluster the call did not take
        rec.update(seq_alt_start=FAR[0], seq_alt_end=FAR[1])
    if rule.endswith("_chose_alt"):           # the call took the second cluster; the first lies away
        rec.update(seq_start=FAR[0], seq_end=FAR[1], seq_alt_start="1000001", seq_alt_end="1500000")
    for p in PLACES:
        rec[f"{p}_start"], rec[f"{p}_end"] = FAR
    return rec


def test_every_rule_has_a_verdict():
    names = _rule_names()
    assert len(names) >= 30 and {"sequence", "none", "cis_trans_pair", "fam_only"} <= set(names)
    for rule in names:
        v = one(minimal(rule), {}).verdict
        assert v and not sg.BANNED.search(v), (rule, v)


def test_labels_have_no_banned_words():
    for s in sg.SIGNALS.values():
        assert not sg.BANNED.search(s.label("CRM")), s.label("CRM")


def test_distance_text():
    assert sg.distance_text(350_000) == "350 kb"
    assert sg.distance_text(2_270_000) == "2.3 Mb"


def test_signal_tests_tomato():
    tests = {t.key: t for t in sg.signal_tests(TOMATO, 12, "Tekay")}
    s, te, fl, fam = tests["sat"], tests["te"], tests["flank"], tests["fam"]
    assert (s.k, s.n, s.chance, s.p, s.verdict) == (7, 12, 8.0, 0.87, "not trusted here")
    assert (te.k, te.n, te.verdict) == (12, 12, "trusted")
    assert (fl.k, fl.n, fl.p, fl.verdict) == (8, 12, 2.4e-05, "trusted")
    assert fam.verdict == "a lead only"
    # Hi-C's and CpG's rows count the chromosomes whose own test passed, and say so
    assert tests["hic"].question == "Contacts the other centromeres: its test passes on 0 of 12 chromosomes"
    assert tests["hic"].per_chrom and not tests["te"].per_chrom
    cpg = {t.key: t for t in sg.signal_tests(dict(TOMATO, cpg_gate_pass="3"), 12, "Tekay")}["cpg"]
    assert (cpg.k, cpg.n) == (3, 12) and cpg.question.endswith(": its test passes on 3 of 12 chromosomes")
    assert [t.key for t in tests.values()] == [k for k in sg.SIGNALS if k in tests]       # registry order
    for t in tests.values():
        assert not sg.BANNED.search(f"{t.label} {t.question}"), t


def test_rows_without_summary():
    rec = dict(FASTA_ONLY, support="satellite", against="te", te_start="70000001", te_end="71000000",
               at_start="60000001", at_end="61000000", hic_gate="pass", hic_anch_core_start="80000001",
               hic_anch_core_end="80250000")
    for summary in (None, {}):
        r = one(rec, summary)
        assert (r.roles["sat"], r.roles["te"]) == ("used", "elsewhere")
        assert r.roles["at"] == r.roles["hic"] == "not"           # far places, but no genome test to trust them by
        tests = sg.signal_tests(summary, 1, "LTR")
        assert tests and {t.verdict for t in tests} == {"not reported"}


AXIS = dict(POA, fold_gate="pass", align_gate="pass")       # both passed; Hi-C's own test passed too, so neither acted
PA2A = dict(chrom="Pa2A", chrom_len="266166041", start="167613549", end="169150825", call_source="sequence",
            relation="concordant", confidence="high", support="satellite,te,hic", rule="hic_agrees",
            seq_start="167613549", seq_end="169150825", seq_evidence="satellite", seq_flags="te_concordant",
            te_clade="CRM", te_start="166619917", te_end="169281576", hic_gate="pass", hic_anch_core_start="166500001",
            hic_anch_core_end="169000000", hic_anch_level="high", fold_core_start="129250001",
            fold_core_end="161750000", align_core_start="167000001", align_core_end="184000000")
PA4B = dict(chrom="Pa4B", chrom_len="72579290", start="30760327", end="32776876", call_source="sequence",
            relation="concordant", confidence="high", support="satellite,te,hic", rule="hic_agrees",
            seq_start="30760327", seq_end="32776876", seq_evidence="satellite", seq_flags="te_concordant",
            te_clade="CRM", te_start="30773497", te_end="32805708", hic_gate="pass", hic_anch_core_start="32500001",
            hic_anch_core_end="33500000", hic_anch_level="medium", fold_core_start="29000001",
            fold_core_end="39000000", align_core_start="20000001", align_core_end="23000000")


FOLD_PAIR = dict(chrom="c2", start="1000001", end="1500000", call_source="hic", relation="functional_only",
                 confidence="medium", support="hic,fold", rule="fold_pair", flags="fold_pair", fold_core_start="1000001",
                 fold_core_end="1500000")                     # the fold and the hub placed this call


def test_hic_readouts_point_elsewhere_only_where_cenfuse_counted_them():
    """The arm fold and alignment passed their tests in this genome but lie 5.9 and 7.8 Mb from the calls: not used,
    even where the fold placed a call on another chromosome; where cenfuse counted the fold against a call, it points
    elsewhere."""
    for r in sg.rows(calls(PA2A, PA4B), AXIS, 500_000):
        assert r.roles["fold"] == r.roles["align"] == "not"
        assert r.verdict == "Satellite, CRM and Hi-C agree"
    r = sg.rows(calls(PA2A, FOLD_PAIR), AXIS, 500_000)[0]        # the fold acted elsewhere in this genome, not here
    assert r.roles["fold"] == "not" and r.verdict == "Satellite, CRM and Hi-C agree"
    r = one(dict(PA2A, against="fold", flags="fold_offset"), AXIS, 500_000)
    assert r.roles["fold"] == "elsewhere"
    assert r.verdict == "Satellite, CRM and Hi-C agree; Hi-C arm fold 5.9 Mb away disagrees"


TRUSTING = dict(POA, cis_line="peak at the core on 12/14 chromosomes (P 2.4e-05): pass", at_gate="pass",
                rep_gate="pass")                             # every signal here passed its genome test
FAR_HIC = dict(FASTA_ONLY, hic_gate="pass", hic_anch_core_start="80000001", hic_anch_core_end="80250000")


@pytest.mark.parametrize("peak", [dict(hic_anch_level="low"), dict(hic_anch_level="high", hic_anch_flags="telomeric"),
                                  dict(hic_anch_level="none"), dict()], ids=["weak", "telomeric", "none", "missing"])
def test_a_far_hic_peak_cenfuse_did_not_count_says_nothing(peak):
    """A Hi-C peak 38 Mb from the call that cenfuse did not count against it (a weak or telomeric peak, or no peak
    recorded at all), and trusted flank contact, AT-rich trough and repeat density far from the call: not used, no
    tail. Without a recorded peak Hi-C has no place (no bar)."""
    rec = dict(FAR_HIC, **peak, cis_peak_start="60000001", cis_peak_end="60250000", at_start="70000001",
               at_end="70100000", rep_start="90000001", rep_end="91000000")
    r = one(rec, TRUSTING, 500_000)
    assert {k: r.roles[k] for k in ("hic", "flank", "at", "rep")} == dict.fromkeys(("hic", "flank", "at", "rep"), "not")
    assert r.verdict == "Satellite array alone"
    peaked = peak.get("hic_anch_level") in ("low", "high")
    assert r.places["hic"] == ((80_000_000, 80_250_000) if peaked else None)


def test_a_hic_peak_cenfuse_counted_against_the_call_points_elsewhere():
    r = one(dict(FAR_HIC, against="hic", flags="hic_offset", rule="hic_offset", hic_anch_level="high"), POA, 500_000)
    assert r.roles["hic"] == "elsewhere"
    assert r.verdict == "Satellite array alone; Hi-C peak 38.0 Mb away disagrees"


def test_present_tracks():
    """Tracks top to bottom in TRACK_ORDER: the Hi-C readouts (the trace first), methylation, the LTR family, the
    satellite arrays, the composition signals."""
    fasta = sg.rows(calls(FASTA_ONLY), {"lines": "satellite"}, None)
    assert sg.present(fasta, {"lines": "satellite"}, {}) == ["sat"]
    tomato = sg.rows(calls(T067, T068, T073, T074), TOMATO, 250_000)
    assert sg.present(tomato, TOMATO, {"te": True, "hic": True}) == ["hic", "flank", "te", "sat"]
    has = {"te": True, "hic": True, "meth": True}
    poa = sg.rows(calls(PA2A, PA4B), AXIS, 500_000)               # fold and alignment passed, never consulted
    assert sg.present(poa, AXIS, has) == ["hic", "meth", "te", "sat"]
    folded = sg.rows(calls(PA2A, PA4B, FOLD_PAIR), AXIS, 500_000)
    assert sg.present(folded, AXIS, has) == ["hic", "fold", "hub", "meth", "te", "sat"]
    hub_only = sg.rows(calls(FASTA_ONLY), {"lines": "satellite,hic"}, 500_000)    # no anchored Hi-C place anywhere
    assert sg.present(hub_only, {"lines": "satellite,hic"}, {"hic": True}) == ["hub", "sat"]   # the hub is the trace


FAM_ONLY = dict(chrom="c3", start="20000001", end="21000000", call_source="composition", relation="functional_only",
                confidence="low", support="families", rule="fam_only", flags="fam_only", fam_start="20000001",
                fam_end="21000000")


def test_repeat_family_cluster_track_only_where_it_acted():
    """The repeat-family cluster's test passes often by chance (a lead only): like the secondary Hi-C readouts, it gets
    a track only in a genome where it placed or opposed a call. Tomato: fam_gate pass, a cluster on every chromosome,
    never used -> no track (it drew a 'not used' bar on all 12 tiles)."""
    has = {"te": True, "hic": True}
    tomato = sg.rows(calls(T067, T068, T073, T074), TOMATO, 250_000)
    assert sg.present(tomato, TOMATO, has) == ["hic", "flank", "te", "sat"]
    assert all(r.roles["fam"] == "not" for r in tomato)
    acted = sg.rows(calls(T067, T068, FAM_ONLY), TOMATO, 250_000)
    assert acted[2].roles["fam"] == "used" and "fam" in acted[2].consulted
    assert sg.present(acted, TOMATO, has) == ["hic", "flank", "te", "sat", "fam"]


def test_placed_by_names_the_signal_that_placed_the_call():
    """The signals that placed a call, not those that agree with it: Hi-C agrees with every Poa annua call but placed
    none; tomato's Hi-C calls were placed by the flank contact; a long tandem repeat is no satellite array, so no
    signal of the figure placed a call kept on one."""
    assert one(PA1B, POA, 500_000).placed == {"sat"}
    assert one(PA1A, POA, 500_000).placed == {"te"}
    assert one(T068, TOMATO, 250_000).placed == {"flank"}
    assert one(T067, TOMATO, 250_000).placed == frozenset()
    assert sg.rows(calls(T067, FAM_ONLY), TOMATO, 250_000)[1].placed == {"fam"}
    assert one(NO_CALL, {}).placed == frozenset()


def test_the_family_line_is_a_lead_through_its_pair_test_only_where_that_test_acts():
    """0.19: where the family test failed but its windows passed the test against the repeat-density peaks, the figure shows
    the family line as a lead, with the pair test's question and count, only where that test lets the line act on some
    chromosome (fam_rep_acts); a genome-wide pass that leaving one chromosome out undoes everywhere (S. pombe) stays as the
    family test said. A call placed through the pair is headed as two readouts agreeing, not as a cluster alone."""
    base = dict(lines="satellite", fam_gate="fail", fam_p="0.0149")
    assert sg._status("fam", dict(base, fam_rep="7/7 pass", fam_rep_acts="7")) == "lead"
    assert sg._status("fam", dict(base, fam_rep="3/3 pass", fam_rep_acts="0")) == "fail"
    assert sg._status("fam", dict(base, fam_rep="9/25 support", fam_rep_acts="0")) == "fail"
    assert sg._status("fam", dict(base, fam_gate="pass")) == "lead"
    row = next(r for r in sg.signal_tests(dict(base, fam_rep="7/7 pass", fam_rep_acts="7"), 7, "CRM") if r.key == "fam")
    assert row.verdict == "a lead only" and "densest repeats" in row.question and (row.k, row.n) == (7, 7)
    rec = dict(call_source="composition", relation="sequence_only", confidence="low", support="families,repeats",
               flags="fam_over_low,fam_rep_pair", rule="fam_over_low", fam_start="10000001", fam_end="16000000",
               rep_start="12000001", rep_end="13000000", start="10000001", end="16000000")
    r = one(rec, dict(base, fam_rep="7/7 pass", fam_rep_acts="7"))
    assert sg.verdict(rec, r).startswith("Repeat-family cluster and densest repeats agree (a lead)")
