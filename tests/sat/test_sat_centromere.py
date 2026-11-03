import pandas as pd
import pysam

from cenfuse.sat.annotate import ARR_COLS
from cenfuse.sat.centromere import finalize, load_regions, longunit_call, satellite_call, write_regions
from cenfuse.sat.longunit import LU_COLS
from cenfuse.sat.params import Params


def test_satellite_call_main_cluster_and_competing():
    a = pd.DataFrame([(0, "c", 5_000_000, 5_200_000, "+", "A", 1100, 0.9, 180),
                      (1, "c", 5_300_000, 5_350_000, "+", "A", 270, 0.9, 180),
                      (2, "c", 8_000_000, 8_120_000, "+", "A", 660, 0.9, 180)], columns=ARR_COLS)
    call = satellite_call(a, {"A"}, 10_000_000, Params())
    assert (call["cen_start"], call["cen_end"]) == (5_000_000, 5_350_000)
    assert call["sat_bp"] == 250_000 and call["second_cluster_bp"] == 120_000 and call["evidence"] == "satellite"


def test_satellite_call_weak_and_absent():
    p = Params()
    small = pd.DataFrame([(0, "c", 100_000, 108_000, "+", "A", 44, 0.9, 180)], columns=ARR_COLS)
    assert satellite_call(small, {"A"}, 10_000_000, p)["evidence"] == "weak_satellite"
    assert satellite_call(small, {"B"}, 10_000_000, p) is None


def test_finalize_flags_gaps_and_near_end(tmp_path):
    seq = "A" * 1000 + "N" * 20 + "A" * 1000 + "N" * 5 + "A" * 2000
    fa = tmp_path / "x.fa"
    fa.write_text(">c\n" + seq + "\n")
    pysam.faidx(str(fa))
    p = Params(subtel_frac_len=0.25)
    call = dict(cen_start=100, cen_end=3000, evidence="satellite", main_family="A", main_family_bp=2000,
                sat_bp=2000, second_cluster_bp=0, cen_families="A")
    with pysam.FastaFile(str(fa)) as h:
        row = finalize("c", len(seq), call, h, None, None, p)
    assert row["n_gaps"] == 1 and row["gap_bp"] == 20
    assert row["flags"] == "low_satellite,near_end,gaps:1" and row["confidence"] == "low"


def test_none_row_when_no_call():
    row = finalize("c", 1000, None, None, None, None, Params())
    assert row["evidence"] == "none" and row["confidence"] == "none" and row["cen_start"] == "."


def test_longunit_call_requires_concentration():
    p = Params()
    lu = pd.DataFrame([("c", 4_000_000, 4_600_000, 2700, 222.2, 0.95)], columns=LU_COLS)
    call = longunit_call(lu, 10_000_000, None, None, p)
    assert call["evidence"] == "long_unit_tandem" and call["cen_start"] == 4_000_000
    scattered = pd.DataFrame([("c", s, s + 30_000, 2700, 11.1, 0.9) for s in range(1_000_000, 9_000_000, 1_000_000)],
                             columns=LU_COLS)
    assert longunit_call(scattered, 10_000_000, None, None, p) is None


def test_regions_roundtrip(tmp_path):
    rows = {"c1": finalize("c1", 1000, None, None, None, None, Params())}
    rows["c2"] = dict(rows["c1"], chrom="c2", cen_start=10, cen_end=500, evidence="satellite", confidence="high")
    rows["c3"] = dict(rows["c2"], chrom="c3", evidence="retrotransposon", overruled_start=100, overruled_end=300,
                      overruled_family="fam001")
    write_regions(rows, tmp_path / "r.tsv")
    back = load_regions(tmp_path / "r.tsv")
    assert back["c1"]["evidence"] == "none" and back["c2"]["cen_end"] == 500 and back["c2"]["chrom_len"] == 1000
    over = ("overruled_start", "overruled_end", "overruled_family")
    assert [back["c3"][k] for k in over] == ["100", "300", "fam001"]            # the replaced array survives
    assert all(back[c][k] == "." for c in ("c1", "c2") for k in over)           # nothing replaced: absent
    # a regions file written before these columns existed (0.8) loads, and its rows read as absent
    older = pd.read_csv(tmp_path / "r.tsv", sep="\t", dtype=str).drop(columns=list(over))
    older.to_csv(tmp_path / "older.tsv", sep="\t", index=False)
    back = load_regions(tmp_path / "older.tsv")
    assert back["c3"]["cen_end"] == 500 and all(back["c3"].get(k, ".") == "." for k in over)


def test_diffuse_region_is_flagged_and_not_high(tmp_path):
    """Real data (Coffea): arrays scattered through pericentromeric heterochromatin chain into 5-22 Mb regions
    that are only 6-9% satellite; such a region must not be reported as a high-confidence centromere."""
    fa = tmp_path / "x.fa"
    fa.write_text(">c\n" + "A" * 20_000_000 + "\n")
    pysam.faidx(str(fa))
    call = dict(cen_start=5_000_000, cen_end=12_000_000, evidence="satellite", main_family="A",
                main_family_bp=700_000, sat_bp=700_000, second_cluster_bp=0, cen_families="A")
    with pysam.FastaFile(str(fa)) as h:
        row = finalize("c", 20_000_000, call, h, None, None, Params())
    assert "diffuse" in row["flags"] and row["confidence"] == "medium"


def _contigs(tmp_path):
    from cenfuse.sat.genomes import Genome
    fa = tmp_path / "ctg.fa"
    fa.write_text("".join(f">{c}\n{'A' * 200_000}\n" for c in ("c1", "c2", "c3")))
    pysam.faidx(str(fa))
    arrays = pd.DataFrame([(0, "c1", 0, 60_000, "+", "A", 333, 0.9, 180),
                           (1, "c2", 50_000, 58_000, "+", "A", 44, 0.9, 180)], columns=ARR_COLS)
    return Genome("ctg", fa, {"c1": 200_000, "c2": 200_000, "c3": 200_000}), arrays


def test_borrowed_calls_need_satellite_evidence(tmp_path, monkeypatch):
    """On a borrowed (contig) genome only full satellite calls count: a contig without the centromeric
    satellite is no evidence of a satellite-free centromere, so no weak-satellite or long-unit fallback, and
    the long-unit scan covers call flanks only."""
    from cenfuse.sat import centromere, longunit
    from cenfuse.sat.seeds import COLS
    g, arrays = _contigs(tmp_path)
    seen = []
    lu = pd.DataFrame([("c3", 80_000, 140_000, 3000, 20.0, 0.95)], columns=LU_COLS)
    monkeypatch.setattr(longunit, "detect", lambda fa, targets, *a: seen.append(targets) or lu)
    p = Params()
    trf = pd.DataFrame(columns=COLS)
    rows, _ = centromere.call_regions(g, arrays, {"A"}, None, trf, tmp_path / "t1", 1, p)
    assert [rows[c]["evidence"] for c in ("c1", "c2", "c3")] == ["satellite", "weak_satellite", "long_unit_tandem"]
    rows, _ = centromere.call_regions(g, arrays, {"A"}, None, trf, tmp_path / "t2", 1, p, borrowed=True)
    assert [rows[c]["evidence"] for c in ("c1", "c2", "c3")] == ["satellite", "none", "none"]
    assert seen[-1] == [("c1", 0, 200_000)]


def test_borrowed_contig_end_flag_does_not_lower_confidence(tmp_path):
    g, _ = _contigs(tmp_path)
    call = dict(cen_start=0, cen_end=60_000, evidence="satellite", main_family="A", main_family_bp=60_000,
                sat_bp=60_000, second_cluster_bp=0, cen_families="A")
    with pysam.FastaFile(str(g.fasta)) as h:
        own = finalize("c1", 200_000, call, h, None, None, Params())
        bor = finalize("c1", 200_000, call, h, None, None, Params(), borrowed=True)
    assert own["flags"] == "near_end" and own["confidence"] == "medium"
    assert bor["flags"] == "contig_end" and bor["confidence"] == "high"


def test_without_a_dominant_satellite_or_te_evidence_calls_are_capped_low(tmp_path):
    """Aegilops: satellite-free centromeres, but interstitial satellite blocks pass every per-family test. With no
    dominant family and no retrotransposon evidence, satellite evidence cannot tell them apart: confidence low."""
    g, _ = _contigs(tmp_path)
    call = dict(cen_start=60_000, cen_end=120_000, evidence="satellite", main_family="A", main_family_bp=60_000,
                sat_bp=60_000, second_cluster_bp=0, cen_families="A")
    with pysam.FastaFile(str(g.fasta)) as h:
        p = Params(low_sat_bp=10_000)
        free = finalize("c1", 200_000, call, h, None, None, p)
        capped = finalize("c1", 200_000, call, h, None, None, p, capped=True)
        disc = finalize("c1", 200_000, call, h, None, None, p, te_flag="te_discordant")
    assert free["confidence"] == "high" and free["flags"] == "."
    assert capped["confidence"] == "low" and capped["flags"] == "no_dominant_satellite"
    assert disc["confidence"] == "medium" and disc["flags"] == "te_discordant"
    weak = dict(call, evidence="weak_satellite", sat_bp=8_000)
    with pysam.FastaFile(str(g.fasta)) as h:
        conc = finalize("c1", 200_000, weak, h, None, None, p, te_flag="te_concordant")
    assert conc["confidence"] == "medium" and "te_concordant" in conc["flags"]          # was low on its own



def test_a_flags_only_marker_lifts_the_cap_only_for_calls_it_supports(tmp_path):
    """No dominant satellite: the cap holds unless retrotransposon evidence speaks for the call. A non-decisive marker
    (Athila, Tekay: flags only) supports only te_concordant calls; a call it contradicts (te_discordant) or does not
    reach keeps the cap. A decisive (CRM) marker lifts it genome-wide, as before."""
    from cenfuse.sat.centromere import honesty_cap
    assert honesty_cap(None, False, True, False) and not honesty_cap(None, False, True, True)       # no TE evidence
    assert honesty_cap(("Athila", {}, False), False, True, False)                                     # flags only
    assert not honesty_cap(("CRM", {}, True), False, True, False)                                     # decisive
    assert not honesty_cap(None, True, True, False) and not honesty_cap(None, False, False, False)
    g, _ = _contigs(tmp_path)
    call = dict(cen_start=60_000, cen_end=120_000, evidence="satellite", main_family="A", main_family_bp=60_000,
                sat_bp=60_000, second_cluster_bp=0, cen_families="A")
    with pysam.FastaFile(str(g.fasta)) as h:
        p = Params(low_sat_bp=10_000)
        kw = dict(capped=True, te_decisive=False)
        conc = finalize("c1", 200_000, call, h, None, None, p, te_flag="te_concordant", **kw)
        disc = finalize("c1", 200_000, call, h, None, None, p, te_flag="te_discordant", **kw)
        none = finalize("c1", 200_000, call, h, None, None, p, **kw)
    assert conc["confidence"] == "high" and "no_dominant_satellite" not in conc["flags"]
    assert disc["confidence"] == "low" and "no_dominant_satellite" in disc["flags"] and "te_discordant" in disc["flags"]
    assert none["confidence"] == "low" and none["flags"] == "no_dominant_satellite"


def test_a_decisive_core_replaces_a_weak_array_beside_it(tmp_path, monkeypatch):
    """A weak satellite call that agrees with a decisive (CRM) core without overlapping it is replaced by the core
    (S. verrucosum chr01, Ae. geniculata chr07U, potato chr11: the core on CENH3, the weak array 0.8-5 Mb away). A
    strong array beside the core is kept here (te_concordant: these cores carry no element midpoints, so the
    head-to-head test of the next tests cannot run), and a flags-only marker replaces nothing. The replaced array is
    kept as overruled_*; a chromosome that had no call takes its core and overrules nothing."""
    from cenfuse.sat import centromere, longunit
    from cenfuse.sat.seeds import COLS
    g, arrays = _contigs(tmp_path)
    monkeypatch.setattr(longunit, "detect", lambda *a: pd.DataFrame(columns=LU_COLS))
    cores = {"c1": dict(start=70_000, end=100_000), "c2": dict(start=70_000, end=120_000),
             "c3": dict(start=70_000, end=100_000)}                       # c3: a core and no call at all
    trf = pd.DataFrame(columns=COLS)
    rows, _ = centromere.call_regions(g, arrays, {"A"}, None, trf, tmp_path / "t1", 1, Params(),
                                      te=("CRM", cores, True))
    assert rows["c1"]["evidence"] == "satellite" and "te_concordant" in rows["c1"]["flags"]
    assert rows["c3"]["evidence"] == "retrotransposon" and rows["c3"]["flags"] == "."
    assert rows["c2"]["evidence"] == "retrotransposon" and "weak_satellite_beside_core" in rows["c2"]["flags"]
    assert (rows["c2"]["cen_start"], rows["c2"]["cen_end"]) == (70_000, 120_000) and rows["c2"]["confidence"] == "medium"
    over = ("overruled_start", "overruled_end", "overruled_family")
    assert [rows["c2"][k] for k in over] == [50_000, 58_000, "A"]               # the weak call the core replaced
    # never in alt_*: fusion lets validated Hi-C move a call to the alternative cluster, which would change calls
    assert rows["c2"]["alt_start"] == rows["c2"]["alt_end"] == "."
    assert all(rows[c][k] == "." for c in ("c1", "c3") for k in over)           # a kept call; the core of no call
    rows, _ = centromere.call_regions(g, arrays, {"A"}, None, trf, tmp_path / "t2", 1, Params(),
                                      te=("Athila", cores, False))
    assert rows["c2"]["evidence"] == "weak_satellite" and "te_concordant" in rows["c2"]["flags"]
    assert all(rows[c][k] == "." for c in rows for k in over)                   # a flags-only marker replaces nothing



def test_a_decisive_core_replaces_a_confident_call_only_if_the_clade_marks_it_more(tmp_path, monkeypatch):
    """A distant decisive (CRM) core replaces a confident satellite call only when the clade's elements sit in the
    core significantly more often than in the call (te.contradicts); otherwise the call stays, te_discordant, one
    level lower (B. rapa CC1/CC2 A09 and soybean X26 chr1/chr6: CENH3 on the satellite, the core window in a broad
    CRM plateau; P. annua Pa1A and Ae. geniculata 7M/2U: the core on the centromere)."""
    import numpy as np
    from cenfuse.sat import centromere, longunit
    from cenfuse.sat.seeds import COLS
    from cenfuse.sat.te import Mids
    g, arrays = _contigs(tmp_path)
    monkeypatch.setattr(longunit, "detect", lambda *a: pd.DataFrame(columns=LU_COLS))
    trf = pd.DataFrame(columns=COLS)
    in_call = list(np.linspace(5_000, 55_000, 10).astype(int))
    alone = centromere.call_regions(g, arrays, {"A"}, None, trf, tmp_path / "t0", 1, Params())[0]["c1"]["confidence"]
    level = {"high": "medium", "medium": "low", "low": "low"}
    plateau = dict(start=150_000, end=170_000, k=12, mids=Mids(sorted(in_call + list(range(151_000, 170_000, 1_600))), 22))
    rows, _ = centromere.call_regions(g, arrays, {"A"}, None, trf, tmp_path / "t1", 1, Params(),
                                      te=("CRM", {"c1": plateau}, True))
    assert rows["c1"]["evidence"] == "satellite" and "te_discordant" in rows["c1"]["flags"]
    assert rows["c1"]["confidence"] == level[alone]
    over = ("overruled_start", "overruled_end", "overruled_family")
    assert all(rows["c1"][k] == "." for k in over)                      # kept (te_discordant): nothing overruled
    marked = dict(start=150_000, end=170_000, k=40, mids=Mids(sorted(in_call[:2] + list(range(150_200, 170_000, 500))), 42))
    rows, _ = centromere.call_regions(g, arrays, {"A"}, None, trf, tmp_path / "t2", 1, Params(),
                                      te=("CRM", {"c1": marked}, True))
    assert rows["c1"]["evidence"] == "retrotransposon" and "satellite_away_from_core" in rows["c1"]["flags"]
    assert (rows["c1"]["cen_start"], rows["c1"]["cen_end"]) == (150_000, 170_000)
    assert [rows["c1"][k] for k in over] == [0, 60_000, "A"]                    # replaced: the satellite call
    rows, _ = centromere.call_regions(g, arrays, {"A"}, None, trf, tmp_path / "t3", 1, Params(),
                                      te=("Athila", {"c1": marked}, False))
    assert "te_discordant" in rows["c1"]["flags"] and all(rows["c1"][k] == "." for k in over)   # flags only: kept


def test_a_decisive_core_replaces_a_distant_call_that_is_low_on_its_own(tmp_path, monkeypatch):
    """The head-to-head test guards confident calls only: a distant decisive core replaces a call that is `low` on its
    own (here: less centromeric satellite than low_sat_bp) even where the clade does not mark the core more."""
    import numpy as np
    from cenfuse.sat import centromere, longunit
    from cenfuse.sat.seeds import COLS
    from cenfuse.sat.te import Mids
    g, arrays = _contigs(tmp_path)
    monkeypatch.setattr(longunit, "detect", lambda *a: pd.DataFrame(columns=LU_COLS))
    trf = pd.DataFrame(columns=COLS)
    p = Params(low_sat_bp=100_000)                     # the 60-kb array is low on its own (low_satellite)
    assert centromere.call_regions(g, arrays, {"A"}, None, trf, tmp_path / "t0", 1, p)[0]["c1"]["confidence"] == "low"
    in_call = list(np.linspace(5_000, 55_000, 10).astype(int))
    plateau = dict(start=150_000, end=170_000, k=12,
                   mids=Mids(sorted(in_call + list(range(151_000, 170_000, 1_600))), 22))   # 12 vs 10: P 0.42
    rows, _ = centromere.call_regions(g, arrays, {"A"}, None, trf, tmp_path / "t1", 1, p,
                                      te=("CRM", {"c1": plateau}, True))
    assert rows["c1"]["evidence"] == "retrotransposon" and "satellite_away_from_core" in rows["c1"]["flags"]
    assert [rows["c1"][k] for k in ("overruled_start", "overruled_end", "overruled_family")] == [0, 60_000, "A"]


def test_a_confident_call_beside_a_decisive_core_moves_only_if_the_clade_marks_the_core(tmp_path, monkeypatch):
    """A strong satellite call that agrees with a decisive (CRM) core without overlapping it moves to the core when
    the clade's elements sit there significantly more (satellite_beside_core; Ae. tauschii Chr2/Chr4: 332 vs 11 and
    292 vs 6 CRM elements, CENH3 inside the core); otherwise it stays, te_concordant (B. rapa A06: 21-23 vs 43-51)."""
    import numpy as np
    from cenfuse.sat import centromere, longunit
    from cenfuse.sat.seeds import COLS
    from cenfuse.sat.te import Mids
    g, arrays = _contigs(tmp_path)
    monkeypatch.setattr(longunit, "detect", lambda *a: pd.DataFrame(columns=LU_COLS))
    trf = pd.DataFrame(columns=COLS)
    in_call = list(np.linspace(5_000, 55_000, 10).astype(int))
    plain = dict(start=70_000, end=100_000, k=12,
                 mids=Mids(sorted(in_call + list(np.linspace(71_000, 99_000, 12).astype(int))), 22))
    rows, _ = centromere.call_regions(g, arrays, {"A"}, None, trf, tmp_path / "t1", 1, Params(),
                                      te=("CRM", {"c1": plain}, True))
    assert rows["c1"]["evidence"] == "satellite" and "te_concordant" in rows["c1"]["flags"]
    over = ("overruled_start", "overruled_end", "overruled_family")
    assert all(rows["c1"][k] == "." for k in over)                              # kept: nothing overruled
    rich = dict(start=70_000, end=100_000, k=60,
                mids=Mids(sorted(in_call[:2] + list(np.linspace(70_500, 99_500, 60).astype(int))), 62))
    rows, _ = centromere.call_regions(g, arrays, {"A"}, None, trf, tmp_path / "t2", 1, Params(),
                                      te=("CRM", {"c1": rich}, True))
    assert rows["c1"]["evidence"] == "retrotransposon" and "satellite_beside_core" in rows["c1"]["flags"]
    assert (rows["c1"]["cen_start"], rows["c1"]["cen_end"]) == (70_000, 100_000)
    assert [rows["c1"][k] for k in over] == [0, 60_000, "A"]                    # replaced: the satellite call beside it


def test_a_long_unit_call_replaced_by_a_core_is_overruled(tmp_path, monkeypatch):
    """A long-unit tandem call that a decisive (CRM) core replaces is recorded as overruled like a satellite call (the
    earlier audit recorded satellite arrays only): one beside the core (weak_satellite_beside_core), and one far from
    it, which the core replaces without a flag. A long-unit call that overlaps the core is kept (te_concordant)."""
    from cenfuse.sat import centromere, longunit
    from cenfuse.sat.seeds import COLS
    g, _ = _contigs(tmp_path)
    lu = pd.DataFrame([("c1", 80_000, 140_000, 3000, 20.0, 0.95), ("c2", 20_000, 60_000, 3000, 13.0, 0.95),
                       ("c3", 80_000, 140_000, 3000, 20.0, 0.95)], columns=LU_COLS)
    monkeypatch.setattr(longunit, "detect", lambda *a: lu)
    cores = {"c1": dict(start=150_000, end=170_000), "c2": dict(start=150_000, end=170_000),
             "c3": dict(start=100_000, end=120_000)}
    rows, _ = centromere.call_regions(g, pd.DataFrame(columns=ARR_COLS), {"A"}, None, pd.DataFrame(columns=COLS),
                                      tmp_path / "t1", 1, Params(), te=("CRM", cores, True))
    over = ("overruled_start", "overruled_end", "overruled_family")
    assert rows["c1"]["evidence"] == "retrotransposon" and "weak_satellite_beside_core" in rows["c1"]["flags"]
    assert [rows["c1"][k] for k in over] == [80_000, 140_000, "long_unit"]
    assert rows["c2"]["evidence"] == "retrotransposon" and (rows["c2"]["cen_start"], rows["c2"]["cen_end"]) == \
        (150_000, 170_000)
    assert [rows["c2"][k] for k in over] == [20_000, 60_000, "long_unit"]
    assert rows["c3"]["evidence"] == "long_unit_tandem" and "te_concordant" in rows["c3"]["flags"]
    assert all(rows["c3"][k] == "." for k in over)
    assert all(rows[c]["alt_start"] == rows[c]["alt_end"] == "." for c in rows)
    rows, _ = centromere.call_regions(g, pd.DataFrame(columns=ARR_COLS), {"A"}, None, pd.DataFrame(columns=COLS),
                                      tmp_path / "t2", 1, Params(), te=("CRM", cores, True), borrowed=True)
    assert all(rows[c][k] == "." for c in rows for k in over)                   # borrowed: no long-unit call at all
