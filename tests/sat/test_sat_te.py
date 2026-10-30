import math

import numpy as np
import pytest

from cenfuse.sat.te import Mids, agrees, binom_sf, centromeric_clade, contradicts, core, load_te_gff, marker

L = 10_000_000


def test_binom_sf_exact_values():
    assert binom_sf(10, 10, 0.5) == pytest.approx(1 / 1024)
    assert binom_sf(0, 10, 0.3) == 1.0
    assert binom_sf(3, 4, 0.5) == pytest.approx(5 / 16)


def test_load_te_gff_reads_clades_from_ltrquest_and_edta_style_lines(tmp_path):
    gff = tmp_path / "g.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "c1\tLTRquest\tLTR_retrotransposon\t101\t200\t.\t+\t.\tID=a;superfamily=Gypsy;clade=CRM;family=f1\n"
        "c1\tLTRquest\tlong_terminal_repeat\t101\t120\t.\t+\t.\tParent=a\n"
        "c1\tEDTA\tGypsy_LTR_retrotransposon\t1001\t2000\t.\t-\t.\tID=b;Clade=CRM\n"
        "c1\tEDTA\tCopia_LTR_retrotransposon\t5001\t6000\t.\t-\t.\tID=c;Clade=Ale\n"
        "c1\tLTRquest\tLTR_retrotransposon\t7001\t8000\t.\t+\t.\tID=d;clade=unknown\n"
        "c1\tLTRquest\tLTR_retrotransposon\t9001\t9100\t.\t+\t.\tID=e\n"
        "c9\tLTRquest\tLTR_retrotransposon\t1\t100\t.\t+\t.\tID=f;clade=CRM\n")
    t = load_te_gff(gff, {"c1": 20_000})
    assert set(t) == {"CRM", "Ale"} and t["CRM"]["c1"].tolist() == [150, 1500] and t["Ale"]["c1"].tolist() == [5500]


DANTE = ("##gff-version 3\n"
         "c1\tdante_ltr\ttransposable_element\t101\t200\t.\t+\t.\tID=TE_1;Final_Classification=Class_I|LTR|Ty3/gypsy|"
         "chromovirus|CRM;Rank=DLTP\n"
         "c1\tdante_ltr\tlong_terminal_repeat\t101\t120\t.\t+\t.\tFinal_Classification=Class_I|LTR|Ty3/gypsy|chromovirus|"
         "CRM;Parent=TE_1\n"
         "c1\tdante\tprotein_domain\t130\t160\t5\t+\t.\tFinal_Classification=Class_I|LTR|Ty3/gypsy|chromovirus|CRM;"
         "Parent=TE_1\n"
         "c1\tdante_ltr\ttransposable_element\t1001\t2000\t.\t-\t.\tID=TE_2;Final_Classification=Class_I|LTR|Ty3/gypsy|"
         "non-chromovirus|OTA|Athila;Rank=D\n"
         "c1\tdante_ltr\ttransposable_element\t3001\t4000\t.\t-\t.\tID=TE_3;Final_Classification=Class_I|LTR|Ty3/gypsy|"
         "chromovirus\n"
         "c1\tdante_ltr\ttransposable_element\t5001\t6000\t.\t+\t.\tID=TE_4;Final_Classification=Class_I|LTR|Ty1/copia|"
         "Ivana\n")


def test_dante_ltr_elements_take_their_lineage(tmp_path):
    from cenfuse.methfuse import read_tes
    from cenfuse.sat.te import dante_lineage, load_te_intervals
    assert dante_lineage("Class_I|LTR|Ty3/gypsy|chromovirus|CRM") == ("CRM", "Gypsy")
    assert dante_lineage("Class_I|LTR|Ty3/gypsy|non-chromovirus|OTA|Tat|Retand") == ("Retand", "Gypsy")
    assert dante_lineage("Class_I|LTR|Ty3/gypsy|chromovirus") == ("", "Gypsy")        # classified to a rank only
    gff = tmp_path / "d.gff3"
    gff.write_text(DANTE)
    t = load_te_gff(gff, {"c1": 20_000})        # parent elements only; the unresolved one and the rank-D Athila (protein
    assert {c: v["c1"].tolist() for c, v in t.items()} == {"CRM": [150], "Ivana": [5500]}    # domains, no LTRs) left out
    s, e, cl = load_te_intervals(gff, {"c1": 20_000})["c1"]
    assert s.tolist() == [100, 5000] and cl == ["CRM", "Ivana"]
    dl = tmp_path / "dl.gff3"
    dl.write_text(DANTE.replace("Athila;Rank=D", "Athila;Rank=DL"))
    assert load_te_gff(dl, {"c1": 20_000})["Athila"]["c1"].tolist() == [1500]
    m = read_tes(gff)                                         # the CHG line's classes: lineage, else superfamily
    assert m.cls.tolist() == ["LTR:CRM", "LTR:Athila", "LTR:Gypsy_unclassified", "LTR:Ivana"]
    assert m.depth.tolist() == [0, 0, 0, 0]


def test_core_finds_a_concentrated_cluster_and_ignores_uniform_spread():
    rng = np.random.default_rng(1)
    clustered = np.sort(np.concatenate([rng.integers(4_000_000, 4_400_000, 60), rng.integers(0, L, 40)]))
    c = core(clustered, L)
    assert c["significant"] and c["start"] <= 4_200_000 <= c["end"] and c["end"] - c["start"] <= 1_000_000
    assert c["k"] >= 55
    assert not core(np.sort(rng.integers(0, L, 100)), L)["significant"]
    assert core(clustered[:2], L) is None                               # too few sites to estimate a peak
    # Sphagnum-like: 8 of 9 elements of a clade in one 0.2-Mb cluster (moss genomes carry few intact LTR-RTs)
    sparse = np.sort(np.concatenate([rng.integers(6_000_000, 6_200_000, 8), [1_000_000]]))
    assert core(sparse, L)["significant"]



def test_core_floor_counts_independent_insertion_sites(tmp_path):
    """Nested (overlapping) elements are one insertion site: one nested cluster cannot make a core on its own, three
    independent insertions at one locus can (Chlorella NS4-2: 2-8 Osser elements per chromosome, all at CENH3)."""
    gff = tmp_path / "g.gff3"
    rows = [(5_000_000 + i * 50, 5_008_000 + i * 50) for i in range(6)] + [(1_000_000, 1_005_000)]
    gff.write_text("##gff-version 3\n" + "".join(
        f"c1\tLTRquest\tLTR_retrotransposon\t{s}\t{e}\t.\t+\t.\tID=e{i};clade=Osser\n" for i, (s, e) in enumerate(rows)))
    m = load_te_gff(gff, {"c1": L})["Osser"]["c1"]
    assert isinstance(m, Mids) and len(m) == 7 and m.sites == 2
    assert core(m, L) is None                                           # 7 elements but 2 sites: below the floor
    c = core(Mids([5_000_000, 5_020_000, 5_040_000], 3), L)             # 3 independent insertions at one locus
    assert c is not None and c["start"] <= 5_020_000 <= c["end"] and c["k"] == 3
    assert core(np.array([5_000_000, 5_020_000]), L) is None            # plain arrays: elements count


def test_a_core_contradicts_a_call_only_if_the_clade_marks_it_significantly_more():
    """Elements in the core vs in the call, binomial at 1/2: B. rapa A09-like (core window 18, call 16) is no
    contradiction; P. annua Pa1A-like (core 62, call on a second cluster of 19) is; no element list: as before."""
    call = (40_000_000, 45_000_000)
    inside = list(np.linspace(40_100_000, 44_900_000, 16).astype(int))
    a09 = dict(start=25_000_000, end=26_000_000, k=18,
               mids=Mids(sorted(inside + list(np.linspace(25_050_000, 25_950_000, 18).astype(int))), 34))
    assert not contradicts(*call, a09)
    inside = list(np.linspace(40_100_000, 44_900_000, 19).astype(int))
    pa1a = dict(start=60_000_000, end=66_000_000, k=62,
                mids=Mids(sorted(inside + list(np.linspace(60_100_000, 65_900_000, 62).astype(int))), 81))
    assert contradicts(*call, pa1a)
    assert contradicts(*call, {k: v for k, v in a09.items() if k != "mids"})


def test_centromeric_clade_is_the_concentrated_one_on_most_chromosomes():
    rng = np.random.default_rng(2)
    chroms = {f"c{i}": L for i in range(5)}
    centro = {c: np.sort(np.concatenate([rng.integers(i * 1_000_000 + 2_000_000, i * 1_000_000 + 2_300_000, 50),
                                         rng.integers(0, L, 30)])) for i, c in enumerate(chroms)}
    spread = {c: np.sort(rng.integers(0, L, 200)) for c in chroms}
    clade, cores = centromeric_clade({"CRM": centro, "Ale": spread}, chroms)
    assert clade == "CRM" and all(cores[c]["significant"] for c in chroms)
    assert centromeric_clade({"Ale": spread}, chroms) == (None, {})
    two = {c: v for c, v in list(centro.items())[:2]}                     # concentrated on only 2 of 5 chromosomes
    assert centromeric_clade({"CRM": two}, chroms)[0] is None


def test_agrees_within_the_wider_of_call_and_core():
    c = dict(start=1_000_000, end=1_500_000)                             # 0.5-Mb core
    assert agrees(1_200_000, 1_300_000, c)                               # inside
    assert agrees(1_900_000, 2_000_000, c)                               # 0.4 Mb away: within the core width
    assert agrees(3_000_000, 6_000_000, c)                               # 1.5 Mb away: within the 3-Mb call width
    assert not agrees(2_100_000, 2_200_000, c)                           # 0.6 Mb away, both narrower


def test_marker_prefers_the_centromere_specific_crm_lineage():
    """CRM (CR chromoviruses: CRM maize, CRR rice, CRB Brassica, CRS sugarcane) is centromere-specific in plants;
    other clades are often pericentromeric. A qualifying CRM is the marker even if another clade is more concentrated,
    and only CRM evidence is decisive."""
    a, b = {"c1": dict(start=1, end=2)}, {"c1": dict(start=5, end=6)}
    assert marker([("Ale", a), ("CRM", b)]) == ("CRM", b, True)
    assert marker([("Retand", a), ("Athila", b)]) == ("Retand", a, False)
    assert marker([]) is None


def test_ltr_overlap_reports_tandem_repeats_inside_retrotransposons(tmp_path):
    """Musa: a 62-bp repeat lies inside SIRE elements; a true satellite array lies outside them."""
    import pandas as pd
    from cenfuse.sat.te import load_te_intervals, ltr_overlap
    gff = tmp_path / "g.gff3"
    gff.write_text("##gff-version 3\n"
                   "c1\tx\tLTR_retrotransposon\t1001\t9000\t.\t+\t.\tID=a;clade=SIRE\n"
                   "c1\tx\tLTR_retrotransposon\t20001\t28000\t.\t+\t.\tID=b;clade=SIRE\n"
                   "c1\tx\tLTR_retrotransposon\t24001\t26000\t.\t+\t.\tID=c\n")       # nested, no clade
    iv = load_te_intervals(gff, {"c1": 100_000})
    assert iv["c1"][2] == ["SIRE", "SIRE"]                                      # unclassified element left out
    arrays = pd.DataFrame(dict(chrom=["c1"] * 3, start=[2000, 21000, 50000], end=[2500, 21500, 60000],
                               family=["T", "T", "S"]))
    ov = ltr_overlap(arrays, iv, longest=8000)
    assert ov["T"] == (1.0, "SIRE") and ov["S"] == (0.0, ".")


def test_census_tells_unclassified_and_misnamed_te_files(tmp_path):
    """An EDTA GFF3 without TEsorter lineages, or a file for another naming, gives no clade: the census says why."""
    from cenfuse.sat.te import te_census
    gff = tmp_path / "edta.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "c1\tEDTA\tGypsy_LTR_retrotransposon\t1001\t2000\t.\t-\t.\tID=b;Classification=LTR/Gypsy\n"
        "c1\tEDTA\tLTR_retrotransposon\t3001\t4000\t.\t-\t.\tID=c;Classification=LTR/unknown\n"
        "c1\tEDTA\tlong_terminal_repeat\t3001\t3100\t.\t-\t.\tParent=c\n"
        "Chr2\tEDTA\tCopia_LTR_retrotransposon\t5001\t6000\t.\t-\t.\tID=d;clade=Ale\n"
        "c1\tEDTA\thelitron\t7001\t8000\t.\t-\t.\tID=e\n")
    assert te_census(gff, {"c1": 20_000, "c2": 20_000}) == dict(ltr=2, classified=0, foreign=1)


def test_te_line_warns_when_no_element_has_a_lineage(tmp_path):
    import io
    import logging
    from cenfuse.teline import te_line
    gff = tmp_path / "edta.gff3"
    gff.write_text("c1\tEDTA\tGypsy_LTR_retrotransposon\t1001\t2000\t.\t-\t.\tID=b;Classification=LTR/Gypsy\n")
    buf, lg = io.StringIO(), logging.getLogger("cenfuse")
    h = logging.StreamHandler(buf)
    lg.addHandler(h)
    try:
        out = te_line(gff, {"c1": 20_000})
    finally:
        lg.removeHandler(h)
    assert out["marker"] is None and "none classified to a lineage" in buf.getvalue()


def test_te_annotation_named_gff_is_found(tmp_path):
    from cenfuse.sat.te import te_gff
    (tmp_path / "tes").mkdir()
    (tmp_path / "tes" / "g.gff").write_text("c1\tx\tLTR_retrotransposon\t1\t100\t.\t+\t.\tID=a;clade=CRM\n")
    assert te_gff(tmp_path, "g").name == "g.gff"


def test_crm_cores_take_precedence_on_their_chromosomes():
    """With a less sensitive annotation CRM concentrates on only some chromosomes (DANTE_LTR rice: 5/12, all on
    CENH3) while a pericentromeric lineage qualifies on most (Tekay 11/12, all off CENH3). The marker stays Tekay, but on
    the chromosomes where CRM has a significant core, that core is used (flags only); a CRM majority stays decisive."""
    from cenfuse.sat.te import candidate_clades, choose_marker
    rng = np.random.default_rng(3)
    chroms = {f"c{i:02d}": L for i in range(12)}
    peak = lambda a: np.sort(np.concatenate([rng.integers(a, a + 300_000, 40), rng.integers(0, L, 20)]))
    tekay = {c: peak(6_000_000) for c in chroms}                         # pericentromeric, every chromosome
    crm3 = {c: (peak(3_000_000) if i < 3 else np.sort(rng.integers(0, L, 30))) for i, c in enumerate(chroms)}
    tes = {"Tekay": tekay, "CRM": crm3}
    cands = candidate_clades(tes, chroms)
    assert [c for c, _ in cands] == ["Tekay"]                            # CRM fails the majority rule
    m = choose_marker(tes, chroms, cands)
    assert m[0] == "Tekay" and not m[2] and len(m[1]) == 12
    assert [c for c, x in m[1].items() if x.get("clade") == "CRM"] == ["c00", "c01", "c02"]
    assert all(m[1][c]["start"] < 3_300_000 and m[1][c]["end"] > 3_000_000 for c in ("c00", "c01", "c02"))
    assert all(m[1][c]["start"] < 6_300_000 and m[1][c]["end"] > 6_000_000 for c in ("c03", "c11"))
    crm_all = {c: peak(3_000_000) for c in chroms}
    t2 = {"Tekay": tekay, "CRM": crm_all}
    m = choose_marker(t2, chroms, candidate_clades(t2, chroms))
    assert m[0] == "CRM" and m[2] and len(m[1]) == 12 and not any("clade" in x for x in m[1].values())
    only = {"CRM": crm3}                                                 # no clade qualifies: the CRM cores alone
    m = choose_marker(only, chroms, candidate_clades(only, chroms))
    assert m[0] == "CRM" and not m[2] and sorted(m[1]) == ["c00", "c01", "c02"]
    assert choose_marker({"Tekay": tekay}, chroms, candidate_clades({"Tekay": tekay}, chroms))[0] == "Tekay"
    assert choose_marker({}, chroms, []) is None
