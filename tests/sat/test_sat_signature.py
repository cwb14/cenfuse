"""The LTR signature (cenfuse.sat.signature): families and lineages as seeds, leave-one-chromosome-out recruitment, what a
signature may decide. Synthetic GFF3 only."""
import numpy as np
import pytest

from cenfuse.sat import signature as S

L = 20_000_000
CHROMS = {f"c{i}": L for i in range(1, 9)}
CEN = {f"c{i}": 2_000_000 * i for i in range(1, 9)}        # a different place on each chromosome


def line(c, s, clade, family, age=None):
    a = f"ID=x;clade={clade}" + (f";family={family}" if family else "") + (f";K2P_T={age}" if age is not None else "")
    return f"{c}\tLTRquest\tLTR_retrotransposon\t{s + 1}\t{s + 5000}\t.\t+\t.\t{a}\n"


def clustered(c, n, clade, family, step=20_000, age=None):
    """n elements at the centromere of c."""
    return [line(c, CEN[c] + i * step, clade, family, age) for i in range(n)]


def dispersed(c, n, clade, family, rng, age=None):
    return [line(c, int(s), clade, family, age) for s in rng.integers(0, L - 6000, n)]


def write(tmp_path, rows):
    p = tmp_path / "g.gff3"
    p.write_text("##gff-version 3\n" + "".join(rows))
    return p


def genome(tmp_path, seed_clade="Tekay"):
    """A centromeric family (fam1) of a broad lineage, a smaller centromeric family (fam2, too few elements to qualify alone),
    and a large dispersed family of the same lineage, on every chromosome."""
    rng = np.random.default_rng(1)
    rows = []
    for c in CHROMS:
        rows += clustered(c, 12, seed_clade, "fam1")
        rows += clustered(c, 2, seed_clade, "fam2", step=30_000)
        rows += dispersed(c, 60, seed_clade, "fam3", rng)
    return write(tmp_path, rows)


def test_load_table_reads_families_and_ages(tmp_path):
    p = write(tmp_path, [line("c1", 1000, "CRM", "f1", 50000), line("c1", 9000, "unknown", "f1"), line("c1", 90000, "Ale", "f2"),
                         line("zz", 5, "CRM", "f1")])
    t = S.load_table(p, CHROMS)
    assert len(t.mid) == 3 and t.families == ["f1", "f2"] and sorted(t.lineages) == ["Ale", "CRM"]
    assert t.age[0] == 50000 and np.isnan(t.age[2])
    assert t.lineages[t.family_lineage[0]] == "CRM"          # the family's lineage: the most common one of its elements


def test_no_family_column_means_no_signature(tmp_path):
    p = write(tmp_path, [line("c1", 1000 + 7000 * i, "CRM", "") for i in range(40)])
    assert S.load_table(p, CHROMS) is None and S.build(p, CHROMS) is None


def test_family_is_more_specific_than_its_lineage(tmp_path):
    t = S.load_table(genome(tmp_path), CHROMS)
    by = {(z["level"], z["name"]): z for z in S.candidates(t)}
    assert by[("family", "fam1")]["specificity"] > 2 * by[("lineage", "Tekay")]["specificity"]
    assert S.candidates(t)[0]["name"] == "fam1"
    assert ("family", "fam3") not in by                       # dispersed: no core on most chromosomes


def test_seed_prefers_crm_then_cores_on_most_chromosomes():
    z = lambda n, lin, k, sp: dict(name=n, lineage=lin, cores=dict.fromkeys(range(k)), specificity=sp)
    tek, few, most = z("t", "Tekay", 8, 50.0), z("a", "CRM", 3, 40.0), z("b", "CRM", 6, 20.0)
    assert S.seed_of([tek, few, most], 8)["name"] == "b"
    assert S.seed_of([tek, few], 8)["name"] == "a"
    assert S.seed_of([tek], 8)["name"] == "t"
    assert S.seed_of([], 8) is None


def test_recruitment_is_leave_one_chromosome_out(tmp_path):
    t = S.load_table(genome(tmp_path), CHROMS)
    seed = S.seed_of(S.candidates(t), len(t.names))
    rounds = S.recruit(t, seed["cores"])
    names = {k: {t.families[f] for f in v[1]} for k, v in rounds.items()}
    assert all({"fam1", "fam2"} <= x and "fam3" not in x for x in names.values())
    # a family concentrated on one chromosome only is not recruited there: the other chromosomes do not hold it
    rows = open(tmp_path / "g.gff3").readlines()[1:] + clustered("c1", 30, "Tekay", "only_c1", step=5_000)
    t2 = S.load_table(write(tmp_path, rows), CHROMS)
    seed = S.seed_of(S.candidates(t2), len(t2.names))
    c1 = t2.names.index("c1")
    assert "only_c1" not in {t2.families[f] for f in S.recruit(t2, seed["cores"])[c1][1]}


def test_signature_cores_lie_at_the_centromeres(tmp_path):
    s = S.build(genome(tmp_path), CHROMS)
    assert s.level == "family" and s.seed == "fam1" and s.lineage == "Tekay" and not s.crm
    assert set(s.significant()) == set(CHROMS)
    for c, x in s.significant().items():
        assert x["start"] <= CEN[c] + 100_000 <= x["end"] and x["strong"]
        assert x["source"] == "signature" and "fam1" in x["families"]


def test_crm_decides_on_most_chromosomes(tmp_path):
    s = S.build(genome(tmp_path, "CRM"), CHROMS)
    assert s.crm and s.decisive
    assert "decisive" in S.decisive(s, True, "pass").why


def test_another_lineage_needs_specificity(tmp_path):
    s = S.build(genome(tmp_path), CHROMS)
    assert not s.decisive
    s.specificity = S.TE_SPECIFIC_FOLD - 1
    assert not S.decisive(s).decisive and "fold" in s.why
    s = S.build(genome(tmp_path), CHROMS)
    s.specificity = S.TE_SPECIFIC_FOLD
    assert S.decisive(s).decisive


def test_another_lineage_does_not_stand_against_a_satellite_that_shows_itself(tmp_path):
    for status in ("pass", "untested"):
        s = S.decisive(S.build(genome(tmp_path), CHROMS), True, status)
        assert not s.decisive and "one per chromosome" in s.why
    assert S.decisive(S.build(genome(tmp_path), CHROMS), True, "fail").decisive         # tomato
    assert S.decisive(S.build(genome(tmp_path), CHROMS), False, "pass").decisive        # moss: no dominant satellite


def test_core_is_strong_only_with_most_of_the_elements(tmp_path):
    rng = np.random.default_rng(2)
    rows = []
    for c in CHROMS:
        rows += clustered(c, 12, "Tekay", "fam1")
    rows += dispersed("c1", 40, "Tekay", "fam1", rng)              # c1: most of the family's elements lie outside its core
    s = S.build(write(tmp_path, rows), CHROMS)
    sig = s.significant()
    assert sig["c2"]["strong"] and not sig["c1"]["strong"]
    assert sig["c1"]["start"] <= CEN["c1"] + 100_000 <= sig["c1"]["end"]
    s.lineage = "CRM"                                               # CRM cores are always strong: the known lineage
    assert s.significant()["c1"]["strong"]


def test_seed_core_stands_where_recruitment_finds_nothing(tmp_path):
    """A lineage split into families too small to recruit (Chlorella's Osser: 54 elements): its own cores stand."""
    rows = []
    for i, c in enumerate(CHROMS):
        rows += [line(c, CEN[c] + j * 20_000, "Osser", f"f{i}_{j}") for j in range(6)]
    s = S.build(write(tmp_path, rows), CHROMS)
    assert s.level == "lineage" and s.seed == "Osser" and len(s.significant()) == len(CHROMS)
    assert {x["source"] for x in s.cores.values()} == {"seed"}


def test_table_round_trip(tmp_path):
    s = S.decisive(S.build(genome(tmp_path), CHROMS))
    S.write_table(s, tmp_path / "t.tsv")
    assert S.read_verdict(tmp_path / "t.tsv") == (True, s.why)
    assert S.read_verdict(tmp_path / "none.tsv") == (None, "")
    head, first = open(tmp_path / "t.tsv").readlines()[:2]
    r = dict(zip(head.rstrip("\n").split("\t"), first.rstrip("\n").split("\t")))
    assert r["seed"] == "fam1" and r["marker_lineage"] == "Tekay" and r["significant"] == "True"


def test_younger_elements_in_the_core_are_reported(tmp_path):
    rng = np.random.default_rng(3)
    rows = []
    for c in CHROMS:
        rows += clustered(c, 20, "CRM", "fam1", age=10_000)
        rows += dispersed(c, 15, "CRM", "fam1", rng, age=900_000)
    s = S.build(write(tmp_path, rows), CHROMS)
    assert s.age_core == 10_000 and s.age_out == 900_000 and s.young_old > 1
