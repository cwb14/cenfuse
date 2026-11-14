"""The LTR signature (cenfuse.sat.signature): families and lineages as seeds, leave-one-chromosome-out recruitment, what a
signature may decide. Synthetic GFF3 only."""
import numpy as np

from types import SimpleNamespace

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


def _two_peaks(tmp_path, young_age, old_age=1_000_000, young_n=6, old_n=10, same=False):
    """fam1 on every chromosome: an old cluster at the centromere and a younger cluster 6 Mb away."""
    rows = []
    for c in CHROMS:
        rows += clustered(c, old_n, "Tekay", "fam1", age=old_age)
        rows += [line(c, CEN[c] + 6_000_000 + i * 20_000, "Tekay", "fam1", young_age if not same else old_age) for i in range(young_n)]
        rows += dispersed(c, 30, "Tekay", "fam3", np.random.default_rng(int(c[1:])), age=500_000)
    return write(tmp_path, rows)


def test_a_comparable_younger_peak_is_chosen(tmp_path):
    """Where the signature's copies form two comparable peaks, the significantly younger one marks the active centromere
    (age_choice), and the core records where the density maximum was."""
    s = S.build(_two_peaks(tmp_path, young_age=1_000, young_n=6), CHROMS)
    moved = 0
    for c, x in s.cores.items():
        if "age_from" in x:
            moved += 1
            assert abs(x["start"] - (CEN[c] + 6_000_000)) < 200_000      # the young cluster
            assert abs(x["age_from"][0] - CEN[c]) < 200_000              # the old density maximum
    assert moved >= 6                                                     # on most chromosomes


def test_the_highest_peak_stands_without_significantly_younger_copies(tmp_path):
    s = S.build(_two_peaks(tmp_path, young_age=1_000, same=True), CHROMS)    # same ages at both peaks
    assert all("age_from" not in x for x in s.cores.values())
    s2 = S.build(_two_peaks(tmp_path, young_age=1_000, young_n=3), CHROMS)   # the younger peak under half as high
    assert all("age_from" not in x for x in s2.cores.values())


def test_a_thin_family_concentrated_across_chromosomes_seeds_a_pooled_signature(tmp_path):
    """pooled_candidates: a family with a few copies per chromosome (too few to be significant on most one by one) that lies at one
    locus on most chromosomes seeds the signature when no other candidate reaches TE_SPECIFIC_FOLD; where another candidate is
    specific enough it keeps its place."""
    rng = np.random.default_rng(2)
    rows = []
    chroms = list(CHROMS)
    for j, c in enumerate(chroms):
        n = 7 if j < 3 else 3                                 # significant cores on 3 chromosomes only
        rows += clustered(c, n, "Tekay", "famT", step=40_000)
        rows += dispersed(c, 30, "SIRE", "famS", rng)
        rows += dispersed(c, 80, "Athila", "famA", rng)
        rows += dispersed(c, 120, "Tekay", "famD", rng)       # the broad part of the lineage, as in real genomes
    p = write(tmp_path, rows)
    t = S.load_table(p, CHROMS)
    pooled = S.pooled_candidates(t)
    assert pooled and pooled[0]["name"] == "famT" and len(pooled[0]["cores"]) >= 3


def _peaks_table(tmp_path, peaks, scatter=8):
    """One chromosome's family: peaks [(position, copies, age)] and `scatter` dispersed old copies. (Table, the elements)."""
    rng = np.random.default_rng(5)
    rows = []
    for pos, n, age in peaks:
        rows += [line("c1", pos + i * 20_000, "Tekay", "fam1", age) for i in range(n)]
    rows += [line("c1", int(x), "Tekay", "fam1", 500_000) for x in rng.integers(0, L - 6000, scatter)]
    t = S.load_table(write(tmp_path, rows), CHROMS)
    return t, np.flatnonzero(t.chrom == 0)


def test_a_younger_peak_must_be_a_significant_core_to_take_over(tmp_path):
    """age_choice: a comparable younger peak that fails te.core's test stays beside the density maximum (the move would swap a
    tested core for an untested one); with enough copies to pass the same test it takes over."""
    t, sel = _peaks_table(tmp_path, [(2_000_000, 8, 1_000_000), (10_000_000, 4, 1_000)])
    second = S._peaks(t, sel, 0)[1]
    assert not S._peak_core(t, sel, 0, second)["significant"]
    assert S.age_choice(t, sel, 0) is None
    t, sel = _peaks_table(tmp_path, [(2_000_000, 10, 1_000_000), (10_000_000, 7, 1_000)])
    moved = S.age_choice(t, sel, 0)
    assert moved is not None and moved["significant"] and abs(moved["start"] - 10_000_000) < 200_000
    assert moved["k"] == 7 and moved["n"] == len(sel)


def test_the_age_test_is_corrected_for_the_peaks_compared(tmp_path, monkeypatch):
    """A P of 0.006 moves the core where one peak was compared, and not where two were (0.012 after Bonferroni)."""
    monkeypatch.setattr(S, "mannwhitneyu", lambda a, b, alternative: SimpleNamespace(pvalue=0.006))
    t, sel = _peaks_table(tmp_path, [(2_000_000, 10, 1_000_000), (10_000_000, 7, 1_000)], scatter=0)
    assert S.age_choice(t, sel, 0) is not None
    t, sel = _peaks_table(tmp_path, [(2_000_000, 10, 1_000_000), (10_000_000, 7, 1_000), (15_000_000, 7, 1_000)], scatter=0)
    assert S.age_choice(t, sel, 0) is None


def test_a_few_strong_chromosomes_do_not_make_a_pooled_family(tmp_path):
    """Fisher's P can be carried by a few chromosomes: a family clustered on 3 of 8 and dispersed on the rest has a strong combined
    P, a core on most chromosomes and 3 significant ones, and does not pool; concentrated on most chromosomes, it does."""
    rng = np.random.default_rng(3)
    few = []
    for j, c in enumerate(CHROMS):
        few += clustered(c, 12, "Tekay", "famT", step=40_000) if j < 3 else dispersed(c, 12, "Tekay", "famT", rng)
        few += dispersed(c, 100, "Tekay", "famD", rng)
    t = S.load_table(write(tmp_path, few), CHROMS)
    assert not any(z["name"] == "famT" for z in S.pooled_candidates(t))
    broad = []
    for j, c in enumerate(CHROMS):
        broad += clustered(c, 12 if j < 3 else 5, "Tekay", "famT", step=40_000)
        broad += dispersed(c, 100, "Tekay", "famD", rng)
    t = S.load_table(write(tmp_path, broad), CHROMS)
    assert any(z["name"] == "famT" for z in S.pooled_candidates(t))


def test_seed_falls_back_to_a_pooled_family():
    """A pooled family seeds only where no CRM candidate exists and the most specific candidate is under TE_SPECIFIC_FOLD."""
    weak = [dict(name="SIRE", lineage="SIRE", level="lineage", specificity=9.5, cores={i: {} for i in range(7)})]
    pooled = [dict(name="famT", lineage="Tekay", level="family", specificity=28.0, cores={0: {}, 3: {}, 4: {}})]
    s = S.seed_of(weak, 12, pooled)
    assert s["name"] == "famT" and s["pooled"]
    assert S.seed_of(weak, 12, []) is weak[0]                                   # no pooled family: as before
    strong = [dict(name="fam1", lineage="Tekay", level="family", specificity=19.0, cores={i: {} for i in range(12)})]
    s = S.seed_of(strong, 12, pooled)
    assert s["name"] == "fam1" and not s.get("pooled")                          # a specific candidate keeps its place
    crm = [dict(name="CRM", lineage="CRM", level="lineage", specificity=7.0, cores={i: {} for i in range(3)})]
    assert S.seed_of(crm, 12, pooled)["name"] == "CRM"                           # CRM precedes it


def test_a_pooled_signature_decides_without_most_significant_cores():
    sig = S.Signature(lineage="Tekay", level="family", seed="famT", specificity=28.0, seed_cores={}, n_chroms=12, pooled=True,
                      cores={f"c{i}": dict(significant=True, share=0.6) for i in range(1, 5)}, candidates=[])
    S.decisive(sig)
    assert sig.decisive
    S.decisive(S.Signature(lineage="Tekay", level="family", seed="famT", specificity=28.0, seed_cores={}, n_chroms=12, pooled=True,
                           cores={f"c{i}": dict(significant=True, share=0.6) for i in range(1, 5)}, candidates=[]), dominant=True,
               satellite_status="pass")                                         # a satellite that shows itself still stands
    weak = S.Signature(lineage="Tekay", level="family", seed="famT", specificity=8.0, seed_cores={}, n_chroms=12, pooled=True,
                       cores={f"c{i}": dict(significant=True, share=0.6) for i in range(1, 5)}, candidates=[])
    S.decisive(weak)
    assert not weak.decisive                                                    # below 10-fold it only flags


def test_a_second_comparable_peak_is_reported(tmp_path):
    """second_peak: where the signature's elements form a second peak at least half as high as the core's, it is reported beside the core (reported only)."""
    s = S.build(_two_peaks(tmp_path, young_age=1_000, same=True), CHROMS)    # same ages: the core stays at the highest peak
    alts = [(c, x["alt"]) for c, x in s.cores.items() if x.get("alt")]
    assert len(alts) >= 6
    for c, (a, b) in alts:
        assert abs(a - (CEN[c] + 6_000_000)) < 300_000


def test_a_signature_whose_cores_are_not_specific_only_flags():
    """decisive: a lineage other than CRM needs its cores to be specific on their chromosomes (median fold), not only its seed (Arachis:
    the seed family is 11.6-fold enriched, the families recruited from it make cores 1.5-fold enriched)."""
    good = S.Signature(lineage="Tekay", level="family", seed="fam1", specificity=20.0, seed_cores={}, n_chroms=4,
                       cores={f"c{i}": dict(significant=True, share=0.6, fold=20.0) for i in range(1, 5)}, candidates=[])
    S.decisive(good)
    assert good.decisive
    broad = S.Signature(lineage="Tekay", level="family", seed="fam1", specificity=11.6, seed_cores={}, n_chroms=4,
                        cores={f"c{i}": dict(significant=True, share=0.2, fold=1.5) for i in range(1, 5)}, candidates=[])
    S.decisive(broad)
    assert not broad.decisive and "1.5-fold on their chromosomes" in broad.why
