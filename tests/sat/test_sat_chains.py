"""Families that are chains of partial hits, not tandem satellites, and what else the satellite engine reports about a
call: the symmetric period test, the families the CRM line rejected, element counts, the array's own spacing, the
clusters left out by the family cap."""
import numpy as np
import pandas as pd

from cenfuse.sat import centromere as cen
from cenfuse.sat import classify as cl
from cenfuse.sat import families, output, te
from cenfuse.sat.annotate import ARR_COLS
from cenfuse.sat.params import Params

P = Params()


def arr(fam, chrom, start, n, mlen=180):
    return (0, chrom, start, start + n * mlen, "+", fam, n, 0.9, mlen)


# ---- #17: the period test

def test_period_test_is_symmetric():
    """155 vs 180 bp passed one way and failed the other (155/180 - 1 = -0.14; 180/155 - 1 = 0.16), so maize's CentC
    chromosomes counted for the knob repeat, and not the reverse."""
    info = {"cent": dict(superfamily="c", length=155), "knob": dict(superfamily="k", length=180),
            "near": dict(superfamily="n", length=165)}
    assert not cl.related("cent", "knob", info) and not cl.related("knob", "cent", info)
    assert cl.related("cent", "near", info) and cl.related("near", "cent", info)


def test_same_superfamily_is_related_at_any_period():
    info = {"a": dict(superfamily="s", length=100), "b": dict(superfamily="s", length=400)}
    assert cl.related("a", "b", info) and cl.related("b", "a", info)


def test_a_knob_family_no_longer_borrows_centc_chromosomes_for_its_dominance():
    """The knob-like family is strong and on most chromosomes, but the largest only on 2 of 6: CentC's 4 no longer count
    for it through the one-way period test."""
    chroms = {f"c{i}": 50_000_000 for i in range(6)}
    rows = []
    for i in range(6):
        rows.append(arr("cent", f"c{i}", 20_000_000, 1200 if i < 4 else 300, 155))
        rows.append(arr("knob", f"c{i}", 40_000_000, 800, 180))
    a = pd.DataFrame(rows, columns=ARR_COLS)
    info = {"cent": dict(superfamily="c", length=155, rdna=".", chain_group="cent"),
            "knob": dict(superfamily="k", length=180, rdna=".", chain_group="knob")}
    fams, mains = cl.family_metrics(a, chroms, None, P, info)
    fams, _ = cl.choose_centromeric(fams, mains, chroms, P, info, a)
    f = fams.set_index("family")
    assert bool(f.loc["cent", "dominant"]) and not bool(f.loc["knob", "dominant"])


# ---- what a chain is

def test_a_chain_has_neither_tandem_arrays_nor_whole_units():
    row = lambda tf, wf: pd.Series(dict(tandem_frac=tf, whole_frac=wf, locus_fold=9.0))
    assert cl.partial_chain(row(0.1, 0.01))
    assert not cl.partial_chain(row(0.6, 0.01))              # in tandem at some period: Cuscuta 1.9 kb
    assert not cl.partial_chain(row(0.1, 0.9))               # whole units, not found by the tandem finder
    assert not cl.partial_chain(row(float("nan"), float("nan")))           # not judged


def test_a_dispersed_repeat_is_a_chain_spread_along_the_chromosome():
    a, b = pd.Series(dict(tandem_frac=0.1, whole_frac=0.1, locus_fold=1.2)), pd.Series(dict(tandem_frac=0.1, whole_frac=0.1, locus_fold=9.0))
    assert cl.dispersed_repeat(a) and cl.partial_chain(a)
    assert not cl.dispersed_repeat(b) and cl.partial_chain(b)


def chain_genome():
    """Every chromosome: a 1.8-kb chain family with the largest cluster, and a real satellite in a smaller one."""
    chroms = {f"c{i}": 50_000_000 for i in range(6)}
    rows = []
    for i in range(6):
        rows.append(arr("chain", f"c{i}", 20_000_000, 200, 1800))
        rows.append(arr("sat", f"c{i}", 30_000_000, 800, 170))
    a = pd.DataFrame(rows, columns=ARR_COLS)
    info = {"chain": dict(superfamily="chain", length=1800, rdna=".", chain_group="chain"),
            "sat": dict(superfamily="sat", length=170, rdna=".", chain_group="sat")}
    return a, chroms, info


def test_dominance_for_the_honesty_cap_is_judged_without_chains():
    """The satellite leads on every chromosome once the chain is left out, and only then: dominant_satellite says so,
    while `dominant` stays as it was, so the choice per chromosome does not change."""
    a, chroms, info = chain_genome()
    fams, mains = cl.family_metrics(a, chroms, None, P, info, tandem={"chain": (0.05, 0.02), "sat": (0.95, 0.95)})
    fams, per_chrom = cl.choose_centromeric(fams, mains, chroms, P, info, a)
    f = fams.set_index("family")
    assert bool(f.loc["chain", "dominant"]) and not bool(f.loc["chain", "dominant_satellite"])
    assert bool(f.loc["sat", "dominant_satellite"])
    assert all(v[0] == "chain" for v in per_chrom.values())


def test_a_genome_where_only_a_chain_dominates_has_no_dominant_satellite():
    chroms = {f"c{i}": 50_000_000 for i in range(6)}
    a = pd.DataFrame([arr("chain", f"c{i}", 20_000_000, 200, 1800) for i in range(6)], columns=ARR_COLS)
    info = {"chain": dict(superfamily="chain", length=1800, rdna=".", chain_group="chain")}
    fams, mains = cl.family_metrics(a, chroms, None, P, info, tandem={"chain": (0.05, 0.02)})
    fams, _ = cl.choose_centromeric(fams, mains, chroms, P, info, a)
    assert bool(fams.dominant.iat[0]) and not bool(fams.dominant_satellite.iat[0])


def test_without_chains_both_dominance_columns_agree():
    a, chroms, info = chain_genome()
    fams, mains = cl.family_metrics(a, chroms, None, P, info, tandem={"chain": (0.95, 0.95), "sat": (0.95, 0.95)})
    fams, _ = cl.choose_centromeric(fams, mains, chroms, P, info, a)
    assert fams.dominant.astype(bool).tolist() == fams.dominant_satellite.astype(bool).tolist()


def test_ltr_overlap_credits_every_element_overlapping_an_array():
    """An array inside two small SIRE elements and beside a larger Tekay overlap: the clade with most bp, not the clade
    of the single largest overlap."""
    arrays = pd.DataFrame([arr("f", "c", 1_000, 10, 100)], columns=ARR_COLS)              # 1000-2000
    iv = {"c": (np.array([1_000, 1_300, 1_650]), np.array([1_250, 1_600, 2_000]), ["SIRE", "SIRE", "Tekay"])}
    share, clade = te.ltr_overlap(arrays, iv, 1_000)["f"]
    assert clade == "SIRE" and round(share, 2) == 0.90


# ---- the cap

def row(fam="chain", conf="high", evidence="satellite", core=None, clade="."):
    r = {k: "." for k in cen.REG_COLS}
    r.update(chrom="c", chrom_len=50_000_000, cen_start=20_000_000, cen_end=22_000_000, evidence=evidence, confidence=conf, flags=".",
             main_family=fam, sat_bp=1_000_000)
    if core:
        r.update(te_clade=clade, te_core_start=core[0], te_core_end=core[1])
    return r


ARRAYS = pd.DataFrame([arr("chain", "c", 20_000_000, 500, 1800)], columns=ARR_COLS)


def cap(r, arrays=ARRAYS, primary=("chain",), chains=("chain",), derived=None, decisive=False):
    rows = {"c": r}
    cen.cap_chain_calls(rows, arrays, set(primary), set(chains), derived or {}, decisive, P)
    return rows["c"]


def test_a_chain_call_with_no_core_is_low():
    r = cap(row())
    assert r["confidence"] == "low" and r["flags"] == "partial_chain"


def test_a_chain_call_on_a_decisive_core_is_at_most_medium():
    r = cap(row(core=(21_000_000, 21_500_000), clade="CRM"), decisive=True)
    assert r["confidence"] == "medium" and "partial_chain" in r["flags"]


def test_a_core_merely_near_a_chain_call_does_not_lift_it():
    """`agrees` allows a gap as wide as the call, so a 10-Mb call agrees with most places on a chromosome: only overlap
    counts here."""
    r = cap(row(core=(23_000_000, 23_500_000), clade="CRM"), decisive=True)
    assert r["confidence"] == "low"


def test_a_flags_only_core_counts_only_if_the_family_is_not_made_of_its_clade():
    near = dict(core=(21_000_000, 21_500_000), clade="Athila")
    assert cap(row(**near))["confidence"] == "medium"
    assert cap(row(**near), derived={"chain": "Athila"})["confidence"] == "low"           # one piece of evidence
    assert cap(row(**near), derived={"chain": "SIRE"})["confidence"] == "medium"


def test_a_chain_call_is_never_raised():
    assert cap(row(conf="low", core=(21_000_000, 21_500_000), clade="CRM"), decisive=True)["confidence"] == "low"


def test_a_decisive_core_inside_the_call_narrows_it():
    r = cap(row(core=(21_000_000, 21_500_000), clade="CRM"), decisive=True)
    assert (r["cen_start"], r["cen_end"]) == (21_000_000, 21_500_000) and "chain_to_core" in r["flags"]


def test_a_core_that_may_not_act_does_not_narrow():
    r = row(core=(21_000_000, 21_500_000), clade="Tekay")
    rows = {"c": r}
    cen.cap_chain_calls(rows, ARRAYS, {"chain"}, {"chain"}, {}, True, P, cores={"c": dict(strong=False)})
    assert (r["cen_start"], r["cen_end"]) == (20_000_000, 22_000_000) and "chain_to_core" not in r["flags"]


def test_a_narrower_core_only_touching_the_call_does_not_narrow():
    r = cap(row(core=(21_900_000, 22_400_000), clade="CRM"), decisive=True)
    assert (r["cen_start"], r["cen_end"]) == (20_000_000, 22_000_000) and "chain_to_core" not in r["flags"]


def test_a_narrowed_call_is_measured_again(tmp_path):
    import pysam
    fa = tmp_path / "g.fa"
    seq = ["A"] * 300_000
    seq[10_000:10_050] = ["N"] * 50                                      # a gap outside the core
    seq[120_000:120_020] = ["N"] * 20                                    # a gap inside it
    fa.write_text(">c\n" + "".join(seq) + "\n")
    pysam.faidx(str(fa))
    arrays = pd.DataFrame([(0, "c", 5_000, 6_000, "+", "chain", 5, 0.9, 200), (0, "c", 110_000, 112_000, "+", "chain", 10, 0.9, 200)], columns=ARR_COLS)
    r = {k: "." for k in cen.REG_COLS}
    r.update(chrom="c", chrom_len=300_000, cen_start=5_000, cen_end=200_000, evidence="satellite", confidence="high", flags="near_end,gaps:2,gene_rich",
             main_family="chain", sat_bp=3_000, main_family_bp=3_000, te_clade="CRM", te_core_start=100_000, te_core_end=150_000, te_core_elements=9,
             n_gaps=2, gap_bp=70, n_genes=3, gene_ratio=2.0)
    genes = {"c": np.array([10_000, 120_000, 140_000, 250_000])}
    with pysam.FastaFile(str(fa)) as f:
        cen.cap_chain_calls({"c": r}, arrays, {"chain"}, {"chain"}, {}, True, P, fa=f, genes=genes)
    assert (r["cen_start"], r["cen_end"]) == (100_000, 150_000) and "chain_to_core" in r["flags"]
    assert (r["n_gaps"], r["gap_bp"]) == (1, 20) and r["sat_bp"] == 2_000 and r["main_family_bp"] == 2_000
    assert r["n_genes"] == 2 and "gaps:1" in r["flags"] and "gaps:2" not in r["flags"] and r["te_call_elements"] == 9
    assert "gene_rich" in r["flags"] and r["gene_ratio"] == 3.0                                       # 2 genes in 50 kb, 4 in 300 kb


def test_a_chain_call_otherwise_never_moves():
    for r in (cap(row(core=(21_000_000, 21_500_000), clade="Athila")),                      # not decisive
              cap(row(core=(23_000_000, 23_500_000), clade="CRM"), decisive=True),          # no overlap
              cap(row(core=(19_000_000, 23_000_000), clade="CRM"), decisive=True)):         # the core is wider
        assert (r["cen_start"], r["cen_end"]) == (20_000_000, 22_000_000) and "chain_to_core" not in r["flags"]


def test_a_real_satellite_layer_in_the_call_spares_it():
    a = pd.DataFrame([arr("chain", "c", 20_000_000, 500, 1800), arr("sat", "c", 20_100_000, 200, 170)], columns=ARR_COLS)  # 34 kb
    r = cap(row(), arrays=a, primary=("chain", "sat"))
    assert r["confidence"] == "high" and r["flags"] == "."


def test_a_thin_layer_does_not_spare_it():
    a = pd.DataFrame([arr("chain", "c", 20_000_000, 500, 1800), arr("sat", "c", 20_100_000, 50, 170)], columns=ARR_COLS)   # 8.5 kb
    assert cap(row(), arrays=a, primary=("chain", "sat"))["confidence"] == "low"


def test_a_real_satellite_outside_the_call_does_not_spare_it():
    a = pd.DataFrame([arr("chain", "c", 20_000_000, 500, 1800), arr("sat", "c", 30_000_000, 800, 170)], columns=ARR_COLS)
    assert cap(row(), arrays=a, primary=("chain", "sat"))["confidence"] == "low"


def test_other_calls_are_left_alone():
    assert cap(row(fam="sat"))["confidence"] == "high"
    for ev in ("retrotransposon", "long_unit_tandem"):
        r = cap(row(evidence=ev))
        assert r["confidence"] == "high" and r["flags"] == "."
    assert cap(dict(row(), evidence="none", confidence="none", cen_start=".", cen_end="."))["confidence"] == "none"


def test_a_weak_chain_call_is_flagged_too():
    r = cap(row(evidence="weak_satellite", conf="medium"))
    assert r["confidence"] == "low" and r["flags"] == "partial_chain"


def test_flags_are_added_to_those_of_the_row():
    r = row(core=(21_000_000, 21_500_000), clade="CRM")
    r["flags"] = "te_concordant"
    assert cap(r, decisive=True)["flags"] == "te_concordant,partial_chain,chain_to_core"


def test_chain_families_lists_only_chains():
    t = pd.DataFrame(dict(family=["a", "b", "c"], tandem_frac=[0.1, 0.9, float("nan")], whole_frac=[0.1, 0.1, float("nan")], locus_fold=[9.0, 9.0, 9.0]))
    assert cl.chain_families(t) == {"a"}
    assert cl.chain_families(pd.DataFrame(dict(family=["a"]))) == set()                  # a table of an older version


def test_derived_families_are_those_mostly_inside_elements():
    t = pd.DataFrame(dict(family=["a", "b"], in_ltr=[0.65, 0.2], ltr_clade=["SIRE", "CRM"]))
    assert cen.derived_families(t) == {"a": "SIRE"}
    assert cen.derived_families(pd.DataFrame(dict(family=["a"]))) == {}


# ---- families the CRM line rejected, element counts

def te_arrays():
    """Family X sits 10 Mb from the CRM cores on every chromosome; Y on them."""
    chroms = {f"c{i}": 50_000_000 for i in range(4)}
    rows = []
    for i in range(4):
        rows.append(arr("X", f"c{i}", 10_000_000, 1000, 170))
        rows.append(arr("Y", f"c{i}", 30_000_000, 300, 170))
    a = pd.DataFrame(rows, columns=ARR_COLS)
    info = {"X": dict(superfamily="X", length=170, rdna=".", chain_group="X"), "Y": dict(superfamily="Y", length=170, rdna=".", chain_group="Y")}
    cores = {c: dict(start=29_900_000, end=30_200_000, k=20, n=30) for c in chroms}
    return a, chroms, info, cores


def test_a_family_rejected_by_the_crm_line_is_marked_in_its_table():
    a, chroms, info, cores = te_arrays()
    fams, mains = cl.family_metrics(a, chroms, None, P, info)
    out, _, rejected = cl.choose_with_te(fams, mains, chroms, P, info, ("CRM", cores, True), a)
    f = out.set_index("family")
    assert "X" in rejected and bool(f.loc["X", "te_rejected"]) and not bool(f.loc["Y", "te_rejected"])


def test_a_flags_only_marker_rejects_nothing():
    a, chroms, info, cores = te_arrays()
    fams, mains = cl.family_metrics(a, chroms, None, P, info)
    out, _, rejected = cl.choose_with_te(fams, mains, chroms, P, info, ("Athila", cores, False), a)
    assert not rejected and not out.te_rejected.astype(bool).any()


def test_rejected_cluster_names_the_array_a_family_would_have_given():
    a, _, _, _ = te_arrays()
    x = cen.rejected_cluster(a[a.chrom == "c0"], {"X"}, 30_000_000, 30_051_000, P)
    assert x == ("X", 10_000_000, 10_170_000)


def test_rejected_cluster_inside_the_call_or_small_is_not_named():
    a, _, _, _ = te_arrays()
    c0 = a[a.chrom == "c0"]
    assert cen.rejected_cluster(c0, {"X"}, 9_000_000, 11_000_000, P) is None                # the call holds it
    assert cen.rejected_cluster(c0, {"none"}, 30_000_000, 30_051_000, P) is None
    small = pd.DataFrame([arr("X", "c0", 10_000_000, 10, 170)], columns=ARR_COLS)
    assert cen.rejected_cluster(small, {"X"}, 30_000_000, 30_051_000, P) is None


def test_elements_are_counted_in_the_call_as_the_head_to_head_test_counts_them():
    mids = np.array([100, 200, 300, 5_000, 5_100], dtype=np.int64)
    c = dict(start=4_000, end=6_000, k=2, mids=mids)
    assert te.elements_in(50, 350, c) == 3
    assert te.elements_in(0, 50, c) == 0
    assert te.elements_in(0, 50, dict(c, mids=None)) is None


# ---- the array's own spacing

def test_spacing_is_the_median_distance_between_neighbors_in_their_array():
    """Rice: the 165-bp CentO profile fits 155-bp copies with a 10-bp gap; the arrays are spaced 155 bp."""
    m = pd.DataFrame(dict(array_id=[1] * 6 + [2] * 4, family="f", start=[0, 155, 310, 465, 620, 775, 10_000, 10_155, 10_310, 10_465]))
    assert output.copy_spacing(m, "f") == 155


def test_spacing_needs_copies_in_a_row():
    assert output.copy_spacing(pd.DataFrame(dict(array_id=[1, 1], family="f", start=[0, 155])), "f") is None
    assert output.copy_spacing(pd.DataFrame(dict(array_id=[1, 2, 3], family="f", start=[0, 155, 310])), "f") is None
    assert output.copy_spacing(pd.DataFrame(dict(array_id=[], family=[], start=[])), "f") is None


# ---- the clusters left out by the cap

def test_dropped_clusters_are_listed():
    seeds = pd.DataFrame(dict(rep=["r1", "r1", "r2"], genome="g", chrom=["c1", "c2", "c1"], start=[0, 100, 500], end=[30_000, 20_100, 21_000],
                              period=[180, 180, 1_900], copies=[160.0, 110.0, 10.0], consensus=["A" * 180, "A" * 180, "C" * 1_900]))
    t = families.dropped_table(seeds, [(50_000, "r1"), (20_500, "r2")])
    assert t.columns.tolist() == families.DROPPED_COLS
    assert t.rep_seed.tolist() == ["r1", "r2"] and t.max_genome_span.tolist() == [50_000, 20_500]
    assert t.period.tolist() == [180, 1_900] and t.n_seeds.tolist() == [2, 1] and t.n_chroms.tolist() == [2, 1]
    assert t.largest_seed.iat[0] == "c1:1-30000" and t.consensus.iat[1] == "C" * 1_900


def test_no_dropped_cluster_gives_an_empty_table_with_its_columns():
    t = families.dropped_table(pd.DataFrame(columns=["rep", "genome", "chrom", "start", "end", "period", "copies", "consensus"]), [])
    assert t.empty and t.columns.tolist() == families.DROPPED_COLS


# ---- the stage: what it reads from the classify table and writes to the regions table

def stage_case(tmp_path, monkeypatch, tandem, tev=None):
    """A real classify table of chain_genome (one chain, one satellite on every chromosome) and a satellite engine whose
    centromere stage calls the chain's cluster on every chromosome. Returns (the pipeline, what call_regions was given)."""
    from cenfuse.sat import pipeline as pl
    a, chroms, info = chain_genome()
    d = tmp_path / "g"
    d.mkdir()
    (d / "g.fa").write_text("".join(f">{c}\nACGT\n" for c in chroms))
    (d / "g.fa.fai").write_text("".join(f"{c}\t{n}\t{i}\t4\t5\n" for i, (c, n) in enumerate(chroms.items())))
    p = pl.Pipeline(d, 1, tmp_path / "cache", Params(), out_dir=tmp_path / "out")
    for st in ("annotate", "classify", "centromeres", "families"):
        (p.work / st).mkdir(parents=True)
    a.to_csv(p.work / "annotate" / "g.arrays.tsv", sep="\t", index=False)
    fams, mains = cl.family_metrics(a, chroms, None, P, info, tandem=tandem)
    fams, per_chrom = cl.choose_centromeric(fams, mains, chroms, P, info, a)
    fams.to_csv(p.work / "classify" / "g.families.tsv", sep="\t", index=False)
    pd.DataFrame(dict(fam_id=["chain", "sat"], superfamily=["chain", "sat"], consensus_len=[1800, 170], rdna=[".", "."],
                      chain_group=["chain", "sat"])).to_csv(p.work / "families" / "families.tsv", sep="\t", index=False)
    seen = {}

    def calls(genome, arrays, primary, genes, seed_rows, tmpdir, threads, params, **kw):
        seen.update(kw)
        rows = {}
        for c, n in chroms.items():
            r = {k: "." for k in cen.REG_COLS}
            r.update(chrom=c, chrom_len=n, cen_start=20_000_000, cen_end=20_360_000, evidence="satellite", confidence="high",
                     flags=".", main_family="chain", sat_bp=360_000, n_gaps=0, gap_bp=0)
            rows[c] = r
        from cenfuse.sat.longunit import LU_COLS
        return rows, pd.DataFrame(columns=LU_COLS)

    monkeypatch.setattr(cen, "call_regions", calls)
    monkeypatch.setattr(p, "seed_table", lambda g: pd.DataFrame(columns=["chrom", "start", "end"]))
    monkeypatch.setattr(p, "te_evidence", lambda g: tev)
    return p, seen


def test_the_stage_flags_chain_calls_and_caps_genomes_with_no_real_dominant(tmp_path, monkeypatch):
    p, seen = stage_case(tmp_path, monkeypatch, {"chain": (0.05, 0.02), "sat": (0.95, 0.95)})
    p.stage_centromeres()
    r = cen.load_regions(p.work / "centromeres" / "g.regions.tsv")
    assert all(x["confidence"] == "low" and "partial_chain" in x["flags"] for x in r.values())
    assert seen["capped"] is False                                       # the satellite leads once the chain is left out


def test_the_stage_leaves_a_real_satellite_family_alone(tmp_path, monkeypatch):
    p, seen = stage_case(tmp_path, monkeypatch, {"chain": (0.95, 0.95), "sat": (0.95, 0.95)})
    p.stage_centromeres()
    r = cen.load_regions(p.work / "centromeres" / "g.regions.tsv")
    assert all(x["confidence"] == "high" and "partial_chain" not in x["flags"] for x in r.values())


def test_the_stage_reads_an_older_classify_table(tmp_path, monkeypatch):
    """Without the new columns: `dominant` decides the cap, as before, and nothing is flagged."""
    p, seen = stage_case(tmp_path, monkeypatch, {"chain": (0.05, 0.02), "sat": (0.95, 0.95)})
    f = p.work / "classify" / "g.families.tsv"
    pd.read_csv(f, sep="\t").drop(columns=["dominant_satellite", "te_rejected", "tandem_frac", "whole_frac", "locus_fold"]).to_csv(f, sep="\t", index=False)
    p.stage_centromeres()
    r = cen.load_regions(p.work / "centromeres" / "g.regions.tsv")
    assert all(x["confidence"] == "high" for x in r.values()) and seen["capped"] is False
