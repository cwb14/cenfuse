"""End to end on the synthetic genome of fuse_synth.py (every planted case, the planted plants), the permuted-Hi-C
negative control, a resume, and one full run with the satellite engine on cenrep's synthetic complex."""
import shutil
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent))
import fuse_synth  # noqa: E402

from cenfuse.pipeline import Run  # noqa: E402


def _calls(out, g="Syn"):
    return pd.read_csv(out / f"{g}.centromeres.tsv", sep="\t", dtype=str).set_index("chrom")


@pytest.fixture(scope="module")
def full(tmp_path_factory):
    root = tmp_path_factory.mktemp("fuse")
    truth = fuse_synth.make(root)
    out = root / "out"
    Run(truth["dir"], out, 4, satellite_dir=truth["sat"]).run()
    return truth, out


def test_planted_cases(full):
    truth, out = full
    t = _calls(out)
    for c in ("chr1", "chr2", "chr3", "chr4", "chr5", "chr6"):
        r = t.loc[c]
        assert (r.relation, r.confidence, r.support) == ("concordant", "high", "satellite,te,hic,meth"), c
    r = t.loc["chr7"]                                     # moved: methylation at the Hi-C core is the second line
    assert r.relation == "offset" and r.confidence == "high" and r["flags"] == "moved_with_meth" and "te" in r.against
    assert int(r.start) <= fuse_synth.bp(fuse_synth.F["chr7"]) <= int(r.end)
    r = t.loc["chr8"]                                     # moved: the CRM core at the Hi-C core is the second line
    assert r.relation == "offset" and r.confidence == "high" and "moved_with_te" in r["flags"]
    r = t.loc["chr9"]                                     # no second line: kept, contested, one level down
    assert (r.relation, r.confidence, r.against, r["flags"]) == ("contested", "medium", "hic", "hic_offset")
    assert int(r.start) - 1 == fuse_synth.S["chr9"] * fuse_synth.BIN
    r = t.loc["chr10"]                                    # no satellite call: placed by Hi-C, TE and methylation
    assert r.relation == "functional_only" and r.confidence == "high" and r.seq_start == "."
    assert int(r.start) <= fuse_synth.bp(fuse_synth.F["chr10"]) <= int(r.end)
    assert (t.hic_gate == "pass").all()


def test_loci_place_the_other_features(full):
    _, out = full
    lt = pd.read_csv(out / "Syn.loci.tsv", sep="\t")
    k = lambda c, pos: lt[(lt.chrom == c) & (lt.start <= pos) & (lt.end >= pos)].kind.tolist()
    assert k("chr5", 5_050_000) == ["rDNA"]
    assert k("chr3", 4_100_000) == ["interstitial_block"]
    assert k("chr1", 30_000) == ["subtelomeric_block"]
    assert k("chr7", fuse_synth.S["chr7"] * fuse_synth.BIN + 100_000) == ["sequence_centromere"]
    assert k("chr9", fuse_synth.bp(fuse_synth.F["chr9"])) == ["alt_centromere_candidate"]
    gap = lt[(lt.chrom == "chr4") & (lt.kind == "centromere")]
    assert gap.gap_bp.iloc[0] == 50_000 and gap.n_lines.iloc[0] == 4
    feats = (out / "Syn.features.gff3").read_text().splitlines()
    assert all(len(x.split("\t")) == 9 for x in feats if not x.startswith("#"))


def test_feature_classes_in_3d(full):
    _, out = full
    fc = pd.read_csv(out / "Syn.feature_contacts.tsv", sep="\t").set_index("feature_class")
    assert fc.loc["centromeres (sequence calls)", "p_enriched"] <= 0.01 and fc.loc["centromeres (sequence calls)", "fold"] > 1.5
    assert fc.loc["chromosome ends", "p_enriched"] <= 0.01           # the planted telomere pole
    assert fc.loc["satellite dispR", "p_enriched"] > 0.05             # blocks dropped at random places: no clustering
    sig = pd.read_csv(out / "Syn.feature_signatures.tsv", sep="\t").set_index("kind")
    assert sig.loc["centromere", "hic_anch_max"] > 2 * sig.loc["subtelomeric_block", "hic_anch_max"]
    assert sig.loc["subtelomeric_block", "hic_telomere_max"] > sig.loc["centromere", "hic_telomere_max"]
    assert sig.loc["centromere", "meth_pooled_extreme"] < sig.loc["subtelomeric_block", "meth_pooled_extreme"]


def test_outputs_and_resume(full):
    truth, out = full
    for f in ("Syn.centromeres.bed", "Syn.centromeres.gff3", "Syn.provenance.tsv", "Syn.plants.tsv",
              "Syn.plant_sites.tsv", "Syn.hic_anchored.bedgraph", "cenfuse_summary.tsv", "hic/Syn.centromeres.tsv",
              "hic/Syn.P3.centromeres.tsv", "methylation/Syn.candidates.tsv", "methylation/Syn.loco_signatures.tsv"):
        assert (out / f).exists(), f
    s = pd.read_csv(out / "cenfuse_summary.tsv", sep="\t").iloc[0]
    assert s.meth_signature == "CHG:LTR:CRM:down" and s.hic_gate_pass == 10 and s.offset == 2
    assert s.meth_deciding == "P1" and s.hic_library == "hic" and s.plants_shifted == "P3:chr3"
    before = {f: (out / f).read_text() for f in ("Syn.centromeres.tsv", "Syn.plants.tsv", "Syn.provenance.tsv")}
    Run(truth["dir"], out, 4, satellite_dir=truth["sat"]).run()          # resumes from every checkpoint
    assert all((out / f).read_text() == v for f, v in before.items())
    sites = pd.read_csv(out / "Syn.plant_sites.tsv", sep="\t")
    assert set(sites["sample"]) == set(fuse_synth.SAMPLES) and sites.groupby("site")["sample"].nunique().eq(3).all()


def test_provenance_of_every_input(full):
    """S1: P1 is the reference plant, P2 and P3 are not (P2 also by its missing segment), P3's Hi-C library is not
    by name, the reference library is of unknown origin and decides; only P1's methylation decides."""
    _, out = full
    pv = pd.read_csv(out / "Syn.provenance.tsv", sep="\t", na_values=".", keep_default_na=False)
    pv = pv.set_index(["kind", "input"])
    assert pv.loc[("meth", "P1"), "verdict"] == "reference" and bool(pv.loc[("meth", "P1"), "decides"])
    for p in ("P2", "P3"):
        assert pv.loc[("meth", p), "verdict"] == "other" and not bool(pv.loc[("meth", p), "decides"])
        assert pv.loc[("meth", p), "fixed_per_kb"] > 20 * pv.loc[("meth", "P1"), "fixed_per_kb"]
    seg = pv.loc[("meth", "P2"), "missing_segments"]
    s, e = (int(x) for x in seg.split(":")[1].split("-"))
    assert seg.startswith(f"{fuse_synth.MISSING[1]}:") and abs(s - fuse_synth.MISSING[2]) <= 1_000_001 \
        and abs(e - fuse_synth.MISSING[3]) <= 1_000_000
    assert pd.isna(pv.loc[("meth", "P1"), "missing_segments"]) and (pv.on_motif.dropna() > 0.99).all()
    assert pv.loc[("hic", "P3"), "verdict"] == "other" and "by name" in pv.loc[("hic", "P3"), "reason"]
    assert pv.loc[("hic", "hic"), "verdict"] == "unknown" and bool(pv.loc[("hic", "hic"), "decides"])
    t = _calls(out)
    assert (t.meth_samples == "P1").all() and (t.hic_provenance == "unknown").all()


def test_plants_are_scored_one_by_one(full):
    """S3: P3's chr3 centromere is shifted to its planted place, by methylation and by its own Hi-C library; no other
    plant x chromosome is shifted."""
    _, out = full
    pl = pd.read_csv(out / "Syn.plants.tsv", sep="\t")
    plant, chrom, cbin = fuse_synth.SHIFT
    r = pl[(pl.plant == plant) & (pl.chrom == chrom)].iloc[0]
    target = fuse_synth.bp(cbin)
    assert r.verdict == "shifted" and r.meth_verdict == "shifted" and r.hic_verdict == "shifted"
    assert r.alt_start <= target <= r.alt_end and r.hic_peak_start <= target + fuse_synth.BIN
    assert "satellite:satA" in r.alt_origins and (r.loss_q <= 0.05 or r.gain_q <= 0.05)
    rest = pl[~((pl.plant == plant) & (pl.chrom == chrom))]
    assert not rest.verdict.isin(["shifted", "conflict", "call_absent"]).any(), rest[rest.verdict == "shifted"]
    assert set(pl[pl.chrom == chrom].set_index("plant").verdict.drop(plant)) == {"shared"}
    t = _calls(out)
    assert t.loc[chrom, "plants_shifted"] == "1" and t.loc["chr1", "plants_shared"] == "3"


def test_permuted_hic_changes_nothing(tmp_path):
    truth = fuse_synth.make(tmp_path, permute=True, meth=False)
    out = tmp_path / "out"
    Run(truth["dir"], out, 4, satellite_dir=truth["sat"]).run()
    t = _calls(out)
    assert (t.hic_gate == "fail").all()
    called = t[t.seq_start != "."]
    assert (called.relation == "sequence_only").all() and (called.start == called.seq_start).all()
    assert t.loc["chr10", "relation"] == "none"


TOOLS = ("FasTAN", "FAtoGDB", "ONEview", "makeblastdb", "blastn", "mafft", "cd-hit-est", "minimap2", "cmsearch")


@pytest.mark.skipif(not all(shutil.which(x) for x in TOOLS), reason="satellite-engine tools missing")
def test_full_pipeline_with_the_satellite_engine(tmp_path):
    sys.path.insert(0, str(Path(__file__).parents[1] / "sat"))
    from sat_common import TEST_PARAMS
    from sat_synth import build_complex
    cdir = tmp_path / "Synth"
    build_complex(cdir)
    fa = cdir / "Pa_A.fa"                                  # unplaced scaffolds: left out by the chromosome set
    with open(fa, "a") as f:
        f.write(">scaffold_9\n" + "ACGT" * 500 + "\n")
    import pysam
    pysam.faidx(str(fa))
    (cdir / "hic").mkdir()                                  # chromosomes far too short for cenhic: the line is skipped
    import gzip
    with gzip.open(cdir / "hic" / "Pa_A.pairs.gz", "wt") as f:
        f.write("Pa_A01\t100\t+\tPa_A02\t200\t-\n" * 50)
    Run(cdir, tmp_path / "out", 4, sat_params=TEST_PARAMS).run()
    out = tmp_path / "out"
    t = _calls(out, "Pa_A")
    assert "scaffold_9" not in t.index and len(t) == 4
    summ = pd.read_csv(out / "satellite" / "Pa_A_centromere_summary.tsv", sep="\t", dtype=str).set_index("chrom")
    assert list(summ.index) == list(t.index)
    called = t[t.seq_start != "."]
    assert len(called) == 3 and (called.relation == "sequence_only").all()
    assert (called.seq_start == summ.loc[called.index, "cen_start"]).all()
    assert (out / "Pa_A.loci.tsv").exists() and (out / "cenfuse_summary.tsv").exists()
    assert (t.hic_gate == ".").all()                        # Hi-C skipped with a warning; the run completed
