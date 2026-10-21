import gzip
import shutil

import pandas as pd
import pysam
import pytest
import meth_synth as synth

from cenfuse.meth import bins, cli


def _args(syn, out, *extra):
    return ["run", "--samples", str(syn["samples"]), "--fai", str(syn["fai"]), "--ltrquest", str(syn["ltrquest"]),
            "--satellites", str(syn["satellites"]), "--genes", str(syn["genes"]), "-o", str(out), "-t", "2", *extra]


def test_end_to_end_with_anchors_and_resume(syn, prep, tmp_path, capsys):
    out = tmp_path / "r" / "x"
    assert cli.main(_args(syn, out, "--anchors", str(syn["anchors"]), "--exclude", str(syn["exclude"]))) == 0
    for suf in ("candidates.bed", "candidates.tsv", "anchors.tsv", "signature.tsv", "qc.tsv", "overview.pdf",
                "S1.adjusted.bedgraph", "S2.adjusted.bedgraph", "S3.adjusted.bedgraph"):
        assert (tmp_path / "r" / f"x.{suf}").stat().st_size > 0, suf
    head = (tmp_path / "r" / "x.candidates.bed").read_text().splitlines()[0]
    assert head.startswith("#cenfuse-meth") and "anchors=hic" in head and "CHG:LTR:CRM:down" in head
    cand = pd.read_csv(tmp_path / "r" / "x.candidates.tsv", sep="\t")
    assert (cand.label == "centromere_like").sum() >= 15
    assert set(cand.anchor_source) == {"hic"}
    assert all(len(ln.split("\t")) == 9 for ln in (tmp_path / "r" / "x.candidates.bed").read_text().splitlines()[1:])
    key = bins.fingerprint([syn["ltrquest"], syn["satellites"], syn["genes"]], "|".join(prep["cmap"].names))
    assert [p.name for p in (tmp_path / "r" / "meth_cache").iterdir()] == [key]     # keyed by the class list too
    capsys.readouterr()
    other = tmp_path / "r" / "y"                         # another prefix in the same folder shares the cache
    assert cli.main(_args(syn, other, "--anchors", str(syn["anchors"]), "--exclude", str(syn["exclude"]))) == 0
    assert "(24 from cache)" in capsys.readouterr().err


def test_without_anchors(syn, tmp_path, capsys):
    out = tmp_path / "n" / "x"
    assert cli.main(_args(syn, out)) == 0
    assert "uncalibrated" in capsys.readouterr().err
    assert not (tmp_path / "n" / "x.anchors.tsv").exists() and not (tmp_path / "n" / "x.signature.tsv").exists()
    cand = pd.read_csv(tmp_path / "n" / "x.candidates.tsv", sep="\t")
    assert set(cand.label) == {"dip"} and set(cand.anchor_source) == {"none"}


def test_edta_input(syn, tmp_path, capsys):
    g = tmp_path / "x.TEanno.gff3"
    with open(syn["ltrquest"]) as f, open(g, "w") as o:
        for ln in f:
            p = ln.split("\t")
            if len(p) == 9 and p[2] == "LTR_retrotransposon":
                p[2], p[8] = "Gypsy_LTR_retrotransposon", "ID=x;Classification=LTR/Gypsy\n"
                o.write("\t".join(p))
    args = _args(syn, tmp_path / "e" / "x")
    i = args.index("--ltrquest")
    args[i:i + 2] = ["--edta", str(g)]
    assert cli.main(args) == 0
    assert "coarser" in capsys.readouterr().err


def test_missing_tabix_index_fails(syn, tmp_path, capsys):
    d = tmp_path / "s"
    d.mkdir()
    shutil.copy(syn["dir"] / "S1.CHG.bedmethyl.gz", d / "S1.CHG.bedmethyl.gz")
    (d / "samples.tsv").write_text("sample\tCHG\nS1\tS1.CHG.bedmethyl.gz\n")
    args = _args(syn, tmp_path / "o" / "x")
    args[args.index("--samples") + 1] = str(d / "samples.tsv")
    assert cli.main(args) == 2
    assert "tabix index" in capsys.readouterr().err


def test_renamed_chromosomes_fail_fast(syn, tmp_path, capsys):
    fai = tmp_path / "g.fai"
    fai.write_text("".join(ln.replace("chr", "") + "\n" for ln in syn["fai"].read_text().splitlines()))
    args = _args(syn, tmp_path / "o" / "x")
    args[args.index("--fai") + 1] = str(fai)
    assert cli.main(args) == 2
    err = capsys.readouterr().err
    assert "same assembly?" in err and "chr1" in err


def test_anchor_on_unknown_sequence_fails(syn, tmp_path, capsys):
    a = tmp_path / "a.bed"
    a.write_text("chrUn\t0\t1000\n")
    assert cli.main(_args(syn, tmp_path / "o" / "x", "--anchors", str(a))) == 2
    assert "absent from the .fai" in capsys.readouterr().err


def test_sample_missing_a_chromosome_and_low_depth_sample(syn, tmp_path, capsys):
    d = tmp_path / "s"
    d.mkdir()
    src = syn["dir"] / "S2.CHG.bedmethyl.gz"
    with gzip.open(src, "rt") as f, open(d / "S2.bed", "w") as o:          # S2 without chr6
        o.writelines(ln for ln in f if not ln.startswith("chr6\t"))
    pysam.tabix_index(str(d / "S2.bed"), preset="bed", force=True)
    with gzip.open(syn["dir"] / "S3.CHG.bedmethyl.gz", "rt") as f, open(d / "L.bed", "w") as o:   # ~3x sample
        for i, ln in enumerate(f):
            p = ln.split("\t")
            nv = 3 if i % 3 else 0
            if nv:
                p[9], p[11], p[12] = "3", str(min(int(p[11]), 3)), str(3 - min(int(p[11]), 3))
                o.write("\t".join(p))
    pysam.tabix_index(str(d / "L.bed"), preset="bed", force=True)
    with gzip.open(syn["dir"] / "S3.CHG.bedmethyl.gz", "rt") as f, open(d / "Z.bed", "w") as o:   # ~1 in 25 sites
        for i, ln in enumerate(f):
            if i % 25 == 0:
                p = ln.split("\t")
                p[9], p[11], p[12] = "3", str(min(int(p[11]), 3)), str(3 - min(int(p[11]), 3))
                o.write("\t".join(p))
    pysam.tabix_index(str(d / "Z.bed"), preset="bed", force=True)
    (d / "samples.tsv").write_text(f"sample\tCHG\nS1\t{syn['dir'] / 'S1.CHG.bedmethyl.gz'}\nS2\tS2.bed.gz\n"
                                   f"L\tL.bed.gz\nZ\tZ.bed.gz\n")
    args = _args(syn, tmp_path / "o" / "x", "--anchors", str(syn["anchors"]))
    args[args.index("--samples") + 1] = str(d / "samples.tsv")
    assert cli.main(args) == 0
    err = capsys.readouterr().err
    assert "sample L: median depth 3.0x" in err
    assert "sample Z: no usable window" in err and "sample L: no usable window" not in err
    cand = pd.read_csv(tmp_path / "o" / "x.candidates.tsv", sep="\t")
    assert cand[(cand["sample"] == "S2") & (cand.chrom == "chr6")].empty
    assert cand[cand["sample"] == "Z"].empty
    s1 = cand[(cand["sample"] == "S1") & (cand.chrom == "chr6")]     # S2 has no chr6, Z no data: L is compared
    assert ((s1.plant_specific_n >= 1) & s1.plant_specific.notna()).any()


def test_fragmented_assembly_fails(syn, tmp_path, capsys):
    fai = tmp_path / "contigs.fai"
    fai.write_text("".join(f"ctg{i}\t1000000\t0\t60\t61\n" for i in range(2000)))
    args = _args(syn, tmp_path / "o" / "x")
    args[args.index("--fai") + 1] = str(fai)
    assert cli.main(args) == 2
    assert "under two" in capsys.readouterr().err


def test_anchor_without_data_warns(syn, tmp_path, capsys):
    ex = tmp_path / "ex.bed"
    ex.write_text("chr1\t18000000\t22000000\n")                 # blank sample column = every sample
    assert cli.main(_args(syn, tmp_path / "o" / "x", "--anchors", str(syn["anchors"]), "--exclude", str(ex))) == 0
    err = capsys.readouterr().err
    assert "anchors without data" in err and "chr1_cen" in err


def test_exclusion_naming_no_sample_fails(syn, tmp_path, capsys):
    ex = tmp_path / "ex.bed"
    ex.write_text("chr1\t15000000\t25000000\tinv1\nchr2\t0\t1000\tS2\nchr3\t0\t1000\t.\n")
    assert cli.main(_args(syn, tmp_path / "o" / "x", "--exclude", str(ex))) == 2
    err = capsys.readouterr().err
    assert "inv1" in err and "leave column 4 blank (or '.') to exclude a region in every sample" in err


def test_calibration_messages(syn, tmp_path, capsys):
    o = tmp_path / "o"                                        # one folder: the runs share the bin cache
    scaf = tmp_path / "scaf.bed"
    scaf.write_text("scaf1\t100000\t200000\tsc\n")            # in the .fai but too short to analyse
    assert cli.main(_args(syn, o / "a", "--anchors", str(scaf))) == 0
    err = capsys.readouterr().err
    assert "no anchor on an analysed chromosome: uncalibrated run" in err
    two = tmp_path / "two.bed"
    two.write_text("chr1\t19750000\t20250000\tc1\nchr2\t27750000\t28250000\tc2\n")
    assert cli.main(_args(syn, o / "b", "--anchors", str(two))) == 0
    err = capsys.readouterr().err
    assert "too few anchors with data to calibrate (need 3)" in err and "no consistent" not in err
    assert cli.main(_args(syn, o / "c", "--anchors", str(syn["anchors_wrong"]))) == 0
    err = capsys.readouterr().err
    assert "no consistent methylation signature" in err and "too few anchors" not in err


@pytest.fixture(scope="module")
def up_run(tmp_path_factory):
    """Synthetic genome whose centromeric CRMs gain CHG (0.45 -> 0.75), run with the true anchors."""
    d = tmp_path_factory.mktemp("up")
    s = synth.make(d / "data", crm_domain=0.75)
    assert cli.main(_args(s, d / "r" / "x", "--anchors", str(s["anchors"]))) == 0
    return s, d / "r" / "x"


def test_end_to_end_up(up_run):
    s, out = up_run
    head = open(f"{out}.candidates.bed").readline()
    assert "direction=up(anchors)" in head and "CHG:LTR:CRM:up" in head
    sig = pd.read_csv(f"{out}.signature.tsv", sep="\t")
    assert ((sig.context == "CHG") & (sig.cls == "LTR:CRM") & (sig.direction == "up") & sig.in_signature).any()
    cand = pd.read_csv(f"{out}.candidates.tsv", sep="\t")
    top = cand[cand["rank"] == 1].set_index(["sample", "chrom"])
    for smp in synth.SAMPLES:
        for c, x in s["cen"].items():
            want = s["shift"][2] if (smp, c) == s["shift"][:2] else x
            r = top.loc[(smp, c)]
            assert abs(r.peak - want) <= 1e6 and r.label == "centromere_like", (smp, c, r.peak, r.label)
