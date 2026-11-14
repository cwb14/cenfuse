import gzip
import logging
import shutil
import struct
from pathlib import Path

import numpy as np
import pysam
import pytest

import meth_synth as synth
from cenfuse.meth import io

ROWS = [("c1", 10, "m", 5, 2, 1), ("c1", 10, "h", 5, 1, 2), ("c1", 20, "m", 2, 1, 0), ("c1", 30, "m", 4, 0, 0)]


def _bedmethyl(path, rows, space=False):
    with open(path, "w") as f:
        for c, p, code, nv, nm, nh in rows:
            head = f"{c}\t{p}\t{p + 1}\t{code}\t{nv}\t+\t{p}\t{p + 1}\t255,0,0"
            tail = [nv, f"{100 * nm / nv:.2f}", nm, nv - nm - nh, nh, 0, 0, 0, 0]
            f.write(head + ("\t" + " ".join(map(str, tail)) if space else "\t" + "\t".join(map(str, tail))) + "\n")
    pysam.tabix_index(str(path), preset="bed", force=True)
    return f"{path}.gz"


@pytest.mark.parametrize("space", [False, True])
def test_bedmethyl_tab_and_space_columns_skip_h_rows(tmp_path, space):
    f = _bedmethyl(tmp_path / "x.bed", ROWS, space)
    with pysam.TabixFile(f) as tb:
        pos, v, m, bad = io.read_bedmethyl(tb, "c1", 0, 100, 3)
        assert io.read_bedmethyl(tb, "absent", 0, 100, 3)[0].size == 0
    assert pos.tolist() == [10, 30]          # the 'h' row and the 2-call site are dropped
    assert v.tolist() == [5, 4]
    assert m.tolist() == [3, 0]              # N_mod + N_other_mod of the 'm' row
    assert bad == 0


def test_encode_bedmethyl_reads_like_modkit(tmp_path):
    """ENCODE bedMethyl (coverage, percent methylated; no 'h' rows) gives the modkit fixture's sites and calls."""
    with open(tmp_path / "e.bed", "w") as f:
        for c, p, code, nv, nm, nh in ROWS:
            if code == "m":
                f.write(f"{c}\t{p}\t{p + 1}\t.\t{nv}\t+\t{p}\t{p + 1}\t0,0,0\t{nv}\t{100 * (nm + nh) / nv:.0f}\n")
    pysam.tabix_index(str(tmp_path / "e.bed"), preset="bed", force=True)
    with pysam.TabixFile(f"{tmp_path / 'e.bed'}.gz") as tb:
        pos, v, m, bad = io.read_bedmethyl(tb, "c1", 0, 100, 3)
    assert pos.tolist() == [10, 30] and v.tolist() == [5, 4] and m.tolist() == [3, 0] and bad == 0


def test_other_tables_stop_with_the_formats_read(tmp_path):
    (tmp_path / "cx.bed").write_text("c1\t10\t11\t+\t3\t2\tCG\n")       # e.g. a Bismark CX report
    pysam.tabix_index(str(tmp_path / "cx.bed"), preset="bed", force=True)
    with pysam.TabixFile(f"{tmp_path / 'cx.bed'}.gz") as tb, pytest.raises(io.InputError, match="ENCODE bedMethyl"):
        io.read_bedmethyl(tb, "c1", 0, 100, 3)


def test_bedmethyl_region_is_half_open(tmp_path):
    f = _bedmethyl(tmp_path / "x.bed", ROWS)
    with pysam.TabixFile(f) as tb:
        assert io.read_bedmethyl(tb, "c1", 0, 30, 3)[0].tolist() == [10]
        assert io.read_bedmethyl(tb, "c1", 30, 31, 3)[0].tolist() == [30]


def test_fai(tmp_path):
    (tmp_path / "g.fai").write_text("a\t100\t0\t60\t61\nb\t50\t0\t60\t61\n")
    assert io.read_fai(tmp_path / "g.fai") == {"a": 100, "b": 50}
    with pytest.raises(io.InputError):
        io.read_fai(tmp_path / "missing.fai")


def test_samples_table(tmp_path):
    f = _bedmethyl(tmp_path / "s1.CHG.bed", ROWS)
    (tmp_path / "s.tsv").write_text(f"sample\tCHG\tCG\nS1\t{f.split('/')[-1]}\t\n")
    s = io.read_samples(tmp_path / "s.tsv")
    assert s[0].name == "S1" and list(s[0].files) == ["CHG"] and s[0].files["CHG"].exists()
    (tmp_path / "bad.tsv").write_text("sample\tCG\nS1\tx\n")
    with pytest.raises(io.InputError, match="CHG"):
        io.read_samples(tmp_path / "bad.tsv")
    (tmp_path / "noidx.bed.gz").write_bytes(gzip.compress(b"c1\t1\t2\n"))
    (tmp_path / "n.tsv").write_text("sample\tCHG\nS1\tnoidx.bed.gz\n")
    with pytest.raises(io.InputError, match="tabix index"):
        io.read_samples(tmp_path / "n.tsv")
    (tmp_path / "d.tsv").write_text(f"sample\tCHG\nS1\t{f}\nS1\t{f}\n")
    with pytest.raises(io.InputError, match="duplicate"):
        io.read_samples(tmp_path / "d.tsv")


@pytest.mark.parametrize("name", ["", "a/b", "S 1", "S1 ", "\u00a0S1"])
def test_samples_table_rejects_bad_names(tmp_path, name):
    f = _bedmethyl(tmp_path / "s1.CHG.bed", ROWS)
    (tmp_path / "s.tsv").write_text(f"sample\tCHG\n{name}\t{f}\n")
    with pytest.raises(io.InputError, match="sample name"):
        io.read_samples(tmp_path / "s.tsv")


def test_gff_attributes_keep_hash_and_fasta_is_cut(tmp_path):
    g = tmp_path / "x.gff3"
    g.write_text("##gff-version 3\n"
                 "c1\tLTRquest\tLTR_retrotransposon\t11\t20\t.\t+\t.\tID=e1;family_rep=c1:1-2#LTR/Gypsy/CRM;"
                 "superfamily=Gypsy;clade=CRM;depth=2\n"
                 "c1\tLTRquest\tlong_terminal_repeat\t11\t13\t.\t+\t.\tID=e1.l;Parent=e1\n"
                 "c1\tLTRquest\tLTR_retrotransposon\t31\t40\t.\t+\t.\tID=e2;superfamily=Copia;clade=unknown;depth=0\n"
                 "c1\tLTRquest\tLTR_retrotransposon\t51\t60\t.\t+\t.\tID=e3;superfamily=unknown;clade=mixture\n"
                 "##FASTA\n>c1\nACGT\n")
    t = io.read_ltrquest(g)
    assert t.start.tolist() == [10, 30, 50] and t.end.tolist() == [20, 40, 60]
    assert t.cls.tolist() == ["LTR:CRM", "LTR:Copia_unclassified", "LTR:unclassified"]
    assert t.depth.tolist() == [2, 0, 0]


def test_ltrquest_without_elements_fails(tmp_path):
    g = tmp_path / "x.gff3"
    g.write_text("c1\tx\tgene\t1\t5\t.\t+\t.\tID=g\n")
    with pytest.raises(io.InputError, match="LTR_retrotransposon"):
        io.read_ltrquest(g)


def test_edta_classes(tmp_path):
    g = tmp_path / "e.TEanno.gff3"
    g.write_text("c1\tEDTA\trepeat_region\t1\t100\t.\t+\t.\tID=r1\n"
                 "c1\tEDTA\tGypsy_LTR_retrotransposon\t5\t95\t.\t+\t.\tID=a;Classification=LTR/Gypsy\n"
                 "c1\tEDTA\tlong_terminal_repeat\t5\t15\t.\t+\t.\tID=b;Classification=LTR/Gypsy\n"
                 "c1\tEDTA\thelitron\t200\t300\t.\t+\t.\tID=c;Classification=DNA/Helitron\n")
    t = io.read_edta(g)
    assert t.cls.tolist() == ["TE:LTR/Gypsy", "TE:DNA/Helitron"] and t.start.tolist() == [4, 199]


def test_satellites_engine_gff_and_bed(tmp_path):
    g = tmp_path / "g_satellite_arrays.gff"
    g.write_text("##gff-version 3\nc1\tcenfuse\tsatellite_DNA\t1\t100\t.\t+\t.\tID=s1;family=fam1;rdna=true\n"
                 "c1\tcenfuse\tsatellite_monomer\t1\t50\t.\t+\t.\tID=m1;family=fam1\n")
    assert io.read_satellites(g).cls.tolist() == ["sat:fam1"]
    b = tmp_path / "s.bed"
    b.write_text("track name=x\nc1\t0\t100\tCentO\nc1\t200\t300\t\n")
    assert io.read_satellites(b).cls.tolist() == ["sat:CentO"]
    (tmp_path / "n.bed").write_text("c1\t0\t100\n")
    with pytest.raises(io.InputError, match="family"):
        io.read_satellites(tmp_path / "n.bed")


def test_genes_only_gene_features(tmp_path):
    g = tmp_path / "g.gff3"
    g.write_text("c1\tx\tgene\t1\t10\t.\t+\t.\tID=g\nc1\tx\tmRNA\t1\t10\t.\t+\t.\tID=t;Parent=g\n")
    assert len(io.read_genes(g)) == 1


def test_anchor_sources_and_midpoints(tmp_path):
    h = tmp_path / "a.centromeres.bed"
    h.write_text("#cenfuse-hic 0.1.1 bin_size=500000\nc1\t1000\t5000\tc1_cen\t300\t.\t1000\t1500\t0,114,178\n")
    d, src = io.read_anchors(h)
    assert src == "hic" and d.mid.tolist() == [1250]
    r = tmp_path / "x_centromeric_repeat.gff"
    r.write_text("c1\tcenfuse\tcentromere\t101\t300\t.\t.\t.\tID=cen_c1;confidence=high\n"
                 "c1\tcenfuse\tsatellite_DNA\t101\t300\t.\t+\t.\tID=s;family=f\n")
    d, src = io.read_anchors(r)
    assert src == "satellite" and len(d) == 1 and d.mid.tolist() == [200]
    u = tmp_path / "chip.bed"
    u.write_text("c1\t0\t1000\n")
    d, src = io.read_anchors(u)
    assert src == "user" and d.mid.tolist() == [500]


def test_bed_with_text_in_thick_columns_is_accepted(tmp_path):
    b = tmp_path / "peaks.bed"
    b.write_text("c1\t0\t1000\tp1\t0\t+\tfoo\tbar\nc1\t2000\t3000\tp2\t0\t+\t2100\t2300\n")
    d = io.read_intervals(b)
    assert np.isnan(d.thick_start[0]) and np.isnan(d.thick_end[0]) and d.thick_start[1] == 2100
    a, _ = io.read_anchors(b)
    assert a.mid.tolist() == [500, 2200]


def test_exclusions_sample_column(tmp_path):
    e = tmp_path / "e.bed"
    e.write_text("c1\t0\t10\tS1\nc1\t20\t30\t.\nc1\t40\t50\n")
    assert io.read_exclusions(e)["sample"].tolist() == ["S1", "", ""]
    assert np.all(io.read_exclusions(e).end.to_numpy() == [10, 30, 50])


def test_anchors_missing_file_is_input_error(tmp_path):
    with pytest.raises(io.InputError):
        io.read_anchors(tmp_path / "absent.centromeres.bed")


def test_anchors_gzipped_hic_matches_plain(tmp_path):
    line = "#cenfuse-hic 0.1.1 bin_size=500000\nc1\t1000\t5000\tc1_cen\t300\t.\t1000\t1500\t0,114,178\n"
    plain = tmp_path / "a.centromeres.bed"
    plain.write_text(line)
    gz = tmp_path / "a.centromeres.bed.gz"
    gz.write_bytes(gzip.compress(line.encode()))
    d_plain, src_plain = io.read_anchors(plain)
    d_gz, src_gz = io.read_anchors(gz)
    assert src_gz == "hic" == src_plain
    assert d_gz.mid.tolist() == d_plain.mid.tolist() == [1250]


def test_csi_only_track_opens_and_is_read(tmp_path):
    """Sequences over 512 Mb need a .csi index (tabix -C); pysam looks for .tbi only unless told."""
    p = tmp_path / "big.bed"
    with open(p, "w") as f:
        for c, s, code, nv, nm, nh in ROWS:
            f.write(f"{c}\t{s}\t{s + 1}\t{code}\t{nv}\t+\t{s}\t{s + 1}\t255,0,0\t"
                    + "\t".join(map(str, [nv, f"{100 * nm / nv:.2f}", nm, nv - nm - nh, nh, 0, 0, 0, 0])) + "\n")
    pysam.tabix_index(str(p), preset="bed", force=True, csi=True)
    f = f"{p}.gz"
    assert not Path(f + ".tbi").exists() and Path(f + ".csi").exists()
    assert io.tabix_contigs(f) == ["c1"]
    assert io.table_layout(f)[0] == "modkit"
    with io.open_tabix(f) as tb:
        assert io.read_bedmethyl(tb, "c1", 0, 100, 3)[0].tolist() == [10, 30]


def test_genes_from_gtf_and_from_transcripts_without_gene_lines(tmp_path):
    """GTF has GFF's columns (1-based closed); a file without `gene` lines (e.g. potato DM v6.1) gives its mRNA or
    transcript lines."""
    gtf = tmp_path / "g.gtf"
    gtf.write_text('c1\tsrc\tgene\t101\t200\t.\t+\t.\tgene_id "g1";\n'
                   'c1\tsrc\ttranscript\t101\t180\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n')
    assert io.read_genes(gtf)[["start", "end"]].values.tolist() == [[100, 200]]
    mrna = tmp_path / "m.gff3"
    mrna.write_text("c1\tx\tmRNA\t1\t10\t.\t+\t.\tID=t1\nc1\tx\texon\t1\t5\t.\t+\t.\tParent=t1\n"
                    "c1\tx\tmRNA\t21\t30\t.\t+\t.\tID=t2\n")
    assert io.read_genes(mrna)[["start", "end"]].values.tolist() == [[0, 10], [20, 30]]
    gz = tmp_path / "t.gtf.gz"
    gz.write_bytes(gzip.compress(b'c1\tsrc\ttranscript\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'))
    assert io.read_genes(gz)[["start", "end"]].values.tolist() == [[100, 200]]


def _pbcpg(path, chrom="c1"):
    """pb-CpG-tools BED (9 columns), bgzip + tabix: cenfuse's CpG line reads it, the methylation engine does not."""
    path.write_text("".join(f"{chrom}\t{p}\t{p + 2}\t50.0\tTotal\t20\t10\t10\t50.0\n" for p in range(100, 1000, 100)))
    pysam.tabix_index(str(path), preset="bed", force=True)
    return f"{path}.gz"


def test_engine_drops_cg_files_that_are_not_bedmethyl(tmp_path, caplog, monkeypatch):
    """The engine reads every file as bedMethyl: a pb-CpG-tools CG next to a modkit CHG is dropped with a warning;
    read for cenfuse's CpG line (require_chg False), it stays."""
    chg, cg = _bedmethyl(tmp_path / "s1.CHG.bed", ROWS), _pbcpg(tmp_path / "pb.bed")
    (tmp_path / "s.tsv").write_text(f"sample\tCHG\tCG\nS1\t{chg}\t{cg}\n")
    monkeypatch.setattr(io.log, "propagate", True)
    with caplog.at_level(logging.WARNING, logger="cenfuse-meth"):
        s = io.read_samples(tmp_path / "s.tsv")
    assert list(s[0].files) == ["CHG"] and "pb.bed.gz" in caplog.text and "9 columns" in caplog.text
    assert list(io.read_samples(tmp_path / "s.tsv", require_chg=False)[0].files) == ["CHG", "CG"]


def test_methylation_line_takes_gtf_genes_and_a_pbcpg_cg(syn, tmp_path):
    """cenfuse's methylation line (methfuse.prepare) with two inputs cenfuse accepts elsewhere: genes as GTF and a
    pb-CpG-tools CG file beside modkit CHG. Neither switches the line off."""
    from types import SimpleNamespace

    from cenfuse import methfuse
    gtf = tmp_path / "genes.gtf"
    with open(syn["genes"]) as f, open(gtf, "w") as o:
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) == 9 and p[2] == "gene":
                o.write("\t".join(p[:8] + [f'gene_id "{p[8][3:]}";']) + "\n")
    d, pb = syn["dir"], _pbcpg(tmp_path / "pb.bed", "chr1")
    (tmp_path / "s.tsv").write_text(f"sample\tCHG\tCG\nS1\t{d / 'S1.CHG.bedmethyl.gz'}\t{pb}\n"
                                    f"S2\t{d / 'S2.CHG.bedmethyl.gz'}\t\nS3\t{d / 'S3.CHG.bedmethyl.gz'}\t\n")
    g = SimpleNamespace(name="g", fasta=str(syn["fai"])[:-4], chroms=list(synth.CHROMS), meth=tmp_path / "s.tsv",
                        tes=syn["ltrquest"], genes=gtf, exclude=None)
    prep = methfuse.prepare(g, syn["satellites"], tmp_path / "cache", 1)
    assert "gene" in prep.cmap.names and [list(s.files) for s in prep.samples] == [["CHG"]] * 3


def _cut_blocks(src, dst, frac):
    """The first frac of src's BGZF blocks: what a killed `modkit pileup | bgzip` leaves."""
    data, offs, o = Path(src).read_bytes(), [], 0
    while o < len(data):
        offs.append(o)
        o += struct.unpack("<H", data[o + 16:o + 18])[0] + 1          # BSIZE: block size - 1
    Path(dst).write_bytes(data[:offs[int(len(offs) * frac)]])


MANY = [("c1", p, "m", 10, 5, 0) for p in range(0, 200_000, 10)]     # 20,000 rows: about 20 BGZF blocks


def test_samples_table_rejects_a_truncated_bgzip_file(tmp_path):
    """Cut at a block boundary and indexed, a file reads in part without an error from htslib: the missing
    end-of-file block gives it away."""
    full = _bedmethyl(tmp_path / "full.bed", list(MANY))
    _cut_blocks(full, tmp_path / "cut.bed.gz", 0.75)
    pysam.tabix_index(str(tmp_path / "cut.bed.gz"), preset="bed", force=True)
    (tmp_path / "ok.tsv").write_text("sample\tCHG\nA\tfull.bed.gz\n")
    assert io.read_samples(tmp_path / "ok.tsv")[0].files["CHG"].name == "full.bed.gz"
    (tmp_path / "s.tsv").write_text("sample\tCHG\nA\tfull.bed.gz\nB\tcut.bed.gz\n")
    with pytest.raises(io.InputError, match="cut.bed.gz.*truncated|truncated.*cut.bed.gz"):
        io.read_samples(tmp_path / "s.tsv")


def test_read_error_is_an_input_error_naming_file_and_region(tmp_path):
    """Cut inside a block, its old index beside it: htslib's read error names neither file nor region."""
    from cenfuse.meth import bins
    full = Path(_bedmethyl(tmp_path / "full.bed", list(MANY)))
    cut = tmp_path / "cut.bed.gz"
    cut.write_bytes(full.read_bytes()[:full.stat().st_size // 2])
    shutil.copy(f"{full}.tbi", f"{cut}.tbi")
    with pytest.raises(io.InputError, match=r"cut\.bed\.gz.*c1:1-200000"):
        bins.bin_chromosome(cut, "c1", 200_000, None, 0, 1)


def test_missing_index_message_says_to_bgzip_first(tmp_path):
    (tmp_path / "plain.bed.gz").write_bytes(gzip.compress(b"c1\t1\t2\n"))
    (tmp_path / "s.tsv").write_text("sample\tCHG\nS1\tplain.bed.gz\n")
    with pytest.raises(io.InputError, match="bgzip"):
        io.read_samples(tmp_path / "s.tsv")


def _lines(rows):
    return [f"{c}\t{p}\t{p + 1}\t{code}\t{nv}\t+\t{p}\t{p + 1}\t255,0,0\t"
            + "\t".join(map(str, [nv, f"{100 * nm / nv:.2f}", nm, nv - nm - nh, nh, 0, 0, 0, 0]))
            for c, p, code, nv, nm, nh in rows]


def _tabix(path, lines):
    path.write_text("".join(x + "\n" for x in lines))
    pysam.tabix_index(str(path), preset="bed", force=True)
    return pysam.TabixFile(f"{path}.gz")


def test_short_first_row_of_a_region_is_one_malformed_row(tmp_path):
    rows = _lines([("c1", p, "m", 10, 5, 0) for p in range(100, 1100, 10)])
    with _tabix(tmp_path / "x.bed", ["\t".join(rows[0].split("\t")[:6])] + rows[1:]) as tb:
        pos, v, m, bad = io.read_bedmethyl(tb, "c1", 0, 10_000, 3)
    assert (len(pos), int(v.sum()), int(m.sum()), bad) == (99, 990, 495, 1)


@pytest.mark.parametrize("code", ["h", "76792"])
def test_5hmc_rows_are_not_extra_cytosines(tmp_path, code):
    """modkit m + 5hmC rows (code h or ChEBI 76792), whole or cut to ENCODE's 11 columns: one site per cytosine."""
    rows = []
    for p in range(100, 1100, 10):
        rows += [("c1", p, "m", 10, 5, 1), ("c1", p, code, 10, 1, 5)]
    with _tabix(tmp_path / "x18.bed", _lines(rows)) as tb:
        pos, v, m, _ = io.read_bedmethyl(tb, "c1", 0, 10_000, 3)
    assert (len(pos), int(v.sum()), int(m.sum())) == (100, 1000, 600)          # 5mC + 5hmC of the m row
    with _tabix(tmp_path / "x11.bed", ["\t".join(x.split("\t")[:11]) for x in _lines(rows)]) as tb:
        pos, v, m, _ = io.read_bedmethyl(tb, "c1", 0, 10_000, 3)
    assert (len(pos), int(v.sum()), int(m.sum())) == (100, 1000, 500)          # the m row's percent: 5mC


def test_bed_with_more_than_twelve_columns(tmp_path):
    (tmp_path / "b13.bed").write_text("c1\t100\t200\tCEN\t0\t+\t100\t200\t0\t1\t100,\t0,\textra\n")
    assert io.read_satellites(tmp_path / "b13.bed")[["start", "end", "cls"]].values.tolist() == [[100, 200, "sat:CEN"]]
    (tmp_path / "b14.bed").write_text("c1\t100\t200\t5\t0\t1\t100\t200\t0\t1\t100\t0\t7\t8\n")
    assert io.read_intervals(tmp_path / "b14.bed")[["start", "end"]].values.tolist() == [[100, 200]]


@pytest.mark.parametrize("name", ["NA", "None", "nan", "null", "NaN"])
def test_samples_named_like_missing_values_keep_their_name(tmp_path, name):
    f = _bedmethyl(tmp_path / "s1.CHG.bed", ROWS)
    (tmp_path / "s.tsv").write_text(f"sample\tCHG\tCG\n{name}\t{f}\t\n")
    s = io.read_samples(tmp_path / "s.tsv")
    assert s[0].name == name and list(s[0].files) == ["CHG"]


def test_damaged_index_is_an_input_error(tmp_path):
    """The engine's checks name the file (cenfuse's CpG line keeps open_tabix's OSError, which it catches)."""
    f, cg = _bedmethyl(tmp_path / "x.bed", ROWS), _bedmethyl(tmp_path / "cg.bed", ROWS)
    Path(f"{f}.tbi").write_bytes(b"garbage" * 10)
    with pytest.raises(io.InputError, match=r"x\.bed\.gz"):
        io.tabix_contigs(f)
    with pytest.raises(OSError):
        io.open_tabix(f)
    (tmp_path / "s.tsv").write_text(f"sample\tCHG\tCG\nS1\t{cg}\t{f}\n")
    with pytest.raises(io.InputError, match=r"S1: CG file unreadable.*x\.bed\.gz"):
        io.read_samples(tmp_path / "s.tsv")


def test_exclusions_naming_no_sample_stop_the_methylation_line(syn, tmp_path):
    """An exclusion BED whose column 4 names no sample of the table excludes nothing; cenfuse's methylation line stops with
    a message naming it (as the standalone engine does) instead of running without the exclusion. A blank or known name
    is fine."""
    from types import SimpleNamespace

    from cenfuse import methfuse
    from cenfuse.util import CenfuseError
    d = syn["dir"]
    (tmp_path / "s.tsv").write_text("".join(f"{s}\t{d / f'{s}.CHG.bedmethyl.gz'}\n" for s in ("S1", "S2", "S3"))
                                    .join(["sample\tCHG\n", ""]))
    bad = tmp_path / "bad.bed"
    bad.write_text("chr5\t50000000\t52000000\tS9\n")
    g = SimpleNamespace(name="g", fasta=str(syn["fai"])[:-4], chroms=list(synth.CHROMS), meth=tmp_path / "s.tsv",
                        tes=syn["ltrquest"], genes=syn["genes"], exclude=bad)
    with pytest.raises(CenfuseError, match="names no sample"):
        methfuse.prepare(g, syn["satellites"], tmp_path / "cache", 1)
    g.exclude = syn["exclude"]                                          # column 4 = S2, a sample of the table
    assert methfuse.prepare(g, syn["satellites"], tmp_path / "cache2", 1) is not None
