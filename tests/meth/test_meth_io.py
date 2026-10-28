import gzip

import numpy as np
import pysam
import pytest

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
