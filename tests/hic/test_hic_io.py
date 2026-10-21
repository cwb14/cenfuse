from __future__ import annotations

import gzip
import warnings

import numpy as np
import pytest
from pandas.errors import ParserWarning

from cenfuse.hic import params
from cenfuse.hic.io import InputError, PairStats, iter_pairs, pairs_format, read_fai, read_intervals


def write(path, text, gz=False):
    op = gzip.open if gz else open
    with op(path, "wt") as f:
        f.write(text)
    return path


def collect(path, chroms):
    st = PairStats()
    out = [np.concatenate(x) for x in zip(*iter_pairs(path, chroms, st))] or [np.array([])] * 4
    return out, st


def test_read_fai_keeps_file_order(tmp_path):
    p = write(tmp_path / "g.fai", "chrB\t500\t0\t60\t61\nchrA\t300\t600\t60\t61\n")
    assert list(read_fai(p).items()) == [("chrB", 500), ("chrA", 300)]


def test_simple6_zero_based_and_chrom_filter(tmp_path):
    p = write(tmp_path / "p.txt", "c1\t10\t+\tc2\t20\t-\nc1\t5\t+\tcX\t7\t+\n")
    (ci1, p1, ci2, p2), st = collect(p, {"c1": 100, "c2": 100})
    assert pairs_format(p).kind == "simple6"
    assert list(ci1) == [0] and list(p1) == [10] and list(ci2) == [1] and list(p2) == [20]
    assert st.read == 2 and st.kept == 1 and st.dropped_chrom == 1 and not st.mapq_checked
    assert "cX" in st.unknown_chroms


def test_4dn_one_based_mapq_filter_gz(tmp_path):
    txt = ("## pairs format v1.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type mapq1 mapq2\n"
           "r1\tc1\t11\tc2\t21\t+\t-\tUU\t60\t60\n"
           "r2\tc1\t11\tc2\t21\t+\t-\tUU\t60\t3\n"
           "r3\tc1\t31\tc1\t41\t+\t-\tUU\t40\t30\n")
    p = write(tmp_path / "p.pairs.gz", txt, gz=True)
    (ci1, p1, ci2, p2), st = collect(p, {"c1": 100, "c2": 100})
    assert list(p1) == [10, 30] and list(p2) == [20, 40]
    assert st.mapq_checked and st.dropped_mapq == 1 and st.kept == 2 and st.cis == 1


def test_4dn_without_columns_line_fails(tmp_path):
    p = write(tmp_path / "bad.pairs", "## pairs format v1.0\nr1\tc1\t1\tc2\t2\t+\t-\n")
    with pytest.raises(InputError, match="#columns"):
        pairs_format(p)


def test_malformed_positions_counted(tmp_path):
    p = write(tmp_path / "p.txt", "c1\tx\t+\tc2\t20\t-\nc1\t5\t+\tc2\t7\t+\n")
    _, st = collect(p, {"c1": 100, "c2": 100})
    assert st.malformed == 1 and st.kept == 1


def test_wrong_width_line_counted_and_no_warning_leaks(tmp_path):
    p = write(tmp_path / "w.txt", "c1\t10\t+\tc2\t20\t-\nc1\t5\t+\tc2\t7\t+\textra\n")
    st = PairStats()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        out = list(iter_pairs(p, {"c1": 100, "c2": 100}, st))
    assert st.malformed == 1 and st.kept == 1 and st.read == 2
    assert len(out) == 1 and list(out[0][1]) == [10]
    assert not [x for x in w if issubclass(x.category, ParserWarning)]


def test_read_intervals_bed_and_gff(tmp_path):
    b = write(tmp_path / "a.bed", "c1\t10\t20\tsat1\nc1\t30\t40\n")
    g = write(tmp_path / "a.gff3", "##gff-version 3\nc1\tcenfuse\tcentromere\t101\t200\t.\t.\t.\tID=cen_c1\n"
                                   "c1\tcenfuse\tsatellite_DNA\t11\t20\t.\t+\t.\tID=a;family=F;rdna=45S_spacer\n")
    tb, tg = read_intervals(b), read_intervals(g)
    # bracket access: DataFrame.attrs is pandas' own metadata dict, not the "attrs" column
    assert list(tb.start) == [10, 30] and tb["attrs"][0]["Name"] == "sat1"
    assert list(tg.start) == [100, 10] and list(tg.end) == [200, 20]
    assert list(tg.type) == ["centromere", "satellite_DNA"] and tg["attrs"][1]["rdna"] == "45S_spacer"
    tn = read_intervals(g, attrs=False)                  # coordinates only (large repeat annotations)
    assert list(tn.start) == [100, 10] and tn["attrs"].isna().all()


def test_read_fai_rejects_a_fasta_naming_file_and_line(tmp_path):
    fa = write(tmp_path / "g.fa", ">chr1\nACGT\n")
    with pytest.raises(InputError, match=r"g\.fa: line 1"):
        read_fai(fa)
    bad = write(tmp_path / "bad.fai", "chr1\t500\t0\t60\t61\nchr2\tabc\t600\t60\t61\n")
    with pytest.raises(InputError, match=r"bad\.fai: line 2"):
        read_fai(bad)
    neg = write(tmp_path / "neg.fai", "chr1\t-5\t0\t60\t61\n")
    with pytest.raises(InputError, match=r"neg\.fai: line 1"):
        read_fai(neg)
    gz = tmp_path / "g.fa.gz"
    gz.write_bytes(gzip.compress(b">chr1\nACGT\n" * 100))
    with pytest.raises(InputError, match="g.fa.gz"):
        read_fai(gz)


def test_read_intervals_bad_coordinate_names_file_and_line(tmp_path):
    g = write(tmp_path / "x.gff3", "##gff-version 3\nc1\tsrc\tcentromere\t101\t200\t.\t.\t.\tID=a\n"
                                   "c1\tsrc\tcentromere\t1.5e3\t2000\t.\t.\t.\tID=b\n")
    with pytest.raises(InputError, match=r"x\.gff3: line 3"):
        read_intervals(g)
    b = write(tmp_path / "x.bed", "c1\t10\t20\nc1\tstart\tend\n")
    with pytest.raises(InputError, match=r"x\.bed: line 2"):
        read_intervals(b)
    binary = tmp_path / "y.bed"
    binary.write_bytes(gzip.compress(b"c1\t10\t20\n" * 100))
    with pytest.raises(InputError, match="y.bed"):
        read_intervals(binary)


def test_fasta_as_annotation_yields_no_intervals(tmp_path):
    fa = write(tmp_path / "g.fa", ">chr1\nACGTACGT\nACGT\n")
    assert read_intervals(fa).empty


def test_params_values():
    assert (params.MAPQ, params.MIN_EXPECTED_CEN, params.HALF_WIN, params.END_WIN, params.MIN_BINS) == (30, 25, 4, 8, 20)
    assert params.BIN_LADDER == (250_000, 500_000, 1_000_000, 2_000_000)
    # simulate: 0.5% substitutions, mates >= 1 kb apart, reads <= 150 bp
    assert (params.SIM_ERR, params.MIN_CIS_DIST, params.SIM_MAXLEN, params.SIM_DEFAULT_LEN) == (0.005, 1000, 150, 150)
    assert params.SIM_MIN_CHROM == 1000 and params.RABL_MIN_BINS == 5 and params.CYCLE_SPREAD == 1
