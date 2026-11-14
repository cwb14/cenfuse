"""End-to-end: the full Hi-C engine `call` pipeline on a synthetic genome with planted centromeres."""
from __future__ import annotations

import bz2
import gzip
import logging
import os
from pathlib import Path

import pandas as pd
import pytest

import cenfuse.hic.cli as cli
from cenfuse.hic.cli import main
from cenfuse.hic.io import InputError
from hic_synth import make_genome


@pytest.fixture(scope="module")
def done(tmp_path_factory):
    d = tmp_path_factory.mktemp("e2e")
    truth = make_genome(d / "g", gz=True)
    prefix = d / "out" / "x"
    main(["call", "--pairs", str(truth["pairs"]), "--fai", str(truth["fai"]), "-o", str(prefix),
          "--satellites", str(truth["sat"]), "--repeats", str(truth["sat"]), "-t", "4"])
    return truth, prefix


def test_outputs_exist(done):
    truth, prefix = done
    for suf in (".centromeres.bed", ".centromeres.tsv", ".score.bedgraph", ".telomere.bedgraph", ".qc.tsv",
                ".overview.pdf"):
        assert (prefix.parent / (prefix.name + suf)).stat().st_size > 0, suf
    assert "# read_class_split" in open(f"{prefix}.qc.tsv").read()


def test_calls_match_truth(done):
    truth, prefix = done
    t = pd.read_csv(f"{prefix}.centromeres.tsv", sep="\t").set_index("chrom")
    size = truth["size"]
    for k, nm in enumerate(truth["names"]):
        if nm in truth["unattached"]:
            assert t.loc[nm, "confidence"] in ("none", "low")
            continue
        assert abs(t.loc[nm, "peak_start"] // size - truth["cens"][k]) <= 1, nm
        assert t.loc[nm, "bin_size"] == 250_000
    assert t.loc["chr2", "peak_minus_sat_Mb"] == pytest.approx(-5.0, abs=0.5)    # functional vs shifted satellite
    assert "rDNA" in str(t.loc[truth["rdna_chrom"], "flags"])
    assert (t.confidence == "high").sum() >= 9


def test_cycle_unstable_flag_only_on_the_unattached_chromosome(done):
    """The default synthetic genome (seed 1) locates in a cycle in which only chr11 (no cluster) jumps far."""
    truth, prefix = done
    t = pd.read_csv(f"{prefix}.centromeres.tsv", sep="\t", keep_default_na=False).set_index("chrom")
    for nm in truth["names"]:
        assert ("cycle_unstable" in t.loc[nm, "flags"].split(",")) == (nm in truth["unattached"]), nm


def test_bed_is_the_anchor(done):
    truth, prefix = done
    lines = open(f"{prefix}.centromeres.bed").read().splitlines()
    assert lines[0].startswith("#cenfuse-hic") and all(len(l.split("\t")) == 9 for l in lines[1:])


def test_rerun_uses_cache(done, caplog):
    truth, prefix = done
    caplog.set_level("INFO", logger="cenfuse-hic")
    main(["call", "--pairs", str(truth["pairs"]), "--fai", str(truth["fai"]), "-o", str(prefix), "-t", "4"])
    assert any("cache hit" in r.message for r in caplog.records)


def test_name_mismatch_fails_fast(tmp_path):
    truth = make_genome(tmp_path / "g", n_trans=20_000, n_cis=20_000)
    fai = tmp_path / "renamed.fai"
    fai.write_text("".join("X" + l for l in open(truth["fai"])))
    with pytest.raises(SystemExit) as e:
        main(["call", "--pairs", str(truth["pairs"]), "--fai", str(fai), "-o", str(tmp_path / "o" / "x")])
    assert e.value.code == 2


def test_too_few_sequences_fails(tmp_path):
    fai = tmp_path / "g.fai"
    fai.write_text("a\t1000000\t0\t60\t61\n")
    p = tmp_path / "p.txt"
    p.write_text("a\t10\t+\ta\t500000\t-\n")
    with pytest.raises(SystemExit) as e:
        main(["call", "--pairs", str(p), "--fai", str(fai), "-o", str(tmp_path / "o" / "x")])
    assert e.value.code == 2


def test_truncated_gz_exits_2(tmp_path, capsys):
    truth = make_genome(tmp_path / "g", n_trans=20_000, n_cis=20_000, gz=True)
    bad = tmp_path / "bad.pairs.gz"
    data = open(truth["pairs"], "rb").read()
    bad.write_bytes(data[: len(data) // 2])
    with pytest.raises(SystemExit) as e:
        main(["call", "--pairs", str(bad), "--fai", str(truth["fai"]), "-o", str(tmp_path / "o" / "x")])
    assert e.value.code == 2
    assert "cenfuse-hic: error:" in capsys.readouterr().err


def _pairs_read_too_early(*a, **k):
    raise AssertionError("pairs were read before every input was checked")


def _stop_after_startup(*a, **k):
    raise InputError("stopped after the startup checks")


def _fai3(d, L=25_000_000):
    fai = d / "g.fai"
    fai.write_text("".join(f"chr{i}\t{L}\t0\t60\t61\n" for i in (1, 2, 3)))
    return fai


def test_compressed_or_binary_pairs_exit_2(tmp_path, capsys):
    """gzip data in a file named .pairs, a .bz2 file and a BAM: a clear error, not a traceback."""
    fai, line = _fai3(tmp_path), b"chr1\t10\t+\tchr2\t20\t-\n"
    for name, data in (("z.pairs", gzip.compress(line)), ("z.pairs.bz2", bz2.compress(line)),
                       ("z.bam", gzip.compress(b"BAM\x01" + bytes(range(256)) * 8))):
        (tmp_path / name).write_bytes(data)
        with pytest.raises(SystemExit) as e:
            main(["call", "--pairs", str(tmp_path / name), "--fai", str(fai), "-o", str(tmp_path / "o" / "x")])
        assert e.value.code == 2 and f"{name}: " in capsys.readouterr().err


@pytest.mark.parametrize("second", ["same", "empty", "header_only", "bad_layout"])
def test_every_pairs_file_is_checked_before_any_is_read(tmp_path, monkeypatch, capsys, second):
    """A file given twice (its contacts would count twice), an empty or header-only library, or an unknown layout among
    the replicates fails in seconds, before the first library is read."""
    fai, good = _fai3(tmp_path), tmp_path / "good.txt"
    good.write_text("chr1\t10\t+\tchr2\t20\t-\n")
    f2 = tmp_path / "second.pairs"
    f2.write_text({"same": "", "empty": "", "bad_layout": "r1 chr1 11 chr2 21 + - UU\n",
                   "header_only": "## pairs format v1.0\n#columns: readID chrom1 pos1 chrom2 pos2\n"}[second])
    monkeypatch.setattr(cli, "_base_matrix", _pairs_read_too_early)
    with pytest.raises(SystemExit) as e:
        main(["call", "--pairs", str(good), str(good if second == "same" else f2), "--fai", str(fai), "-o",
              str(tmp_path / "o" / "x")])
    assert e.value.code == 2 and "cenfuse-hic: error:" in capsys.readouterr().err


def test_no_usable_pairs_gives_the_counts(tmp_path, capsys):
    """Every pair below MAPQ 30 (e.g. Juicer's placeholder MAPQ 0): the error names the file and the MAPQ filter, not
    a sequence-name mismatch."""
    p = tmp_path / "q0.pairs"
    p.write_text("## pairs format v1.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 mapq1 mapq2\n"
                 + "".join(f"r{i}\tchr1\t{100 + i}\tchr2\t{500 + i}\t+\t-\t0\t0\n" for i in range(3)))
    with pytest.raises(SystemExit) as e:
        main(["call", "--pairs", str(p), "--fai", str(_fai3(tmp_path)), "-o", str(tmp_path / "o" / "x")])
    err = capsys.readouterr().err
    assert e.value.code == 2 and "q0.pairs: no usable pairs" in err and "3 with MAPQ < 30" in err


def test_cache_key_covers_the_middle_of_a_small_file(tmp_path):
    """A file of 1-2 MiB is fingerprinted whole: a same-size edit past its first MiB changes the key."""
    p = tmp_path / "p.txt"
    data = bytearray(b"a" * (3 << 19))                  # 1.5 MiB
    p.write_bytes(bytes(data))
    st = p.stat()
    k1 = cli.fast_checksum(p)
    data[5 << 18] = ord("b")                            # at 1.25 MiB
    p.write_bytes(bytes(data))
    os.utime(p, ns=(st.st_atime_ns, st.st_mtime_ns))
    assert cli.fast_checksum(p) != k1


@pytest.mark.parametrize("opt", ["--satellites", "--repeats"])
def test_missing_annotation_fails_before_reading_pairs(done, tmp_path, monkeypatch, capsys, opt):
    truth, _ = done
    monkeypatch.setattr(cli, "_base_matrix", _pairs_read_too_early)
    with pytest.raises(SystemExit) as e:
        main(["call", "--pairs", str(truth["pairs"]), "--fai", str(truth["fai"]), "-o", str(tmp_path / "x"),
              opt, str(truth["sat"]), str(tmp_path / "missing.gff")])
    assert e.value.code == 2
    err = capsys.readouterr().err
    assert "cenfuse-hic: error:" in err and "missing.gff" in err


def test_bad_gff_coordinate_fails_before_reading_pairs(done, tmp_path, monkeypatch, capsys):
    truth, _ = done
    g = tmp_path / "sat.gff3"
    g.write_text("##gff-version 3\nchr1\tsrc\tcentromere\tx\t200\t.\t.\t.\tID=a\n")
    monkeypatch.setattr(cli, "_base_matrix", _pairs_read_too_early)
    with pytest.raises(SystemExit) as e:
        main(["call", "--pairs", str(truth["pairs"]), "--fai", str(truth["fai"]), "-o", str(tmp_path / "x"),
              "--satellites", str(g)])
    assert e.value.code == 2
    assert "sat.gff3: line 2" in capsys.readouterr().err


def test_fasta_as_fai_exits_2_naming_the_line(tmp_path, capsys):
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\nACGT\n")
    p = tmp_path / "p.txt"
    p.write_text("chr1\t10\t+\tchr1\t500\t-\n")
    with pytest.raises(SystemExit) as e:
        main(["call", "--pairs", str(p), "--fai", str(fa), "-o", str(tmp_path / "x")])
    assert e.value.code == 2
    err = capsys.readouterr().err
    assert "cenfuse-hic: error:" in err and "g.fa: line 1" in err


def test_annotation_without_intervals_is_warned(done, tmp_path, monkeypatch, caplog):
    truth, _ = done
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\nACGTACGT\n")
    monkeypatch.setattr(cli, "_base_matrix", _stop_after_startup)
    caplog.set_level("INFO", logger="cenfuse-hic")
    with pytest.raises(SystemExit):
        main(["call", "--pairs", str(truth["pairs"]), "--fai", str(truth["fai"]), "-o", str(tmp_path / "x"),
              "--repeats", str(fa)])
    assert any("g.fa" in r.message and "no BED/GFF3 intervals" in r.message and r.levelname == "WARNING"
               for r in caplog.records)


def test_repeats_with_cool_input_is_warned(done, tmp_path, monkeypatch, caplog):
    truth, _ = done
    cool = tmp_path / "x.cool"
    cool.write_bytes(b"")                       # never opened: the run stops after the startup checks
    rep = tmp_path / "r.bed"
    rep.write_text("chr1\t10\t20\n")
    monkeypatch.setattr(cli, "_base_matrix", _stop_after_startup)
    caplog.set_level("INFO", logger="cenfuse-hic")
    with pytest.raises(SystemExit):
        main(["call", "--pairs", str(cool), "--fai", str(truth["fai"]), "-o", str(tmp_path / "x"),
              "--repeats", str(rep)])
    assert any("--repeats" in r.message and ".cool" in r.message and r.levelname == "WARNING"
               for r in caplog.records)


def test_zero_threads_exits_2(done):
    truth, prefix = done
    with pytest.raises(SystemExit) as e:
        main(["call", "--pairs", str(truth["pairs"]), "--fai", str(truth["fai"]), "-o", str(prefix), "-t", "0"])
    assert e.value.code == 2


def test_corrupt_cache_is_rebuilt(done, tmp_path, caplog):
    truth, prefix = done
    import shutil
    out = tmp_path / "copy" / "x"
    shutil.copytree(str(prefix) + ".cache", str(out) + ".cache")
    for f in Path(str(out) + ".cache").glob("*.npz"):
        f.write_bytes(b"not a zip")
    caplog.set_level("INFO", logger="cenfuse-hic")
    main(["call", "--pairs", str(truth["pairs"]), "--fai", str(truth["fai"]), "-o", str(out), "-t", "4"])
    msgs = [r.message for r in caplog.records]
    assert any("unreadable" in m and "rebuilding" in m for m in msgs)
    assert not any("cache hit" in m for m in msgs)
    assert Path(f"{out}.centromeres.tsv").stat().st_size > 0


def test_uncalled_row_is_none_with_empty_coordinates(done):
    truth, prefix = done
    t = pd.read_csv(f"{prefix}.centromeres.tsv", sep="\t").set_index("chrom")
    nm = truth["unattached"][0]
    assert t.loc[nm, "confidence"] == "none"
    assert pd.isna(t.loc[nm, "peak_start"]) and pd.isna(t.loc[nm, "core_start"])
    assert t.confidence.isin(["high", "medium", "low", "none"]).all()


def test_bootstrap_interval_rounds_outwards():
    import numpy as np

    from cenfuse.hic.cli import _boot_interval
    assert _boot_interval(np.array([3, 4])) == (3, 4)                   # percentiles 3.025, 3.975 (int: 3, 3)
    assert _boot_interval(np.array([10] * 97 + [11, 12, 13])) == (10, 11)   # 97.5th percentile 10.525 (int: 10)
    assert _boot_interval(np.array([7] * 50)) == (7, 7)


def test_over_median_skips_empty_and_nonpositive_chromosomes():
    import warnings

    import numpy as np

    from cenfuse.hic.bins import Bins
    from cenfuse.hic.cli import _over_median
    b = Bins({"a": 400, "b": 400, "c": 400}, 100)
    S = np.array([1, 2, 2, 4, np.nan, np.nan, np.nan, np.nan, 0, 0, 0, 1.0])     # a: median 2; b: all NaN; c: median 0
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        r = _over_median(S, b, np.ones(3, bool))
    assert np.allclose(r[:4], [0.5, 1, 1, 2]) and np.isnan(r[4:]).all()


@pytest.fixture(scope="module")
def plus(done, tmp_path_factory):
    """One run with two extra .fai sequences: a 1-Mb contig (the standalone engine's chromosome-set rule leaves it out
    of the bin rule, and it is shorter than 20 x 250 kb) and a 30-Mb chrX with no pairs.
    Log records are collected on the root logger (the cli replaces the cenfuse-hic logger's own handlers)."""
    truth, _ = done
    d = tmp_path_factory.mktemp("plus")
    fai = d / "plus.fai"
    fai.write_text(open(truth["fai"]).read() + "tiny_contig\t1000000\t0\t60\t61\nchrX\t30000000\t0\t60\t61\n")
    records = []
    h = logging.Handler(logging.INFO)
    h.emit = records.append
    root = logging.getLogger()
    root.addHandler(h)
    try:
        main(["call", "--pairs", str(truth["pairs"]), "--fai", str(fai), "-o", str(d / "o" / "x"), "-t", "4"])
    finally:
        root.removeHandler(h)
    return truth, d / "o" / "x", records


def test_short_sequence_is_warned(plus):
    truth, prefix, records = plus
    assert any("tiny_contig" in r.getMessage() and r.levelname == "WARNING" for r in records)


def test_chromosome_without_contacts_is_warned_and_flagged(plus):
    truth, prefix, records = plus
    assert any("chrX" in r.getMessage() and "no_contacts" in r.getMessage() and r.levelname == "WARNING"
               for r in records)
    t = pd.read_csv(f"{prefix}.centromeres.tsv", sep="\t", keep_default_na=False).set_index("chrom")
    assert t.loc["chrX", "flags"] == "no_contacts" and t.loc["chrX", "confidence"] == "none"
    assert t.loc["chrX", "peak_start"] == "" and "tiny_contig" not in t.index
    assert int(t.loc["chr1", "bin_size"]) == 250_000          # as in 0.9
    assert not any("no_contacts" in t.loc[nm, "flags"] for nm in truth["names"])
    assert (t.loc[truth["names"], "confidence"] == "high").sum() >= 9


def test_fewer_than_three_chromosomes_with_contacts_fails(done, tmp_path, capsys):
    truth, _ = done
    lines = open(truth["fai"]).read().splitlines(keepends=True)[:2]          # chr1, chr2 have pairs
    fai = tmp_path / "two.fai"
    fai.write_text("".join(lines) + "chrX\t30000000\t0\t60\t61\nchrY\t30000000\t0\t60\t61\n")
    with pytest.raises(SystemExit) as e:
        main(["call", "--pairs", str(truth["pairs"]), "--fai", str(fai), "-o", str(tmp_path / "o" / "x"), "-t", "2"])
    assert e.value.code == 2
    err = capsys.readouterr().err
    assert "cenfuse-hic: error:" in err and "contacts" in err


def test_replicate_libraries_are_merged(tmp_path):
    """0.10.1: two halves of one library given as replicates reproduce the full library's calls (one summed matrix)."""
    t = make_genome(tmp_path / "g", n_trans=150_000, n_cis=150_000)
    lines = Path(t["pairs"]).read_text().splitlines(keepends=True)
    h1, h2 = tmp_path / "h1.txt", tmp_path / "h2.txt"
    h1.write_text("".join(lines[0::2]))
    h2.write_text("".join(lines[1::2]))
    main(["call", "--pairs", str(t["pairs"]), "--fai", str(t["fai"]), "-o", str(tmp_path / "full" / "x"), "-t", "2"])
    main(["call", "--pairs", str(h1), str(h2), "--fai", str(t["fai"]), "-o", str(tmp_path / "rep" / "x"), "-t", "2"])
    a = pd.read_csv(tmp_path / "full" / "x.centromeres.tsv", sep="\t")
    b = pd.read_csv(tmp_path / "rep" / "x.centromeres.tsv", sep="\t")
    pd.testing.assert_frame_equal(a, b)
    qc = (tmp_path / "rep" / "x.qc.tsv").read_text()
    assert "files" in qc and "files" not in (tmp_path / "full" / "x.qc.tsv").read_text()


def test_replicates_in_different_pair_layouts(tmp_path, caplog):
    """A 4DN file with MAPQ columns and a simple 6-column file as replicates: each is read with its own layout; MAPQ
    is reported unchecked (one file has none) and nothing crashes."""
    t = make_genome(tmp_path / "g", n_trans=60_000, n_cis=60_000)
    lines = Path(t["pairs"]).read_text().splitlines()
    four = tmp_path / "a.pairs"
    with open(four, "w") as f:
        f.write("## pairs format v1.0.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type mapq1 mapq2\n")
        for i, line in enumerate(lines[0::2]):
            c1, p1, s1, c2, p2, s2 = line.split("\t")
            f.write(f"r{i}\t{c1}\t{p1}\t{c2}\t{p2}\t{s1}\t{s2}\tUU\t60\t60\n")
    six = tmp_path / "b.txt"
    six.write_text("\n".join(lines[1::2]) + "\n")
    with caplog.at_level(logging.WARNING, logger="cenfuse-hic"):
        main(["call", "--pairs", str(four), str(six), "--fai", str(t["fai"]), "-o", str(tmp_path / "o" / "x"), "-t", "1"])
    assert "uniqueness (MAPQ >= 30) was not checked" in caplog.text
    assert (tmp_path / "o" / "x.centromeres.tsv").exists()
