"""GENOME_DIR discovery (inputs.py): accepted annotation names, Hi-C formats cenfuse cannot read failing clearly
instead of being ignored, and the chromosome set (the length rule, the 0.8 rule of --classify-from targets, and the
user's <genome>.chromosomes.txt)."""
import logging

import numpy as np
import pytest

from cenfuse import inputs, methfuse
from cenfuse.pipeline import Run
from cenfuse.util import CenfuseError, log


def _genome(d, chroms=("c1", "c2")):
    (d / "g.fa").write_text("".join(f">{c}\nACGT\n" for c in chroms))
    (d / "g.fa.fai").write_text("".join(f"{c}\t4\t{i}\t4\t5\n" for i, c in enumerate(chroms)))


def test_gff_and_gtf_annotations_are_found(tmp_path):
    _genome(tmp_path)
    for sub, name in (("genes", "g.gtf.gz"), ("tes", "g.gff")):
        (tmp_path / sub).mkdir()
        (tmp_path / sub / name).write_text("")
    g = inputs.discover(tmp_path)[0]
    assert g.genes.name == "g.gtf.gz" and g.tes.name == "g.gff"


def test_two_gene_annotations_for_one_genome_fail(tmp_path):
    _genome(tmp_path)
    (tmp_path / "genes").mkdir()
    for n in ("g.gff3", "g.gff"):
        (tmp_path / "genes" / n).write_text("")
    with pytest.raises(CenfuseError, match="several inputs"):
        inputs.discover(tmp_path)


@pytest.mark.parametrize("name, hint", [("g.hic", "hictk convert"), ("g.bam", "pairtools parse"),
                                        ("g.validPairs", "allValidPairs")])
def test_unreadable_hic_formats_fail_with_a_recipe(tmp_path, name, hint):
    _genome(tmp_path)
    (tmp_path / "hic").mkdir()
    (tmp_path / "hic" / name).write_text("")
    with pytest.raises(CenfuseError, match=hint):
        inputs.discover(tmp_path)


def _assembly(d, name, lengths):
    """d/name.fa (a stub: discover reads only the .fai) and d/name.fa.fai with these sequence lengths."""
    (d / f"{name}.fa").write_text("".join(f">{c}\nA\n" for c in lengths))
    (d / f"{name}.fa.fai").write_text("".join(f"{c}\t{n}\t{i}\t{n}\t{n + 1}\n"
                                              for i, (c, n) in enumerate(lengths.items())))


ASM = {"c1": 300_000, "c2": 300_000, "c3": 300_000, "s1": 5_000, "s2": 5_000}
# 24 chromosomes and a 7.7-Mb scaffold, then scaffolds falling gradually to 1 kb: 0.8 takes all 669 sequences, 0.9 25
GRADUAL = {f"s{i}": int(n) for i, n in enumerate(np.concatenate([np.linspace(60e6, 15e6, 24), [7.68e6],
                                                                   np.geomspace(1.6e6, 1e3, 644)]))}


def test_chromosome_list_file_overrides_the_rule(tmp_path):
    _assembly(tmp_path, "asm", ASM)
    assert list(inputs.discover(tmp_path)[0].chroms) == ["c1", "c2", "c3"]          # the rule's set
    (tmp_path / "asm.chromosomes.txt").write_bytes(b"# mine\r\nc3 \r\n\r\nc1\r\n")   # a Windows editor's file
    g = inputs.discover(tmp_path)[0]
    assert list(g.chroms.items()) == [("c1", 300_000), ("c3", 300_000)]             # .fai order, not the file's
    assert g.chrom_source == "asm.chromosomes.txt"
    note = next(n for n in g.notes if "asm.chromosomes.txt" in n)
    assert "2 of 5 sequences" in note and "other 3 sequences (310,000 bp)" in note


def test_the_notes_say_which_source_set_the_chromosomes(tmp_path):
    _assembly(tmp_path, "asm", ASM)
    g = inputs.discover(tmp_path)[0]
    assert g.chrom_source == "length rule"
    assert "chromosome set by the length rule" in g.notes
    assert "2 shorter sequences (unplaced scaffolds; 10,000 bp) are not analyzed" in g.notes


def test_the_notes_say_when_the_set_ends_at_a_gradual_tail(tmp_path):
    """Where the tail test cut 0.8's set, the notes (log and report) say that the set ends at a gradual tail, how many
    sequences were left out and the length of the shortest one kept, and how to bring back a chromosome cut with the
    tail. They do not call the sequences left out unplaced scaffolds: dot chromosomes can be among them. A clean gap
    (the 0.8 set) and a --classify-from target (0.8's rule) keep the 0.8 wording."""
    _assembly(tmp_path, "asm", GRADUAL)
    for name in ("clean", "T"):
        (tmp_path / name).mkdir()
    _assembly(tmp_path / "clean", "asm", ASM)
    _assembly(tmp_path / "T", "T", GRADUAL)
    g = inputs.discover(tmp_path)[0]
    assert len(g.chroms) == 25 and g.notes[0] == "chromosome set by the length rule" and len(g.notes) == 2
    assert g.notes[1] == ("the set ends at a gradual tail of short sequences: 644 shorter sequences "
                          f"({sum(GRADUAL.values()) - sum(g.chroms.values()):,} bp) are not analyzed, and the shortest "
                          "one kept is 7,680,000 bp; if chromosomes were left out, list every chromosome in "
                          "asm.chromosomes.txt")
    assert inputs.discover(tmp_path / "clean")[0].notes == [
        "chromosome set by the length rule", "2 shorter sequences (unplaced scaffolds; 10,000 bp) are not analyzed"]
    assert inputs.discover(tmp_path / "T", legacy=frozenset({"T"}))[0].notes == [
        "chromosome set by the length rule of cenfuse 0.8 (--classify-from target)"]


@pytest.mark.parametrize("text, msg", [
    (b"c1\nnope\n", "'nope'"),                                       # not a sequence of the .fai
    (b"c1\nc2\nc1\n", "more than once: 'c1'"),
    (b"# only comments\n\n   \n", "no sequence names"),
    (b"", "no sequence names"),
    (b"a\nb\nc\nd\ne\nf\ng\n", "'a', 'b', 'c', 'd', 'e' and 2 more"),
    ("c1\n".encode("utf-16"), "not UTF-8 text"),                     # PowerShell's default encoding
])
def test_chromosome_list_file_errors(tmp_path, text, msg):
    _assembly(tmp_path, "asm", ASM)
    (tmp_path / "asm.chromosomes.txt").write_bytes(text)
    with pytest.raises(CenfuseError, match=msg) as e:
        inputs.discover(tmp_path)
    assert "asm.chromosomes.txt" in str(e.value)


def test_chromosome_list_file_accepts_hand_written_files(tmp_path):
    fai = {"HG002#1#chr1": 10, "c2": 20, "c3": 30}
    f = tmp_path / "x.chromosomes.txt"
    # a byte-order mark and CRLF (Notepad), a lone CR, blanks and tabs around a name, an indented comment, a '#' inside
    # a name (pangenome naming), no newline at the end
    f.write_bytes("\ufeff# header\r\n  HG002#1#chr1\t \r\n\r   # skip c2\r\nc3".encode())
    assert inputs.read_chromosome_list(f, fai) == {"HG002#1#chr1": 10, "c3": 30}


def test_a_chromosome_list_for_no_genome_is_reported(tmp_path, monkeypatch, caplog):
    _assembly(tmp_path, "asm", ASM)
    (tmp_path / "asn.chromosomes.txt").write_text("c1\n")                          # a typo in the genome's name
    monkeypatch.setattr(log, "propagate", True)
    with caplog.at_level(logging.WARNING, logger="cenfuse"):
        g = inputs.discover(tmp_path)[0]
    assert g.chrom_source == "length rule"
    assert any("asn.chromosomes.txt" in r.getMessage() and "no genome" in r.getMessage() for r in caplog.records)


def test_classify_from_targets_keep_the_08_rule(tmp_path):
    for name in ("S", "T"):
        _assembly(tmp_path, name, GRADUAL)
    s, t = inputs.discover(tmp_path, legacy=frozenset({"T"}))
    assert (len(s.chroms), s.chrom_source) == (25, "length rule")
    assert (len(t.chroms), t.chrom_source) == (669, "length rule of cenfuse 0.8 (--classify-from target)")
    assert "chromosome set by the length rule of cenfuse 0.8 (--classify-from target)" in t.notes
    assert len(inputs.discover(tmp_path)[1].chroms) == 25                           # without legacy: the 0.9 set
    assert len(inputs.chromosome_set(GRADUAL)) == 25 and len(inputs.chromosome_set(GRADUAL, legacy=True)) == 669


def test_a_chromosome_list_also_overrides_a_classify_from_target(tmp_path):
    _assembly(tmp_path, "T", GRADUAL)
    (tmp_path / "T.chromosomes.txt").write_text("s0\ns1\n")
    t = inputs.discover(tmp_path, legacy=frozenset({"T"}))[0]
    assert list(t.chroms) == ["s0", "s1"] and t.chrom_source == "T.chromosomes.txt"


def test_run_gives_classify_from_targets_the_08_rule(tmp_path):
    for name in ("S", "T"):
        _assembly(tmp_path, name, GRADUAL)
    run = Run(tmp_path, tmp_path / "out", classify_from="T=S")
    assert {g.name: (len(g.chroms), g.chrom_source) for g in run.genomes} == {
        "S": (25, "length rule"), "T": (669, "length rule of cenfuse 0.8 (--classify-from target)")}
    assert {g.name: len(g.chroms) for g in Run(tmp_path, tmp_path / "out").genomes} == {"S": 25, "T": 25}


def test_run_reports_a_malformed_classify_from_as_a_cenfuse_error(tmp_path):
    _assembly(tmp_path, "T", ASM)
    with pytest.raises(CenfuseError, match="--classify-from"):
        Run(tmp_path, tmp_path / "out", classify_from="T")


def test_methylation_D_counts_the_chromosomes_cenfuse_analyzes(tmp_path, monkeypatch, caplog):
    """10 chromosomes of 25 Mb and 590 scaffolds of 1 kb: the rule counts 10 (D = 270 kb); a list naming all 600
    sequences makes the count 600 and D 5 kb, under two 50-kb bins, which stops the methylation line before it reads a
    sample (the genome has none)."""
    lengths = {**{f"chr{i}": 25_000_000 for i in range(10)}, **{f"s{i}": 1_000 for i in range(590)}}
    _assembly(tmp_path, "asm", lengths)
    (tmp_path / "asm.chromosomes.txt").write_text("\n".join(lengths) + "\n")
    g = inputs.discover(tmp_path)[0]
    assert len(g.chroms) == 600
    monkeypatch.setattr(log, "propagate", True)
    with caplog.at_level(logging.WARNING, logger="cenfuse"):
        assert methfuse.prepare(g, tmp_path / "arrays.gff", tmp_path / "cache", 1) is None
    assert "expected centromere size D = 5 kb is under two 50-kb bins" in caplog.text


@pytest.mark.parametrize("name", ["g.allValidPairs", "g.allValidPairs.gz", "g.merged_nodups.txt",
                                  "g.merged_nodups.txt.gz"])
def test_hicpro_and_juicer_files_are_libraries(tmp_path, name):
    _genome(tmp_path)
    (tmp_path / "hic").mkdir()
    (tmp_path / "hic" / name).write_text("")
    g = inputs.discover(tmp_path)[0]
    assert [(x.name, x.path.name) for x in g.hic] == [("hic", name)]



def test_pairs_are_preferred_over_a_juicer_file_beside_them(tmp_path, caplog):
    _genome(tmp_path)
    (tmp_path / "hic").mkdir()
    for n in ("g.pairs.gz", "g.merged_nodups.txt.gz"):
        (tmp_path / "hic" / n).write_text("")
    g = inputs.discover(tmp_path)[0]
    assert [x.path.name for x in g.hic] == ["g.pairs.gz"]
    assert "g.merged_nodups.txt.gz not used" in caplog.text


@pytest.mark.parametrize("reader", ["cenfuse", "satellite engine"])
def test_a_plain_gzip_fasta_gets_a_recipe_samtools_accepts(tmp_path, reader):
    """A plain-gzip genome stops at once with a recipe that works: `samtools faidx` refuses plain gzip ('Cannot index files
    compressed with gzip, please use bgzip'), so the recipe recompresses with bgzip first. It is run here as given."""
    import gzip
    import re
    import shutil
    import subprocess
    from cenfuse.sat import genomes
    from cenfuse.sat.util import SatelliteError
    f = tmp_path / "g.fa.gz"
    with gzip.open(f, "wt") as h:
        h.write(">c1\nACGTACGTAC\nGGCC\n>c2\nTTTTAAAA\n")
    read = (lambda: inputs.discover(tmp_path)) if reader == "cenfuse" else (lambda: genomes.discover_genomes(tmp_path))
    with pytest.raises(CenfuseError if reader == "cenfuse" else SatelliteError) as e:
        read()
    recipe = re.search(r"recompress and index: (.+?)\)", str(e.value)).group(1)
    assert recipe == f"gzip -dc {f} | bgzip > {f}.tmp && mv {f}.tmp {f} && samtools faidx {f}"
    if not (shutil.which("bgzip") and shutil.which("samtools")):
        pytest.skip("bgzip or samtools not on PATH")
    subprocess.run(recipe, shell=True, check=True, cwd=tmp_path, capture_output=True)
    g = read()[0]
    assert dict(g.chroms) == {"c1": 14, "c2": 8}


@pytest.mark.parametrize("sample", ["NA", "None", "nan", "null", "NULL"])
def test_a_hic_sample_named_like_a_missing_value_keeps_its_name(tmp_path, sample):
    """pandas reads NA, None, nan, null and NULL as missing by default: the samples table turned such a plant's name into '' and
    stopped with "bad sample name ''"."""
    d = tmp_path / "hic"
    d.mkdir()
    (d / "lib.pairs.gz").write_text("")
    (d / "g.tsv").write_text(f"sample\tfile\n{sample}\tlib.pairs.gz\n")
    assert [(x.name, x.path.name) for x in inputs.hic_libraries(d, "g")] == [(sample, "lib.pairs.gz")]
