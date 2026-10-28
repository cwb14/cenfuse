"""GENOME_DIR discovery (inputs.py): accepted annotation names, and Hi-C formats cenfuse cannot read failing clearly
instead of being ignored."""
import pytest

from cenfuse import inputs
from cenfuse.util import CenfuseError


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
                                        ("g.merged_nodups.txt.gz", "pairs")])
def test_unreadable_hic_formats_fail_with_a_recipe(tmp_path, name, hint):
    _genome(tmp_path)
    (tmp_path / "hic").mkdir()
    (tmp_path / "hic" / name).write_text("")
    with pytest.raises(CenfuseError, match=hint):
        inputs.discover(tmp_path)
