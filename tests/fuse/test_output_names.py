"""Output files keep chromosome names as given, including names that hold ':' (locus IDs are 'genome:chrom:Ln')."""
import pandas as pd

from cenfuse import output
from cenfuse.loci import LOCUS_COLS


def test_member_lines_keep_a_chromosome_name_with_a_colon(tmp_path):
    lt = pd.DataFrame([dict({c: "." for c in LOCUS_COLS}, locus="G:chr1:A:L1", chrom="chr1:A", start=11, end=20,
                            length=10, kind="telomere", lines="telomere", n_lines=1)], columns=LOCUS_COLS)
    feats = [("G:chr1:A:L1", (10, 20, "telomere_array", "telomere", "AACCCT", dict(side="start", bp=10)))]
    output.write_loci(lt, feats, "G", {"chr1:A": 1000}, tmp_path)
    rows = [l.split("\t") for l in (tmp_path / "G.features.gff3").read_text().splitlines() if not l.startswith("#")]
    assert [r[0] for r in rows] == ["chr1:A", "chr1:A"] and rows[1][2] == "telomere_array"
