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
    assert [r[0] for r in rows] == ["chr1:A", "chr1:A"]
    assert rows[0][2] == "telomere" and "kind=telomere" in rows[0][8]
    assert rows[1][2] == "telomeric_repeat" and "member=telomere_array" in rows[1][8]
    assert "locus=G:chr1:A:L1" in rows[1][8] and "Parent=" not in rows[1][8]


def test_every_type_is_a_sequence_ontology_term_and_the_kind_is_kept(tmp_path):
    """GFF3 types are SO terms; cenfuse's own kind and member type stay as attributes."""
    from cenfuse.loci import KIND_ORDER
    so = set(output.SO_KIND.values()) | set(output.SO_MEMBER.values())
    assert so <= {"centromere", "telomere", "telomeric_repeat", "tandem_repeat", "satellite_DNA", "repeat_region",
                  "epigenetically_modified_region", "biological_region"}
    assert set(KIND_ORDER) <= set(output.SO_KIND)
    members = ("rdna_array", "rdna5s_array", "telomere_array", "telomere_interstitial", "integrated_call",
               "sequence_call", "satellite_cluster", "te_core", "hic_anchored_core", "hic_denovo_core", "hic_hub",
               "meth_candidate")
    assert set(members) <= set(output.SO_MEMBER)
