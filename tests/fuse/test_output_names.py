"""Output files keep chromosome names as given, including names that hold ':' (locus IDs are 'genome:chrom:Ln'); GFF3 column 1
percent-encodes what the GFF3 spec does not allow in a seqid."""
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


def test_gff3_seqids_escape_what_the_spec_does_not_allow(tmp_path):
    """GFF3 column 1 and ##sequence-region percent-encode every character outside [a-zA-Z0-9.:^*$@!+_?-|], as the GFF3 spec
    requires of a seqid; the TSV and BED keep the name as given (0.16.0 wrote 'chr1|x=1;y', 'chr#3%2C' and 'chr,4*' raw)."""
    from urllib.parse import unquote
    names = {"chr1|x=1;y": "chr1|x%3D1%3By", "chr#3%2C": "chr%233%252C", "chr,4*": "chr%2C4*", "chr 5~": "chr%205%7E",
             "Chr1.2:A_b-c": "Chr1.2:A_b-c"}
    t = pd.DataFrame([dict(chrom=c, chrom_len=1000, start="11", end="20", relation="sequence_only", confidence="high",
                           call_source="sequence", support="satellite", against=".", flags=".", rule="sequence",
                           seq_start="11", seq_end="20", seq_evidence="satellite", main_family="f1", offset_bp=".",
                           cpg_kt_start="13", cpg_kt_end="15", cpg_dip_depth=".", cpg_dip_z=".", cpg_pos_p=".",
                           cpg_kt_meth=".", cpg_array_meth=".") for c in names])
    output.write_calls(t, "G", {c: 1000 for c in names}, tmp_path)
    gff = (tmp_path / "G.centromeres.gff3").read_text().splitlines()
    rows = [x.split("\t") for x in gff if not x.startswith("#")]
    assert [r[0] for r in rows] == [e for e in names.values() for _ in range(2)]        # the call and its kinetochore
    assert [x.split()[1] for x in gff if x.startswith("##sequence-region")] == list(names.values())
    assert [unquote(r[0]) for r in rows[::2]] == list(names)
    bed = [x.split("\t")[0] for x in (tmp_path / "G.centromeres.bed").read_text().splitlines() if not x.startswith("#")]
    assert bed == list(names)
    lt = pd.DataFrame([dict({c: "." for c in LOCUS_COLS}, locus="G:chr1|x=1;y:L1", chrom="chr1|x=1;y", start=11, end=20,
                            length=10, kind="telomere", lines="telomere", n_lines=1)], columns=LOCUS_COLS)
    feats = [("G:chr1|x=1;y:L1", (10, 20, "telomere_array", "telomere", "AACCCT", dict(side="start", bp=10)))]
    output.write_loci(lt, feats, "G", {"chr1|x=1;y": 1000}, tmp_path)
    rows = [x.split("\t") for x in (tmp_path / "G.features.gff3").read_text().splitlines() if not x.startswith("#")]
    assert [r[0] for r in rows] == ["chr1|x%3D1%3By", "chr1|x%3D1%3By"]
