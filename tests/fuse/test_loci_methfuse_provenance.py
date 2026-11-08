"""Unit tests of loci.py, methfuse.py and provenance.py (future.md #65): until 0.14 they were covered end to end only."""
import numpy as np

from cenfuse import loci, methfuse, provenance


def test_runs_of_a_mask():
    assert loci._runs(np.array([0, 1, 1, 0, 1, 0, 0, 1], dtype=bool)) == [(1, 2), (4, 4), (7, 7)]
    assert loci._runs(np.zeros(4, dtype=bool)) == []


def test_members_merge_into_loci_with_hic_padded_one_bin():
    m = [(100, 200, "satellite_cluster", "satellite", "x", {"role": "centromeric"}),
         (190, 300, "te_core", "te", "CRM", {}),
         (520, 600, "hic_anchored_core", "hic", "peak", {}),             # 220 bp away: joins only with the 250-bp pad
         (2000, 2100, "telomere_array", "telomere", "TTTAGGG", {})]
    groups = loci.merge(m, tol=250)
    assert [len(g) for g in groups] == [3, 1]
    assert [len(g) for g in loci.merge(m, tol=0)] == [2, 1, 1]


def test_the_kind_of_a_locus_follows_its_strongest_member():
    sat = lambda role, rdna=".": (0, 1, "satellite_cluster", "satellite", "f", {"role": role, "rdna": rdna})
    assert loci.kind([sat("centromeric"), (0, 1, "integrated_call", "fused", "x", {})]) == "centromere"
    assert loci.kind([(0, 1, "telomere_array", "telomere", "u", {})]) == "telomere"
    assert loci.kind([sat("other", rdna="45S")]) == "rDNA"
    assert loci.kind([sat("subtelomeric")]) == "subtelomeric_block"
    assert loci.kind([(0, 1, "te_core", "te", "CRM", {"marker": True})]) == "alt_centromere_candidate"
    assert loci.kind([(0, 1, "te_core", "te", "Tekay", {})]) == "te_core"
    assert loci.kind([sat("other")]) == "satellite_block"


def test_methylation_sites_merge_within_half_a_window():
    raw = [("sequence_call", 1000, 3000), ("hic_peak", 2500, 2700), ("te_core", 50_000, 52_000), ("provisional", -1, -1)]
    s = methfuse.merge_sites(raw, W=10_000, L=60_000)
    assert [x["origins"] for x in s] == [["sequence_call", "hic_peak"], ["te_core"]]
    assert s[0]["mid"] == 2000 and (s[0]["start"], s[0]["end"]) == (0, 7000)
    assert s[1]["end"] == 56_000


def test_chg_sites_and_motif_checks():
    seq = "AACAGTTCTGAACGTT"                    # CAG at 2, CTG at 7
    assert list(provenance.symmetric_chg(seq, 100)) == [102, 107]
    assert len(provenance.symmetric_chg("CA", 0)) == 0
    pos = np.array([102, 104, 110])            # + strand CAG C; the G of CAG on - strand; CG at 12 (C at 112 is the CG)
    strand = np.array(["+", "-", "+"])
    hits = provenance.on_motif("NN" + seq + "NN", 98, pos, strand, "CHG")
    assert list(hits) == [True, True, False]
    assert list(provenance.on_motif("NN" + seq + "NN", 98, np.array([112]), np.array(["+"]), "CG")) == [True]
