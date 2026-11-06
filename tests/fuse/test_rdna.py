"""45S arrays from the universal 18S motif (rdna.py) and their use as NOR flags on Hi-C peaks."""
import random

import numpy as np
import pandas as pd
import pysam
import pytest

from cenfuse import hicfuse, rdna


RC = str.maketrans("ACGT", "TGCA")
NUCLEAR = rdna.MOTIF + "GCTCCTACCGATTGAATGGT"                # nuclear 18S, as in most plants
PLASTID = rdna.MOTIF + "ACACTATGGGAGCTGGCCAT"                # plastid 16S (rice chr4 9.23-9.37 Mb)
BACTERIAL = rdna.MOTIF + "ACACCATGGGAGTTGGTTTT"              # bacterial 16S, plant mitochondrial 18S
NUMT = rdna.MOTIF + "ACCCTCCTCAAGTATACTTC"                   # animal mitochondrial 12S
MUTATED = rdna.MOTIF + "ACTCCTACCGATTGAATGGT"                # nuclear, G -> A at the first base (barley 7H: 3 of 4 hits)


def revcomp(s: str) -> str:
    return s.translate(RC)[::-1]


def faidx(tmp_path, seq: str, name: str = "g.fa"):
    fa = tmp_path / name
    fa.write_text(">c\n" + seq + "\n")
    pysam.faidx(str(fa))
    return fa


def bases(seed: int):
    rnd = random.Random(seed)
    return lambda n: "".join(rnd.choice("ACGT") for _ in range(n))


def test_tandem_motif_runs_are_nors_and_single_hits_are_not(tmp_path):
    bg = bases(1)
    unit = bg(9_000)
    seq = bg(200_000) + (unit[:4000] + NUCLEAR + unit[4000:]) * 5 + bg(300_000) + revcomp(NUCLEAR) + bg(200_000)
    fa = faidx(tmp_path, seq)
    pos = rdna.hits(fa, "c", len(seq))
    assert len(pos) == 6                                     # 5 tandem copies + 1 isolated reverse-strand copy
    nor = rdna.detect(fa, {"c": len(seq)})
    assert len(nor) == 1 and nor[0][3] == 5 and nor[0][4] == 9_036
    assert nor[0][1] < 200_000 + 4000 < 200_000 + 4 * 9036 + 4000 < nor[0][2]


@pytest.mark.parametrize("organelle", [PLASTID, BACTERIAL, NUMT])
def test_organelle_motif_hits_make_no_array(tmp_path, organelle):
    """Rice chr4 9.23-9.37 Mb, reported as a NOR, is a chloroplast insertion. Organelle SSU rRNA has the motif followed by
    A where nuclear 18S has G, and differs from nuclear 18S at most of the bases after it, on either strand."""
    bg = bases(2)
    spacer = bg(9_000)
    seq = bg(50_000) + (organelle + spacer + revcomp(organelle) + spacer) * 3 + bg(50_000)
    fa = faidx(tmp_path, seq)
    assert len(rdna.hits(fa, "c", len(seq))) == 0
    assert rdna.detect(fa, {"c": len(seq)}) == []


def test_nuclear_copies_with_a_mutated_first_base_stay_an_array(tmp_path):
    """Barley 7H, 10 5.8S genes in 120 kb: 3 of its 4 motif hits have A at that base, and the next 19 bases are nuclear."""
    bg = bases(3)
    spacer = bg(9_000)
    seq = bg(50_000) + (MUTATED + spacer) * 2 + revcomp(MUTATED) + spacer + NUCLEAR + bg(50_000)
    fa = faidx(tmp_path, seq)
    nor = rdna.detect(fa, {"c": len(seq)})
    assert len(nor) == 1 and nor[0][3] == 4


def test_a_nor_keeps_its_nuclear_units_and_loses_an_organelle_hit_inside_it(tmp_path):
    bg = bases(4)
    spacer = bg(9_000)
    seq = bg(50_000) + (NUCLEAR + spacer) * 2 + PLASTID + spacer + (NUCLEAR + spacer) * 2 + bg(50_000)
    fa = faidx(tmp_path, seq)
    nor = rdna.detect(fa, {"c": len(seq)})
    assert len(nor) == 1 and nor[0][3] == 4


@pytest.mark.parametrize("after", ["G", "C", "T", "N"])
def test_only_a_following_a_can_mark_an_organelle_hit(tmp_path, after):
    """Nuclear 18S of some lineages differs from the usual context; another first base is no evidence for an organelle,
    however different the bases after it."""
    bg = bases(5)
    seq = bg(100) + rdna.MOTIF + after + "ACTATGGGAGCTGGCCAT" + bg(100)
    assert len(rdna.hits(faidx(tmp_path, seq), "c", len(seq))) == 1


def test_motif_at_the_sequence_ends_with_too_little_context_is_kept(tmp_path):
    """A hit with no or too few following bases on the sequence has no evidence against it, on either strand."""
    body = bases(6)(500)
    seq = revcomp(rdna.MOTIF) + body + rdna.MOTIF + "A"
    assert list(rdna.hits(faidx(tmp_path, seq), "c", len(seq))) == [0, 500 + len(rdna.MOTIF)]


def test_a_truncated_context_of_most_of_its_length_is_still_judged(tmp_path):
    seq = bases(7)(500) + PLASTID[:len(rdna.MOTIF) + 15]
    assert len(rdna.hits(faidx(tmp_path, seq), "c", len(seq))) == 0


@pytest.mark.parametrize("gene,keep", [(PLASTID, 0), (NUCLEAR, 1), (MUTATED, 1)])
def test_context_is_read_across_a_chunk_edge(tmp_path, monkeypatch, gene, keep):
    """Forward strand: the bases after the motif sit in the next chunk."""
    monkeypatch.setattr(rdna, "CHUNK", 1_000)
    bg = bases(8)
    seq = bg(1_000 - len(rdna.MOTIF)) + gene + bg(500)
    assert len(rdna.hits(faidx(tmp_path, seq), "c", len(seq))) == keep


@pytest.mark.parametrize("gene,keep", [(PLASTID, 0), (NUCLEAR, 1), (MUTATED, 1)])
def test_context_of_a_reverse_hit_is_read_across_a_chunk_edge(tmp_path, monkeypatch, gene, keep):
    """Reverse strand: the bases after the motif, read along the gene, lie before the hit, in the previous chunk."""
    monkeypatch.setattr(rdna, "CHUNK", 1_000)
    bg = bases(9)
    extra = len(gene) - len(rdna.MOTIF)
    seq = bg(1_005 - extra) + revcomp(gene) + bg(500)
    assert len(rdna.hits(faidx(tmp_path, seq), "c", len(seq))) == keep


def test_a_hit_on_a_chunk_edge_is_counted_once(tmp_path, monkeypatch):
    monkeypatch.setattr(rdna, "CHUNK", 1_000)
    bg = bases(10)
    seq = bg(995) + NUCLEAR + bg(500)
    assert list(rdna.hits(faidx(tmp_path, seq), "c", len(seq))) == [995]


def test_peak_beside_an_unmappable_nor_cannot_oppose():
    t = pd.DataFrame([dict(chrom="c", bin_size=250_000, anch_peak_start=4_500_000, anch_peak_end=4_750_000,
                           anch_strong=True, anch_level="high", anch_flags="", seq_start=9_000_000,
                           seq_end=12_000_000, masked=((2_000_000, 4_500_000),))])
    nors = [("c", 2_100_000, 3_900_000, 20, 44_300)]        # 2.4 blind bins between the NOR and the peak
    out = hicfuse.nor_peaks(t, nors)
    assert out.anch_flags.iat[0] == "rDNA" and not out.anch_strong.iat[0]
    far = hicfuse.nor_peaks(t, [("c", 500_000, 900_000, 9, 44_300)])     # mappable bins between: not this NOR
    assert far.anch_strong.iat[0]
    merged = hicfuse.add_nors({"c": (np.array([1]), np.array([2]))}, nors)
    assert list(merged["c"][0]) == [1, 2_100_000]


def test_5s_arrays_are_tandem_runs_of_whole_genes(tmp_path):
    """A tandem run of whole 5S genes is an array; a lone gene and a tandem run of half genes (pseudogene / SINE
    fragments) are not."""
    rnd = random.Random(2)
    bg = lambda n: "".join(rnd.choice("ACGT") for _ in range(n))
    gene = rdna.five_s_hmm().consensus.upper()
    assert len(gene) == 120
    rc = gene.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    spacer = bg(380)
    seq = (bg(300_000) + (gene + spacer) * 6 + bg(400_000) + rc + bg(400_000) + (gene[:60] + bg(440)) * 6
           + bg(300_000))
    fa = tmp_path / "g.fa"
    fa.write_text(">c\n" + seq + "\n")
    pysam.faidx(str(fa))
    pos = rdna.genes_5s(fa, {"c": len(seq)}, 2)
    assert len(pos["c"]) == 7                                # 6 tandem + 1 lone reverse-strand gene; halves not
    arr = rdna.detect_5s(fa, {"c": len(seq)}, 2)
    assert len(arr) == 1 and arr[0][3] == 6 and arr[0][4] == 500
    assert arr[0][1] <= 300_000 < 300_000 + 6 * 500 <= arr[0][2]
