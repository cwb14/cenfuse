"""45S arrays from the universal 18S motif (rdna.py) and their use as NOR flags on Hi-C peaks."""
import random

import numpy as np
import pandas as pd
import pysam

from cenfuse import hicfuse, rdna


def test_tandem_motif_runs_are_nors_and_single_hits_are_not(tmp_path):
    rnd = random.Random(1)
    bg = lambda n: "".join(rnd.choice("ACGT") for _ in range(n))
    unit = bg(9_000)
    rc = rdna.MOTIF.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    seq = bg(200_000) + (unit[:4000] + rdna.MOTIF + unit[4000:]) * 5 + bg(300_000) + rc + bg(200_000)
    fa = tmp_path / "g.fa"
    fa.write_text(">c\n" + seq + "\n")
    pysam.faidx(str(fa))
    pos = rdna.hits(fa, "c", len(seq))
    assert len(pos) == 6                                     # 5 tandem copies + 1 isolated reverse-strand copy
    nor = rdna.detect(fa, {"c": len(seq)})
    assert len(nor) == 1 and nor[0][3] == 5 and nor[0][4] == 9_016
    assert nor[0][1] < 200_000 + 4000 < 200_000 + 4 * 9016 + 4000 < nor[0][2]


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
