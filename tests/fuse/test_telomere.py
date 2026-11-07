"""Telomere rules for the ends the strict rule misses: stretches under half unit copies that reach the sequence end
(soybean chr3: 48 kb, 20.2 kb of copies) and units of 13-29 bp (C. tropicalis: 23 bp)."""
import random

import pysam

from cenfuse import telomere as tm

UNIT = "TTTAGGG"
LONG_UNIT = "ACGGATGTCACGATCATTGGTCA"        # 23 bp, primitive, as in Candida


def bg(rnd, n: int) -> str:
    return "".join(rnd.choice("ACGT") for _ in range(n))


def loose(rnd, unit: str, copies: int, spacer: int, blocks: int) -> str:
    """A stretch of blocks of `copies` units, each followed by `spacer` random bases: density copies*k / (copies*k + spacer)."""
    return "".join(unit * copies + bg(rnd, spacer) for _ in range(blocks))


def write(tmp_path, seqs: dict):
    fa = tmp_path / "g.fa"
    fa.write_text("".join(f">{c}\n{s}\n" for c, s in seqs.items()))
    pysam.faidx(str(fa))
    return fa, {c: len(s) for c, s in seqs.items()}


def side_of(rows, chrom, side):
    return next((r for r in rows if r[0] == chrom and r[1] == side), None)


def genome(tmp_path, rnd, left_end: str, chrom_len: int = 700_000):
    """Three chromosomes with strict telomere arrays at most ends; c1's start holds `left_end`."""
    tel = UNIT * 200
    rc = tel.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    seqs = {"c1": left_end + bg(rnd, chrom_len) + tel,
            "c2": rc + bg(rnd, chrom_len) + tel,
            "c3": rc + bg(rnd, chrom_len) + tel}
    return write(tmp_path, seqs)


def test_a_loose_stretch_reaching_the_end_is_an_array(tmp_path):
    rnd = random.Random(1)
    stretch = loose(rnd, UNIT, 10, 140, 200)                  # 42 kb at one third density, as soybean chr3
    assert tm.end_array(stretch + bg(rnd, 60_000), tm.canonical(UNIT), "start") is None   # the strict rule misses it
    fa, chroms = genome(tmp_path, rnd, stretch)
    r = tm.detect(fa, chroms)
    assert r["unit"] == tm.canonical(UNIT)
    row = side_of(r["ends"], "c1", "start")
    assert row is not None and row[2] == 0 and row[4] >= 10_000
    assert len(r["ends"]) == 6


def test_a_loose_stretch_away_from_the_end_is_not_an_array(tmp_path):
    rnd = random.Random(2)
    stretch = bg(rnd, 5_000) + loose(rnd, UNIT, 10, 140, 200)
    fa, chroms = genome(tmp_path, rnd, stretch)
    r = tm.detect(fa, chroms)
    assert side_of(r["ends"], "c1", "start") is None and len(r["ends"]) == 5


def test_scattered_copies_below_the_floor_are_not_an_array(tmp_path):
    rnd = random.Random(3)
    stretch = loose(rnd, UNIT, 5, 980, 90)                    # copies every kb: 3.4% of 91 kb, chained by CHAIN_GAP
    assert tm.chained(stretch, tm.canonical(UNIT))[0][2] >= tm.MIN_ARRAY
    fa, chroms = genome(tmp_path, rnd, stretch)
    r = tm.detect(fa, chroms)
    assert side_of(r["ends"], "c1", "start") is None


def test_a_strict_array_is_kept_when_a_loose_stretch_holds_more_copies(tmp_path):
    """The exemption only fills ends without an array: an end the strict rule answers does not change."""
    rnd = random.Random(4)
    strict = UNIT * 100
    seq = strict + bg(rnd, 2_000) + loose(rnd, UNIT, 10, 140, 200)
    fa, chroms = genome(tmp_path, rnd, seq)
    r = tm.detect(fa, chroms)
    row = side_of(r["ends"], "c1", "start")
    assert row[2] == 0 and row[3] < len(strict) + 50


def test_loose_stretches_do_not_choose_the_unit(tmp_path):
    """A decoy unit in loose stretches at every end covers more bp than the real unit but has no strict array."""
    rnd = random.Random(5)
    decoy = "TTGGACC"
    tel = UNIT * 60
    seqs = {c: loose(rnd, decoy, 10, 140, 60) + tel + bg(rnd, 700_000) + tel + loose(rnd, decoy, 10, 140, 60)
            for c in ("c1", "c2", "c3")}
    fa, chroms = write(tmp_path, seqs)
    assert tm.detect(fa, chroms)["unit"] == tm.canonical(UNIT)


def long_unit_genome(tmp_path, rnd, ends_with: int, n: int = 4):
    tel = LONG_UNIT * 30
    seqs, left = {}, ends_with
    for i in range(n):
        a = tel if left > 0 else ""
        left -= 1
        b = tel if left > 0 else ""
        left -= 1
        seqs[f"c{i}"] = a + bg(rnd, 700_000) + b
    return write(tmp_path, seqs)


def test_a_long_unit_is_found_where_no_short_unit_qualifies(tmp_path):
    fa, chroms = long_unit_genome(tmp_path, random.Random(6), 7)    # 7 of 8 ends, as C. tropicalis's 13 of 14
    r = tm.detect(fa, chroms)
    assert r["unit"] == tm.canonical(LONG_UNIT) and len(r["unit"]) == 23
    assert len(r["ends"]) == 7


def test_a_long_unit_at_a_few_ends_is_not_a_telomere(tmp_path):
    """A subtelomeric minisatellite at two of eight ends: it is enriched there, but a telomere reaches most ends."""
    fa, chroms = long_unit_genome(tmp_path, random.Random(7), 2)
    assert tm.detect(fa, chroms)["unit"] is None


def test_a_short_unit_outranks_a_long_one_that_covers_more_bp(tmp_path):
    rnd = random.Random(8)
    tel, sub = UNIT * 60, LONG_UNIT * 60
    seqs = {c: sub + tel + bg(rnd, 700_000) + tel + sub for c in ("c1", "c2", "c3")}
    fa, chroms = write(tmp_path, seqs)
    assert tm.detect(fa, chroms)["unit"] == tm.canonical(UNIT)


def test_left_out_sequences_with_telomeres_at_both_ends_are_flagged(tmp_path):
    """future.md #8: a complete molecule the length rule cut (a dot chromosome) is named; a scaffold with one telomere is not."""
    import random
    from cenfuse import inputs
    rng = random.Random(5)
    rnd = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    tel = "TTTAGGG" * 300
    rev = "CCCTAAA" * 300
    seqs = {"chr1": tel + rnd(400_000) + rev, "dot": tel + rnd(250_000) + rev, "scaf": tel + rnd(250_000), "tiny": tel + rnd(1_000) + rev}
    fa = tmp_path / "g.fa"
    fa.write_text("".join(f">{c}\n{s}\n" for c, s in seqs.items()))
    import pysam
    pysam.faidx(str(fa))
    fai = {c: len(s) for c, s in seqs.items()}
    assert inputs.left_out(fa, fai, {"chr1": fai["chr1"]}, "AAACCCT") == [("dot", fai["dot"])]
    assert inputs.left_out(fa, fai, {"chr1": fai["chr1"]}, None) == []
    assert inputs.left_out(fa, fai, fai, "AAACCCT") == []
