"""Long models with a low-complexity tract search with those positions neutral (background emissions). Such a tract matches
every low-complexity stretch of the genome, so the model's filters pass nearly everything and each window costs model x
window DP (Lespedeza: 25-50 s per 2-Mb chunk with no hit). Neutral, the tract triggers nothing. A copy usually aligns through
it; where it does not, its pieces are joined, and a copy is whole by the positions that can score."""
from types import SimpleNamespace

import numpy as np
import pysam
import pytest

import cenfuse.sat.annotate as annotate
from cenfuse.sat.families import single_hmm
from cenfuse.sat.params import FULL_MONOMER_COV
from sat_common import dna

TRACT = "AT" * 60                                                    # 120 bp of low complexity


def mutate(rng, s: str, rate: float) -> str:
    b = list(s)
    for i in np.flatnonzero(rng.random(len(b)) < rate):
        b[i] = "ACGT"[(("ACGT".index(b[i])) + int(rng.integers(1, 4))) % 4]
    return "".join(b)


def write(tmp_path, seq: str, name: str = "c"):
    fa = tmp_path / "g.fa"
    fa.write_text(f">{name}\n{seq}\n")
    pysam.faidx(str(fa))
    return SimpleNamespace(name="g", fasta=fa, chroms={name: len(seq)})


def search(hmm, g, plain: bool = False, threads: int = 2, **kw):
    """search_genome's hits for one model; plain = with every model too short to be neutralized."""
    c = next(iter(g.chroms))
    if plain:
        with pytest.MonkeyPatch.context() as m:
            m.setattr(annotate, "TRACT_MIN_MODEL", 10 ** 9)
            return annotate.search_genome([hmm], g, threads, 1e-5, chunk=400_000, overlap=10_000, **kw)[c]
    return annotate.search_genome([hmm], g, threads, 1e-5, chunk=400_000, overlap=10_000, **kw)[c]


def whole(hits, M: int):
    return hits[(hits.hmm_to - hits.hmm_from + 1) >= FULL_MONOMER_COV * M].sort_values("start")


def same_copies(a, b, tol: int = 5):
    """The same hits, in the same places (within tol bp at each end)."""
    assert len(a) == len(b), (len(a), len(b))
    assert (abs(a.start.to_numpy() - b.start.to_numpy()) <= tol).all()
    assert (abs(a.end.to_numpy() - b.end.to_numpy()) <= tol).all()


# ---- which models, which stretches

def test_dustmasker_finds_a_simple_repeat_and_little_in_random_sequence():
    rng = np.random.default_rng(1)
    a, b = dna(rng, 2_000), dna(rng, 2_000)
    iv = annotate.low_complexity({"x": a + TRACT + b, "y": dna(rng, 4_000)})
    s, e = max(iv["x"].tolist(), key=lambda t: t[1] - t[0])
    assert abs(s - 2_000) <= 3 and abs(e - 2_120) <= 3
    assert (iv["y"][:, 1] - iv["y"][:, 0]).sum() < 100           # random sequence: only chance stretches of a few bases


def _models(**seqs):
    return [single_hmm(k, v) for k, v in seqs.items()]


def test_only_long_models_with_a_moderate_tract_change():
    rng = np.random.default_rng(2)
    models = _models(ok=dna(rng, 400) + TRACT + dna(rng, 400),
                     short=dna(rng, 150) + TRACT + dna(rng, 150),                # under TRACT_MIN_MODEL: stays exact
                     no_tract=dna(rng, 1_700),
                     little=dna(rng, 800) + "AT" * 12 + dna(rng, 800),           # 24 bp: under one DUST window
                     mostly=dna(rng, 300) + "AT" * 700 + dna(rng, 300))          # a low-complexity satellite
    out, tracts = annotate.neutralize(models)
    assert list(tracts) == ["ok"] and [h.name for h in out] == [h.name for h in models]
    assert [o is h for o, h in zip(out, models)] == [False, True, True, True, True]


def test_no_long_model_runs_no_dustmasker(monkeypatch):
    monkeypatch.setattr(annotate, "low_complexity", lambda seqs: pytest.fail("dustmasker ran"))
    stand_in = SimpleNamespace(name="m", M=100)
    out, tracts = annotate.neutralize([stand_in])
    assert out == [stand_in] and tracts == {}


def test_tract_positions_emit_the_background_and_others_do_not_change():
    rng = np.random.default_rng(3)
    h = single_hmm("fam", dna(rng, 400) + TRACT + dna(rng, 400))
    before = np.array([list(h.match_emissions[k]) for k in range(1, h.M + 1)])
    (n,), tracts = annotate.neutralize([h])
    assert np.array_equal(before, np.array([list(h.match_emissions[k]) for k in range(1, h.M + 1)]))     # a copy
    after = np.array([list(n.match_emissions[k]) for k in range(1, n.M + 1)])
    assert np.allclose(after[410:510], 0.25, atol=0.02)          # inside the tract
    assert np.array_equal(after[:380], before[:380]) and np.array_equal(after[-380:], before[-380:])
    assert n.consensus == h.consensus and n.M == h.M and n.name == h.name
    assert any(abs(a - 400) <= 3 and abs(b - 520) <= 3 for a, b in tracts["fam"].tolist())


def row(start, end, strand, hmm_from, hmm_to, M=1_000, bits=500.0, evalue=1e-30, identity=0.9):
    return ("c", start, end, strand, "f", bits, evalue, hmm_from, hmm_to, M, identity)


TRACT_400 = np.array([[400, 500]])                               # model positions 401-500 of 1,000: 900 can score


def span(r):
    return r[7:9]


# ---- a hit is whole by the positions that can score

def test_a_hit_over_90pct_of_the_scoring_positions_is_whole():
    near = lambda rows, tr: annotate.whole_by_scoring_positions(rows, tr)
    start = np.array([[0, 120]])                                 # a terminal tract: alignments stop at its inner side
    assert span(near([row(0, 880, "+", 121, 1_000)], start)[0]) == (1, 1_000)
    assert span(near([row(0, 880, "+", 300, 1_000)], start)[0]) == (300, 1_000)                # a fragment stays one
    end = np.array([[880, 1_000]])
    assert span(near([row(0, 880, "+", 1, 880)], end)[0]) == (1, 1_000)
    assert span(near([row(0, 500, "+", 1, 500)], end)[0]) == (1, 500)
    inside = np.array([[60, 340]])                               # Chlamydomonas' 924-bp family: 60 scoring bp beyond the tract
    assert span(near([row(0, 660, "+", 340, 1_000)], inside)[0]) == (1, 1_000)                 # 660 of 720
    assert span(near([row(0, 950, "+", 1, 950)], TRACT_400)[0]) == (1, 1_000)                  # 850 of 900
    assert span(near([row(0, 900, "+", 1, 900)], TRACT_400)[0]) == (1, 900)                    # 800 of 900
    assert span(near([row(0, 700, "+", 300, 1_000)], TRACT_400)[0]) == (300, 1_000)            # 600 of 900


def test_only_the_hmm_positions_of_a_whole_hit_change():
    r = row(100, 980, "-", 121, 1_000, bits=321.0, evalue=1e-9, identity=0.77)
    assert annotate.whole_by_scoring_positions([r], np.array([[0, 120]]))[0] == row(100, 980, "-", 1, 1_000, bits=321.0, evalue=1e-9, identity=0.77)


# ---- pieces of one copy

def test_pieces_either_side_of_a_tract_are_one_copy():
    """Potato's 1,362-bp family at 99% identity: 9-514 and 608-1360. Where the model holds an internal duplicate, a jump to it
    outscores a walk over a tract that earns nothing, so the copy comes out in pieces."""
    a, b = row(100, 500, "+", 1, 400, bits=400.0, evalue=1e-40, identity=0.9), row(600, 1_100, "+", 501, 1_000, bits=600.0, evalue=1e-60, identity=0.8)
    joined, = annotate.join_split_copies([a, b], TRACT_400)
    assert joined[:4] == ("c", 100, 1_100, "+") and joined[4] == "f" and span(joined) == (1, 1_000) and joined[9] == 1_000
    assert joined[5] == 1_000.0 and joined[6] == 1e-60           # the scores add up; the best E-value stands
    assert joined[10] == pytest.approx((400 * 0.9 + 500 * 0.8) / 900)


def test_pieces_join_on_the_minus_strand_too():
    a, b = row(100, 600, "-", 501, 1_000), row(700, 1_100, "-", 1, 400)        # along the target the model runs down
    joined, = annotate.join_split_copies([a, b], TRACT_400)
    assert (joined[1], joined[2]) == (100, 1_100) and span(joined) == (1, 1_000)


@pytest.mark.parametrize("a, b, why", [
    (row(100, 500, "+", 1, 400), row(900, 1_400, "+", 501, 1_000), "a target gap beyond the tract's own length + 20 bp"),
    (row(100, 500, "+", 1, 300), row(800, 1_100, "+", 701, 1_000), "gap holds 300 positions that could score: over 10% of 900"),
    (row(100, 500, "+", 501, 900), row(510, 700, "+", 901, 1_000), "no tract between them"),
    (row(100, 500, "+", 501, 1_000), row(600, 1_000, "+", 1, 400), "the next copy of an array: not colinear"),
    (row(100, 500, "+", 1, 400), row(600, 1_100, "-", 501, 1_000), "other strand"),
    (row(100, 500, "-", 1, 400), row(600, 1_100, "-", 501, 1_000), "on the minus strand the model must run down"),
])
def test_other_neighbors_stay_apart(a, b, why):
    assert sorted(annotate.join_split_copies([a, b], TRACT_400)) == sorted([a, b]), why


def test_a_tract_may_differ_in_length_between_copies():
    """Simple sequence varies in length between copies; the rest of a model does not."""
    a = row(100, 500, "+", 1, 400)
    assert len(annotate.join_split_copies([a, row(500 + 60, 1_060, "+", 501, 1_000)], TRACT_400)) == 1      # 40 bp shorter
    assert len(annotate.join_split_copies([a, row(500 + 215, 1_215, "+", 501, 1_000)], TRACT_400)) == 1     # 115 bp longer
    assert len(annotate.join_split_copies([a, row(500 + 225, 1_225, "+", 501, 1_000)], TRACT_400)) == 2     # beyond 100 + 20


def test_a_copy_is_joined_past_a_hit_on_its_internal_duplicate():
    """The hit on the duplicate overlaps both pieces in the target, and stays a hit of its own."""
    a, dup, b = row(100, 500, "+", 1, 400), row(420, 600, "+", 700, 880, bits=100.0), row(600, 1_100, "+", 501, 1_000)
    out = annotate.join_split_copies([a, dup, b], TRACT_400)
    assert sorted(span(r) for r in out) == [(1, 1_000), (700, 880)]


def test_three_pieces_chain_through_two_tracts():
    tr = np.array([[300, 340], [600, 640]])
    pieces = [row(0, 300, "+", 1, 300), row(340, 600, "+", 341, 600), row(640, 1_000, "+", 641, 1_000)]
    joined, = annotate.join_split_copies(pieces, tr)
    assert span(joined) == (1, 1_000) and (joined[1], joined[2]) == (0, 1_000) and joined[5] == 1_500.0


def test_rows_without_a_split_come_back_as_they_were():
    rows = [row(100, 500, "+", 1, 400), row(5_000, 5_100, "+", 700, 800), row(9_000, 9_500, "-", 501, 1_000)]
    assert annotate.join_split_copies(rows, TRACT_400) == rows


# ---- copies come out whole and in one piece

def tandem_case(tmp_path, seed: int, copy, n: int = 10, spacer: int = 0, rate: float = 0.05):
    rng = np.random.default_rng(seed)
    decoys = "".join(dna(rng, 3_000) + "AT" * 80 for _ in range(10))
    body = "".join(mutate(rng, copy, rate) + dna(rng, spacer) for _ in range(n))
    return write(tmp_path, dna(rng, 20_000) + decoys + body + decoys + dna(rng, 20_000))


def test_copies_with_a_tract_inside_are_whole_as_in_a_plain_search(tmp_path):
    rng = np.random.default_rng(5)
    cons = dna(rng, 800) + TRACT + dna(rng, 800)
    hmm = single_hmm("fam", cons)
    g = tandem_case(tmp_path, 6, cons)
    plain, new = search(hmm, g, plain=True), search(hmm, g)
    assert len(whole(plain, hmm.M)) >= 9
    same_copies(whole(plain, hmm.M), whole(new, hmm.M))


def test_a_divergent_stretch_beside_the_tract_does_not_split_a_copy(tmp_path):
    """Masking the target split such a copy: an alignment against N stops short of the tract. A neutral position costs
    nothing to cross, so the copy stays one hit."""
    rng = np.random.default_rng(7)
    left, right = dna(rng, 800), dna(rng, 800)
    hmm = single_hmm("fam", left + TRACT + right)
    copies = [mutate(rng, left[:-40], 0.03) + dna(rng, 40) + TRACT + dna(rng, 40) + mutate(rng, right[40:], 0.03) for _ in range(10)]
    g = write(tmp_path, dna(rng, 20_000) + "".join(copies) + dna(rng, 20_000))
    plain, new = search(hmm, g, plain=True), search(hmm, g)
    assert len(whole(plain, hmm.M)) >= 8
    same_copies(whole(plain, hmm.M), whole(new, hmm.M))
    assert len(new) == len(plain)


@pytest.mark.parametrize("left, right", [(0, 1_500), (1_500, 0)])
def test_a_tract_at_a_model_end_leaves_copies_whole(tmp_path, left, right):
    """Sheep's 1,174-bp family: a neutral position scores nothing, so a local alignment stops before a terminal tract, and
    with a long tract the hit would no longer count as whole."""
    rng = np.random.default_rng(8)
    tract = "AT" * 120                                           # 240 bp, 14% of the model: more than whole hits may miss
    cons = dna(rng, left) + tract + dna(rng, right)
    hmm = single_hmm("fam", cons)
    g = tandem_case(tmp_path, 9, cons, spacer=300)
    plain, new = search(hmm, g, plain=True), search(hmm, g)
    assert len(whole(plain, hmm.M)) >= 8
    assert len(whole(new, hmm.M)) == len(whole(plain, hmm.M))


def _gc_rich(rng, n: int) -> str:
    return "".join("ACGT"[i] for i in rng.choice(4, size=n, p=[0.17, 0.35, 0.35, 0.13]))


@pytest.mark.parametrize("n", [1, 10])
def test_a_long_biased_tract_near_the_model_start_does_not_cost_a_copy(tmp_path, n):
    """Chlamydomonas' 924-bp family: 59 bp, a 279-bp G-rich tract, 586 bp. nhmmer's composition correction (null2) charges a
    hit for the biased residues it holds, even though neutral positions earned nothing on them, so hits stop at the tract and
    the copy comes out as 337-924 (in an array also 1-61). Turning that correction off is no answer: 45% more hits in potato,
    none of them seen by the plain model. So: over 90% of the positions that can score, and joined where both pieces show."""
    rng = np.random.default_rng(21)
    tract = ("GGCGCAGGCGACGGC" * 20)[:279]
    cons = _gc_rich(rng, 59) + tract + _gc_rich(rng, 586)
    hmm = single_hmm("fam", cons)
    g = write(tmp_path, dna(rng, 5_000) + "".join(mutate(rng, cons, 0.02) for _ in range(n)) + dna(rng, 5_000))
    plain, new = search(hmm, g, plain=True), search(hmm, g)
    assert len(whole(plain, hmm.M)) == n and len(whole(new, hmm.M)) == n
    a, b = whole(plain, hmm.M), whole(new, hmm.M)
    assert (np.minimum(a.end.to_numpy(), b.end.to_numpy()) - np.maximum(a.start.to_numpy(), b.start.to_numpy()) >= 0.6 * hmm.M).all()


def test_the_search_joins_the_pieces_of_a_copy(tmp_path, monkeypatch):
    """The whole path: a real split needs a model with an internal duplicate, so the search reports two pieces here."""
    rng = np.random.default_rng(23)
    hmm = single_hmm("fam", dna(rng, 400) + TRACT + dna(rng, 400))
    pieces = [("c", 1_000, 1_400, "+", "fam", 400.0, 1e-30, 1, 400, hmm.M, 0.9),
              ("c", 1_520, 1_920, "+", "fam", 400.0, 1e-30, 521, hmm.M, hmm.M, 0.9)]
    monkeypatch.setattr(annotate, "_rows", lambda *a: list(pieces))
    hits = search(hmm, write(tmp_path, dna(rng, 30_000)))
    assert hits[["start", "end", "hmm_from", "hmm_to"]].values.tolist() == [[1_000, 1_920, 1, hmm.M]]


def test_a_plain_model_keeps_its_pieces(tmp_path, monkeypatch):
    rng = np.random.default_rng(24)
    hmm = single_hmm("fam", dna(rng, 920))                       # no tract
    pieces = [("c", 1_000, 1_400, "+", "fam", 400.0, 1e-30, 1, 400, hmm.M, 0.9),
              ("c", 1_520, 1_920, "+", "fam", 400.0, 1e-30, 521, hmm.M, hmm.M, 0.9)]
    monkeypatch.setattr(annotate, "_rows", lambda *a: list(pieces))
    assert len(search(hmm, write(tmp_path, dna(rng, 30_000)))) == 2


def test_hits_of_the_tract_alone_are_dropped(tmp_path):
    """What goes on purpose: a stretch of the tract with nothing of the model around it marks no copy."""
    rng = np.random.default_rng(10)
    hmm = single_hmm("fam", dna(rng, 800) + TRACT + dna(rng, 800))
    g = write(tmp_path, "".join(dna(rng, 4_000) + "AT" * 200 for _ in range(30)))
    assert search(hmm, g).empty


def test_a_model_without_a_tract_is_searched_exactly_as_before(tmp_path):
    rng = np.random.default_rng(11)
    cons = dna(rng, 1_600)
    hmm = single_hmm("fam", cons)
    g = tandem_case(tmp_path, 12, cons)
    assert search(hmm, g, plain=True).equals(search(hmm, g))


def test_thread_count_does_not_change_the_hits(tmp_path):
    rng = np.random.default_rng(13)
    cons = dna(rng, 800) + TRACT + dna(rng, 800)
    hmm = single_hmm("fam", cons)
    g = tandem_case(tmp_path, 14, cons)
    assert search(hmm, g, threads=1).equals(search(hmm, g, threads=4))


def test_e_values_stay_near_those_of_the_plain_model(tmp_path):
    rng = np.random.default_rng(15)
    cons = dna(rng, 800) + TRACT + dna(rng, 800)
    hmm = single_hmm("fam", cons)
    g = tandem_case(tmp_path, 16, cons, rate=0.25)               # weak enough for E-values above the floor of 0
    plain, new = whole(search(hmm, g, plain=True), hmm.M), whole(search(hmm, g), hmm.M)
    same_copies(plain, new, tol=20)
    shown = (plain.evalue.to_numpy() > 0) & (new.evalue.to_numpy() > 0)
    assert shown.any() and (new.bitscore.to_numpy()[shown] <= plain.bitscore.to_numpy()[shown] + 1).all()
    assert (new.bitscore.to_numpy()[shown] > 47).all()           # each stays significant on its own (E 1e-5 on 500 Mb)


def test_cost_table_counts_one_search_per_job(tmp_path):
    rng = np.random.default_rng(17)
    cons = dna(rng, 800) + TRACT + dna(rng, 800)
    hmm = single_hmm("fam", cons)
    g = tandem_case(tmp_path, 18, cons)
    cost = {}
    search(hmm, g, cost=cost)
    assert cost["fam"]["searches"] == 1 and cost["fam"]["cpu_s"] > 0 and cost["fam"]["hits"] >= 10


def test_a_neutral_model_is_faster_on_low_complexity_rich_sequence(tmp_path):
    """The point of it, in nhmmer's own CPU time: a model with a tract, sequence full of such stretches, no copies."""
    rng = np.random.default_rng(19)
    hmm = single_hmm("fam", dna(rng, 800) + "AAT" * 60 + dna(rng, 800))
    g = write(tmp_path, "".join(dna(rng, 800) + "AAT" * 50 for _ in range(300)))
    plain, new = {}, {}
    search(hmm, g, plain=True, cost=plain)
    search(hmm, g, cost=new)
    assert new["fam"]["cpu_s"] < plain["fam"]["cpu_s"]
