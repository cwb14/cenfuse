"""The AT-trough line (cenfuse.atline) and its rules in fuse.decide. Synthetic FASTA only."""
import numpy as np
import pysam
import pytest

from cenfuse import atline, fuse

RNG = np.random.default_rng(11)
NT = np.array(list("ACGT"))
BIN = 10_000


def seq(n, gc):
    p = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return RNG.choice(NT, size=n, p=p)


def fasta(tmp_path, chroms, troughs=None, n_run=None):
    """One FASTA at 40% GC; troughs: {chrom: (start, end, gc)} planted; n_run: {chrom: (start, end)} set to N."""
    f = tmp_path / "g.fa"
    with open(f, "w") as o:
        for c, L in chroms.items():
            a = seq(L, 0.40)
            if troughs and c in troughs:
                s, e, g = troughs[c]
                a[s:e] = seq(e - s, g)
            if n_run and c in n_run:
                a[n_run[c][0]:n_run[c][1]] = "N"
            o.write(f">{c}\n{''.join(a)}\n")
    pysam.faidx(str(f))
    return f


CHROMS = {f"c{i}": 300_000 + 20_000 * i for i in range(1, 11)}


@pytest.fixture(scope="module")
def planted(tmp_path_factory):
    tr = {c: (100_000 + 1000 * i, 104_000 + 1000 * i, 0.10) for i, c in enumerate(CHROMS)}
    f = fasta(tmp_path_factory.mktemp("at"), CHROMS, tr)
    return f, tr, atline.scan(atline.gc_bins(f, CHROMS), CHROMS)


def test_gc_bins_counts_gc_and_called_bases(tmp_path):
    f = tmp_path / "t.fa"
    f.write_text(">c\n" + "G" * 1000 + "A" * 500 + "N" * 500 + "\n")
    pysam.faidx(str(f))
    g, m = atline.gc_bins(f, {"c": 2000})["c"]
    assert g.tolist()[:2] == [1000, 0] and m.tolist()[:2] == [1000, 500]


def test_windows_skip_mostly_unknown_sequence():
    g, m = np.array([100, 100, 0, 100]), np.array([1000, 1000, 100, 1000])
    f = atline._windows(g, m, 1)
    assert f[0] == pytest.approx(0.1) and np.isinf(f[2])


def test_troughs_are_independent():
    f = np.full(1000, 0.4)
    f[100:104] = 0.1                                         # one regional centromere, several low windows
    f[500] = 0.2
    t = atline.troughs(f, 2, 3)
    assert [j for _, j in t][:2] == [100, 500] and abs(t[2][1] - 100) > atline.APART * 2


def test_one_trough_per_chromosome_passes(planted):
    f, tr, r = planted
    assert r["status"] == "pass" and r["best"] <= 5_000
    x = r["scales"][r["best"]]
    assert x["covered"] == len(CHROMS)
    for c, (s, e, _) in tr.items():
        assert x["calls"][c]["start"] <= (s + e) // 2 <= x["calls"][c]["end"] and x["calls"][c]["first"]


def test_troughs_piled_on_one_chromosome_fail(tmp_path):
    chroms = dict(list(CHROMS.items()))
    many = [f for f in chroms]
    f = tmp_path / "g.fa"
    with open(f, "w") as o:
        for c, L in chroms.items():
            a = seq(L, 0.40)
            if c == many[0]:
                for s in range(30_000, L - 10_000, 30_000):
                    a[s:s + 4000] = seq(4000, 0.10)
            o.write(f">{c}\n{''.join(a)}\n")
    pysam.faidx(str(f))
    r = atline.scan(atline.gc_bins(f, chroms), chroms)
    assert r["status"] == "fail"


def test_no_trough_at_all_fails(tmp_path):
    f = fasta(tmp_path, CHROMS)
    assert atline.scan(atline.gc_bins(f, CHROMS), CHROMS)["status"] == "fail"


def test_too_few_chromosomes_are_untested(tmp_path):
    ch = dict(list(CHROMS.items())[:3])
    tr = {c: (100_000, 104_000, 0.10) for c in ch}
    f = fasta(tmp_path, ch, tr)
    assert atline.scan(atline.gc_bins(f, ch), ch)["status"] == "untested"


def _synth_bins(lengths, troughs=(), gc=0.45, seed=3):
    """GC bins of 1 kb with noise: {chrom: (GC counts, called)}; troughs: (chrom, start bin, bins, gc)."""
    rng = np.random.default_rng(seed)
    out = {}
    for c, L in lengths.items():
        nb = L // atline.BIN + 1
        g = rng.binomial(atline.BIN, gc, nb).astype(np.int32)
        for t_c, s, n, t_gc in troughs:
            if t_c == c:
                g[s:s + n] = rng.binomial(atline.BIN, t_gc, n)
        out[c] = (g, np.full(nb, atline.BIN, dtype=np.int32))
    return out


def test_deepest_class_is_the_largest_step_in_depth():
    d = np.array([0.1] * 8 + [0.30, 0.31, 0.32, 0.33])
    assert atline.deepest_class(d, 6) == 8                    # the step lies after the eighth deepest, between N and 1.5 N
    d = np.array([0.10, 0.11, 0.12, 0.13, 0.14, 0.30, 0.31, 0.32, 0.33, 0.34])
    assert atline.deepest_class(d, 5) == 5                    # already at N
    assert atline.deepest_class(np.array([0.1, 0.2]), 5) == 2 # no more troughs than that


def test_a_chromosome_never_holds_more_deepest_troughs_than_it_has_room_for():
    cap = np.array([1, 1, 4])
    n = atline._counts(atline._deal(np.array([10, 10, 40]), cap, 6), 3)
    assert (n <= cap).all() and (n.sum(axis=1) == 6).all()


def test_more_than_the_chromosomes_of_a_karyotype_are_not_tested():
    lengths = {f"s{i}": 200_000 for i in range(atline.MAX_CHROMS + 1)}
    r = atline.scan(_synth_bins(lengths), lengths)
    assert r["status"] == "untested" and r["scales"] == {}


def test_contigs_are_not_tested_at_the_windows_only_some_of_them_reach():
    """The count of sequences decides: 300 contigs, of which 180 reach ten windows at the smaller sizes, are contigs."""
    lengths = {f"s{i}": 200_000 if i < 180 else 3_000 for i in range(300)}
    r = atline.scan(_synth_bins(lengths), lengths)
    assert r["status"] == "untested" and r["scales"] == {}


def test_evenness_ranks_even_counts_below_piled_ones_and_shares_by_length():
    share = np.ones(5) / 5
    even, piled = np.array([2, 2, 2, 2, 2]), np.array([10, 0, 0, 0, 0])
    assert atline.evenness(even, 10, share) < atline.evenness(piled, 10, share)
    # one centromere on each of 8 chromosomes, the longest with its extra troughs as its length allows
    share = np.array([9.3, 4.3, 4.2, 4.2, 4.1, 3.5, 3.4, 3.3])
    share = share / share.sum()
    fair, clumped = np.array([3, 1, 1, 1, 1, 1, 1, 1]), np.array([4, 4, 1, 1, 0, 0, 0, 0])
    assert atline.evenness(fair, 10, share) < atline.evenness(clumped, 10, share)


def test_a_chromosome_with_extra_deep_troughs_still_passes_if_its_extra_share_is_fair():
    """V. dahliae-like: 8 chromosomes, the longest holds three deep AT troughs, the others one each. The N deepest cover 6 of 8
    (the former test failed), but the deepest class falls more evenly than in proportion to length."""
    lengths = {"a": 9_300_000, "b": 4_300_000, "c": 4_200_000, "d": 4_150_000, "e": 4_100_000, "f": 3_500_000, "g": 3_400_000, "h": 3_300_000}
    tr = [("a", 2900, 100, 0.20), ("a", 7000, 100, 0.20), ("a", 500, 100, 0.20)] + [(c, 1800, 100, 0.20) for c in lengths if c != "a"]
    r = atline.scan(_synth_bins(lengths, tr), lengths)
    assert r["status"] == "pass" and r["scales"][r["best"]]["m"] >= 8


def test_troughs_that_fall_at_random_rarely_pass():
    lengths = {f"c{i}": 300_000 + 20_000 * i for i in range(1, 11)}
    passes = sum(atline.scan(_synth_bins(lengths, seed=200 + i), lengths)["status"] == "pass" for i in range(20))
    assert passes <= 4              # about the level, in pure noise: the calibration is measured on the benchmark's GC shuffles


def _peaks(calls, shift=0):
    return {c: (x["start"] + shift, x["start"] + shift + BIN) for c, x in calls.items()}


def test_validate_against_hic_peaks(planted):
    _, _, r = planted
    calls = r["scales"][r["best"]]["calls"]
    v = atline.validate(r, _peaks(calls), dict.fromkeys(calls, BIN), CHROMS)
    assert v["status"] == "pass" and v["k"] == len(CHROMS)
    far = atline.validate(r, _peaks(calls, 150_000), dict.fromkeys(calls, BIN), CHROMS)
    assert far["status"] == "fail"


def test_validate_beyond_chance_on_a_minority_supports(planted):
    _, _, r = planted
    calls = r["scales"][r["best"]]["calls"]
    pk = {c: (x["start"], x["start"] + 1000) for c, x in calls.items()}
    for c in list(pk)[4:]:
        pk[c] = (pk[c][0] + 150_000, pk[c][1] + 150_000)
    v = atline.validate(r, pk, dict.fromkeys(calls, 1000), CHROMS)
    assert v["status"] == "support" and v["k"] == 4


def test_validate_needs_chromosomes(planted):
    _, _, r = planted
    assert atline.validate(r, {}, {}, CHROMS)["status"] == "untested"


def test_on_chromosome(planted):
    _, tr, r = planted
    res = dict(r, window=r["best"], hic=None)
    x = atline.on_chromosome(res, "c1")
    assert x["status"] == "pass" and not x["validated"] and not x["agrees"]
    k, s0 = r["best"], res["scales"][r["best"]]["calls"]["c1"]                   # the deepest window, inside the half-depth run
    assert s0["start"] <= x["trough_start"] and x["trough_end"] <= s0["end"] and x["trough_end"] - x["trough_start"] == k
    assert atline.on_chromosome(None, "c1") is None and atline.on_chromosome(res, "zz") is None


# ---------------------------------------------------------------- the rules
HB = 10_000


def ev(at=None, denovo=None, seq=None):
    s = dict(start=seq[0] if seq else -1, end=seq[1] if seq else -1, evidence="satellite" if seq else "none",
             confidence="high" if seq else "none", flags=".", main_family="A" if seq else ".")
    h = None
    if denovo is not None:
        h = dict(callable=True, bin_size=HB, gate_pass=False, gate_support=False, gate_genome_untested=True, anch_strong=False,
                 anch_level="none", anch_peak_start=-1, anch_peak_end=-1, anch_core_start=-1, anch_core_end=-1, masked=(),
                 denovo_peak_start=denovo, denovo_peak_end=denovo + HB, denovo_core_start=denovo, denovo_core_end=denovo + 2 * HB,
                 denovo_confidence="high", denovo_flags="")
    return fuse.Evidence("c", 1_000_000, s, None, h, None, lambda a, b: 0, at=at)


TROUGH = dict(start=300_000, end=303_000, gc=0.1, status="pass", window=1000, validated=False, support=False, agrees=False)


def test_a_passing_trough_places_a_low_call():
    d = fuse.decide(ev(at=TROUGH))
    assert d["call"] == (300_000, 303_000) and d["confidence"] == "low" and d["source"] == "composition"
    assert "at_only" in d["flags"] and d["support"] == ["at"]


def test_a_failing_trough_places_nothing():
    d = fuse.decide(ev(at=dict(TROUGH, status="fail")))
    assert d["call"] is None and d["relation"] == "none"


def test_a_trough_validated_by_hic_places_a_medium_call():
    d = fuse.decide(ev(at=dict(TROUGH, status="fail", validated=True)))
    assert d["call"] == (300_000, 303_000) and d["confidence"] == "medium"


def test_a_trough_never_replaces_a_sequence_call():
    d = fuse.decide(ev(at=dict(TROUGH, validated=True, agrees=True), seq=(700_000, 900_000)))
    assert d["call"] == (700_000, 900_000) and "at" not in d["support"] and "at_only" not in d["flags"]


def _low_alt(alt=(100_000, 900_000)):
    """Evidence of a low sequence call at 2-5 Mb with a competing cluster alt."""
    return dict(start=2_000_000, end=5_000_000, evidence="satellite", confidence="low", flags="competing_cluster", main_family="A",
                alt_start=alt[0], alt_end=alt[1])


def test_a_passing_trough_chooses_the_competing_cluster_of_a_low_call():
    e = fuse.Evidence("c", 10_000_000, _low_alt(), None, None, None, lambda a, b: 0, at=dict(TROUGH, start=300_000, end=500_000))
    d = fuse.decide(e)
    assert d["call"] == (100_000, 900_000) and d["confidence"] == "low" and d["rule"] == "at_chose_alt"
    assert "at_chose_alt" in d["flags"] and "at" in d["support"]


def test_the_deepest_window_decides_where_the_half_depth_run_touches_both_clusters():
    at = dict(TROUGH, start=0, end=4_000_000, trough_start=300_000, trough_end=500_000)  # a broad run over both, its deepest window on the alt
    d = fuse.decide(fuse.Evidence("c", 10_000_000, _low_alt(), None, None, None, lambda a, b: 0, at=at))
    assert d["call"] == (100_000, 900_000) and d["rule"] == "at_chose_alt"
    at = dict(TROUGH, start=0, end=4_000_000, trough_start=3_000_000, trough_end=3_200_000)                # its deepest window on the call
    d = fuse.decide(fuse.Evidence("c", 10_000_000, _low_alt(), None, None, None, lambda a, b: 0, at=at))
    assert d["call"] == (2_000_000, 5_000_000) and "at_chose_alt" not in d["flags"]


def test_a_failing_or_elsewhere_trough_leaves_a_low_call_as_it_is():
    fail = fuse.Evidence("c", 10_000_000, _low_alt(), None, None, None, lambda a, b: 0,
                         at=dict(TROUGH, status="fail", start=300_000, end=500_000))
    assert fuse.decide(fail)["call"] == (2_000_000, 5_000_000)
    elsewhere = fuse.Evidence("c", 10_000_000, _low_alt(), None, None, None, lambda a, b: 0,
                              at=dict(TROUGH, start=8_000_000, end=8_100_000))
    assert fuse.decide(elsewhere)["call"] == (2_000_000, 5_000_000)
    on_call = fuse.Evidence("c", 10_000_000, _low_alt(), None, None, None, lambda a, b: 0,
                            at=dict(TROUGH, start=3_000_000, end=3_100_000))
    assert fuse.decide(on_call)["call"] == (2_000_000, 5_000_000)


def _validated_hic(peak):
    """A strong validated Hi-C peak at `peak`, as the anchored line reports one."""
    return dict(callable=True, bin_size=HB, gate_pass=True, gate_support=True, gate_genome_untested=False, anch_strong=True, anch_level="high",
                anch_peak_start=peak, anch_peak_end=peak + HB, anch_core_start=peak, anch_core_end=peak + HB, masked=(), denovo_peak_start=-1,
                denovo_peak_end=-1, denovo_core_start=-1, denovo_core_end=-1, denovo_confidence="none", denovo_flags="")


def test_a_validated_peak_beside_a_low_call_keeps_it_from_the_trough():
    h = _validated_hic(5_000_000 + 2 * HB)                                       # beyond the agreement tolerance, not far
    e = fuse.Evidence("c", 10_000_000, _low_alt(), None, h, None, lambda a, b: 0, at=dict(TROUGH, start=300_000, end=500_000))
    d = fuse.decide(e)
    assert "hic_near" in d["flags"] and "at_chose_alt" not in d["flags"] and d["call"] == (2_000_000, 5_000_000)


def test_the_trough_still_chooses_where_hic_cannot_tell_the_two_clusters_apart():
    """A peak on the call, with no mappable bin between the call and its competing cluster: Hi-C has no say (mouse: T2T arrays)."""
    h = dict(_validated_hic(2_000_000), masked=[(900_000, 2_000_000)])
    at = dict(TROUGH, start=300_000, end=500_000)
    d = fuse.decide(fuse.Evidence("c", 10_000_000, _low_alt(), None, h, None, lambda a, b: 0, at=at))
    assert "hic_alt_ambiguous" in d["flags"] and d["rule"] == "at_chose_alt" and d["call"] == (100_000, 900_000)


def test_methylation_is_read_at_the_cluster_the_trough_chose():
    """The verdict on the abandoned call says nothing about the chosen cluster: the provisional item (the call as it stands) does."""
    at = dict(TROUGH, start=300_000, end=500_000)
    like, empty = dict(label="centromere_like", score=0.9), dict(label="", score=float("nan"))
    d = fuse.decide(fuse.Evidence("c", 10_000_000, _low_alt(), None, None, dict(sequence_call=like, provisional=empty), lambda a, b: 0, at=at))
    assert d["rule"] == "at_chose_alt" and "meth" not in d["support"] and d["confidence"] == "low" and d["meth_call"]["label"] == ""
    d = fuse.decide(fuse.Evidence("c", 10_000_000, _low_alt(), None, None, dict(sequence_call=empty, provisional=like), lambda a, b: 0, at=at))
    assert d["rule"] == "at_chose_alt" and "meth" in d["support"] and d["confidence"] == "medium"


def test_methylation_is_read_at_the_cluster_hic_chose():
    h = _validated_hic(400_000)                                                   # on the competing cluster
    like, empty = dict(label="centromere_like", score=0.9), dict(label="", score=float("nan"))
    d = fuse.decide(fuse.Evidence("c", 10_000_000, _low_alt(), None, h, dict(sequence_call=like, provisional=empty), lambda a, b: 0))
    assert "hic_chose_alt" in d["flags"] and d["call"] == (100_000, 900_000) and "meth" not in d["support"]


def test_a_trough_never_moves_a_confident_call_to_its_competing_cluster():
    high = dict(_low_alt(), confidence="medium")
    e = fuse.Evidence("c", 10_000_000, high, None, None, None, lambda a, b: 0, at=dict(TROUGH, start=300_000, end=500_000))
    assert fuse.decide(e)["call"] == (2_000_000, 5_000_000)


def test_hic_alone_with_an_agreeing_trough_rises_a_level_and_takes_its_interval():
    d = fuse.decide(ev(at=dict(TROUGH, validated=True, agrees=True), denovo=300_000))
    assert d["call"] == (300_000, 303_000) and d["confidence"] == "medium" and "at_pair" in d["flags"] and "at" in d["support"]
    assert "hic" in d["support"]
    d = fuse.decide(ev(at=dict(TROUGH, support=True, agrees=True), denovo=300_000))     # the support tier does the same
    assert d["confidence"] == "medium" and "at_pair" in d["flags"]


def test_hic_alone_against_a_validated_trough_is_flagged_not_moved():
    d = fuse.decide(ev(at=dict(TROUGH, validated=True), denovo=600_000))
    assert d["call"] == (600_000, 620_000) and d["confidence"] == "low" and "at_elsewhere" in d["flags"] and "at" in d["against"]


def test_hic_alone_with_an_untested_trough_changes_nothing():
    d = fuse.decide(ev(at=TROUGH, denovo=600_000))
    assert d["call"] == (600_000, 620_000) and d["confidence"] == "low" and not {"at_pair", "at_elsewhere"} & set(d["flags"])


# ---------------------------------------------------------------- what the code review found
def _ltr_lines(te_pair=None, cis=None, at=None):
    h = dict(callable=True, bin_size=HB, gate_pass=False, gate_support=False, gate_genome_untested=False, anch_strong=False, anch_level="none",
             anch_peak_start=-1, anch_peak_end=-1, anch_core_start=-1, anch_core_end=-1, masked=())
    if te_pair:
        h["te_pair"] = te_pair
    if cis:
        h["cis"] = cis
    s = dict(start=-1, end=-1, evidence="none", confidence="none", flags=".", main_family=".")
    return fuse.decide(fuse.Evidence("c", 1_000_000, s, None, h, None, lambda a, b: 0, at=at))


def test_a_trough_places_only_what_nothing_else_places():
    d = _ltr_lines(te_pair=dict(core=(600_000, 620_000), status="pass", clade="CRM"), at=TROUGH)
    assert d["call"] == (600_000, 620_000) and d["confidence"] == "medium" and "te_hic_pair" in d["flags"] and "at_only" not in d["flags"]
    d = _ltr_lines(cis=dict(peak=(600_000, 610_000), status="pass"), at=TROUGH)
    assert d["call"] == (600_000, 610_000) and "cis_only" in d["flags"] and "at_only" not in d["flags"]


def test_a_trough_far_from_the_hic_call_does_not_take_it():
    """`agrees` says the trough lies by the prior-free peak, not by this call."""
    far = dict(TROUGH, support=True, agrees=True)                                  # the trough at 300 kb, the cis peak at 600 kb
    d = _ltr_lines(cis=dict(peak=(600_000, 610_000), status="pass"), at=far)
    assert d["call"] == (600_000, 610_000) and "at_pair" not in d["flags"] and d["confidence"] == "low"
    d = _ltr_lines(cis=dict(peak=(300_000, 310_000), status="pass"), at=far)       # by the call: the pair stands
    assert d["call"] == (300_000, 303_000) and "at_pair" in d["flags"]


def test_a_validated_trough_near_the_call_without_a_peak_there_is_not_elsewhere():
    near = dict(TROUGH, validated=True, agrees=False)                              # no prior-free peak on this chromosome
    d = _ltr_lines(cis=dict(peak=(300_000, 310_000), status="pass"), at=near)
    assert "at_elsewhere" not in d["flags"] and "at_pair" not in d["flags"]
    d = _ltr_lines(cis=dict(peak=(600_000, 610_000), status="pass"), at=near)
    assert "at_elsewhere" in d["flags"]


def _two_windows():
    """A FASTA test that passed at one window, and Hi-C agreeing at another."""
    calls = lambda s: {f"c{i}": dict(start=s + i * 1000, end=s + i * 1000 + 500, gc=0.1, first=True) for i in range(6)}
    return dict(scales={1000: dict(calls=calls(100_000), p=0.01), 50_000: dict(calls=calls(400_000), p=0.9)}, best=1000, p=0.02, status="pass")


def test_hic_support_does_not_switch_a_passing_fasta_test_to_another_window():
    res = _two_windows()
    lengths = {f"c{i}": 1_000_000 for i in range(6)}
    peaks = {f"c{i}": (400_000 + i * 1000, 400_000 + i * 1000 + 500) for i in range(3)}      # 3 of 6 at the untested window
    peaks.update({f"c{i}": (800_000, 800_500) for i in range(3, 6)})
    v = atline.validate(res, peaks, {c: 100 for c in peaks}, lengths)
    assert v["status"] == "support" and v["window"] == 50_000 and set(v["by_window"]) == {1000, 50_000}
    out = dict(res, window=res["best"], hic=v)
    r = atline.on_chromosome(out, "c0")
    assert r["start"] == 100_000 and r["window"] == 1000 and r["support"] and not r["agrees"]     # read at the window that passed


def test_with_hic_follows_a_window_only_where_that_says_something():
    res = _two_windows()
    lengths = {f"c{i}": 1_000_000 for i in range(6)}
    table = lambda ps: __import__("pandas").DataFrame([dict(chrom=c, callable=True, denovo_peak_start=a, denovo_peak_end=b, bin_size=100) for c, (a, b) in ps.items()])
    onw = {f"c{i}": (400_000 + i * 1000, 400_000 + i * 1000 + 500) for i in range(6)}          # all at the second window
    assert atline.with_hic(res, dict(table=table(onw)), lengths)["window"] == 50_000              # validated there
    some = {f"c{i}": (400_000 + i * 1000, 400_000 + i * 1000 + 500) for i in range(3)}
    some.update({f"c{i}": (800_000, 800_500) for i in range(3, 6)})
    assert atline.with_hic(res, dict(table=table(some)), lengths)["window"] == 1000               # support: stays
    fail = dict(res, status="fail")
    assert atline.with_hic(fail, dict(table=table(some)), lengths)["window"] == 50_000            # the FASTA test said nothing


# ---------------------------------------------------------------- 0.16: the trough as second line for the cis readout (#33)
def test_troughs_are_tested_against_the_cis_extremes_where_the_prior_free_peaks_cannot():
    """C. neoformans in G1: the prior-free peaks sit at the clustered telomeres, the cis extremes at the centromeres."""
    res = dict(_two_windows(), status="fail")
    lengths = {f"c{i}": 1_000_000 for i in range(6)}
    ends = {f"c{i}": (0, 100) for i in range(6)}                                                   # prior-free peaks at the ends
    tab = __import__("pandas").DataFrame([dict(chrom=c, callable=True, denovo_peak_start=a, denovo_peak_end=b, bin_size=100)
                                          for c, (a, b) in ends.items()])
    cis = dict(status="pass", tested_by="trans", peaks={f"c{i}": (400_000 + i * 1000, 400_000 + i * 1000 + 100) for i in range(6)})
    out = atline.with_hic(res, dict(table=tab, cis=cis), lengths)
    assert out["hic"]["status"] == "fail" and out["cis"]["status"] == "pass" and out["window"] == 50_000
    r = atline.on_chromosome(out, "c2")
    assert r["cis_validated"] and r["cis_agrees"] and not r["validated"] and not r["agrees"]
    out = atline.with_hic(res, dict(table=tab, cis=dict(cis, status="support")), lengths)     # a trans-tested line acts at pass only
    assert out["cis"] is None


def test_the_cis_test_is_not_run_where_the_prior_free_peaks_already_test_the_troughs():
    res = _two_windows()
    lengths = {f"c{i}": 1_000_000 for i in range(6)}
    onw = {f"c{i}": (400_000 + i * 1000, 400_000 + i * 1000 + 500) for i in range(6)}
    tab = __import__("pandas").DataFrame([dict(chrom=c, callable=True, denovo_peak_start=a, denovo_peak_end=b, bin_size=100)
                                          for c, (a, b) in onw.items()])
    cis = dict(status="pass", tested_by="trans", peaks={f"c{i}": (800_000, 800_100) for i in range(6)})
    out = atline.with_hic(res, dict(table=tab, cis=cis), lengths)
    assert out["hic"]["status"] == "pass" and out["cis"] is None


def test_a_cis_call_with_a_trough_tested_against_the_cis_extremes_is_medium():
    tested = dict(TROUGH, cis_validated=True, cis_agrees=True)
    d = _ltr_lines(cis=dict(peak=(300_000, 310_000), status="pass"), at=tested)
    assert d["call"] == (300_000, 303_000) and "at_pair" in d["flags"] and "cis_only" in d["flags"]
    assert d["confidence"] == "medium" and {"hic", "at"} <= set(d["support"])
    d = _ltr_lines(cis=dict(peak=(300_000, 310_000), status="pass"), at=dict(TROUGH, cis_support=True, cis_agrees=True))
    assert d["confidence"] == "medium" and "at_pair" in d["flags"]                       # support: where the two agree


def test_a_cis_call_keeps_low_without_a_tested_or_agreeing_trough():
    d = _ltr_lines(cis=dict(peak=(300_000, 310_000), status="pass"), at=dict(TROUGH, cis_validated=True, cis_agrees=False))
    assert d["call"] == (300_000, 310_000) and "at_pair" not in d["flags"] and d["confidence"] == "low"
    d = _ltr_lines(cis=dict(peak=(300_000, 310_000), status="pass"), at=dict(TROUGH, cis_agrees=True))       # untested trough
    assert "at_pair" not in d["flags"] and d["confidence"] == "low"
    d = _ltr_lines(cis=dict(peak=(600_000, 610_000), status="pass"), at=dict(TROUGH, cis_validated=True, cis_agrees=True))
    assert "at_pair" not in d["flags"] and d["confidence"] == "low"                      # the trough is not by this call


def test_the_cis_test_of_the_troughs_does_not_pair_a_prior_free_call():
    """A trough tested against the cis extremes says nothing about a call the prior-free peak placed."""
    d = fuse.decide(ev(at=dict(TROUGH, cis_validated=True, cis_agrees=True), denovo=300_000))
    assert "at_pair" not in d["flags"] and d["confidence"] == "low"
