"""Hi-C and fusion rules added after the T2T mammal benchmark: distances in mappable bins, partner windows on the
mappable flanks of an unmappable array, the gate's `untested` state, Hi-C choosing between two competing sequence
clusters, and low-confidence calls from clean medium peaks where the sequence line has none."""
import numpy as np

from cenfuse import fuse, hicfuse
from cenfuse.hic.bins import Bins

BIN = 250_000


def masked_chrom(lo=40, hi=60, nb=100):
    """One 25-Mb chromosome of 100 bins whose bins [lo, hi) are unmappable (an array of near-identical units)."""
    bins = Bins({"c": nb * BIN}, BIN)
    good = np.ones(bins.N, bool)
    good[lo:hi] = False
    return bins, good


def test_visible_distance_skips_unmappable_bins():
    bins, good = masked_chrom()
    d = hicfuse.visible_dist(bins, good, 0, 45 * BIN, 50 * BIN)          # call inside the array
    assert d[39] == 0 and d[60] == 0 and d[38] == 1 and d[61] == 1 and d[47] == 0
    dc = hicfuse.visible_dist(bins, good, 0, 45 * BIN, 50 * BIN, cap=4)  # 5 unmappable bins left, 10 right of it
    assert dc[39] == 1 and dc[60] == 6
    assert hicfuse._agrees_bin(bins, good, 0, 39, 45 * BIN, 50 * BIN)   # the array's near edge agrees
    assert not hicfuse._agrees_bin(bins, good, 0, 60, 45 * BIN, 50 * BIN)   # across 10 blind bins it does not
    assert hicfuse._agrees_bin(bins, good, 0, 60, 45 * BIN, 58 * BIN)   # a call spanning the array: its edge agrees


def test_partner_window_is_the_mappable_flank_of_an_unmappable_call():
    bins, good = masked_chrom()
    w = hicfuse.visible_window(bins, good, 0, 45 * BIN, 50 * BIN)
    assert len(w) == 9 and good[w].all() and set(w) <= set(range(35, 40)) | set(range(60, 65))
    allgood = np.ones(bins.N, bool)                                      # a mappable call: midpoint +/- 4
    w2 = hicfuse.visible_window(bins, allgood, 0, 45 * BIN, 50 * BIN)
    assert list(w2) == list(range(43, 52))


def test_visible_gap_and_chance_use_the_same_geometry():
    bins, good = masked_chrom()
    runs = hicfuse.masked_runs(bins, good, 0)
    assert runs == ((40 * BIN, 60 * BIN),)
    assert hicfuse.visible_gap(60 * BIN, 61 * BIN, 45 * BIN, 50 * BIN, runs) == 0
    assert hicfuse.visible_gap(70 * BIN, 71 * BIN, 45 * BIN, 50 * BIN, runs) == 10 * BIN
    assert hicfuse.visible_gap(60 * BIN, 61 * BIN, 45 * BIN, 50 * BIN, runs, 4 * BIN) == 6 * BIN
    S = np.where(good, 1.0, np.nan)
    # agreeing bins (blind bins free up to 4): 39 only on the left (5 blind bins), none on the right (10)
    assert hicfuse._chance(S, bins, good, 0, 45 * BIN, 50 * BIN) == 1 / 80
    # a call spanning the array (40-58): 38, 39 and 60, 61 agree
    assert hicfuse._chance(S, bins, good, 0, 40 * BIN, 58 * BIN) == 4 / 80


def ev(seq=(11_250_000, 12_500_000), conf="low", alt=None, peak=None, level="high", gate=True, strong=True,
       flags="", masked=(), gate_n=10):
    s = dict(start=seq[0] if seq else -1, end=seq[1] if seq else -1, evidence="satellite" if seq else "none",
             confidence=conf if seq else "none", flags="competing_cluster" if alt else ".", main_family="A",
             alt_start=alt[0] if alt else -1, alt_end=alt[1] if alt else -1)
    h = None
    if peak is not None:
        h = dict(bin_size=BIN, gate_pass=gate, gate_n=gate_n, anch_strong=strong, anch_level=level,
                 anch_flags=flags, anch_peak_start=peak, anch_peak_end=peak + BIN, anch_core_start=peak,
                 anch_core_end=peak + BIN, masked=masked)
    sat = (lambda a, b: 10_000) if seq else (lambda a, b: 0)
    return fuse.Evidence("c", 25_000_000, s, None, h, None, sat)


def test_peak_at_the_edge_of_the_unmappable_array_supports_the_call():
    runs = ((10_000_000, 15_000_000),)
    d = fuse.decide(ev(seq=(11_250_000, 14_500_000), peak=15_000_000, masked=runs))   # 2 blind bins between
    assert d["relation"] == "concordant" and "hic" in d["support"] and d["confidence"] == "medium"
    d = fuse.decide(ev(seq=(11_250_000, 14_500_000), peak=15_000_000))      # no mask: 2 mappable bins away
    assert d["relation"] != "concordant"
    d = fuse.decide(ev(peak=15_000_000, masked=runs))     # call 10 blind bins away: Hi-C cannot place it
    assert d["relation"] != "concordant" and "hic" not in d["support"]


def test_validated_hic_chooses_the_competing_cluster():
    d = fuse.decide(ev(alt=(20_000_000, 20_500_000), peak=20_250_000))
    assert d["call"] == (20_000_000, 20_500_000) and "hic_chose_alt" in d["flags"]
    assert d["relation"] == "concordant" and d["confidence"] == "medium"
    weak = fuse.decide(ev(alt=(20_000_000, 20_500_000), peak=20_250_000, level="medium", strong=False))
    assert weak["call"] == (11_250_000, 12_500_000) and "hic_chose_alt" not in weak["flags"]   # a move needs strong
    row = fuse.row(ev(alt=(20_000_000, 20_500_000), peak=20_250_000), d)
    assert row["offset_bp"] == "." and row["seq_alt_start"] == 20_000_001
    d = fuse.decide(ev(alt=(20_000_000, 20_500_000), peak=20_250_000, gate=False))     # not validated: no choice
    assert d["call"] == (11_250_000, 12_500_000) and "hic_chose_alt" not in d["flags"]


def test_clean_medium_peak_places_a_low_call_where_the_sequence_line_has_none():
    d = fuse.decide(ev(seq=None, peak=8_000_000, level="medium", strong=False))
    assert d["relation"] == "functional_only" and d["confidence"] == "low" and "hic_medium_peak" in d["flags"]
    for kw in (dict(flags="hub"), dict(flags="telomeric"), dict(gate=False), dict(level="low")):
        d = fuse.decide(ev(seq=None, peak=8_000_000, strong=False, **{"level": "medium", **kw}))
        assert d["call"] is None


def test_gate_status_untested_without_references():
    assert fuse.gate_status(dict(gate_pass=True, gate_n=5)) == "pass"
    assert fuse.gate_status(dict(gate_pass=False, gate_n=10)) == "fail"
    assert fuse.gate_status(dict(gate_pass=False, gate_n=0)) == "untested"
    assert fuse.gate_status(dict(gate_pass=False, gate_n=2)) == "untested"


def test_lopo_fast_path_equals_the_engine_rule():
    """lopo_all restricts the leave-one-partner-out score to the chromosome scored; same peaks as cluster_score."""
    from cenfuse.hic.score import cluster_score
    rng = np.random.default_rng(3)
    bins = Bins({f"c{i}": 40 * BIN for i in range(6)}, BIN)
    good = rng.random(bins.N) > 0.1
    members = np.ones(6, bool)
    oe = {y: np.where(rng.random(bins.N) > 0.2, rng.random(bins.N), np.nan) for y in range(6)}
    for y in oe:
        oe[y][bins.sl(y)] = np.nan
    fast = hicfuse.lopo_all(oe, bins, good, members)
    for x in range(6):
        ref = []
        for y in oe:
            if y == x:
                continue
            v = cluster_score(oe, bins, good, members, exclude={x: {y}})[bins.sl(x)]
            ref.append(int(np.nanargmax(v)) if np.isfinite(v).any() else -1)
        assert list(fast[x]) == ref


def test_gate_needs_a_majority_of_the_whole_reference_set():
    """2 of 4 references agree: leaving out a dissenter gives 2 of 3, but there is no majority, so no gate."""
    import pandas as pd
    rows = [dict(chrom=f"c{i}", seq_start=0, seq_end=1, seq_confidence="high", anch_peak_start=0, chance=0.01,
                 anch_agrees_seq=i < 2) for i in range(4)]
    assert not hicfuse.gate(pd.DataFrame(rows)).gate_pass.any()
    rows[2]["anch_agrees_seq"] = True                                   # 3 of 4: a majority
    assert hicfuse.gate(pd.DataFrame(rows)).gate_pass.all()            # every leave-one-out set: >= 2 of 3, P 3e-4
