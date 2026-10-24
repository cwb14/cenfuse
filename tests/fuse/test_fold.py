"""Arm fold (fold.py): the anti-diagonal ridge of a chromosome's cis map, its validation against the de novo Hi-C
peaks, and the fusion rule that lets the two 3D lines move or lower a sequence call."""
import numpy as np
import pandas as pd

from cenfuse import fold, fuse

BIN = 250_000


def cis(nb=120, c=50.0, ridge=3.0):
    """Distance-decay contacts with a fold: pairs (c - d, c + d) at d >= S_MIN / 2 are `ridge` x enriched."""
    i, j = np.indices((nb, nb))
    W = 1.0 / (1.0 + np.abs(i - j))
    on = (np.abs(i + j - 2 * c) <= 1) & (np.abs(i - j) >= fold.S_MIN)
    return np.where(on, ridge * W, W)


def test_fold_point_on_the_ridge():
    W, good = cis(), np.ones(120, bool)
    F, n = fold.profile(W, good)
    f = fold.call(F, BIN, 120 * BIN)
    assert abs(f["point"] - 50.5 * BIN) <= BIN and f["core_start"] <= 50 * BIN < f["core_end"]
    assert f["core_end"] - f["core_start"] < 10 * BIN


def test_fold_point_inside_an_unmappable_array():
    """The array's bins are masked: the center is read from its mappable flanks."""
    good = np.ones(120, bool)
    good[45:56] = False
    F, _ = fold.profile(cis(), good)
    f = fold.call(F, BIN, 120 * BIN)
    assert 45 * BIN <= f["point"] <= 56 * BIN


def test_telocentric_centers_are_not_scanned():
    """A ridge at 7% of the chromosome (below CI_MIN) is not a fold point; the scan stays inside 12.5-87.5%."""
    F, _ = fold.profile(cis(c=8.0), np.ones(120, bool))
    f = fold.call(F, BIN, 120 * BIN)
    assert f is None or fold.CI_MIN * 120 * BIN <= f["point"] <= (1 - fold.CI_MIN) * 120 * BIN


def table(peaks):
    return pd.DataFrame([dict(chrom=f"c{i}", callable=True, denovo_confidence="medium", denovo_peak_start=p,
                              denovo_peak_end=p + BIN) for i, p in enumerate(peaks)])


def folds(n):
    return {f"c{i}": dict(point=50 * BIN, core_start=48 * BIN, core_end=53 * BIN, scanned=90 * BIN) for i in range(n)}


def test_validation_needs_a_majority_beyond_chance():
    v = fold.validate(folds(6), table([50 * BIN] * 5 + [10 * BIN]))
    assert v["status"] == "pass" and v["k"] == 5 and v["agrees"]["c5"] is False
    assert fold.validate(folds(6), table([50 * BIN] * 3 + [10 * BIN] * 3))["status"] == "fail"     # no majority
    assert fold.validate(folds(2), table([50 * BIN] * 2))["status"] == "untested"


def ev(seq=(500_000, 1_000_000), conf="low", fold_status="pass", agrees=True, dn=12_000_000, gate=False):
    s = dict(start=seq[0] if seq else -1, end=seq[1] if seq else -1, evidence="satellite" if seq else "none",
             confidence=conf if seq else "none", flags=".", main_family="A", alt_start=-1, alt_end=-1)
    h = dict(bin_size=BIN, gate_pass=gate, gate_n=5, anch_strong=False, anch_level="low", anch_flags="",
             anch_peak_start=-1, anch_peak_end=-1, anch_core_start=-1, anch_core_end=-1, masked=(),
             denovo_peak_start=dn, denovo_peak_end=dn + BIN, denovo_core_start=dn - BIN, denovo_core_end=dn + 2 * BIN)
    f = dict(point=12_500_000, core_start=11_000_000, core_end=14_000_000, scanned=18_000_000, status=fold_status,
             agrees=agrees)
    return fuse.Evidence("c", 25_000_000, s, None, h, None, lambda a, b: 0, f)


def test_two_3d_lines_move_a_call_at_neither():
    d = fuse.decide(ev())
    assert d["call"] == (11_750_000, 12_500_000) and d["relation"] == "offset" and d["confidence"] == "medium"
    assert "fold_pair" in d["flags"] and {"hic", "fold"} <= set(d["support"])
    d = fuse.decide(ev(conf="high"))                       # the fold is the second line at the peak, not at the call
    assert d["call"] == (11_750_000, 12_500_000) and "moved_with_fold" in d["flags"] and d["confidence"] == "medium"
    row = fuse.row(ev(), d)
    assert row["fold_gate"] == "pass" and row["fold_core_start"] == 11_000_001 and row["offset_bp"] != "."
    d = fuse.decide(ev(seq=None))
    assert d["relation"] == "functional_only" and d["confidence"] == "medium"


def test_disagreeing_3d_lines_lower_the_call_one_level_each():
    d = fuse.decide(ev(conf="high", agrees=False, dn=20_000_000))
    assert d["call"] == (500_000, 1_000_000) and d["relation"] == "contested" and d["confidence"] == "low"
    assert {"fold_offset", "hic_denovo_offset"} <= set(d["flags"]) and {"fold", "hic"} <= set(d["against"])
    d = fuse.decide(ev(conf="high", agrees=False, dn=700_000))    # de novo at the call, fold elsewhere: one level
    assert d["confidence"] == "medium" and d["flags"] == ["fold_offset"]


def test_a_call_inside_the_fold_core_is_kept_not_raised():
    d = fuse.decide(ev(seq=(12_000_000, 12_250_000), conf="medium"))
    assert d["call"] == (12_000_000, 12_250_000) and d["confidence"] == "medium" and d["relation"] == "sequence_only"
    assert "fold" in d["support"] and not d["against"]


def test_fold_acts_only_without_the_anchored_gate_and_after_its_own_test():
    for kw in (dict(gate=True), dict(fold_status="fail"), dict(fold_status="untested")):
        d = fuse.decide(ev(conf="high", **kw))
        assert d["call"] == (500_000, 1_000_000) and d["confidence"] == "high" and "fold_pair" not in d["flags"]
