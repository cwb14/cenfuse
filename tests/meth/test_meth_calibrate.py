import logging

import numpy as np
import pandas as pd
import pytest
from meth_common import run_scan

from cenfuse.meth import calibrate, params
from cenfuse.meth.scan import Pooled


def test_anchor_table_merges_and_drops(caplog):
    a = pd.DataFrame(dict(chrom=["c1", "c1", "c2", "scaf"], start=[1_000_000, 1_400_000, 5_000_000, 0],
                          end=[1_200_000, 1_600_000, 6_000_000, 10], mid=[1_100_000, 1_500_000, 5_500_000, 5],
                          name=["p1", "p2", "c2_cen", "s"]))
    with caplog.at_level(logging.WARNING, logger="cenfuse-meth"):
        t = calibrate.anchor_table(a, {"c1": 50_000_000, "c2": 50_000_000})
    assert set(t) == {"c1", "c2"}
    assert (t["c1"]["start"], t["c1"]["end"], t["c1"]["mid"]) == (1_000_000, 1_600_000, 1_300_000)
    assert t["c2"]["bin"] == 5_500_000 // params.BIN and t["c2"]["name"] == "c2_cen"
    assert "merged" in caplog.text and "not analysed" in caplog.text


def test_permutation_on_constructed_pools():
    rng = np.random.default_rng(3)
    pooled = {"CHG": {}}
    for i in range(6):
        anchor_num = np.array([-200.0, 0.0, 0.0, -200.0])         # class 0 drops 0.2 at every anchor; 'all' too
        anchor_vok = np.array([1000.0, 1000.0, 1000.0, 3000.0])
        noise = rng.normal(0, 10, 300)
        else_ = {0: (noise, np.full(300, 1000.0)), 1: (noise, np.full(300, 1000.0)),
                 3: (noise, np.full(300, 3000.0))}
        if i < 2:
            else_[2] = (noise, np.full(300, 1000.0))               # class 2 has data elsewhere on 2 chromosomes only
        pooled["CHG"][f"c{i}"] = Pooled(anchor_num, anchor_vok, else_)
    totals = {"CHG": np.array([100.0, 1000.0, 1000.0])}            # class 'a' is rare genome-wide: concentrated
    cal = calibrate.calibrate(pooled, totals, ["a", "b", "c"]).set_index("cls")
    assert cal.loc["a", "p"] == pytest.approx(1 / (1 + params.N_PERM)) and cal.loc["a", "in_signature"]
    assert cal.loc["a", "anchor_share"] == pytest.approx(1 / 3) and cal.loc["a", "genome_share"] == pytest.approx(1 / 21)
    assert cal.loc["b", "p"] > 0.05 and not cal.loc["b", "in_signature"]
    assert np.isnan(cal.loc["c", "p"]) and cal.loc["c", "chromosomes_tested"] == 2
    assert cal.loc["all", "q"] <= params.FDR and not cal.loc["all", "in_signature"]
    spread = {"CHG": np.array([1000.0, 1000.0, 1000.0])}            # 'a' as common as elsewhere: not a marker
    cal = calibrate.calibrate(pooled, spread, ["a", "b", "c"]).set_index("cls")
    assert cal.loc["a", "q"] <= params.FDR and not cal.loc["a", "enriched"] and not cal.loc["a", "in_signature"]


def _pools(rng, dev=(-200.0, 0.0), noise_sd=10.0, n_chrom=6):
    """Constructed pools: names a, b, c; a and b deviate by dev/1000 at every anchor; c has no data at the anchors."""
    pooled = {"CHG": {}}
    for i in range(n_chrom):
        num = np.array([dev[0], dev[1], 0.0, dev[0] + dev[1]])
        vok = np.array([1000.0, 1000.0, 0.0, 2000.0])
        noise = rng.normal(0, noise_sd, 300)
        else_ = {0: (noise, np.full(300, 1000.0)), 1: (noise, np.full(300, 1000.0)), 3: (noise, np.full(300, 2000.0))}
        pooled["CHG"][f"c{i}"] = Pooled(num, vok, else_)
    return pooled, {"CHG": np.array([10.0, 10.0, 10_000.0])}       # a, b rare genome-wide: concentrated


def test_perm_floor():
    m = int(params.FDR * (params.N_PERM + 1))
    assert calibrate.perm_floor_ok(m) and not calibrate.perm_floor_ok(m + 1)
    assert m >= 500                                              # ~500 class x context tests still allowed


def test_calibrate_warns_when_permutations_cannot_reach_the_fdr(monkeypatch, caplog):
    pooled, totals = _pools(np.random.default_rng(1))
    with caplog.at_level(logging.WARNING, logger="cenfuse-meth"):
        calibrate.calibrate(pooled, totals, ["a", "b", "c"])
    assert "too many class x context tests" not in caplog.text
    monkeypatch.setattr(params, "N_PERM", 20)                   # 3 tests: 3 / 21 > FDR
    with caplog.at_level(logging.WARNING, logger="cenfuse-meth"):
        calibrate.calibrate(pooled, totals, ["a", "b", "c"])
    assert "too many class x context tests" in caplog.text


def test_weak_entry_is_not_in_the_signature():
    pooled, totals = _pools(np.random.default_rng(2), dev=(-10.0, -200.0), noise_sd=0.5)   # a: -0.01, b: -0.2
    cal = calibrate.calibrate(pooled, totals, ["a", "b", "c"]).set_index("cls")
    a, b = cal.loc["a"], cal.loc["b"]
    assert a.mean_dev == pytest.approx(-0.01) and a.q <= params.FDR and a.enriched
    assert not a.detectable and not a.in_signature
    assert b.detectable and b.in_signature


def test_bh():
    q = calibrate.bh(np.array([0.01, 0.04, np.nan, 0.03]))
    assert np.allclose(q[[0, 1, 3]], [0.03, 0.04, 0.04]) and np.isnan(q[2])


def test_true_anchors_give_crm_chg_down_signature(prep, syn):
    _, _, cal = run_scan(prep, syn["anchors"], None)
    sig = calibrate.signature(cal)
    assert ("CHG", "LTR:CRM", "down") in set(zip(sig.context, sig.cls, sig.direction))
    assert "sat:satA" not in set(sig.cls) and "LTR:Tekay" not in set(sig.cls)
    assert calibrate.learned_direction(cal) == ("down", "anchors")


def test_wrong_anchors_give_no_signature(prep, syn):
    _, _, cal = run_scan(prep, syn["anchors_wrong"], None)
    assert calibrate.signature(cal).empty
    assert calibrate.learned_direction(cal) == ("down", "default")


def test_calibration_table_columns(prep, syn):
    _, _, cal = run_scan(prep, syn["anchors"], None)
    assert list(cal.columns) == ["context", "cls", "anchors_with_data", "mean_dev", "chromosomes_tested", "p",
                                 "anchor_share", "genome_share", "q", "direction", "enriched", "detectable",
                                 "in_signature"]
    assert set(cal.context) == {"CHG", "CG"}                 # S1 also has CG
    assert cal[cal.chromosomes_tested < params.MIN_ANCHORS].p.isna().all()
    assert (cal.p.dropna() >= 1 / (1 + params.N_PERM)).all()
