import numpy as np
import pandas as pd
import pytest

from cenfuse.meth import candidates, params, track
from cenfuse.meth.scan import median_of_others, pooled_at, scan


def test_pick_separation_cap_and_direction():
    geo = track.geometry(0.5e6)                      # sep = 40 bins
    y = np.zeros(1000)
    y[[100, 120, 300, 500, 700, 800, 900]] = [-0.5, -0.6, -0.4, -0.3, -0.2, -0.1, -0.05]
    assert candidates.pick(y, geo, "down") == [120, 300, 500, 700, 800]
    assert candidates.pick(-y, geo, "up") == [120, 300, 500, 700, 800]
    assert candidates.pick(np.full(100, 0.1), geo, "down") == []
    assert candidates.pick(np.full(100, np.nan), geo, "down") == []


def test_half_width_and_robust_z():
    y = np.array([0, -0.02, -0.06, -0.1, -0.07, -0.04, np.nan, -0.1, 0])
    assert candidates.half_width(y, 3, "down") == (2, 4)
    z = np.random.default_rng(1).normal(0, 0.01, 500)
    z[250] = -0.1
    assert candidates.robust_z(z, 250) < -8


def test_signature_score_one_class():
    names = ["LTR:CRM", "other"]
    sig = pd.DataFrame(dict(context=["CHG"], cls=["LTR:CRM"], mean_dev=[-0.2], direction=["down"]))
    vok = {"CHG": np.array([1000.0, 1000.0])}
    wd, mt, sc = candidates.signature_match({"CHG": np.array([-150.0, 0])}, vok, sig, names)
    assert (wd, mt) == (1, 1) and sc == pytest.approx(0.75)
    assert candidates.signature_match({"CHG": np.array([-50.0, 0])}, vok, sig, names)[1:] == (0, pytest.approx(0.25))
    assert candidates.signature_match({"CHG": np.array([80.0, 0])}, vok, sig, names)[2] == pytest.approx(-0.4)
    low = {"CHG": np.array([params.MIN_CALLS - 1.0, 1000.0])}
    wd, mt, sc = candidates.signature_match({"CHG": np.array([-300.0, 0])}, low, sig, names)
    assert (wd, mt) == (0, 0) and np.isnan(sc)


def test_signature_score_weights_classes_by_calls():
    names = ["LTR:CRM", "gene"]
    sig = pd.DataFrame(dict(context=["CHG", "CHG"], cls=["LTR:CRM", "gene"], mean_dev=[-0.2, 0.06],
                            direction=["down", "up"]))
    wd, mt, sc = candidates.signature_match({"CHG": np.array([-150.0, 0.0])}, {"CHG": np.array([1000.0, 1000.0])},
                                            sig, names)
    assert (wd, mt) == (2, 1) and sc == pytest.approx(0.03 / 0.0436)      # gene alone fails, the pattern holds


def test_median_of_others():
    s = np.array([[1.0, np.nan, 1.0, np.nan], [3.0, 2.0, np.nan, np.nan], [5.0, np.nan, np.nan, np.nan]])
    m, k = median_of_others(s)
    assert m[:3].tolist() == [3.0, 2.0, 1.0] and np.isnan(m[3]) and k.tolist() == [3, 1, 1, 0]


def test_pooled_at_keeps_windows_with_data_away_from_anchor_and_ends():
    num = np.zeros((2, 100))
    vok = np.full((2, 100), float(params.MIN_CALLS))
    vok[1, :] = 0
    vok[1, [5, 30, 50, 55, 90]] = params.MIN_CALLS
    p = pooled_at(num, vok, a=50, sep=10)
    assert p.anchor_vok.tolist() == [params.MIN_CALLS, params.MIN_CALLS, 2 * params.MIN_CALLS]
    assert len(p.elsewhere[0][1]) == len([i for i in range(10, 90) if abs(i - 50) >= 10])
    assert p.elsewhere[1][1].tolist() == [params.MIN_CALLS]              # only bin 30 qualifies


def test_fewer_than_three_samples_have_no_plant_specific_score(prep):
    two = prep["samples"][:2]
    sc = scan(prep["bins"], two, prep["chroms"], len(prep["cmap"].names), prep["geo"], {}, "down")
    assert sc.picks and all(np.isnan(p.ps) and np.isnan(p.ps_p) for p in sc.picks)


def test_scan_without_direction_collects_tracks_only(prep):
    sc = scan(prep["bins"], prep["samples"], prep["chroms"], len(prep["cmap"].names), prep["geo"], {}, None)
    assert sc.picks == [] and sc.at_anchor == [] and sc.pooled == {}
    y = sc.adj["S1"]["chr1"]
    assert len(y) == 60_000_000 // params.BIN and np.isfinite(y).mean() > 0.95
