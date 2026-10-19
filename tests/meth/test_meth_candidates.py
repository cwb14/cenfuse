import pytest
from meth_common import run_scan

from cenfuse.meth import calibrate, candidates, output, params, track
from cenfuse.meth.scan import Pick


@pytest.fixture(scope="module")
def calibrated(prep, syn):
    anchors, _, cal = run_scan(prep, syn["anchors"], None)
    direction, _ = calibrate.learned_direction(cal)
    _, sc, _ = run_scan(prep, syn["anchors"], direction)
    sig = calibrate.signature(cal)
    cand = candidates.table(sc.picks, sig, prep["cmap"].names, anchors, prep["geo"], prep["D"], prep["chroms"])
    return anchors, sc, sig, cand


def test_pick_in_last_bin_is_clipped_to_the_chromosome(tmp_path):
    b = params.BIN
    L = 10 * b + 7_000                                   # last bin (10) is 7 kb long
    geo = track.geometry(0.65e6)
    p = Pick("S1", "c1", 1, 10, -0.1, 8, 10, -3.0)
    cand = candidates.table([p], calibrate.signature(None), ["other"], {}, geo, 0.65e6, {"c1": L})
    r = cand.iloc[0]
    assert r.window_end <= L and r.dip_end <= L and r.peak < L
    assert r.window_end == L and r.dip_end == L and r.peak == (10 * b + L) // 2
    assert r.width == r.dip_end - r.dip_start
    out = tmp_path / "c.bed"
    output.write_bed(cand, "#h", out)
    f = out.read_text().splitlines()[1].split("\t")
    s, e, ts, te = int(f[1]), int(f[2]), int(f[6]), int(f[7])
    assert 0 <= s <= ts < te <= e <= L
    assert bool(r.near_end)


def test_near_end_flags_picks_within_sep_of_either_end():
    geo = track.geometry(0.65e6)
    nb = 4 * geo.sep                                     # chromosome of nb full bins
    bins_ = [geo.sep - 1, geo.sep, nb - 1 - geo.sep, nb - geo.sep]
    picks = [Pick("S1", "c1", i + 1, j, -0.1, j, j, -3.0) for i, j in enumerate(bins_)]
    cand = candidates.table(picks, calibrate.signature(None), ["other"], {}, geo, 0.65e6, {"c1": nb * params.BIN})
    assert cand.near_end.tolist() == [True, False, False, True]


def test_planted_centromeres_rank_first_and_centromere_like(calibrated, syn):
    _, _, _, cand = calibrated
    top = cand[cand["rank"] == 1].set_index(["sample", "chrom"])
    for s in ("S1", "S2", "S3"):
        for c, x in syn["cen"].items():
            want = syn["shift"][2] if (s, c) == syn["shift"][:2] else x
            r = top.loc[(s, c)]
            assert abs(r.peak - want) <= 1e6, (s, c, r.peak)
            assert r.label == "centromere_like", (s, c)


def test_true_negative_and_decoy_are_not_centromere_like(calibrated, syn):
    _, _, _, cand = calibrated
    assert (cand[cand.chrom == "chr6"].label == "dip").all()
    c, a, b = syn["decoy"]
    dec = cand[(cand.chrom == c) & (cand.dip_end > a) & (cand.dip_start < b)]
    assert dec.empty


def test_shifted_centromere_is_plant_specific(calibrated, syn):
    _, _, _, cand = calibrated
    s, c, x = syn["shift"]
    r = cand[(cand["sample"] == s) & (cand.chrom == c) & ((cand.peak - x).abs() <= 1e6)].iloc[0]
    assert r.label == "centromere_like" and r.sample_specific and r.plant_specific < 0


def test_duplicated_segment_is_filtered(calibrated, syn):
    _, _, _, cand = calibrated
    s, c, a, b = syn["dup"]
    near = cand[(cand["sample"] == s) & (cand.chrom == c) & (cand.peak > a - 5e5) & (cand.peak < b + 5e5)]
    assert (near.value > -0.05).all()


def test_anchored_report(calibrated, prep):
    anchors, sc, sig, cand = calibrated
    rep = candidates.anchored(sc, anchors, cand, "down", sig, prep["cmap"].names)
    assert len(rep) == len(anchors) * 3
    r = rep[(rep["sample"] == "S1") & (rep.chrom == "chr1")].iloc[0]
    assert r.candidate_rank == 1 and r.candidate_label == "centromere_like" and r.window_rank <= 5
    assert rep.sites_vs_best.between(0, 1).all()


def test_uncalibrated_labels_are_dip(prep):
    _, sc, _ = run_scan(prep, None, "down")
    cand = candidates.table(sc.picks, calibrate.signature(None), prep["cmap"].names, {}, prep["geo"], prep["D"],
                            prep["chroms"])
    assert set(cand.label) == {"dip"} and cand.anchor.eq("").all()
