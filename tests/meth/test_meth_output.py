import numpy as np
import pandas as pd

from cenfuse.meth import bins, output, params, plot


def _cand():
    return pd.DataFrame(dict(sample=["S1", "S2"], chrom=["c1", "c1"], rank=[1, 2],
                             label=["centromere_like", "dip"], peak=[1_025_000, 3_025_000],
                             window_start=[950_000, 2_950_000], window_end=[1_100_000, 3_100_000],
                             dip_start=[1_000_000, 2_900_000], dip_end=[1_050_000, 3_150_000],
                             value=[-0.12, -2.0]))


def test_bed9_is_valid(tmp_path):
    p = tmp_path / "x.candidates.bed"
    output.write_bed(_cand(), "#cenmeth test", p)
    lines = p.read_text().splitlines()
    assert lines[0] == "#cenmeth test"
    for ln in lines[1:]:
        f = ln.split("\t")
        assert len(f) == 9
        s, e, sc, ts, te = int(f[1]), int(f[2]), int(f[4]), int(f[6]), int(f[7])
        assert s < e and s <= ts < te <= e and 0 <= sc <= 1000
        assert f[3].count(":") == 3 and f[5] == "."
    assert lines[1].split("\t")[8] == "0,114,178" and lines[2].split("\t")[4] == "1000"


def test_empty_candidates_write_header_only(tmp_path):
    p = tmp_path / "e.bed"
    output.write_bed(pd.DataFrame(), "#h", p)
    assert p.read_text() == "#h\n"


def test_bedgraph_skips_nan_and_clips_end(tmp_path):
    p = tmp_path / "x.bedgraph"
    output.write_bedgraph({"c1": np.array([0.1, np.nan, -0.2])}, {"c1": 2 * params.BIN + 10}, p)
    rows = [ln.split("\t") for ln in p.read_text().splitlines()]
    assert [r[1] for r in rows] == ["0", str(2 * params.BIN)] and rows[-1][2] == str(2 * params.BIN + 10)


def test_qc_rows():
    st = [bins.DepthStats("S1", 20.0, 100, 2, 3), bins.DepthStats("S2", 4.0, 100, 0, 0)]
    q = output.qc_rows({"D": 1}, st, {"S1": {"satellite": 10, "transposon": 50}, "S2": {"satellite": 5}})
    d = {(r.section, r.key): r.value for r in q.itertuples()}
    assert d[("run", "D")] == 1 and d[("sample:S2", "satellite_sites_vs_best")] == 0.5
    assert d[("sample:S2", "low_depth")] is True and d[("sample:S1", "bins_dropped_excluded")] == 3


def test_overview_writes_pdf(tmp_path):
    adj = {"S1": {"c1": np.linspace(-0.1, 0.1, 80)}, "S2": {"c1": np.zeros(80)}}
    p = tmp_path / "o.pdf"
    plot.overview(adj, {"c1": 80 * params.BIN}, _cand(), {"c1": dict(start=900_000, end=1_100_000, mid=1_000_000)},
                  "footer", p)
    assert p.read_bytes()[:4] == b"%PDF"


def test_bed_sorted_on_written_start_and_integer_coords(tmp_path):
    c = pd.DataFrame(dict(sample=["S1", "S1"], chrom=["c1", "c1"], rank=[1, 2],
                          label=["dip", "dip"], peak=[1.0, 1.0],
                          window_start=[1_000_000.0, 900_000.0], window_end=[1_100_000.0, 1_250_000.0],
                          dip_start=[1_000_000.0, 1_050_000.0], dip_end=[1_100_000.0, 1_200_000.0],
                          value=[-1.0, -1.0]))
    p = tmp_path / "order.bed"
    output.write_bed(c, "#h", p)
    rows = [ln.split("\t") for ln in p.read_text().splitlines()[1:]]
    starts = [int(r[1]) for r in rows]
    assert starts == sorted(starts)
    assert all("." not in r[1] and "." not in r[2] and "." not in r[6] for r in rows)


def test_overview_paginates_and_embeds_truetype(tmp_path):
    adj = {"S1": {f"c{i}": np.linspace(-0.1, 0.1, 80) for i in range(9)}}
    p = tmp_path / "many.pdf"
    plot.overview(adj, {f"c{i}": 80 * params.BIN for i in range(9)}, pd.DataFrame(), {}, "footer", p)
    raw = p.read_bytes()
    assert raw.count(b"/Type /Page") - raw.count(b"/Type /Pages") == 2
    assert b"/CIDFontType2" in raw and b"/FontFile2" in raw


def test_all_nan_track_writes_no_figure(tmp_path):
    p = tmp_path / "nan.pdf"
    plot.overview({"S1": {"c1": np.full(10, np.nan)}}, {"c1": 10 * params.BIN}, pd.DataFrame(), {}, "f", p)
    assert not p.exists()


def test_tsv_writes_large_floats_without_exponent(tmp_path):
    p = tmp_path / "t.tsv"
    output.write_tsv(pd.DataFrame(dict(sample=["S1", "S2"], anchor_distance=[22_350_000.0, np.nan])), p)
    assert p.read_text().splitlines()[1:] == ["S1\t22350000", "S2\t"]
