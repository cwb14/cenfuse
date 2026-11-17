"""Close-up pages (panel E, cenfuse.closeup): one tile per chromosome at one shared width, the centromere as a band
through every track, each signal's own place as a bar (an arrow at the tile edge when it lies beyond the window), its
role in the call, and the verdict. Synthetic genomes as in test_plot._fig: chromosomes c1..cn, each with a satellite
call and one satellite array under it, LTR elements and Hi-C where a test needs them."""
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text

from cenfuse import closeup
from cenfuse import signals as sg

MM = 1 / 25.4


def _calls(n=1, L=10_000_000, **by_chrom):
    """A centromeres table of n chromosomes c1..cn of length L, each called on its satellite array at 4-5 Mb (1-based
    closed, as the table writes it); by_chrom['c2'] = {...} overrides that chromosome's columns."""
    recs = []
    for i in range(1, n + 1):
        rec = dict(chrom=f"c{i}", chrom_len=str(L), start="4000001", end="5000000", seq_start="4000001",
                   seq_end="5000000", seq_evidence="satellite", relation="sequence_only", confidence="high",
                   call_source="sequence", support="satellite", against=".", rule="sequence", flags=".")
        rec.update(by_chrom.get(f"c{i}", {}))
        recs.append(rec)
    cols = set().union(*recs)
    return pd.DataFrame([{**{c: "." for c in cols}, **r} for r in recs])


def _f(calls, te=False, hic=False):
    """The figure dict `plot.overview` receives, for a centromeres table: one satellite array under each satellite
    place; with te, 40 CRM elements around each call and 20 spread over the chromosome; with hic, a 100-kb Hi-C score
    peaking at each call."""
    rows = calls.to_dict("records")
    arrays = pd.DataFrame([dict(chrom=r["chrom"], start=int(r["seq_start"]) - 1, end=int(r["seq_end"]),
                                role="centromeric") for r in rows if r["seq_start"] != "."],
                          columns=["chrom", "start", "end", "role"])
    f = dict(genome=SimpleNamespace(name="g", lines=["satellite"]), seq=SimpleNamespace(arrays=arrays), tel={},
             hic=None, meth=None, calls=calls, telo={}, plants=None, out="x", rdna={})
    lens = {r["chrom"]: int(r["chrom_len"]) for r in rows}
    if te:
        mids = {}
        for r in rows:
            L = lens[r["chrom"]]
            c = (int(r["start"]) + int(r["end"])) // 2 if r["start"] != "." else L // 2
            m = np.concatenate([np.linspace(c - 400_000, c + 400_000, 40), np.linspace(0, L - 1, 20)])
            mids[r["chrom"]] = np.sort(np.clip(m, 0, L - 1).astype(np.int64))
        f["tel"] = dict(marker=("CRM", {}, True), mids={"CRM": mids}, signature=None,
                        density={c: (100_000, np.bincount(m // 100_000, minlength=lens[c] // 100_000 + 1))
                                 for c, m in mids.items()})
    if hic:
        from cenfuse.hic.bins import Bins
        size = 100_000
        b = Bins(lens, size)
        S = np.ones(b.N)
        for k, r in enumerate(rows):
            if r["start"] != ".":
                S[b.off[k] + (int(r["start"]) // size)] = 4.0
        f["hic"] = dict(bins=b, members=np.ones(len(lens), bool), S=S, denovo_S=S.copy(), bin_size=size)
    return f


def _texts(figs):
    return [t for fig in figs for t in fig.findobj(Text) if t.get_text()]


def _pages(f, rows, keys):
    figs = list(closeup.pages(f, rows, keys))
    return figs


@pytest.fixture(autouse=True)
def _close():
    yield
    plt.close("all")


def test_window_is_shared_and_holds_used_places():
    """One width for every tile, a multiple of 0.5 Mb, wide enough for the call and the places of the signals used or
    pointing elsewhere (10% margin); each window holds them."""
    calls = _calls(3, L=60_000_000,
                   c2=dict(support="satellite,te", te_start="9000001", te_end="9500000"),     # used, 4 Mb right
                   c3=dict(start="20000001", end="20500000", seq_start="20000001", seq_end="20500000",
                           against="te", te_start="15000001", te_end="16000000"))            # elsewhere, 4 Mb left
    rows = sg.rows(calls, None)
    assert rows[1].roles["te"] == "used" and rows[2].roles["te"] == "elsewhere"
    W, wins = closeup.window(rows, ["sat", "te"], [r.L for r in rows])
    assert W >= 1.1 * (9_500_000 - 4_000_000) and W % 500_000 == 0
    assert all(hi - lo == pytest.approx(W) for lo, hi in wins)
    for r, (lo, hi) in zip(rows, wins):
        for iv in (r.call, r.places["te"]):
            if iv:
                assert lo <= iv[0] and iv[1] <= hi, (r.chrom, iv, lo, hi)


def test_window_shrinks_for_small_chromosomes():
    """Yeast-sized chromosomes (1 Mb): the window shrinks to at most 1 Mb, positions read in kb, the scale bar is
    100 kb, and nothing is drawn past a chromosome end."""
    calls = _calls(8, L=1_000_000, c1=dict(start="40001", end="80000", seq_start="40001", seq_end="80000"),
                   **{f"c{i}": dict(start="400001", end="450000", seq_start="400001", seq_end="450000")
                      for i in range(2, 9)})
    f = _f(calls, te=True, hic=True)
    rows = sg.rows(calls, None, 100_000)
    W, wins = closeup.window(rows, ["sat", "te", "hic"], [r.L for r in rows])
    assert W <= 1_000_000
    figs = _pages(f, rows, ["sat", "te", "hic"])
    texts = _texts(figs)
    win = [t.get_text() for t in texts if str(t.get_gid()).startswith("window:")]
    assert len(win) == 8 and all(x.endswith(" kb") for x in win), win
    assert [t.get_text() for t in texts if t.get_gid() == "scale"] == ["100 kb (all close-ups)"]
    # c1's window starts left of the chromosome (its call sits 40-80 kb from the start): nothing drawn below 0 or past L
    ax = next(a for fig in figs for a in fig.axes if a.get_gid() == "tile:c1")
    assert ax.get_xlim()[0] < 0
    L = 1_000_000
    for art in ax.get_children():
        if isinstance(art, Line2D) and len(art.get_xdata()):
            x = np.asarray(art.get_xdata(), float)
            assert x.min() >= 0 and x.max() <= L, art
        elif isinstance(art, Rectangle) and art is not ax.patch:
            assert art.get_x() >= 0 and art.get_x() + art.get_width() <= L, art
        elif isinstance(art, PolyCollection):
            for p in art.get_paths():
                assert p.vertices[:, 0].min() >= 0 and p.vertices[:, 0].max() <= L


def test_far_place_becomes_edge_arrow():
    """A signal's place 30 Mb from the call is an arrow at the tile edge with its position, in the signal's color; the
    window keeps the width it has without that place. A missing place draws no bar, but its role still shows."""
    far = dict(against="te", te_start="35000001", te_end="36000000")
    calls = _calls(2, L=100_000_000, c1=far, c2=dict(support="satellite,te"))      # c2: te used, its place '.'
    rows = sg.rows(calls, None)
    assert rows[0].roles["te"] == "elsewhere" and rows[1].roles["te"] == "used" and rows[1].places["te"] is None
    W, _ = closeup.window(rows, ["sat", "te"], [r.L for r in rows])
    W0, _ = closeup.window(sg.rows(_calls(2, L=100_000_000), None), ["sat", "te"], [100_000_000] * 2)
    assert W == W0
    figs = _pages(_f(calls, te=True), rows, ["sat", "te"])
    texts = _texts(figs)
    arrow = [t for t in texts if t.get_gid() == "arrow:c1:te"]
    assert len(arrow) == 1 and arrow[0].get_text() == "35.5 Mb"
    assert matplotlib.colors.to_hex(arrow[0].get_color()) == sg.SIGNALS["te"].color.lower()
    bars = {a.get_gid() for fig in figs for a in fig.findobj(Line2D) if str(a.get_gid()).startswith("place:")}
    assert "place:c1:sat" in bars and "place:c2:sat" in bars
    assert "place:c1:te" not in bars and "place:c2:te" not in bars             # off the window / no place
    roles = [t.get_text() for t in texts if t.get_gid() == "role:c2:te"]
    assert roles == ["used"]


def test_pages_stay_under_240_mm():
    """60 chromosomes: every close-up page is at most 240 mm tall and 180 mm wide, and every chromosome has one tile."""
    calls = _calls(60)
    rows = sg.rows(calls, None)
    figs = _pages(_f(calls, te=True), rows, ["sat", "te"])
    assert len(figs) > 1
    for fig in figs:
        w, h = fig.get_size_inches() / MM
        assert h <= 240 + 1e-6 and w == pytest.approx(180)
    names = [t.get_text() for t in _texts(figs) if str(t.get_gid()).startswith("name:")]
    assert sorted(names) == sorted(f"c{i}" for i in range(1, 61))


def test_chromosome_without_call():
    """A chromosome without a call (relation none): its tile says so and draws no band; its signals are still drawn."""
    calls = _calls(2, c2=dict(start=".", end=".", relation="none", confidence="none", call_source="none",
                              rule="none", support=".", seq_start="2000001", seq_end="2500000"))
    rows = sg.rows(calls, None)
    figs = _pages(_f(calls, te=True), rows, ["sat", "te"])
    verdicts = {t.get_gid(): t.get_text() for t in _texts(figs) if str(t.get_gid()).startswith("verdict:")}
    assert verdicts["verdict:c2"] == "No centromere called"
    bands = {a.get_gid() for fig in figs for a in fig.findobj(Rectangle) if str(a.get_gid()).startswith("band:")}
    assert "band:c1" in bands and "band:c2" not in bands
    bars = {a.get_gid() for fig in figs for a in fig.findobj(Line2D)}
    assert "place:c2:sat" in bars


def test_fasta_only_has_one_track():
    """A FASTA-only run: the satellite arrays are the one track."""
    calls = _calls(3)
    rows = sg.rows(calls, None)
    keys = sg.present(rows, None, {})
    assert keys == ["sat"]
    figs = _pages(_f(calls), rows, keys)
    tracks = {t.get_text() for t in _texts(figs) if str(t.get_gid()).startswith("track:")}
    assert tracks == {"satellite arrays"}
    assert all(r.verdict == "Satellite array alone" for r in rows)


def test_no_banned_words():
    """No word of cenfuse's jargon on a close-up page: every label, role, verdict and key entry, with every signal's
    track drawn."""
    calls = _calls(4, L=40_000_000,
                   c1=dict(support="satellite,te,hic", rule="hic_agrees", hic_gate="pass", te_clade="CRM",
                           te_start="3500001", te_end="5500000", hic_anch_core_start="4000001",
                           hic_anch_core_end="4500000"),
                   c2=dict(call_source="hic", rule="cis_only", flags="cis_only", support="hic", confidence="low",
                           start="12000001", end="12250000", cis_peak_start="12000001", cis_peak_end="12250000",
                           seq_start="12000001", seq_end="12250000", seq_evidence="long_unit_tandem",
                           against="te", te_start="30000001", te_end="31000000"),
                   c3=dict(seq_alt_start="7000001", seq_alt_end="7500000", at_start="4200001", at_end="4300000"),
                   c4=dict(start=".", end=".", relation="none", confidence="none", call_source="none", rule="none",
                           support="."))
    rows = sg.rows(calls, None, 100_000)
    figs = _pages(_f(calls, te=True, hic=True), rows, list(sg.SIGNALS))
    texts = [t.get_text() for t in _texts(figs)]
    assert any("long tandem repeat" in x for x in texts) and any("second satellite cluster" in x for x in texts)
    hits = [x for x in texts if sg.BANNED.search(x)]
    assert not hits, hits


def test_hic_scale_is_genome_wide_and_clips_with_a_caret():
    """One Hi-C scale for the genome, the 99.5th percentile of the trace over all chromosomes (at least 2x), as on the
    genome pages: a spike above it is clipped and marked by a caret at the track's top, and the key states the top."""
    calls = _calls(3)
    f = _f(calls, hic=True)
    b = f["hic"]["bins"]
    f["hic"]["S"][b.off[b.names.index("c1")] + 45] = 40.0        # one spike beside c1's call (4-5 Mb)
    rows = sg.rows(calls, None, 100_000)
    figs = _pages(f, rows, ["hic", "sat"])
    texts = [t.get_text() for t in _texts(figs)]
    assert any("top 4.0×" in x for x in texts), texts                  # not the spike's 40x
    carets = [a.get_gid() for fig in figs for a in fig.findobj(Line2D) if str(a.get_gid()).startswith("clip:")]
    assert carets == ["clip:c1"]


def test_a_long_chromosome_name_never_reaches_the_confidence_chip():
    """A 60-character chromosome name is shortened in the middle (an ellipsis) so that it and the window's span end
    left of the confidence chip; a short name is drawn whole."""
    long = "scaffold_0001_RagTag_HiC_pilon_polished_haplotype1_v2_chr01A"
    assert len(long) == 60
    calls = _calls(2, c1=dict(confidence="medium"))
    calls.loc[0, "chrom"] = long
    f = _f(calls)
    rows = sg.rows(calls, None)
    figs = _pages(f, rows, ["sat"])
    fig = figs[0]
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for c in (long, "c2"):
        name = next(t for t in _texts(figs) if t.get_gid() == f"name:{c}")
        chip = next(a for a in fig.findobj(Line2D) if a.get_gid() == f"chip:{c}")
        x_mid = fig.transFigure.transform((chip.get_xdata()[0], chip.get_ydata()[0]))[0]
        x_chip = x_mid - chip.get_markersize() / 2 * fig.dpi / 72              # the chip's left edge, pixels
        right = max([name.get_window_extent(renderer).x1] +
                    [t.get_window_extent(renderer).x1 for t in _texts(figs) if t.get_gid() == f"window:{c}"])
        assert right < x_chip, (c, right, x_chip)
        if c == long:
            assert "\u2026" in name.get_text() and name.get_text()[:5] == long[:5] and name.get_text()[-3:] == long[-3:]
        else:
            assert name.get_text() == "c2"
