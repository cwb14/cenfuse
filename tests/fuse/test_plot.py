"""Overview plot: panel C (each line's call relative to the sequence call) scales to the calls drawn."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cenfuse import plot


def test_offsets_axis_covers_sequence_only_calls():
    """FASTA-only run with every chromosome called: only rectangle patches are drawn, which do not rescale the axis
    by themselves (mouse and maize kept 0-1 Mb, showing the upper half of each call)."""
    rows = [dict(chrom=f"c{i}", seq=(1_000_000, 9_000_000), call=(1_000_000, 9_000_000), te=None, hic=None)
            for i in range(3)]
    fig, ax = plt.subplots()
    plot._offsets(ax, rows)
    fig.canvas.draw()
    lo, hi = ax.get_ylim()
    plt.close(fig)
    assert lo <= -4 and hi >= 4


def _fig(rdna):
    import pandas as pd
    from types import SimpleNamespace
    calls = pd.DataFrame([dict(chrom="c1", chrom_len=10_000_000, start="4000001", end="5000000", seq_start="4000001",
                               seq_end="5000000", relation="sequence_only", confidence="high", call_source="sequence",
                               support="satellite", against=".", te_start=".", te_end=".", hic_anch_core_start=".",
                               hic_anch_core_end=".", hic_anch_peak_start=".", hic_anch_peak_end=".")])
    arrays = pd.DataFrame(dict(chrom=["c1"], start=[4_100_000], end=[4_900_000], role=["centromeric"]))
    return dict(genome=SimpleNamespace(name="g", lines=["satellite"]), seq=SimpleNamespace(arrays=arrays), tel={},
                hic=None, meth=None, calls=calls, telo={}, plants=None, out="x", rdna=rdna)


def test_overview_marks_45s_and_5s_rdna_arrays():
    """45S arrays (NORs) as filled and 5S arrays as open green diamonds at their midpoints, in the key and the legend;
    a genome without either says so."""
    from matplotlib.colors import to_hex
    f = _fig(dict(nors=[("c1", 100_000, 400_000, 30, 9_100)], five_s=[("c1", 8_000_000, 8_200_000, 80, 500)]))
    fig = plot._pages(f)[0]
    marks = sorted((float(x.get_xdata()[0]), to_hex(x.get_markerfacecolor())) for x in fig.axes[0].lines
                   if x.get_marker() == "D")
    assert marks == [(250_000.0, to_hex(plot.ROLE_COLOR["rDNA"])), (8_100_000.0, "#ffffff")]
    labels = [t.get_text() for t in fig.legends[0].get_texts()]
    assert "45S rDNA array (NOR)" in labels and "5S rDNA array" in labels
    plt.close(fig)
    assert "1 of 45S" in plot.legends(f)[1] and "1 of 5S" in plot.legends(f)[1]
    f = _fig(dict(nors=[], five_s=[]))
    fig = plot._pages(f)[0]
    assert not [x for x in fig.axes[0].lines if x.get_marker() == "D"]
    assert not {"45S rDNA array (NOR)", "5S rDNA array"} & {t.get_text() for t in fig.legends[0].get_texts()}
    plt.close(fig)
    assert "No 45S or 5S rDNA array was found" in plot.legends(f)[1]
