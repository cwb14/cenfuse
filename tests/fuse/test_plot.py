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
