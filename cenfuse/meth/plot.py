"""Overview PDF: every analysed chromosome to scale, one adjusted CHG track per sample, candidates as boxes
(blue = centromere_like, grey outline = dip), anchors as a light band with a triangle. Several pages when needed."""
from __future__ import annotations

import logging
import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

from . import params  # noqa: E402

BLUE, GREY, INK, BAND = "#0072B2", "#8C8C8C", "#1A1A1A", "#E3EEF7"
PER_PAGE = 8
log = logging.getLogger("cenmeth")


def overview(adj: dict, lengths: dict, cand, anchors: dict, footer: str, path) -> None:
    """adj: sample -> chrom -> track; cand: candidates table; anchors: chrom -> dict(start, end, mid)."""
    samples = list(adj)
    chroms = list(lengths)
    vals = np.concatenate([adj[s][c] for s in samples for c in chroms]) if samples and chroms else np.zeros(0)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        log.warning("overview: no finite adjusted values; no figure written")
        return
    scale = 0.45 / max(np.quantile(np.abs(vals), 0.999), 1e-6)
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
                         "font.size": 6, "pdf.fonttype": 42})
    Lmax = max(lengths.values()) / 1e6
    ns = len(samples)
    with PdfPages(path) as pdf:
        for p0 in range(0, len(chroms), PER_PAGE):
            page = chroms[p0:p0 + PER_PAGE]
            fig, axs = plt.subplots(len(page), 1, figsize=(7.1, 0.16 * ns * len(page) + 0.35 * len(page) + 0.6),
                                    squeeze=False)
            for ax, c in zip(axs[:, 0], page):
                x = (np.arange(len(adj[samples[0]][c])) + 0.5) * params.BIN / 1e6
                a = anchors.get(c)
                if a is not None:
                    ax.axvspan(a["start"] / 1e6, a["end"] / 1e6, color=BAND, lw=0, zorder=0)
                    ax.plot([a["mid"] / 1e6], [0.6], marker="v", color=INK, ms=3, clip_on=False)
                for i, s in enumerate(samples):
                    base = -i
                    ax.axhline(base, color="#D9D9D9", lw=0.3, zorder=1)
                    ax.plot(x, base + scale * adj[s][c], color=INK, lw=0.4, zorder=2)
                    if ax is axs[0, 0]:
                        ax.text(Lmax * 1.005, base, s, fontsize=5, va="center")
                    if cand is None or cand.empty:
                        continue
                    for r in cand[(cand["sample"] == s) & (cand.chrom == c)].itertuples():
                        like = r.label == "centromere_like"
                        ax.add_patch(Rectangle((r.dip_start / 1e6, base - 0.45), (r.dip_end - r.dip_start) / 1e6, 0.9,
                                               facecolor=BLUE if like else "none", edgecolor=BLUE if like else GREY,
                                               lw=0.5, alpha=0.6 if like else 1.0, zorder=3))
                ax.set_xlim(0, Lmax)
                ax.set_ylim(-ns + 0.4, 0.7)
                ax.set_yticks([])
                ax.set_ylabel(c, rotation=0, ha="right", va="center", fontsize=6)
                ax.tick_params(labelsize=5, length=1.5)
                for sp in ("top", "right", "left"):
                    ax.spines[sp].set_visible(False)
            handles = [Patch(facecolor=BLUE, edgecolor=BLUE, alpha=0.6, label="centromere_like"),
                       Patch(facecolor="none", edgecolor=GREY, label="dip"),
                       Patch(facecolor=BAND, edgecolor="none", label="anchor"),
                       Line2D([], [], color=INK, lw=0.6, label="adjusted CHG (one row per sample)")]
            axs[0, 0].legend(handles=handles, loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=4, frameon=False,
                             fontsize=5, borderaxespad=0.2)
            axs[-1, 0].set_xlabel("Position (Mb)", fontsize=6)
            fig.text(0.01, 0.005, textwrap.fill(footer, 190), fontsize=5, va="bottom")
            fig.tight_layout(rect=(0, 0.045, 0.95, 1))
            pdf.savefig(fig)
            plt.close(fig)
