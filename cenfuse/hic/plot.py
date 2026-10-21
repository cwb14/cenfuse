"""Overview PDF: every callable chromosome to scale with its centromere-cluster score, core and satellite call."""
from __future__ import annotations

import logging

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .output import as_finite  # noqa: E402

BLUE, GREY, INK = "#0072B2", "#8C8C8C", "#1A1A1A"
log = logging.getLogger("cenfuse-hic")


def overview(S_ratio, bins, members, rows, sat, rabl, path):
    ks = [int(k) for k in np.where(members)[0]]
    if not ks:
        log.warning("overview: no member chromosomes; no figure written")
        return
    member_vals = np.concatenate([S_ratio[bins.sl(k)] for k in ks]).astype(float)
    finite = member_vals[np.isfinite(member_vals)]
    if finite.size == 0:
        log.warning("overview: no finite score on member chromosomes; no figure written")
        return
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
                         "font.size": 6, "pdf.fonttype": 42})
    Lmax = bins.lengths[ks].max() / 1e6
    fig, axs = plt.subplots(len(ks), 1, figsize=(7.1, 0.55 * len(ks) + 0.8), squeeze=False)
    byname = {r["chrom"]: r for r in rows}
    ymax = finite.max()
    for ax, k in zip(axs[:, 0], ks):
        sl = bins.sl(k)
        x = (bins.start[sl] + bins.end[sl]) / 2e6
        ax.plot(x, S_ratio[sl], color=BLUE, lw=0.7)
        ax.axhline(1, color="#D0D0D0", lw=0.4)
        r = byname.get(bins.names[k])
        if r and as_finite(r["core_start"]) is not None and as_finite(r["core_end"]) is not None:
            ax.axvspan(float(r["core_start"]) / 1e6, float(r["core_end"]) / 1e6, color="#E3EEF7", lw=0, zorder=0)
            ax.text(Lmax * 1.005, ymax * 0.5, f"{r['confidence']} {r['flags']}".strip(), fontsize=5, va="center")
        if k in sat:
            ax.plot([sat[k] / 1e6], [0], marker="^", color=INK, ms=3, clip_on=False)
        ax.set_xlim(0, Lmax)
        ax.set_ylim(0, ymax * 1.05)
        ax.set_ylabel(bins.names[k], rotation=0, ha="right", va="center", fontsize=6)
        ax.tick_params(labelsize=5, length=1.5)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axs[-1, 0].set_xlabel("Position (Mb); line = centromere-cluster score (x chromosome median); band = core; "
                          "triangle = satellite centromere", fontsize=6)
    fig.text(0.01, 0.995, f"Rabl test: {rabl['negative']}/{rabl['n']} chromosomes, P = {rabl['p']:.2g}"
             f"{'' if rabl['passed'] else ' (not passed)'}", fontsize=6, va="top")
    fig.tight_layout(rect=(0, 0, 0.93, 0.98))
    fig.savefig(path)
    plt.close(fig)
