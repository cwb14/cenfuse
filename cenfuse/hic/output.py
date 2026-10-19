"""Writers: BED9 calls (also the anchor input of the companion tool cenmeth), per-chromosome TSV, bedGraph tracks,
QC table."""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import __version__

CALL_COLS = ["chrom", "length", "callable", "bin_size", "peak_start", "peak_end", "core_start", "core_end",
             "peak_over_median", "boot_lo", "boot_hi", "boot_in_core", "lopo_min", "lopo_max", "lopo_stable",
             "null_p", "confidence", "flags", "sat_mid", "peak_minus_sat_Mb", "rounds", "converged"]
RGB = {"high": "0,114,178", "medium": "86,180,233", "low": "190,190,190"}
# Statistic columns: written with 4 significant digits. Coordinate columns: written as integers (bp).
STAT_COLS = ["peak_over_median", "boot_in_core", "null_p", "peak_minus_sat_Mb"]
COORD_COLS = ["length", "bin_size", "peak_start", "peak_end", "core_start", "core_end", "boot_lo", "boot_hi",
              "lopo_min", "lopo_max", "sat_mid"]


def as_finite(v):
    """Float value of v, or None when v is missing (None, empty, NaN, inf or non-numeric)."""
    if v is None or (isinstance(v, str) and v.strip() == ""):
        return None
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if np.isfinite(x) else None


def _int_text(v):
    x = as_finite(v)
    return "" if x is None else str(int(round(x)))


def _stat_text(v):
    x = as_finite(v)
    return "" if x is None else f"{x:.4g}"


def write_bed(rows, path, bin_size):
    with open(path, "w") as f:
        f.write(f"#cenhic {__version__} bin_size={bin_size}\n")
        for r in rows:
            label = r["confidence"]
            if not isinstance(label, str) or label not in RGB:
                continue
            coords = [as_finite(r[c]) for c in ("core_start", "core_end", "peak_start", "peak_end")]
            if any(c is None for c in coords):
                continue
            ratio = as_finite(r["peak_over_median"])
            score = 0 if ratio is None else min(1000, max(0, int(round(100 * ratio))))
            cs, ce, ps, pe = (int(round(c)) for c in coords)
            f.write("\t".join(str(x) for x in (r["chrom"], cs, ce, f"{r['chrom']}_cen", score,
                                                ".", ps, pe, RGB[label])) + "\n")


def write_tsv(rows, path):
    df = pd.DataFrame(rows, columns=CALL_COLS)
    for c in COORD_COLS:
        df[c] = df[c].map(_int_text)
    for c in STAT_COLS:
        df[c] = df[c].map(_stat_text)
    df.to_csv(path, sep="\t", index=False)


def write_bedgraph(values, bins, members, path, name):
    with open(path, "w") as f:
        f.write(f'track type=bedGraph name="{name}"\n')
        for k in np.where(members)[0]:
            sl = bins.sl(int(k))
            for s, e, v in zip(bins.start[sl], bins.end[sl], values[sl]):
                if np.isfinite(v):
                    f.write(f"{bins.names[int(k)]}\t{s}\t{e}\t{v:.4g}\n")


def write_qc(sections, path):
    with open(path, "w") as f:
        for name, df in sections.items():
            f.write(f"# {name}\n")
            df.to_csv(f, sep="\t", index=False, float_format="%.4g")
