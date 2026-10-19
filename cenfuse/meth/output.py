"""Output files: candidates BED9 and TSV, anchored report, signature, per-sample adjusted bedGraph, QC table."""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import params

RGB = {"centromere_like": "0,114,178", "dip": "140,140,140"}


def write_tsv(df: pd.DataFrame, path) -> None:
    df.to_csv(path, sep="\t", index=False, float_format="%.10g", na_rep="")


def write_bed(cand: pd.DataFrame, header: str, path) -> None:
    """BED9: span = half-depth dip united with the window; name sample:chrom:rank:label; score = |value| x 1000
    (capped at 1000); thick = sliding window at the extreme; colour by label."""
    with open(path, "w") as f:
        f.write(f"{header}\n")
        if cand.empty:
            return
        t = cand.assign(bed_start=np.minimum(cand.dip_start, cand.window_start).astype(np.int64),
                        bed_end=np.maximum(cand.dip_end, cand.window_end).astype(np.int64))
        for r in t.sort_values(["chrom", "bed_start", "sample"]).itertuples():
            score = int(min(1000, round(abs(r.value) * 1000)))
            f.write(f"{r.chrom}\t{r.bed_start}\t{r.bed_end}\t{r.sample}:{r.chrom}:{r.rank}:{r.label}\t{score}\t.\t"
                    f"{int(r.window_start)}\t{int(r.window_end)}\t{RGB[r.label]}\n")


def write_bedgraph(adj: dict, lengths: dict, path) -> None:
    """One sample's adjusted CHG track: bin start, end (clipped to the chromosome), value; NaN bins left out."""
    with open(path, "w") as f:
        for c, y in adj.items():
            j = np.flatnonzero(np.isfinite(y))
            s = j * params.BIN
            e = np.minimum(s + params.BIN, lengths[c])
            for a, b, v in zip(s, e, y[j]):
                f.write(f"{c}\t{a}\t{b}\t{v:.5g}\n")


def qc_rows(run: dict, depth_stats, tiers: dict) -> pd.DataFrame:
    """Long table (section, key, value): run-level settings, then per sample depth, dropped bins and observed
    cytosines per class tier (with the share of the best-covered sample)."""
    rows = [("run", k, v) for k, v in run.items()]
    best = {t: max(tiers[s].get(t, 0) for s in tiers) for t in ("satellite", "transposon", "gene", "other")}
    for d in depth_stats:
        sec = f"sample:{d.sample}"
        rows += [(sec, "median_depth", round(d.median_depth, 2)), (sec, "bins_with_data", d.bins_with_data),
                 (sec, "bins_dropped_depth", d.dropped_depth), (sec, "bins_dropped_excluded", d.dropped_excluded),
                 (sec, "low_depth", bool(d.median_depth < params.LOW_DEPTH))]
        for t, b in best.items():
            n = tiers[d.sample].get(t, 0)
            rows += [(sec, f"{t}_sites", int(n)), (sec, f"{t}_sites_vs_best", round(n / b, 4) if b else "")]
    return pd.DataFrame(rows, columns=["section", "key", "value"])
