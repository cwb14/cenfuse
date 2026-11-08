#!/usr/bin/env python3
"""How close the right calls are to their truth domain (future.md #37: overlap is lenient for wide calls).

  python -m benchmark.precision CHROM.tsv OUT.md [--tsv OUT.tsv]
CHROM.tsv: the scorer's per-chromosome table (benchmark/score.py: release<TAG>.chrom.tsv). For every right call of the real arms
(the randomized-Hi-C controls NEG, NEG2 and TNEG left out; truth by agreement with another caller and the holocentric controls left
out, since neither is a domain): the call's width, the domain's width, the distance between their midpoints, whether the call's
midpoint lies inside the domain, and the share of the call that lies on the domain (0 for a call right within the tolerance
without overlap). Writes the markdown table (one row per truth class, all right calls and the confident ones) and, with --tsv,
one row per call."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

CONTROLS = {"NEG", "NEG2", "TNEG"}
CLASSES = [("animals: CENP-A, CenSat active array, minor satellite", {"CENP-A", "CENP-A windows", "CenSat active", "minor satellite"}),
           ("plants and algae: CENH3", {"CENH3", "CENH3 (published)"}),
           ("fungi and *Plasmodium*: mapped centromeres", {"CEN"})]


def calls(chrom: pd.DataFrame) -> pd.DataFrame:
    kinds = set().union(*(k for _, k in CLASSES))                       # agreement and holocentric: no domain
    d = chrom[~chrom.arm.isin(CONTROLS) & chrom.kind.isin(kinds) & chrom.outcome.astype(str).str.startswith("right")].copy()
    d = d[d.truth.astype(str).str.contains("-")]
    s, e = d.start.astype(int) - 1, d.end.astype(int)                  # calls: 1-based inclusive
    ts = d.truth.str.split("-").str[0].astype(int)
    te = d.truth.str.split("-").str[1].astype(int)                       # truth: half-open, as the scorer wrote it
    d["call_bp"], d["truth_bp"] = e - s, te - ts
    d["mid_dist"] = ((s + e) / 2 - (ts + te) / 2).abs()
    d["mid_inside"] = ((s + e) / 2 >= ts) & ((s + e) / 2 <= te)
    d["on_domain"] = (np.minimum(e, te) - np.maximum(s, ts)).clip(lower=0) / (e - s)
    d["confident"] = d.confidence.isin(["high", "medium"])
    return d


def _mb(x: float) -> str:
    return f"{x / 1e6:.2f} Mb" if x >= 1e5 else f"{x / 1e3:.0f} kb"


def table(d: pd.DataFrame) -> str:
    out = ["| truth | calls | right calls | median call width | median domain width | median midpoint distance (90th percentile) | "
           "midpoint inside the domain | median share of the call on the domain |",
           "|---|---|---|---|---|---|---|---|"]
    for name, kinds in CLASSES:
        x = d[d.kind.isin(kinds)]
        for label, y in (("all right", x), ("`high`, `medium`", x[x.confident])):
            if not len(y):
                continue
            out.append(f"| {name} | {label} | {len(y):,} | {_mb(y.call_bp.median())} | {_mb(y.truth_bp.median())} | "
                       f"{_mb(y.mid_dist.median())} ({_mb(y.mid_dist.quantile(0.9))}) | {100 * y.mid_inside.mean():.0f}% | "
                       f"{100 * y.on_domain.median():.0f}% |")
    return "\n".join(out) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("chrom", type=Path, help="the scorer's per-chromosome table (release<TAG>.chrom.tsv)")
    ap.add_argument("out", type=Path, help="markdown table to write")
    ap.add_argument("--tsv", type=Path, help="also one row per right call")
    a = ap.parse_args(argv)
    d = calls(pd.read_csv(a.chrom, sep="\t", dtype=str, keep_default_na=False))
    a.out.write_text(table(d))
    if a.tsv:
        d.to_csv(a.tsv, sep="\t", index=False)
    print(a.out.read_text(), end="")


if __name__ == "__main__":
    main()
