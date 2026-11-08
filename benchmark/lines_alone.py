#!/usr/bin/env python3
"""Each line alone vs the integrated call, on the genomes with truth and Hi-C.

  python -m benchmark.lines_alone OUT.tsv [--runs ROOT] [--manifest truth.tsv]
Every line is scored within 4 Hi-C bins of the truth (v3 bench/score.py's `right_near`, NEAR_BINS = 4: Hi-C's
resolution), and by >= 1-bp overlap (`right_overlap`, the mammal benchmark's rule). Intervals: integrated (start/end),
satellite line (seq_start/seq_end), prior-free Hi-C (hic_denovo_peak_*), anchored Hi-C (hic_anch_peak_*). Truth: the
manifest (mammals: tier FW). Runs: the FULL arm of the 9 T2T mammals with CENP-A truth, the SH arm of rice, barley,
potato and horse EquCab3. Port of v12 bench/rel/lines_alone.py."""
import argparse
import sys
from pathlib import Path

import pandas as pd

from benchmark import paths
from benchmark import score as S

MAM_TRUTH = ["chm13", "mPanTro3", "mPanPan1", "mGorGor1", "mPonAbe1", "mPonPyg2", "MFA8", "Sscrofa", "horse1"]
JOBS = [("FULL", g) for g in MAM_TRUTH] + [("SH", g) for g in ("rice", "barley", "potato", "horse")]
LINES = {"integrated": ("start", "end"), "satellite line": ("seq_start", "seq_end"),
         "prior-free Hi-C": ("hic_denovo_peak_start", "hic_denovo_peak_end"),
         "anchored Hi-C": ("hic_anch_peak_start", "hic_anch_peak_end")}


def iv(r, a, b):
    try:
        s, e = int(float(r[a])), int(float(r[b]))
    except (ValueError, KeyError):
        return None
    return (s - 1, e) if s > 0 and e > 0 else None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", help="output table")
    ap.add_argument("--runs", type=Path, default=paths.RUNS, help="release root [paths.RUNS]")
    ap.add_argument("--manifest", type=Path, default=paths.MANIFEST, help="truth manifest [benchmark/truth.tsv]")
    a = ap.parse_args(argv)
    T, warned = S.load_manifest(a.manifest), set()
    rows = []
    for arm, g in JOBS:
        if not S.usable(T[g], warned):
            continue
        run = a.runs / arm / g
        t = pd.read_csv(run / f"{g}.centromeres.tsv", sep="\t", dtype=str, keep_default_na=False)
        summ = pd.read_csv(run / "cenfuse_summary.tsv", sep="\t", dtype=str).iloc[0]
        bin_bp = int(float(summ["hic_bin"]))
        truth = S.load_truth(T[g])
        for line, (x, y) in LINES.items():
            n = near = strict = 0
            for r in t.to_dict("records"):
                if r["chrom"] not in truth:
                    continue
                n += 1
                c = iv(r, x, y)
                if c is None:
                    continue
                gaps = [max(0, max(c[0], ts) - min(c[1], te)) for ts, te in truth[r["chrom"]]]
                near += min(gaps) <= 4 * bin_bp
                strict += any(max(c[0], ts) < min(c[1], te) for ts, te in truth[r["chrom"]])
            rows.append(dict(genome=g, arm=arm, line=line, chromosomes=n, right_near=near, right_overlap=strict,
                             hic_bin=bin_bp))
    if not rows:
        sys.exit("lines_alone: no run with truth")
    r = pd.DataFrame(rows)
    r.to_csv(a.out, sep="\t", index=False)
    for col in ("right_near", "right_overlap"):
        p = r.pivot_table(index="genome", columns="line", values=col, aggfunc="first")
        p["chromosomes"] = r.groupby("genome").chromosomes.first()
        print(col)
        print(p.to_string())
        print(r[r.genome.isin(MAM_TRUTH)].groupby("line")[[col, "chromosomes"]].sum().to_string())


if __name__ == "__main__":
    main()
