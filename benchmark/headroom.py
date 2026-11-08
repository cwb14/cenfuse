#!/usr/bin/env python3
"""Headroom of the integration: on every arm with truth, is each line's own locus right where the final call is not
(wrong, low or none)?

  python -m benchmark.headroom ROOT JOBS OUT.tsv [--manifest truth.tsv]
Reads ROOT/<arm>/<id>/<g>.centromeres.tsv of every job (holocentric controls skipped); truth and tolerance by the
frozen scorer (benchmark/score.py). Writes one row per chromosome with R/W/. per line, and prints the counts.
Port of v12 bench/headroom.py."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from benchmark import paths
from benchmark import score as S

LINES = dict(final=("start", "end"), seq=("seq_start", "seq_end"), alt=("seq_alt_start", "seq_alt_end"),
             te=("te_start", "te_end"), at=("at_start", "at_end"), anch=("hic_anch_peak_start", "hic_anch_peak_end"),
             denovo=("hic_denovo_peak_start", "hic_denovo_peak_end"), fold=("fold_core_start", "fold_core_end"),
             align=("align_core_start", "align_core_end"), cis=("cis_peak_start", "cis_peak_end"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path, help="release root: ROOT/<arm>/<id>")
    ap.add_argument("jobs", help="job table (arm, id, ...)")
    ap.add_argument("out", help="output table")
    ap.add_argument("--manifest", type=Path, default=paths.MANIFEST, help="truth manifest [benchmark/truth.tsv]")
    a = ap.parse_args(argv)
    T, warned = S.load_manifest(a.manifest), set()
    out = []
    for p in S.jobs_of(a.jobs):
        arm, gid = p[0], p[1]
        run = a.root / arm / gid
        for cen in sorted(run.glob("*.centromeres.tsv")):
            g = cen.name[: -len(".centromeres.tsv")]
            if g not in T or T[g].rule == "holo" or not S.usable(T[g], warned):
                continue
            t = T[g]
            truth = S.load_truth(t)
            tol = S.tolerance(t, S.hic_bin_of(S.run_summary(run, g)))
            strict = t.rule in S.STRICT
            c = pd.read_csv(cen, sep="\t", dtype=str, keep_default_na=False)
            for r in c.to_dict("records"):
                ch = r["chrom"]
                if ch not in truth:
                    continue
                row = dict(arm=arm, genome=g, chrom=ch, conf=r["confidence"], flags=r["flags"],
                           seqconf=r.get("seq_confidence", "."), hic_gate=r.get("hic_gate", "."))
                for ln, (x, y) in LINES.items():
                    s, e = r.get(x, "."), r.get(y, ".")
                    if s in (".", "") or e in (".", ""):
                        row[ln] = "."
                        continue
                    s, e = int(s) - 1, int(e)
                    hit = False
                    for ts, te in truth[ch]:
                        gap = max(0, max(s, ts) - min(e, te))
                        hit |= (max(s, ts) < min(e, te)) if strict else gap <= tol
                    row[ln] = "R" if hit else "W"
                out.append(row)
    D = pd.DataFrame(out)
    D.to_csv(a.out, sep="\t", index=False)
    D["fin"] = np.where(D.final == "R", np.where(D.conf.isin(["high", "medium"]), "R_conf", "R_low"),
                        np.where(D.final == "W", np.where(D.conf.isin(["high", "medium"]), "W_conf", "W_low"), "none"))
    print(D.groupby("fin").size())
    bad = D[D.fin.isin(["W_conf", "W_low", "none"])]
    print("\nnot-right final calls where a line is right:")
    for ln in LINES:
        if ln != "final":
            x = bad[bad[ln] == "R"]
            print(f"{ln:7s} {len(x):4d}  ", x.groupby(["arm", "genome"]).size().to_dict())
    print("\nnot-right final calls where NO line is right:",
          int((bad[[x for x in LINES if x != "final"]] != "R").all(axis=1).sum()))
    low = D[D.fin == "R_low"]
    print("\nright low calls: which other lines also right")
    for ln in LINES:
        if ln != "final":
            x = low[low[ln] == "R"]
            print(f"{ln:7s} {len(x):4d}  ", x.groupby(["arm", "genome"]).size().to_dict())


if __name__ == "__main__":
    main()
