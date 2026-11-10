#!/usr/bin/env python3
"""The Rabl-axis table for docs/benchmark.md: every run with Hi-C (release arms, randomized controls, P. annua) of a job
table, grouped by nucleus type, how many pass the trans arm alignment, with k/n and P.

  python -m benchmark.rabl_table OUT.md [--runs ROOT] [--jobs JOBS] [--canonical FILE]
Reads ROOT/<arm>/<id>/cenfuse_summary.tsv; writes OUT.md and RABL_RUNS.md (the run count) beside it.
JOBS defaults to the 220 arms of 0.12's own run (paths.JOBS_RABL), as the docs since 0.13 use.
Port of v12 bench/rabl_table.py."""
import argparse
from pathlib import Path

import pandas as pd

from benchmark import paths
from benchmark.md_tables import canonical

SHORT = {"scer": "yeast", "pfal": "*Plasmodium*", "ncra": "*Neurospora*", "spom": "*S. pombe*", "vdah": "*V. dahliae*",
         "ztri": "*Z. tritici*", "cneoA": "*C. neoformans* anaphase", "cneoG1": "*C. neoformans* G1",
         "msym": "*M. sympodialis*", "ctro": "*C. tropicalis*", "celegans": "*C. elegans*", "cnig": "*C. nigoni*",
         "bmori": "*B. mori*", "rbrev": "*R. breviuscula*", "chm13": "CHM13", "annua": "*P. annua*",
         "infirma": "*P. infirma*", "barley": "barley", "colcc": "*Arabidopsis*"}
FUNGI = ["scer", "spom", "pfal", "ncra", "vdah", "ztri", "cneoA", "msym", "ctro"]
CTRL = {"NEG", "NEG2", "TNEG"}
GROUPS = [("fungi and *Plasmodium*", lambda a, g, L: a not in CTRL and g in FUNGI),
          ("*C. neoformans* in G1 (centromeres unclustered)", lambda a, g, L: a not in CTRL and g == "cneoG1"),
          ("barley", lambda a, g, L: a not in CTRL and g == "barley"),
          ("*P. annua* and *P. infirma* (no ChIP truth)", lambda a, g, L: a not in CTRL and g in ("annua", "infirma")),
          ("territorial mammals", lambda a, g, L: a not in CTRL and L == "mammals"),
          ("other plants and algae, chicken", lambda a, g, L: a not in CTRL and L not in ("mammals", "fungi", "protists",
                                                                                         "holocentric")
           and g not in ("barley", "annua", "infirma")),
          ("holocentric", lambda a, g, L: a not in CTRL and L == "holocentric"),
          ("randomized Hi-C", lambda a, g, L: a in CTRL)]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path, help="output markdown file")
    ap.add_argument("--runs", type=Path, default=paths.RUNS, help="release root [paths.RUNS]")
    ap.add_argument("--jobs", default=paths.JOBS_RABL, help="job table of the runs to read [paths.JOBS_RABL]")
    ap.add_argument("--canonical", type=Path, default=paths.CANONICAL, help="genome table [benchmark/canonical.tsv]")
    a = ap.parse_args(argv)
    can = canonical(a.canonical)
    rows = []
    for arm, gid in (x.rstrip("\n").split("\t")[:2] for x in open(a.jobs)):
        f = a.runs / arm / gid / "cenfuse_summary.tsv"
        if not f.exists():
            continue
        for _, r in pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False).iterrows():
            if r.get("align_gate", "") in ("", ".") and r.get("fold_gate", "") in ("", "."):
                continue                                   # no Hi-C in this run
            g = r.genome
            L = ("holocentric" if g in ("celegans", "cnig", "bmori", "rbrev") else
                 can.at[g, "lineage"].split(" ")[0] if g in can.index else "plants")
            rows.append(dict(arm=arm, genome=g, lineage=L, aln=r.get("align_gate", ""), k_n=r.get("align_on_denovo", ""),
                             p=r.get("align_p", "")))
    d = pd.DataFrame(rows)
    used = set()
    out = ["| nuclei | runs | alignment passes | examples (k/n, P) |", "|---|---|---|---|"]
    for name, test in GROUPS:
        x = d[[test(a_, g, l) and i not in used for i, (a_, g, l) in enumerate(zip(d.arm, d.genome, d.lineage))]]
        used |= set(x.index)
        if not len(x):
            continue
        ok = x[x.aln == "pass"]
        ex = ok if len(ok) else x[x.k_n.str.contains("/")].assign(
            _k=lambda y: y.k_n.str.split("/").str[0].astype(int) / y.k_n.str.split("/").str[1].astype(int).clip(lower=1)
        ).sort_values("_k", ascending=False).head(1)
        ex = ex.drop_duplicates("genome")
        txt = ", ".join(f"{SHORT.get(r.genome, can.at[r.genome, 'label'] if r.genome in can.index else r.genome)} "
                        f"{r.k_n} ({float(r.p):.2g})" for r in ex.itertuples() if r.k_n not in ("", "."))
        out.append(f"| {name} | {len(x)} | {len(ok)} | {('best: ' if not len(ok) and txt else '') + txt} |")
    assert len(used) == len(d), "a run fell in no group"
    a.out.write_text("\n".join(out) + "\n")
    (a.out.parent / "RABL_RUNS.md").write_text(f"{len(d)}\n")
    print("\n".join(out))
    print("runs", len(d))


if __name__ == "__main__":
    main()
