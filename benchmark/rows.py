#!/usr/bin/env python3
"""The docs' result rows from a score.py per-chromosome table, one line per row; with a second prefix, both side by
side ('**' marks a row that differs). Columns: chromosomes | right (high/medium/low) | confidently wrong | wrong low |
no call.

  python -m benchmark.rows PREFIX [PREFIX2]
Port of v12 bench/rel/rows.py."""
import argparse

import pandas as pd

MAM9 = ["chm13", "mPanTro3", "mPanPan1", "mGorGor1", "mPonAbe1", "mPonPyg2", "MFA8", "Sscrofa", "horse1"]
AGREE = ["mSymSyn1", "MMU8", "sheep1", "goat1", "mGirTip1"]
BRAPA = ["CC1", "CC2", "CC3", "CCA03", "PC1", "TC1", "TU1", "OI1", "CX1", "MIZ1", "WTC1"]
BRAPA_T = ["CC1", "CC2", "CC3", "CCA03"]
SMALL = ["scer", "pfal", "spom", "ncra", "vdah", "ztri", "cneoA", "cneoG1", "msym", "ctro"]
HOLO = ["celegans", "bmori", "rbrev", "cnig"]
ROWS = [("9 T2T mammals", MAM9, ["S", "FULL", "CPG", "NEG", "NEG2"]),
        ("5 mammals (satellome agreement)", AGREE, ["S", "FULL", "NEG", "NEG2"]),
        *[(g, [g], None) for g in ["rice", "barley", "potato", "horse", "colcc", "ggswu", "mo17", "mmT2T", "ppat",
                                   "chlamy", "brachy"]],
        ("B. rapa 11", BRAPA, ["S"]), ("B. rapa 4", BRAPA_T, ["T", "S", "SH"]),
        *[(g, [g], None) for g in ["gmX26", "sbBTx623", "aetAL878", "csNS42", "tomato", "hg002"] + SMALL + HOLO]]


def tab(prefix):
    """(row label, arm) -> the row's text; a group row only where every genome of the group was scored."""
    d = pd.read_csv(f"{prefix}.chrom.tsv", sep="\t", dtype={"start": str})
    out = {}
    for label, gs, arms in ROWS:
        x = d[d.genome.isin(gs)]
        for arm in (arms or sorted(x.arm.unique())):
            y = x[x.arm == arm]
            if y.empty or (len(gs) > 1 and y.genome.nunique() < len(gs)):
                continue
            o = y.outcome.value_counts()
            r = int((y.ok == "R").sum())
            out[(label, arm)] = (f"{len(y)} | {r} ({o.get('right_high', 0)}/{o.get('right_medium', 0)}/"
                                 f"{o.get('right_low', 0)}) | cw {o.get('conf_wrong', 0)} | wl {o.get('wrong_low', 0)}"
                                 f" | none {o.get('none', 0)}")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefix", help="score.py output prefix (PREFIX.chrom.tsv)")
    ap.add_argument("prefix2", nargs="?", help="a second prefix, shown beside the first (an earlier release)")
    x = ap.parse_args(argv)
    a = tab(x.prefix)
    b = tab(x.prefix2) if x.prefix2 else {}
    for k in list(a) + [k for k in b if k not in a]:
        va, vb = a.get(k, "-"), b.get(k, "-")
        mark = "" if not b else ("   " if va == vb else " **")
        print(f"{mark}{k[0]:34s} {k[1]:5s} {va:48s} {vb if b else ''}")


if __name__ == "__main__":
    main()
