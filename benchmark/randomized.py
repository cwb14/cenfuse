#!/usr/bin/env python3
"""Randomized Hi-C controls of a release: did any test pass, and does any call differ from the same genome's run
without Hi-C?

  python -m benchmark.randomized JOBS ROOT OUT.tsv
Reads ROOT/<NEG|NEG2|TNEG>/<id>/cenfuse_summary.tsv and <g>.centromeres.tsv; the partner is the same genome in arm SG
if it was run, else S; for TNEG, arm T (the same inputs apart from Hi-C). Port of v12 bench/rel/randomized.py."""
import argparse
from pathlib import Path

import pandas as pd

CALL = ["start", "end", "confidence"]


def table(f, index=None):
    t = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False)
    return t.set_index(index) if index else t


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("jobs", help="job table (arm, id, ...)")
    ap.add_argument("root", type=Path, help="release root: ROOT/<arm>/<id>")
    ap.add_argument("out", help="output table")
    a = ap.parse_args(argv)
    done = {tuple(line.split("\t")[:2]) for line in open(a.jobs)}
    rows = []
    for arm, gid in sorted(done):
        if arm not in ("NEG", "NEG2", "TNEG"):
            continue
        d = a.root / arm / gid
        s = table(d / "cenfuse_summary.tsv")
        partner = "T" if arm == "TNEG" else next((x for x in ("SG", "S") if (x, gid) in done), None)
        for r in s.itertuples():
            g = r.genome
            c = table(d / f"{g}.centromeres.tsv", "chrom")
            changed = None
            if partner:
                b = table(a.root / partner / gid / f"{g}.centromeres.tsv", "chrom")
                changed = int(sum(any(c.at[ch, k] != b.at[ch, k] for k in CALL) for ch in c.index))
            rows.append(dict(arm=arm, id=gid, genome=g, gate=r.hic_gate_pass, support=r.hic_gate_support, fold=r.fold_gate,
                             align=r.align_gate, partner=partner or ".",
                             calls_changed=changed if changed is not None else ".",
                             hic_flags=int(c["flags"].str.contains("hic", regex=False).sum()),
                             hic_support=int(c["support"].str.contains("hic", regex=False).sum()),
                             hic_end=int(c["flags"].str.contains("hic_end", regex=False).sum()),
                             te_hic=getattr(r, "te_hic", "."), cis_line=getattr(r, "cis_line", ".")))
    d = pd.DataFrame(rows)
    d.to_csv(a.out, sep="\t", index=False)
    num = lambda c: pd.to_numeric(d[c], errors="coerce").fillna(0)      # noqa: E731
    print(f"{len(d)} runs: gate pass {int((num('gate') > 0).sum())}, support {int((num('support') > 0).sum())}, "
          f"fold passes {int((d.fold == 'pass').sum())}, align passes {int((d['align'] == 'pass').sum())}; "
          f"runs with a call differing from the partner {int((num('calls_changed') > 0).sum())}; "
          f"with a hic flag {int((d.hic_flags > 0).sum())}; with hic support {int((d.hic_support > 0).sum())}; "
          f"hic_end {int((d.hic_end > 0).sum())}; "
          f"te_hic passes or supports {int(d.te_hic.astype(str).str.contains('pass|support').sum())}; "
          f"cis passes or supports {int(d.cis_line.astype(str).str.contains('pass|support').sum())}")


if __name__ == "__main__":
    main()
