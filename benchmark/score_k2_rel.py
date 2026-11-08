#!/usr/bin/env python3
"""The morphology check K2 (v7 bench/score_k.py, frozen) on a release's runs: acrocentric vs bi-armed from the
integrated call's centromeric index (acrocentric below 0.125) against the literature karyotype, on right calls
(within 1 Mb of truth).

  python -m benchmark.score_k2_rel OUT.tsv [--runs ROOT] [--manifest truth.tsv]
Reads ROOT/FULL/<g> for goat, sheep, horse T2T and pig and ROOT/S/mmT2T for mouse (the arms v7 scored), the truth
files of the manifest, and NCBI assembly reports (paths.ASM_REPORTS) or the mouse name map for chromosome names.
Port of v12 bench/rel/score_k2_rel.py with v7's k2() inlined."""
import argparse
from pathlib import Path

import pandas as pd

from benchmark import paths
from benchmark import score as S

KARYO = {   # species: (arm, genome, GenBank accession, acrocentric chromosomes (Y excluded)); mouse: all acrocentric
    "goat": ("FULL", "goat1", "GCA_040806595.1", {str(i) for i in range(1, 30)} | {"X"}),
    "sheep": ("FULL", "sheep1", "GCA_040805955.1", {str(i) for i in range(4, 27)} | {"X"}),
    "horse": ("FULL", "horse1", "GCA_052818215.1", {str(i) for i in range(14, 32)}),
    "pig": ("FULL", "Sscrofa", "GCA_000003025.7", {str(i) for i in range(13, 19)}),
    "mouse": ("S", "mmT2T", None, None),
}


def report_map(gca, reports=None):
    """accession -> chromosome name (assembled molecules) from an NCBI assembly report."""
    m = {}
    for line in open(Path(reports or paths.ASM_REPORTS) / f"{gca}_assembly_report.txt"):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) > 4 and p[1] == "assembled-molecule":
            m[p[4]] = p[2]
    return m


def k2(rel, T, mouse_names=None):
    rows = []
    names = dict(x.split()[::-1][:2] for x in open(mouse_names or paths.MOUSE_NAMES))
    for sp, (arm, g, gca, acro) in KARYO.items():
        t = pd.read_csv(rel / arm / g / f"{g}.centromeres.tsv", sep="\t", dtype=str, keep_default_na=False)
        if gca:
            m = report_map(gca)
        else:
            m = {acc: c.replace("chr", "") for acc, c in names.items()}
            acro = set(m.values()) - {"Y"}
        tr = {}
        for f in T[g].files:                      # v7's reader: every line counts (no tier FW for horse and pig)
            S.read_truth(f, tr=tr)
        for r in t.itertuples():
            c = m.get(r.chrom)
            if c is None or c in ("Y", "MT") or r.chrom not in tr:
                continue
            if r.start == ".":
                rows.append(dict(genome=sp, chrom=r.chrom, chr=c, right=False, ci=".", expected="acro" if c in acro
                                 else "biarmed", observed="none", agree=False))
                continue
            s, e = int(r.start) - 1, int(r.end)
            gap = min(max(0, max(s, a) - min(e, b)) for a, b in tr[r.chrom])
            obs = "acro" if float(r.ci) < 0.125 else "biarmed"
            exp = "acro" if c in acro else "biarmed"
            rows.append(dict(genome=sp, chrom=r.chrom, chr=c, right=gap <= 1_000_000, ci=r.ci, morphology=r.morphology,
                             expected=exp, observed=obs, agree=obs == exp, confidence=r.confidence))
    return pd.DataFrame(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", help="output table")
    ap.add_argument("--runs", type=Path, default=paths.RUNS, help="release root [paths.RUNS]")
    ap.add_argument("--manifest", type=Path, default=paths.MANIFEST, help="truth manifest [benchmark/truth.tsv]")
    a = ap.parse_args(argv)
    T = S.load_manifest(a.manifest)
    lost = [g for _, g, _, _ in KARYO.values() if T[g].missing]
    if lost:
        raise SystemExit(f"score_k2_rel: truth missing for {', '.join(lost)}")
    d = k2(a.runs, T)
    d.to_csv(a.out, sep="\t", index=False)
    r = d[d.right]
    print(f"right calls: {len(r)}; agree {int(r.agree.sum())}/{len(r)}")
    print(r.groupby("genome").agree.agg(["sum", "count"]).to_string())
    print(r[~r.agree][["genome", "chrom", "chr", "ci", "expected", "observed", "confidence"]].to_string(index=False))


if __name__ == "__main__":
    main()
