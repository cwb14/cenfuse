#!/usr/bin/env python3
"""The frozen benchmark scorer: each chromosome of a set of cenfuse runs, scored by the rule its benchmark used.

  python -m benchmark.score JOBS ROOT|ref OUT_PREFIX [--manifest truth.tsv] [-v]

A port of v12 bench/rel/score_rel.py (the scorer of 0.13's docs, identical to v11's) whose truth table is read from
benchmark/truth.tsv, with truth paths relative to the data root (benchmark/paths.py). Same rules, columns, row order
and numbers. JOBS: one arm per line (arm, id, inputs, threads, satellite dir, reference run). ROOT: a release root
holding ROOT/<arm>/<id>; 'ref' scores each job's reference run (column 6) instead. Writes OUT_PREFIX.chrom.tsv (one
row per chromosome with truth, or of a holocentric control) and OUT_PREFIX.runs.tsv (per arm and genome: right
high/medium/low, confidently wrong, wrong low, no call, chance, P, Hi-C summary), and prints the runs table.

Rules (manifest columns rule, tol_bp, hic_bins; tol = the distance allowed between call and truth):
  mammal   >= 1 bp overlap with tier FW: the truth files in order (W, then F over it), one interval per chromosome,
           the last line wins (v3 bench/score.py)
  overlap  >= 1 bp overlap with a truth interval (satellome agreement; Arabidopsis)
  active   >= 1 bp overlap with a CenSat row named active_hor* (HG002)
  gap      gap <= tol, tol = hic_bins x the run's Hi-C bin when both exist, else tol_bp (gap 0: overlap or touching)
  bin1     the same with 1 bin: 1 Hi-C bin of the run, 250 kb without Hi-C (small genomes)
  holo     no truth: every call is wrong; a confident one is confidently wrong
Outcome: right_<confidence>, conf_wrong, wrong_low or none (no call); confident = high or medium. Chance of a row =
(call + truth + 2 tol) / chromosome length, at most 1 (v6 score.py); P = Poisson-binomial tail of the right count.
A genome whose truth file is missing is skipped with a warning."""
import argparse
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from benchmark import paths
from cenfuse.util import poisson_binomial_sf

RULES = ("mammal", "overlap", "active", "gap", "bin1", "holo")
STRICT = ("mammal", "overlap", "active")          # right = >= 1 bp overlap; the other rules: gap <= tol
HIC = ("hic_bin", "hic_gate_pass", "hic_gate_support", "hic_gate_untested", "fold_gate", "align_gate", "hic_axis",
       "hic_architecture", "cpg_gate_pass", "hic_sequence_conflict")
RUNS_SHOWN = ["arm", "genome", "kind", "chromosomes", "right", "right_high", "right_medium", "right_low", "conf_wrong",
              "wrong_low", "none", "chance", "hic_gate_pass", "fold_gate", "align_gate"]


@dataclass(frozen=True)
class Truth:
    genome: str
    kind: str
    rule: str
    tol: int                 # bp allowed between call and truth (without Hi-C, when hic_bins is set)
    hic_bins: int | None     # Hi-C bins allowed when the run has Hi-C
    files: tuple             # truth BEDs (none for a holocentric control)
    source: str = "."
    notes: str = "."

    @property
    def missing(self):
        return [f for f in self.files if not f.exists()]


def load_manifest(path=None, data=None):
    """genome -> Truth, from the manifest; truth paths are relative to data (an absolute path is kept)."""
    path, data = Path(path or paths.MANIFEST), Path(data or paths.DATA)
    m = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    lack = [c for c in ("genome", "kind", "rule", "tol_bp", "hic_bins", "truth") if c not in m.columns]
    if lack:
        sys.exit(f"score: {path}: no column {', '.join(lack)}")
    T = {}
    for r in m.to_dict("records"):
        g = r["genome"]
        files = () if r["truth"] in (".", "") else tuple(data / f for f in r["truth"].split(";"))
        if g in T or r["rule"] not in RULES or (r["rule"] == "holo") != (not files):
            sys.exit(f"score: {path}: {g}: duplicate genome, unknown rule {r['rule']!r}, or truth files that do not "
                     "fit the rule (holo has none, every other rule needs one)")
        T[g] = Truth(g, r["kind"], r["rule"], int(r["tol_bp"]),
                     int(r["hic_bins"]) if r["hic_bins"] not in (".", "") else None, files,
                     r.get("source", "."), r.get("notes", "."))
    return T


def mammal_truth(files):
    """Tier FW (v3 bench/score.py load_truth): the files in order, one interval per chromosome, the last line wins."""
    d = {}
    for f in files:
        for line in open(f):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3:
                d[p[0]] = [(int(p[1]), int(p[2]))]
    return d


def read_truth(path, active=False, tr=None):
    """chromosome -> [(start, end)] of a BED; active: only rows named active_hor* (CenSat)."""
    tr = {} if tr is None else tr
    for line in open(path):
        p = line.split()
        if len(p) < 3 or line.startswith(("#", "track")):
            continue
        if active and not (len(p) > 3 and p[3].startswith("active_hor")):
            continue
        tr.setdefault(p[0], []).append((int(p[1]), int(p[2])))
    return tr


def load_truth(t):
    """chromosome -> [(start, end)] (0-based, half-open), or None for a holocentric control."""
    if t.rule == "holo":
        return None
    if t.rule == "mammal":
        return mammal_truth(t.files)
    tr = {}
    for f in t.files:
        read_truth(f, active=t.rule == "active", tr=tr)
    return tr


def tolerance(t, hic_bin):
    """bp allowed between call and truth: hic_bins x the run's Hi-C bin when both exist, else tol_bp."""
    return t.hic_bins * hic_bin if t.hic_bins and hic_bin else t.tol


def run_summary(run, gid):
    """The run's cenfuse_summary.tsv row of genome gid (the first row if no row names it)."""
    s = pd.read_csv(Path(run) / "cenfuse_summary.tsv", sep="\t", dtype=str, keep_default_na=False)
    return s[s.genome == gid].iloc[0] if "genome" in s and (s.genome == gid).any() else s.iloc[0]


def hic_bin_of(s):
    hb = s.get("hic_bin", ".")
    return int(float(hb)) if hb not in (".", "", None) else None


def score_genome(arm, gid, run, cen, t):
    """One row per chromosome of cen (a <gid>.centromeres.tsv) with truth, and the run's Hi-C summary fields."""
    c = pd.read_csv(cen, sep="\t", dtype=str, keep_default_na=False)
    s = run_summary(run, gid)
    hic_bin = hic_bin_of(s)
    truth = load_truth(t)
    tol = tolerance(t, hic_bin)
    strict = t.rule in STRICT
    rows = []
    for r in c.to_dict("records"):
        ch = r["chrom"]
        if truth is not None and ch not in truth:
            continue
        L = int(r["chrom_len"]) if r.get("chrom_len", ".") not in (".", "") else np.nan
        base = dict(arm=arm, genome=gid, chrom=ch, chrom_len=L, kind=t.kind, rule=t.rule, tol=tol, hic_bin=hic_bin or ".",
                    start=r["start"], end=r["end"], confidence=r["confidence"], relation=r["relation"],
                    flags=r["flags"], support=r["support"], against=r.get("against", "."),
                    seq_start=r.get("seq_start", "."), seq_end=r.get("seq_end", "."),
                    seq_confidence=r.get("seq_confidence", "."), cpg_kt_start=r.get("cpg_kt_start", "."),
                    cpg_kt_end=r.get("cpg_kt_end", "."), morphology=r.get("morphology", "."),
                    ci=r.get("ci", "."), telomeres=r.get("telomeres", "."), t2t=r.get("t2t", "."))
        conf = r["confidence"] in ("high", "medium")
        if r["start"] in (".", ""):
            rows.append(dict(base, ok="none", gap=np.nan, truth=".", chance=0.0, outcome="none"))
            continue
        a, z = int(r["start"]) - 1, int(r["end"])
        if truth is None:                                   # holocentric: no single centromere exists
            rows.append(dict(base, ok="W", gap=np.nan, truth=".", chance=0.0,
                             outcome="conf_wrong" if conf else "wrong_low"))
            continue
        best = None
        for ts, te in truth[ch]:                            # the nearest interval, a hit first
            gap = max(0, max(a, ts) - min(z, te))
            hit = (max(a, ts) < min(z, te)) if strict else gap <= tol
            key = (not hit, gap)
            if best is None or key < best[0]:
                best = (key, ts, te, gap, hit)
        _, ts, te, gap, hit = best
        chance = min(1.0, ((z - a) + (te - ts) + 2 * tol) / L) if L == L else np.nan
        ok = "R" if hit else "W"
        outcome = ("right_" + r["confidence"]) if hit else ("conf_wrong" if conf else "wrong_low")
        rows.append(dict(base, ok=ok, gap=gap, truth=f"{ts}-{te}", chance=chance, outcome=outcome))
    return pd.DataFrame(rows), {k: s.get(k, ".") for k in HIC}


def summarize(d, arm, gid, kind, hic):
    n = len(d)
    k = int((d.ok == "R").sum())
    out = dict(arm=arm, genome=gid, kind=kind, chromosomes=n, right=k,
               right_high=int((d.outcome == "right_high").sum()), right_medium=int((d.outcome == "right_medium").sum()),
               right_low=int((d.outcome == "right_low").sum()), conf_wrong=int((d.outcome == "conf_wrong").sum()),
               wrong_low=int((d.outcome == "wrong_low").sum()), none=int((d.outcome == "none").sum()),
               chance=round(float(d.chance.sum()), 2) if n else 0.0,
               P=float(poisson_binomial_sf(k, list(d.chance))) if n and kind != "holocentric" else np.nan)
    out.update(hic)
    return out


def done_ok(run, arm, gid):
    """True if the run finished with exit 0; reports a pending or failed run on stderr."""
    if not (run / "DONE").exists():
        print(f"pending: {arm}/{gid} ({run})", file=sys.stderr)
        return False
    done = (run / "DONE").read_text().split()
    if done and done[0].isdigit() and done[0] != "0":
        print(f"failed: {arm}/{gid} exit {done[0]}", file=sys.stderr)
        return False
    return True


def jobs_of(path):
    """The job table's rows (tab-separated fields); blank and '#' lines skipped."""
    return [p for p in (line.rstrip("\n").split("\t") for line in open(path)) if p[0] and not p[0].startswith("#")]


def usable(t, warned):
    """False (with one warning per genome) if a truth file of t is missing."""
    if t.missing:
        if t.genome not in warned:
            print(f"warning: truth of {t.genome} missing, its runs skipped: {', '.join(map(str, t.missing))}",
                  file=sys.stderr)
            warned.add(t.genome)
        return False
    return True


def score_jobs(jobs, root, T, verbose=False):
    """Per-chromosome tables and per-run summaries of every finished run of the job table."""
    chrom, runs, warned = [], [], set()
    for p in jobs_of(jobs):
        arm, gid = p[0], p[1]
        run = paths.local(p[5]) if root == "ref" else Path(root) / arm / gid
        if not done_ok(run, arm, gid):
            continue
        for cen in sorted(run.glob("*.centromeres.tsv")):
            g = cen.name[: -len(".centromeres.tsv")]
            if g not in T or not usable(T[g], warned):
                continue
            d, hic = score_genome(arm, g, run, cen, T[g])
            if d.empty:
                print(f"no chromosome with truth: {arm}/{g} ({run})", file=sys.stderr)
                continue
            chrom.append(d)
            runs.append(summarize(d, arm, g, T[g].kind, hic))
            if verbose:
                print(f"scored {arm}/{g}: {len(d)} chromosomes, {runs[-1]['right']} right", file=sys.stderr)
    return chrom, runs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("jobs", help="job table (tab-separated: arm, id, inputs, threads, satellite dir, reference run)")
    ap.add_argument("root", help="release root holding ROOT/<arm>/<id>, or 'ref': each job's reference run")
    ap.add_argument("prefix", help="output prefix: PREFIX.chrom.tsv, PREFIX.runs.tsv")
    ap.add_argument("--manifest", type=Path, default=paths.MANIFEST, help="truth manifest [benchmark/truth.tsv]")
    ap.add_argument("-v", "--verbose", action="store_true", help="log every scored run on stderr")
    a = ap.parse_args(argv)
    chrom, runs = score_jobs(a.jobs, a.root, load_manifest(a.manifest), a.verbose)
    if not chrom:
        sys.exit("score: no finished run with truth")
    pd.concat(chrom, ignore_index=True).to_csv(f"{a.prefix}.chrom.tsv.tmp", sep="\t", index=False)
    r = pd.DataFrame(runs)
    r.to_csv(f"{a.prefix}.runs.tsv.tmp", sep="\t", index=False)
    os.replace(f"{a.prefix}.chrom.tsv.tmp", f"{a.prefix}.chrom.tsv")     # atomic: readers never see a partial table
    os.replace(f"{a.prefix}.runs.tsv.tmp", f"{a.prefix}.runs.tsv")
    with pd.option_context("display.width", 250, "display.max_rows", 400, "display.max_columns", 30):
        print(r[RUNS_SHOWN].to_string(index=False))


if __name__ == "__main__":
    main()
