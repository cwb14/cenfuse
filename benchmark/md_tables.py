#!/usr/bin/env python3
"""Markdown tables and numbers for docs/benchmark.md (and the gallery's rows) from a release's scores.

  python -m benchmark.md_tables PREFIX OUTDIR [--runs ROOT] [--canonical FILE]
PREFIX: score.py output (PREFIX.chrom.tsv, PREFIX.runs.tsv); ROOT: the release root the scores came from (for caveats,
and the FASTA-only runs' summaries: telomeres and rDNA). Writes OUTDIR/<NAME>.md for each placeholder NAME of
benchmark.md computed here, and OUTDIR/gallery_rows.tsv. Every number is computed from the tables; nothing is typed in.
Port of v12 bench/rel/md_tables.py."""
import argparse
import math
import re
from pathlib import Path

import pandas as pd

from benchmark import paths

ITAL = {"Brachypodium distachyon", "Aegilops tauschii", "Arabidopsis thaliana", "Physcomitrium patens",
        "Chlamydomonas reinhardtii", "Chlorella sorokiniana", "Saccharomyces cerevisiae", "Schizosaccharomyces pombe",
        "Neurospora crassa", "Verticillium dahliae", "Zymoseptoria tritici", "Cryptococcus neoformans",
        "Malassezia sympodialis", "Candida tropicalis", "Plasmodium falciparum", "Caenorhabditis elegans",
        "Caenorhabditis nigoni", "Bombyx mori", "Rhynchospora breviuscula", "Brassica rapa"}
TRUTH = {"CENP-A": "CENP-A", "agreement": "none (satellome locus)", "CENH3": "CENH3", "CENH3 (published)":
         "CENH3 (published domains)", "CENP-A windows": "CENP-A windows", "minor satellite": "minor-satellite array",
         "CEN": "CEN / CENP-A / CenH3", "holocentric": "none (holocentric)", "CenSat active": "CenSat active array"}
ARMTXT = {"FULL": "genes, Hi-C", "TGH": "`tes/`, genes, Hi-C", "TH": "`tes/`, Hi-C", "THn": "`tes/`, Hi-C",
          "TG": "`tes/`, genes", "M": "`tes/`, ONT methylation", "ONT": "ONT CpG methylation", "CPG": "genes, Hi-C, CpG"}


def canonical(path=None):
    """genome -> label, lineage, centromere_type, arm_fasta, arm_all (the arm each table shows as 'all inputs')."""
    return pd.read_csv(path or paths.CANONICAL, sep="\t").set_index("genome")


def label(g, can, other=None):
    """The genome's label, its species in italics."""
    s = can.at[g, "label"] if g in can.index else (other or {}).get(g, g)
    for sp in sorted(ITAL, key=len, reverse=True):
        if sp in s:
            return s.replace(sp, f"*{sp}*")
    if s.startswith("moss (Physcomitrium"):
        return "moss (*Physcomitrium patens*)"
    return s


def cell(r):
    """right/chromosomes (confidently wrong); for a holocentric control: confident calls / chromosomes."""
    if r is None:
        return ""
    if r["kind"] == "holocentric":
        return f"{r['conf_wrong']} confident / {r['chromosomes']}"
    s = f"{r['right']}/{r['chromosomes']}"
    if r["conf_wrong"]:
        s += f" ({r['conf_wrong']})"
    return s


def caveat(rel, g, arm, ch):
    f = rel / arm / g / f"{g}.centromeres.tsv"
    if not f.exists():
        return "."
    t = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False).set_index("chrom")
    c = t.at[ch, "caveat"] if ch in t.index else "."
    return c.replace("|", "/")


def pct(k, n):
    return f"{100 * k / n:.1f}%" if n else ""


def calibration(tr):
    """How often each confidence level is right: one sentence and a table, over the chromosomes with truth."""
    n, right = len(tr), int((tr.ok == "R").sum())
    lines = [f"{tr.genome.nunique()} genomes, each with all its inputs: {n} chromosomes with truth, {right} right "
             f"({pct(right, n)}), {int((tr.outcome == 'conf_wrong').sum())} confidently wrong.", "",
             "| confidence | calls | right | wrong | share right |", "|---|---|---|---|---|"]
    called = 0
    for lev in ("high", "medium", "low"):
        y = tr[(tr.confidence == lev) & (tr.ok != "none")]
        k = int((y.ok == "R").sum())
        called += len(y)
        lines.append(f"| `{lev}` | {len(y)} | {k} | {len(y) - k} | {pct(k, len(y))} |")
    lines.append(f"| no call | {n - called} | | | |")
    return "\n".join(lines)


def telomeres(rel, order, lab):
    """Telomere and 45S recall per genome, from the summaries of the FASTA-only runs (ROOT/S/<g>): the unit found,
    chromosome ends with a terminal array (of 2 per chromosome), complete chromosomes (an array at both ends, no gap),
    interstitial arrays, 45S arrays and the chromosomes holding them. One row per assembly (cneoA, cneoG1 share one);
    the assemblies with an array at every end as a sentence. Returns (table, sentence)."""
    lines = ["| genome | chromosomes | telomere unit | ends with an array | complete chromosomes | interstitial arrays "
             "| 45S arrays (chromosomes) |", "|---|---|---|---|---|---|---|"]
    seen, full = set(), []
    for g in order:
        f = rel / "S" / g / "cenfuse_summary.tsv"
        if not f.exists():
            continue
        name = re.sub(r" \((anaphase|G1) Hi-C\)$", "", lab(g))
        if name in seen:
            continue
        seen.add(name)
        r = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False).iloc[0]
        n = int(r["chromosomes"])
        unit = r["telomere_unit"] if r["telomere_unit"] not in ("", ".") else "none found"
        nor = len([x for x in r["rdna_chromosomes"].split(",") if x not in (".", "")])
        if int(r["telomere_ends"]) == 2 * n:
            full.append(f"{name} ({r.get('t2t_chromosomes', '.')} complete)")
            continue
        lines.append(f"| {name} | {n} | {unit} | {r['telomere_ends']}/{2 * n} | {r.get('t2t_chromosomes', '.')} | "
                     f"{r['telomere_interstitial']} | {r['rdna_arrays']} ({nor}) |")
    return "\n".join(lines), (f"An array at every chromosome end in {len(full)} assemblies: " + ", ".join(full) + "."
                               if full else "No assembly has an array at every chromosome end.")


def write_tables(prefix, out, rel, can):
    out.mkdir(parents=True, exist_ok=True)
    runs = pd.read_csv(f"{prefix}.runs.tsv", sep="\t", dtype={"hic_gate_pass": str})
    chrom = pd.read_csv(f"{prefix}.chrom.tsv", sep="\t", dtype={"start": str, "end": str})
    R = {(r.genome, r.arm): r._asdict() for r in runs.itertuples(index=False)}

    def get(g, a):
        return R.get((g, a))

    def w(name, text):
        (out / f"{name}.md").write_text(text.rstrip("\n") + "\n")

    def lab(g):
        return label(g, can)

    # ---- every genome
    lines = ["| genome | truth | FASTA only | + genes | + `tes/` | + Hi-C | all inputs | chance (FASTA only) |",
             "|---|---|---|---|---|---|---|---|"]
    order = list(can.index)
    lineage_prev = None
    for g in order:
        if not any(k[0] == g for k in R):
            continue
        lin = can.at[g, "lineage"]
        if lin != lineage_prev:
            lines.append(f"| **{lin}** | | | | | | | |")
            lineage_prev = lin
        s, genes, tes, hic = get(g, "S"), get(g, "SG") or get(g, "G"), get(g, "T") or get(g, "HT"), get(g, "SH")
        allarm = can.at[g, "arm_all"]
        al = get(g, allarm)
        kind = (s or al)["kind"]
        alltxt = cell(al) + f" ({ARMTXT.get(allarm, allarm)})" if al and allarm not in ("S", "T", "SH", "SG", "G") else ""
        ch = f"{s['chance']:.1f}" if s and kind not in ("holocentric",) else ""
        lines.append(f"| {lab(g)} | {TRUTH.get(kind, kind)} | {cell(s)} | {cell(genes)} | {cell(tes)} | {cell(hic)} "
                     f"| {alltxt} | {ch} |")
    w("TABLE_EVERY_GENOME", "\n".join(lines))

    # ---- mammals
    mam = [g for g in order if str(can.at[g, "lineage"]).startswith("mammals") and get(g, "FULL")]
    ml = []
    for g in mam:
        s, f = get(g, "S"), get(g, "FULL")
        gate = f["hic_gate_pass"]
        try:
            sup = int(float(f.get("hic_gate_support") or 0))
        except ValueError:
            sup = 0
        extra = [f"support tier on {sup}"] if str(gate) in ("0", "0.0") and sup else []
        gtxt = f"{gate}" + (f" ({', '.join(extra)})" if extra else "")
        ml.append(f"| {lab(g)} | {TRUTH.get(f['kind'], f['kind'])} | {cell(s)} | {cell(f)} | {gtxt} |")
    w("TABLE_MAMMALS", "\n".join(ml))
    m9 = chrom[(chrom.arm == "FULL") & (chrom.kind == "CENP-A") & (chrom.genome.isin(mam))]
    hi = m9[m9.confidence == "high"]
    w("MAMMAL_HIGH", f"{int((hi.ok == 'R').sum())} of {len(hi)} times")

    # ---- small genomes
    sl = []
    for g in ["scer", "pfal", "spom", "ncra", "vdah", "ztri", "cneoA", "cneoG1", "msym", "ctro"]:
        s, h = get(g, "S"), get(g, "SH")
        if not h:
            continue
        b = h.get("hic_bin")
        bt = (f"{int(float(b)) // 1000} kb" if b not in (None, ".", "") and not (isinstance(b, float) and math.isnan(b))
              else ".")
        hh = f"{cell(h)} [{h['right_high']}/{h['right_medium']}/{h['right_low']}]" if h["right"] else cell(h)
        sl.append(f"| {lab(g)} | {can.at[g, 'centromere_type']} | {bt} | {cell(s)} | {hh} |")
    w("TABLE_SMALL", "| genome | centromeres | Hi-C bin | FASTA only | + Hi-C |\n|---|---|---|---|---|\n"
      + "\n".join(sl))

    # ---- randomized
    neg = runs[runs.arm.isin(["NEG", "NEG2"])]
    w("RANDOM_RUNS", f"{neg[['arm', 'genome']].drop_duplicates().shape[0]}")

    # ---- confident misses (canonical arms)
    rows = []
    for g in order:
        a = can.at[g, "arm_all"]
        x = chrom[(chrom.genome == g) & (chrom.arm == a) & (chrom.outcome == "conf_wrong")]
        for r in x.itertuples():
            gap = f"{r.gap / 1e6:.1f} Mb" if r.gap == r.gap else "."
            rows.append(f"| {lab(g)} ({ARMTXT.get(a, a)}) | {r.chrom} | `{r.confidence}` | {gap} | "
                        f"{caveat(rel, g, a, r.chrom)} |")
    w("TABLE_CONF_MISSES", "| genome (inputs besides the FASTA) | chromosome | confidence | distance from truth | "
      "caveat (what would settle it) |\n"
      "|---|---|---|---|---|\n" + "\n".join(rows))

    # ---- calibration over canonical arms
    cn = chrom[chrom.apply(lambda r: r.genome in can.index and can.at[r.genome, "arm_all"] == r.arm, axis=1)]
    tr = cn[~cn.kind.isin(["agreement", "holocentric"])]
    w("CALIBRATION", calibration(tr))

    # ---- support tier
    sup = chrom[chrom["flags"].fillna("").astype(str).str.contains("hic_support")]
    sup = sup.drop_duplicates(["genome", "chrom", "arm"])
    w("SUPPORT_TIER", f"{len(sup)} calls flagged hic_support in {sup[['genome', 'arm']].drop_duplicates().shape[0]} runs; "
      f"outcomes: {sup.outcome.value_counts().to_dict()}; by run: "
      + ", ".join(f"{a}/{g} {n}" for (g, a), n in sup.groupby(['genome', 'arm']).size().items()))

    # ---- telomeres and 45S rDNA: the FASTA-only run of every genome (both read only the FASTA)
    tab, full = telomeres(rel, order, lab)
    w("TABLE_TELOMERES", tab)
    w("TELOMERES_FULL", full)

    # ---- gallery rows
    gal = []
    for g in order:
        a = can.at[g, "arm_all"]
        r = get(g, a)
        if r:
            gal.append(f"{g}\t{a}\t{cell(r)}\t{r['right_high']}/{r['right_medium']}/{r['right_low']}\t{r['kind']}")
    (out / "gallery_rows.tsv").write_text("\n".join(gal) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefix", help="score.py output prefix")
    ap.add_argument("out", type=Path, help="output folder")
    ap.add_argument("--runs", type=Path, default=paths.RUNS, help="release root of the scored runs [paths.RUNS]")
    ap.add_argument("--canonical", type=Path, default=paths.CANONICAL, help="genome table [benchmark/canonical.tsv]")
    a = ap.parse_args(argv)
    write_tables(a.prefix, a.out, a.runs, canonical(a.canonical))


if __name__ == "__main__":
    main()
