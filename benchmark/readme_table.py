#!/usr/bin/env python3
"""The README's results table, from a release's runs table (every arm).

  python -m benchmark.readme_table OUT.md [PREFIX]
PREFIX: score.py output prefix [paths.SCORES]. Cells are computed, but the row texts are typed in, a few of them
results ("a kinetochore on 40", "all 20 at the right end", "unchanged"): check them at each release.
Port of v12 bench/readme_table.py."""
import argparse
from pathlib import Path

import pandas as pd

from benchmark import paths

MAM = ["chm13", "mPanTro3", "mPanPan1", "mGorGor1", "mPonAbe1", "mPonPyg2", "MFA8", "Sscrofa", "horse1"]
BRAPA = ["CC1", "CC2", "CC3", "CCA03", "PC1", "TC1", "TU1", "OI1", "CX1", "MIZ1", "WTC1"]
SMALL = ["scer", "pfal", "spom", "ncra"]
FUNGI = ["vdah", "ztri", "cneoA", "cneoG1", "msym", "ctro"]
HOLO = ["celegans", "cnig", "bmori", "rbrev"]


def table(r):
    def tot(arms, genomes):
        x = r[r.arm.isin(arms) & r.genome.isin(genomes)]
        return int(x.right.sum()), int(x.chromosomes.sum()), int(x.conf_wrong.sum()), int(len(x))

    def add(a, b):
        return tuple(x + y for x, y in zip(a, b))

    def cell(t, low=False):
        right, n, cw, runs = t
        return f"{right}/{n}, {cw}" + (" (all `low`)" if low else "") if runs else ""

    rows = [
        ("9 T2T mammals (human, great apes, crab-eating macaque, pig, horse)", "satellites", cell(tot(["S"], MAM)),
         "+ genes + Hi-C: " + cell(tot(["FULL"], MAM))),
        ("HG002, diploid human T2T", "alpha satellite", cell(tot(["S"], ["hg002"])),
         "+ CpG: a kinetochore on 40, all in the active array"),
        ("rice", "CentO satellite", cell(tot(["S"], ["rice"])), "+ Hi-C, `tes/`, genes: " + cell(tot(["TGH"], ["rice"]))),
        ("*Brassica rapa*, 11 accessions", "satellite", cell(tot(["S"], BRAPA)),
         f"+ `tes/` (4 accessions): {cell(tot(['T'], BRAPA))}; {tot(['HT'], BRAPA)[3]} more accessions: "
         f"{cell(tot(['HT'], BRAPA))}"),
        ("soybean", "satellite", cell(tot(["S"], ["gmX26"])),
         "+ Hi-C: " + cell(tot(["SH"], ["gmX26"])) + "; + `tes/`: " + cell(tot(["T"], ["gmX26"]))),
        ("sorghum, *Brachypodium*, *Arabidopsis*", "satellite", cell(tot(["S"], ["sbBTx623", "brachy", "colcc"])),
         "+ `tes/`, Hi-C, methylation: " + cell(add(tot(["THn"], ["sbBTx623", "brachy"]), tot(["M"], ["colcc"])))),
        ("grape", "satellite", cell(tot(["HS"], ["vvPN"])),
         "; ".join(f"+ {w}: " + cell(tot([a], ["vvPN"])) for a, w in (("HT", "`tes/`"), ("HSH", "Hi-C"), ("HTH", "both"))
                   if tot([a], ["vvPN"])[3])),
        ("maize (Mo17), with knobs", "CentC satellite", cell(tot(["S"], ["mo17"])),
         "+ `tes/`, genes: " + cell(tot(["TG"], ["mo17"]))),
        ("chicken", "chromosome-specific satellites", cell(tot(["S"], ["ggswu"]), True),
         "+ Hi-C: " + cell(tot(["SH"], ["ggswu"]), True)),
        ("mouse, all acrocentric", "minor satellite beside the major satellite",
         cell(tot(["S"], ["mmT2T"])) + " (`low`; since 0.15 the AT trough chooses the minor satellite array on most chromosomes)", "+ Hi-C: unchanged"),
        ("barley", "CRM retrotransposons", cell(tot(["S"], ["barley"])),
         "+ `tes/`: " + cell(tot(["T"], ["barley"])) + "; + Hi-C: " + cell(tot(["SH"], ["barley"])) + "; both: "
         + cell(tot(["TH"], ["barley"]))),
        ("potato", "satellite on half, none on half", cell(tot(["S"], ["potato"])),
         "+ `tes/`: " + cell(tot(["T"], ["potato"])) + "; + Hi-C too: " + cell(tot(["TH"], ["potato"]))),
        ("*Aegilops tauschii*", "CRM beside satellites", cell(tot(["S"], ["aetAL878"])),
         "+ `tes/`: " + cell(tot(["T"], ["aetAL878"]))),
        ("moss (*Physcomitrium*)", "Copia (Bryco) retrotransposon clusters", cell(tot(["S"], ["ppat"]), True),
         "+ `tes/`: " + cell(tot(["T"], ["ppat"])) + "; + Hi-C: " + cell(tot(["SH"], ["ppat"])) + "; both: "
         + cell(tot(["TH"], ["ppat"]))),
        ("*Chlorella sorokiniana*", "2-7 Copia (Osser) insertions", cell(tot(["S"], ["csNS42"])),
         "+ `tes/`: " + cell(tot(["T"], ["csNS42"]))),
        ("*Chlamydomonas*", "LINE clusters", cell(tot(["S"], ["chlamy"]), True), "+ `tes/`: unchanged"),
        ("tomato", "no repeat, beside a dominant satellite", cell(tot(["S"], ["tomato"])),
         "+ Hi-C: " + cell(tot(["SH"], ["tomato"])) + "; + `tes/`: " + cell(tot(["T"], ["tomato"])) + "; both: "
         + cell(tot(["TH"], ["tomato"]))),
        ("horse EquCab3 (older assembly)", "satellites; chr11 none", cell(tot(["S"], ["horse"]), True), "+ Hi-C: unchanged"),
        ("budding and fission yeast, *Plasmodium*, *Neurospora*", "point, regional, AT-rich", cell(tot(["S"], SMALL)),
         "+ Hi-C: " + cell(tot(["SH"], SMALL))),
        ("5 more fungi (6 Hi-C libraries)", "regional", cell(tot(["S"], FUNGI)),
         "+ Hi-C: " + cell(tot(["SH"], FUNGI))),
    ]
    L = ["| genomes | centromeres made of | FASTA only: right, confidently wrong | more inputs: right, confidently wrong |",
         "|---|---|---|---|"]
    L += [f"| {a} | {b} | {c} | {d} |" for a, b, c, d in rows]
    x = r[r.arm.isin(["S"]) & r.genome.isin(HOLO)]
    y = r[r.arm.isin(["SH"]) & r.genome.isin(HOLO)]
    L.append(f"| 4 holocentric genomes (any confident call is wrong) | no single centromere | 0 confident in "
             f"{int(x.chromosomes.sum())}" if int(x.conf_wrong.sum()) == 0 else
             f"| 4 holocentric genomes | no single centromere | {int(x.conf_wrong.sum())} confident in "
             f"{int(x.chromosomes.sum())}")
    L[-1] += f" | + Hi-C: {int(y.conf_wrong.sum())} confident |"
    return L


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path, help="output markdown file")
    ap.add_argument("prefix", nargs="?", default=paths.SCORES, help="score.py output prefix [paths.SCORES]")
    a = ap.parse_args(argv)
    L = table(pd.read_csv(f"{a.prefix}.runs.tsv", sep="\t", dtype={"hic_gate_pass": str}))
    a.out.write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
