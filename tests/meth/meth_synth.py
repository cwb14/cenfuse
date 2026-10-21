"""Tiny synthetic genome for the tests: six 60-Mb chromosomes plus two short unplaced scaffolds, one CHG site per
kb, three samples (S1 also has CG), LTRquest-style GFF3, satellite-engine-style arrays GFF3, genes GFF3, anchors and
exclusions. Planted:
- chr1-chr5: a centromere (satellite arrays + a CRM-rich pericentromere) whose CHG dip is carried by CRM cytosines
  inside the domain (CRM 0.45 -> 0.20); the satellite stays CHG-high (0.68).
- chr6: true negative (CRM-rich region, no dip).
- chr2 48-50 Mb: gene-dense decoy (genes 70% of sequence at CHG 0.10): plain CHG dips, the adjusted track does not.
- S1 chr4 10.0-10.1 Mb: duplicated segment (3x depth, CHG 0.15): removed by the depth filter.
- S3 chr3: its centromere moved from 35 to 38 Mb (plant-specific).
- Some CRMs nest inside Tekay elements (depth 1): their cytosines must be classed as CRM.
Truth is returned as a dict."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pysam

L = 60_000_000
CHROMS = [f"chr{i}" for i in range(1, 7)]
CEN = {"chr1": 20e6, "chr2": 28e6, "chr3": 35e6, "chr4": 41e6, "chr5": 25e6}
CRM_CLUSTER = {**CEN, "chr6": 30e6}
SHIFT = ("S3", "chr3", 38e6)
DECOY = ("chr2", 48e6, 50e6)
DUP = ("S1", "chr4", 10_000_000, 10_100_000)
DOMAIN = 650_000          # dip width (bp), ~ D for this genome
P = {"sat": 0.68, "CRM": 0.45, "Tekay": 0.50, "gene": 0.10, "other": 0.45}
P_CRM_DOMAIN = 0.20
SAMPLES = ("S1", "S2", "S3")


def annotation(rng):
    """Interval lists per class (0-based half-open). Depth-0 elements sit in separate 10-kb slots (no overlaps), so
    the only overlaps are the deliberate ones: nested CRMs (depth 1) and satellite arrays over everything."""
    rows = []
    slot = np.arange(0, L, 10_000)
    for c in CHROMS:
        u = rng.random(len(slot))
        near = np.abs(slot - CRM_CLUSTER[c]) <= 6e6
        crm = u < np.where(near, 0.35, 0.01)
        tek = ~crm & (u < np.where(near, 0.35, 0.01) + 0.30)
        gene = ~crm & ~tek & (u > 0.75)
        if c == DECOY[0]:
            dec = (slot >= DECOY[1]) & (slot < DECOY[2])
            crm, tek, gene = crm & ~dec, tek & ~dec, gene | dec
        rows += [(c, a, a + 8000, "CRM", 0) for a in slot[crm]]
        rows += [(c, a, a + 8000, "Tekay", 0) for a in slot[tek]]
        glen = np.where((c == DECOY[0]) & (slot >= DECOY[1]) & (slot < DECOY[2]), 7000, 3000)
        rows += [(c, a, a + g, "gene", 0) for a, g in zip(slot[gene], glen[gene])]
        nest = slot[tek & near][:40]                                   # CRMs nested inside Tekay (depth 1)
        rows += [(c, a + 2000, a + 5000, "CRM", 1) for a in nest]
        if c in CEN:                                                   # satellite arrays, centromere +/- 0.3 Mb
            rows += [(c, a, a + 20_000, "satA", 0) for a in range(int(CEN[c] - 3e5), int(CEN[c] + 3e5), 40_000)]
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "cls", "depth"])


def site_classes(ann, chrom, pos):
    """Truth class per site, the methylation engine's priority: satellite > nested CRM (depth 1) > depth-0 TE > gene > other."""
    a = ann[ann.chrom == chrom]
    out = np.array(["other"] * len(pos), dtype=object)
    for cls, depth in (("gene", 0), ("Tekay", 0), ("CRM", 0), ("CRM", 1), ("satA", 0)):
        x = a[(a.cls == cls) & (a.depth == depth)]
        for s, e in zip(x.start.to_numpy(), x.end.to_numpy()):
            i0, i1 = np.searchsorted(pos, [s, e])
            out[i0:i1] = cls
    return out


def _bedmethyl(path, rows):
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(path, "w") as f:
        for c, p, nv, nm in rows:
            f.write(f"{c}\t{p}\t{p + 1}\tm\t{nv}\t+\t{p}\t{p + 1}\t255,0,0\t{nv}\t{100 * nm / nv:.2f}\t{nm}\t"
                    f"{nv - nm}\t0\t0\t0\t0\t0\n")
    pysam.tabix_index(str(path), preset="bed", force=True)   # writes path.gz + path.gz.tbi, removes path
    return Path(f"{path}.gz")


def write_gff_ltr(path, ann):
    with open(path, "w") as f:
        f.write("##gff-version 3\n")
        for i, r in enumerate(ann[ann.cls.isin(["CRM", "Tekay"])].itertuples()):
            sf = "Gypsy"
            clade = r.cls if i % 25 else "unknown"      # 4% unclassified
            f.write(f"{r.chrom}\tLTRquest\tLTR_retrotransposon\t{r.start + 1}\t{r.end}\t.\t+\t.\tID=e{i};"
                    f"Name={r.chrom}:{r.start + 1}-{r.end};family_rep={r.chrom}:1-2#LTR/{sf}/{r.cls};"
                    f"superfamily={sf};clade={clade};depth={r.depth}\n")
            f.write(f"{r.chrom}\tLTRquest\tlong_terminal_repeat\t{r.start + 1}\t{r.start + 300}\t.\t+\t.\t"
                    f"ID=e{i}.lLTR;Parent=e{i}\n")


def write_arrays(path, ann):
    with open(path, "w") as f:
        f.write("##gff-version 3\n")
        for i, r in enumerate(ann[ann.cls == "satA"].itertuples()):
            f.write(f"{r.chrom}\tcenfuse\tsatellite_DNA\t{r.start + 1}\t{r.end}\t.\t+\t.\tID=sat{i};family=satA;"
                    f"role=centromeric;monomer_len=150;n_monomers=133;mean_identity=0.95;in_centromere=true\n")


def write_genes(path, ann):
    with open(path, "w") as f:
        f.write("##gff-version 3\n")
        for i, r in enumerate(ann[ann.cls == "gene"].itertuples()):
            f.write(f"{r.chrom}\tsim\tgene\t{r.start + 1}\t{r.end}\t.\t+\t.\tID=g{i}\n")
            f.write(f"{r.chrom}\tsim\tmRNA\t{r.start + 1}\t{r.end}\t.\t+\t.\tID=g{i}.t1;Parent=g{i}\n")


def make(d, seed: int = 1, depth: float = 40.0, crm_domain: float = P_CRM_DOMAIN) -> dict:
    """Write the synthetic dataset into folder d; return paths and truth. crm_domain = CHG of CRM cytosines inside
    the centromeric domain (default 0.20: a loss; above 0.45 a gain)."""
    d = Path(d)
    d.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    fai = d / "genome.fa.fai"
    with open(fai, "w") as f:
        for c in CHROMS:
            f.write(f"{c}\t{L}\t0\t60\t61\n")
        f.write("scaf1\t500000\t0\t60\t61\nscaf2\t400000\t0\t60\t61\n")
    ann = annotation(rng)
    write_gff_ltr(d / "ltrquest.gff3", ann)
    write_arrays(d / "g_satellite_arrays.gff", ann)
    write_genes(d / "genes.gff3", ann)
    pos = {c: np.unique(rng.integers(0, L, L // 1000)) for c in CHROMS}
    cls = {c: site_classes(ann, c, pos[c]) for c in CHROMS}
    files = {}
    for smp in SAMPLES:
        for ctx in (("CHG", "CG") if smp == "S1" else ("CHG",)):
            rows = []
            for c in CHROMS:
                p = pos[c]
                k = cls[c]
                prob = np.array([P["sat" if x == "satA" else x] for x in k], float)
                cen = SHIFT[2] if (smp, c) == SHIFT[:2] else CEN.get(c)
                if cen is not None:
                    prob[(k == "CRM") & (np.abs(p - cen) <= DOMAIN / 2)] = crm_domain
                if ctx == "CG":
                    prob = np.clip(prob + 0.4, 0, 0.95)
                nv = rng.poisson(depth, len(p))
                if (smp, c) == DUP[:2]:
                    dup = (p >= DUP[2]) & (p < DUP[3])
                    nv[dup] = rng.poisson(3 * depth, int(dup.sum()))
                    prob[dup] = 0.15
                nm = rng.binomial(nv, prob)
                keep = nv > 0
                rows += list(zip([c] * int(keep.sum()), p[keep].tolist(), nv[keep].tolist(), nm[keep].tolist()))
            files[(smp, ctx)] = _bedmethyl(d / f"{smp}.{ctx}.bedmethyl", rows)
    with open(d / "samples.tsv", "w") as f:
        f.write("sample\tCHG\tCG\n")
        for smp in SAMPLES:
            cg = files.get((smp, "CG"))
            f.write(f"{smp}\t{files[(smp, 'CHG')].name}\t{cg.name if cg else ''}\n")
    with open(d / "anchors.bed", "w") as f:
        f.write("#cenfuse-hic 0.1.1 bin_size=500000\n")
        for c, x in CEN.items():
            x = int(x)
            f.write(f"{c}\t{x - 750000}\t{x + 750000}\t{c}_cen\t500\t.\t{x - 250000}\t{x + 250000}\t0,114,178\n")
    with open(d / "anchors_wrong.bed", "w") as f:
        for c, x in CEN.items():
            y = int(x + 15e6 if x < 30e6 else x - 15e6)
            f.write(f"{c}\t{y - 250000}\t{y + 250000}\t{c}_wrong\n")
    with open(d / "exclude.bed", "w") as f:
        f.write("chr5\t50000000\t52000000\tS2\n")
    return dict(dir=d, fai=fai, samples=d / "samples.tsv", ltrquest=d / "ltrquest.gff3",
                satellites=d / "g_satellite_arrays.gff", genes=d / "genes.gff3", anchors=d / "anchors.bed",
                anchors_wrong=d / "anchors_wrong.bed", exclude=d / "exclude.bed", ann=ann, cen=CEN, shift=SHIFT,
                decoy=DECOY, dup=DUP)
