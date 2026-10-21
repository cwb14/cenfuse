"""Synthetic GENOME_DIR for cenfuse at the engines' real scales (Hi-C bins 250 kb, 100 per chromosome, enough for
the Hi-C engine's null; methylation D ~ 270 kb): ten 25-Mb chromosomes, satellite-engine-format sequence calls (imported with --satellite-dir), LTRquest-style TEs, genes,
Hi-C pairs from a Rabl nucleus, and CHG bedMethyl for three plants.

Functional centromere F (Hi-C cluster; CRM elements with CHG loss) vs satellite call S:
- chr1-chr6: S = F (concordant).
- chr7: S at 3.75 Mb (high), F at 18.75 Mb. CRMs are dense at S (the TE core) and sparser at F, CHG loss only at F:
  methylation is the second line, so the call moves.
- chr8: S at 5 Mb (high), F at 18.75 Mb, the CRM core at F: TE is the second line, so the call moves.
- chr9: S at 6.25 Mb (high), F at 20 Mb, the CRM core at S and no CRM at F: no second line, so `contested`.
- chr10: no satellite call, F at 12.5 Mb with its CRM core: `functional_only`.
Satellite calls are 400 kb of 20-kb arrays every 40 kb (CRMs fall in the gaps, as in real arrays). Other
features: a subtelomeric family at every chromosome end, an interstitial knob on chr3, an rDNA-tagged family on chr5,
a Hi-C contact hub on chr2, and a 50-kb N gap inside chr4's centromere.
permute=True moves every read end to a uniform random place on its own chromosome (negative control).

Plants (docs/plants.md). The sequence is random A/T with a CAG at every methylation site (a CHG cytosine on each
strand, rows on both strands every 300 bp), so reads that show another base remove rows as modkit does:
- P1 is the reference plant: no fixed differences.
- P2 and P3 are other plants: fixed differences at 1.2% of sites each, plus 0.4% shared (related plants).
- P2 lacks chr6 15-21 Mb (no rows: a missing segment).
- P3's chr3 centromere sits at 20 Mb (SHIFT), on a second satA array where the reference's CRMs are sparse; its
  CHG loss and its own Hi-C library (hic/Syn.tsv: libraries 'hic', the reference plant's, and 'P3') follow it.
Every row is also lost alone with probability 0.02 (failed calls, read errors)."""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pysam

L = 25_000_000
BIN = 250_000
NB = L // BIN
CHROMS = [f"chr{i}" for i in range(1, 11)]
F = {"chr1": 25, "chr2": 38, "chr3": 50, "chr4": 62, "chr5": 75, "chr6": 30, "chr7": 75, "chr8": 75, "chr9": 80,
     "chr10": 50}                                           # functional centromere bin
S = {**{c: F[c] for c in CHROMS[:6]}, "chr7": 15, "chr8": 20, "chr9": 25}  # satellite call bin (none on chr10)
CRM_CORE = {**F, "chr7": 15, "chr9": 25}                   # densest CRM cluster
CRM_EXTRA = {"chr7": 75}                                   # sparser CRMs (still hypomethylated) at F
HUB = ("chr2", 90)
KNOB = ("chr3", 4_000_000, 4_200_000)
RDNA = ("chr5", 5_000_000, 5_150_000)
DISP = [("chr2", 3_000_000), ("chr4", 21_000_000), ("chr6", 18_000_000), ("chr8", 11_000_000)]   # no 3D clustering
GAP = ("chr4", 15_600_000, 15_650_000)
SAMPLES = ("P1", "P2", "P3")
SHIFT = ("P3", "chr3", 80)                                 # plant, chromosome, its centromere bin (20 Mb)
MISSING = ("P2", "chr6", 15_000_000, 21_000_000)           # a segment the plant lacks
SITE_STEP = 300                                            # CHG sites (CAG) every 300 bp, rows on both strands
FIXED = {"P2": 0.012, "P3": 0.012}                         # private fixed differences per site
SHARED_FIXED = 0.004                                       # ... shared by P2 and P3
ROW_LOSS = 0.02                                            # a row lost alone (noise)
P = {"sat": 0.68, "CRM": 0.45, "CRM_cen": 0.15, "Tekay": 0.50, "gene": 0.10, "other": 0.45}
DOMAIN = 300_000                                           # hypomethylated domain width around F (bp)


def bp(b: int) -> int:
    return b * BIN + BIN // 2


def sites() -> np.ndarray:
    """+ strand C of every CAG (the - strand C is 2 bp further)."""
    return np.arange(150, L - 3, SITE_STEP)


def write_fasta(path: Path, rng) -> None:
    """Random A/T with a CAG at every methylation site: no other C or G, so the only CHG cytosines are the sites."""
    pos = sites()
    with open(path, "wb") as f, open(f"{path}.fai", "w") as fai:
        off = 0
        for c in CHROMS:
            seq = np.where(rng.random(L) < 0.5, ord("A"), ord("T")).astype(np.uint8)
            seq[pos], seq[pos + 1], seq[pos + 2] = ord("C"), ord("A"), ord("G")
            if c == GAP[0]:
                seq[GAP[1]:GAP[2]] = ord("N")
            head = f">{c}\n".encode()
            off += len(head)
            fai.write(f"{c}\t{L}\t{off}\t{L}\t{L + 1}\n")
            f.write(head + seq.tobytes() + b"\n")
            off += L + 1


def write_satellite(sat: Path, genome: str, plants: bool = False) -> None:
    """Satellite-engine outputs: summary (1-based), arrays GFF (roles, rdna tag), families TSV. plants: a second satA array at
    the shifted plant's centromere (200 kb at chr3 20 Mb)."""
    sat.mkdir(parents=True, exist_ok=True)
    cols = ["chrom", "chrom_len", "cen_start", "cen_end", "cen_len", "evidence", "confidence", "flags", "main_family",
            "main_family_bp", "cen_sat_bp", "sat_fraction", "n_arrays", "largest_array_bp", "n_monomers",
            "mean_identity", "hor_monomers", "hor_fraction", "top_tandem_hor", "n_gaps", "gap_bp", "n_genes",
            "gene_ratio", "other_tandem_bp", "te_clade", "te_core_start", "te_core_end"]
    arrays = []
    with open(sat / f"{genome}_centromere_summary.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for c in CHROMS:
            if c in S:
                s, e = S[c] * BIN, S[c] * BIN + 400_000
                row = [c, L, s + 1, e, e - s, "satellite", "high", ".", "satA", 200_000, 200_000, 0.5, 10,
                       20_000, 1300, 0.9, 0, 0, ".", 0, 0, ".", ".", 0, ".", ".", "."]
                arrays += [(c, a, a + 20_000, "satA", "centromeric", ".") for a in range(s, e, 40_000)]
            else:
                row = [c, L] + ["."] * 3 + ["none", "none"] + ["."] * 20
            f.write("\t".join(map(str, row)) + "\n")
            arrays += [(c, 0, 60_000, "subT", "subtelomeric", "."), (c, L - 60_000, L, "subT", "subtelomeric", ".")]
    if plants:
        s0 = SHIFT[2] * BIN
        arrays += [(SHIFT[1], a, a + 20_000, "satA", "centromeric", ".") for a in range(s0, s0 + 200_000, 40_000)]
    arrays.append((KNOB[0], KNOB[1], KNOB[2], "knob", "interstitial", "."))
    arrays.append((RDNA[0], RDNA[1], RDNA[2], "rdnaR", "rDNA", "45S_spacer"))
    arrays += [(c, s, s + 60_000, "dispR", "dispersed", ".") for c, s in DISP]
    with open(sat / f"{genome}_satellite_arrays.gff", "w") as f:
        f.write("##gff-version 3\n")
        for i, (c, s, e, fam, role, rd) in enumerate(sorted(arrays), 1):
            extra = f";rdna={rd}" if rd != "." else ""
            f.write(f"{c}\tcenfuse\tsatellite_DNA\t{s + 1}\t{e}\t.\t+\t.\tID=sat_{i};family={fam};role={role}{extra};"
                    f"monomer_len=150;n_monomers={(e - s) // 150};mean_identity=0.9;in_centromere=false\n")
    with open(sat / "satellite_families.tsv", "w") as f:
        f.write("family\tgenome\tclass\trole\n")
        for fam, role in (("satA", "centromeric"), ("subT", "subtelomeric"), ("knob", "interstitial"),
                          ("rdnaR", "rDNA"), ("dispR", "dispersed")):
            f.write(f"{fam}\t{genome}\t{'centromeric' if role == 'centromeric' else 'other'}\t{role}\n")


def te_elements(rng, plants: bool = False) -> list:
    """(chrom, start, end, clade): CRM clusters (core + extra), sparse CRM, uniform Tekay. plants: extra CRMs at the
    shifted plant's centromere."""
    out = []
    extra = {**CRM_EXTRA, SHIFT[1]: SHIFT[2]} if plants else CRM_EXTRA
    for c in CHROMS:
        for center, n in ((CRM_CORE.get(c), 40), (extra.get(c), 14)):
            if center is not None:
                for m in rng.normal(bp(center), 120_000, n):
                    m = int(np.clip(m, 5_000, L - 5_000))
                    out.append((c, m - 4_000, m + 4_000, "CRM"))
        ring = rng.integers(0, L - 8_000, 100)                                                 # ring references
        if c == "chr9":                                    # by design no CRM at chr9's F: no methylation second line
            ring = ring[np.abs(ring - bp(F[c])) > 2 * DOMAIN]
        out += [(c, int(m), int(m) + 8_000, "CRM") for m in ring]
        out += [(c, int(m), int(m) + 8_000, "Tekay") for m in rng.integers(0, L - 8_000, 300)]
    return out


def write_tes(path: Path, elems: list) -> None:
    with open(path, "w") as f:
        f.write("##gff-version 3\n")
        for i, (c, s, e, clade) in enumerate(sorted(elems), 1):
            f.write(f"{c}\tLTRquest\tLTR_retrotransposon\t{s + 1}\t{e}\t.\t+\t.\tID=te{i};clade={clade};"
                    f"superfamily={'Gypsy'};depth=0\n")


def write_genes(path: Path, rng) -> list:
    genes = []
    with open(path, "w") as f:
        f.write("##gff-version 3\n")
        for c in CHROMS:
            for i, m in enumerate(sorted(rng.integers(100_000, L - 100_000, 350))):
                if any(abs(m - bp(b)) < 600_000 for b in {F[c], S.get(c, F[c])}):
                    continue
                genes.append((c, int(m), int(m) + 3_000))
                f.write(f"{c}\tsynth\tgene\t{m + 1}\t{m + 3000}\t.\t+\t.\tID={c}_g{i}\n")
    return genes


def write_pairs(path: Path, rng, n_trans=400_000, n_cis=400_000, permute=False, cen_bins=None) -> None:
    """The Hi-C engine's Rabl model (tests/hic/hic_synth.py) at the functional centromeres F (or cen_bins), plus a hub."""
    names = CHROMS
    N = NB * len(names)
    chrom_of = np.repeat(np.arange(len(names)), NB)
    local = np.arange(N) % NB
    cen = np.array([(cen_bins or F)[c] for c in names])[chrom_of]
    arm = np.where(local <= cen, np.maximum(cen, 1), np.maximum(NB - 1 - cen, 1))
    r = np.abs(local - cen) / arm
    d = np.abs(local - cen).astype(float)
    fc, ft, kc = np.exp(-r / 0.3), np.exp(-(1 - r) / 0.3), np.exp(-d ** 2 / 8)
    w = 1 + 8 * fc[:, None] * fc[None, :] + 4 * ft[:, None] * ft[None, :] + 10 * kc[:, None] * kc[None, :]
    hub = names.index(HUB[0]) * NB + HUB[1]
    w[hub, :] *= 4
    w[:, hub] *= 4
    iu = np.triu_indices(N, 1)
    p = np.where(chrom_of[iu[0]] == chrom_of[iu[1]], 0.0, w[iu])
    p /= p.sum()
    pick = rng.choice(len(p), size=n_trans, p=p)
    bi, bj = iu[0][pick], iu[1][pick]
    ci = rng.choice(N, size=n_cis)
    dd = rng.geometric(0.3, size=n_cis) - 1
    cj = np.clip(local[ci] + np.where(rng.random(n_cis) < 0.5, -1, 1) * dd, 0, NB - 1) + chrom_of[ci] * NB
    b1, b2 = np.concatenate([bi, ci]), np.concatenate([bj, cj])
    p1 = local[b1] * BIN + rng.integers(0, BIN, len(b1))
    p2 = local[b2] * BIN + rng.integers(0, BIN, len(b2))
    if permute:
        p1, p2 = rng.integers(0, L, len(b1)), rng.integers(0, L, len(b2))
    c = np.array(names)
    with gzip.open(path, "wt") as f:
        for c1, q1, c2, q2 in zip(c[chrom_of[b1]], p1, c[chrom_of[b2]], p2):
            f.write(f"{c1}\t{q1}\t+\t{c2}\t{q2}\t-\n")


def site_class(c: str, pos: np.ndarray, sat_iv, crm_iv, gene_iv) -> np.ndarray:
    out = np.full(len(pos), "other", dtype=object)
    for cls, ivs in (("gene", gene_iv), ("CRM", crm_iv), ("sat", sat_iv)):
        for s, e in ivs.get(c, []):
            i0, i1 = np.searchsorted(pos, [s, e])
            out[i0:i1] = cls
    return out


def write_meth(mdir: Path, rng, elems, genes, plants: bool = False) -> Path:
    """CHG bedMethyl per plant, rows on both strands of every CAG: CRM cytosines lose CHG inside the domain around the
    plant's centromere; the satellite stays CHG-high. plants: genotypes, P2's missing segment and P3's shift (module
    doc). The table lists the bgzip+tabix files relative to meth/."""
    mdir.mkdir(parents=True, exist_ok=True)
    sat_iv = {}
    for c in CHROMS:
        if c in S:
            sat_iv.setdefault(c, []).extend((a, a + 20_000) for a in range(S[c] * BIN, S[c] * BIN + 400_000, 40_000))
    crm_iv, gene_iv = {}, {}
    for c, s, e, clade in elems:
        if clade == "CRM":
            crm_iv.setdefault(c, []).append((s, e))
    for c, s, e in genes:
        gene_iv.setdefault(c, []).append((s, e))
    pos = sites()
    shared = {c: rng.random(len(pos)) < SHARED_FIXED for c in CHROMS}
    rows = ["sample\tCHG"]
    for smp in SAMPLES:
        path = mdir / f"{smp}.CHG.bedmethyl"
        with open(path, "w") as f:
            for c in CHROMS:
                cls = site_class(c, pos, sat_iv, crm_iv, gene_iv)
                pm = np.array([P[k] for k in cls])
                cen = SHIFT[2] if plants and (smp, c) == SHIFT[:2] else F[c]
                dom = np.abs(pos - bp(cen)) <= DOMAIN // 2
                if c in CRM_EXTRA:
                    dom |= np.abs(pos - bp(CRM_EXTRA[c])) <= DOMAIN // 2
                pm = np.where((cls == "CRM") & dom, P["CRM_cen"], pm)
                keep = np.ones(len(pos), bool)
                if c == GAP[0]:
                    keep &= (pos + 2 < GAP[1]) | (pos >= GAP[2])
                if plants and smp in FIXED:
                    keep &= ~(shared[c] | (rng.random(len(pos)) < FIXED[smp]))
                if plants and (smp, c) == MISSING[:2]:
                    keep &= (pos + 2 < MISSING[2]) | (pos >= MISSING[3])
                out = []
                for strand, q in (("+", pos), ("-", pos + 2)):
                    ok = keep & (rng.random(len(pos)) >= ROW_LOSS)
                    nv = rng.poisson(14, len(pos)) + 3
                    nm = rng.binomial(nv, pm)
                    out += [(int(x), strand, int(a), int(b)) for x, a, b in zip(q[ok], nv[ok], nm[ok])]
                out.sort()
                f.write("".join(f"{c}\t{x}\t{x + 1}\tm\t{a}\t{s}\t{x}\t{x + 1}\t255,0,0\t{a}\t{100 * b / a:.2f}\t"
                                f"{b}\t{a - b}\t0\t0\t0\t0\t0\n" for x, s, a, b in out))
        pysam.tabix_index(str(path), preset="bed", force=True)
        rows.append(f"{smp}\t{path.name}.gz")
    (mdir / "Syn.tsv").write_text("\n".join(rows) + "\n")
    return mdir / "Syn.tsv"


def make(root, permute: bool = False, meth: bool = True, seed: int = 3, plants: bool = True) -> dict:
    """plants (with meth): the planted plants of the module doc and P3's Hi-C library."""
    root = Path(root)
    g = root / "Syn"
    for sub in ("tes", "genes", "hic", "meth"):
        (g / sub).mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    plants = plants and meth
    write_fasta(g / "Syn.fa", rng)
    write_satellite(root / "sat", "Syn", plants)
    elems = te_elements(rng, plants)
    write_tes(g / "tes" / "Syn.gff3", elems)
    genes = write_genes(g / "genes" / "Syn.gff3", rng)
    write_pairs(g / "hic" / "Syn.pairs.gz", rng, permute=permute)
    if plants:
        write_pairs(g / "hic" / "Syn.P3.pairs.gz", rng, cen_bins={**F, SHIFT[1]: SHIFT[2]})
        (g / "hic" / "Syn.tsv").write_text("sample\tfile\nhic\tSyn.pairs.gz\nP3\tSyn.P3.pairs.gz\n")
    if meth:
        write_meth(g / "meth", rng, elems, genes, plants)
    return dict(dir=g, sat=root / "sat", F=F, S=S, bin=BIN)
