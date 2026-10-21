"""Tiny simulated Rabl nucleus: .fai + simple 6-column pairs (+ satellite GFF) with known centromeres.

Rabl = two poles plus a tight centromere cluster. Trans weight between bins i (chrom x) and j (chrom y != x):
  1 + Bc * Fc(r_i) * Fc(r_j)      centromere pole (broad):  Fc(r) = exp(-r / 0.3)
    + Bt * Ft(r_i) * Ft(r_j)      telomere pole (broad):    Ft(r) = exp(-(1 - r) / 0.3)
    + Bk * Kc(d_i) * Kc(d_j)      centromere cluster (tight): Kc(d) = exp(-d^2 / 8), d = bins from the centromere
r = distance from the centromere / that arm's length (0 at the centromere, 1 at the end). The broad pole terms make
the iteration converge from far starts (as in real Hi-C data); the tight term gives a sharp core (peaks 3.4-4.1x
the median, as observed in grass Hi-C). Chromosomes in `unattached` (and all, when rabl=False) have weight 1. `hub` =
(chrom index, bin): that bin's trans weights x4. Cis pairs: geometric distance in bins."""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np

SIZE = 250_000
LENGTHS = (100, 120, 90, 110, 96, 104, 88, 116, 92, 108, 100)
CENS = (40, 30, 60, 55, 20, 70, 44, 58, 30, 80, 50)


def polar(nb, cen):
    i = np.arange(nb)
    arm = np.where(i <= cen, max(cen, 1), max(nb - 1 - cen, 1))
    return np.abs(i - cen) / arm


def make_genome(outdir, lengths=LENGTHS, cens=CENS, n_trans=300_000, n_cis=300_000, rabl=True,
                unattached=(10,), hub=(2, 10), sat_shift=(1, 20), rdna_on=3, seed=1, gz=False):
    rng = np.random.default_rng(seed)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    names = [f"chr{k + 1}" for k in range(len(lengths))]
    nb = np.array(lengths)
    off = np.concatenate([[0], np.cumsum(nb)[:-1]])
    N = int(nb.sum())
    chrom_of = np.repeat(np.arange(len(nb)), nb)
    local = np.arange(N) - off[chrom_of]
    r = np.concatenate([polar(n, c) for n, c in zip(nb, cens)])
    d = np.abs(local - np.array(cens)[chrom_of]).astype(float)
    Bc, Bt, Bk = (8.0, 4.0, 10.0) if rabl else (0.0, 0.0, 0.0)
    fc, ft, kc = np.exp(-r / 0.3), np.exp(-(1 - r) / 0.3), np.exp(-d ** 2 / 8)
    w = 1 + Bc * fc[:, None] * fc[None, :] + Bt * ft[:, None] * ft[None, :] + Bk * kc[:, None] * kc[None, :]
    flat = np.isin(chrom_of, list(unattached))
    w[flat, :] = 1.0
    w[:, flat] = 1.0
    hub_glob = None
    if hub is not None:
        hub_glob = int(off[hub[0]] + hub[1])
        w[hub_glob, :] *= 4
        w[:, hub_glob] *= 4
    iu = np.triu_indices(N, 1)
    p = np.where(chrom_of[iu[0]] == chrom_of[iu[1]], 0.0, w[iu])
    p /= p.sum()
    pick = rng.choice(len(p), size=n_trans, p=p)
    bi, bj = iu[0][pick], iu[1][pick]
    ci = rng.choice(N, size=n_cis)
    dd = rng.geometric(0.3, size=n_cis) - 1
    cj = np.clip(local[ci] + np.where(rng.random(n_cis) < 0.5, -1, 1) * dd, 0, nb[chrom_of[ci]] - 1) + off[chrom_of[ci]]
    b1, b2 = np.concatenate([bi, ci]), np.concatenate([bj, cj])
    p1 = local[b1] * SIZE + rng.integers(0, SIZE, len(b1))
    p2 = local[b2] * SIZE + rng.integers(0, SIZE, len(b2))
    lengths_bp = nb * SIZE
    with open(out / "genome.fa.fai", "w") as f:
        o = 0
        for nm, L in zip(names, lengths_bp):
            f.write(f"{nm}\t{L}\t{o}\t60\t61\n")
            o += int(L + L // 60 + 1)
    pairs = out / ("pairs.txt.gz" if gz else "pairs.txt")
    c = np.array(names)
    with (gzip.open if gz else open)(pairs, "wt") as f:
        for c1, q1, c2, q2 in zip(c[chrom_of[b1]], p1, c[chrom_of[b2]], p2):
            f.write(f"{c1}\t{q1}\t+\t{c2}\t{q2}\t-\n")
    sat_bin = list(cens)
    if sat_shift is not None:
        sat_bin[sat_shift[0]] += sat_shift[1]
    with open(out / "sat.gff", "w") as f:
        f.write("##gff-version 3\n")
        for k, nm in enumerate(names):
            s = sat_bin[k] * SIZE
            f.write(f"{nm}\tcenfuse\tcentromere\t{s + 1}\t{s + SIZE}\t.\t.\t.\tID=cen_{nm}\n")
        if rdna_on is not None:
            s = cens[rdna_on] * SIZE
            f.write(f"{names[rdna_on]}\tcenfuse\tsatellite_DNA\t{s + 1}\t{s + 50_000}\t.\t+\t.\t"
                    f"ID=sat_rdna;family=R1;rdna=45S_spacer\n")
    return dict(names=names, cens=list(cens), size=SIZE, unattached=[names[k] for k in unattached],
                hub=hub_glob, sat_bin=sat_bin, rdna_chrom=names[rdna_on] if rdna_on is not None else None,
                pairs=pairs, fai=out / "genome.fa.fai", sat=out / "sat.gff", dir=out)
