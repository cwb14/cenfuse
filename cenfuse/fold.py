"""Arm folding: a second Hi-C signature of the centromere, read from each chromosome's own (cis) contact map.

In a Rabl nucleus both arms of a chromosome run from the centromere pole to the telomere pole, so the bin d away
from the centromere on one arm lies beside the bin d away on the other: the cis map carries a ridge along the
anti-diagonal through the centromere (the 'X' of grass Hi-C maps). Trans clustering (the Hi-C engine: centromeres of
different chromosomes meet) and cis folding (the two arms of one chromosome meet) are separate readouts of the
nucleus, so where the two place the centromere at one locus on most chromosomes, beyond chance, each validates the
other without the sequence line. That third line settles the genome-wide dispute fuse.sequence_conflict reports
(barley: satellite calls at the chromosome ends, prior-free Hi-C at the CENH3 domains, fold with Hi-C).

Benchmark (docs/benchmark.md): the fold point lies at the CENH3 domain in rice (11/12 at 500 kb) and
barley (4/7 within 8 Mb) and agrees with the de novo Hi-C peak in barley 6/7 and Poa annua 8/14; it is at chance in
potato and in every T2T mammal (no Rabl fold in those nuclei) and in all randomized Hi-C controls. Its resolution is
the arm scale (core ~5-10% of the chromosome): it can contradict a call at the wrong end of a chromosome, not place
one inside an array.

Rules:
- Fold score of center c (in half bins): mean observed/expected (ICE-balanced, divided by the chromosome's mean at
  that distance) over the anti-diagonal pairs (c - d, c + d) with 2d >= S_MIN, both bins mappable. A center can lie
  inside an unmappable array: its pairs are the array's mappable flanks.
- Centers with a centromeric index below CI_MIN are not scanned: below 12.5% a chromosome is telocentric (Levan et
  al. 1964) and its short arm cannot fold. Without that bound the few short-range pairs near the ends win the
  maximum on noise (barley at 250 kb-1 Mb).
- Fold point = the top center; core = the unbroken run of centers >= median + (peak - median) / 2 (the Hi-C
  engine's half-maximum core rule).
- Validation: a chromosome agrees when its de novo Hi-C peak overlaps the fold core; chance = (core + peak length)
  / scanned length. The fold passes at Poisson-binomial P <= GATE_P with a majority of the tested chromosomes
  agreeing and >= GATE_MIN_REF tested (the Hi-C gate's rules).

Trans arm alignment (0.10.1), the same Rabl axis read from inter-chromosomal contacts: a locus d bins from its
centromere sits about d from the centromere pole on every chromosome, so trans contacts are enriched between loci at
equal distances from their own centromeres (the trans 'X' of Rabl Hi-C maps). With the other chromosomes' centers
at their prior-free Hi-C peaks, chromosome x's alignment score of center c is the mean, over x's mappable bins i at
|i - c| > HALF_WIN, of the row-centered O/E between i and the partner bins at radial distance |i - c|; the
centromere cluster itself is left out on both sides (partner windows, radius <= HALF_WIN, and prior-free cores).
It pools every partner's contacts, so it works where a chromosome's own cis map is too small for the fold (yeasts,
Plasmodium: fold 6/16 and 3/14; alignment 13/16 and 12/14 hold the prior-free peak; randomized Hi-C 0). Center,
core, CI_MIN and validation are the fold's. At a chromosome end centromere- and telomere-centered radial geometry
cannot be told apart (moss: every prior-free peak at a telomere, and without CI_MIN the alignment 'agreed' with them
on 22/26), the reason CI_MIN also applies here.
- Each test passes at P <= GATE_P / 2 (Bonferroni: two readouts of one axis); a peak agrees when it overlaps or
  touches the core, so the chance is (core + 2 peak) / scanned. The genome's axis line is the fold if it passes,
  else the alignment where the fold leans the same way (fuse.axis_line)."""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from .hic.params import HALF_WIN
from .params import GATE_MIN_REF, GATE_P
from .util import gap, log, poisson_binomial_sf

S_MIN = 2 * HALF_WIN + 1   # bins: a pair spans beyond the +/- HALF_WIN centromeric domain, so it measures arm contact
CI_MIN = 0.125             # Levan et al. 1964: centromeric index < 12.5% = telocentric ('t'), no second arm to fold
N_MIN = 3                  # heuristic: fewest anti-diagonal pairs whose mean scores a center


def profile(W: np.ndarray, good: np.ndarray):
    """(F, n) per anti-diagonal t = i + j (center t / 2 in bins) of a dense balanced cis block W (nb x nb)."""
    nb = len(good)
    tot, n = np.zeros(2 * nb - 1), np.zeros(2 * nb - 1)
    for s in range(S_MIN, nb):
        i = np.arange(nb - s)
        ok = good[i] & good[i + s]
        if not ok.any():
            continue
        v = W[i[ok], i[ok] + s]
        e = v.mean()                       # expected at distance s: the chromosome's mean over mappable pairs
        if e > 0:
            t = 2 * i[ok] + s
            np.add.at(tot, t, v / e)
            np.add.at(n, t, 1)
    F = np.full(2 * nb - 1, np.nan)
    np.divide(tot, n, out=F, where=n >= N_MIN)
    return F, n


def call(F: np.ndarray, size: int, length: int) -> dict | None:
    """Fold point and half-maximum core (bp, 0-based half-open) over the centers with CI >= CI_MIN, or None."""
    pos = (np.arange(len(F)) / 2 + 0.5) * size        # center bp: the midpoint of the pair's two bin centers
    return _core(F, pos, size, length, size / 2)


def scan(W, bins, good, members) -> dict:
    """{chromosome: call()} for every callable chromosome; W = the genome's ICE-balanced sparse matrix."""
    out = {}
    for k in np.where(members)[0]:
        sl = bins.sl(int(k))
        F, _ = profile(W[sl, sl].toarray(), good[sl])
        c = call(F, bins.size, int(bins.lengths[k]))
        if c is not None:
            out[bins.names[int(k)]] = c
    return out


def _core(F: np.ndarray, pos: np.ndarray, size: int, length: int, per: float) -> dict | None:
    """Top center and half-maximum core (bp, 0-based half-open) of profile F at centers pos (bp), over the centers
    with CI >= CI_MIN; per = bp each center stands for (scanned length). None when nothing is scanned."""
    F = F.copy()
    F[(pos < CI_MIN * length) | (pos > (1 - CI_MIN) * length)] = np.nan
    ok = np.isfinite(F)
    if not ok.any():
        return None
    j = int(np.nanargmax(F))
    med = float(np.nanmedian(F))
    half = med + (F[j] - med) / 2
    lo = hi = j
    while lo - 1 >= 0 and ok[lo - 1] and F[lo - 1] >= half:
        lo -= 1
    while hi + 1 < len(F) and ok[hi + 1] and F[hi + 1] >= half:
        hi += 1
    return dict(point=int(pos[j]), core_start=max(0, int(pos[lo] - size / 2)),
                core_end=min(length, int(pos[hi] + size / 2)), scanned=int(ok.sum() * per),
                height=float(F[j]), median=med)


def align_profile(W, bins, good, mean, members, centers: dict, k: int, cores: dict | None = None) -> np.ndarray:
    """Alignment score A(c) of every bin c of chromosome k (module doc). W: the genome's balanced sparse matrix (CSR),
    mean: hic.score.block_means, centers: {chromosome index: center bin} of the partners (k's own is ignored),
    cores: {chromosome index: (first, last) local bin of its prior-free core}. The centromere cluster itself is left
    out on both sides, so the score reads arm contacts only: partner bins within HALF_WIN of their center or inside
    their prior-free core, and bins of k within HALF_WIN of c (a cluster alone, with no aligned arms, otherwise
    passes the test: a synthetic +/-4-bin centromere domain did, 15/16)."""
    nb = int(bins.nb[k])
    rad = np.full(bins.N, -1, dtype=np.int64)
    for y, c in centers.items():
        y = int(y)
        if y != k and members[y]:
            loc = np.arange(int(bins.nb[y]))
            r = np.abs(loc - int(c))
            r[r <= HALF_WIN] = -1                                  # the partner window: the cluster, not the arms
            if cores and y in cores:
                lo, hi = cores[y]
                r[(loc >= lo) & (loc <= hi)] = -1                  # ... and the partner's prior-free core
            rad[bins.sl(y)] = r
    m = mean[k, bins.chrom_of]
    valid = (rad >= 0) & good & np.isfinite(m) & (m > 0)
    A = np.full(nb, np.nan)
    if not valid.any():
        return A
    cols = np.flatnonzero(valid)
    R = int(rad[cols].max()) + 1
    ind = sp.csr_matrix((1.0 / m[cols], (np.arange(len(cols)), rad[cols])), shape=(len(cols), R))
    S = np.asarray((W[bins.sl(k)][:, cols] @ ind).todense())       # nb x R: summed O/E by radial distance
    n_r = np.bincount(rad[cols], minlength=R).astype(float)
    P = np.full((nb, R), np.nan)
    np.divide(S, n_r, out=P, where=n_r > 0)
    P[:, :HALF_WIN + 1] = np.nan
    P[~good[bins.sl(k)]] = np.nan
    rows = np.isfinite(P).any(axis=1)
    P[rows] -= np.nanmean(P[rows], axis=1, keepdims=True)
    i = np.arange(nb)
    for c in range(nb):
        d = np.abs(i - c)
        sel = (d > HALF_WIN) & (d < R)
        v = P[i[sel], d[sel]]
        v = v[np.isfinite(v)]
        if len(v) >= N_MIN:
            A[c] = v.mean()
    return A


def align_centers(table, bins) -> tuple:
    """({chromosome index: prior-free peak bin}, {chromosome index: (first, last) local bin of its prior-free core})
    of the callable chromosomes with a prior-free call (hicfuse table): the partners of the alignment."""
    idx = {n: i for i, n in enumerate(bins.names)}
    cen, cores = {}, {}
    for r in table.itertuples(index=False):
        if r.chrom in idx and r.callable and r.denovo_confidence != "none" and r.denovo_peak_start >= 0:
            k = idx[r.chrom]
            cen[k] = int(r.denovo_peak_start) // bins.size
            if r.denovo_core_start >= 0:
                cores[k] = (int(r.denovo_core_start) // bins.size, max(int(r.denovo_core_start),
                                                                         int(r.denovo_core_end) - 1) // bins.size)
    return cen, cores


def align_scan(W, bins, good, mean, members, centers: dict, cores: dict | None = None) -> dict:
    """{chromosome: alignment call} (the fold's call fields) for every callable chromosome with >= 2 centered
    partners; centers and cores as in align_profile."""
    out = {}
    W = W.tocsr()
    for k in np.where(members)[0]:
        k = int(k)
        if sum(1 for y in centers if int(y) != k and members[int(y)]) < 2:
            continue
        A = align_profile(W, bins, good, mean, members, centers, k, cores)
        c = _core(A, (np.arange(len(A)) + 0.5) * bins.size, bins.size, int(bins.lengths[k]), bins.size)
        if c is not None:
            out[bins.names[k]] = c
    return out


def validate(folds: dict, table, label: str = "arm fold") -> dict:
    """Agreement of the fold (or alignment) cores with the de novo Hi-C peaks (module doc). table: the hicfuse
    table. Passes at P <= GATE_P / 2 (Bonferroni over the fold and the alignment)."""
    rows, chance = {}, []
    for r in table.itertuples(index=False):
        f = folds.get(r.chrom)
        if f is None or not r.callable or r.denovo_confidence == "none" or r.denovo_peak_start < 0:
            continue
        agree = gap(r.denovo_peak_start, r.denovo_peak_end, f["core_start"], f["core_end"]) == 0
        rows[r.chrom] = agree
        core, peak = f["core_end"] - f["core_start"], r.denovo_peak_end - r.denovo_peak_start
        # agreement = overlap or touching (gap 0): a peak bin agrees from one bin before the core to one after it
        chance.append(min(1.0, (core + 2 * peak) / max(1, f["scanned"])))
    n, k = len(rows), int(sum(rows.values()))
    p = poisson_binomial_sf(k, chance) if n else 1.0
    status = "untested" if n < GATE_MIN_REF else "pass" if p <= GATE_P / 2 and 2 * k > n else "fail"
    log.info(f"Hi-C {label}: {len(folds)} chromosomes scanned; its core holds the de novo peak on {k}/{n} "
             f"(P {p:.2g}): {status}")
    return dict(agrees=rows, k=k, n=n, p=p, status=status)
