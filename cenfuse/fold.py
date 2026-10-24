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
  agreeing and >= GATE_MIN_REF tested (the Hi-C gate's rules)."""
from __future__ import annotations

import numpy as np

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
    F = F.copy()
    pos = (np.arange(len(F)) / 2 + 0.5) * size        # center bp: the midpoint of the pair's two bin centers
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
                core_end=min(length, int(pos[hi] + size / 2)), scanned=int(ok.sum() * size / 2),
                height=float(F[j]), median=med)


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


def validate(folds: dict, table) -> dict:
    """Agreement of the fold cores with the de novo Hi-C peaks (module doc). table: the hicfuse table."""
    rows, chance = {}, []
    for r in table.itertuples(index=False):
        f = folds.get(r.chrom)
        if f is None or not r.callable or r.denovo_confidence == "none" or r.denovo_peak_start < 0:
            continue
        agree = gap(r.denovo_peak_start, r.denovo_peak_end, f["core_start"], f["core_end"]) == 0
        rows[r.chrom] = agree
        core, peak = f["core_end"] - f["core_start"], r.denovo_peak_end - r.denovo_peak_start
        chance.append(min(1.0, (core + peak) / max(1, f["scanned"])))
    n, k = len(rows), int(sum(rows.values()))
    p = poisson_binomial_sf(k, chance) if n else 1.0
    status = "untested" if n < GATE_MIN_REF else "pass" if p <= GATE_P and 2 * k > n else "fail"
    log.info(f"Hi-C arm fold: {len(folds)} chromosomes scanned; fold core holds the de novo peak on {k}/{n} "
             f"(P {p:.2g}): {status}")
    return dict(agrees=rows, k=k, n=n, p=p, status=status)
