"""Expected centromere size D, the analyzed chromosomes, and the composition-adjusted track: each class is compared
with its own methylation in a ring 1.5-4 D away, deviations are summed over classes and a 0.4-D sliding window."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import params


def d_formula(genome_size: float, n_chrom: int) -> float:
    """Expected CENH3 domain per centromere (bp): published scaling of CENH3 area with genome size, anchored on
    maize Mo17: D = 22.2 Mb x (G / 2.18 Gb)^0.97 / n."""
    return params.D_REF_BP * (genome_size / params.D_REF_G) ** params.D_SLOPE / n_chrom


def _half_index(L: np.ndarray) -> int:
    """Index of the sequence that takes the running total of L (sorted descending) to half the assembly."""
    return int(np.searchsorted(np.cumsum(L), L.sum() / 2))


def tenfold_count(lengths) -> int:
    """Chromosomes in a chromosome-level assembly by the rule of cenfuse 0.8: the sequences above the first >= GAP_FOLD
    drop in sorted length once half the assembly is accounted for (chromosomes hold most of it; unplaced scaffolds are
    far shorter). Every sequence when there is no such drop. Kept for the --classify-from targets and for old output
    folders; a gradual scaffold tail has no such drop, so it takes the whole tail (chromosome_count)."""
    L = np.sort(np.asarray(list(lengths), float))[::-1]
    half = _half_index(L)
    drop = L[:-1] >= params.GAP_FOLD * L[1:]
    drop[:half] = False
    hit = np.flatnonzero(drop)
    return int(hit[0]) + 1 if hit.size else len(L)


def _median_cut(L: np.ndarray) -> int:
    """First index k, at or after the half index (and not 0), whose length is under 1/GAP_FOLD of the median of the k
    longer sequences; len(L) when there is none. L is sorted descending, so each median is read off two positions
    ((L[(k-1)//2] + L[k//2]) / 2, what np.median gives): one pass for assemblies of 10^6 sequences."""
    k = np.arange(max(_half_index(L), 1), len(L))
    hit = np.flatnonzero(params.GAP_FOLD * L[k] < (L[(k - 1) // 2] + L[k // 2]) / 2)
    return int(k[hit[0]]) if hit.size else len(L)


def chromosome_count(lengths) -> int:
    """Chromosomes in a chromosome-level assembly (cenfuse 0.9): the tenfold_count set, unless its small members form a
    gradual scaffold tail. Small = under 1/GAP_FOLD of the median length of the sequences before them (the first is at
    _median_cut); they are a tail when they span a GAP_FOLD range or more (no clean gap separates them), and the set
    then ends at the first of them. A few small sequences of similar size before a clean gap stay (dot chromosomes:
    D. melanogaster chr4, chicken microchromosomes; the median cut alone drops them).
    Why: 383 assemblies: the 0.8 set on all 234 benchmark and batch genomes; on 121 polyploid assemblies with project
    chromosome lists, extra sequences 5,614 -> 31, none missing."""
    L = np.sort(np.asarray(list(lengths), float))[::-1]
    n0, n2 = tenfold_count(L), _median_cut(L)
    if n2 >= n0:
        return n0
    small = L[n2:n0]
    return n0 if small.max() < params.GAP_FOLD * small.min() else n2


def chromosome_set(fai: dict[str, int], n: int | None = None) -> tuple[list[str], float, float]:
    """Analyzed chromosomes = sequences >= MIN_LEN_D x D, where D = d_formula(G, n), G = all sequences in the .fai
    and n = chromosome_count unless given (cenfuse passes the size of its own chromosome set: a --classify-from
    target's 0.8 set, a set the user listed). Returns (names in .fai order, G, D); names is empty when nothing is
    long enough."""
    G = float(sum(fai.values()))
    D = d_formula(G, chromosome_count(fai.values()) if n is None else n)
    return [c for c, L in fai.items() if L >= params.MIN_LEN_D * D], G, D


@dataclass(frozen=True)
class Geometry:
    D: float     # expected domain size (bp)
    w: int       # sliding window (bins)
    lo: int      # ring inner edge (bins)
    hi: int      # ring outer edge (bins)
    sep: int     # SEP_D x D (bins)


def geometry(D: float) -> Geometry:
    b = params.BIN
    w = max(1, round(params.WIN_D * D / b))
    lo = max(w // 2 + 1, round(params.RING_LO_D * D / b))
    hi = max(lo + 1, round(params.RING_HI_D * D / b))
    return Geometry(D, w, lo, hi, max(1, round(params.SEP_D * D / b)))


def _cs(x: np.ndarray) -> np.ndarray:
    return np.concatenate([np.zeros(x.shape[:-1] + (1,)), np.cumsum(x, axis=-1)], axis=-1)


def ring_sum(x: np.ndarray, lo: int, hi: int) -> np.ndarray:
    """Along the last axis: sum of x over j with lo <= |i - j| <= hi."""
    n = x.shape[-1]
    cs, i = _cs(x), np.arange(n)
    right = cs[..., np.clip(i + hi + 1, 0, n)] - cs[..., np.clip(i + lo, 0, n)]
    left = cs[..., np.clip(i - lo + 1, 0, n)] - cs[..., np.clip(i - hi, 0, n)]
    return left + right


def slide(x: np.ndarray, w: int) -> np.ndarray:
    """Along the last axis: sum over the w-bin window i - (w-1)//2 .. i + w//2."""
    n = x.shape[-1]
    cs, i = _cs(x), np.arange(n)
    return cs[..., np.clip(i + w // 2 + 1, 0, n)] - cs[..., np.clip(i - (w - 1) // 2, 0, n)]


def class_sums(v: np.ndarray, m: np.ndarray, geo: Geometry) -> tuple[np.ndarray, np.ndarray]:
    """v, m: (classes, bins). Per class: windowed sum of (m - v x ring reference) and of the calls that had a
    reference. Classes without data in the ring contribute nothing."""
    rv, rm = ring_sum(v, geo.lo, geo.hi), ring_sum(m, geo.lo, geo.hi)
    ok = rv > 0
    ref = np.divide(rm, rv, out=np.zeros_like(rm), where=ok)
    dev = np.where(ok, m - v * ref, 0.0)
    return slide(dev, geo.w), slide(np.where(ok, v, 0.0), geo.w)


def ratio(num: np.ndarray, vok: np.ndarray) -> np.ndarray:
    """Deviation (methylated share above or below the class references) where >= MIN_CALLS calls, else NaN."""
    return np.divide(num, vok, out=np.full(np.shape(num), np.nan), where=np.asarray(vok) >= params.MIN_CALLS)


def adjusted(num: np.ndarray, vok: np.ndarray) -> np.ndarray:
    """Composition-adjusted track = sum over classes of deviations / sum of calls (classes x bins -> bins)."""
    return ratio(num.sum(axis=0), vok.sum(axis=0))
