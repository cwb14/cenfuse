"""One locus per chromosome: a test of whether a kind of feature behaves like centromeres, in the genome at hand.

A chromosome has one centromere; anything else falls roughly in proportion to length, with long chromosomes holding several and
short ones none. So take the N largest (or deepest) features of the kind, N = the number of chromosomes, and count the chromosomes
they cover: N features that fall in proportion to length cover far fewer than N. This asks nothing about any one locus, needs no
known centromere, and does not care whether the satellite differs from chromosome to chromosome.

On the benchmark, for the N largest clusters of the families called centromeric (cenfuse's satellite line, FASTA only): every
genome with satellite centromeres but the siamang passes (P <= 0.05), down to chromosome-specific satellites (chicken: 36 of 40, null 17), while
tomato, potato, barley and Ae. tauschii fail (P 0.24-0.91), which hold 16 of the FASTA-only confident errors, and so do the
holocentric genomes. The same test on GC troughs is the AT line's (cenfuse.atline)."""
from __future__ import annotations

import numpy as np

from .classify import make_clusters

ALPHA = 0.05            # conventional level
N_NULL = 10_000         # null draws: P resolution 1e-4
SEED = 7                # fixed: reruns give identical P values
MAX_SEQUENCES = 5_000   # more sequences than this are contigs, not chromosomes: no test (and the null grows with them)


def null_covered(lengths, n: int, seed: int = SEED, draws: int = N_NULL) -> np.ndarray:
    """Chromosomes covered by n features that fall in proportion to length, draws times."""
    lengths = np.asarray(lengths, dtype=float)
    p = lengths / lengths.sum()
    x = np.random.default_rng(seed).choice(len(p), size=(draws, n), p=p)
    x.sort(axis=1)
    return 1 + (np.diff(x, axis=1) != 0).sum(axis=1) if n > 1 else np.ones(draws, dtype=np.int64)


def covered_test(top: list, lengths, alpha: float = ALPHA) -> dict:
    """top: the chromosomes of the len(lengths) largest features, largest first (fewer if there are not so many: the null
    then draws that many). Returns dict(n, features, covered, null_median, p, best_possible, status): pass; fail; untested
    (even all chromosomes covered could not pass, or too many sequences to speak of chromosomes)."""
    n, m = len(lengths), len(top)
    covered = len(set(top))
    if n > MAX_SEQUENCES:
        return dict(n=n, features=m, covered=covered, null_median=float("nan"), p=1.0, best_possible=1.0, status="untested")
    full = null_covered(lengths, n)
    null = full if m == n else null_covered(lengths, m) if m else np.zeros(len(full), dtype=np.int64)
    p = (1 + int((null >= covered).sum())) / (1 + len(null))
    best = (1 + int((full >= n).sum())) / (1 + len(full))
    status = "pass" if p <= alpha else "untested" if best > alpha else "fail"
    return dict(n=n, features=m, covered=covered, null_median=float(np.median(null)), p=p, best_possible=best, status=status)


def satellite_test(arrays, families: set, chroms: dict, p) -> dict:
    """The test on the clusters of the given families (pooled per chromosome, as a call's cluster is, with the engine's cluster
    gap). Without such arrays: status 'none'."""
    a = arrays[arrays.family.isin(families)]
    if a.empty:
        return dict(n=len(chroms), features=0, covered=0, null_median=float("nan"), p=1.0, best_possible=1.0, status="none")
    clusters = []
    for c, d in a.groupby("chrom"):
        if c in chroms:
            clusters += [(x.bp, c) for x in make_clusters(list(zip(d.start, d.end)), p.cluster_gap)]
    clusters.sort(key=lambda t: (-t[0], t[1]))
    return covered_test([c for _, c in clusters[:len(chroms)]], list(chroms.values()))


def sentence(t: dict) -> str:
    """'37/40 chromosomes hold the 40 largest clusters (null 17, P 0.0001): pass'"""
    if not t or t.get("status") == "none":
        return "."
    if not np.isfinite(t["null_median"]):
        return f"not tested on {t['n']} sequences"
    m = t.get("features", t["n"])
    what = f"the {m} largest clusters" if m == t["n"] else f"all {m} clusters"
    return f"{t['covered']}/{t['n']} chromosomes hold {what} (null {t['null_median']:.0f}, P {t['p']:.2g}): {t['status']}"
