"""Telomere arrays, found de novo: no motif is given.

The genome's telomere unit is the primitive 5-12-bp unit that most often repeats in tandem in the terminal
TERM_SCAN of the chromosome ends, counted on both strands as one canonical unit. Any telomerase RNA template variant
(TTTAGGG in most plants, TTAGGG in vertebrates, TTAGG in insects) is found this way. Each end then gets its array:
the array holding the most unit copies (either strand, short breaks tolerated) in the terminal TERM_SEARCH; its
length, and the gap from the array to the sequence end. Units shorter than 30 bp are outside cenrep's satellite search, so the two
never report the same array."""
from __future__ import annotations

import random
import re
from collections import Counter

import pysam

from .util import log

TERM_SCAN = 20_000         # bp at each end used to learn the unit (telomere arrays are kb-long at assembled ends)
TERM_SEARCH = 100_000      # bp at each end searched for the array (assemblies can leave a short non-telomeric tail)
K_RANGE = range(5, 13)     # telomerase templates make 5-12-bp units (TTAGG 5 ... TTTTAGGGG 9; some algae 8-12)
MAX_BREAK = 30             # copies of the unit closer than this are one run (sequencing errors, variant copies)
CHAIN_GAP = 1_000          # runs closer than this form one array (variant-repeat stretches split real arrays)
MIN_ARRAY = 300            # bp of unit copies in an array: ~40 copies of a 7-bp unit, beyond chance tandem matches
MIN_COPIES = 4             # runs shorter than this (chance matches every ~100 bp in random sequence) are not chained
MIN_DENSITY = 0.5          # an array is mostly unit copies: covered bp >= half of its span
MIN_ENDS = 2               # a telomere unit forms arrays at two or more chromosome ends ...
END_ENRICH = 10            # ... and is >= 10-fold denser in tandem at the ends than inside the chromosomes
INTERIOR = 5               # interior windows (TERM_SCAN bp each) sampled per chromosome for that comparison
SEED = 7


def canonical(u: str) -> str:
    """Smallest rotation of the unit or of its reverse complement (one name for both strands and all phases)."""
    rc = u.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    return min(min(x[i:] + x[:i] for i in range(len(x))) for x in (u, rc))


def primitive(u: str) -> bool:
    """False for units that are themselves repeats of a shorter unit (ATATAT, AAAAA)."""
    k = len(u)
    return not any(k % p == 0 and u == u[:p] * (k // p) for p in range(1, k))


def tandem_units(seq: str) -> Counter:
    """Canonical primitive units counted at every place where a k-mer is followed by the same k-mer."""
    c = Counter()
    for k in K_RANGE:
        for i in range(0, len(seq) - 2 * k + 1):
            u = seq[i:i + k]
            if u == seq[i + k:i + 2 * k] and "N" not in u and primitive(u):
                c[(k, canonical(u))] += 1
    return c


def _ends(fa, chroms: dict, size: int):
    for c, L in chroms.items():
        w = min(size, L // 2)
        yield c, "start", fa.fetch(c, 0, w).upper(), 0
        yield c, "end", fa.fetch(c, L - w, L).upper(), L - w


def _interior(fa, chroms: dict):
    rng = random.Random(SEED)
    for c, L in chroms.items():
        if L > 4 * TERM_SEARCH:
            for _ in range(INTERIOR):
                s = rng.randrange(TERM_SEARCH, L - TERM_SEARCH - TERM_SCAN)
                yield fa.fetch(c, s, s + TERM_SCAN).upper()


def runs(seq: str, unit: str):
    """[start, end, bp covered] runs of the unit (either strand, any phase), copies within MAX_BREAK merged. Every
    rotation is searched, so a tandem array matches at each of its bases; coverage is the union of the matches."""
    k = len(unit)
    rc = unit.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    rots = {x[i:] + x[:i] for x in (unit, rc) for i in range(k)}
    out = []
    for m in re.finditer("(?=(" + "|".join(sorted(rots)) + "))", seq):
        h = m.start()
        if out and h - out[-1][1] <= MAX_BREAK:
            if h > out[-1][1]:
                out[-1][2] += k                       # a new copy after a break
            else:
                out[-1][2] += max(0, h + k - out[-1][1])   # overlapping match: only the new bases
            out[-1][1] = max(out[-1][1], h + k)
        else:
            out.append([h, h + k, k])
    return out


def end_array(seq: str, unit: str, side: str):
    """The telomere array of one end: runs of >= MIN_COPIES chained while gaps <= CHAIN_GAP into arrays (>= MIN_ARRAY
    bp of unit copies, >= MIN_DENSITY of their span), and the array holding the most unit bp kept (a short decoy run
    nearer the end must not hide the real array). (start, end, unit bp) in seq coordinates, or None."""
    k = len(unit)
    rs = [r for r in runs(seq, unit) if r[2] >= MIN_COPIES * k]
    if not rs:
        return None
    arrays, cur = [], list(rs[0])
    for s2, e2, b2 in rs[1:]:
        if s2 - cur[1] > CHAIN_GAP:
            arrays.append(cur)
            cur = [s2, e2, b2]
        else:
            cur[1], cur[2] = max(cur[1], e2), cur[2] + b2
    arrays.append(cur)
    ok = [a for a in arrays if a[2] >= MIN_ARRAY and a[2] >= MIN_DENSITY * (a[1] - a[0])]
    if not ok:
        return None
    best = max(ok, key=lambda a: (a[2], a[1] if side == "end" else -a[0]))   # most unit bp; ties: nearer the end
    return tuple(best)


def detect(fasta, chroms: dict) -> dict:
    """{'unit': canonical unit or None, 'enrichment': fold, 'ends': [(chrom, side, start, end, unit_bp,
    gap_to_end)]} (0-based)."""
    with pysam.FastaFile(str(fasta)) as fa:
        term = Counter()
        for _, _, s, _ in _ends(fa, chroms, TERM_SCAN):
            term.update(tandem_units(s))
        if not term:
            return dict(unit=None, enrichment=float("nan"), ends=[])
        inner, n_inner = Counter(), 0
        for s in _interior(fa, chroms):
            inner.update(tandem_units(s))
            n_inner += 1
        n_term = 2 * len(chroms)
        ranked = sorted(term, key=lambda x: -term[x] * x[0])      # bp covered in tandem: long units not penalised
        for key in ranked[:10]:
            k, unit = key
            if term[key] * k < MIN_ARRAY:                      # too little tandem sequence at the ends to judge
                continue
            enrich = (term[key] / n_term) / max(inner[key] / max(n_inner, 1), 1 / max(n_inner, 1))
            if enrich < END_ENRICH:
                continue
            rows = []
            for c, side, s, off in _ends(fa, chroms, TERM_SEARCH):
                a = end_array(s, unit, side)
                if a:
                    gap_end = off + a[0] if side == "start" else chroms[c] - (off + a[1])
                    rows.append((c, side, off + a[0], off + a[1], a[2], gap_end))
            if len(rows) >= MIN_ENDS:
                log.info(f"telomeres: unit {unit} ({enrich:.0f}-fold denser at the ends); arrays at {len(rows)}/"
                         f"{n_term} chromosome ends")
                return dict(unit=unit, enrichment=round(enrich, 1), ends=rows)
    log.info("telomeres: no 5-12-bp unit is enriched in tandem arrays at the chromosome ends")
    return dict(unit=None, enrichment=float("nan"), ends=[])
