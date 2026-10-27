"""Telomere arrays, found de novo: no motif is given.

The genome's telomere unit is the primitive 5-12-bp unit that most often repeats in tandem in the terminal
TERM_SCAN of the chromosome ends, counted on both strands as one canonical unit. Any telomerase RNA template variant
(TTTAGGG in most plants, TTAGGG in vertebrates, TTAGG in insects) is found this way. Each end then gets its array:
the array holding the most unit copies (either strand, short breaks tolerated) in the terminal TERM_SEARCH; its
length, and the gap from the array to the sequence end. Units shorter than 30 bp are outside the satellite engine's search, so the two
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
MIN_TERMINAL = 100         # ... or this many bp when the array reaches the sequence end (within TERMINAL): a short telomere
                           # (Neurospora NC12: 157-207 bp at 13/14 ends) or one the assembler trimmed. The unit is already
                           # the genome's telomere unit (END_ENRICH), and 100 bp of one exact 5-12-bp unit is no chance match
TERMINAL = 1_000           # an array this close to the sequence end reaches it
MIN_COPIES = 4             # runs shorter than this (chance matches every ~100 bp in random sequence) are not chained
SCAN_CHUNK = 10_000_000    # bp read at once in the interstitial scan (memory)
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


def chained(seq: str, unit: str) -> list:
    """[start, end, unit bp] of every array: runs of >= MIN_COPIES copies chained while gaps <= CHAIN_GAP."""
    rs = [r for r in runs(seq, unit) if r[2] >= MIN_COPIES * len(unit)]
    if not rs:
        return []
    arrays, cur = [], list(rs[0])
    for s2, e2, b2 in rs[1:]:
        if s2 - cur[1] > CHAIN_GAP:
            arrays.append(cur)
            cur = [s2, e2, b2]
        else:
            cur[1], cur[2] = max(cur[1], e2), cur[2] + b2
    arrays.append(cur)
    return arrays


def end_array(seq: str, unit: str, side: str):
    """The telomere array of one end: chained runs (>= MIN_ARRAY bp of unit copies, or >= MIN_TERMINAL bp for an array
    within TERMINAL of the sequence end; >= MIN_DENSITY of their span), and the array holding the most unit bp kept (a
    short decoy run nearer the end must not hide the real array). (start, end, unit bp) in seq coordinates, or None."""
    to_end = (lambda a: a[0]) if side == "start" else (lambda a: len(seq) - a[1])
    ok = [a for a in chained(seq, unit) if a[2] >= MIN_DENSITY * (a[1] - a[0])
          and (a[2] >= MIN_ARRAY or (a[2] >= MIN_TERMINAL and to_end(a) <= TERMINAL))]
    if not ok:
        return None
    best = max(ok, key=lambda a: (a[2], a[1] if side == "end" else -a[0]))   # most unit bp; ties: nearer the end
    return tuple(best)


def interstitial(fasta, chroms: dict, unit: str) -> list:
    """[(chrom, start, end, unit bp)] of interstitial telomeric sequences (ITS): arrays of the genome's telomere unit
    inside the chromosomes (beyond the terminal TERM_SEARCH), by the end-array rules (>= MIN_ARRAY bp of unit copies at
    >= MIN_DENSITY). Relics of chromosome fusions (human 2q13) and the cores of some centromeres (chicken: 78-210 kb at
    the CENP-A domain of several chromosomes). Only stretches with two tandem copies in one phase are chained (any
    three copies, in any phase, contain them), so a 3-Gb genome takes about a minute."""
    k = len(unit)
    rc = unit.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    seeds = (unit * 2, rc * 2)
    out = []
    with pysam.FastaFile(str(fasta)) as fa:
        for c, L in chroms.items():
            lo_c, hi_c = TERM_SEARCH, L - TERM_SEARCH
            hits = []
            for s in range(lo_c, max(lo_c, hi_c), SCAN_CHUNK):
                seq = fa.fetch(c, s, min(hi_c, s + SCAN_CHUNK + 2 * k)).upper()
                hits += [s + i for x in seeds for i in _find_all(seq, x)]
            hits.sort()
            stretches = []                                     # seeds closer than CHAIN_GAP: one stretch
            for h in hits:
                if stretches and h - stretches[-1][1] <= CHAIN_GAP:
                    stretches[-1][1] = h + 2 * k
                else:
                    stretches.append([h, h + 2 * k])
            for a0, a1 in stretches:
                lo, hi = max(lo_c, a0 - CHAIN_GAP), min(hi_c, a1 + CHAIN_GAP)
                for a in chained(fa.fetch(c, lo, hi).upper(), unit):
                    if a[2] >= MIN_ARRAY and a[2] >= MIN_DENSITY * (a[1] - a[0]):
                        out.append((c, lo + a[0], lo + a[1], a[2]))
    return out


def _find_all(seq: str, x: str):
    i = seq.find(x)
    while i >= 0:
        yield i
        i = seq.find(x, i + 1)


def detect(fasta, chroms: dict) -> dict:
    """{'unit': canonical unit or None, 'enrichment': fold, 'ends': [(chrom, side, start, end, unit_bp,
    gap_to_end)], 'interstitial': [(chrom, start, end, unit_bp)]} (0-based)."""
    with pysam.FastaFile(str(fasta)) as fa:
        term = Counter()
        for _, _, s, _ in _ends(fa, chroms, TERM_SCAN):
            term.update(tandem_units(s))
        if not term:
            return dict(unit=None, enrichment=float("nan"), ends=[], interstitial=[])
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
                its = interstitial(fasta, chroms, unit)
                log.info(f"telomeres: unit {unit} ({enrich:.0f}-fold denser at the ends); arrays at {len(rows)}/"
                         f"{n_term} chromosome ends; {len(its)} interstitial arrays")
                return dict(unit=unit, enrichment=round(enrich, 1), ends=rows, interstitial=its)
    log.info("telomeres: no 5-12-bp unit is enriched in tandem arrays at the chromosome ends")
    return dict(unit=None, enrichment=float("nan"), ends=[], interstitial=[])
