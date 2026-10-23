"""45S rDNA arrays (nucleolus organizers, NORs), found with the universal 18S rRNA motif.

Every 45S unit holds one 18S gene, and the 18S "universal" region (GTACACACCGCCCGTC, E. coli 16S positions
1390-1405) is identical across eukaryotes. A NOR is a tandem run of 45S units, so it shows as a run of motif hits at
the unit spacing (CHM13 chr15: 44.3 kb; rice chr9: 8.4 kb). Single hits are pseudogenes or organelle insertions
(mitochondrial 12S holds the motif too) and are ignored. An exact scan costs seconds per Gb, where an HMM search of
the 18S model costs ~40 min per 3 Gb and stalls on low-complexity sequence; the satellite engine's own rDNA tags
see only spacer repeats next to satellite arrays."""
from __future__ import annotations

import re

import numpy as np
import pysam

from .util import log

MOTIF = "GTACACACCGCCCGTC"       # 18S universal region (one per 45S unit), both strands searched
UNIT_MAX = 100_000              # prior: 45S units are ~8-45 kb (plants, fungi ~8-11; mammals 43-45); hits farther apart
                                # than twice the longest unit are separate loci
MIN_UNITS = 3                   # a tandem array has >= 3 copies (the satellite engine's min_array_monomers)
CHUNK = 10_000_000              # bp fetched at once (bounds memory on Gb chromosomes)
_RC = MOTIF.translate(str.maketrans("ACGT", "TGCA"))[::-1]
_PAT = re.compile(f"{MOTIF}|{_RC}")


def hits(fasta, chrom: str, length: int) -> np.ndarray:
    """Start positions of the motif (either strand) on one chromosome."""
    out, k = [], len(MOTIF)
    with pysam.FastaFile(str(fasta)) as fa:
        for s in range(0, length, CHUNK):
            seq = fa.fetch(chrom, s, min(length, s + CHUNK + k - 1)).upper()
            out.extend(s + m.start() for m in _PAT.finditer(seq))
    return np.unique(np.array(out, dtype=np.int64))


def arrays(pos: np.ndarray, length: int) -> list:
    """[(start, end, units, unit_bp)]: runs of >= MIN_UNITS hits spaced <= UNIT_MAX apart, padded by one unit."""
    out = []
    if len(pos) < MIN_UNITS:
        return out
    cuts = np.flatnonzero(np.diff(pos) > UNIT_MAX) + 1
    for run in np.split(pos, cuts):
        if len(run) >= MIN_UNITS:
            unit = int(np.median(np.diff(run)))
            out.append((max(0, int(run[0]) - unit), min(length, int(run[-1]) + unit), len(run), unit))
    return out


def detect(fasta, chroms: dict) -> list:
    """[(chrom, start, end, units, unit_bp)] for every NOR of the genome."""
    out = []
    for c, L in chroms.items():
        out.extend((c, *a) for a in arrays(hits(fasta, c, L), L))
    if out:
        log.info(f"rDNA: {len(out)} 45S arrays (NORs) on {len({x[0] for x in out})} chromosomes: "
                 + ", ".join(f"{c}:{s / 1e6:.1f}-{e / 1e6:.1f} Mb ({n} units, {u / 1e3:.1f} kb)"
                             for c, s, e, n, u in out[:6]) + (" ..." if len(out) > 6 else ""))
    else:
        log.info("rDNA: no 45S array (NOR) in the assembly")
    return out
