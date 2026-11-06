"""rDNA arrays: 45S (nucleolus organizers, NORs) found with the universal 18S rRNA motif, and 5S found with the Rfam
5S rRNA model.

Every 45S unit holds one 18S gene, and the 18S "universal" region (GTACACACCGCCCGTC, E. coli 16S positions
1390-1405) is identical across eukaryotes. A NOR is a tandem run of 45S units, so it shows as a run of motif hits at
the unit spacing (CHM13 chr15: 44.3 kb; rice chr9: 8.4 kb). Single hits are pseudogenes and are ignored. Organelle
SSU rRNA (plastid and bacterial 16S, mitochondrial 12S and 18S) holds the motif too, followed by A where nuclear 18S has
G, and its next 19 bases differ from nuclear 18S at 11-17 positions, so such a hit is not counted: insertions of
organelle DNA otherwise form arrays (rice chr4 9.23-9.37 Mb: 4 hits, 99.8% identical to the plastid genome; 145 'arrays'
on one Salicornia chromosome). The A alone would not do: nuclear copies with that one base mutated keep the next 19
bases (at most 5 differ; barley 7H, with 10 5.8S genes in 120 kb). Strand tells nothing: large real NORs hold both
orientations (Coffea arabica: 658 units, 336 and 322). An exact scan costs seconds per Gb, where an HMM search of the
18S model costs ~40 min per 3 Gb and stalls on low-complexity sequence; the satellite engine's own rDNA tags see only
spacer repeats next to satellite arrays.

5S genes (120 bp, in their own tandem arrays) are less conserved than the 18S core, so they are found with the profile
HMM that Rfam embeds in its 5S covariance model (RF00001, in sat/data/rfam_rrna.cm), searched like the satellite
engine's monomers: nhmmer per chunk on a thread pool, E-values scaled to the chromosome. A 5S array is a tandem run of
whole genes; pseudogenes and 5S-derived fragments are mostly partial (CHM13: all 126 hits of its 1q42 array whole, 108
of 434 scattered hits). The satellite engine sees 5S arrays only as one of its satellite families, with units up to
2 kb: the 2.2-kb units of mammals and birds are beyond that (none of 14 T2T mammals had a 5S family)."""
from __future__ import annotations

import io
import re
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pysam

from .util import log

MOTIF = "GTACACACCGCCCGTC"       # 18S universal region (one per 45S unit), both strands searched
ORGANELLE_NEXT = "A"             # the base after the motif in organelle SSU rRNA (nuclear 18S: G; other bases say nothing)
NUCLEAR_AFTER = "CTCCTACCGATTGGATGGT"   # the 19 bases after that one in nuclear 18S (human, plants, yeast differ at 1-2)
UNIT_MAX = 100_000              # prior: 45S units are ~8-45 kb (plants, fungi ~8-11; mammals 43-45); hits farther apart
                                # than twice the longest unit are separate loci
MIN_UNITS = 3                   # a tandem array has >= 3 copies (the satellite engine's min_array_monomers)
CHUNK = 10_000_000              # bp fetched at once (bounds memory on Gb chromosomes)
FIVE_S_E = 1e-5                 # per-gene E-value scaled to the chromosome (the satellite engine's nhmmer_evalue)
FIVE_S_UNIT_MAX = 20_000        # prior: 5S units are 0.2-2.5 kb (plants 0.3-1 kb; mammals and birds 2.2 kb), and a 5S gene
                                # sits in each 9.1-kb rDNA unit of budding yeast; genes farther apart than twice that are
                                # separate loci
_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_RC = MOTIF.translate(_COMPLEMENT)[::-1]
_PAT = re.compile(f"{MOTIF}|{_RC}")


def organelle_type(seq: str, i: int, forward: bool) -> bool:
    """Whether the motif at seq[i:], read along its gene, is an organelle SSU rRNA copy: followed by ORGANELLE_NEXT, and
    the bases after that differ from NUCLEAR_AFTER at most positions (majority rule; 1,100 organelle hits differ at 11-17
    of 19, 110 nuclear hits with the first base mutated at 1-5, 4 hits between, over 143 genomes). On the reverse strand
    those bases lie before the hit, reverse complemented. Less than half the context on the sequence: not judged."""
    k, n = len(MOTIF), 1 + len(NUCLEAR_AFTER)
    ctx = seq[i + k:i + k + n] if forward else seq[max(0, i - n):i].translate(_COMPLEMENT)[::-1]
    rest = ctx[1:]
    if ctx[:1] != ORGANELLE_NEXT or 2 * len(rest) < len(NUCLEAR_AFTER):
        return False
    return 2 * sum(a != b for a, b in zip(rest, NUCLEAR_AFTER)) > len(rest)


def hits(fasta, chrom: str, length: int) -> np.ndarray:
    """Start positions of the motif (either strand) on one chromosome, organelle copies left out. Each chunk is read with
    the length of the context on each side, so the bases after a hit are seen across a chunk edge."""
    out, k, n = [], len(MOTIF), 1 + len(NUCLEAR_AFTER)
    with pysam.FastaFile(str(fasta)) as fa:
        for s in range(0, length, CHUNK):
            lo = max(0, s - n)
            seq = fa.fetch(chrom, lo, min(length, s + CHUNK + k + n)).upper()
            for m in _PAT.finditer(seq):
                pos = lo + m.start()
                if s <= pos < s + CHUNK and not organelle_type(seq, m.start(), m.group() == MOTIF):
                    out.append(pos)
    return np.array(sorted(out), dtype=np.int64)


def arrays(pos: np.ndarray, length: int, unit_max: int = UNIT_MAX) -> list:
    """[(start, end, units, unit_bp)]: runs of >= MIN_UNITS hits spaced <= unit_max apart, padded by one unit."""
    out = []
    if len(pos) < MIN_UNITS:
        return out
    cuts = np.flatnonzero(np.diff(pos) > unit_max) + 1
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


_HMM = {}


def five_s_hmm():
    """The Rfam 5S rRNA profile HMM (the filter HMM inside the bundled covariance model), read as DNA."""
    if "5S" not in _HMM:
        from pyhmmer.plan7 import HMMFile
        cm = (Path(__file__).parent / "sat" / "data" / "rfam_rrna.cm").read_text()
        i = cm.index("HMMER3/f", cm.index("NAME     5S_rRNA"))
        block = cm[i:cm.index("//", i) + 2].replace("ALPH  RNA", "ALPH  DNA") + "\n"
        with HMMFile(io.BytesIO(block.encode())) as f:
            _HMM["5S"] = f.read()
    return _HMM["5S"]


def genes_5s(fasta, chroms: dict, threads: int) -> dict:
    """{chrom: sorted 0-based starts of whole 5S genes} (either strand; E <= FIVE_S_E for the chromosome; >=
    FULL_MONOMER_COV of the model). Chunks of ANNOT_CHUNK overlap by two genes, so a gene cut by one chunk edge is whole
    in the next; duplicates in overlaps are merged."""
    from pyhmmer.easel import DigitalSequenceBlock, TextSequence

    from .sat.annotate import ABC, release_pipeline, search_one
    from .sat.params import ANNOT_CHUNK, FULL_MONOMER_COV
    hmm = five_s_hmm()
    jobs = [(c, s, min(L, s + ANNOT_CHUNK + 2 * hmm.M)) for c, L in chroms.items() for s in range(0, L, ANNOT_CHUNK)]
    local = threading.local()

    def work(job):
        c, s, e = job
        if not hasattr(local, "fa"):
            local.fa = pysam.FastaFile(str(fasta))
        blk = DigitalSequenceBlock(ABC, [TextSequence(name=c, sequence=local.fa.fetch(c, s, e).upper()).digitize(ABC)])
        scale, out = chroms[c] / (e - s), []
        for h in search_one(hmm, blk):
            a = h.best_domain.alignment
            if h.evalue * scale <= FIVE_S_E and a.hmm_to - a.hmm_from + 1 >= FULL_MONOMER_COV * hmm.M:
                out.append(s + min(a.target_from, a.target_to) - 1)
        return c, out

    pos = {}
    with ThreadPoolExecutor(max(1, threads)) as ex:
        for c, p in ex.map(work, jobs):
            pos.setdefault(c, set()).update(p)
    release_pipeline()
    return {c: np.array(sorted(v), dtype=np.int64) for c, v in pos.items() if v}


def detect_5s(fasta, chroms: dict, threads: int = 1) -> list:
    """[(chrom, start, end, units, unit_bp)] for every 5S rDNA array of the genome: tandem runs of >= MIN_UNITS whole
    genes spaced <= FIVE_S_UNIT_MAX apart."""
    pos = genes_5s(fasta, chroms, threads)
    out = []
    for c, L in chroms.items():
        out.extend((c, *a) for a in arrays(pos.get(c, np.zeros(0, np.int64)), L, FIVE_S_UNIT_MAX))
    if out:
        log.info(f"rDNA: {len(out)} 5S arrays on {len({x[0] for x in out})} chromosomes: "
                 + ", ".join(f"{c}:{s / 1e6:.2f}-{e / 1e6:.2f} Mb ({n} copies, {u}-bp unit)"
                             for c, s, e, n, u in sorted(out, key=lambda x: -x[3])[:6])
                 + (" ..." if len(out) > 6 else ""))
    else:
        log.info("rDNA: no 5S array in the assembly")
    return out
