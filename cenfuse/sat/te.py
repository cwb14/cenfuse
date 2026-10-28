"""Retrotransposon evidence: the LTR-RT clade that concentrates at one locus per chromosome marks the centromere
(Athila in A. thaliana, CRM in grasses). Optional input GENOME_DIR/tes/<genome>.gff3[.gz]; design:
docs/2026-10-05-retrotransposon-evidence.md."""
import gzip
import math
from pathlib import Path

import numpy as np

from .params import TE_ALPHA, TE_BINS, TE_CENTROMERIC_CLADE, TE_MIN_BIN, TE_MIN_ELEMENTS, TE_SMOOTH
from .util import log

IGNORED_CLADES = {"unknown", "mixture", "na", "none", "."}
# DANTE_LTR (REXdb) classification ranks above the lineage: an element classified only this far has no clade
DANTE_RANKS = {"class_i", "ltr", "ty1/copia", "ty3/gypsy", "chromovirus", "non-chromovirus", "ota", "tat"}


def dante_lineage(path: str) -> tuple:
    """(clade, superfamily) of a DANTE_LTR Final_Classification (Class_I|LTR|Ty3/gypsy|chromovirus|CRM -> CRM, Gypsy):
    the lineage is the last field unless that is a rank above it; superfamily Copia or Gypsy from Ty1/copia, Ty3/gypsy."""
    f = path.split("|")
    sf = next((x.split("/")[1].capitalize() for x in f if x.lower() in ("ty1/copia", "ty3/gypsy")), "")
    return ("" if f[-1].lower() in DANTE_RANKS else f[-1]), sf


def _clade(ftype: str, attrs: dict) -> str:
    """The clade of one GFF3 feature that is a classified LTR retrotransposon, else ''. Two formats: a type holding
    'LTR_retrotransposon' with clade= (LTRquest; EDTA + TEsorter), or a DANTE_LTR 'transposable_element' whose
    Final_Classification is an LTR path (dante_lineage) and whose Rank holds L (LTRs found: a rank-D domain cluster is
    a fragment, not an element; potato: 195 of DANTE_LTR's 223 CRM, with cores up to 15 Mb off CENH3). Keys are
    lower-cased."""
    if "LTR_retrotransposon" in ftype:
        clade = attrs.get("clade", "")
    elif ftype == "transposable_element" and "|LTR|" in attrs.get("final_classification", ""):
        if "L" not in attrs.get("rank", "DL"):          # DANTE_LTR rank D: protein domains, no LTRs found
            return ""
        clade = dante_lineage(attrs["final_classification"])[0]
    else:
        return ""
    return "" if clade.lower() in IGNORED_CLADES else clade


def _te_rows(gff, chroms: dict):
    """(chrom, start, end, clade) of every classified LTR retrotransposon on the genome's sequences (0-based)."""
    with (gzip.open(gff, "rt") if str(gff).endswith(".gz") else open(gff)) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\r\n").split("\t")
            if len(p) < 9 or p[0] not in chroms or not ("LTR_retrotransposon" in p[2] or p[2] == "transposable_element"):
                continue
            clade = _clade(p[2], {k.lower(): v for k, _, v in (kv.partition("=") for kv in p[8].split(";")) if v})
            if clade:
                yield p[0], int(p[3]) - 1, int(p[4]), clade


def te_census(gff, chroms: dict) -> dict:
    """Counts of a tes/ file's LTR retrotransposon features (types holding 'LTR_retrotransposon', DANTE_LTR elements of
    an LTR path): on the genome's sequences (ltr), those with a lineage there (classified, as _clade reads them), and
    on other sequence names (foreign). Tells an unclassified file (EDTA without TEsorter: superfamily only) or one for
    another naming from one whose clades are merely not concentrated."""
    n = dict(ltr=0, classified=0, foreign=0)
    with (gzip.open(gff, "rt") if str(gff).endswith(".gz") else open(gff)) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\r\n").split("\t")
            if len(p) < 9:
                continue
            attrs = {k.lower(): v for k, _, v in (kv.partition("=") for kv in p[8].split(";")) if v}
            if not ("LTR_retrotransposon" in p[2] or (p[2] == "transposable_element"
                                                      and "|LTR|" in attrs.get("final_classification", ""))):
                continue
            if p[0] not in chroms:
                n["foreign"] += 1
                continue
            n["ltr"] += 1
            n["classified"] += bool(_clade(p[2], attrs))
    return n


TE_EXT = (".gff3", ".gff3.gz", ".gff", ".gff.gz")


def te_gff(complex_dir, genome: str):
    for ext in TE_EXT:
        f = Path(complex_dir) / "tes" / f"{genome}{ext}"
        if f.exists():
            return f
    return None


def load_te_gff(gff, chroms: dict) -> dict:
    """{clade: {chrom: sorted midpoints}} of classified LTR retrotransposons (_clade) on the genome's sequences."""
    out = {}
    for c, s, e, clade in _te_rows(gff, chroms):
        out.setdefault(clade, {}).setdefault(c, []).append((s + e) // 2)
    return {c: {ch: np.sort(np.array(v, dtype=np.int64)) for ch, v in d.items()} for c, d in out.items()}


def load_te_intervals(gff, chroms: dict) -> dict:
    """{chrom: (starts, ends, clades)} of clade-classified LTR retrotransposons, sorted by start, 0-based half-open.
    Elements without a clade are left out: an array of kb-long tandem units looks like an LTR pair to structural
    detectors, so unclassified calls on tandem arrays are weak evidence of a retrotransposon (Carassius 1.6-2.0-kb
    families: 50-85% inside 'unknown' elements)."""
    out = {}
    for c, s, e, clade in _te_rows(gff, chroms):
        out.setdefault(c, []).append((s, e, clade))
    res = {}
    for c, v in out.items():
        v.sort()
        res[c] = (np.array([x[0] for x in v], dtype=np.int64), np.array([x[1] for x in v], dtype=np.int64),
                  [x[2] for x in v])
    return res


def ltr_overlap(arrays, intervals: dict, longest: int) -> dict:
    """{family: (share of array bp inside LTR retrotransposons, clade holding most of that overlap)}. A tandem repeat
    that lies inside retrotransposons is part of them (Musa: a 62-bp repeat inside SIRE elements), or derived from
    them (Coffea CRM / Tekay kb units); reported, not judged. longest: the longest element, bounds the search."""
    out = {}
    for fam, d in arrays.groupby("family"):
        bp = inside = 0
        clade_bp = {}
        for c, x in d.groupby("chrom"):
            if c not in intervals:
                bp += int((x.end - x.start).sum())
                continue
            st, en, cl = intervals[c]
            for s, e in zip(x.start.to_numpy(), x.end.to_numpy()):
                bp += int(e - s)
                lo, hi = np.searchsorted(st, s - longest), np.searchsorted(st, e)
                if lo >= hi:
                    continue
                os_ = np.maximum(0, np.minimum(en[lo:hi], e) - np.maximum(st[lo:hi], s))
                if not os_.any():
                    continue
                cov = np.zeros(int(e - s), dtype=bool)                  # union of overlapping elements
                for a, b in zip(np.maximum(st[lo:hi], s) - s, np.minimum(en[lo:hi], e) - s):
                    if b > a:
                        cov[a:b] = True
                inside += int(cov.sum())
                k = lo + int(np.argmax(os_))
                clade_bp[cl[k]] = clade_bp.get(cl[k], 0) + int(os_.max())
        out[fam] = (inside / bp if bp else 0.0, max(clade_bp, key=clade_bp.get) if clade_bp else ".")
    return out


def binom_sf(k: int, n: int, q: float) -> float:
    """P(X >= k), X ~ Binomial(n, q); exact sum in log space."""
    if k <= 0:
        return 1.0
    if q <= 0:
        return 0.0
    if q >= 1:
        return 1.0
    lq, l1q = math.log(q), math.log1p(-q)
    terms = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * lq + (n - i) * l1q
             for i in range(k, n + 1)]
    m = max(terms)
    return min(1.0, math.exp(m) * sum(math.exp(t - m) for t in terms))


def core(mids: np.ndarray, L: int) -> dict:
    """Peak region of one clade on one chromosome: bins of L/TE_BINS (>= TE_MIN_BIN), smoothed over TE_SMOOTH bins;
    core = bins around the maximum at >= half the peak (full width at half maximum). Significant if the binomial
    tail of the core count, Bonferroni-corrected for L/core-length positions, is below TE_ALPHA."""
    n = len(mids)
    if n < TE_MIN_ELEMENTS:
        return None
    b = max(TE_MIN_BIN, L // TE_BINS)
    nb = L // b + 1
    h = np.bincount(np.minimum(mids // b, nb - 1), minlength=nb).astype(float)
    sm = np.convolve(h, np.ones(TE_SMOOTH), mode="same")
    pk = int(np.argmax(sm))
    half = sm[pk] / 2
    lo = pk
    while lo > 0 and sm[lo - 1] >= half:
        lo -= 1
    hi = pk
    while hi < nb - 1 and sm[hi + 1] >= half:
        hi += 1
    s, e = lo * b, min(L, (hi + 1) * b)
    k = int(np.searchsorted(mids, e) - np.searchsorted(mids, s))
    q = (e - s) / L
    pval = min(1.0, binom_sf(k, n, q) * (L / (e - s)))
    return dict(start=s, end=e, k=k, n=n, share=k / n, p=pval, significant=pval < TE_ALPHA)


def candidate_clades(tes: dict, chroms: dict) -> list:
    """Clades that concentrate at one locus per chromosome: a significant core on the majority of chromosomes.
    Returns [(clade, {chrom: significant core})], most concentrated (median core share) first. Several clades can
    qualify (pericentromeric ones too), and the most concentrated is not always the centromeric one (B. juncea A),
    so callers keep them all."""
    out = []
    for clade, d in tes.items():
        cores = {c: core(d[c], L) for c, L in chroms.items() if c in d}
        sig = {c: x for c, x in cores.items() if x and x["significant"]}
        if len(sig) * 2 > len(chroms):
            out.append((float(np.median([x["share"] for x in sig.values()])), clade, sig))
    out.sort(key=lambda t: (-t[0], t[1]))
    for med, clade, sig in out:
        log.debug(f"te: concentrated clade {clade} (median core share {med:.2f}, {len(sig)} sequences)")
    return [(clade, sig) for _, clade, sig in out]


def marker(cands: list):
    """(clade, {chrom: core}, decisive) of the centromere marker among the candidate clades (most concentrated first),
    or None. A qualifying CRM clade is the marker and its evidence is decisive (it may reject a satellite and sets
    confidence); otherwise the most concentrated clade is the marker and only flags disagreement (B. juncea A: only
    the pericentromeric Retand/Athila qualify, their cores sit beside CentBr)."""
    for clade, cores in cands:
        if clade.lower() == TE_CENTROMERIC_CLADE.lower():
            return clade, cores, True
    return (cands[0][0], cands[0][1], False) if cands else None


def choose_marker(tes: dict, chroms: dict, cands: list):
    """The marker among the candidate clades (marker), with CRM cores taking precedence chromosome by chromosome: when the
    decisive marker is not CRM, a significant core of a CRM clade replaces the marker's core on its chromosome (core
    'clade' = CRM; flags only). CRM is the prior centromere-specific lineage, so its own significant concentration on a
    chromosome outranks another clade's there (DANTE_LTR rice: CRM cores on 5/12 chromosomes, all on CENH3, short of the
    majority; the marker Tekay's cores on 11/12, all 0.2-6.2 Mb away). Without any qualifying clade, the CRM cores alone
    are the (flags-only) marker."""
    m = marker(cands)
    if m is not None and m[2]:
        return m
    for clade, d in tes.items():
        if clade.lower() != TE_CENTROMERIC_CLADE.lower():
            continue
        cores = {c: core(d[c], L) for c, L in chroms.items() if c in d}
        sig = {c: dict(x, clade=clade) for c, x in cores.items() if x and x["significant"]}
        if sig:
            log.debug(f"te: {clade} significant on {len(sig)}/{len(chroms)} sequences: its cores there")
            return (m[0], {**m[1], **sig}, False) if m is not None else (clade, sig, False)
    return m


def centromeric_clade(tes: dict, chroms: dict):
    """(clade, {chrom: significant core}) of the most concentrated candidate clade, or (None, {})."""
    c = candidate_clades(tes, chroms)
    return c[0] if c else (None, {})


def gap(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, max(a0, b0) - min(a1, b1))


def agrees(call_start: int, call_end: int, c: dict) -> bool:
    """Two estimates of one locus agree if the gap between them is at most the wider of the two (centromeric
    retrotransposons are partly pericentromeric, so a core can sit beside the satellite array)."""
    return gap(call_start, call_end, c["start"], c["end"]) <= max(c["end"] - c["start"], call_end - call_start)
