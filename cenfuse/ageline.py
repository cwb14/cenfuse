"""Insertion age as a tie-breaker between two sites (docs/future.md #31).

Centromeric retrotransposons insert into active centromeric chromatin, so the youngest copies of a genome's centromeric LTR
families mark where each centromere is now, and older copies where it was (maize: recently active CRM subgroups colonize the
present centromeres, older ones their flanks; Wolfgruber et al. 2009, PLoS Genet 5:e1000743). A centromere that moved leaves its
satellite and its accumulated LTR copies behind; the two sequence lines then agree on the old site, and only validated Hi-C sees
the new one (P. annua Pa1A: Hi-C 21 Mb from the satellite; CRM copies at the Hi-C site have a median age of about 0 years, at the
satellite site 325,000; one family alone 0 vs 157,000, P 2e-7).

As a locator on its own, the core of the young copies does worse than the density core (on CENH3 for 128 vs 155 of 172
chromosomes with truth), so age never places or moves a call by itself. It only breaks the tie the gate cannot: a strong peak of
validated Hi-C far from a `high` sequence call moves it only with a second line at the peak (fuse), and copies of the chromosome's
signature families that are younger at the peak than at the call are one (flag moved_with_age).

Like every line, it must first show in the genome at hand that its premise holds (premise): on the OTHER chromosomes, the younger
half of the signature's dated copies lies in the cores more often than the older half. That holds for CRM and Bryco, whose young
copies concentrate in the cores, and not for lineages that accumulate near centromeres with age (Athila, Retand, most Tekay:
docs/design.md, the age of the elements); there age says nothing about where a centromere is now.

Ages are LTRquest's K2P_T (insertion time from the divergence of the two LTRs); a file without them gives no test."""
from __future__ import annotations

import numpy as np
from scipy.stats import fisher_exact, mannwhitneyu

from .params import HIC_FAR_BINS, HIC_TOL_BINS
from .sat import signature as sigmod
from .sat.params import TE_CONTRA_ALPHA, TE_MIN_SITES

AGE_ALPHA = TE_CONTRA_ALPHA     # the LTR line's level for overruling a call (0.01, conventional), one-sided
AGE_MIN = TE_MIN_SITES          # dated copies needed at each site, as the three insertions that make an LTR core (copies, not sites)


def load(gff, chroms: dict):
    """The element table of a tes/ file (signature.load_table), or None: no families, or no element with an age."""
    if gff is None:
        return None
    t = sigmod.load_table(gff, chroms)
    if t is None or not np.isfinite(t.age).any():
        return None
    return t


def family_mask(t, sig, chrom: str) -> np.ndarray:
    """Elements of chromosome chrom with an age that belong to its signature: the families learned for it from the other
    chromosomes (signature.build, leave one chromosome out), else the seed (a family, or with a lineage-level seed its lineage)."""
    names = {c: i for i, c in enumerate(t.names)}
    if chrom not in names:
        return np.zeros(len(t.mid), dtype=bool)
    on = (t.chrom == names[chrom]) & np.isfinite(t.age)
    core = sig.cores.get(chrom) or {}
    fams = core.get("families") or [sig.seed]
    if sig.level == "lineage" and list(fams) == [sig.seed]:
        if sig.lineage not in t.lineages:
            return np.zeros(len(t.mid), dtype=bool)
        return on & (t.lineage == t.lineages.index(sig.lineage))
    codes = [t.families.index(f) for f in fams if f in t.families]
    return on & np.isin(t.family, codes)


def compare(t, mask: np.ndarray, peak: tuple, call: tuple, margin: int) -> dict:
    """Ages of the masked elements whose midpoints lie at the peak and at the call, each widened by margin bp on both sides
    (elements at both are dropped). p: one-sided Mann-Whitney, the copies at the peak younger; NaN with fewer than AGE_MIN
    dated copies at either site."""
    def zone(iv):
        return (t.mid >= iv[0] - margin) & (t.mid < iv[1] + margin)
    a, b = zone(peak), zone(call)
    x, y = t.age[mask & a & ~b], t.age[mask & b & ~a]
    out = dict(n_peak=int(len(x)), n_call=int(len(y)), med_peak=float(np.median(x)) if len(x) else float("nan"),
               med_call=float(np.median(y)) if len(y) else float("nan"), p=float("nan"))
    if len(x) >= AGE_MIN and len(y) >= AGE_MIN:
        out["p"] = float(mannwhitneyu(x, y, alternative="less").pvalue)
    return out


def premise(t, sig, chrom: str) -> float:
    """One-sided P (Fisher's exact test) that, on the chromosomes other than chrom with a significant signature core, the younger
    half of the signature's dated copies (split at their pooled median age; copies at the median left out) lies in the cores more
    often than the older half: young copies mark where this genome's centromeres are now. NaN without dated copies on both sides."""
    names = {c: i for i, c in enumerate(t.names)}
    keep = np.zeros(len(t.mid), dtype=bool)
    inside = np.zeros(len(t.mid), dtype=bool)
    for c, core in sig.cores.items():
        if c == chrom or c not in names or not core.get("significant"):
            continue
        m = family_mask(t, sig, c)
        keep |= m
        inside |= m & (t.mid >= int(core["start"])) & (t.mid < int(core["end"]))
    if keep.sum() < 2 * AGE_MIN:
        return float("nan")
    age = t.age[keep]
    med = np.median(age)
    young, old, ins = age < med, age > med, inside[keep]
    table = [[int((young & ins).sum()), int((young & ~ins).sum())], [int((old & ins).sum()), int((old & ~ins).sum())]]
    if not (young.any() and old.any()) or not (ins.any() and (~ins).any()):
        return float("nan")
    return float(fisher_exact(table, alternative="greater")[1])


def for_genome(gff, chroms: dict, sig, hic_table, calls) -> dict:
    """{chrom: compare(...)} for the chromosomes where age can break a tie: a `high` sequence call, and a strong peak of validated
    Hi-C (gate passed) whose core lies more than a bin from it (fuse decides whether it is far enough to oppose); {} without such a
    chromosome, a signature, Hi-C or dated elements. The margin is the Hi-C line's resolution (HIC_FAR_BINS bins), on both sites
    alike."""
    if sig is None or hic_table is None or not len(hic_table):
        return {}
    rows = {r["chrom"]: r for r in hic_table.to_dict("records")}
    need = []
    for r in calls.itertuples(index=False):
        h = rows.get(r.chrom)
        if h is None or int(r.start) < 0 or r.confidence != "high" or int(h.get("anch_core_start", -1)) < 0 \
                or not h.get("gate_pass") or not h.get("anch_strong"):
            continue
        core, call, b = (int(h["anch_core_start"]), int(h["anch_core_end"])), (int(r.start), int(r.end)), int(h["bin_size"])
        if max(0, max(core[0], call[0]) - min(core[1], call[1])) > HIC_TOL_BINS * b:
            need.append((r.chrom, core, call, HIC_FAR_BINS * b))
    if not need:
        return {}
    t = load(gff, chroms)
    if t is None:
        return {}
    return {c: dict(compare(t, family_mask(t, sig, c), pk, cl, m), premise=premise(t, sig, c)) for c, pk, cl, m in need}
