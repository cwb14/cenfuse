"""The LTR signature: which LTR retrotransposon families mark this genome's centromeres.

A lineage is a weak proxy for centromere targeting: only part of the CRM clade is centromeric (Neumann et al. 2011), and
other lineages do it too (Bryco in moss, Osser in Chlorella, a Tekay family in tomato). So the signature is found at two
levels, lineage and family, and then learned for each chromosome from the other chromosomes:

1. Candidates. A lineage or a family whose elements concentrate at one locus on most chromosomes (te.core, the same
   test as te.candidate_clades). A CRM family or the CRM lineage is a candidate from TE_SEED_MIN_CORES cores on, since CRM is
   the known centromeric lineage (every CRM core short of a majority lay on CENH3, docs/design.md).
2. Specificity of a candidate: the share of its elements inside its cores / the share of those chromosomes the cores
   cover. It is a fold enrichment, comparable between a family of 50 elements and a lineage of 50,000.
3. The seed: the most specific CRM candidate, else the most specific candidate. Only classified lineages seed: a family
   without one can be a tandem array read as LTR pairs (te.load_te_intervals), which would not be independent of the
   satellite line.
4. Recruitment, leave one chromosome out. For chromosome x, the families enriched in the current cores of the OTHER
   chromosomes (binomial against the share of those chromosomes the cores cover, Bonferroni over the families tested,
   P <= TE_ALPHA) make x's signature, and x's core is their densest locus (te.core). x never shapes its own signature, so
   a family recruited elsewhere that also concentrates on x is independent evidence there. The first round starts from
   the seed's cores, each next one from the significant cores of the round before, until they stop changing.

What a signature may do is decided in decisive(): a CRM seed decides as it always has; another lineage decides when it is
specific, unless a dominant satellite falls one per chromosome (cenfuse.sat.karyotype). A file without family= (not LTRquest) has
no signature: the lineage rules of te.py stand alone."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import binom, chi2, mannwhitneyu

from ..runtime import atomic_open
from .params import (TE_AGE_ALPHA, TE_AGE_MIN, TE_ALPHA, TE_BINS, TE_CENTROMERIC_CLADE, TE_MAJOR_SHARE, TE_MIN_BIN, TE_MIN_SITES,
                     TE_PEAK_MAX, TE_PEAK_SHARE, TE_POOLED_BREADTH, TE_POOLED_MIN_CORES, TE_RECRUIT_MIN, TE_RECRUIT_ROUNDS, TE_SEED_MIN_CORES,
                     TE_SMOOTH, TE_SPECIFIC_FOLD)
from .te import IGNORED_CLADES, Coords, Mids, binom_sf, core
from .util import log, open_text

NO_FAMILY = {"", ".", "na", "none", "unknown"}
FAMILY_MIN_SHARE = 0.5     # a file holds families when most of its LTR retrotransposons name one (majority rule)


@dataclass
class Table:
    """The LTR retrotransposons of one genome as columns (one row per element; -1 = none, NaN = unknown)."""
    names: list
    lengths: np.ndarray
    chrom: np.ndarray
    start: np.ndarray
    end: np.ndarray
    mid: np.ndarray
    lineage: np.ndarray
    family: np.ndarray
    age: np.ndarray
    lineages: list
    families: list
    family_lineage: np.ndarray     # family -> the most common lineage of its elements, -1 if none


@dataclass
class Signature:
    lineage: str                   # the seed's lineage
    level: str                     # 'lineage' or 'family'
    seed: str                      # the seed's name
    specificity: float             # the seed's fold enrichment in its cores
    seed_cores: dict               # {chrom: te.core} significant cores of the seed
    cores: dict                    # {chrom: te.core + source ('signature'/'seed') and families (names)}; significant or not
    candidates: list               # the candidates, most specific first: dicts without their cores
    pooled: bool = False           # seeded by a family's concentration pooled over the chromosomes (seed_of)
    age_core: float = float("nan")  # median K2P time of the signature's elements inside / outside its cores
    age_out: float = float("nan")
    young_old: float = float("nan")   # share of the younger half in the cores / share of the older half
    n_chroms: int = 0
    why: str = ""
    decisive: bool = False

    @property
    def crm(self) -> bool:
        return self.lineage.lower() == TE_CENTROMERIC_CLADE.lower()

    def significant(self) -> dict:
        """The significant cores; each knows whether it is strong (may replace a call: CRM's, or another lineage's holding
        most of the signature's elements on its chromosome)."""
        return {c: dict(x, strong=self.crm or x["share"] >= TE_MAJOR_SHARE) for c, x in self.cores.items() if x["significant"]}

    def marker(self):
        """The (clade, {chrom: significant core}, decisive) triple the satellite engine reads (te.marker)."""
        return self.lineage, self.significant(), self.decisive


def _attrs(col: str) -> dict:
    return {k.lower(): v for k, _, v in (kv.partition("=") for kv in col.split(";")) if v}


def load_table(gff, chroms: dict):
    """Table of the LTR_retrotransposon features on the genome's chromosomes, or None when most name no family (not an
    LTRquest annotation) or there are none."""
    names = list(chroms)
    idx = {c: i for i, c in enumerate(names)}
    rows, lin_code, fam_code = [], {}, {}
    xy = Coords(gff)
    with open_text(gff) as f:
        for n, line in enumerate(f, 1):
            if line.startswith("#"):
                continue
            p = line.rstrip("\r\n").split("\t")
            if len(p) < 9 or p[0] not in idx or "LTR_retrotransposon" not in p[2]:
                continue
            iv = xy(p, n)
            if iv is None:
                continue
            a = _attrs(p[8])
            cl, fa = a.get("clade", ""), a.get("family", "")
            lin = -1 if cl.lower() in IGNORED_CLADES or not cl else lin_code.setdefault(cl, len(lin_code))
            fam = -1 if fa.lower() in NO_FAMILY else fam_code.setdefault(fa, len(fam_code))
            try:
                age = float(a.get("k2p_t", "nan"))
            except ValueError:
                age = float("nan")
            rows.append((idx[p[0]], iv[0], iv[1], lin, fam, age))
    xy.done()
    if not rows:
        return None
    x = np.array(rows, dtype=np.float64)
    fam = x[:, 4].astype(np.int64)
    if (fam >= 0).mean() < FAMILY_MIN_SHARE:
        return None
    lin = x[:, 3].astype(np.int64)
    fl = np.full(len(fam_code), -1, dtype=np.int64)
    both = (fam >= 0) & (lin >= 0)
    if both.any():
        order = np.lexsort((lin[both], fam[both]))
        pairs, counts = np.unique(np.column_stack((fam[both][order], lin[both][order])), axis=0, return_counts=True)
        best = {}
        for (f, l), n in zip(pairs.tolist(), counts.tolist()):
            if n > best.get(f, (0, -1))[0] or (n == best.get(f, (0, -1))[0] and l < best[f][1]):
                best[f] = (n, l)
        for f, (_, l) in best.items():
            fl[f] = l
    st, en = x[:, 1].astype(np.int64), x[:, 2].astype(np.int64)
    return Table(names, np.array([chroms[c] for c in names], dtype=np.int64), x[:, 0].astype(np.int64), st, en,
                 (st + en) // 2, lin, fam, x[:, 5], list(lin_code), list(fam_code), fl)


def _sites(s: np.ndarray, e: np.ndarray) -> int:
    """Independent insertion sites: overlapping or abutting elements (nested insertions) count once."""
    if not len(s):
        return 0
    o = np.argsort(s, kind="stable")
    s, e = s[o], e[o]
    return int(1 + (s[1:] > np.maximum.accumulate(e)[:-1]).sum())


def _core_of(t: Table, sel: np.ndarray, k: int):
    """te.core of the elements sel (indices) on chromosome index k, or None."""
    x = sel[t.chrom[sel] == k]
    if not len(x):
        return None
    return core(Mids(np.sort(t.mid[x]), _sites(t.start[x], t.end[x])), int(t.lengths[k]))


def _cores(t: Table, sel: np.ndarray) -> dict:
    """{chromosome index: significant core} of the elements sel (indices)."""
    out = {}
    for k in np.unique(t.chrom[sel]):
        c = _core_of(t, sel, int(k))
        if c and c["significant"]:
            out[int(k)] = c
    return out


def _in_cores(t: Table, sel: np.ndarray, cores: dict) -> np.ndarray:
    """Which of the elements sel lie inside the core of their chromosome."""
    ins = np.zeros(len(sel), dtype=bool)
    ch = t.chrom[sel]
    for k, c in cores.items():
        on = ch == k
        ins[on] = (t.mid[sel[on]] >= c["start"]) & (t.mid[sel[on]] < c["end"])
    return ins


def specificity(t: Table, sel: np.ndarray, cores: dict) -> tuple:
    """(fold enrichment, share inside the cores, share of the chromosomes covered), over the chromosomes with a core."""
    if not cores:
        return 0.0, 0.0, 0.0
    on = sel[np.isin(t.chrom[sel], list(cores))]
    share = float(_in_cores(t, on, cores).mean()) if len(on) else 0.0
    width = sum(c["end"] - c["start"] for c in cores.values()) / float(t.lengths[list(cores)].sum())
    return (share / width if width else 0.0), share, width


def candidates(t: Table) -> list:
    """Lineages and families that qualify as seeds (module doc), most specific first (ties by name, for a reproducible
    choice). Each: level, name, lineage, n, cores {chromosome index: core}, specificity, share, width."""
    n_chr = len(t.names)
    most = n_chr // 2 + 1
    out = []

    def add(level, name, lineage, sel):
        crm = lineage.lower() == TE_CENTROMERIC_CLADE.lower()
        need = TE_SEED_MIN_CORES if crm else most
        if len(sel) < need * TE_MIN_SITES:
            return
        cores = _cores(t, sel)
        if len(cores) >= need:
            sp, share, width = specificity(t, sel, cores)
            out.append(dict(level=level, name=name, lineage=lineage, n=len(sel), cores=cores, specificity=sp,
                            share=share, width=width))

    for code, name in enumerate(t.lineages):
        add("lineage", name, name, np.flatnonzero(t.lineage == code))
    order = np.argsort(t.family, kind="stable")
    fam_sorted = t.family[order]
    for code, name in enumerate(t.families):
        lin = int(t.family_lineage[code])
        if lin < 0:
            continue
        lo, hi = np.searchsorted(fam_sorted, code), np.searchsorted(fam_sorted, code, side="right")
        add("family", name, t.lineages[lin], order[lo:hi])
    out.sort(key=lambda z: (-z["specificity"], z["name"]))
    return out


def pooled_candidates(t: Table) -> list:
    """Families (with a lineage) too thin to be significant on most chromosomes one by one, but significantly concentrated at one
    locus per chromosome taken together: more than half of the chromosomes show the family's core at the conventional level by
    themselves (TE_POOLED_BREADTH: a few strong chromosomes do not make it, as they could for Fisher's P alone), Fisher's combination
    of the core P of every chromosome holding one (te.core, TE_MIN_SITES sites) is at P <= TE_ALPHA after Bonferroni over the
    families tested, TE_POOLED_MIN_CORES cores are significant to recruit from, and the family is at least TE_SPECIFIC_FOLD enriched
    in its cores. A centromeric family of 120 elements on 12 chromosomes has about 10 copies per chromosome: too few for te.core to
    be significant on most (tomato LA1028), whereas a family's concentration across chromosomes has power to spare. The fold is an
    effect size, not a test: core windows fitted to few copies reach 10-fold by chance (uniform positions, 10 copies per
    chromosome, all reach it); the tests are the others (none of 2,500 such uniform families passed them).
    Each: as candidates(), plus n_cores_any (chromosomes with a core, significant or not)."""
    n_chr, n_fam = len(t.names), max(1, len(t.families))
    most = n_chr // 2 + 1
    order = np.argsort(t.family, kind="stable")
    fam_sorted = t.family[order]
    out = []
    for code, name in enumerate(t.families):
        lin = int(t.family_lineage[code])
        if lin < 0:
            continue
        lo, hi = np.searchsorted(fam_sorted, code), np.searchsorted(fam_sorted, code, side="right")
        sel = order[lo:hi]
        if len(sel) < most * TE_MIN_SITES:
            continue
        every = {}
        for k in np.unique(t.chrom[sel]):
            c = _core_of(t, sel, int(k))
            if c:
                every[int(k)] = c
        if 2 * sum(c["p"] <= TE_POOLED_BREADTH for c in every.values()) <= n_chr:
            continue
        stat = -2 * sum(np.log(max(c["p"], 1e-300)) for c in every.values())
        pval = float(chi2.sf(stat, 2 * len(every))) * n_fam
        sig = {k: c for k, c in every.items() if c["significant"]}
        sp, share, width = specificity(t, sel, every)
        if pval <= TE_ALPHA and sp >= TE_SPECIFIC_FOLD and len(sig) >= TE_POOLED_MIN_CORES:
            out.append(dict(level="family", name=name, lineage=t.lineages[lin], n=len(sel), cores=sig, specificity=sp,
                            share=share, width=width, pooled_p=pval, n_cores_any=len(every)))
    out.sort(key=lambda z: (-len(z["cores"]), -z["specificity"], z["name"]))
    return out


def seed_of(cands: list, n_chroms: int, pooled: list | None = None):
    """The most specific CRM candidate with cores on most chromosomes, else the most specific CRM candidate, else the most
    specific candidate, else None. A small CRM family on a few chromosomes can be more specific than the CRM lineage on
    most of them (rice: 26-fold on 3, 16-fold on 8), and the seed's cores stand wherever recruitment finds nothing.
    pooled (0.15): where no CRM candidate exists and the most specific candidate is under TE_SPECIFIC_FOLD (it could not decide) or
    there is none, the pooled family with the most significant cores (pooled_candidates; its Signature is marked pooled)."""
    crm = [z for z in cands if z["lineage"].lower() == TE_CENTROMERIC_CLADE.lower()]
    most = [z for z in crm if 2 * len(z["cores"]) > n_chroms]
    if most or crm:
        return (most or crm)[0]
    if pooled and (not cands or cands[0]["specificity"] < TE_SPECIFIC_FOLD):
        return dict(pooled[0], pooled=True)
    return cands[0] if cands else None


def recruit(t: Table, anchors: dict) -> dict:
    """{chromosome index: (core or None, recruited family codes)} for every chromosome (module doc, step 4). anchors:
    {chromosome index: core}."""
    nf, nk = len(t.families), len(t.names)
    has = t.family >= 0
    ins = np.zeros(len(t.mid), dtype=bool)
    ins[has] = _in_cores(t, np.flatnonzero(has), anchors)
    anch = np.zeros(nk, dtype=bool)
    anch[list(anchors)] = True
    on = has & anch[t.chrom]
    total = np.bincount(t.family[on] * nk + t.chrom[on], minlength=nf * nk).reshape(nf, nk)
    inside = np.bincount(t.family[on & ins] * nk + t.chrom[on & ins], minlength=nf * nk).reshape(nf, nk)
    tot_all, in_all = total.sum(axis=1), inside.sum(axis=1)
    width = np.zeros(nk)
    for k, c in anchors.items():
        width[k] = c["end"] - c["start"]
    length = np.where(anch, t.lengths, 0).astype(float)
    order = np.argsort(t.family, kind="stable")
    fam_sorted = t.family[order]
    out = {}
    for k in range(nk):
        lens = length.sum() - length[k]
        if anch.sum() - anch[k] < 2 or lens <= 0:
            out[k] = (None, np.zeros(0, dtype=np.int64))
            continue
        q = (width.sum() - width[k]) / lens
        tot, inn = tot_all - total[:, k], in_all - inside[:, k]
        test = np.flatnonzero(tot >= TE_RECRUIT_MIN)
        if not len(test):
            out[k] = (None, test)
            continue
        p = binom.sf(inn[test] - 1, tot[test], q) * len(test)
        fams = test[p <= TE_ALPHA]
        if not len(fams):
            out[k] = (None, fams)
            continue
        sel = np.concatenate([order[np.searchsorted(fam_sorted, f):np.searchsorted(fam_sorted, f, side="right")]
                              for f in fams])
        out[k] = (_core_of(t, sel, k), fams)
    return out


def _peaks(t: Table, sel: np.ndarray, k: int) -> list:
    """The density peaks of elements sel on chromosome index k, highest first: te.core's profile and half-height run, a
    peak's bins taken from those of the higher peaks (so two peaks never share a bin). Each: start, end, height, sites,
    ages (the finite K2P times inside it)."""
    L = int(t.lengths[k])
    b = max(TE_MIN_BIN, L // TE_BINS)
    nb = L // b + 1
    mids = t.mid[sel]
    h = np.bincount(np.minimum(mids // b, nb - 1), minlength=nb).astype(float)
    sm = np.convolve(h, np.ones(TE_SMOOTH), mode="same")
    used, out = np.zeros(nb, dtype=bool), []
    for i in np.argsort(-sm, kind="stable"):
        if sm[i] <= 0 or len(out) >= TE_PEAK_MAX:
            break
        if used[max(0, i - TE_SMOOTH // 2 - 1):min(nb, i + TE_SMOOTH // 2 + 2)].any():
            continue                                    # within reach of a higher peak: its shoulder
        half, lo, hi = sm[i] / 2, int(i), int(i)
        while lo > 0 and sm[lo - 1] >= half and not used[lo - 1]:
            lo -= 1
        while hi < nb - 1 and sm[hi + 1] >= half and not used[hi + 1]:
            hi += 1
        used[lo:hi + 1] = True
        s, e = lo * b, min(L, (hi + 1) * b)
        ins = (mids >= s) & (mids < e)
        out.append(dict(start=s, end=e, height=float(sm[i]), sites=_sites(t.start[sel][ins], t.end[sel][ins]),
                        ages=t.age[sel][ins][np.isfinite(t.age[sel][ins])]))
    return out


def _peak_core(t: Table, sel: np.ndarray, k: int, p: dict) -> dict:
    """The density peak p of elements sel on chromosome index k as a core: te.core's fields and te.core's test."""
    mids = np.sort(t.mid[sel])
    kin = int(np.searchsorted(mids, p["end"]) - np.searchsorted(mids, p["start"]))
    q = (p["end"] - p["start"]) / float(t.lengths[k])
    pval = min(1.0, binom_sf(kin, len(mids), q) / q)
    return dict(start=p["start"], end=p["end"], k=kin, n=len(mids), share=kin / len(mids), p=pval, significant=pval < TE_ALPHA, mids=mids)


def age_choice(t: Table, sel: np.ndarray, k: int) -> dict | None:
    """The peak to use instead of the highest, or None to keep the highest (te.core's own choice). Where the signature's copies
    form several comparable peaks (at least half the height of the highest, TE_PEAK_SHARE, with TE_MIN_SITES sites), the youngest
    significantly younger than the highest marks where the centromere is now: copies insert into the active kinetochore, so
    the older peak is where it was (young CRM subgroups colonize the present maize centromeres, older ones their flanks,
    Wolfgruber et al. 2009; P. annua Pa1A, future.md #31). A peak must itself be a significant core (te.core's test), or the move
    would swap a tested core for an untested one. One-sided Mann-Whitney P <= TE_AGE_ALPHA, Bonferroni over the peaks compared on
    this chromosome, TE_AGE_MIN dated copies at each peak, the same test as ageline's. Returns the chosen peak as a core (te.core's
    fields). Unlike ageline it needs no premise that young copies concentrate in the cores: within a genome whose young copies are
    spread out (tomato, young_old 1.0) this picks the younger of two clusters, and the density maximum is then the old cluster.
    Limits: recruited families are pooled, and families burst at different times, so a peak of a younger family reads as younger."""
    ps = _peaks(t, sel, k)
    if len(ps) < 2:
        return None
    top = ps[0]
    if len(top["ages"]) < TE_AGE_MIN:
        return None
    tested = []
    for p in ps[1:]:
        if p["height"] >= TE_PEAK_SHARE * top["height"] and p["sites"] >= TE_MIN_SITES and len(p["ages"]) >= TE_AGE_MIN:
            c = _peak_core(t, sel, k, p)
            if c["significant"]:
                tested.append((p, c))
    best, best_age = None, np.median(top["ages"])
    for p, c in tested:
        a = float(np.median(p["ages"]))
        if a < best_age and mannwhitneyu(p["ages"], top["ages"], alternative="less").pvalue * len(tested) <= TE_AGE_ALPHA:
            best, best_age = c, a
    return best


def second_peak(t: Table, sel: np.ndarray, k: int, pick: dict):
    """(start, end) of the highest density peak of elements sel on chromosome index k other than the one the core is, when at
    least half as high as it (TE_PEAK_SHARE) and with TE_MIN_SITES sites: reported only, a second place to look (a chromosome's
    signature elements can lie at the old and the present centromere), else None."""
    if not len(sel):
        return None
    ps = _peaks(t, sel, k)
    if not ps:
        return None
    here = max(ps, key=lambda p: min(p["end"], pick["end"]) - max(p["start"], pick["start"]))
    other = [p for p in ps if p is not here and p["height"] >= TE_PEAK_SHARE * here["height"] and p["sites"] >= TE_MIN_SITES]
    return (int(other[0]["start"]), int(other[0]["end"])) if other else None


def _age_columns(t: Table, cores: dict, fams: dict) -> tuple:
    """(median K2P time inside the cores, outside them, share of the younger half inside / share of the older half) of
    the signature's elements on the chromosomes with a core."""
    keep = np.zeros(len(t.mid), dtype=bool)
    for k in cores:
        keep |= (t.chrom == k) & np.isin(t.family, fams.get(k, ()))
    sel = np.flatnonzero(keep & np.isfinite(t.age))
    if len(sel) < 2:
        return float("nan"), float("nan"), float("nan")
    ins = _in_cores(t, sel, cores)
    age = t.age[sel]
    inner = float(np.median(age[ins])) if ins.any() else float("nan")
    outer = float(np.median(age[~ins])) if (~ins).any() else float("nan")
    med = np.median(age)
    y, o = age <= med, age > med
    if not y.any() or not o.any() or not ins[o].any():
        return inner, outer, float("nan")
    return inner, outer, float(ins[y].mean() / ins[o].mean())


def build(gff, chroms: dict):
    """The genome's Signature, or None (no family column, or no candidate). Its decisive field holds what is known
    without the satellite line: a CRM seed with significant cores on most chromosomes. decisive() completes it."""
    t = load_table(gff, chroms)
    if t is None:
        return None
    cands = candidates(t)
    seed = seed_of(cands, len(t.names), pooled_candidates(t))
    if seed is None:
        return None
    seed_elements = (np.flatnonzero(t.family == t.families.index(seed["name"])) if seed["level"] == "family"
                     else np.flatnonzero(t.lineage == t.lineages.index(seed["name"])))
    anchors, rounds, seen = seed["cores"], {}, None
    for _ in range(TE_RECRUIT_ROUNDS):
        rounds = recruit(t, anchors)
        sig = {k: c for k, (c, _) in rounds.items() if c and c["significant"]}
        now = sorted((k, c["start"], c["end"]) for k, c in sig.items())
        if now == seen:
            break
        seen = now
        if len(sig) >= TE_SEED_MIN_CORES:
            anchors = sig
    cores, fams = {}, {}
    for k, name in enumerate(t.names):
        c, f = rounds.get(k, (None, np.zeros(0, dtype=np.int64)))
        s = seed["cores"].get(k)
        if c is not None and c["significant"]:
            pick, src = c, "signature"
        elif s is not None:
            pick, src = s, "seed"
        elif c is not None:
            pick, src = c, "signature"
        else:
            continue
        # the copies this core is read from, whose ages may choose another of their density peaks
        sel = np.flatnonzero((t.chrom == k) & np.isin(t.family, f)) if src == "signature" and len(f) else seed_elements[t.chrom[seed_elements] == k]
        if len(sel) and (moved := age_choice(t, sel, k)) is not None:
            pick = dict(moved, age_from=(int(pick["start"]), int(pick["end"])))
        alt = second_peak(t, sel, k, pick)
        cores[name] = dict(pick, source=src, families=[t.families[i] for i in f] if src == "signature" else [seed["name"]],
                           fold=pick["share"] / ((pick["end"] - pick["start"]) / float(t.lengths[k])), alt=alt)
        fams[k] = f
    sig = {k: dict(cores[t.names[k]]) for k in range(len(t.names)) if t.names[k] in cores and cores[t.names[k]]["significant"]}
    ages = _age_columns(t, sig, fams) if sig else (float("nan"),) * 3
    s = Signature(lineage=seed["lineage"], level=seed["level"], seed=seed["name"], specificity=float(seed["specificity"]),
                  seed_cores={t.names[k]: c for k, c in seed["cores"].items()}, cores=cores, pooled=bool(seed.get("pooled")),
                  candidates=[{k: v for k, v in z.items() if k != "cores"} | {"n_cores": len(z["cores"])} for z in cands],
                  age_core=ages[0], age_out=ages[1], young_old=ages[2], n_chroms=len(t.names))
    s.decisive = bool(s.crm and 2 * len(s.significant()) > s.n_chroms)
    s.why = (f"{'CRM' if s.crm else s.lineage} seed {s.seed}, {s.specificity:.1f}-fold in its cores, significant signature "
             f"cores on {len(s.significant())}/{s.n_chroms} chromosomes" + ("; a family concentrated across chromosomes, too few "
             "copies per chromosome to be significant on most" if s.pooled else ""))
    log.debug(f"te: signature: {s.why}; {len(cands)} candidates")
    return s


def rejecting_cores(cores: dict) -> dict:
    """The cores specific enough to speak against a satellite family (classify.te_rejected): rejecting a family says more than
    replacing one call."""
    return {c: x for c, x in cores.items() if x.get("fold", TE_SPECIFIC_FOLD) >= TE_SPECIFIC_FOLD and x.get("strong", True)}


def decisive(s: Signature, dominant: bool = False, satellite_status: str = "none") -> Signature:
    """Complete s.decisive. CRM decides on most chromosomes (build). Another lineage needs a majority of significant cores (a pooled
    seed has shown its majority at TE_POOLED_BREADTH instead, pooled_candidates: it is the case where most chromosomes hold too few
    copies to be significant at TE_ALPHA) and
    TE_SPECIFIC_FOLD, and must not stand against a satellite line that has shown itself: a dominant satellite whose largest clusters
    fall one per chromosome (karyotype.satellite_test, status pass, or too few chromosomes to tell). Tomato's dominant satellite
    fails that test and its Tekay family decides; Arabidopsis' CEN180 passes and Athila only flags; moss has no dominant satellite."""
    if s.crm:
        s.why += "; decisive" if s.decisive else "; flags only (cores on too few chromosomes)"
        return s
    sig = s.significant()
    med = float(np.median([x.get("fold", TE_SPECIFIC_FOLD) for x in sig.values()])) if sig else 0.0
    if 2 * len(sig) <= s.n_chroms and not s.pooled:
        s.why += "; flags only (cores on too few chromosomes)"
    elif s.specificity < TE_SPECIFIC_FOLD:
        s.why += f"; flags only (under {TE_SPECIFIC_FOLD:.0f}-fold)"
    elif med < TE_SPECIFIC_FOLD:                       # the seed is specific but the recruited families' cores are not (Arachis:
        s.why += f"; flags only (its cores are {med:.1f}-fold on their chromosomes, under {TE_SPECIFIC_FOLD:.0f})"  # 11.6-fold seed, cores 1.5-fold)
    elif dominant and satellite_status in ("pass", "untested"):
        s.why += "; flags only (a dominant satellite falls one per chromosome)"
    else:
        s.decisive = True
        s.why += "; decisive"
    return s


def _age_from(x: dict) -> str:
    """Where the density maximum lay when age_choice moved the core (1-based inclusive), else '.'."""
    a = x.get("age_from")
    return f"{int(a[0]) + 1}-{int(a[1])}" if a else "."


def _alt(x: dict) -> str:
    """The second comparable density peak (second_peak), 1-based inclusive, else '.'."""
    a = x.get("alt")
    return f"{int(a[0]) + 1}-{int(a[1])}" if a else "."


def write_table(s: Signature, path) -> None:
    """One row per chromosome: the seed and what the signature found there (0-based half-open coordinates become 1-based
    inclusive, as in the satellite engine's summary)."""
    cols = ["chrom", "marker_lineage", "seed_level", "seed", "seed_specificity", "decisive", "core_start", "core_end",
            "core_elements", "chrom_elements", "core_p", "significant", "source", "n_families", "families", "why",
            "age_in_cores", "age_outside", "young_old", "core_fold", "core_age_from", "core_alt"]
    fmt = lambda x: "." if x != x else f"{x:.0f}" if abs(x) >= 100 else f"{x:.3g}"
    with atomic_open(path) as f:
        f.write("\t".join(cols) + "\n")
        for c, x in s.cores.items():
            fams = x.get("families", [])
            row = [c, s.lineage, s.level, s.seed, f"{s.specificity:.1f}", str(s.decisive), str(int(x["start"]) + 1), str(int(x["end"])),
                   str(int(x["k"])), str(int(x["n"])), f"{x['p']:.3g}", str(bool(x["significant"])), x.get("source", "."),
                   str(len(fams)), ",".join(fams[:20]) or ".", s.why, fmt(s.age_core), fmt(s.age_out), fmt(s.young_old),
                   f"{x.get('fold', float('nan')):.1f}", _age_from(x), _alt(x)]
            f.write("\t".join(row) + "\n")


def read_verdict(path) -> tuple:
    """(decisive, why) from a table of write_table; (None, '') when absent or empty."""
    try:
        with open(path) as f:
            head = f.readline().rstrip("\n").split("\t")
            first = f.readline().rstrip("\n").split("\t")
    except OSError:
        return None, ""
    if len(first) != len(head):
        return None, ""
    r = dict(zip(head, first))
    return r.get("decisive") == "True", r.get("why", "")
