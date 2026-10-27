"""S4: per-family metrics and the per-chromosome centromeric-family choice."""
import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .genomes import count_in
from .te import agrees
from .params import ACRO_ALPHA, FULL_MONOMER_COV, LOCUS_FOLD, RELATED_LEN_TOL, TANDEM_FRAC
from .util import clipped_bp, log, merge_intervals, union_bp

CLS_COLS = ["family", "bp", "n_arrays", "prevalence", "concentration", "subtel_frac", "gene_ratio",
            "median_main_bp", "eligible", "cls", "strength", "defines", "dominant", "reason", "subtelomeric",
            "rdna_tag", "role", "in_ltr", "ltr_clade", "tandem_frac", "whole_frac", "locus_fold", "acrocentric"]
# what a family is in one genome, first match wins (assign_roles); the tool's purpose is the centromere, the rest is
# what its tests found along the way
ROLES = ["rDNA", "centromeric", "cen_associated", "subtelomeric", "interstitial", "dispersed", "minor"]


@dataclass
class Cluster:
    start: int
    end: int
    bp: int


def make_clusters(iv, gap: int) -> list:
    """Arrays closer than gap form one cluster; bp = union bp of member arrays. Sorted by bp desc, start."""
    out, members, run_end = [], [], -1
    for s, e in sorted(iv):
        if members and s - run_end > gap:
            out.append(Cluster(members[0][0], run_end, union_bp(members)))
            members, run_end = [], -1
        members.append((s, e))
        run_end = max(run_end, e)
    if members:
        out.append(Cluster(members[0][0], run_end, union_bp(members)))
    return sorted(out, key=lambda c: (-c.bp, c.start))


def end_zone(L: int, p) -> int:
    """The terminal subtel_frac_len of a chromosome at each end. No absolute floor: a 1-Mb floor made every chromosome
    <= 2 Mb all end zone (Cryptococcus) and covered a CENH3 centromere 1.0 Mb from a chromosome end (Fragaria x
    ananassa Chr1D, Type II satellite; Jin 2025)."""
    return int(p.subtel_frac_len * L)


def terminal_bp(iv, L: int, p, mids=None) -> int:
    """bp of iv in the subtelomeres: the stretch between each chromosome end and its outermost gene (the textbook
    subtelomere). A centromere has a gene-bearing arm on both sides even next to an end (Fragaria x ananassa Chr1D:
    CENH3-bound Type II array 1.0 Mb from the end, 5 genes beyond it); subtelomeric repeats have none beyond them,
    whether at one end or both (rice TrsC, Solanum 182-bp blocks). Without >= 2 gene loci on the sequence: the
    terminal subtel_frac_len at each end."""
    if mids is not None and len(mids) >= 2 and mids[0] < mids[-1]:
        lo, hi = int(mids[0]), int(mids[-1])
    else:
        z = end_zone(L, p)
        lo, hi = z, L - z
    return clipped_bp(iv, 0, lo) + clipped_bp(iv, hi, L)


def near_end_sides(iv, L: int, p, mids=None) -> tuple:
    """(left, right) bp of iv near each chromosome end: in the subtelomere (terminal_bp) or within the terminal end
    zone. Gene models can reach into a subtelomere (Solanum acaule SG3 183-bp blocks), so a gene-bounded array
    this close to an end is ambiguous position evidence, not direct evidence of an interior locus."""
    z = end_zone(L, p)
    lo, hi = z, L - z
    if mids is not None and len(mids) >= 2 and mids[0] < mids[-1]:
        lo, hi = max(lo, int(mids[0])), min(hi, int(mids[-1]))
    return clipped_bp(iv, 0, lo), clipped_bp(iv, max(lo, hi), L)


def near_end_bp(iv, L: int, p, mids=None) -> int:
    return min(union_bp(iv), sum(near_end_sides(iv, L, p, mids)))


def interval_gap(a, b) -> int:
    return max(0, max(a.start, b.start) - min(a.end, b.end))


def related(f: str, g: str, fam_info) -> bool:
    """Families that may share a centromere as layers: same superfamily or period within 15%
    (CentBr1/2, CentO 155/165, CEN180/AthCEN159). Without family info every pair counts as related."""
    if fam_info is None:
        return True
    a, b = fam_info[f], fam_info[g]
    return a["superfamily"] == b["superfamily"] or abs(a["length"] / b["length"] - 1) <= RELATED_LEN_TOL


def same_sequence(f: str, g: str, fam_info) -> bool:
    """Families related by sequence: same superfamily (>= RELATED_IDENT) or chain group (close subfamilies, CentO
    155/165). Period alone is not kinship: the Populus 148-bp knob repeat and the 144/156-bp centromeric satellites
    share no sequence. Without family info every pair counts as related."""
    if fam_info is None:
        return True
    a, b = fam_info[f], fam_info[g]
    return a["superfamily"] == b["superfamily"] or (a.get("chain_group", f) == b.get("chain_group", g))


def terminal_share(arrays: pd.DataFrame, chroms: dict, genes, p) -> dict:
    """{family: (bp in subtelomeres, total bp)} of one genome (terminal_bp)."""
    out = {}
    for fam, g in arrays.groupby("family"):
        s = b = 0
        for c, d in g.groupby("chrom"):
            iv = list(zip(d.start, d.end))
            s += terminal_bp(iv, chroms[c], p, genes.get(c) if genes else None)
            b += union_bp(iv)
        out[fam] = (s, b)
    return out


def family_metrics(arrays: pd.DataFrame, chroms: dict, genes, p, fam_info=None, subtel_vote: dict = None,
                   tandem: dict = None):
    """fam_info: {fam: {superfamily, length, rdna}} from the families stage (optional). A family with most of its bp in
    subtelomeres (terminal_bp) is subtelomeric. subtel_vote: {fam: (genomes where it is subtelomeric, genomes that
    carry it)} over the complex. Position evidence has three grades. Interior (most bp away from both ends: neither
    beyond the outermost genes nor in the terminal end zone) is direct evidence and is never overruled. Terminal
    (beyond the outermost genes) and near-end (gene-bounded but in the end zone, near_end_bp) are ambiguous: a
    truncated assembly end or gene models reaching into a subtelomere. Such a family is subtelomeric here only if it
    is terminal in at least half of the genomes. Cases: Populus clone 741 chr14A/14D (terminal in 1 genome: kept);
    the RRCC radish subgenome, whose centromeres all end their chromosomes, does not make CL1 subtelomeric in the
    interior R. sativus; S. acaule SG3 183-bp blocks (near-end, terminal elsewhere: rejected); TrsC and the Solanum
    182-bp repeat (terminal in most genomes: rejected). A near-end or terminal family with blocks at an end of most
    of this genome's chromosomes is subtelomeric here whatever the vote: a truncated assembly end is one chromosome's
    (Populus: 1 end each), a subtelomeric family sits at the ends of most (M. sookensis G 687 bp on 10/14 while
    Cent728 holds the centres; Brassica CL25-like ~350 bp, Triticeae pSc119.2-like 118-120 bp, Nicotiana HRS60-like
    180 bp).

    Two tests added after the T2T mammal benchmark (docs/design.md):
    - tandem: {family: (tandem_frac, whole_frac)} from tandem_support. A family with neither most of its array bp
      inside de novo tandem repeats nor most of its hit bp in whole units, and not concentrated at a locus either
      (locus_fold: its main cluster's density over its mean density on the chromosome, median over chromosomes,
      below LOCUS_FOLD), is a dispersed repeat, never centromeric. Retrotransposon-derived arrays that cluster in
      centromeres (Coffea, Glycine, Rosa) are not tandem either, but they are concentrated (fold >= 2.5).
    - acrocentric: a family that would be subtelomeric is the centromere of an acrocentric karyotype when it is
      strong, sits at an end of most chromosomes, and at both ends of so few that the two ends cannot be independent
      (ACRO_ALPHA), and no strong interior family competes. Subtelomeric repeats occupy ends independently of each
      other; a centromere occurs once per chromosome (goat satellite I: one end of 29/31 chromosomes, both ends of
      0)."""
    total = sum(chroms.values())
    n_genes = sum(len(v) for v in genes.values()) if genes else 0
    density = n_genes / total if n_genes else None
    rows, mains = [], {}
    for fam, g in arrays.groupby("family"):
        bp = subtel = near = end_chroms = both_chroms = gcount = gspan = 0
        shares, main_bps, folds = [], [], []
        for c, d in g.groupby("chrom"):
            iv = list(zip(d.start, d.end))
            bpc = union_bp(iv)
            bp += bpc
            gc = genes.get(c) if genes else None
            subtel += terminal_bp(iv, chroms[c], p, gc)
            sides = near_end_sides(iv, chroms[c], p, gc)
            near += min(bpc, sum(sides))
            end_chroms += max(sides) >= p.chrom_min_bp            # a terminal block at either end
            both_chroms += min(sides) >= p.chrom_min_bp           # ... and at both ends
            m = make_clusters(iv, p.cluster_gap)[0]
            mains[(fam, c)] = m
            if bpc >= p.chrom_min_bp:
                shares.append(m.bp / bpc)
                main_bps.append(m.bp)
                folds.append((m.bp / bpc) / (max(1, m.end - m.start) / chroms[c]))
                if density:
                    gcount += count_in(genes.get(c), m.start, m.end)
                    gspan += m.end - m.start
        prevalence = len(shares) / len(chroms)
        conc = float(np.median(shares)) if shares else float("nan")
        subtel_frac = subtel / bp if bp else float("nan")
        gene_ratio = (gcount / gspan) / density if (density and gspan) else float("nan")
        fails = []
        rdna = fam_info[fam]["rdna"] if fam_info and fam in fam_info else "."
        if rdna != ".":
            fails.append(f"rDNA ({rdna})")
        if bp < p.chrom_min_bp:
            fails.append(f"bp={bp}<{p.chrom_min_bp}")
        if not conc >= p.min_concentration:
            fails.append(f"concentration={conc:.2f}<{p.min_concentration}")
        vote = (subtel_vote or {}).get(fam)
        terminal = not subtel_frac < p.max_subtel_frac
        near_frac = near / bp if bp else float("nan")
        at_ends = terminal or near_frac >= p.max_subtel_frac
        most_ends = 2 * end_chroms > len(chroms)                 # terminal blocks on most chromosomes of this genome
        if vote:
            subtelomeric = at_ends and (2 * vote[0] >= vote[1] or most_ends)
        else:
            subtelomeric = terminal or (at_ends and most_ends)
        subtel_why = None
        if subtelomeric:
            why = (f"subtelomeric in {vote[0]}/{vote[1]} genomes" if vote and 2 * vote[0] >= vote[1] else
                   f"terminal blocks on {end_chroms}/{len(chroms)} sequences" if most_ends else f">={p.max_subtel_frac}")
            subtel_why = (f"subtel_frac={subtel_frac:.2f}" + ("" if terminal else f", near ends {near_frac:.2f}")
                          + f", {why}")
        if not math.isnan(gene_ratio) and gene_ratio >= p.max_gene_ratio:
            fails.append(f"gene_ratio={gene_ratio:.2f}>={p.max_gene_ratio}")
        tf, wf = (tandem or {}).get(fam, (float("nan"), float("nan")))
        fold = float(np.median(folds)) if folds else float("nan")
        if tf < TANDEM_FRAC and wf < TANDEM_FRAC and fold < LOCUS_FOLD:    # NaN (no seeds/monomers) compares False
            fails.append(f"dispersed repeat: {tf:.0%} of bp in tandem repeats, {wf:.0%} whole-unit hits "
                         f"(< {TANDEM_FRAC:.0%} each), main cluster {fold:.1f}x its mean density (< {LOCUS_FOLD:g}x)")
        rows.append(dict(family=fam, bp=bp, n_arrays=len(g), prevalence=prevalence, concentration=conc,
                         subtel_frac=subtel_frac, gene_ratio=gene_ratio,
                         median_main_bp=float(np.median(main_bps)) if main_bps else 0.0,
                         eligible=False, cls="other", strength=".", defines=".", dominant=False, reason=".",
                         subtelomeric=bool(subtelomeric), rdna_tag=rdna, role=".", in_ltr=float("nan"), ltr_clade=".",
                         tandem_frac=tf, whole_frac=wf, locus_fold=fold, acrocentric=False,
                         _fails=fails, _subtel=subtel_why, _ends=(end_chroms, both_chroms, most_ends), _near=near))
    strong = lambda r: r["prevalence"] >= p.strong_prevalence and r["median_main_bp"] >= p.strong_cluster_bp
    # an interior rival must be interior as a satellite kind (subfamilies pooled, as for dominance): a subfamily at the
    # inner edge of a terminal block is part of that block (T2T mouse: a 58-bp major-satellite subfamily 2-13 Mb from
    # the end of 10/20 chromosomes lies beyond the end zone while its kind is 61% terminal)
    known = sorted(r["family"] for r in rows if fam_info is not None and r["family"] in fam_info)
    kind_of = {f: i for i, k in enumerate(kinds(known, fam_info)) for f in k}
    kbp, knear = {}, {}
    for r in rows:
        k = kind_of.setdefault(r["family"], r["family"])
        kbp[k] = kbp.get(k, 0) + r["bp"]
        knear[k] = knear.get(k, 0) + r["_near"]
    interior = any(strong(r) and not r["_fails"] and not r["_subtel"]
                   and knear[kind_of[r["family"]]] < p.max_subtel_frac * kbp[kind_of[r["family"]]] for r in rows)
    n = len(chroms)
    for r in rows:
        e, b, most = r["_ends"]
        if r["_subtel"] and not r["_fails"] and not interior and strong(r) and most:
            q = (e + b) / (2 * n)                               # per-end occupancy (one-ended + 2 x both-ended ends)
            pv = binom_cdf(b, n, q * q)                         # both ends occupied by chance, if ends independent
            if pv <= ACRO_ALPHA:
                r.update(subtelomeric=False, acrocentric=True, _subtel=None,
                         _acro=f"acrocentric: terminal at one end of {e - b}/{n} sequences, both ends of {b} "
                               f"(independent ends: {n * q * q:.1f} expected, P={pv:.1e}); no strong interior family")
        fails = r.pop("_fails") + ([r["_subtel"]] if r["_subtel"] else [])
        r["eligible"] = not fails
        r["reason"] = "; ".join(fails) if fails else "eligible" + (f"; {r['_acro']}" if r.get("_acro") else "")
        for k in ("_subtel", "_ends", "_acro", "_near"):
            r.pop(k, None)
    return pd.DataFrame(rows, columns=CLS_COLS), mains


def binom_cdf(k: int, n: int, q: float) -> float:
    """P(X <= k), X ~ Binomial(n, q)."""
    if q <= 0:
        return 1.0
    if q >= 1:
        return 1.0 if k >= n else 0.0
    lq, l1q = math.log(q), math.log1p(-q)
    return min(1.0, sum(math.exp(math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * lq
                                 + (n - i) * l1q) for i in range(0, min(k, n) + 1)))


def _merged(iv):
    """Disjoint sorted intervals: (starts, ends, bp before each)."""
    m = np.array(merge_intervals(iv), dtype=np.int64).reshape(-1, 2)
    return m[:, 0], m[:, 1], np.concatenate([[0], np.cumsum(m[:, 1] - m[:, 0])[:-1]])


def _covered(S, E, before, s, e):
    """bp of [s, e) (arrays) inside the disjoint intervals (S, E, before)."""
    def upto(x):                                          # covered bp in [0, x)
        i = np.searchsorted(S, x, side="right") - 1
        out = np.zeros(len(x), np.int64)
        ok = i >= 0
        out[ok] = before[i[ok]] + np.minimum(x[ok], E[i[ok]]) - S[i[ok]]
        return out
    return upto(e) - upto(s)


def tandem_support(arrays: pd.DataFrame, monomers: pd.DataFrame, seeds: pd.DataFrame, p) -> dict:
    """{family: (tandem_frac, whole_frac)}. tandem_frac: share of the family's array bp inside de novo tandem repeats
    (seed intervals of any period: satellites built of short subunits, such as HSat2/3, are found at their subunit
    period). whole_frac: share of its monomer-hit bp in hits that cover
    >= FULL_MONOMER_COV of the model. A satellite is tandem: a tandem-repeat finder finds its arrays, and its hits are
    whole units head to tail. One tandem copy of an interspersed element (LINE, SINE or LTR fragment) can seed a
    family that nhmmer then finds genome-wide, chained into 'arrays' of partial hits: the T2T chimpanzee 1,523-bp
    family (719 Mb on every chromosome) has 7% of its bp in tandem repeats and 0.3% whole hits; every CENP-A-validated
    centromeric family of the mammal benchmark has >= 45% and >= 90%."""
    out = {}
    if arrays.empty:
        return out
    sm = {c: _merged(list(zip(d.start, d.end))) for c, d in seeds.groupby("chrom")} if len(seeds) else {}
    tot, cov = {}, {}
    for (fam, c), d in arrays.groupby(["family", "chrom"]):
        S, E, _ = _merged(list(zip(d.start, d.end)))
        tot[fam] = tot.get(fam, 0) + int((E - S).sum())
        if c in sm:
            cov[fam] = cov.get(fam, 0) + int(_covered(*sm[c], S, E).sum())
    whole = {}
    if monomers is not None and len(monomers):
        bp = (monomers.end - monomers.start).to_numpy()
        full = (monomers.hmm_to - monomers.hmm_from + 1).to_numpy() >= FULL_MONOMER_COV * monomers.model_len.to_numpy()
        g = pd.DataFrame(dict(family=monomers.family.to_numpy(), bp=bp, wbp=np.where(full, bp, 0)))
        s = g.groupby("family")[["bp", "wbp"]].sum()
        whole = {f: (r.wbp / r.bp if r.bp else 0.0) for f, r in s.iterrows()}
    for fam, b in tot.items():
        out[fam] = (cov.get(fam, 0) / b if b else 0.0, whole.get(fam, float("nan")))
    return out


def kinds(families, fam_info) -> list:
    """Satellite kinds: connected components of `same_sequence` (subfamilies of one satellite, e.g. the 6 alpha
    subfamilies of a primate genome, or CentO 155/165). Sorted lists, in first-seen order."""
    out = []
    for f in families:
        hit = [k for k in out if any(same_sequence(f, x, fam_info) for x in k)]
        out = [k for k in out if k not in hit] + [sorted({f, *(x for k in hit for x in k)})]
    return out


def acro_rivals(acro, primary, fam_info) -> frozenset:
    """Acrocentric families unrelated by sequence to every centromeric family (`primary`): another satellite that
    also occurs once per chromosome at its centromeric end, so the sequence line cannot tell which one carries the
    kinetochore (centromere.call_regions). A family missing from fam_info counts as related (no rival)."""
    info = fam_info or {}
    kin = lambda f, x: f not in info or x not in info or same_sequence(f, x, info)
    return frozenset(f for f in set(acro) - set(primary) if not any(kin(f, x) for x in primary))


def largest_kinds(arrays: pd.DataFrame, elig: set, chroms: dict, p, fam_info) -> list:
    """Per chromosome, the satellite kind with the largest main cluster (arrays of all its eligible families pooled;
    >= chrom_min_bp), as a list of family sets."""
    a = arrays[arrays.family.isin(elig)]
    kind_of = {f: i for i, k in enumerate(kinds(sorted(elig), fam_info)) for f in k}
    members = {}
    for f, i in kind_of.items():
        members.setdefault(i, set()).add(f)
    out = []
    for c, d in a.groupby("chrom"):
        if c not in chroms:
            continue
        best = None
        for i, x in d.groupby(d.family.map(kind_of)):
            m = make_clusters(list(zip(x.start, x.end)), p.cluster_gap)[0]
            if m.bp >= p.chrom_min_bp and (best is None or (m.bp, -i) > (best[0], -best[1])):
                best = (m.bp, i)
        if best:
            out.append(members[best[1]])
    return out


def choose_centromeric(fams: pd.DataFrame, mains: dict, chroms: dict, p, fam_info=None, arrays=None):
    """arrays (with fam_info): the largest cluster per chromosome is judged per satellite kind (largest_kinds), so a
    satellite the families stage split into subfamilies is not outcounted by one pericentric family (T2T bonobo:
    alpha in 6 subfamilies lost dominance to a 43-bp satellite, which then took 12/25 calls 1-5 Mb from CENP-A).
    Without arrays: per family, as before."""
    elig = set(fams.loc[fams.eligible.astype(bool), "family"])
    by_chrom = {c: sorted(((mains[(f, c)].bp, f) for f in elig
                           if (f, c) in mains and mains[(f, c)].bp >= p.chrom_min_bp), reverse=True) for c in chroms}
    # genome-wide dominant family: strong AND its kind is the largest cluster on >= strong_prevalence of
    # chromosomes; presence alone let a widespread pericentromeric satellite (B. napus A06, 228 bp) take precedence
    if arrays is not None and fam_info is not None:
        largest = largest_kinds(arrays, elig, chroms, p, fam_info)
    else:
        largest = [{cands[0][1]} for cands in by_chrom.values() if cands]
    dominant = {r.family for r in fams.itertuples(index=False) if r.family in elig
                and r.prevalence >= p.strong_prevalence and r.median_main_bp >= p.strong_cluster_bp
                and sum(any(related(w, r.family, fam_info) for w in k) for k in largest)
                >= p.strong_prevalence * len(chroms)}
    per_chrom = {}
    for c, cands in by_chrom.items():
        if not cands:
            continue
        # a family related to a genome-wide dominant family wins where it has a qualifying cluster;
        # otherwise a larger pericentromeric satellite (e.g. TR238 next to CentBr) would define the centromere
        pref = [bf for bf in cands if any(related(bf[1], s, fam_info) for s in dominant)]
        best = (pref or cands)[0][1]
        bc = mains[(best, c)]
        per_chrom[c] = [best] + [f for _, f in cands if f != best and interval_gap(mains[(f, c)], bc) <= p.cluster_gap
                                 and related(f, best, fam_info)]
    defines = {}
    for c, fl in per_chrom.items():
        for f in fl:
            defines.setdefault(f, []).append(c)
    fams = fams.copy()
    fams["dominant"] = fams.family.isin(dominant)
    for i, r in fams.iterrows():
        if r.family in defines:
            strong = r.prevalence >= p.strong_prevalence and r.median_main_bp >= p.strong_cluster_bp
            fams.loc[i, ["cls", "strength", "defines", "reason"]] = [
                "centromeric", "strong" if strong else "weak", ",".join(defines[r.family]),
                f"eligible; centromeric family on {len(defines[r.family])} chromosome(s)"]
        elif r.eligible:
            fams.loc[i, "reason"] = "eligible; never the main or co-located cluster on a chromosome"
    return fams, per_chrom


def apply_rdna(fams: pd.DataFrame, rdna: dict, p, elsewhere: dict = None) -> pd.DataFrame:
    """Families whose largest arrays sit among 45S rRNA genes (rdna_neighbourhood) are rDNA spacer repeats (NORs):
    never centromeric. elsewhere: {fam: (genome, k, n)} tagged so in another genome of the complex. rRNA genes are
    often missing from assemblies, so their presence next to the arrays in one genome is the evidence, and their
    absence in another is not (Lemna: one family flanked by 45S genes at 5/9 and 6/10 arrays in two genomes, 1/6 in a
    third, where it was called centromeric)."""
    fams = fams.copy()
    if "rdna_tag" not in fams:
        fams["rdna_tag"] = "."
    for i, r in fams.iterrows():
        n, k = rdna.get(r.family, (0, 0))
        tag = r.rdna_tag if r.rdna_tag != "." else "45S_spacer"
        if n and k / n >= p.rdna_array_frac:
            fams.loc[i, ["eligible", "reason", "rdna_tag"]] = [
                False, f"rDNA-associated: 45S rRNA gene within {p.rdna_flank} bp of {k}/{n} largest arrays "
                       f"(rDNA spacer repeat)", tag]
        elif elsewhere and r.family in elsewhere:
            fams.loc[i, "rdna_tag"] = tag                    # one family = one sequence: the tag holds everywhere
            if r.eligible:
                g, k2, n2 = elsewhere[r.family]
                fams.loc[i, ["eligible", "reason"]] = [
                    False, f"rDNA-associated in {g}: 45S rRNA gene within {p.rdna_flank} bp of {k2}/{n2} largest "
                           f"arrays there"]
    return fams


def borrow_centromeric(fams: pd.DataFrame, mains: dict, sources: dict, p):
    """Contig assemblies: the per-chromosome tests (prevalence, end zones, concentration) are meaningless on
    contigs, so the target takes the union of its source genomes' centromeric families (strongest strength);
    its own metrics (eligible, reason) stay in the table for reference only. sources: {genome: classify table}.
    defines = target sequences where the family has a main cluster >= chrom_min_bp."""
    rank = {"strong": 2, "weak": 1}
    cen = {}
    for t in sources.values():
        for r in t[t.cls == "centromeric"].itertuples(index=False):
            if rank.get(r.strength, 0) >= rank.get(cen.get(r.family), 0):
                cen[r.family] = r.strength
    src = "+".join(sources)
    missing = sorted(set(cen) - set(fams.family))
    if missing:
        log.info(f"classify: centromeric in {src} but without arrays in the target: {', '.join(missing)}")
    per_chrom = {}
    fams = fams.copy()
    for i, r in fams.iterrows():
        if r.family in cen:
            chroms = [c for (f, c), m in mains.items() if f == r.family and m.bp >= p.chrom_min_bp]
            for c in chroms:
                per_chrom.setdefault(c, []).append(r.family)
            fams.loc[i, ["cls", "strength", "defines", "reason"]] = [
                "centromeric", cen[r.family], ",".join(chroms) or ".",
                f"borrowed from {src} ({cen[r.family]} there); own contig metrics: {r.reason}"]
        else:
            fams.loc[i, "reason"] = f"not centromeric in {src}; own contig metrics: {r.reason}"
    return fams, per_chrom


def assign_cen_associated(fams: pd.DataFrame, arrays: pd.DataFrame, regions: dict, p) -> pd.DataFrame:
    reg = {c: (r["cen_start"], r["cen_end"]) for c, r in regions.items() if r["evidence"] != "none"}
    fams = fams.copy()
    for i, r in fams.iterrows():
        if r.cls == "centromeric" or not r.bp:
            continue
        d = arrays[arrays.family == r.family]
        inside = sum(clipped_bp(list(zip(x.start, x.end)), *reg[c]) for c, x in d.groupby("chrom") if c in reg)
        if inside / r.bp >= p.cen_assoc_frac:
            fams.loc[i, ["cls", "reason"]] = ["cen_associated",
                                              f"{r.reason}; {inside / r.bp:.0%} of bp inside centromeres"]
    return fams


def assign_roles(fams: pd.DataFrame, p, borrowed: bool = False) -> pd.DataFrame:
    """One role per family in this genome (ROLES, first match wins), from the tests classify already ran:
    rDNA (5S rRNA by Rfam on the consensus, or 45S spacer: arrays among 45S rRNA genes), centromeric,
    cen_associated (most bp inside centromere calls), subtelomeric (the family_metrics verdict), interstitial (concentrated at one
    locus per chromosome, not the centromere: knobs and other heterochromatic blocks land here; the engine cannot tell
    them apart), dispersed (spread along chromosomes), minor (< chrom_min_bp in this genome: too little to judge).
    Borrowed genomes (contigs): per-chromosome tests are meaningless, so only rDNA / centromeric / cen_associated
    are assigned and the rest is 'unassessed'."""
    fams = fams.copy()
    roles = []
    for r in fams.itertuples(index=False):
        tag = getattr(r, "rdna_tag", ".")
        if isinstance(tag, str) and tag not in (".", ""):
            roles.append("rDNA")
        elif r.cls in ("centromeric", "cen_associated"):
            roles.append(r.cls)
        elif borrowed:
            roles.append("unassessed")
        elif bool(getattr(r, "subtelomeric", False)):
            roles.append("subtelomeric")
        elif r.bp < p.chrom_min_bp:
            roles.append("minor")
        elif r.concentration >= p.min_concentration:
            roles.append("interstitial")
        else:
            roles.append("dispersed")
    fams["role"] = roles
    return fams


def te_rejected(fams: pd.DataFrame, mains: dict, cores, fam_info=None) -> dict:
    """Centromeric family kinds (sequence-related families: subfamilies, CentBr1/2, CentO 155/165) whose main clusters
    disagree with the retrotransposon core on the majority of the chromosomes they define, for EVERY concentrated
    clade (cores: one {chrom: core} per candidate clade, or a single dict). A chromosome agrees if any member of the
    kind agrees there; judged per kind, so a subfamily on a few chromosomes is not rejected alone; judged against
    every clade, because which concentrated clade is centromeric is not known. Returns {family: (disagreeing,
    tested)} (counts for the most supportive clade) for every member of each rejected kind."""
    cands = [cores] if isinstance(cores, dict) else list(cores)
    cen = [r for r in fams[fams.cls == "centromeric"].itertuples(index=False)]
    kinds = []                                           # connected components of `same_sequence`
    for r in cen:
        hit = [k for k in kinds if any(same_sequence(r.family, x.family, fam_info) for x in k)]
        merged = [r] + [x for k in hit for x in k]
        kinds = [k for k in kinds if k not in hit] + [merged]
    out = {}
    for k in kinds:
        tallies = []
        for cc in cands:
            votes = {}
            for r in k:
                for c in str(r.defines).split(","):
                    m, core = mains.get((r.family, c)), cc.get(c)
                    if core and m:
                        votes[c] = votes.get(c, False) or agrees(m.start, m.end, core)
            if votes:
                tallies.append((sum(not v for v in votes.values()), len(votes)))
        if tallies and all(bad * 2 > n for bad, n in tallies):
            best = min(tallies, key=lambda t: t[0] / t[1])
            out.update({r.family: best for r in k})
    return out


def choose_with_te(fams: pd.DataFrame, mains: dict, chroms: dict, p, fam_info, te, arrays=None):
    """choose_centromeric, then drop families away from the retrotransposon centromere and choose again, until no
    centromeric family disagrees: a family that was outcompeted by a rejected one (e.g. the CentT550-like satellite
    beside interstitial blocks in Aegilops) gets its chance. te: (clade, {chrom: significant core}, decisive) or
    None; only decisive evidence (CRM, te.marker) rejects. Returns (fams, per_chrom, {rejected family: (disagreeing, tested)})."""
    rejected = {}
    while True:
        out, per_chrom = choose_centromeric(fams, mains, chroms, p, fam_info, arrays)
        rej = te_rejected(out, mains, [te[1]], fam_info) if te and te[2] else {}
        if not rej:
            return out, per_chrom, rejected
        rejected.update(rej)
        fams = fams.copy()
        for i, r in fams.iterrows():
            if r.family in rej:
                k, n = rej[r.family]
                fams.loc[i, ["eligible", "reason"]] = [
                    False, f"main clusters away from the {te[0]} retrotransposon centromere on {k}/{n} chromosomes"]
