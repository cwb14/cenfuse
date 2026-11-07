"""S5: per-chromosome centromere region call: satellite clusters first, then tandem-only fallbacks."""
import re

import numpy as np
import pandas as pd
import pysam

from . import longunit
from .classify import end_zone, make_clusters, terminal_bp
from .genomes import count_in
from .params import TE_CENTROMERIC_CLADE, TE_DERIVED_FRAC
from .te import agrees, elements_in, gap, overrules
from .util import log, union_bp

REG_COLS = ["chrom", "chrom_len", "cen_start", "cen_end", "evidence", "confidence", "flags", "main_family",
            "main_family_bp", "sat_bp", "second_cluster_bp", "cen_families", "n_gaps", "gap_bp", "n_genes",
            "gene_ratio", "te_clade", "te_core_start", "te_core_end", "alt_start", "alt_end", "overruled_start",
            "overruled_end", "overruled_family", "te_core_elements", "te_call_elements", "rejected_family",
            "rejected_start", "rejected_end"]
LOWER = {"high": "medium", "medium": "low", "low": "low"}
RANK = {"low": 0, "medium": 1, "high": 2}
EMPTY = np.zeros(0, np.int64)


def satellite_call(arrays_c: pd.DataFrame, primary: set, L: int, p, chosen: set = None, core: dict = None):
    """The cluster of centromeric-family arrays (all families pooled, so layered satellites share one call) that holds
    the main cluster of the families classify chose on this chromosome (`chosen`), or else the largest. A decisive
    retrotransposon core (`core`) outranks the choice: when only the largest cluster agrees with it, the call is the
    largest (S. verrucosum chr07, CENH3 on the largest cluster at the CRM core). The largest other cluster is the
    competing candidate."""
    a = arrays_c[arrays_c.family.isin(primary)]
    if a.empty:
        return None
    cl = make_clusters(list(zip(a.start, a.end)), p.cluster_gap)
    m = cl[0]
    ch = a[a.family.isin(chosen)] if chosen else a.iloc[:0]
    if len(ch):
        anc = make_clusters(list(zip(ch.start, ch.end)), p.cluster_gap)[0]
        m = next(x for x in cl if x.start <= anc.start and x.end >= anc.end)
    if core is not None and m is not cl[0] and not agrees(m.start, m.end, core) and agrees(cl[0].start, cl[0].end, core):
        m = cl[0]
    cl = [m] + [x for x in cl if x is not m]
    if m.bp >= p.chrom_min_bp:
        ev = "satellite"
    elif m.bp >= p.fallback_min_bp:
        ev = "weak_satellite"
    else:
        return None
    inside = a[(a.start >= m.start) & (a.end <= m.end)]
    fam_bp = {f: union_bp(zip(g.start, g.end)) for f, g in inside.groupby("family")}
    main = max(fam_bp, key=lambda f: (fam_bp[f], f))
    out = dict(cen_start=int(m.start), cen_end=int(m.end), evidence=ev, main_family=main,
               main_family_bp=fam_bp[main], sat_bp=m.bp, second_cluster_bp=cl[1].bp if len(cl) > 1 else 0,
               cen_families=",".join(sorted(fam_bp)))
    if len(cl) > 1 and cl[1].bp >= p.competing_ratio * m.bp:      # the competing candidate, reported for other lines
        out.update(alt_start=int(cl[1].start), alt_end=int(cl[1].end))
    return out


def longunit_call(lu_c: pd.DataFrame, L: int, genes_c, gene_density, p):
    """Fallback 2: long-unit arrays that pass the S4 eligibility tests define a low-confidence region."""
    if lu_c.empty:
        return None
    iv = list(zip(lu_c.start, lu_c.end))
    bpc = union_bp(iv)
    cl = make_clusters(iv, p.cluster_gap)
    m = cl[0]
    if m.bp < p.chrom_min_bp or m.bp / bpc < p.min_concentration or terminal_bp(iv, L, p, genes_c) / bpc >= p.max_subtel_frac:
        return None
    if genes_c is not None and gene_density:
        if (count_in(genes_c, m.start, m.end) / (m.end - m.start)) / gene_density >= p.max_gene_ratio:
            return None
    return dict(cen_start=int(m.start), cen_end=int(m.end), evidence="long_unit_tandem", main_family="long_unit",
                main_family_bp=m.bp, sat_bp=m.bp, second_cluster_bp=cl[1].bp if len(cl) > 1 else 0,
                cen_families="long_unit")


def weak_call(call: dict, p) -> bool:
    """A call of weak satellite evidence: a weak_satellite or long_unit_tandem fallback, or a cluster below low_sat_bp
    (the low_satellite flag)."""
    return call["evidence"] in ("weak_satellite", "long_unit_tandem") or call["sat_bp"] < p.low_sat_bp


def te_call(core: dict, replaced: dict = None) -> dict:
    """Centromere placed by the centromeric retrotransposon core alone (no satellite or long-unit call there).
    replaced: the call the core replaces; it is kept as overruled_start/overruled_end/overruled_family, because the
    outputs showed no trace of the array a core overruled. Never in alt_start/alt_end: fusion lets validated Hi-C move
    a call to the alternative cluster, so the replaced array there would change calls."""
    out = dict(cen_start=int(core["start"]), cen_end=int(core["end"]), evidence="retrotransposon", main_family=".",
               main_family_bp=0, sat_bp=0, second_cluster_bp=0, cen_families=".")
    if replaced:
        out.update(overruled_start=int(replaced["cen_start"]), overruled_end=int(replaced["cen_end"]),
                   overruled_family=replaced["main_family"])
    return out


def honesty_cap(te, borrowed: bool, primary: bool, dominant: bool) -> bool:
    """Whether a genome's satellite calls are capped at low: centromeric families exist but none dominates, so
    satellite evidence alone cannot tell a centromeric array from an interstitial one. Retrotransposon evidence lifts the
    cap genome-wide only when a CRM marker decides; any other marker lifts it per call, for te_concordant calls (finalize),
    because a line that contradicts a call, or says nothing about it, is no evidence for it (in the 185-genome batch a
    flags-only marker had lifted 47 te_discordant calls to medium/high). Borrowed genomes are never capped."""
    return not borrowed and primary and not dominant and not (te is not None and te[2] and
                                                              te[0].lower() == TE_CENTROMERIC_CLADE.lower())


def finalize(chrom: str, L: int, call, fa, genes_c, chrom_density, p, borrowed: bool = False, te_flag=None,
             capped: bool = False, te_decisive: bool = True, acro: frozenset = frozenset(),
             not_one: bool = False) -> dict:
    """borrowed (contig assemblies): near_end is replaced by contig_end (call within contig_end_bp of a contig
    end: the array may continue in another contig), which does not lower confidence. te_flag: te_concordant (at
    least medium: two independent lines agree) / te_discordant (lowers confidence one level); both only with
    te_decisive evidence (CRM), otherwise flags only, and a retrotransposon-only call is low. capped (honesty_cap): no
    dominant satellite and no decisive retrotransposon evidence, so satellite evidence alone cannot tell a centromeric
    array from an interstitial one: confidence low, unless a flags-only marker supports this call (te_concordant).
    A competing cluster (a second cluster at least half the main one) costs two levels: the satellite line cannot
    tell which array is the active one (T2T mammals: 2/6 such calls on the CENP-A domain), so the call is low and the
    second cluster is reported (alt_start/alt_end) for the other lines to choose. A rival satellite's cluster
    (call_regions) competes the same way, whatever its size. acro: acrocentric families (classify); a call on one is
    flagged acrocentric. not_one: the largest clusters of the centromeric families do not fall one per chromosome
    (karyotype.satellite_test: fail), so satellite evidence alone cannot say those arrays are centromeres: a satellite call is
    low unless a retrotransposon core agrees with it (flag not_one_per_chromosome)."""
    row = {k: "." for k in REG_COLS}
    row.update(chrom=chrom, chrom_len=L, evidence="none", confidence="none")
    if call is None:
        return row
    s, e = call["cen_start"], call["cen_end"]
    runs = [m.end() - m.start() for m in re.finditer("[Nn]{%d,}" % p.min_gap_n, fa.fetch(chrom, s, e))]
    flags = []
    if call["evidence"] == "retrotransposon":
        pass                                             # satellite-amount flags do not apply
    elif call["sat_bp"] < p.low_sat_bp:
        flags.append("low_satellite")
    if call["evidence"] != "retrotransposon" and (call.get("rival_bp") or call["second_cluster_bp"] > 0 and
                                                  call["second_cluster_bp"] >= p.competing_ratio * call["sat_bp"]):
        flags.append("competing_cluster")
    z = end_zone(L, p)
    if not borrowed and (s < z or e > L - z):
        flags.append("near_end")
    if call["evidence"] != "retrotransposon" and e - s >= p.diffuse_min_len and \
            call["sat_bp"] < p.diffuse_max_frac * (e - s):
        flags.append("diffuse")                        # broad, satellite-sparse block (e.g. Coffea heterochromatin)
    core = len(flags) + ("competing_cluster" in flags)
    if borrowed and (s < p.contig_end_bp or e > L - p.contig_end_bp):
        flags.append("contig_end")
    n_genes, gene_ratio = ".", "."
    if genes_c is not None:
        n_genes = count_in(genes_c, s, e)
        if chrom_density:
            gene_ratio = round((n_genes / (e - s)) / chrom_density, 3)
            if gene_ratio > 1:
                flags.append("gene_rich")
    if runs:
        flags.append(f"gaps:{len(runs)}")
    if call["evidence"] == "satellite":
        conf = ["high", "medium", "low", "low", "low"][core]
    elif call["evidence"] == "retrotransposon":
        conf = "medium" if te_decisive else "low"
    else:
        conf = "low"
    if call.get("main_family") in acro:
        flags.append("acrocentric")
    if te_flag:
        flags.append(te_flag)
    if te_flag and te_decisive:
        if te_flag == "te_discordant":
            conf = LOWER[conf]
        elif te_flag == "te_concordant" and RANK[conf] < RANK["medium"]:
            conf = "medium"                  # two independent lines agree: at least as sure as the core alone
    if capped and call["evidence"] != "retrotransposon" and te_flag != "te_concordant":
        flags.append("no_dominant_satellite")
        conf = "low"
    if not_one and call["evidence"] in ("satellite", "weak_satellite") and te_flag != "te_concordant":
        flags.append("not_one_per_chromosome")
        conf = "low"
    row.update(call)
    row.update(confidence=conf, flags=",".join(flags) or ".", n_gaps=len(runs), gap_bp=sum(runs),
               n_genes=n_genes, gene_ratio=gene_ratio)
    return row


def derived_families(fams: pd.DataFrame) -> dict:
    """{family: clade} of the families mostly inside LTR retrotransposons (TE_DERIVED_FRAC; only known with tes/)."""
    if not {"in_ltr", "ltr_clade"} <= set(fams.columns):
        return {}
    share = pd.to_numeric(fams.in_ltr, errors="coerce")
    return {f: c for f, c, s in zip(fams.family, fams.ltr_clade, share) if s >= TE_DERIVED_FRAC}


def remeasure(r: dict, fa, genes_c, arrays_c: pd.DataFrame, primary: set, p) -> None:
    """A call's interval changed after finalize: its columns and position flags describe the new one. Confidence stays."""
    s, e, L = r["cen_start"], r["cen_end"], r["chrom_len"]
    runs = [m.end() - m.start() for m in re.finditer("[Nn]{%d,}" % p.min_gap_n, fa.fetch(r["chrom"], s, e))]
    flags = [f for f in str(r["flags"]).split(",") if f not in (".", "", "near_end", "gene_rich") and not f.startswith("gaps:")]
    z = end_zone(L, p)
    if s < z or e > L - z:
        flags.append("near_end")
    n_genes, gene_ratio = ".", "."
    if genes_c is not None:
        n_genes = count_in(genes_c, s, e)
        if len(genes_c):
            gene_ratio = round((n_genes / (e - s)) / (len(genes_c) / L), 3)
            if gene_ratio > 1:
                flags.append("gene_rich")
    if runs:
        flags.append(f"gaps:{len(runs)}")
    a = arrays_c[arrays_c.family.isin(primary) & (arrays_c.start < e) & (arrays_c.end > s)]
    sat = union_bp((max(x, s), min(y, e)) for x, y in zip(a.start, a.end))
    main = a[a.family == r["main_family"]]
    r.update(flags=",".join(flags) or ".", n_gaps=len(runs), gap_bp=sum(runs), n_genes=n_genes, gene_ratio=gene_ratio, sat_bp=sat,
             main_family_bp=union_bp((max(x, s), min(y, e)) for x, y in zip(main.start, main.end)))


def cap_chain_calls(rows: dict, arrays: pd.DataFrame, primary: set, chains: set, derived: dict, decisive: bool, p,
                    cores: dict = None, fa=None, genes: dict = None) -> int:
    """Calls that stand on a chain of partial hits (classify.partial_chain) are not satellite calls: flag partial_chain,
    and at most `medium` where a retrotransposon core overlaps the call, `low` otherwise. Only CRM may decide among
    retrotransposons, so a core of another clade counts only as a second line, not if the family is made of that clade
    (derived: one piece of evidence twice). Overlap, not te.agrees: its gap may be as wide as the call, and a call
    spanning 23% of a chromosome agrees with about 70% of random places. A call holding chrom_min_bp of real satellite is
    a layered satellite call and stays as it is. Applied after finalize, so that the call a core replaces stays what it was:
    this lowers confidence and, where a decisive core inside the call is narrower, narrows the call to it (chain_to_core:
    such cores are more gene-poor than the chains, with the same homoeolog agreement; 7 of 7 right on the genomes with truth,
    before and after). A core that may not act on its chromosome (cores: its marker's, strong False) does not narrow. fa, genes:
    to measure a narrowed call again (remeasure). Returns the calls flagged."""
    real = arrays[arrays.family.isin(set(primary) - set(chains))]
    n = 0
    for c, r in rows.items():
        if r["evidence"] not in ("satellite", "weak_satellite") or r["main_family"] not in chains:
            continue
        s, e = r["cen_start"], r["cen_end"]
        x = real[(real.chrom == c) & (real.start >= s) & (real.end <= e)]
        if union_bp(zip(x.start, x.end)) >= p.chrom_min_bp:
            continue
        clade = str(r.get("te_clade", "."))
        on_core = r.get("te_core_start", ".") != "." and gap(s, e, int(r["te_core_start"]), int(r["te_core_end"])) == 0
        independent = decisive or derived.get(r["main_family"], "").lower() != clade.lower()
        most = "medium" if on_core and independent else "low"
        if RANK[r["confidence"]] > RANK[most]:
            r["confidence"] = most
        r["flags"] = "partial_chain" if r["flags"] in (".", "") else r["flags"] + ",partial_chain"
        inside = on_core and s <= int(r["te_core_start"]) and int(r["te_core_end"]) <= e
        if inside and decisive and (cores or {}).get(c, {}).get("strong", True) and int(r["te_core_end"]) - int(r["te_core_start"]) < e - s:
            r["cen_start"], r["cen_end"] = int(r["te_core_start"]), int(r["te_core_end"])    # a chain says where the repeat
            if fa is not None:                                                                # gathers, the core how narrowly
                remeasure(r, fa, (genes or {}).get(c) if genes else None, arrays[arrays.chrom == c], primary, p)
            r["flags"] += ",chain_to_core"
            if r.get("te_core_elements", ".") != ".":
                r["te_call_elements"] = r["te_core_elements"]
        n += 1
    return n


def rejected_cluster(arrays_c: pd.DataFrame, rejected: set, start: int, end: int, p):
    """(family, start, end) of the largest cluster, on one chromosome, of the families the CRM line rejected, when it is
    centromere-sized and lies outside the call [start, end): the array the call would otherwise have stood on. Else None."""
    a = arrays_c[arrays_c.family.isin(rejected)]
    if a.empty:
        return None
    m = make_clusters(list(zip(a.start, a.end)), p.cluster_gap)[0]
    if m.bp < p.chrom_min_bp or gap(m.start, m.end, start, end) == 0:
        return None
    inside = a[(a.start >= m.start) & (a.end <= m.end)]
    fam_bp = {f: union_bp(zip(g.start, g.end)) for f, g in inside.groupby("family")}
    return max(fam_bp, key=lambda f: (fam_bp[f], f)), int(m.start), int(m.end)


def call_regions(genome, arrays, primary: set, genes, seed_rows, tmpdir, threads: int, p, borrowed: bool = False,
                 te=None, capped: bool = False, acro: frozenset = frozenset(), per_chrom: dict = None,
                 rivals: frozenset = frozenset(), not_one: bool = False, weak_cores: dict = None):
    """borrowed (contig assemblies with families from --classify-from): only full satellite calls count. A contig
    without the centromeric satellite is no evidence of a satellite-free centromere, so there is no weak-satellite
    or long-unit fallback, and the long-unit scan covers call flanks only.
    te: (clade, {chrom: significant core}, decisive) of the centromere-marking LTR clade, or None. A call that
    agrees with its chromosome's core is te_concordant; one that does not is te_discordant, except a long-unit call,
    which the core replaces, and, with decisive evidence, a satellite call, which the core also replaces
    (satellite_away_from_core; Ae. geniculata 2U/7M: CENH3 on the CRM core, not on the distant array), unless the call
    is high or medium on its own and the core does not contradict it (te.contradicts: the clade does not mark the
    core significantly more than the call; B. rapa CC1 A09: a CRM plateau over 22-42 Mb, CentBr on CENH3): te_discordant;
    a call that agrees with a decisive core without overlapping it is replaced by the core when the core contradicts
    it in the same way (satellite_beside_core; Ae. tauschii Chr2/Chr4: arrays 5-7 Mb beside a 292-332-element CRM
    core holding CENH3; B. rapa A06: CentBr call holds more CRB than the core, kept); a chromosome
    without a call takes the core; so does a weak satellite call (weak_call) beside a decisive core, agreeing with it but
    not overlapping it (weak_satellite_beside_core: S. verrucosum chr01, Ae. geniculata chr07U, potato chr11 had the core
    on CENH3 and the weak array 0.8-5 Mb away). A call (satellite, weak satellite or long-unit) the core replaces is
    kept in overruled_start/overruled_end/overruled_family (te_call), never in alt_*; a kept call (te_concordant,
    te_discordant) and a chromosome without a call record nothing. Decisive evidence (CRM) sets confidence;
    other clades only flag. Without satellite calls the long-unit scan covers the core and its flanks instead of the
    whole chromosome. per_chrom: {chrom: the families classify chose there (its `defines`)}; the call is the cluster
    holding them (satellite_call), so a family chosen on a few chromosomes cannot take the call elsewhere where its
    cluster is merely larger (T2T siamang: an 80-bp satellite chosen on 2 chromosomes took 11 calls from alpha; CHM13
    chr21: the call took an inactive alpha array 2 Mb from the active one). Chromosomes without a choice take the
    largest cluster (weak-satellite fallback).
    rivals: acrocentric families unrelated to every centromeric family (classify). A satellite call competes with a
    rival's terminal block outside it (a cluster >= chrom_min_bp reaching the end zone, the blocks that made the
    rival acrocentric; reported in alt_start/alt_end): two satellites both occur once per chromosome at its
    centromeric end, and which carries the kinetochore is not in the sequence (T2T mouse: the 234-bp major
    satellite out-sizes the 120-bp minor satellite that carries CENP-A).
    not_one: the centromeric satellites failed the one-per-chromosome test (finalize). weak_cores: {chrom: core} of a decisive
    signature where its core is too thin to be significant: such a core is the better candidate where nothing else stands, on a
    chromosome without a call or with a satellite call of a family that failed the test. The call stays low (weak_core; potato:
    right on chr04, chr09 and chr12, whose satellite calls were not)."""
    clade, cores, decisive = te if te else (None, {}, False)
    per_chrom = per_chrom or {}
    tmpdir.mkdir(parents=True, exist_ok=True)
    total = sum(genome.chroms.values())
    gd = (sum(len(v) for v in genes.values()) / total) if genes else None
    by_chrom = dict(tuple(arrays.groupby("chrom"))) if len(arrays) else {}
    sat = {c: satellite_call(by_chrom[c], primary, L, p, per_chrom.get(c), cores.get(c) if decisive else None)
           if c in by_chrom else None for c, L in genome.chroms.items()}
    if borrowed:
        sat = {c: x if x and x["evidence"] == "satellite" else None for c, x in sat.items()}
    elif rivals:
        for c, x in sat.items():
            r = by_chrom[c][by_chrom[c].family.isin(rivals)] if x else None
            if r is None or r.empty:
                continue
            L = genome.chroms[c]
            z = end_zone(L, p)
            ends = [m for side in (r[r.end <= x["cen_start"]], r[r.start >= x["cen_end"]]) if len(side)
                    for m in make_clusters(list(zip(side.start, side.end)), p.cluster_gap)
                    if m.bp >= p.chrom_min_bp and (m.start < z or m.end > L - z)]     # its terminal block
            if ends:
                m = max(ends, key=lambda k: (k.bp, -k.start))
                x["rival_bp"] = m.bp
                if "alt_start" not in x:
                    x.update(alt_start=int(m.start), alt_end=int(m.end))
    def target(c, L):
        x = sat[c] or (dict(cen_start=cores[c]["start"], cen_end=cores[c]["end"]) if c in cores else None)
        return (c, max(0, x["cen_start"] - p.lu_flank), min(L, x["cen_end"] + p.lu_flank)) if x else (c, 0, L)
    targets = [target(c, L) for c, L in genome.chroms.items() if sat[c] or not borrowed]
    lu = longunit.detect(genome.fasta, targets, longunit.build_masks([arrays, seed_rows]), tmpdir, threads, p)
    rows = {}
    lu_by = dict(tuple(lu.groupby("chrom"))) if len(lu) else {}
    with pysam.FastaFile(str(genome.fasta)) as fa:
        for c, L in genome.chroms.items():
            gc = None if genes is None else genes.get(c, EMPTY)
            call = sat[c]
            if call is None and not borrowed:
                call = longunit_call(lu_by.get(c, lu.iloc[:0]), L, gc, gd, p)
            te_flag, core = None, cores.get(c)
            decisive_all, decisive = decisive, bool(decisive and (core or {}).get("strong", True))
            if core:
                if call and decisive and not borrowed and weak_call(call, p) and \
                        agrees(call["cen_start"], call["cen_end"], core) and \
                        gap(call["cen_start"], call["cen_end"], core["start"], core["end"]) > 0:
                    te_flag, call = "weak_satellite_beside_core", te_call(core, call)   # the core outranks a weak array
                elif call and agrees(call["cen_start"], call["cen_end"], core):
                    te_flag = "te_concordant"
                    if decisive and not borrowed and call["evidence"] in ("satellite", "weak_satellite") and \
                            gap(call["cen_start"], call["cen_end"], core["start"], core["end"]) > 0 and \
                            core.get("mids") is not None and overrules(call["cen_start"], call["cen_end"], core):
                        # the array lies beside the CRM centromere (Ae. tauschii Chr2, Chr4)
                        te_flag, call = "satellite_beside_core", te_call(core, call)
                elif call and call["evidence"] != "long_unit_tandem" and (not decisive or borrowed):
                    te_flag = "te_discordant"
                elif not borrowed:
                    if call and call["evidence"] in ("satellite", "weak_satellite"):
                        alone = finalize(c, L, call, fa, gc, (len(gc) / L) if gc is not None else None, p, borrowed,
                                         None, capped, decisive, acro, not_one)["confidence"]
                        if alone in ("high", "medium") and \
                                not overrules(call["cen_start"], call["cen_end"], core):
                            te_flag = "te_discordant"            # the clade marks the call's locus about as much as the core
                        else:
                            te_flag = "satellite_away_from_core"  # decisive core outranks a distant satellite array
                    if te_flag != "te_discordant":
                        call = te_call(core, call)                # a chromosome without a call overrules nothing
            thin = (weak_cores or {}).get(c) if decisive_all and not core and not borrowed else None
            if thin and (call is None or (not_one and call["evidence"] in ("satellite", "weak_satellite"))):
                te_flag, call, core, decisive = "weak_core", te_call(thin, call) if call else te_call(thin), thin, False
            rows[c] = finalize(c, L, call, fa, gc, (len(gc) / L) if gc is not None else None, p, borrowed,
                               te_flag, capped, decisive, acro, not_one)
            if core:
                rows[c].update(te_clade=core.get("clade", clade), te_core_start=int(core["start"]),
                               te_core_end=int(core["end"]))
                if rows[c]["evidence"] != "none" and core.get("mids") is not None:
                    # the two counts of the head-to-head test (te.contradicts): the core, and the call or what it replaced
                    over = rows[c].get("overruled_start", ".")
                    iv = (over, rows[c]["overruled_end"]) if over != "." else (rows[c]["cen_start"], rows[c]["cen_end"])
                    rows[c].update(te_core_elements=int(core["k"]), te_call_elements=elements_in(*iv, core))
            decisive = decisive_all
    log.debug(f"centromeres: {genome.name}: {len(lu)} long-unit arrays")
    return rows, lu


def write_regions(rows: dict, path) -> None:
    pd.DataFrame(list(rows.values()), columns=REG_COLS).to_csv(path, sep="\t", index=False)


def load_regions(path) -> dict:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    out = {}
    for r in df.to_dict("records"):
        r["chrom_len"] = int(r["chrom_len"])
        if r["evidence"] != "none":
            r["cen_start"], r["cen_end"] = int(r["cen_start"]), int(r["cen_end"])
        out[r["chrom"]] = r
    return out

