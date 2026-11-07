"""Interpretation (README, 'Reading a call'): what each call stands on and why it may be wrong, in words, and the run
report (cenfuse_report.md). Both are read from the tables cenfuse writes, so they say nothing the tables do not; nothing
here changes a call.

- basis: the sequence evidence (satellite family and unit, its bp in the call; a retrotransposon core), a move or
  placement by validated 3D lines and its second line, each supporting line, the CpG kinetochore.
- caveat: the array a TE core overruled, each confidence-lowering flag and opposing line, gaps in the call, and the
  input that would decide where cenfuse's own tests say which one can (CENP-A/CENH3 data for two arrays of one
  satellite: neither Hi-C nor CpG methylation chooses between them when no mappable bins separate them)."""
from __future__ import annotations

from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__
from .morphology import levan
from .params import AXIS_LEAN_P
from .sat.params import LEN_TOL
from .seqcalls import bp_by_role

CEN_ROLES = ("centromeric", "cen_associated")
LINE = {"te": "LTR retrotransposons (tes/)", "hic": "Hi-C", "meth": "methylation", "genes": "genes"}


def _toks(x) -> list:
    return [t for t in str(x).split(",") if t not in ("", ".", "nan")]


def _int(x):
    try:
        return int(x)
    except (TypeError, ValueError):
        return None


def size(bp) -> str:
    """bp as kb below 1 Mb, else Mb."""
    bp = int(bp)
    return f"{bp / 1e6:.2f} Mb" if bp >= 1_000_000 else f"{max(1, round(bp / 1e3)) if bp else 0} kb"


def interval(s, e) -> str:
    """1-based inclusive start, end -> 'a-b Mb'."""
    return f"{(int(s) - 1) / 1e6:.2f}-{int(e) / 1e6:.2f} Mb"


def context(seq, tel: dict, hic_given: bool, tes_given: bool, telo: dict | None = None) -> dict:
    """What the sentences need besides a row: family unit lengths, centromeric satellite bp in an interval, the LTR
    marker, which inputs were given, and the interstitial telomeric arrays."""
    unit = {}
    if len(seq.families) and {"family", "consensus_len"} <= set(seq.families.columns):
        unit = {r.family: _int(r.consensus_len) for r in seq.families.itertuples(index=False)}

    def sat_bp(chrom, s, e):
        return sum(v for k, v in bp_by_role(seq.arrays, chrom, s, e).items() if k in CEN_ROLES)
    m = tel.get("marker")
    sig = tel.get("signature")
    return dict(unit=unit, sat_bp=sat_bp, marker=m[0] if m else None, decisive=bool(m and m[2]), hic=hic_given,
                tes=tes_given, its=list((telo or {}).get("interstitial", [])), te_why=sig.why if sig is not None else "")


def _years(x) -> str:
    """'1234 (56)' (median years, copies; fuse._age) -> '1,234 years, 56 copies'."""
    try:
        y, n = str(x).split(" (")
        return f"{int(float(y)):,} years, {n.rstrip(')')} copies"
    except ValueError:
        return "?"


def _cis_trans(r) -> bool:
    """The row's cis signal was tested by Hi-C's trans contacts (hicfuse.cis_trans_line), not by LTR cores."""
    return str(getattr(r, "cis_test", ".")) == "trans"


def basis(r, ctx: dict) -> str:
    """What the call stands on, in words ('.' without a call)."""
    if r.relation == "none" or _int(r.start) is None:
        return "."
    flags, sflags, sup = set(_toks(r.flags)), set(_toks(r.seq_flags)), _toks(r.support)
    s, e = int(r.start) - 1, int(r.end)
    over = _int(getattr(r, "seq_overruled_start", None)), _int(getattr(r, "seq_overruled_end", None))
    where = f" ({interval(*over)})" if None not in over else ""    # where the replaced array lay (older tables: absent)
    out = []
    if r.call_source == "sequence":
        ev = r.seq_evidence
        if ev in ("satellite", "weak_satellite") and "partial_chain" in sflags:
            unit = ctx["unit"].get(r.main_family)
            out.append(f"repeat family {r.main_family}" + (f" ({unit}-bp model)" if unit else "") + ", chains of partial "
                       f"copies, {size(ctx['sat_bp'](r.chrom, s, e))} in the call")
        elif ev in ("satellite", "weak_satellite"):
            unit = ctx["unit"].get(r.main_family)
            sp = _int(getattr(r, "seq_spacing", None))               # the arrays' own spacing, where a profile fits
            spaced = (f", copies spaced {sp} bp here" if sp and unit and abs(sp - unit) > LEN_TOL * unit else "")  # others
            out.append(f"{'weak ' if ev == 'weak_satellite' else ''}satellite {r.main_family}"
                       + (f" ({unit}-bp unit{spaced})" if unit else "") + f", {size(ctx['sat_bp'](r.chrom, s, e))} of "
                       f"centromeric satellite in the call")
        elif ev == "long_unit_tandem":
            out.append(f"long-unit tandem array, {size(ctx['sat_bp'](r.chrom, s, e))} of satellite in the call")
        elif ev == "retrotransposon":
            replaced = sflags & {"weak_satellite_beside_core", "satellite_away_from_core", "satellite_beside_core"}
            seed = getattr(r, "te_seed", ".")
            k, n = _int(getattr(r, "te_core_elements", None)), _int(getattr(r, "te_chrom_elements", None))
            out.append(f"{r.te_clade} retrotransposon core" + ("" if replaced else " (no satellite call)") +
                       (f", signature seeded by {seed}" if seed not in (".", r.te_clade, None) else "") +
                       (f", {k} of its {n} elements here" if k is not None and n else ""))
        if "weak_satellite_beside_core" in sflags:
            out.append(f"the {r.te_clade} core replaced a weak satellite array beside it" + where)
        if "satellite_away_from_core" in sflags:
            out.append(f"the {r.te_clade} core replaced a satellite array elsewhere" + where)
        if "satellite_beside_core" in sflags:
            out.append(f"the {r.te_clade} core replaced a satellite array beside it that holds significantly fewer "
                       f"{r.te_clade} elements" + where)
        if "hic_chose_alt" in flags:
            out.append("the second satellite cluster, chosen by a strong Hi-C peak (Hi-C finds the other chromosomes' centromeres "
                       "beyond chance, but not on most)" if "hic_support" in flags else
                       "the second satellite cluster, chosen by validated Hi-C")
        if "acrocentric" in sflags:
            out.append("the family sits at one end of most chromosomes (acrocentric karyotype)")
    elif r.call_source == "composition":
        out.append(f"the deepest AT trough of this chromosome (GC {float(r.at_gc):.0%}), in a genome whose deepest troughs "
                   f"fall one per chromosome" if str(r.at_gate) == "pass" else
                   f"the deepest AT trough of this chromosome (GC {float(r.at_gc):.0%}), where Hi-C's prior-free peak agrees, "
                   f"in a genome where the two agree beyond chance" if "at_denovo_pair" in flags else
                   f"the deepest AT trough of this chromosome (GC {float(r.at_gc):.0%}), in a genome whose troughs agree "
                   f"with Hi-C's prior-free peaks")
    else:
        if "te_hic_pair" in flags:
            how = f"the {r.te_clade} core and Hi-C's peak, which agree here and beyond chance genome-wide"
        elif "cis_only" in flags:
            how = ("the extreme of the contact crossing each bin, a cis signal that Hi-C's trans contacts confirm on most chromosomes"
                   if _cis_trans(r) else
                   "the extreme of the contact crossing each bin, a cis signal that agrees with the LTR cores beyond chance genome-wide")
        elif "fold_pair" in flags:
            how = "the arm fold and Hi-C's prior-free peak, which agree here"
        elif "align_pair" in flags:
            how = "the trans arm alignment and Hi-C's prior-free peak, which agree here"
        elif "hic_medium_peak" in flags:
            how = "a clean medium validated Hi-C peak"
        elif "hic_denovo_only" in flags:
            axis_ok = getattr(r, "fold_gate", ".") == "pass" or getattr(r, "align_gate", ".") == "pass"
            how = ("Hi-C's prior-free peak alone (too few reference calls in this genome to test Hi-C against" +
                   ("; the Rabl-axis test passed genome-wide but does not agree here)" if axis_ok else ")"))
        else:
            how = "a strong validated Hi-C peak"
        if r.relation == "offset":
            off = _int(r.offset_bp)
            out.append(f"moved {size(off) + ' ' if off is not None else ''}from the sequence call "
                       f"({interval(r.seq_start, r.seq_end)}) by {how}")
        else:
            out.append(f"placed by {how} (no sequence call)")
        second = {"te": f"the {r.te_clade} core at the peak", "meth": "CHG methylation at the peak",
                  "fold": "the arm fold", "align": "the trans arm alignment",
                  "age": f"copies of the centromeric LTR families younger at the peak than at the sequence call (median "
                         f"{_years(getattr(r, 'age_peak', '.'))} vs {_years(getattr(r, 'age_call', '.'))}, P "
                         f"{getattr(r, 'age_p', '.')}): the centromere moved"}
        named = set()
        for f in sorted(flags):
            if f.startswith("moved_with_"):
                named |= set(f[len("moved_with_"):].split("+"))
                out += [f"second line: {second[x]}" for x in f[len("moved_with_"):].split("+") if x in second]
        if "te" in sup and "te" not in named:            # a call below high moved with the CRM core at the peak
            out.append(f"the {r.te_clade} core at the peak agrees")
        if "at_pair" in flags:
            out.append(f"an AT trough at the peak agrees (GC {float(r.at_gc):.0%}) and gives the interval")
        if "satellite" in sup:
            out.append(f"{size(ctx['sat_bp'](r.chrom, s, e))} of centromeric satellite inside")
    for line in sup:
        if line == "te" and r.call_source == "sequence" and r.seq_evidence != "retrotransposon":
            out.append(f"the {r.te_clade} core agrees")
        elif line == "hic" and r.call_source == "sequence":
            out.append("Hi-C contact is higher at this chromosome end than at the other, as on most chromosomes of this genome"
                       if "hic_end" in flags else
                       "the contact crossing each bin has its extreme here, a cis signal that Hi-C's trans contacts confirm on most "
                       "chromosomes" if "cis_agrees" in flags and _cis_trans(r) else
                       "the contact crossing each bin has its extreme here, a cis signal that agrees with the LTR cores beyond chance"
                       if "cis_agrees" in flags else
                       f"Hi-C's peak lies at the {r.te_clade} core (the two agree beyond chance genome-wide)"
                       if "te_hic_pair" in flags else
                       "Hi-C agrees (support only: its peaks find the other chromosomes' centromeres beyond chance, "
                       "but not on most)" if "hic_support" in flags else "validated Hi-C agrees")
        elif line == "meth":
            out.append("CHG methylation is centromere-like")
        elif line == "fold" and "fold_pair" not in flags:
            out.append("the arm fold agrees")
        elif line == "align" and "align_pair" not in flags:
            out.append("the trans arm alignment agrees")
    ks, ke = _int(getattr(r, "cpg_kt_start", None)), _int(getattr(r, "cpg_kt_end", None))
    if ks is not None and ke is not None:
        out.append(f"CpG kinetochore {interval(ks, ke)} ({size(ke - ks + 1)})")
    return "; ".join(out) or "."


def _could(r, ctx: dict) -> list:
    """Inputs not given that could confirm or move a call: Hi-C, an LTR core (tes/)."""
    out = []
    if not ctx["hic"]:
        out.append("Hi-C")
    if not ctx["tes"]:
        out.append("an LTR core (tes/; plants)")
    return out


def caveat(r, ctx: dict) -> str:
    """Why the call may be wrong, and the input that would decide ('.' when nothing is known against it)."""
    failed = "; Hi-C failed its test here" if ctx["hic"] and str(getattr(r, "hic_gate", ".")) == "fail" else ""
    if r.relation == "none" or _int(r.start) is None:
        need = _could(r, ctx)
        head = (f"no line supports the only candidate, a long-unit tandem array ({interval(r.seq_start, r.seq_end)})"
                if "long_unit_unsupported" in _toks(r.flags) else "no line found a centromere (satellite-free?)")
        return head + failed + (f" -> {' or '.join(need)} can place one" if need else "")
    flags, sflags, against = set(_toks(r.flags)), set(_toks(r.seq_flags)), _toks(r.against)
    out, hint = [], None
    over = _int(getattr(r, "seq_overruled_start", None)), _int(getattr(r, "seq_overruled_end", None))
    if None not in over:      # the call a TE core replaced: a satellite or a long-unit array (older tables: absent)
        out.append(f"the tandem array at {interval(*over)} was overruled by the "
                   f"{r.te_clade if _toks(r.te_clade) else 'TE'} core")
    rej = _int(getattr(r, "seq_rejected_start", None)), _int(getattr(r, "seq_rejected_end", None))
    if None not in rej:       # a family the LTR line rejected at classification, its array elsewhere (older tables: absent)
        cl = r.te_clade if _toks(r.te_clade) else "LTR"
        out.append(f"the {cl} line rejected {r.seq_rejected_family}, whose array at {interval(*rej)} lies away from the "
                   f"{cl} cores on most chromosomes")
    te_iv = (_int(r.te_start), _int(r.te_end))
    alt_iv = (_int(r.seq_alt_start), _int(r.seq_alt_end))
    te_in_alt = None not in te_iv + alt_iv and te_iv[0] <= alt_iv[1] and te_iv[1] >= alt_iv[0]
    if r.call_source == "sequence" and "competing_cluster" in sflags and "hic_chose_alt" not in flags:
        alt = f" ({interval(*alt_iv)})" if alt_iv[0] is not None else ""
        where = f", where the {r.te_clade} core lies" if te_in_alt else ""
        out.append(f"another satellite cluster{alt} competes{where}: the sequence cannot tell which is active")
        hint = "CENP-A/CENH3 data would decide" + ("" if ctx["tes"] else ", or an LTR core (tes/; plants)") + \
            ("" if ctx["hic"] else ", or Hi-C when mappable DNA separates them")
    if r.call_source == "sequence":
        if "partial_chain" in sflags:
            out.append("not a tandem satellite: chains of partial copies of a longer repeat, which mark where that repeat "
                       "concentrates, possibly beside the centromere")
            hint = hint or ("CENP-A/CENH3 data would decide" if ctx["tes"] else
                            "an LTR core (tes/; plants) or CENP-A/CENH3 data would decide")
        if "no_dominant_satellite" in sflags:
            out.append("no satellite family dominates the genome")
            need = _could(r, ctx)
            hint = hint or (f"{' or '.join(need)} would confirm" if need else None)
        if "not_one_per_chromosome" in sflags:
            out.append("the largest arrays of the centromeric satellites do not fall one per chromosome in this genome, so "
                       "they may not be its centromeres")
            need = _could(r, ctx)
            hint = hint or (f"{' or '.join(need)} would decide" if need else "CENP-A/CENH3 data would decide")
        if r.seq_evidence == "retrotransposon" and not ctx["decisive"]:
            out.append(f"placed by a {r.te_clade} core, in a genome where this lineage does not decide"
                       + (f" ({ctx['te_why'].rpartition('; ')[2]})" if ctx.get("te_why") else ""))
        elif r.seq_evidence == "retrotransposon" and "weak_core" in sflags:
            out.append(f"a {r.te_clade} core of too few elements to be significant alone: the best candidate, not a finding")
        elif r.seq_evidence == "retrotransposon" and r.confidence == "low":
            out.append(f"the {r.te_clade} core holds under half of the signature's elements on this chromosome")
        if "low_satellite" in sflags:
            out.append("little centromeric satellite")
        if "diffuse" in sflags:
            out.append("a broad block with little satellite")
        if r.seq_evidence == "long_unit_tandem":
            out.append("a long-unit tandem repeat, not a satellite family")
        if "near_end" in sflags:
            out.append("near a chromosome end")
        if "gene_rich" in sflags:
            out.append("gene-rich for a centromere")
        if "te_discordant" in sflags and _int(r.te_start) is not None and not te_in_alt:
            out.append(f"the {r.te_clade} core lies elsewhere ({interval(r.te_start, r.te_end)})")
    elif r.call_source == "composition":
        out.append("a GC trough says where the sequence is AT-rich, not where the kinetochore sits")
        hint = hint or ("CENP-A/CENH3 data would confirm" if ctx["hic"] else "Hi-C would test it")
    elif r.confidence != "high":
        sup = _toks(r.support)                                 # the evidence basis reports, so the two sentences agree
        here = [x for x, on in ((f"{size(ctx['sat_bp'](r.chrom, int(r.start) - 1, int(r.end)))} of centromeric satellite",
                                 "satellite" in sup), (f"the {r.te_clade} core", "te" in sup),
                                ("CHG methylation", "meth" in sup), ("an AT trough", "at" in sup)) if on]
        out.append(f"rests on the 3D lines; the other evidence at this place is {' and '.join(here)}" if here else
                   "rests on the 3D lines alone: no sequence or LTR evidence at this place")
        if "hic_denovo_only" in flags:
            hint = "CENP-A/CENH3 data would confirm"
    peak = (f" ({interval(r.hic_anch_peak_start, r.hic_anch_peak_end)})"
            if _int(getattr(r, "hic_anch_peak_start", None)) is not None else "")
    if "hic_alt_ambiguous" in flags:
        out.append(f"Hi-C's peak agrees with this array and with the competing cluster "
                   f"({interval(r.seq_alt_start, r.seq_alt_end)}) alike: no mappable bin separates them, so Hi-C "
                   f"cannot choose")
        hint = hint or "CENP-A/CENH3 data would decide"
    if "hic_offset" in flags:
        out.append(f"validated Hi-C points elsewhere{peak}, with no second line there")
        hint = hint or "an LTR core or CHG methylation at either place (plants), or CENP-A/CENH3 data, would decide"
    if "fold_offset" in flags:
        out.append(f"the arm fold points elsewhere ({interval(r.fold_core_start, r.fold_core_end)})")
    if "align_offset" in flags:
        out.append(f"the trans arm alignment points elsewhere ({interval(r.align_core_start, r.align_core_end)})")
    if "hic_denovo_offset" in flags:
        out.append(f"Hi-C's prior-free peak points elsewhere "
                   f"({interval(r.hic_denovo_peak_start, r.hic_denovo_peak_end)})")
    if "hic_denovo_elsewhere" in flags:
        out.append(f"Hi-C's prior-free peak, untested here, lies elsewhere "
                   f"({interval(r.hic_denovo_peak_start, r.hic_denovo_peak_end)})")
        hint = hint or "CENP-A/CENH3 data would decide"
    if "hic_disagrees_weak" in flags:
        out.append(f"Hi-C's peak lies elsewhere{peak} but is not strong enough to oppose")
    if "te" in against and not ("te_discordant" in sflags and r.call_source == "sequence") and not te_in_alt:
        out.append(f"the {r.te_clade} core lies elsewhere" + (f" ({interval(r.te_start, r.te_end)})"
                                                              if _int(r.te_start) is not None else ""))
    if "hic" in against and "hic_offset" not in flags and "hic_denovo_offset" not in flags and \
            not ("cis_elsewhere" in flags and str(getattr(r, "hic_gate", ".")) != "pass"):     # the cis line said it above
        out.append(f"validated Hi-C points elsewhere{peak}")
    if "at_elsewhere" in flags:
        out.append(f"the AT trough lies elsewhere ({interval(r.at_start, r.at_end)})")
    if "cis_elsewhere" in flags:
        out.append("the cis Hi-C signal has its extreme elsewhere" + (f" ({interval(r.cis_peak_start, r.cis_peak_end)})"
                                                                      if _int(getattr(r, "cis_peak_start", None)) else ""))
    if "cis_blind" in flags:
        out.append(f"the cis Hi-C signal, which passed its test on most chromosomes, has its extreme elsewhere "
                   f"({interval(r.cis_peak_start, r.cis_peak_end)}), but it is read outside the chromosome ends and cannot judge "
                   f"a call this close to one")
    if "cis_only" in flags:
        out.append("a cis Hi-C signal alone")
        hint = hint or "CENP-A/CENH3 data would decide"
    n = _int(getattr(r, "call_gaps", None))
    if n:
        out.append(f"{n} assembly gap{'s' if n > 1 else ''} inside the call")
    if not out:
        return "."
    return "; ".join(out) + (failed if r.confidence != "high" else "") + (f" -> {hint}" if hint else "")


def columns(t: pd.DataFrame, ctx: dict) -> pd.DataFrame:
    """t with basis and caveat appended."""
    t = t.copy()
    rows = list(t.itertuples(index=False))
    t["basis"] = [basis(r, ctx) for r in rows]
    t["caveat"] = [caveat(r, ctx) for r in rows]
    return t


# ------------------------------------------------------------------------------------------------ report
def _pval(x) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return 1.0


def _ends_sig(s: str) -> bool:
    """'2.14x P=0.000999' -> P <= 0.01 (the feature-class test's level)."""
    try:
        return float(str(s).split("P=")[1]) <= 0.01
    except (IndexError, ValueError):
        return False


def _n(x) -> int:
    v = _int(x)
    return v if v is not None else 0


def _list(items, k: int = 8) -> str:
    items = list(items)
    return ", ".join(items[:k]) + (f" and {len(items) - k} more" if len(items) > k else "")


def _md(x) -> str:
    return str(x).replace("|", "/")


def _count(n: int, one: str, many: str | None = None) -> str:
    return f"{n} {one if n == 1 else many or one + 's'}"


RANK = {"none": -1, "low": 0, "medium": 1, "high": 2}


def facts(t: pd.DataFrame, summ: dict, lt: pd.DataFrame, fc: pd.DataFrame | None, ctx: dict) -> dict:
    """The counts the report states, from the tables only (bench checker recomputes them from the written files)."""
    called = t[t.start.astype(str) != "."]
    conf = called.confidence.value_counts()
    seq = called[called.call_source == "sequence"]
    sat = seq[seq.seq_evidence.isin(["satellite", "weak_satellite"])]
    fams = sat.main_family.value_counts()
    sup = called.support.astype(str).str.split(",")
    against = called.against.astype(str).str.split(",")
    kinds = lt.kind.value_counts() if len(lt) else pd.Series(dtype=int)
    sat_cols = [c for c in lt.columns if c.startswith("sat_bp_")]
    blk = {k: (int(kinds.get(k, 0)), int(lt.length[lt.kind == k].sum()),
               int(lt.loc[lt.kind == k, sat_cols].apply(pd.to_numeric, errors="coerce").sum().sum()) if len(lt) else 0)
           for k in ("interstitial_block", "subtelomeric_block")}
    f = dict(chromosomes=len(t), called=len(called), high=int(conf.get("high", 0)), medium=int(conf.get("medium", 0)),
             low=int(conf.get("low", 0)), none=len(t) - len(called),
             satellite_calls=len(sat), families=[(fam, int(n), ctx["unit"].get(fam)) for fam, n in fams.items()],
             te_calls=int((seq.seq_evidence == "retrotransposon").sum()),
             at_calls=int((called.call_source == "composition").sum()),
             chain_calls=int(seq.seq_flags.astype(str).str.contains("partial_chain").sum()),
             tandem_calls=int((seq.seq_evidence == "long_unit_tandem").sum()),
             hic_moved=int((called.relation == "offset").sum()),
             hic_placed=int(((called.relation == "functional_only")
                             & ~called["flags"].astype(str).str.contains("hic_denovo_only")).sum()),
             hic_denovo_only=int(called["flags"].astype(str).str.contains("hic_denovo_only").sum()),
             hic_agree=int(sup.apply(lambda x: "hic" in x).sum()),
             hic_support=int(called["flags"].astype(str).str.contains("hic_support").sum()),
             hic_chose_alt=int((called["flags"].astype(str).str.contains("hic_support")
                                & called["flags"].astype(str).str.contains("hic_chose_alt")).sum()),
             hic_support_raised=int((called["flags"].astype(str).str.contains("hic_support")
                                     & (called.confidence.map(RANK) > called.seq_confidence.map(RANK))).sum()),
             hic_against=int(against.apply(lambda x: "hic" in x).sum()),
             te_agree=int(sup.apply(lambda x: "te" in x).sum()), meth_agree=int(sup.apply(lambda x: "meth" in x).sum()),
             contested=int((called.relation == "contested").sum()),
             kinetochores=int((called.cpg_kt_start.astype(str) != ".").sum()) if "cpg_kt_start" in called else 0,
             call_median=float(called.length.astype(int).median()) if len(called) else np.nan,
             blocks=blk, its_in_calls=sum(1 for c, s, e, _ in ctx.get("its", []) for r in called.itertuples()
                                          if r.chrom == c and s < int(r.end) and e > int(r.start) - 1))
    if "morphology" in t:
        k, kl = {}, {}
        for r in called.itertuples(index=False):
            d = k if r.confidence in ("high", "medium") else kl
            cls = levan(float(r.arm_ratio))
            d[cls] = d.get(cls, 0) + 1
        f.update(karyotype=k, karyotype_low=kl, t2t=int((t.t2t == "yes").sum()),
                 gapped_calls=[r.chrom for r in called.itertuples() if _n(r.call_gaps)],
                 short_end=[r.chrom for r in called.itertuples() if r.short_arm == "end" and r.confidence in
                            ("high", "medium") and levan(float(r.arm_ratio)) != "metacentric"],
                 telomere_both=int((t.telomeres == "both").sum()))
    if fc is not None and len(fc):                    # the centromeres always (the positive control), others when
        f["contacts"] = [(r.feature_class, float(r.fold), float(r.p_enriched)) for r in fc.itertuples(index=False)
                         if r.kind == "centromere" or (float(r.p_enriched) <= 0.01 and r.kind in
                                                       ("ends", "telomere", "nor", "interstitial", "subtelomeric",
                                                        "rDNA", "te_core"))]   # self-associated (P <= 0.01)
    return f


def genome_section(name: str, g, t: pd.DataFrame, summ: dict, lt: pd.DataFrame, fc, ctx: dict, lines: pd.DataFrame | None = None) -> str:
    """One genome's part of the report (Markdown)."""
    f = facts(t, summ, lt, fc, ctx)
    given = ["FASTA"] + [LINE[x] for x, ok in (("te", g.tes), ("hic", g.hic), ("meth", g.meth), ("genes", g.genes))
                         if ok]
    missing = [LINE[x] for x, ok in (("te", g.tes), ("hic", g.hic), ("meth", g.meth), ("genes", g.genes)) if not ok]
    L = [f"## {name}", "", f"**{f['chromosomes']} chromosomes. Inputs**: {', '.join(given)}" +
         (f"; not given: {', '.join(missing)}" if missing else "") + ".", ""]
    if g.notes:
        L += [f"**Notes**: {'; '.join(g.notes)}.", ""]
    L += ["### Centromeres", "",
          f"- **{f['called']}/{f['chromosomes']} chromosomes called**: {f['high']} high, {f['medium']} medium, "
          f"{f['low']} low" + (" (uncertain: the best candidate, another place not excluded)" if f["low"] else "") +
          (f"; {f['none']} without a call" if f["none"] else "") + "."]
    made = []
    if f["satellite_calls"]:
        fam = ", ".join(f"{x} ({u} bp) on {n}" if u else f"{x} on {n}" for x, n, u in f["families"][:6])
        made.append(f"satellite arrays on {f['satellite_calls']}: {fam}" +
                    (f" and {len(f['families']) - 6} more families" if len(f["families"]) > 6 else ""))
    if f["chain_calls"]:
        made.append(f"chains of partial copies of a longer repeat, not tandem satellites, on {f['chain_calls']} of those")
    if f["te_calls"]:
        made.append(f"retrotransposon cores on {f['te_calls']}")
    if f["at_calls"]:
        made.append(f"AT troughs on {f['at_calls']}")
    if f["tandem_calls"]:
        made.append(f"long-unit tandem arrays on {f['tandem_calls']}")
    if f["hic_moved"] or f["hic_placed"]:
        made.append(f"the validated 3D lines moved {f['hic_moved']} and placed {f['hic_placed']} calls")
    if f["hic_denovo_only"]:
        made.append(f"Hi-C's prior-free peaks alone placed {f['hic_denovo_only']} `low` calls (too few reference calls "
                    f"to test Hi-C)")
    if made:
        L.append(f"- **Basis**: {'; '.join(made)}.")
    if t.seq_flags.astype(str).str.contains("no_dominant_satellite").any():
        L.append("- **No satellite family dominates the genome**: satellite arrays alone cannot tell a centromere from "
                 "other blocks, so those calls are capped `low` unless another line supports them.")
    if g.hic:
        gate = _n(summ.get("hic_gate_pass"))
        n = _n(summ.get("hic_chromosomes"))
        ax = summ.get("hic_axis")                                         # fuse.genome_axis's decision
        if ax is None:                                                    # a summary written before 0.10.1
            ax = "fold" if summ.get("fold_gate") == "pass" else \
                "align" if summ.get("align_gate") == "pass" and _pval(summ.get("fold_p")) <= AXIS_LEAN_P else "."
        axis = {"fold": "arm fold", "align": "trans arm alignment"}.get(ax)
        lean = ax == "align"
        untested_all = bool(n) and _n(summ.get("hic_gate_untested")) == n
        if gate:
            L.append(f"- **Hi-C** passed its test (it finds the other chromosomes' centromeres) on {gate}/{n} "
                     f"chromosomes: it agrees with {f['hic_agree']} calls" +
                     (f" and points elsewhere on {f['hic_against']} (kept, contested)" if f["hic_against"] else "") +
                     (f"; it moved {f['hic_moved']} and placed {f['hic_placed']}" if f["hic_moved"] or f["hic_placed"]
                      else "") + ".")
        elif _n(summ.get("hic_gate_support")):
            L.append(f"- **Hi-C** found the other chromosomes' centromeres beyond chance but not on most "
                     f"({_n(summ.get('hic_gate_support'))}/{n} chromosomes passed the support-only test): it supports "
                     f"{f['hic_support']} calls where its peak agrees ({f['hic_support_raised']} raised one level)" +
                     (f", {f['hic_chose_alt']} of them on the competing satellite cluster its strong peak chose" if f["hic_chose_alt"]
                      else "") + "; on this test it moves no call.")
        elif str(summ.get("hic_gate_pass", ".")) != "." and axis:
            L.append(f"- **Hi-C** " + (f"could not be tested against reference calls (too few) on any chromosome"
                                         if untested_all else f"failed its anchored test on every chromosome ({n})") +
                     f", but its {axis} passed: the "
                     f"two 3D lines moved {f['hic_moved']} and placed {f['hic_placed']} calls where they agree" +
                     (f"; Hi-C's prior-free peaks alone placed {f['hic_denovo_only']} `low` calls where they do not"
                      if f["hic_denovo_only"] else "") + ".")
        elif (f["hic_moved"] or f["hic_placed"]) and (summ.get("te_hic", ".") not in (".", None) or
                                                       summ.get("cis_line", ".") not in (".", None)):
            L.append(f"- **Hi-C** " + (f"could not be tested against reference calls (too few)" if untested_all else
                                         f"failed its anchored test on every chromosome ({n})") +
                     f", but its readings tested otherwise (below: the LTR cores, or its cis signal) moved {f['hic_moved']} "
                     f"and placed {f['hic_placed']} calls.")
        elif n and _n(summ.get("hic_gate_untested")) == n and not axis:
            L.append(f"- **Hi-C** could not be tested: too few reference calls to test it against" +
                     (f"; its prior-free peaks alone placed {f['hic_denovo_only']} `low` calls on chromosomes without "
                      f"one" if f["hic_denovo_only"] else "; it changes no call") + ".")
        elif str(summ.get("hic_gate_pass", ".")) != ".":
            L.append(f"- **Hi-C** failed its test on every chromosome ({n}): it reports but changes no call.")
        else:
            L.append("- **Hi-C**: no library could decide (see the log).")
        if str(summ.get("hic_sequence_conflict")) == "True":
            L.append("- **Warning**: Hi-C's own peaks and the sequence calls disagree genome-wide; one line is wrong "
                     "everywhere" + (f"; the {axis} sided with Hi-C" if axis else
                                     "; tes/ or methylation would decide") + ".")
        if summ.get("fold_gate") == "pass":
            L.append(f"- **Arm fold** (Rabl nuclei) passed its test ({summ.get('fold_on_denovo')}).")
        arch = {"rabl_axis": "a Rabl axis (centromeres clustered at one pole, arms aligned; Type I nucleus)",
                "centromere_clustering": "centromeres that meet in trans on most chromosomes",
                "weak_centromere_clustering": "centromeres that meet in trans beyond chance but not on most "
                                              "chromosomes (territorial, Type II-like nucleus)",
                "untested": "no test was possible (too few reference calls, no Rabl axis)",
                "not_detected": "no centromere signature (a territorial or holocentric nucleus, or sequence calls "
                                "wrong genome-wide: Hi-C alone cannot tell which)"}.get(summ.get("hic_architecture"))
        cis_l = str(summ.get("cis_line", "."))
        if summ.get("hic_architecture") in ("not_detected", "untested") and ("tested by trans contacts" in cis_l
                                                                             and ": pass," in cis_l):
            arch = ("centromeres that meet in trans when anchored on Hi-C's own cis extremes (the sequence calls failing as "
                    "anchors): see the contact crossing each bin, below")
        if arch:
            ends = summ.get("hic_3d_ends", ".")
            L.append(f"- **3D organization** of this library: {arch}" +
                     (f"; chromosome ends meet in trans ({ends})" if ends not in (".", None) and
                      _ends_sig(ends) else "") + ".")
        if summ.get("align_gate") == "pass" and lean:
            L.append(f"- **Trans arm alignment** (Rabl nuclei: loci at equal distances from their centromeres meet) "
                     f"passed its test ({summ.get('align_on_denovo')})" +
                     (", the arm fold having failed: it acts in the fold's place" if summ.get("fold_gate") != "pass"
                      else "; the arm fold, which passed too, decides") + ".")
    if summ.get("sat_one_per_chromosome", ".") not in (".", None):
        L.append(f"- **Do the centromeric satellites fall one per chromosome?** {summ.get('sat_one_per_chromosome')}."
                 + (" Their calls are capped at `low` unless an LTR core agrees." if
                    str(summ.get("sat_one_per_chromosome")).endswith("fail") else ""))
    if g.tes:
        m = summ.get("te_marker", ".")
        if m not in (".", None):
            seed, sp = summ.get("te_seed", "."), summ.get("te_specificity", ".")
            why = str(summ.get("te_verdict", ".")).rpartition("; ")[2]        # 'flags only (...)'
            L.append(f"- **LTR retrotransposons**: {m} marks centromeres ("
                     + ("decisive" if str(summ.get("te_decisive")) == "True" else why if "flags only" in why else f"flags only: {why}") + ")"
                     + (f": signature seeded by {seed}, {sp}-fold in its cores" if seed not in (".", None) else "")
                     + f"; its core agrees with {f['te_agree']} calls, and {f['te_calls']} calls stand on it.")
        else:
            L.append("- **LTR retrotransposons**: no lineage or family concentrates at one locus per chromosome.")
    if "tested by trans contacts" in str(summ.get("cis_line", ".")):
        L.append(f"- **The contact crossing each bin (cis)**, tested by Hi-C's trans contacts (no LTR cores; the anchored test "
                 f"failed or could not run): {summ.get('cis_line')}.")
    if summ.get("te_hic", ".") not in (".", None):
        L.append(f"- **Hi-C on the LTR cores**: anchored on them (the sequence calls having failed as anchors), Hi-C's peak lies at "
                 f"the core on {summ.get('te_hic')}. The contact crossing each bin (cis): {summ.get('cis_line', '.')}.")
    if summ.get("at_gate", ".") not in (".", None):
        hic_at = str(summ.get("at_hic", "."))
        L.append(f"- **AT troughs** (most chromosomes have no sequence call): the genome's deepest GC troughs "
                 + ("fall one per chromosome" if summ.get("at_gate") == "pass" else
                    "do not fall one per chromosome" if summ.get("at_gate") == "fail" else
                    "cannot be tested from the FASTA (too few chromosomes)")
                 + f" (P {summ.get('at_p')}, windows of {summ.get('at_window')} bp)"
                 + (f"; against Hi-C's prior-free peaks: {hic_at}" if hic_at != "." else "") + ".")
    if g.meth:
        if "meth_signature" in summ:
            sig = str(summ.get("meth_signature"))
            L.append("- **CHG methylation**: " + ("no centromere signature found: the line does not act." if
                                                  sig in ("none", ".", "") else
                                                  f"signature {sig}; supports {f['meth_agree']} calls."))
        if "cpg_gate_pass" in summ:
            gate = _n(summ.get("cpg_gate_pass"))
            L.append(f"- **CpG methylation**: " + (f"its dips mark the other chromosomes' centromeres on {gate} "
                                                   f"chromosomes (its test); kinetochores inside {f['kinetochores']} "
                                                   f"calls." if gate else
                                                   "failed its test (dips do not mark the other chromosomes' "
                                                   "centromeres beyond chance): no kinetochore placed."))
    tested = lines[lines.chromosomes >= 3] if lines is not None and len(lines) else None
    if tested is not None and len(tested) and len(tested) <= 28:
        L.append("- **Do the lines agree with each other?** (`line_tests.tsv`) " + "; ".join(
            f"{r.line_a.replace('_', ' ')} and {r.line_b.replace('_', ' ')} {r.agree}/{r.chromosomes} (expected {r.expected:g})"
            for r in tested[tested.verdict.str.startswith("agree")].itertuples()) + "." if
                 tested.verdict.str.startswith("agree").any() else
                 "- **Do the lines agree with each other?** (`line_tests.tsv`) no two lines agree beyond chance.")
    other = [x for x in str(summ.get("inputs_other", ".")).split(",") if x not in (".", "")]
    if other:
        L.append(f"- **Inputs from other plants** (they score their plant and do not decide the calls): "
                 f"{_list(x.split(':', 1)[1] + ' (' + x.split(':', 1)[0] + ')' for x in other)}.")
    if _n(summ.get("plants")):
        moved = {v: [x for x in str(summ.get(f"plants_{v}", ".")).split(",") if x not in (".", "")]
                 for v in ("shifted", "possible_shift", "call_absent", "conflict")}
        L.append(f"- **Several plants**: {_n(summ.get('plants'))} scored on their own; " +
                 ("; ".join(f"{v.replace('_', ' ')}: {_list(x.replace(':', ' ') for x in xs)}"
                            for v, xs in moved.items() if xs) or "every centromere shared with the assembled plant")
                 + " (`<genome>.plants.tsv`).")
    unc = t[(t.confidence != "high") | (t.relation == "contested")]
    if len(unc):
        L += ["", "**Calls that are not `high`, and why:**", "", "| chromosome | call | confidence | caveat |",
              "|---|---|---|---|"]
        for r in unc.itertuples(index=False):
            call = interval(r.start, r.end) if str(r.start) != "." else "none"
            L.append(f"| {_md(r.chrom)} | {call} | {r.confidence} | {_md(r.caveat)} |")
    if "karyotype" in f:
        k = f["karyotype"]
        order = ("metacentric", "submetacentric", "subtelocentric", "acrocentric")
        kl = f["karyotype_low"]
        L += ["", "### Chromosomes", "",
              "- **Morphology** (arm ratio by sequence, Levan 1964), confident calls: " +
              (", ".join(f"{k[c]} {c}" for c in order if c in k) or "none") +
              (f"; low calls: {', '.join(f'{kl[c]} {c}' for c in order if c in kl)}" if kl else "") + "."]
        L.append(f"- **Assembly**: {f['t2t']}/{f['chromosomes']} chromosomes telomere-to-telomere (telomere arrays at "
                 f"both ends, no gap); {f['telomere_both']} with telomere arrays at both ends; " +
                 (f"calls holding assembly gaps: {_list(f['gapped_calls'])}" if f["gapped_calls"] else
                  "no call holds an assembly gap") + ".")
        if f["short_end"]:
            L.append(f"- **Short arm at the sequence end** (confident calls, not metacentric; the cytogenetic "
                     f"convention puts the short arm first): {_list(f['short_end'])}.")
    L += ["", "### Other features", ""]
    unit = summ.get("telomere_unit", ".")
    if unit not in (".", None):
        its = _n(summ.get("telomere_interstitial"))
        L.append(f"- **Telomeres**: unit {unit}, arrays at {_n(summ.get('telomere_ends'))}/{2 * f['chromosomes']} "
                 f"ends; " + (_count(its, "interstitial telomeric array") +
                              (f", {f['its_in_calls']} inside centromere calls" if f["its_in_calls"] else "")
                              if its else "no interstitial telomeric array") + ".")
    else:
        L.append("- **Telomeres**: no telomere unit found at the chromosome ends.")
    nor, five = _n(summ.get("rdna_arrays")), _n(summ.get("rdna5s_arrays"))
    L.append("- **rDNA**: " + (f"{_count(nor, '45S array')} ({_count(nor, 'NOR').split(' ', 1)[1]}) on "
                               f"{_md(summ.get('rdna_chromosomes'))}" if nor else "no 45S array") + "; " +
             (f"{_count(five, '5S array')} on {_md(summ.get('rdna5s_chromosomes'))}" if five else "no 5S array") + ".")
    blk = f["blocks"]
    kinds = (("interstitial_block", "interstitial", " (knobs, pericentromeric satellites and others)"),
             ("subtelomeric_block", "subtelomeric", ""))
    blocks = [f"{blk[k][0]} {label}, {size(blk[k][1])}, {blk[k][2] / blk[k][1]:.0%} satellite{note}"
              for k, label, note in kinds if blk[k][0]]
    L.append("- **Satellite blocks**: " + ("; ".join(blocks) if blocks else "none outside the calls") + ".")
    if f.get("contacts"):
        L.append("- **3D** (Hi-C contacts among a class's instances on different chromosomes, vs random places): " +
                 "; ".join(f"{c} {x:.1f}x (P {p:.2g})" for c, x, p in f["contacts"]) + ".")
    L += ["", "### Per chromosome", "",
          "| chromosome | length | call | confidence | morphology (CI) | telomeres | basis |",
          "|---|---|---|---|---|---|---|"]
    for r in t.itertuples(index=False):
        call = interval(r.start, r.end) if str(r.start) != "." else "none"
        mo = f"{r.morphology} ({r.ci})" if "morphology" in t and str(r.morphology) != "." else "."
        tel = getattr(r, "telomeres", ".")
        L.append(f"| {_md(r.chrom)} | {size(r.chrom_len)} | {call} | {r.confidence} | {mo} | {tel} | "
                 f"{_md(r.basis)} |")
    return "\n".join(L) + "\n"


def write_report(sections: list, out: Path, genome_dir: Path) -> Path:
    """cenfuse_report.md: one section per genome."""
    head = (f"# cenfuse report\n\ncenfuse {__version__}, {date.today()}, `{genome_dir}`. One row per chromosome: "
            f"`<genome>.centromeres.tsv` (its `basis` and `caveat` columns are the sentences below); every feature: "
            f"`<genome>.loci.tsv`; figure: `cenfuse_overview.pdf`. Confidence: `high`, `medium`, `low` = "
            f"uncertain.\n\n")
    p = Path(out) / "cenfuse_report.md"
    p.write_text(head + "\n".join(sections))
    return p
