"""Interpretation (README, 'Reading a call'): what each call stands on and why it may be wrong, in words, and the run
report (cenfuse_report.md). Both are read from the tables cenfuse writes, so they say nothing the tables do not; nothing
here changes a call.

- basis: the sequence evidence (satellite family and unit, its bp in the call; a retrotransposon core), a move or
  placement by validated 3D lines and its second line, each supporting line, the CpG kinetochore.
- caveat: each confidence-lowering flag and opposing line, gaps in the call, and the input that would decide where
  cenfuse's own tests say which one can (CENP-A/CENH3 data for two arrays of one satellite: neither Hi-C nor CpG
  methylation chooses between them when no mappable bins separate them)."""
from __future__ import annotations

from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__
from .morphology import levan
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
    return dict(unit=unit, sat_bp=sat_bp, marker=m[0] if m else None, decisive=bool(m and m[2]), hic=hic_given,
                tes=tes_given, its=list((telo or {}).get("interstitial", [])))


def basis(r, ctx: dict) -> str:
    """What the call stands on, in words ('.' without a call)."""
    if r.relation == "none" or _int(r.start) is None:
        return "."
    flags, sflags, sup = set(_toks(r.flags)), set(_toks(r.seq_flags)), _toks(r.support)
    s, e = int(r.start) - 1, int(r.end)
    out = []
    if r.call_source == "sequence":
        ev = r.seq_evidence
        if ev in ("satellite", "weak_satellite"):
            unit = ctx["unit"].get(r.main_family)
            out.append(f"{'weak ' if ev == 'weak_satellite' else ''}satellite {r.main_family}"
                       + (f" ({unit}-bp unit)" if unit else "") + f", {size(ctx['sat_bp'](r.chrom, s, e))} of "
                       f"centromeric satellite in the call")
        elif ev == "long_unit_tandem":
            out.append(f"long-unit tandem array, {size(ctx['sat_bp'](r.chrom, s, e))} of satellite in the call")
        elif ev == "retrotransposon":
            replaced = sflags & {"weak_satellite_beside_core", "satellite_away_from_core", "satellite_beside_core"}
            out.append(f"{r.te_clade} retrotransposon core" + ("" if replaced else " (no satellite call)"))
        if "weak_satellite_beside_core" in sflags:
            out.append("the CRM core replaced a weak satellite array beside it")
        if "satellite_away_from_core" in sflags:
            out.append("the CRM core replaced a satellite array elsewhere")
        if "satellite_beside_core" in sflags:
            out.append("the CRM core replaced a satellite array beside it that holds significantly fewer CRM elements")
        if "hic_chose_alt" in flags:
            out.append("the second satellite cluster, chosen by validated Hi-C")
        if "acrocentric" in sflags:
            out.append("the family sits at one end of most chromosomes (acrocentric karyotype)")
    else:
        if "fold_pair" in flags:
            how = "the arm fold and Hi-C's prior-free peak, which agree here"
        elif "hic_medium_peak" in flags:
            how = "a clean medium validated Hi-C peak"
        else:
            how = "a strong validated Hi-C peak"
        if r.relation == "offset":
            off = _int(r.offset_bp)
            out.append(f"moved {size(off) + ' ' if off is not None else ''}from the sequence call "
                       f"({interval(r.seq_start, r.seq_end)}) by {how}")
        else:
            out.append(f"placed by {how} (no sequence call)")
        second = {"te": f"the {r.te_clade} core at the peak", "meth": "CHG methylation at the peak",
                  "fold": "the arm fold"}
        named = set()
        for f in sorted(flags):
            if f.startswith("moved_with_"):
                named |= set(f[len("moved_with_"):].split("+"))
                out += [f"second line: {second[x]}" for x in f[len("moved_with_"):].split("+") if x in second]
        if "te" in sup and "te" not in named:            # a call below high moved with the CRM core at the peak
            out.append(f"the {r.te_clade} core at the peak agrees")
        if "satellite" in sup:
            out.append(f"{size(ctx['sat_bp'](r.chrom, s, e))} of centromeric satellite inside")
    for line in sup:
        if line == "te" and r.call_source == "sequence" and r.seq_evidence != "retrotransposon":
            out.append(f"the {r.te_clade} core agrees")
        elif line == "hic" and r.call_source == "sequence":
            out.append("validated Hi-C agrees")
        elif line == "meth":
            out.append("CHG methylation is centromere-like")
        elif line == "fold" and "fold_pair" not in flags:
            out.append("the arm fold agrees")
    ks, ke = _int(getattr(r, "cpg_kt_start", None)), _int(getattr(r, "cpg_kt_end", None))
    if ks is not None and ke is not None:
        out.append(f"CpG kinetochore {interval(ks, ke)} ({size(ke - ks + 1)})")
    return "; ".join(out) or "."


def _could(r, ctx: dict) -> list:
    """Inputs not given that could confirm or move a call: Hi-C, a CRM core (tes/)."""
    out = []
    if not ctx["hic"]:
        out.append("Hi-C")
    if not ctx["tes"]:
        out.append("a CRM core (tes/; plants)")
    return out


def caveat(r, ctx: dict) -> str:
    """Why the call may be wrong, and the input that would decide ('.' when nothing is known against it)."""
    failed = "; Hi-C failed its test here" if ctx["hic"] and str(getattr(r, "hic_gate", ".")) == "fail" else ""
    if r.relation == "none" or _int(r.start) is None:
        need = _could(r, ctx)
        return ("no line found a centromere (satellite-free?)" + failed +
                (f" -> {' or '.join(need)} can place one" if need else ""))
    flags, sflags, against = set(_toks(r.flags)), set(_toks(r.seq_flags)), _toks(r.against)
    out, hint = [], None
    te_iv = (_int(r.te_start), _int(r.te_end))
    alt_iv = (_int(r.seq_alt_start), _int(r.seq_alt_end))
    te_in_alt = None not in te_iv + alt_iv and te_iv[0] <= alt_iv[1] and te_iv[1] >= alt_iv[0]
    if r.call_source == "sequence" and "competing_cluster" in sflags and "hic_chose_alt" not in flags:
        alt = f" ({interval(*alt_iv)})" if alt_iv[0] is not None else ""
        where = f", where the {r.te_clade} core lies" if te_in_alt else ""
        out.append(f"another satellite cluster{alt} competes{where}: the sequence cannot tell which is active")
        hint = "CENP-A/CENH3 data would decide" + ("" if ctx["tes"] else ", or a CRM core (tes/; plants)") + \
            ("" if ctx["hic"] else ", or Hi-C when mappable DNA separates them")
    if r.call_source == "sequence":
        if "no_dominant_satellite" in sflags:
            out.append("no satellite family dominates the genome")
            need = _could(r, ctx)
            hint = hint or (f"{' or '.join(need)} would confirm" if need else None)
        if r.seq_evidence == "retrotransposon" and not ctx["decisive"]:
            out.append(f"placed by a {r.te_clade} core: only CRM cores may decide, other clades flag")
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
    elif r.confidence != "high":
        out.append("rests on the 3D lines alone: no sequence or LTR evidence at this place")
    peak = (f" ({interval(r.hic_anch_peak_start, r.hic_anch_peak_end)})"
            if _int(getattr(r, "hic_anch_peak_start", None)) is not None else "")
    if "hic_offset" in flags:
        out.append(f"validated Hi-C points elsewhere{peak}, with no second line there")
        hint = hint or "a CRM core or CHG methylation at either place (plants), or CENP-A/CENH3 data, would decide"
    if "fold_offset" in flags:
        out.append(f"the arm fold points elsewhere ({interval(r.fold_core_start, r.fold_core_end)})")
    if "hic_denovo_offset" in flags:
        out.append(f"Hi-C's prior-free peak points elsewhere "
                   f"({interval(r.hic_denovo_peak_start, r.hic_denovo_peak_end)})")
    if "hic_disagrees_weak" in flags:
        out.append(f"Hi-C's peak lies elsewhere{peak} but is not strong enough to oppose")
    if "te" in against and not ("te_discordant" in sflags and r.call_source == "sequence") and not te_in_alt:
        out.append(f"the {r.te_clade} core lies elsewhere" + (f" ({interval(r.te_start, r.te_end)})"
                                                              if _int(r.te_start) is not None else ""))
    if "hic" in against and "hic_offset" not in flags and "hic_denovo_offset" not in flags:
        out.append(f"validated Hi-C points elsewhere{peak}")
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
             tandem_calls=int((seq.seq_evidence == "long_unit_tandem").sum()),
             hic_moved=int((called.relation == "offset").sum()),
             hic_placed=int((called.relation == "functional_only").sum()),
             hic_agree=int(sup.apply(lambda x: "hic" in x).sum()),
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


def genome_section(name: str, g, t: pd.DataFrame, summ: dict, lt: pd.DataFrame, fc, ctx: dict) -> str:
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
    if f["te_calls"]:
        made.append(f"retrotransposon cores on {f['te_calls']}")
    if f["tandem_calls"]:
        made.append(f"long-unit tandem arrays on {f['tandem_calls']}")
    if f["hic_moved"] or f["hic_placed"]:
        made.append(f"the validated 3D lines moved {f['hic_moved']} and placed {f['hic_placed']} calls")
    if made:
        L.append(f"- **Basis**: {'; '.join(made)}.")
    if t.seq_flags.astype(str).str.contains("no_dominant_satellite").any():
        L.append("- **No satellite family dominates the genome**: satellite arrays alone cannot tell a centromere from "
                 "other blocks, so those calls are capped `low` unless another line supports them.")
    if g.hic:
        gate = _n(summ.get("hic_gate_pass"))
        n = _n(summ.get("hic_chromosomes"))
        if gate:
            L.append(f"- **Hi-C** passed its test (it finds the other chromosomes' centromeres) on {gate}/{n} "
                     f"chromosomes: it agrees with {f['hic_agree']} calls" +
                     (f" and points elsewhere on {f['hic_against']} (kept, contested)" if f["hic_against"] else "") +
                     (f"; it moved {f['hic_moved']} and placed {f['hic_placed']}" if f["hic_moved"] or f["hic_placed"]
                      else "") + ".")
        elif str(summ.get("hic_gate_pass", ".")) != "." and summ.get("fold_gate") == "pass":
            L.append(f"- **Hi-C** failed its anchored test on every chromosome ({n}), but its arm fold passed: the "
                     f"two 3D lines moved {f['hic_moved']} and placed {f['hic_placed']} calls where they agree.")
        elif str(summ.get("hic_gate_pass", ".")) != ".":
            L.append(f"- **Hi-C** failed its test on every chromosome ({n}): it reports but changes no call.")
        else:
            L.append("- **Hi-C**: no library could decide (see the log).")
        if str(summ.get("hic_sequence_conflict")) == "True":
            L.append("- **Warning**: Hi-C's own peaks and the sequence calls disagree genome-wide; one line is wrong "
                     "everywhere" + ("; the arm fold sided with Hi-C" if summ.get("fold_gate") == "pass" else
                                     "; tes/ or methylation would decide") + ".")
        if summ.get("fold_gate") == "pass":
            L.append(f"- **Arm fold** (Rabl nuclei) passed its test ({summ.get('fold_on_denovo')}).")
    if g.tes:
        m = summ.get("te_marker", ".")
        if m not in (".", None):
            L.append(f"- **LTR retrotransposons**: {m} marks centromeres "
                     f"({'decisive' if str(summ.get('te_decisive')) == 'True' else 'flags only: only CRM may decide'}"
                     f"); its core agrees with {f['te_agree']} calls.")
        else:
            L.append("- **LTR retrotransposons**: no clade concentrates at one locus per chromosome.")
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
