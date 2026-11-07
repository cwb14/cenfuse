"""Integration rules (docs/design.md, 'Fusion rules').

The integrated call is the sequence call unless a validated functional line moves it:
- Hi-C may act on chromosome x only when its gate passes (hicfuse.gate), and only through a strong anchored peak
  (all three Hi-C engine tests, no hub/rDNA/telomeric flag). A peak supports a call within HIC_TOL_BINS, only when it is a
  call at all (Hi-C level not 'none').
- A strong peak that points to another locus (more than HIC_FAR_BINS from the sequence call, the method's resolution,
  and with a core that does not reach the call) replaces a low or medium sequence call, and places a call on a
  chromosome without one. Closer than that, Hi-C neither supports nor opposes the call (flag hic_near).
- It replaces a high sequence call only with a second line at the Hi-C peak and not at the sequence call: a decisive
  (CRM) TE core within HIC_FAR_BINS of the peak and beyond it from the call, or methylation `centromere_like` at the
  peak with a measured non-centromere-like `dip` at the call (no data is not evidence), or (0.13, ageline) copies of the
  chromosome's signature families significantly younger at the peak than at the call (moved_with_age: the centromere
  moved and left its satellite and old copies behind). Without one the call stays and is `contested`.
- Support tier (hicfuse.gate): where the gate fails only for want of a majority but agreements on the other
  chromosomes are beyond chance (>= GATE_MIN_REF of them at significant peaks, Hi-C level medium or high; P <=
  GATE_P), Hi-C may only support: a kept sequence call whose own significant anchored peak agrees goes up one level
  (relation concordant, flag hic_support). On this tier Hi-C never opposes, moves or places a call; since 0.13 a strong
  peak on the competing cluster chooses it, as at pass (the sequence line's own second candidate; flags hic_chose_alt,
  hic_support). Where a Rabl-axis line passed (the arm fold, or the trans arm alignment), its rules govern Hi-C
  instead.
- Hi-C alone (0.10): in a genome with fewer than GATE_MIN_REF reference calls the gate cannot be tested on any
  chromosome (gate_genome_untested, pipeline.evidences), so nothing can validate Hi-C (point, regional and AT-rich
  centromeres: yeasts, Plasmodium, Neurospora). There, on a
  chromosome without a sequence call where no Rabl-axis line pairs, a clean prior-free Hi-C call (engine level high;
  no hub, rDNA, telomeric, not_converged or cycle_unstable flag) gives a low functional_only call at its core (flag
  hic_denovo_only): the best candidate, untested. Never where a sequence call or a test exists (prior-free high
  calls were right on 45% of T2T mammal chromosomes; with 10-kb bins on 16/16 yeast and 14/14 Plasmodium ones).
  Where such a genome has a sequence call and the clean prior-free peak lies beyond HIC_FAR_BINS of it, the call is
  kept and only flagged (hic_denovo_elsewhere; its caveat says where): nothing validates either line there.
- A call with a competing cluster whose anchored peak agrees with both (no mappable bin between them) gets no Hi-C
  support and no raise (flag hic_alt_ambiguous; 0.10.1): the agreement cannot tell the two apart (mouse: the major-
  and minor-satellite arrays at each chromosome's start), as hic_chose_alt needs mappable bins to choose.
- Methylation never moves a call alone.
- Competing clusters (the satellite engine's second cluster, at least half the main one): a strong validated peak
  that agrees with the second cluster and is separated from the main one by mappable bins chooses it (flag
  hic_chose_alt). Choosing moves the call, so it needs a strong peak like any move (CHM13 chr21: a medium p-arm peak
  chose the inactive array); Hi-C chooses between two sequence candidates and does not place the call.
- No sequence call, gate passed, a clean (no hub/rDNA/telomeric) anchored peak at Hi-C level medium: a low-confidence
  functional_only call (flag hic_medium_peak). Horse chr11, the satellite-free centromere, was found this way and
  discarded by the rule that asked for level high.
- Arm fold (fold.py), only where the anchored gate has not validated Hi-C on the chromosome and the fold passed its
  own test (fold cores hold the de novo Hi-C peaks genome-wide): the fold and the de novo peak are then two validated
  3D lines. Both at one locus and the sequence call at neither (beyond HIC_FAR_BINS mappable bins of the fold core and
  of the de novo peak): the call moves to the de novo core, medium (flag fold_pair; moved_with_fold for a high call,
  the fold being the second line at the peak and not at the call); no sequence call: a functional_only call there.
  Otherwise each validated 3D line that opposes the kept sequence call lowers it one level (contested; flags
  fold_offset, hic_denovo_offset). The fold is arm-scale (core ~5-10% of the chromosome): it never raises a call.
  Barley: satellite calls at the chromosome ends, the fold and de novo Hi-C at the CENH3 domains.
- Rabl axis (0.10.1): the trans arm alignment (fold.py) is a second readout of the same axis; where the fold's test
  fails but leans the same way (P <= AXIS_LEAN_P) and the alignment's test passes, the alignment takes the fold's
  place in every rule above (flags align_pair, moved_with_align, align_offset; support/against 'align'). Yeast and
  Plasmodium: the fold fails in these small nuclei, the alignment holds the prior-free peak on 13/16 and 12/14
  chromosomes. Hi-C alone then gives medium calls where the two agree, and keeps its low call where they do not
  (0.10: no call where the fold had passed). Otherwise the alignment is reported only (axis_line).

- LTR cores and Hi-C (0.12, hicfuse.anchored_on_cores): where the gate fails on most chromosomes because the sequence calls are
  wrong anchors, Hi-C is scored against the LTR signature's cores, weak ones too. If its peaks lie at the cores beyond chance
  (status support) or on most chromosomes as well (pass), a chromosome whose peak lies at its core holds two lines at one locus
  (flag te_hic_pair): a call there is supported, a low call elsewhere moves there (medium), with status pass any call does (a
  high one stays high, the core being its second line), and a chromosome without a call takes the core (medium).
- The cis signal (0.12, hicfuse.cis_line): where the gate fails on most chromosomes, the extreme of the contact crossing each bin,
  its sign and its level tested against the LTR cores of the other chromosomes (status support or pass). At a call it adds a level
  (cis_agrees); within Hi-C's resolution (HIC_FAR_BINS) it neither supports nor opposes (cis_near); beyond, a low call moves to it
  and stays low (cis_only), as does a chromosome without a call; a medium or high call stays, flagged (cis_elsewhere), one level
  down at status pass: one tested line against another, neither more sure.
- Which end (0.12, hicfuse.end_test): a call in the terminal eighth of a chromosome, where Hi-C's peaks cannot count because ends
  cluster in 3D anyway, goes up one level when the anchored score is higher at its end than at the other, in a genome where that
  holds on most of the other chromosomes with such a call (flag hic_end). Only where neither the gate, its support tier nor a
  Rabl-axis line speaks; never with a competing cluster at the same end (hic_alt_ambiguous).
- The cis signal tested by trans contacts (0.13, hicfuse.cis_trans_line): without LTR cores, the same rules where that test passes
  (pass only), read outside the terminal eighths: a call inside one is neither supported, opposed nor moved (blind); where the
  extreme lies beyond HIC_FAR_BINS of it, the call is flagged cis_blind and its caveat says so.
- AT troughs (0.12, atline.py): only on a chromosome without a sequence call. Where the genome's troughs passed the
  one-per-chromosome test, a chromosome nothing else places takes its trough, low (flag at_only). Where they also agree with
  Hi-C's prior-free peaks genome-wide (validated), a Hi-C call in agreement takes the trough's interval and rises a level
  (at_pair: low -> medium, medium -> high), one in disagreement is flagged (at_elsewhere), and a chromosome nothing else
  places takes its trough at medium. Agreement beyond chance but on a minority (atline.validate: support) allows at_pair,
  and (0.13) a chromosome nothing else places takes its trough, low, where a clean prior-free peak that placed no call of its
  own agrees with it (at_denovo_pair). A trough never moves or opposes a satellite or LTR call.

Every distance that involves a Hi-C interval is counted in mappable bins (hicfuse.visible_gap): Hi-C cannot see into
unmappable arrays, so a peak at an array's edge is adjacent to a call inside it. Agreement discounts at most
HIC_FAR_BINS unmappable bins; opposition needs that many mappable bins. Across a longer blind stretch Hi-C cannot
place the centromere on either side, so it neither supports nor opposes.

Confidence uses the satellite engine's levels, changed one step at a time (its te_concordant convention):
- A kept sequence call goes up one level for agreeing validated Hi-C and up one for methylation support, at most
  MAX_RAISE in total. Contested: down one level (one per opposing 3D line under the fold rule), offset by
  methylation support.
- A moved or functional-only call is medium, and high with a second line."""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Callable

import numpy as np
import pandas as pd

from .hicfuse import END_ZONE, SUPPORT_LEVELS, clean_peak, visible_gap
from .ageline import AGE_ALPHA
from .params import AXIS_LEAN_P, GATE_MIN_REF, HIC_FAR_BINS, HIC_TOL_BINS, MAX_RAISE, SAT_SUPPORT_BP
from .sat.te import agrees
from .util import gap, level_shift

FUSED_COLS = ["chrom", "chrom_len", "start", "end", "length", "call_source", "relation", "confidence", "support",
              "against", "flags", "offset_bp", "seq_start", "seq_end", "seq_evidence", "seq_confidence", "seq_flags",
              "seq_alt_start", "seq_alt_end", "seq_overruled_start", "seq_overruled_end", "seq_rejected_family",
              "seq_rejected_start", "seq_rejected_end", "seq_spacing", "main_family", "te_clade", "te_seed",
              "te_start", "te_end", "te_share", "te_chrom_elements", "te_core_elements", "te_call_elements", "at_gate", "at_start", "at_end",
              "at_gc", "hic_gate", "hic_gate_k", "hic_gate_n",
              "hic_gate_p", "hic_gate_refs", "hic_end_k", "hic_end_n",
              "hic_anch_peak_start", "hic_anch_peak_end", "hic_anch_core_start", "hic_anch_core_end", "hic_anch_level",
              "hic_anch_flags", "hic_anch_agrees_seq", "hic_denovo_peak_start", "hic_denovo_peak_end",
              "hic_denovo_confidence", "hic_denovo_flags", "fold_gate", "fold_point", "fold_core_start",
              "fold_core_end", "align_gate", "align_point", "align_core_start", "align_core_end", "cis_peak_start", "cis_peak_end", "cis_test", "age_peak", "age_call", "age_p", "age_premise_p", "meth_seq",
              "meth_seq_score", "meth_hic", "meth_hic_score",
              "meth_call", "meth_call_score"]


@dataclass
class Evidence:
    """Everything known about one chromosome (0-based half-open intervals; -1 = absent)."""
    chrom: str
    chrom_len: int
    seq: dict                 # start, end, evidence, confidence, flags, main_family; alt_start/alt_end (competing
                              # cluster, -1 = none); overruled_start/overruled_end (the call a TE core replaced: with
                              # a decisive CRM core a satellite, weak-satellite or long-unit call, with a core of any
                              # clade a long-unit call away from it; -1 = none): reported only, never a candidate for
                              # the call
    te: dict | None           # clade, start, end, decisive (the marker clade's core)
    hic: dict | None          # one row of the hicfuse table, or None (no Hi-C, or not a callable chromosome)
    meth: dict | None = None  # item -> {'label', 'score'} pooled over samples (methfuse.support), or None
    sat_bp: Callable | None = None   # (start, end) -> bp of centromeric-role satellite arrays inside
    fold: dict | None = None  # fold.call() of this chromosome + status (the genome's fold test) and agrees (its core
                              # holds the de novo peak), or None
    align: dict | None = None  # the same for the trans arm alignment (fold.align_scan), or None
    at: dict | None = None    # the AT-trough line on this chromosome: start, end, gc, status (the genome's FASTA test,
                              # atline.scan), validated (the genome's troughs agree with Hi-C's prior-free peaks) and
                              # agrees (they do here), or None
    age: dict | None = None   # ageline.compare: insertion ages of the signature's copies at the anchored Hi-C core vs at the
                              # sequence call (n_peak, n_call, med_peak, med_call, p one-sided: younger at the peak), or None


def _iv(d, s="start", e="end"):
    return (int(d[s]), int(d[e])) if d and int(d[s]) >= 0 else None


def hic_core(ev: Evidence):
    h = ev.hic
    return (int(h["anch_core_start"]), int(h["anch_core_end"])) if h and int(h["anch_core_start"]) >= 0 else None


def hic_peak(ev: Evidence):
    h = ev.hic
    return (int(h["anch_peak_start"]), int(h["anch_peak_end"])) if h and int(h["anch_peak_start"]) >= 0 else None


def _bins(ev: Evidence, n: int) -> int:
    return n * int(ev.hic["bin_size"])


def _vgap(ev: Evidence, a, b, capped: bool) -> int:
    """bp between intervals a and b that Hi-C can see: unmappable stretches of this chromosome removed (capped: at
    most HIC_FAR_BINS bins of them, hicfuse). Agreement is capped, opposition is not: across a blind stretch longer
    than the method's resolution Hi-C neither supports nor opposes a call."""
    h = ev.hic or {}
    return visible_gap(*a, *b, h.get("masked") or (), HIC_FAR_BINS * int(h.get("bin_size") or 0) if capped else None)


def _peak_agrees(ev: Evidence, iv) -> bool:
    """A significant anchored peak (Hi-C level not 'none') within HIC_TOL_BINS mappable bins of iv; a clean one when the
    gate's references were low calls (hicfuse.gate: there only clean agreements validate Hi-C), or included low calls and this
    chromosome's call is not confident (refs 'mixed')."""
    pk = hic_peak(ev)
    refs = ev.hic.get("gate_refs") if ev.hic else None
    need_clean = refs == "low" or (refs == "mixed" and ev.seq.get("confidence") not in ("high", "medium"))
    return bool(pk and iv and ev.hic.get("anch_level", "none") != "none"
                and (not need_clean or clean_peak(ev.hic.get("anch_flags")))
                and _vgap(ev, pk, iv, True) <= _bins(ev, HIC_TOL_BINS))


def in_end_zone(ev: Evidence, iv) -> bool:
    """The interval's midpoint lies in a terminal eighth of the chromosome (hicfuse.END_ZONE, Levan's telocentric bound)."""
    mid = (iv[0] + iv[1]) / 2
    return min(mid, ev.chrom_len - mid) < END_ZONE * ev.chrom_len


def _peak_far(ev: Evidence, iv) -> bool:
    """The peak points to another locus: farther from iv than the method's resolution (HIC_FAR_BINS mappable bins),
    and its core does not reach iv."""
    pk, core = hic_peak(ev), hic_core(ev)
    return bool(pk and iv and _vgap(ev, pk, iv, False) > _bins(ev, HIC_FAR_BINS)
                and (core is None or _vgap(ev, core, iv, False) > _bins(ev, HIC_TOL_BINS)))


def _te_agrees(ev: Evidence, iv) -> bool:
    """The satellite engine's rule (gap at most the wider interval's width): used for kept sequence calls, as that engine does."""
    return bool(ev.te and iv and agrees(iv[0], iv[1], dict(start=ev.te["start"], end=ev.te["end"])))


def _te_decides(ev: Evidence) -> bool:
    """A core that may act as a second line or oppose: a decisive marker's, and strong (signature.Signature.significant)."""
    return bool(ev.te and ev.te.get("decisive") and ev.te.get("strong", True))


def _te_near(ev: Evidence, iv, n: int) -> bool:
    """The TE core within n mappable Hi-C bins of iv (one fixed tolerance for both places a move compares)."""
    return bool(ev.te and iv and _vgap(ev, (ev.te["start"], ev.te["end"]), iv, False) <= _bins(ev, n))


DENOVO_UNCLEAN = ("not_converged", "cycle_unstable")     # with NOT_STRONG_FLAGS: a prior-free call that is not clean


def _clean_denovo(h) -> bool:
    """A prior-free Hi-C call at engine level high without a hub, rDNA, telomeric, not_converged or cycle_unstable
    flag."""
    fl = str(h.get("denovo_flags") or "")
    return bool(h.get("denovo_confidence") == "high" and clean_peak(fl)
                and not set(fl.split(",")) & set(DENOVO_UNCLEAN))


def _clean_flags(h) -> bool:
    """The prior-free call carries no hub, rDNA, telomeric, not_converged or cycle_unstable flag (any level)."""
    fl = str(h.get("denovo_flags") or "")
    return bool(clean_peak(fl) and not set(fl.split(",")) & set(DENOVO_UNCLEAN))


def _clean_medium(h) -> bool:
    """An anchored peak at Hi-C level medium without a hub, rDNA or telomeric flag."""
    return bool(h and h.get("anch_level") == "medium" and clean_peak(h.get("anch_flags")))


def _fold_far(ev: Evidence, f_core, iv) -> bool:
    """iv lies beyond the fold core by more than HIC_FAR_BINS mappable bins (opposition, uncapped as for peaks)."""
    return _vgap(ev, f_core, iv, False) > _bins(ev, HIC_FAR_BINS)


def _denovo_far(ev: Evidence, peak, core, iv) -> bool:
    """The de novo peak points to another locus (the rule of _peak_far, for the prior-free peak)."""
    return _vgap(ev, peak, iv, False) > _bins(ev, HIC_FAR_BINS) and _vgap(ev, core, iv, False) > _bins(ev, HIC_TOL_BINS)


def _meth(ev: Evidence, item: str) -> dict:
    return (ev.meth or {}).get(item) or {"label": "", "score": np.nan}


def genome_axis(fo: dict | None):
    """The genome's validated Rabl-axis line from the fold stage's result (pipeline.stage_fold): 'fold' if the arm
    fold's test passed; 'align' if the trans arm alignment's passed and the fold leans the same way (P <=
    AXIS_LEAN_P); else None. axis_line applies the same rule per chromosome; the summary reports this one."""
    if not fo:
        return None
    if fo.get("status") == "pass":
        return "fold"
    al = fo.get("align") or {}
    if al.get("status") == "pass" and float(fo.get("p", 1.0)) <= AXIS_LEAN_P:
        return "align"
    return None


def axis_line(ev: Evidence) -> tuple:
    """(the chromosome's entry of the genome's validated Rabl-axis line, its name): the arm fold if its test passed;
    else the trans arm alignment if its test passed and the fold's own test leans the same way (P <= AXIS_LEAN_P);
    else (None, None). Each test at GATE_P / 2 (fold.py). The alignment is read relative to the other chromosomes'
    prior-free peaks, so any arrangement symmetric about them looks Rabl-like (holocentric C. elegans, whose
    arm/center compartments are symmetric about each chromosome's middle: 3/5 agreements, P 0.007; its fold, read
    from each chromosome's own cis map, P 1); a Rabl axis shows in both readouts (yeast fold P 0.005, Plasmodium
    0.02, barley at half depth 0.013)."""
    if ev.fold and ev.fold.get("status") == "pass":
        return ev.fold, "fold"
    a = ev.align
    if a and a.get("status") == "pass" and a.get("fold_status") != "pass" and \
            float(a.get("fold_p", 1.0)) <= AXIS_LEAN_P:  # the genome's axis line: the fold where it passed (even on a
        return a, "align"                                # chromosome it does not call), else the corroborated alignment
    return None, None


def decide(ev: Evidence) -> dict:
    """One chromosome's integrated call (see module doc). With ev.meth None it is the provisional call. A long-unit
    tandem array that stays the call with no support from Hi-C, a retrotransposon core or methylation is reported as a
    feature, not a call (flag long_unit_unsupported): 24 such calls from the FASTA alone were 2 right, 16 wrong, 6 without
    truth. Satellite bp is no support (it says what else lies there), nor is an arm axis, which spans 5-10% of a
    chromosome. The chromosome is then decided as one without a sequence call, so Hi-C can still place one."""
    d = _decide(ev)
    if d["source"] == "sequence" and ev.seq["evidence"] == "long_unit_tandem" and \
            not {"hic", "te", "meth"} & set(d["support"]):
        d = _decide(replace(ev, seq=dict(ev.seq, start=-1, end=-1, evidence="none", confidence="none", alt_start=-1,
                                         alt_end=-1)))
        d["flags"] = d["flags"] + ["long_unit_unsupported"]
    return d


def _decide(ev: Evidence) -> dict:
    seq_iv = _iv(ev.seq)
    sconf = ev.seq["confidence"] if seq_iv else "none"
    alt_iv = _iv(ev.seq, "alt_start", "alt_end") if seq_iv and "alt_start" in ev.seq else None
    h = ev.hic
    valid = bool(h and h.get("gate_pass"))
    ax, axis = axis_line(ev) if h and not valid else (None, None)   # the validated Rabl-axis line (fold, alignment)
    weak = bool(h and not valid and not ax and h.get("gate_support"))   # support tier: Hi-C may only support
    core, peak = hic_core(ev), hic_peak(ev)
    strong = bool(valid and h.get("anch_strong") and core)
    agree = bool(valid and seq_iv and _peak_agrees(ev, seq_iv))
    flags, relation, source = [], "sequence_only", "sequence"
    # Hi-C picks the competing cluster, only if mappable bins separate it from the other (it can tell the two apart), and with a
    # strong peak. On the support tier too (0.13): there Hi-C may support what its peak agrees with, and the competing cluster is
    # the sequence line's own second candidate (moss without tes/: 5 calls, all onto CENH3); the call then rises as support does.
    # Where Hi-C read on the LTR cores (te_pair) or its cis signal acts on this chromosome, those tested readings decide instead
    alt_peak = strong or bool(weak and h.get("anch_strong") and core and not h.get("te_pair") and not h.get("cis"))
    if alt_peak and alt_iv and not agree and _peak_agrees(ev, alt_iv) and peak and \
            _vgap(ev, peak, seq_iv, False) > _bins(ev, HIC_TOL_BINS):
        seq_iv, agree = alt_iv, strong
        flags.append("hic_chose_alt")
    far = bool(seq_iv and _peak_far(ev, seq_iv))
    # where no mappable bin separates the call from its competing cluster, Hi-C cannot tell the two apart (the
    # principle of hic_chose_alt, which needs mappable bins between them to choose): an agreement with the call is
    # no evidence for it over the competitor, so Hi-C neither raises nor supports the call (mouse: the minor- and
    # major-satellite arrays at each chromosome's start, the peak at the inner edge of the blind block)
    amb = bool(h and alt_iv and seq_iv and "hic_chose_alt" not in flags and _vgap(ev, seq_iv, alt_iv, False) == 0)
    fo = ax
    f_core = (fo["core_start"], fo["core_end"]) if fo else None
    # the alignment (newer, conditional on the prior-free positions) pairs and opposes only through a clean prior-free
    # peak (no hub, rDNA, telomeric, not_converged, cycle_unstable flag); the fold keeps 0.10's rule (barley: its 6
    # moves, one at an rDNA-flagged peak, were right)
    dn_ok = axis != "align" or _clean_flags(h)
    dn_peak = _iv(h, "denovo_peak_start", "denovo_peak_end") if fo and dn_ok else None
    dn_core = _iv(h, "denovo_core_start", "denovo_core_end") if fo and dn_ok else None
    pair = bool(fo and fo.get("agrees") and dn_peak and dn_core)
    untested = bool(h and h.get("gate_genome_untested") and not valid and _clean_denovo(h))
    # Hi-C alone: also where a passed axis line does not pair on this chromosome (0.10: no call there)
    dn_only = _iv(h, "denovo_core_start", "denovo_core_end") if untested and not pair and seq_iv is None else None
    m_seq, m_hic = _meth(ev, "sequence_call"), _meth(ev, "hic_peak")
    at_ok = bool(ev.at and (ev.at.get("status") == "pass" or ev.at.get("validated")))   # the genome's troughs are tested
    second = []
    if strong:
        if _te_decides(ev) and _te_near(ev, peak, HIC_FAR_BINS) and not (seq_iv and _te_near(ev, seq_iv, HIC_FAR_BINS)):
            second.append("te")
        if m_hic["label"] == "centromere_like" and (seq_iv is None or m_seq["label"] == "dip"):
            second.append("meth")
        if seq_iv and ev.age and float(ev.age.get("p", np.nan)) <= AGE_ALPHA and \
                float(ev.age.get("premise", np.nan)) <= AGE_ALPHA:     # younger centromeric copies at the peak than at the call,
            second.append("age")                                         # in a genome where young copies mark the cores (#31)
    if pair and (seq_iv is None or (_fold_far(ev, f_core, seq_iv) and _denovo_far(ev, dn_peak, dn_core, seq_iv))):
        call, source = dn_core, "hic"                              # two validated 3D lines at one locus, the
        relation = "functional_only" if seq_iv is None else "offset"   # sequence call at neither
        flags.append(f"{axis}_pair")
        if sconf == "high":
            flags.append(f"moved_with_{axis}")
    elif seq_iv is None and strong:
        call, relation, source = core, "functional_only", "hic"
    elif seq_iv is None and valid and core and _clean_medium(h):
        call, relation, source = core, "functional_only", "hic"
        flags.append("hic_medium_peak")
    elif seq_iv is None and dn_only:
        call, relation, source = dn_only, "functional_only", "hic"
        flags.append("hic_denovo_only")
    elif seq_iv is None:
        call, relation = None, "none"
    elif agree and amb:
        call = seq_iv
        flags.append("hic_alt_ambiguous")
    elif agree:
        call, relation = seq_iv, "concordant"
    elif weak and _peak_agrees(ev, seq_iv) and h.get("anch_level") in SUPPORT_LEVELS:
        call = seq_iv
        if amb:
            flags.append("hic_alt_ambiguous")
        else:
            relation = "concordant"
            flags.append("hic_support")
    elif not valid and not ax and h and h.get("end_pass") and h.get("end_higher") and not amb:
        call, relation = seq_iv, "concordant"          # Hi-C tells the call's end from the other, in a genome where it does so
        flags.append("hic_end")
    elif strong and far and (sconf != "high" or second):
        call, relation, source = core, "offset", "hic"
        if sconf == "high":
            flags.append("moved_with_" + "+".join(second))
    elif strong and far:
        call, relation = seq_iv, "contested"
        flags.append("hic_offset")
    else:
        call = seq_iv
        if valid and peak:
            flags.append("hic_near" if not far else "hic_disagrees_weak")
    pair = (h or {}).get("te_pair")      # Hi-C's peak at this chromosome's LTR core, the two lines tested genome-wide
    # both readings run where the gate fails on most chromosomes; a chromosome whose own gate passed keeps what validated Hi-C said
    if pair and not valid and source == "sequence" and call is not None and relation != "contested":
        if gap(*call, *pair["core"]) == 0:                                           # both at the call: Hi-C supports it
            if not amb:                                                              # no mappable bin between it and its competitor
                relation = "concordant"
                flags.append("te_hic_pair")
            elif "hic_alt_ambiguous" not in flags:
                flags.append("hic_alt_ambiguous")
        elif sconf == "low" or pair["status"] == "pass":     # two lines at one locus against one line, unsure or on a minority
            call, relation, source = pair["core"], "offset", "hic"
            flags.append("te_hic_pair")
            if sconf == "high":
                second.append("te")
                flags.append("moved_with_te")
    elif pair and call is None:
        call, relation, source = pair["core"], "functional_only", "hic"
        flags.append("te_hic_pair")
    cis = (h or {}).get("cis")           # the cis signal's extreme, its sign and level tested against the LTR cores or trans contacts
    # a line read outside the terminal eighths (cis_trans_line) has no say on a call inside one: it cannot see there
    blind = bool(cis and cis.get("interior") and call is not None and in_end_zone(ev, call))
    if cis and not valid and not blind and "te_hic_pair" not in flags and source == "sequence" and call is not None and \
            relation != "contested":
        if _vgap(ev, cis["peak"], call, True) <= _bins(ev, HIC_TOL_BINS):              # at the call: an independent line agrees
            if not amb:
                relation = "concordant"
                flags.append("cis_agrees")
            elif "hic_alt_ambiguous" not in flags:
                flags.append("hic_alt_ambiguous")
        elif _vgap(ev, cis["peak"], call, False) <= _bins(ev, HIC_FAR_BINS):           # within Hi-C's resolution: neither way
            flags.append("cis_near")
        elif sconf == "low":                 # elsewhere, the call unsure: the better candidate, still unsure
            call, relation, source = cis["peak"], "offset", "hic"
            flags.append("cis_only")
        else:
            flags.append("cis_elsewhere")
            if cis["status"] == "pass":      # a line at the LTR cores on most chromosomes lies elsewhere: one line against one
                relation = "contested"
    elif cis and call is None:
        call, relation, source = cis["peak"], "functional_only", "hic"
        flags.append("cis_only")
    elif blind and not valid and source == "sequence" and _vgap(ev, cis["peak"], call, False) > _bins(ev, HIC_FAR_BINS):
        flags.append("cis_blind")            # said in the caveat: the extreme lies elsewhere, but the line cannot judge an end
    # two lines at a chromosome nothing else places (0.13): the trough and a clean prior-free Hi-C peak that placed no call of its
    # own agree here, in a genome where the two agree beyond chance (atline.validate: support); C. neoformans in anaphase: 3
    at_two = bool(ev.at and ev.at.get("support") and ev.at.get("agrees") and h and _clean_flags(h)
                  and h.get("denovo_confidence", "none") != "none")
    if call is None and (at_ok or at_two):           # nothing else places this chromosome
        call, relation, source = (ev.at["start"], ev.at["end"]), "sequence_only", "composition"
        flags.append("at_only" if at_ok else "at_denovo_pair")
    if call is None:
        return dict(call=None, relation="none", confidence="none", source="none", support=[], against=[],
                    flags=flags, second=second, meth_call={"label": "", "score": np.nan})
    if seq_iv is None and source == "hic" and ev.at:
        trough = (ev.at["start"], ev.at["end"])
        at_near = _vgap(ev, trough, call, True) <= _bins(ev, HIC_TOL_BINS)
    else:
        trough, at_near = None, False
    at_pair = bool(at_near and ev.at.get("agrees") and (ev.at.get("validated") or ev.at.get("support")))
    if at_pair:                  # two independent lines on a chromosome without a sequence call: the trough is the finer
        call = trough                                                # (a window of 1-500 kb, a Hi-C bin of 10-250 kb)
        flags.append("at_pair")
    elif trough and not at_near and ev.at.get("validated"):
        flags.append("at_elsewhere")      # said in the caveat; a Hi-C call stands against a trough of the same standing
    if untested and not fo and source == "sequence":
        dpk, dco = _iv(h, "denovo_peak_start", "denovo_peak_end"), _iv(h, "denovo_core_start", "denovo_core_end")
        if dpk and dco and _denovo_far(ev, dpk, dco, call):
            flags.append("hic_denovo_elsewhere")
    f_opp = d_opp = False
    if fo and source == "sequence":
        f_opp, d_opp = _fold_far(ev, f_core, call), bool(dn_peak and _denovo_far(ev, dn_peak, dn_core, call))
        if f_opp or d_opp:
            relation = "contested"
            flags += [x for x, y in ((f"{axis}_offset", f_opp), ("hic_denovo_offset", d_opp)) if y]
    elif fo and "hic_denovo_only" in flags and _fold_far(ev, f_core, call):
        f_opp = True                                       # the axis line passed genome-wide but points elsewhere
        flags.append(f"{axis}_offset")                     # here: said, not hidden (the call is low already)
    at_core = "hic_denovo_only" in flags or (axis and f"{axis}_pair" in flags)   # the call sits at the prior-free
    m_call = _meth(ev, "provisional") if at_core else m_hic if source == "hic" else m_seq   # core: its own item
    support = ["at"] if at_pair or {"at_only", "at_denovo_pair"} & set(flags) else []
    if ev.sat_bp and ev.sat_bp(*call) >= SAT_SUPPORT_BP:
        support.append("satellite")
    if (source == "sequence" and _te_agrees(ev, call)) or (source == "hic" and "te" in second) or "te_hic_pair" in flags:
        support.append("te")
    if source == "hic" and "age" in second:
        support.append("age")
    if valid and "hic_alt_ambiguous" not in flags and (_peak_agrees(ev, call) or
                                                       (source == "hic" and h.get("anch_level") != "none")):
        support.append("hic")
    elif {"hic_support", "hic_denovo_only", "te_hic_pair", "hic_end", "cis_agrees", "cis_only", "at_denovo_pair"} & set(flags):
        support.append("hic")
    if axis and f"{axis}_pair" in flags:
        support += ["hic", axis]
    elif fo and source == "sequence" and gap(*call, *f_core) == 0:
        support.append(axis)
    if m_call["label"] == "centromere_like":
        support.append("meth")
    against = []
    if valid and strong and _peak_far(ev, call):
        against.append("hic")
    if _te_decides(ev) and not (_te_agrees(ev, call) if source == "sequence" else "te" in second or "te_hic_pair" in flags):
        against.append("te")
    if "at_elsewhere" in flags:
        against.append("at")
    if "cis_elsewhere" in flags:
        against.append("hic")
    against += [x for x, y in ((axis, f_opp), ("hic", d_opp)) if y]
    if source == "sequence":
        if relation == "contested":
            conf = level_shift(sconf, -max(1, int(f_opp) + int(d_opp)) + int("meth" in support))
        else:
            raise_by = int("hic" in support) + int("meth" in support)
            conf = level_shift(sconf, min(MAX_RAISE, raise_by))
        if relation == "sequence_only" and "meth" in support:
            relation = "concordant"
    elif source == "composition":      # a trough alone is low; tested against Hi-C's peaks genome-wide, medium
        conf = "medium" if ev.at.get("validated") else "low"
    elif "cis_only" in flags:
        conf = "low"
    elif "hic_medium_peak" in flags or "hic_denovo_only" in flags:
        conf = "medium" if at_pair else "low"
    else:
        conf = "high" if second or at_pair else "medium"
    return dict(call=call, relation=relation, confidence=conf, source=source, support=support, against=against,
                flags=flags, second=second, meth_call=m_call)


def row(ev: Evidence, d: dict) -> dict:
    """The output row of one chromosome (coordinates 1-based inclusive like the satellite engine's summary; '.' when absent)."""
    na = "."
    h, call, seq_iv = ev.hic or {}, d["call"], _iv(ev.seq)
    one = lambda x: int(x) + 1 if x is not None and int(x) >= 0 else na
    end = lambda x: int(x) if x is not None and int(x) >= 0 else na
    m = {k: _meth(ev, k) for k in ("sequence_call", "hic_peak")}
    m["call"] = d["meth_call"]
    off = na
    alt_iv = _iv(ev.seq, "alt_start", "alt_end") if seq_iv and "alt_start" in ev.seq else None
    if call and seq_iv and call != seq_iv and not (alt_iv and call == alt_iv):
        off = gap(*call, *seq_iv)
    over_iv = _iv(ev.seq, "overruled_start", "overruled_end") if seq_iv and "overruled_start" in ev.seq else None
    return dict(
        chrom=ev.chrom, chrom_len=ev.chrom_len, start=one(call[0]) if call else na, end=end(call[1]) if call else na,
        length=call[1] - call[0] if call else na, call_source=d["source"], relation=d["relation"],
        confidence=d["confidence"], support=",".join(d["support"]) or na, against=",".join(d["against"]) or na,
        flags=",".join(d["flags"]) or na, offset_bp=off,
        seq_start=one(seq_iv[0]) if seq_iv else na, seq_end=end(seq_iv[1]) if seq_iv else na,
        seq_evidence=ev.seq["evidence"], seq_confidence=ev.seq["confidence"], seq_flags=ev.seq["flags"],
        seq_alt_start=one(alt_iv[0]) if alt_iv else na, seq_alt_end=end(alt_iv[1]) if alt_iv else na,
        seq_overruled_start=one(over_iv[0]) if over_iv else na, seq_overruled_end=end(over_iv[1]) if over_iv else na,
        seq_rejected_family=ev.seq.get("rejected_family", na), seq_rejected_start=one(ev.seq.get("rejected_start")),
        seq_rejected_end=end(ev.seq.get("rejected_end")), seq_spacing=end(ev.seq.get("spacing")),
        main_family=ev.seq["main_family"], te_clade=ev.te["clade"] if ev.te else na,
        te_seed=(ev.te.get("seed") or na) if ev.te else na,
        te_start=one(ev.te["start"]) if ev.te else na, te_end=end(ev.te["end"]) if ev.te else na,
        te_share=round(ev.te["share"], 3) if ev.te and ev.te.get("share") is not None else na,
        te_chrom_elements=ev.te.get("n", na) if ev.te else na,
        te_core_elements=end(ev.seq.get("te_core_elements")), te_call_elements=end(ev.seq.get("te_call_elements")),
        at_gate=("pass_hic" if ev.at.get("validated") else ev.at["status"]) if ev.at else na,
        at_start=one(ev.at["start"]) if ev.at else na, at_end=end(ev.at["end"]) if ev.at else na,
        at_gc=round(ev.at["gc"], 3) if ev.at else na,
        hic_gate=gate_status(h) if h else na,
        hic_gate_k=h.get("gate_k", na) if h else na, hic_gate_n=h.get("gate_n", na) if h else na,
        hic_gate_p=h.get("gate_p", na) if h else na,
        hic_gate_refs=h.get("gate_refs", na) if h else na,
        hic_end_k=h.get("end_k", na) if h else na, hic_end_n=h.get("end_n", na) if h else na,
        hic_anch_peak_start=one(h.get("anch_peak_start", -1)) if h else na,
        hic_anch_peak_end=end(h.get("anch_peak_end", -1)) if h else na,
        hic_anch_core_start=one(h.get("anch_core_start", -1)) if h else na,
        hic_anch_core_end=end(h.get("anch_core_end", -1)) if h else na,
        hic_anch_level=h.get("anch_level", na) if h else na, hic_anch_flags=h.get("anch_flags", "") or na,
        hic_anch_agrees_seq=h.get("anch_agrees_seq", na) if h else na,
        hic_denovo_peak_start=one(h.get("denovo_peak_start", -1)) if h else na,
        hic_denovo_peak_end=end(h.get("denovo_peak_end", -1)) if h else na,
        hic_denovo_confidence=h.get("denovo_confidence", na) if h else na,
        hic_denovo_flags=h.get("denovo_flags", "") or na,
        fold_gate=ev.fold["status"] if ev.fold else na, fold_point=one(ev.fold["point"]) if ev.fold else na,
        fold_core_start=one(ev.fold["core_start"]) if ev.fold else na,
        fold_core_end=end(ev.fold["core_end"]) if ev.fold else na,
        align_gate=ev.align["status"] if ev.align else na, align_point=one(ev.align["point"]) if ev.align else na,
        align_core_start=one(ev.align["core_start"]) if ev.align else na,
        align_core_end=end(ev.align["core_end"]) if ev.align else na,
        cis_peak_start=one(h["cis"]["peak"][0]) if h.get("cis") else na,
        cis_peak_end=end(h["cis"]["peak"][1]) if h.get("cis") else na,
        cis_test=h["cis"].get("tested_by", "ltr") if h.get("cis") else na,
        age_peak=_age(ev, "med_peak", "n_peak"), age_call=_age(ev, "med_call", "n_call"),
        age_p=float(f"{ev.age['p']:.3g}") if ev.age and np.isfinite(ev.age.get("p", np.nan)) else na,
        age_premise_p=float(f"{ev.age['premise']:.3g}") if ev.age and np.isfinite(ev.age.get("premise", np.nan)) else na,
        meth_seq=m["sequence_call"]["label"] or na, meth_seq_score=m["sequence_call"]["score"],
        meth_hic=m["hic_peak"]["label"] or na, meth_hic_score=m["hic_peak"]["score"],
        meth_call=m["call"]["label"] or na, meth_call_score=m["call"]["score"])


def _age(ev: Evidence, med: str, n: str):
    """'median years (n copies)' of the signature's dated copies at one site (ageline.compare), or '.'."""
    a = ev.age or {}
    if not a.get(n):
        return "."
    return f"{a[med]:.0f} ({a[n]})"


def gate_status(h: dict) -> str:
    """pass / support (the support tier only: hicfuse.gate) / fail / untested (fewer than GATE_MIN_REF reference
    chromosomes: no majority test was possible)."""
    if h.get("gate_pass"):
        return "pass"
    if h.get("gate_support"):
        return "support"
    n = h.get("gate_n")
    return "untested" if n is None or (isinstance(n, (int, float, np.integer)) and int(n) < GATE_MIN_REF) else "fail"


def table(evs: list, decisions: list) -> pd.DataFrame:
    return pd.DataFrame([row(e, d) for e, d in zip(evs, decisions)], columns=FUSED_COLS)


def sequence_conflict(t: pd.DataFrame) -> dict:
    """Genome-level diagnostic. Hi-C's prior-free calls are confident (Hi-C high/medium) on most chromosomes but agree
    with the sequence calls on fewer than half, and the gate fails on most: the two lines disagree genome-wide, so at
    least one is wrong everywhere, and two lines cannot tell which. A third line decides: the arm fold (fold.py, from
    the same Hi-C; Rabl nuclei only), tes/ (CRM cores) or methylation; without one the sequence calls are kept. Seen
    both ways: barley (satellite line wrong, prior-free Hi-C right; the fold sides with Hi-C) and S. orangutan /
    giraffe (the reverse; no fold)."""
    s = t[t.seq_start.ge(0) & t.callable]
    n = len(s)
    conf = int(s.denovo_confidence.isin(["high", "medium"]).sum())
    agree = int(s.denovo_agrees_seq.sum())
    gate = int(s.gate_pass.sum())
    flag = bool(n >= GATE_MIN_REF and 2 * conf > n and 2 * agree < n and 2 * gate < n)
    return dict(conflict=flag, n=n, denovo_confident=conf, denovo_on_seq=agree, gate_pass=gate)
