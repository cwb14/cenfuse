"""Integration rules (docs/design.md, 'Integration').

The integrated call is the sequence call unless a validated functional line moves it:
- Hi-C may act on chromosome x only when its gate passes (hicfuse.gate), and only through a strong anchored peak
  (all three Hi-C engine tests, no hub/rDNA/telomeric flag). A peak supports a call within HIC_TOL_BINS, only when it is a
  call at all (Hi-C level not 'none').
- A strong peak that points to another locus (more than HIC_FAR_BINS from the sequence call, the method's resolution,
  and with a core that does not reach the call) replaces a low or medium sequence call, and places a call on a
  chromosome without one. Closer than that, Hi-C neither supports nor opposes the call (flag hic_near).
- It replaces a high sequence call only with a second line at the Hi-C peak and not at the sequence call: a decisive
  TE core within HIC_FAR_BINS of the peak and beyond it from the call, or methylation `centromere_like` at the
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
  calls were right on 14 of 63 T2T mammal chromosomes; with 10-kb bins on 16/16 yeast and 14/14 Plasmodium ones).
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
  own agrees with it (at_denovo_pair). A trough never moves or opposes a confident satellite or LTR call; a passing genome's
  trough at the competing cluster of a low satellite call (its deepest window there, not at the call) chooses that cluster, still
  low (0.15, at_chose_alt; mouse: the minor satellite, not the pericentromeric major satellite).

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

from dataclasses import dataclass, field, replace
from typing import Callable

import numpy as np
import pandas as pd

from .hicfuse import END_ZONE, SUPPORT_LEVELS, clean_peak, visible_gap
from .ageline import AGE_ALPHA
from .params import AXIS_LEAN_P, GATE_MIN_REF, HIC_FAR_BINS, HIC_TOL_BINS, MAX_RAISE, SAT_SUPPORT_BP
from .sat.te import agrees
from .util import gap, level_shift

FUSED_COLS = ["chrom", "chrom_len", "start", "end", "length", "call_source", "relation", "confidence", "support",
              "against", "flags", "rule", "offset_bp", "seq_start", "seq_end", "seq_evidence", "seq_confidence", "seq_flags",
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
    """A prior-free Hi-C call at engine level high without a hub, rDNA, telomeric, not_converged or cycle_unstable flag. A call
    held at medium only because its null test cannot be run has passed every test it could take, and counts (null_untestable,
    hic.confidence.level; 0.15, C. tropicalis: too few partner chromosomes with bins to resample; the same flag marks the rarer
    degenerate case of a non-finite observed ratio, not separated here: future.md #69)."""
    fl = str(h.get("denovo_flags") or "")
    level = h.get("denovo_confidence")
    tested = level == "high" or (level == "medium" and "null_untestable" in fl.split(","))
    return bool(tested and clean_peak(fl) and not set(fl.split(",")) & set(DENOVO_UNCLEAN))


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


@dataclass
class _Ctx:
    """What the rules read about one chromosome, worked out once from its Evidence (_context)."""
    ev: Evidence
    seq_iv: tuple | None       # the sequence call: its competing cluster when Hi-C chose that
    sconf: str                 # the sequence call's confidence ('none' without one)
    alt_iv: tuple | None       # the competing cluster
    h: dict | None             # this chromosome's Hi-C row
    valid: bool                # Hi-C's test passed here
    ax: dict | None            # the validated Rabl-axis line's entry (fold or alignment), where the test did not pass
    axis: str | None           # its name
    weak: bool                 # the support tier: Hi-C may only support
    core: tuple | None         # the anchored Hi-C core and peak
    peak: tuple | None
    strong: bool               # a strong anchored peak at a validated chromosome
    agree: bool                # it agrees with the sequence call
    chose: bool                # Hi-C chose the competing cluster
    far: bool                  # the strong peak points elsewhere
    amb: bool                  # no mappable bin between the call and its competing cluster
    f_core: tuple | None       # the axis line's core
    dn_peak: tuple | None      # the prior-free peak and core the axis line reads
    dn_core: tuple | None
    pair: bool                 # the axis line and the prior-free peak at one locus
    untested: bool             # Hi-C alone: nothing can test it in this genome, and its prior-free call is clean
    dn_only: tuple | None      # the prior-free core where Hi-C alone places a call
    m_seq: dict                # methylation at the sequence call and at Hi-C's peak
    m_hic: dict
    at_ok: bool                # the genome's AT troughs are tested


@dataclass
class _Call:
    """The call as the rules leave it: place, relation, where it came from, the rule that set it, flags, second lines."""
    call: tuple | None = None
    relation: str = "sequence_only"
    source: str = "sequence"
    rule: str = "sequence"
    flags: list = field(default_factory=list)
    second: list = field(default_factory=list)
    at_pair: bool = False      # the AT trough refined a Hi-C call
    f_opp: bool = False        # the axis line, the prior-free peak, oppose the kept call
    d_opp: bool = False


def _context(ev: Evidence) -> tuple:
    """(_Ctx, the _Call it starts from): every quantity the rules share, and Hi-C's choice of the competing cluster."""
    seq_iv = _iv(ev.seq)
    sconf = ev.seq["confidence"] if seq_iv else "none"
    alt_iv = _iv(ev.seq, "alt_start", "alt_end") if seq_iv and "alt_start" in ev.seq else None
    h = ev.hic
    valid = bool(h and h.get("gate_pass"))
    ax, axis = axis_line(ev) if h and not valid else (None, None)
    weak = bool(h and not valid and not ax and h.get("gate_support"))
    core, peak = hic_core(ev), hic_peak(ev)
    strong = bool(valid and h.get("anch_strong") and core)
    agree = bool(valid and seq_iv and _peak_agrees(ev, seq_iv))
    st = _Call()
    # Hi-C picks the competing cluster, only if mappable bins separate it from the other (it can tell the two apart), and with a
    # strong peak. On the support tier too (0.13): there Hi-C may support what its peak agrees with, and the competing cluster is
    # the sequence line's own second candidate (moss without tes/: 5 calls, all onto CENH3); the call then rises as support does.
    # Where Hi-C read on the LTR cores (te_pair) or its cis signal acts on this chromosome, those tested readings decide instead
    alt_peak = strong or bool(weak and h.get("anch_strong") and core and not h.get("te_pair") and not h.get("cis"))
    if alt_peak and alt_iv and not agree and _peak_agrees(ev, alt_iv) and peak and \
            _vgap(ev, peak, seq_iv, False) > _bins(ev, HIC_TOL_BINS):
        seq_iv, agree = alt_iv, strong
        st.flags.append("hic_chose_alt")
    chose = "hic_chose_alt" in st.flags
    far = bool(seq_iv and _peak_far(ev, seq_iv))
    # where no mappable bin separates the call from its competing cluster, Hi-C cannot tell the two apart (the
    # principle of hic_chose_alt, which needs mappable bins between them to choose): an agreement with the call is
    # no evidence for it over the competitor, so Hi-C neither raises nor supports the call (mouse: the minor- and
    # major-satellite arrays at each chromosome's start, the peak at the inner edge of the blind block)
    amb = bool(h and alt_iv and seq_iv and not chose and _vgap(ev, seq_iv, alt_iv, False) == 0)
    f_core = (ax["core_start"], ax["core_end"]) if ax else None
    # the alignment (newer, conditional on the prior-free positions) pairs and opposes only through a clean prior-free
    # peak (no hub, rDNA, telomeric, not_converged, cycle_unstable flag); the fold keeps 0.10's rule (barley: its 6
    # moves, one at an rDNA-flagged peak, were right)
    dn_ok = axis != "align" or _clean_flags(h)
    dn_peak = _iv(h, "denovo_peak_start", "denovo_peak_end") if ax and dn_ok else None
    dn_core = _iv(h, "denovo_core_start", "denovo_core_end") if ax and dn_ok else None
    pair = bool(ax and ax.get("agrees") and dn_peak and dn_core)
    untested = bool(h and h.get("gate_genome_untested") and not valid and _clean_denovo(h))
    # Hi-C alone: also where a passed axis line does not pair on this chromosome (0.10: no call there)
    dn_only = _iv(h, "denovo_core_start", "denovo_core_end") if untested and not pair and seq_iv is None else None
    m_seq, m_hic = _meth(ev, "sequence_call"), _meth(ev, "hic_peak")
    at_ok = bool(ev.at and (ev.at.get("status") == "pass" or ev.at.get("validated")))
    if strong:                   # a second line at the strong peak and not at the call: what lets Hi-C move a high call
        if _te_decides(ev) and _te_near(ev, peak, HIC_FAR_BINS) and not (seq_iv and _te_near(ev, seq_iv, HIC_FAR_BINS)):
            st.second.append("te")
        if m_hic["label"] == "centromere_like" and (seq_iv is None or m_seq["label"] == "dip"):
            st.second.append("meth")
        if seq_iv and ev.age and float(ev.age.get("p", np.nan)) <= AGE_ALPHA and \
                float(ev.age.get("premise", np.nan)) <= AGE_ALPHA:     # younger centromeric copies at the peak than at the call,
            st.second.append("age")                                      # in a genome where young copies mark the cores (#31)
    c = _Ctx(ev, seq_iv, sconf, alt_iv, h, valid, ax, axis, weak, core, peak, strong, agree, chose, far, amb, f_core, dn_peak,
             dn_core, pair, untested, dn_only, m_seq, m_hic, at_ok)
    return c, st


# ---- the cascade: the first rule that applies sets the call (each returns whether it applied)
def _axis_pair(c: _Ctx, s: _Call) -> bool:
    """Two validated 3D lines at one locus (the axis line's core holds the prior-free peak), the sequence call at neither."""
    if not (c.pair and (c.seq_iv is None or (_fold_far(c.ev, c.f_core, c.seq_iv) and
                                             _denovo_far(c.ev, c.dn_peak, c.dn_core, c.seq_iv)))):
        return False
    s.call, s.source = c.dn_core, "hic"
    s.relation = "functional_only" if c.seq_iv is None else "offset"
    s.flags.append(f"{c.axis}_pair")
    s.rule = f"{c.axis}_pair"
    if c.sconf == "high":
        s.flags.append(f"moved_with_{c.axis}")
    return True


def _hic_strong_peak(c: _Ctx, s: _Call) -> bool:
    """No sequence call: a strong validated peak places one at its core."""
    if not (c.seq_iv is None and c.strong):
        return False
    s.call, s.relation, s.source, s.rule = c.core, "functional_only", "hic", "hic_strong_peak"
    return True


def _hic_medium_peak(c: _Ctx, s: _Call) -> bool:
    """No sequence call, the test passed: a clean medium peak places a low call."""
    if not (c.seq_iv is None and c.valid and c.core and _clean_medium(c.h)):
        return False
    s.call, s.relation, s.source, s.rule = c.core, "functional_only", "hic", "hic_medium_peak"
    s.flags.append("hic_medium_peak")
    return True


def _hic_denovo_only(c: _Ctx, s: _Call) -> bool:
    """Hi-C alone, in a genome where nothing can test it: its clean prior-free core."""
    if not (c.seq_iv is None and c.dn_only):
        return False
    s.call, s.relation, s.source, s.rule = c.dn_only, "functional_only", "hic", "hic_denovo_only"
    s.flags.append("hic_denovo_only")
    return True


def _no_call(c: _Ctx, s: _Call) -> bool:
    if c.seq_iv is not None:
        return False
    s.call, s.relation, s.rule = None, "none", "none"
    return True


def _hic_agrees_ambiguous(c: _Ctx, s: _Call) -> bool:
    """Validated Hi-C agrees, but no mappable bin separates the call from its competing cluster: no support."""
    if not (c.agree and c.amb):
        return False
    s.call, s.rule = c.seq_iv, "hic_chose_alt" if c.chose else "sequence"
    s.flags.append("hic_alt_ambiguous")
    return True


def _hic_agrees(c: _Ctx, s: _Call) -> bool:
    if not c.agree:
        return False
    s.call, s.relation, s.rule = c.seq_iv, "concordant", "hic_chose_alt" if c.chose else "hic_agrees"
    return True


def _hic_support(c: _Ctx, s: _Call) -> bool:
    """The support tier: the call's own significant peak agrees."""
    if not (c.weak and _peak_agrees(c.ev, c.seq_iv) and c.h.get("anch_level") in SUPPORT_LEVELS):
        return False
    s.call, s.rule = c.seq_iv, "hic_chose_alt" if c.chose else "sequence"
    if c.amb:
        s.flags.append("hic_alt_ambiguous")
    else:
        s.relation = "concordant"
        s.flags.append("hic_support")
        s.rule = "hic_chose_alt" if c.chose else "hic_support"
    return True


def _hic_end(c: _Ctx, s: _Call) -> bool:
    """Hi-C tells the call's end from the other, in a genome where it does so on most chromosomes."""
    if not (not c.valid and not c.ax and c.h and c.h.get("end_pass") and c.h.get("end_higher") and not c.amb):
        return False
    s.call, s.relation, s.rule = c.seq_iv, "concordant", "hic_end"
    s.flags.append("hic_end")
    return True


def _hic_moved(c: _Ctx, s: _Call) -> bool:
    """A strong validated peak elsewhere moves a low or medium call, and a high one with a second line there."""
    if not (c.strong and c.far and (c.sconf != "high" or s.second)):
        return False
    s.call, s.relation, s.source, s.rule = c.core, "offset", "hic", "hic_moved"
    if c.sconf == "high":
        s.flags.append("moved_with_" + "+".join(s.second))
    return True


def _hic_offset(c: _Ctx, s: _Call) -> bool:
    """A strong validated peak elsewhere, no second line: a high call is kept and contested."""
    if not (c.strong and c.far):
        return False
    s.call, s.relation, s.rule = c.seq_iv, "contested", "hic_offset"
    s.flags.append("hic_offset")
    return True


def _kept(c: _Ctx, s: _Call) -> bool:
    """The sequence call as it is (a validated peak beside it or too weak to oppose is said in the flags)."""
    s.call, s.rule = c.seq_iv, "hic_chose_alt" if c.chose else "sequence"
    if c.valid and c.peak:
        s.flags.append("hic_near" if not c.far else "hic_disagrees_weak")
    return True


CASCADE = (_axis_pair, _hic_strong_peak, _hic_medium_peak, _hic_denovo_only, _no_call, _hic_agrees_ambiguous, _hic_agrees,
           _hic_support, _hic_end, _hic_moved, _hic_offset, _kept)


# ---- the phases: readings that act where the gate fails, in this order, each on the call the earlier ones left
def _te_hic_pair(c: _Ctx, s: _Call) -> None:
    """Hi-C's peak at this chromosome's LTR core, the two lines tested genome-wide (hicfuse.anchored_on_cores). Both readings run
    where the gate fails on most chromosomes; a chromosome whose own gate passed keeps what validated Hi-C said."""
    pair = (c.h or {}).get("te_pair")
    if pair and not c.valid and s.source == "sequence" and s.call is not None and s.relation != "contested":
        if gap(*s.call, *pair["core"]) == 0:                                      # both at the call: Hi-C supports it
            if not c.amb:                                                         # no mappable bin between it and its competitor
                s.relation, s.rule = "concordant", "te_hic_pair"
                s.flags.append("te_hic_pair")
            elif "hic_alt_ambiguous" not in s.flags:
                s.flags.append("hic_alt_ambiguous")
        elif c.sconf == "low" or pair["status"] == "pass":  # two lines at one locus against one line, unsure or on a minority
            s.call, s.relation, s.source, s.rule = pair["core"], "offset", "hic", "te_hic_pair"
            s.flags.append("te_hic_pair")
            if c.sconf == "high":
                s.second.append("te")
                s.flags.append("moved_with_te")
    elif pair and s.call is None:
        s.call, s.relation, s.source, s.rule = pair["core"], "functional_only", "hic", "te_hic_pair"
        s.flags.append("te_hic_pair")


def _cis(c: _Ctx, s: _Call) -> None:
    """The cis signal's extreme, its sign and level tested against the LTR cores or trans contacts. A line read outside the
    terminal eighths (cis_trans_line) has no say on a call inside one: it cannot see there."""
    cis = (c.h or {}).get("cis")
    blind = bool(cis and cis.get("interior") and s.call is not None and in_end_zone(c.ev, s.call))
    ev = c.ev
    if cis and not c.valid and not blind and "te_hic_pair" not in s.flags and s.source == "sequence" and s.call is not None \
            and s.relation != "contested":
        if _vgap(ev, cis["peak"], s.call, True) <= _bins(ev, HIC_TOL_BINS):            # at the call: an independent line agrees
            if not c.amb:
                s.relation, s.rule = "concordant", "cis_agrees"
                s.flags.append("cis_agrees")
            elif "hic_alt_ambiguous" not in s.flags:
                s.flags.append("hic_alt_ambiguous")
        elif _vgap(ev, cis["peak"], s.call, False) <= _bins(ev, HIC_FAR_BINS):         # within Hi-C's resolution: neither way
            s.flags.append("cis_near")
        elif c.sconf == "low":               # elsewhere, the call unsure: the better candidate, still unsure
            s.call, s.relation, s.source, s.rule = cis["peak"], "offset", "hic", "cis_only"
            s.flags.append("cis_only")
        else:
            s.flags.append("cis_elsewhere")
            if cis["status"] == "pass":      # a line at the LTR cores on most chromosomes lies elsewhere: one line against one
                s.relation, s.rule = "contested", "cis_elsewhere"
    elif cis and s.call is None:
        s.call, s.relation, s.source, s.rule = cis["peak"], "functional_only", "hic", "cis_only"
        s.flags.append("cis_only")
    elif blind and not c.valid and s.source == "sequence" and _vgap(ev, cis["peak"], s.call, False) > _bins(ev, HIC_FAR_BINS):
        s.flags.append("cis_blind")          # said in the caveat: the extreme lies elsewhere, but the line cannot judge an end


def _at_places(c: _Ctx, s: _Call) -> None:
    """A chromosome nothing else places takes its AT trough: where the genome's troughs passed their test (at_only), or where
    the trough and a clean prior-free Hi-C peak that placed no call of its own agree, in a genome where the two agree beyond
    chance (0.13, at_denovo_pair; C. neoformans in anaphase: 3). Where the troughs passed their test and the sequence call
    is low and has a competing cluster, the trough's deepest window lying on that cluster and not on the call chooses it
    (0.15, at_chose_alt, as Hi-C does with hic_chose_alt): a tested line picks between the sequence line's own two candidates.
    Not where validated Hi-C has a peak on this chromosome and can tell the two apart: it has chosen already, or the call stands
    beside its peak. A heuristic, named as such (design.md): two satellites differ in composition whatever their centromeres do,
    so the deepest trough can lie on the inactive one; only T2T mouse exercises the rule, and the call stays low unless another
    line agrees at the chosen cluster."""
    ev, h = c.ev, c.h
    at_two = bool(ev.at and ev.at.get("support") and ev.at.get("agrees") and h and _clean_flags(h)
                  and h.get("denovo_confidence", "none") != "none")
    if s.call is None and (c.at_ok or at_two):
        s.call, s.relation, s.source = (ev.at["start"], ev.at["end"]), "sequence_only", "composition"
        s.flags.append("at_only" if c.at_ok else "at_denovo_pair")
        s.rule = s.flags[-1]
    elif c.at_ok and c.sconf == "low" and c.alt_iv and s.rule == "sequence" and (c.amb or not (c.valid and c.peak)):
        lo, hi = ev.at.get("trough_start", ev.at["start"]), ev.at.get("trough_end", ev.at["end"])    # the deepest window itself
        if gap(lo, hi, *c.alt_iv) == 0 and gap(lo, hi, *s.call) > 0:
            s.call, s.rule = c.alt_iv, "at_chose_alt"
            s.flags.append("at_chose_alt")


PHASES = (_te_hic_pair, _cis, _at_places)


# ---- after a call is placed: what refines or contests it, then its support, what is against it, and its confidence
def _refine(c: _Ctx, s: _Call) -> None:
    """The AT trough's finer interval for a Hi-C call it agrees with (at_pair) or a flag where it lies elsewhere; Hi-C alone's
    prior-free peak elsewhere (said); the Rabl-axis line and the prior-free peak opposing a kept call (contested)."""
    ev, h = c.ev, c.h
    if c.seq_iv is None and s.source == "hic" and ev.at:
        trough = (ev.at["start"], ev.at["end"])
        at_near = _vgap(ev, trough, s.call, True) <= _bins(ev, HIC_TOL_BINS)
    else:
        trough, at_near = None, False
    s.at_pair = bool(at_near and ev.at.get("agrees") and (ev.at.get("validated") or ev.at.get("support")))
    if s.at_pair:              # two independent lines on a chromosome without a sequence call: the trough is the finer
        s.call, s.rule = trough, "at_pair"                           # (a window of 1-500 kb, a Hi-C bin of 10-250 kb)
        s.flags.append("at_pair")
    elif trough and not at_near and ev.at.get("validated"):
        s.flags.append("at_elsewhere")    # said in the caveat; a Hi-C call stands against a trough of the same standing
    if c.untested and not c.ax and s.source == "sequence":
        dpk, dco = _iv(h, "denovo_peak_start", "denovo_peak_end"), _iv(h, "denovo_core_start", "denovo_core_end")
        if dpk and dco and _denovo_far(ev, dpk, dco, s.call):
            s.flags.append("hic_denovo_elsewhere")
    if c.ax and s.source == "sequence":
        s.f_opp = _fold_far(ev, c.f_core, s.call)
        s.d_opp = bool(c.dn_peak and c.dn_core and _denovo_far(ev, c.dn_peak, c.dn_core, s.call))
        if s.f_opp or s.d_opp:
            s.relation, s.rule = "contested", f"{c.axis}_offset" if s.f_opp else "hic_denovo_offset"
            s.flags += [x for x, y in ((f"{c.axis}_offset", s.f_opp), ("hic_denovo_offset", s.d_opp)) if y]
    elif c.ax and "hic_denovo_only" in s.flags and _fold_far(ev, c.f_core, s.call):
        s.f_opp = True                                     # the axis line passed genome-wide but points elsewhere
        s.flags.append(f"{c.axis}_offset")                 # here: said, not hidden (the call is low already)


def _support(c: _Ctx, s: _Call, m_call: dict) -> list:
    """The lines that agree with the call."""
    ev, flags = c.ev, s.flags
    support = ["at"] if s.at_pair or {"at_only", "at_denovo_pair", "at_chose_alt"} & set(flags) else []
    if ev.sat_bp and ev.sat_bp(*s.call) >= SAT_SUPPORT_BP:
        support.append("satellite")
    if (s.source == "sequence" and _te_agrees(ev, s.call)) or (s.source == "hic" and "te" in s.second) or "te_hic_pair" in flags:
        support.append("te")
    if s.source == "hic" and "age" in s.second:
        support.append("age")
    if c.valid and "hic_alt_ambiguous" not in flags and (_peak_agrees(ev, s.call) or
                                                         (s.source == "hic" and c.h.get("anch_level") != "none")):
        support.append("hic")
    elif {"hic_support", "hic_denovo_only", "te_hic_pair", "hic_end", "cis_agrees", "cis_only", "at_denovo_pair"} & set(flags):
        support.append("hic")
    if c.axis and f"{c.axis}_pair" in flags:
        support += ["hic", c.axis]
    elif c.ax and s.source == "sequence" and gap(*s.call, *c.f_core) == 0:
        support.append(c.axis)
    if m_call["label"] == "centromere_like":
        support.append("meth")
    return support


def _against(c: _Ctx, s: _Call) -> list:
    """The lines that point elsewhere."""
    ev, flags = c.ev, s.flags
    against = []
    if c.valid and c.strong and _peak_far(ev, s.call):
        against.append("hic")
    if _te_decides(ev) and not (_te_agrees(ev, s.call) if s.source == "sequence" else "te" in s.second or "te_hic_pair" in flags):
        against.append("te")
    if "at_elsewhere" in flags:
        against.append("at")
    if "cis_elsewhere" in flags:
        against.append("hic")
    against += [x for x, y in ((c.axis, s.f_opp), ("hic", s.d_opp)) if y]
    return against


def _confidence(c: _Ctx, s: _Call, support: list) -> str:
    """One level up per agreeing tested line (at most MAX_RAISE), down for a contest; a moved or functional-only call is
    medium, high with a second line; a trough alone low, medium where Hi-C validated the troughs."""
    if s.source == "sequence":
        if s.relation == "contested":
            conf = level_shift(c.sconf, -max(1, int(s.f_opp) + int(s.d_opp)) + int("meth" in support))
        else:
            conf = level_shift(c.sconf, min(MAX_RAISE, int("hic" in support) + int("meth" in support)))
        if s.relation == "sequence_only" and "meth" in support:
            s.relation = "concordant"
        return conf
    if s.source == "composition":
        return "medium" if c.ev.at.get("validated") else "low"
    if "cis_only" in s.flags:
        return "low"
    if "hic_medium_peak" in s.flags or "hic_denovo_only" in s.flags:
        return "medium" if s.at_pair else "low"
    return "high" if s.second or s.at_pair else "medium"


def _decide(ev: Evidence) -> dict:
    """The rules in their declared order: the first rule of CASCADE that applies sets the call; the PHASES then act on it in
    turn; then the call is refined or contested (_refine), and its support, what is against it and its confidence follow. The
    order is the specification (design.md, Integration): a later reading may override an earlier one, as the LTR-core and cis
    readings override a plain sequence call, and the AT trough places what nothing else did or chooses a low call's competing
    cluster. `rule` names the last rule
    that set the call's place or relation."""
    c, s = _context(ev)
    next(r for r in CASCADE if r(c, s))
    for phase in PHASES:
        phase(c, s)
    if s.call is None:
        return dict(call=None, relation="none", confidence="none", source="none", support=[], against=[],
                    flags=s.flags, rule="none", second=s.second, meth_call={"label": "", "score": np.nan})
    _refine(c, s)
    # methylation is read where the call now is: its own item (the provisional call) where that is not the sequence call or Hi-C's
    # peak, as at the prior-free core or at a competing cluster that was chosen (the verdict on the abandoned call says nothing)
    at_core = "hic_denovo_only" in s.flags or (c.axis and f"{c.axis}_pair" in s.flags) or "at_chose_alt" in s.flags or c.chose
    m_call = _meth(ev, "provisional") if at_core else c.m_hic if s.source == "hic" else c.m_seq
    support = _support(c, s, m_call)
    against = _against(c, s)
    conf = _confidence(c, s, support)
    if s.rule == "sequence" and s.source == "sequence" and "meth" in support and s.relation == "concordant":
        s.rule = "meth_support"                # methylation alone raised or confirmed a kept call
    return dict(call=s.call, relation=s.relation, confidence=conf, source=s.source, support=support, against=against,
                flags=s.flags, rule=s.rule, second=s.second, meth_call=m_call)


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
        flags=",".join(d["flags"]) or na, rule=d.get("rule", na), offset_bp=off,
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
    the same Hi-C; Rabl nuclei only) or tes/ (Hi-C tested on the LTR cores); without one the sequence calls are kept. Seen
    both ways: barley (satellite line wrong, prior-free Hi-C right; the fold sides with Hi-C) and S. orangutan /
    giraffe (the reverse; no fold)."""
    s = t[t.seq_start.ge(0) & t.callable]
    n = len(s)
    conf = int(s.denovo_confidence.isin(["high", "medium"]).sum())
    agree = int(s.denovo_agrees_seq.sum())
    gate = int(s.gate_pass.sum())
    flag = bool(n >= GATE_MIN_REF and 2 * conf > n and 2 * agree < n and 2 * gate < n)
    return dict(conflict=flag, n=n, denovo_confident=conf, denovo_on_seq=agree, gate_pass=gate)
