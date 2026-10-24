"""Integration rules (docs/design.md, 'provisional' and 'fuse').

The integrated call is the sequence call unless a validated functional line moves it:
- Hi-C may act on chromosome x only when its gate passes (hicfuse.gate), and only through a strong anchored peak
  (all three Hi-C engine tests, no hub/rDNA/telomeric flag). A peak supports a call within HIC_TOL_BINS, only when it is a
  call at all (Hi-C level not 'none').
- A strong peak that points to another locus (more than HIC_FAR_BINS from the sequence call, the method's resolution,
  and with a core that does not reach the call) replaces a low or medium sequence call, and places a call on a
  chromosome without one. Closer than that, Hi-C neither supports nor opposes the call (flag hic_near).
- It replaces a high sequence call only with a second line at the Hi-C peak and not at the sequence call: a decisive
  (CRM) TE core within HIC_FAR_BINS of the peak and beyond it from the call, or methylation `centromere_like` at the
  peak with a measured non-centromere-like `dip` at the call (no data is not evidence). Without one the call stays and
  is `contested`.
- Methylation never moves a call alone.
- Competing clusters (the satellite engine's second cluster, at least half the main one): validated Hi-C that agrees
  with the second cluster and not the main one chooses it (flag hic_chose_alt). Hi-C chooses between two sequence
  candidates; it does not place the call.
- No sequence call, gate passed, a clean (no hub/rDNA/telomeric) anchored peak at Hi-C level medium: a low-confidence
  functional_only call (flag hic_medium_peak). Horse chr11, the satellite-free centromere, was found this way and
  discarded by the rule that asked for level high.

Every distance that involves a Hi-C interval is counted in mappable bins (hicfuse.visible_gap): Hi-C cannot see into
unmappable arrays, so a peak at an array's edge is adjacent to a call inside it. Agreement discounts at most
HIC_FAR_BINS unmappable bins; opposition needs that many mappable bins. Across a longer blind stretch Hi-C cannot
place the centromere on either side, so it neither supports nor opposes.

Confidence uses the satellite engine's levels, changed one step at a time (its te_concordant convention):
- A kept sequence call goes up one level for agreeing validated Hi-C and up one for methylation support, at most
  MAX_RAISE in total. Contested: down one level, offset by methylation support.
- A moved or functional-only call is medium, and high with a second line."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
import pandas as pd

from .hicfuse import visible_gap
from .params import GATE_MIN_REF, HIC_FAR_BINS, HIC_TOL_BINS, MAX_RAISE, NOT_STRONG_FLAGS, SAT_SUPPORT_BP
from .sat.te import agrees
from .util import gap, level_shift

FUSED_COLS = ["chrom", "chrom_len", "start", "end", "length", "call_source", "relation", "confidence", "support",
              "against", "flags", "offset_bp", "seq_start", "seq_end", "seq_evidence", "seq_confidence", "seq_flags",
              "seq_alt_start", "seq_alt_end", "main_family", "te_clade", "te_start", "te_end", "hic_gate", "hic_gate_k", "hic_gate_n", "hic_gate_p",
              "hic_anch_peak_start", "hic_anch_peak_end", "hic_anch_core_start", "hic_anch_core_end", "hic_anch_level",
              "hic_anch_flags", "hic_anch_agrees_seq", "hic_denovo_peak_start", "hic_denovo_peak_end",
              "hic_denovo_confidence", "hic_denovo_flags", "meth_seq", "meth_seq_score", "meth_hic", "meth_hic_score",
              "meth_call", "meth_call_score"]


@dataclass
class Evidence:
    """Everything known about one chromosome (0-based half-open intervals; -1 = absent)."""
    chrom: str
    chrom_len: int
    seq: dict                 # start, end, evidence, confidence, flags, main_family; alt_start/alt_end (competing
                              # cluster, -1 = none)
    te: dict | None           # clade, start, end, decisive (the marker clade's core)
    hic: dict | None          # one row of the hicfuse table, or None (no Hi-C, or not a callable chromosome)
    meth: dict | None = None  # item -> {'label', 'score'} pooled over samples (methfuse.support), or None
    sat_bp: Callable | None = None   # (start, end) -> bp of centromeric-role satellite arrays inside


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
    """A significant anchored peak (Hi-C level not 'none') within HIC_TOL_BINS mappable bins of iv."""
    pk = hic_peak(ev)
    return bool(pk and iv and ev.hic.get("anch_level", "none") != "none"
                and _vgap(ev, pk, iv, True) <= _bins(ev, HIC_TOL_BINS))


def _peak_far(ev: Evidence, iv) -> bool:
    """The peak points to another locus: farther from iv than the method's resolution (HIC_FAR_BINS mappable bins),
    and its core does not reach iv."""
    pk, core = hic_peak(ev), hic_core(ev)
    return bool(pk and iv and _vgap(ev, pk, iv, False) > _bins(ev, HIC_FAR_BINS)
                and (core is None or _vgap(ev, core, iv, False) > _bins(ev, HIC_TOL_BINS)))


def _te_agrees(ev: Evidence, iv) -> bool:
    """The satellite engine's rule (gap at most the wider interval's width): used for kept sequence calls, as that engine does."""
    return bool(ev.te and iv and agrees(iv[0], iv[1], dict(start=ev.te["start"], end=ev.te["end"])))


def _te_near(ev: Evidence, iv, n: int) -> bool:
    """The TE core within n mappable Hi-C bins of iv (one fixed tolerance for both places a move compares)."""
    return bool(ev.te and iv and _vgap(ev, (ev.te["start"], ev.te["end"]), iv, False) <= _bins(ev, n))


def _clean_medium(h) -> bool:
    """An anchored peak at Hi-C level medium without a hub, rDNA or telomeric flag."""
    return bool(h and h.get("anch_level") == "medium"
                and not set(str(h.get("anch_flags") or "").split(",")) & set(NOT_STRONG_FLAGS))


def _meth(ev: Evidence, item: str) -> dict:
    return (ev.meth or {}).get(item) or {"label": "", "score": np.nan}


def decide(ev: Evidence) -> dict:
    """One chromosome's integrated call (see module doc). With ev.meth None it is the provisional call."""
    seq_iv = _iv(ev.seq)
    sconf = ev.seq["confidence"] if seq_iv else "none"
    alt_iv = _iv(ev.seq, "alt_start", "alt_end") if seq_iv and "alt_start" in ev.seq else None
    h = ev.hic
    valid = bool(h and h.get("gate_pass"))
    core, peak = hic_core(ev), hic_peak(ev)
    strong = bool(valid and h.get("anch_strong") and core)
    agree = bool(valid and seq_iv and _peak_agrees(ev, seq_iv))
    flags, relation, source = [], "sequence_only", "sequence"
    if valid and alt_iv and not agree and _peak_agrees(ev, alt_iv):     # Hi-C picks the competing cluster
        seq_iv, agree = alt_iv, True
        flags.append("hic_chose_alt")
    far = bool(seq_iv and _peak_far(ev, seq_iv))
    m_seq, m_hic = _meth(ev, "sequence_call"), _meth(ev, "hic_peak")
    second = []
    if strong:
        decisive = bool(ev.te and ev.te.get("decisive"))
        if decisive and _te_near(ev, peak, HIC_FAR_BINS) and not (seq_iv and _te_near(ev, seq_iv, HIC_FAR_BINS)):
            second.append("te")
        if m_hic["label"] == "centromere_like" and (seq_iv is None or m_seq["label"] == "dip"):
            second.append("meth")
    if seq_iv is None and strong:
        call, relation, source = core, "functional_only", "hic"
    elif seq_iv is None and valid and core and _clean_medium(h):
        call, relation, source = core, "functional_only", "hic"
        flags.append("hic_medium_peak")
    elif seq_iv is None:
        call, relation = None, "none"
    elif agree:
        call, relation = seq_iv, "concordant"
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
    if call is None:
        return dict(call=None, relation="none", confidence="none", source="none", support=[], against=[],
                    flags=flags, second=second, meth_call={"label": "", "score": np.nan})
    m_call = m_hic if source == "hic" else m_seq
    support = []
    if ev.sat_bp and ev.sat_bp(*call) >= SAT_SUPPORT_BP:
        support.append("satellite")
    if (source == "sequence" and _te_agrees(ev, call)) or (source == "hic" and "te" in second):
        support.append("te")
    if valid and (_peak_agrees(ev, call) or (source == "hic" and h.get("anch_level") != "none")):
        support.append("hic")
    if m_call["label"] == "centromere_like":
        support.append("meth")
    against = []
    if valid and strong and _peak_far(ev, call):
        against.append("hic")
    if ev.te and ev.te.get("decisive") and not (_te_agrees(ev, call) if source == "sequence" else "te" in second):
        against.append("te")
    if source == "sequence":
        if relation == "contested":
            conf = level_shift(sconf, -1 + int("meth" in support))
        else:
            raise_by = int("hic" in support) + int("meth" in support)
            conf = level_shift(sconf, min(MAX_RAISE, raise_by))
        if relation == "sequence_only" and "meth" in support:
            relation = "concordant"
    elif "hic_medium_peak" in flags:
        conf = "low"
    else:
        conf = "high" if second else "medium"
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
    return dict(
        chrom=ev.chrom, chrom_len=ev.chrom_len, start=one(call[0]) if call else na, end=end(call[1]) if call else na,
        length=call[1] - call[0] if call else na, call_source=d["source"], relation=d["relation"],
        confidence=d["confidence"], support=",".join(d["support"]) or na, against=",".join(d["against"]) or na,
        flags=",".join(d["flags"]) or na, offset_bp=off,
        seq_start=one(seq_iv[0]) if seq_iv else na, seq_end=end(seq_iv[1]) if seq_iv else na,
        seq_evidence=ev.seq["evidence"], seq_confidence=ev.seq["confidence"], seq_flags=ev.seq["flags"],
        seq_alt_start=one(alt_iv[0]) if alt_iv else na, seq_alt_end=end(alt_iv[1]) if alt_iv else na,
        main_family=ev.seq["main_family"], te_clade=ev.te["clade"] if ev.te else na,
        te_start=one(ev.te["start"]) if ev.te else na, te_end=end(ev.te["end"]) if ev.te else na,
        hic_gate=gate_status(h) if h else na,
        hic_gate_k=h.get("gate_k", na) if h else na, hic_gate_n=h.get("gate_n", na) if h else na,
        hic_gate_p=h.get("gate_p", na) if h else na,
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
        meth_seq=m["sequence_call"]["label"] or na, meth_seq_score=m["sequence_call"]["score"],
        meth_hic=m["hic_peak"]["label"] or na, meth_hic_score=m["hic_peak"]["score"],
        meth_call=m["call"]["label"] or na, meth_call_score=m["call"]["score"])


def gate_status(h: dict) -> str:
    """pass / fail / untested (fewer than GATE_MIN_REF reference chromosomes: no majority test was possible)."""
    if h.get("gate_pass"):
        return "pass"
    n = h.get("gate_n")
    return "untested" if n is None or (isinstance(n, (int, float, np.integer)) and int(n) < GATE_MIN_REF) else "fail"


def table(evs: list, decisions: list) -> pd.DataFrame:
    return pd.DataFrame([row(e, d) for e, d in zip(evs, decisions)], columns=FUSED_COLS)


def sequence_conflict(t: pd.DataFrame) -> dict:
    """Genome-level diagnostic. Hi-C's prior-free calls are confident (Hi-C high/medium) on most chromosomes but agree
    with the sequence calls on fewer than half, and the gate fails on most: the two lines disagree genome-wide, so at
    least one is wrong everywhere, and two lines cannot tell which. cenfuse then keeps the sequence calls. A third line
    decides: tes/ (CRM cores) or methylation. Seen in the benchmark both ways: barley (satellite line wrong, prior-free
    Hi-C right) and S. orangutan / giraffe (the reverse)."""
    s = t[t.seq_start.ge(0) & t.callable]
    n = len(s)
    conf = int(s.denovo_confidence.isin(["high", "medium"]).sum())
    agree = int(s.denovo_agrees_seq.sum())
    gate = int(s.gate_pass.sum())
    flag = bool(n >= GATE_MIN_REF and 2 * conf > n and 2 * agree < n and 2 * gate < n)
    return dict(conflict=flag, n=n, denovo_confident=conf, denovo_on_seq=agree, gate_pass=gate)
