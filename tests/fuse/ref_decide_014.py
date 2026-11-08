"""The integration cascade of cenfuse 0.14.0 as one function (fuse._decide before it was split into rule functions; future.md #63),
kept unchanged as the reference that tests/fuse/test_fuse_equivalence.py compares the split rules against."""
import numpy as np

from cenfuse.ageline import AGE_ALPHA
from cenfuse.fuse import (Evidence, _bins, _clean_denovo, _clean_flags, _clean_medium, _denovo_far, _fold_far, _iv, _meth,
                          _peak_agrees, _peak_far, _te_agrees, _te_decides, _te_near, _vgap, axis_line, hic_core, hic_peak,
                          in_end_zone)
from cenfuse.hicfuse import SUPPORT_LEVELS
from cenfuse.params import HIC_FAR_BINS, HIC_TOL_BINS, MAX_RAISE, SAT_SUPPORT_BP
from cenfuse.util import gap, level_shift


def decide_014(ev: Evidence) -> dict:
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
    chose = "hic_chose_alt" in flags
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
        rule = f"{axis}_pair"
        if sconf == "high":
            flags.append(f"moved_with_{axis}")
    elif seq_iv is None and strong:
        call, relation, source, rule = core, "functional_only", "hic", "hic_strong_peak"
    elif seq_iv is None and valid and core and _clean_medium(h):
        call, relation, source, rule = core, "functional_only", "hic", "hic_medium_peak"
        flags.append("hic_medium_peak")
    elif seq_iv is None and dn_only:
        call, relation, source, rule = dn_only, "functional_only", "hic", "hic_denovo_only"
        flags.append("hic_denovo_only")
    elif seq_iv is None:
        call, relation, rule = None, "none", "none"
    elif agree and amb:
        call, rule = seq_iv, "hic_chose_alt" if chose else "sequence"
        flags.append("hic_alt_ambiguous")
    elif agree:
        call, relation, rule = seq_iv, "concordant", "hic_chose_alt" if chose else "hic_agrees"
    elif weak and _peak_agrees(ev, seq_iv) and h.get("anch_level") in SUPPORT_LEVELS:
        call, rule = seq_iv, "hic_chose_alt" if chose else "sequence"
        if amb:
            flags.append("hic_alt_ambiguous")
        else:
            relation = "concordant"
            flags.append("hic_support")
            rule = "hic_chose_alt" if chose else "hic_support"
    elif not valid and not ax and h and h.get("end_pass") and h.get("end_higher") and not amb:
        call, relation, rule = seq_iv, "concordant", "hic_end"   # Hi-C tells the call's end from the other, where it does so
        flags.append("hic_end")
    elif strong and far and (sconf != "high" or second):
        call, relation, source, rule = core, "offset", "hic", "hic_moved"
        if sconf == "high":
            flags.append("moved_with_" + "+".join(second))
    elif strong and far:
        call, relation, rule = seq_iv, "contested", "hic_offset"
        flags.append("hic_offset")
    else:
        call, rule = seq_iv, "hic_chose_alt" if chose else "sequence"
        if valid and peak:
            flags.append("hic_near" if not far else "hic_disagrees_weak")
    pair = (h or {}).get("te_pair")      # Hi-C's peak at this chromosome's LTR core, the two lines tested genome-wide
    # both readings run where the gate fails on most chromosomes; a chromosome whose own gate passed keeps what validated Hi-C said
    if pair and not valid and source == "sequence" and call is not None and relation != "contested":
        if gap(*call, *pair["core"]) == 0:                                           # both at the call: Hi-C supports it
            if not amb:                                                              # no mappable bin between it and its competitor
                relation, rule = "concordant", "te_hic_pair"
                flags.append("te_hic_pair")
            elif "hic_alt_ambiguous" not in flags:
                flags.append("hic_alt_ambiguous")
        elif sconf == "low" or pair["status"] == "pass":     # two lines at one locus against one line, unsure or on a minority
            call, relation, source, rule = pair["core"], "offset", "hic", "te_hic_pair"
            flags.append("te_hic_pair")
            if sconf == "high":
                second.append("te")
                flags.append("moved_with_te")
    elif pair and call is None:
        call, relation, source, rule = pair["core"], "functional_only", "hic", "te_hic_pair"
        flags.append("te_hic_pair")
    cis = (h or {}).get("cis")           # the cis signal's extreme, its sign and level tested against the LTR cores or trans contacts
    # a line read outside the terminal eighths (cis_trans_line) has no say on a call inside one: it cannot see there
    blind = bool(cis and cis.get("interior") and call is not None and in_end_zone(ev, call))
    if cis and not valid and not blind and "te_hic_pair" not in flags and source == "sequence" and call is not None and \
            relation != "contested":
        if _vgap(ev, cis["peak"], call, True) <= _bins(ev, HIC_TOL_BINS):              # at the call: an independent line agrees
            if not amb:
                relation, rule = "concordant", "cis_agrees"
                flags.append("cis_agrees")
            elif "hic_alt_ambiguous" not in flags:
                flags.append("hic_alt_ambiguous")
        elif _vgap(ev, cis["peak"], call, False) <= _bins(ev, HIC_FAR_BINS):           # within Hi-C's resolution: neither way
            flags.append("cis_near")
        elif sconf == "low":                 # elsewhere, the call unsure: the better candidate, still unsure
            call, relation, source, rule = cis["peak"], "offset", "hic", "cis_only"
            flags.append("cis_only")
        else:
            flags.append("cis_elsewhere")
            if cis["status"] == "pass":      # a line at the LTR cores on most chromosomes lies elsewhere: one line against one
                relation, rule = "contested", "cis_elsewhere"
    elif cis and call is None:
        call, relation, source, rule = cis["peak"], "functional_only", "hic", "cis_only"
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
        rule = flags[-1]
    if call is None:
        return dict(call=None, relation="none", confidence="none", source="none", support=[], against=[],
                    flags=flags, rule="none", second=second, meth_call={"label": "", "score": np.nan})
    if seq_iv is None and source == "hic" and ev.at:
        trough = (ev.at["start"], ev.at["end"])
        at_near = _vgap(ev, trough, call, True) <= _bins(ev, HIC_TOL_BINS)
    else:
        trough, at_near = None, False
    at_pair = bool(at_near and ev.at.get("agrees") and (ev.at.get("validated") or ev.at.get("support")))
    if at_pair:                  # two independent lines on a chromosome without a sequence call: the trough is the finer
        call, rule = trough, "at_pair"                               # (a window of 1-500 kb, a Hi-C bin of 10-250 kb)
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
            relation, rule = "contested", f"{axis}_offset" if f_opp else "hic_denovo_offset"
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
    if rule == "sequence" and source == "sequence" and "meth" in support and relation == "concordant":
        rule = "meth_support"                # methylation alone raised or confirmed a kept call
    return dict(call=call, relation=relation, confidence=conf, source=source, support=support, against=against,
                flags=flags, rule=rule, second=second, meth_call=m_call)


