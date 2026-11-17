"""What the overview figure says about each call, in plain words, read from `<genome>.centromeres.tsv` and the genome's
`cenfuse_summary.tsv` row: where each signal alone puts the centromere, its role in the call, one sentence on how the
call was made, and each signal's genome-wide test. Pure logic: nothing here draws or changes a call.

A signal is one readout on its own: the satellite arrays, the LTR family's cluster, each Hi-C readout (contact with the
other centromeres, the contact hub from Hi-C alone, the contact between flanks, the arm fold, the arm alignment), CHG
and CpG methylation, the AT-rich trough, the repeat density and the repeat-family cluster. The contact with the other
centromeres has a place only where cenfuse records a peak (hic_anch_level not none). Its role in a call, first match
wins: 1 listed in `support` -> used; 2 listed in `against` -> points elsewhere; 3 failed its genome test (or
methylation measured in other plants) -> not used; 4 the satellite arrays or the LTR family's cluster, having passed
its test, has a place more than one Hi-C bin (250 kb without Hi-C) from the call -> points elsewhere: the place a moved
or overruled call leaves behind, which cenfuse does not list in `against`; 5 -> not used. Every other signal points
elsewhere only where cenfuse counted it against the call (rule 2): one far from the call that cenfuse did not count (a
weak or telomeric Hi-C peak, say) says nothing against it. Without a summary row the genome tests are unknown: rules
1, 2 and 5 only. The verdict follows the branches of interpret.basis, in its order; the same branches name the signals
that placed the call (Row.placed)."""
from __future__ import annotations

import re
from dataclasses import dataclass

import pandas as pd

NO_HIC_TOL = 250_000            # prior: without Hi-C, a signal whose own place lies farther than this from the call
                                # points elsewhere (with Hi-C: one bin); one bin at the usual Hi-C resolution
SAT_KINDS = ("satellite", "weak_satellite", "long_unit_tandem")
LEFT_BEHIND = frozenset({"sat", "te"})      # rule 4: the places a moved or overruled call leaves behind
NO_PEAK = (".", "none", "", "nan", "None")  # hic_anch_level where cenfuse records no Hi-C peak
KIND = {"satellite": "satellite array", "weak_satellite": "small satellite array",
        "long_unit_tandem": "long tandem repeat"}
# seq_flags of a satellite call an LTR core replaced (the array is in seq_overruled_*)
OVERRULED = frozenset({"satellite_away_from_core", "satellite_beside_core", "minority_core", "weak_satellite_beside_core"})
# the Hi-C readout the token `hic` in `support` names, by rule or flag (fuse._support); `against`: fuse._against
CIS_SUPPORT = frozenset({"cis_agrees", "cis_only", "cis_trans_pair", "cis_end_pair"})
HUB_SUPPORT = frozenset({"hic_denovo_only", "fold_pair", "align_pair", "at_denovo_pair"})
TWO_HIC = frozenset({"cis_trans_pair", "cis_end_pair", "fold_pair", "align_pair"})   # two Hi-C readouts place the call
HIC_KEYS = frozenset({"hic", "hub", "flank", "fold", "align"})
# the secondary Hi-C readouts and the repeat-family cluster (a lead: its genome test is uncalibrated), consulted in some
# genomes only, and the rules and flags that show cenfuse consulted one on a chromosome (besides its token in support or
# against): only those cenfuse consulted in a genome get a close-up track
ACTED = {"hub": frozenset({"hic_denovo_only", "at_denovo_pair", "fold_pair", "align_pair", "hic_denovo_offset",
                           "hic_denovo_elsewhere"}),
         "flank": frozenset({"cis_agrees", "cis_only", "cis_elsewhere", "cis_trans_pair", "cis_end_pair", "cis_near",
                             "cis_blind"}),
         "fold": frozenset({"fold_pair", "fold_offset", "moved_with_fold"}),
         "align": frozenset({"align_pair", "align_offset", "moved_with_align"}),
         "fam": frozenset({"fam_only"})}
# the order of the signals' tracks, top to bottom, in the figure: the Hi-C readouts (the Hi-C trace first), methylation,
# the LTR family, the satellite arrays, then the FASTA's composition signals
TRACK_ORDER = ("hic", "flank", "fold", "align", "hub", "meth", "cpg", "te", "sat", "at", "rep", "fam")
ROLE_WORDS = {"used": "used", "elsewhere": "points elsewhere", "not": "not used"}
ELSE, GRAY = "#D55E00", "#A6A6A6"           # the ink of a role, in every panel: points elsewhere; not used
ROLE_INK = {"used": "black", "elsewhere": ELSE, "not": GRAY}
SAT_QUESTION = "Largest arrays one per chromosome"
VERDICTS = ("trusted", "not trusted here", "a lead only", "not used (other plants)", "not reported")
# words that are cenfuse's jargon, not a biologist's: never on a figure page (whole words, plurals too)
BANNED = re.compile(r"\b(?:sequence call|integrated|anchored|prior-free|de novo|concordant|offset|contested|functional|"
                    r"gate|cis|signature|line|rabl)s?\b", re.I)


@dataclass(frozen=True)
class Signal:
    key: str
    name: str                       # track label; {clade} = the LTR family's lineage
    color: str                      # as plot.py draws the signal
    cols: tuple[str, str] | None    # its own place: start and end columns (1-based closed) of asm.centromeres.tsv
    token: str | None               # its name in `support` / `against`
    test: tuple[str, ...]           # the columns of its genome test (summary row; hic_gate and cpg_gate: the call's row)
    short: str                      # in a verdict's list of agreeing signals
    where: str | None               # in a verdict's tail: "{where} 2.3 Mb away disagrees"

    def label(self, clade: str = "LTR") -> str:
        return self.name.format(clade=clade)


SIGNALS = {s.key: s for s in (
    Signal("sat", "satellite arrays", "#D55E00", ("seq_start", "seq_end"), "satellite", ("sat_one_per_chromosome",),
           "satellite", "satellite array"),
    Signal("te", "{clade} elements", "#E69F00", ("te_start", "te_end"), "te", ("te_decisive", "te_verdict"),
           "{clade}", "{clade} cluster"),
    Signal("hic", "Hi-C contact with other centromeres", "#0072B2", ("hic_anch_core_start", "hic_anch_core_end"),
           "hic", ("hic_gate_pass", "hic_gate_support", "hic_gate_untested", "hic_chromosomes", "te_hic", "hic_gate"),
           "Hi-C", "Hi-C peak"),
    Signal("hub", "Hi-C contact hub, from Hi-C alone", "#0072B2", ("hic_denovo_peak_start", "hic_denovo_peak_end"),
           "hic", (), "Hi-C", "Hi-C hub"),
    Signal("flank", "Hi-C contact between flanks", "#6A3D9A", ("cis_peak_start", "cis_peak_end"), "hic",
           ("cis_line",), "Hi-C", "Hi-C flank contact"),
    Signal("fold", "Hi-C: arms fold back here", "#0072B2", ("fold_core_start", "fold_core_end"), "fold",
           ("fold_gate", "fold_on_denovo", "fold_p"), "Hi-C", "Hi-C arm fold"),
    Signal("align", "Hi-C: arms align", "#0072B2", ("align_core_start", "align_core_end"), "align",
           ("align_gate", "align_on_denovo", "align_p"), "Hi-C", "Hi-C arm alignment"),
    Signal("meth", "CHG methylation", "#303030", None, "meth",
           ("inputs_reference", "inputs_unknown", "inputs_other", "meth_signature"), "CHG methylation", None),
    Signal("cpg", "CpG methylation dip", "#303030", ("cpg_kt_start", "cpg_kt_end"), "cpg", ("cpg_gate_pass", "cpg_gate"),
           "CpG methylation dip", "CpG dip"),
    Signal("at", "AT-rich trough", "#009E73", ("at_start", "at_end"), "at", ("at_gate", "at_p", "at_hic"),
           "AT-rich trough", "AT-rich trough"),
    Signal("rep", "repeat density", "#882255", ("rep_start", "rep_end"), "repeats", ("rep_gate", "rep_p"),
           "repeat density", "repeat-density peak"),
    Signal("fam", "repeat-family cluster", "#882255", ("fam_start", "fam_end"), "families", ("fam_gate", "fam_p"),
           "repeat-family cluster", "repeat-family cluster"),
)}
TOKEN = {s.token: s.key for s in SIGNALS.values() if s.token != "hic"}


@dataclass
class Row:
    """One chromosome: its call (0-based half-open, None without one), each signal's own place and role, the second
    satellite cluster, the verdict, the signals that placed the call (placed_by), the secondary Hi-C readouts and
    repeat-family cluster cenfuse consulted on it (ACTED), and whether Hi-C's contact with the other centromeres was
    scored here, with a peak or not (present: which Hi-C trace the figure draws)."""
    chrom: str
    L: int
    call: tuple[int, int] | None
    conf: str
    places: dict
    roles: dict
    alt: tuple[int, int] | None
    verdict: str = ""
    consulted: frozenset = frozenset()
    placed: frozenset = frozenset()
    hic_scored: bool = False


@dataclass
class TestRow:
    """One signal's genome test: what it asks, in plain words with its numbers; k of n chromosomes, the number chance
    gives, P (each None where the summary does not say); the verdict (one of VERDICTS); per_chrom where k counts the
    chromosomes whose own test passed (Hi-C, CpG), not those the signal marks."""
    __test__ = False                # a data class, not a test class, where a test module imports it
    key: str
    label: str
    question: str
    k: int | None
    n: int | None
    chance: float | None
    p: float | None
    verdict: str
    per_chrom: bool = False


# ------------------------------------------------------------------------------------------------ parsing
def _iv(a, b) -> tuple[int, int] | None:
    """1-based closed (start, end) -> 0-based half-open, as plot._iv (not imported: plot loads matplotlib); None where
    either is missing."""
    try:
        return int(a) - 1, int(b)
    except (TypeError, ValueError):
        return None


def _gap(a, b) -> int:
    """bp between two intervals' nearest edges (0 where they overlap)."""
    return max(0, max(a[0], b[0]) - min(a[1], b[1]))


def _toks(x) -> list:
    return [t for t in str(x).split(",") if t not in ("", ".", "nan", "None")]


def _num(x, kind=float):
    try:
        return kind(x)
    except (TypeError, ValueError):
        return None


def _s(summary: dict, col: str) -> str:
    v = summary.get(col)
    return "." if v is None or str(v) in ("", "nan", "None") else str(v)


_KN = re.compile(r"(\d+)/(\d+)")                      # '7/12 chromosomes', '9/14'
_CHANCE = re.compile(r"\(null ([\d.]+)")              # '(null 8, P 0.87)'
_P = re.compile(r"\bP[ =]([\d.]+(?:e[+-]?\d+)?)")     # 'P 0.87', 'P 2.4e-05', 'P=0.202'
_STATUS = re.compile(r": (\w+)")                      # '...: fail', '...: pass, tested by trans contacts'


def _kn(s: str):
    m = _KN.search(s)
    return (int(m[1]), int(m[2])) if m else (None, None)


def _p(s: str):
    m = _P.search(s)
    return float(m[1]) if m else None


# ------------------------------------------------------------------------------------------------ genome tests
def _status(key: str, s: dict) -> str | None:
    """The genome test of a signal from the summary row: 'pass', 'support', 'fail', 'untested', 'lead' (the
    repeat-family test passed: not calibrated), 'other plants' (methylation from other plants only), or None (the signal
    did not run, or has no test)."""
    if key == "sat":
        x = _s(s, "sat_one_per_chromosome")
        m = _STATUS.search(x)
        return m[1] if m else "untested"                # 'not tested on N sequences', or no centromeric family ('.')
    if key == "te":
        d = _s(s, "te_decisive")
        return "pass" if d == "True" else "fail" if d == "False" or "te" in _toks(s.get("lines")) else None
    if key == "hic":
        n = _num(s.get("hic_chromosomes"), int)
        if n is None:
            return None
        te_hic = _STATUS.search(_s(s, "te_hic"))       # Hi-C anchored on the LTR cores, where the gate failed
        if _num(s.get("hic_gate_pass"), int):
            return "pass"
        if _num(s.get("hic_gate_support"), int):
            return "support"
        if te_hic and te_hic[1] in ("pass", "support"):
            return te_hic[1]
        return "untested" if _num(s.get("hic_gate_untested"), int) == n else "fail"
    if key == "flank":
        m = _STATUS.search(_s(s, "cis_line"))
        return m[1] if m else None                      # 'ambiguous' (both signs pass) is no verdict: untrusted
    if key in ("fold", "align", "rep"):
        x = _s(s, f"{key}_gate")
        return None if x == "." else x
    if key == "fam":
        x = _s(s, "fam_gate")
        return None if x == "." else "lead" if x == "pass" else x
    if key == "at":                                     # the troughs pass alone or where Hi-C validated them (fuse at_ok)
        x = _s(s, "at_gate")
        return None if x == "." else "pass" if x == "pass" or _s(s, "at_hic").endswith(" pass") else x
    if key == "cpg":
        k = _num(s.get("cpg_gate_pass"), int)
        return None if k is None else "pass" if k else "fail"
    if key == "meth":                                   # reference and unknown samples decide (docs/plants.md)
        mine = [t for c in ("inputs_reference", "inputs_unknown") for t in _toks(s.get(c)) if t.startswith("meth:")]
        other = [t for t in _toks(s.get("inputs_other")) if t.startswith("meth:")]
        if mine or (not other and _s(s, "meth_samples") != "."):
            return "fail" if _s(s, "meth_signature") == "none" else "pass"
        return "other plants" if other else None
    return None                                         # hub: no genome test of its own


def _trust(key: str, s: dict, rec: dict) -> bool | None:
    """True where the signal passed its genome test (on this chromosome for Hi-C and CpG, whose tests leave it out),
    False where it failed or comes from other plants, None where untested, untestable or a lead."""
    if key in ("hic", "cpg"):
        x = str(rec.get(f"{key}_gate", "."))
        st = x if x in ("pass", "support", "fail") else None
    else:
        st = _status(key, s)
    return True if st in ("pass", "support") else False if st in ("fail", "other plants") else None


def signal_tests(summary: dict | None, n_chrom: int, clade: str) -> list[TestRow]:
    """Each signal's genome test, in registry order, for the signals that ran. Without a summary row: the satellite
    arrays alone, 'not reported'."""
    s = summary or {}
    if not s:
        return [TestRow("sat", SIGNALS["sat"].label(clade), SAT_QUESTION, None, None, None, None, "not reported")]
    out = []
    for key in SIGNALS:
        st = _status(key, s)
        if st is None and key != "sat":
            continue
        k = n = chance = p = None
        per_chrom = False                                # k counts chromosomes whose own test passed (Hi-C, CpG)
        if key == "sat":
            x = _s(s, "sat_one_per_chromosome")
            (k, n), m, p = _kn(x), _CHANCE.search(x), _p(x)
            chance = float(m[1]) if m else None
            q = SAT_QUESTION
        elif key == "te":
            k, n = _kn(_s(s, "te_verdict"))
            q = f"A dense {clade} cluster"
        elif key == "hic":
            te_hic = _s(s, "te_hic")
            if st in ("pass", "support") and not (_num(s.get("hic_gate_pass"), int) or _num(s.get("hic_gate_support"), int)):
                (k, n), p = _kn(te_hic), _p(te_hic)
                q = f"Contacts peak at the {clade} clusters"
            elif st == "untested":
                q = "Contacts the other centromeres: too few calls to test it against"
            else:
                n = _num(s.get("hic_chromosomes"), int) or n_chrom
                k = _num(s.get("hic_gate_pass"), int) or _num(s.get("hic_gate_support"), int) or 0
                per_chrom = True
                q = ("Contacts the other centromeres, beyond chance but not on most" if st == "support" else
                     "Contacts the other centromeres")
        elif key == "flank":
            x = _s(s, "cis_line")
            (k, n), p = _kn(x), _p(x)
            q = ("Flanks touch most where contacts with other chromosomes peak" if "tested by trans" in x else
                 f"Flanks touch most at the {clade} cluster")
        elif key in ("fold", "align"):
            k, n = _kn(_s(s, f"{key}_on_denovo"))
            p = _num(s.get(f"{key}_p"))
            q = "Arms fold back at Hi-C's own peak" if key == "fold" else \
                "Arms of different chromosomes align at Hi-C's own peak"
        elif key == "meth":
            q = ("Measured in other plants: it does not place this plant's centromeres" if st == "other plants" else
                 "Measured in this plant" + (": no centromere pattern found" if st == "fail" else ""))
        elif key == "cpg":
            k, n, per_chrom = _num(s.get("cpg_gate_pass"), int), n_chrom, True
            q = "CpG dips at the other centromeres"
        else:
            p = _num(s.get(f"{key}_p"))
            q = {"at": "The deepest AT-rich trough falls one per chromosome",
                 "rep": "The densest repeats fall one per chromosome",
                 "fam": "Repeat-family clusters fall one per chromosome, a test that often passes by chance"}[key]
        if k is not None and n:
            q += (f": its test passes on {k} of {n} chromosomes" if per_chrom else f": {k} of {n} chromosomes") + \
                (f" (chance: {chance:g})" if chance is not None else "")
        verdict = {"pass": "trusted", "support": "trusted", "lead": "a lead only",
                   "other plants": "not used (other plants)"}.get(st, "not trusted here")
        out.append(TestRow(key, SIGNALS[key].label(clade), q, k, n, chance, p, verdict, per_chrom))
    return out


# ------------------------------------------------------------------------------------------------ roles
def _place(key: str, rec: dict) -> tuple[int, int] | None:
    """Where the signal alone puts the centromere on this chromosome (0-based half-open), None where it puts none."""
    if key == "sat":                                     # the array the call stands on, or the one an LTR core replaced
        if str(rec.get("seq_evidence", ".")) in SAT_KINDS:
            return _iv(rec.get("seq_start"), rec.get("seq_end"))
        return _iv(rec.get("seq_overruled_start"), rec.get("seq_overruled_end"))
    if key == "hic" and str(rec.get("hic_anch_level", ".")) in NO_PEAK:
        return None                                      # Hi-C scored this chromosome but found no peak
    cols = SIGNALS[key].cols
    p = _iv(rec.get(cols[0]), rec.get(cols[1])) if cols else None
    if p is None and key == "cpg":                       # no dip in the call's arrays: the dip on another array
        p = _iv(rec.get("cpg_alt_start"), rec.get("cpg_alt_end"))
    return p


def _sides(rec: dict) -> tuple[set, set, set]:
    """The call's rule and flags, and the signals its `support` and `against` name: the token `hic` is the Hi-C readout
    that acted (support: the cis extreme for a cis rule, the contact hub for Hi-C alone, an axis pair or an AT pair, else
    the contact with the other centromeres; against: the cis extreme for cis_elsewhere, the hub for hic_denovo_offset,
    else the contact with the other centromeres)."""
    tags = set(_toks(rec.get("rule"))) | set(_toks(rec.get("flags")))
    sup = {TOKEN[t] for t in _toks(rec.get("support")) if t in TOKEN}
    if "hic" in _toks(rec.get("support")):
        sup.add("flank" if tags & CIS_SUPPORT else "hub" if tags & HUB_SUPPORT else "hic")
    ag = {TOKEN[t] for t in _toks(rec.get("against")) if t in TOKEN}
    if "hic" in _toks(rec.get("against")):
        hic = {k for k, f in (("flank", "cis_elsewhere"), ("hub", "hic_denovo_offset")) if f in tags}
        ag |= hic or {"hic"}
    return tags, sup, ag


def _clade(rec: dict, s: dict) -> str:
    for x in (rec.get("te_clade"), s.get("te_marker")):
        if _toks(x):
            return str(x)
    return "LTR"


def clade_of(calls: pd.DataFrame, summary: dict | None) -> str:
    """The LTR family's lineage, for labels: the first te_clade of the table, else the summary's marker, else 'LTR'."""
    first = next((x for x in calls.get("te_clade", pd.Series(dtype=str)).astype(str) if _toks(x)), None)
    return _clade({"te_clade": first}, summary or {})


def role_label(key: str, role: str, summary: dict | None) -> str:
    """A role in words: 'used', 'points elsewhere', 'not used', or 'not used (other plants)' for methylation measured
    in other plants only."""
    if role == "not" and key == "meth" and _status("meth", summary or {}) == "other plants":
        return "not used (other plants)"
    return ROLE_WORDS[role]


def _consulted(tags: set, sup: set, ag: set) -> frozenset:
    """The secondary Hi-C readouts and the repeat-family cluster cenfuse consulted on one chromosome (ACTED): named in
    support or against, or by a rule or flag (moved_with_te+fold names the fold)."""
    moved = {f"moved_with_{x}" for t in tags if t.startswith("moved_with_") for x in t[len("moved_with_"):].split("+")}
    return frozenset(k for k in ACTED if k in sup | ag or ACTED[k] & (tags | moved))


def consulted(rows: list[Row]) -> set:
    """The secondary Hi-C readouts and the repeat-family cluster cenfuse consulted in this genome, on any chromosome:
    only these get a close-up track."""
    return set().union(*(r.consulted for r in rows))


def rows(calls: pd.DataFrame, summary: dict | None, hic_bin: int | None = None) -> list[Row]:
    """One Row per chromosome of `<genome>.centromeres.tsv` (read as strings or not). `summary` is the genome's
    cenfuse_summary.tsv row (None or {}: genome tests unknown); `hic_bin` the Hi-C bin, rule 4's distance."""
    s = summary or {}
    tol = int(hic_bin) if hic_bin else NO_HIC_TOL
    recs = calls.to_dict("records")
    sides = [_sides(rec) for rec in recs]
    out = [Row(str(rec["chrom"]), int(rec["chrom_len"]),
               _iv(rec.get("start"), rec.get("end")) if str(rec.get("relation")) != "none" else None,
               str(rec.get("confidence", ".")), {k: _place(k, rec) for k in SIGNALS}, {},
               _iv(rec.get("seq_alt_start"), rec.get("seq_alt_end")), consulted=_consulted(*x),
               hic_scored=_iv(rec.get("hic_anch_core_start"), rec.get("hic_anch_core_end")) is not None)
           for rec, x in zip(recs, sides)]
    for rec, (tags, sup, ag), r in zip(recs, sides, out):
        for k, p in r.places.items():
            far = k in LEFT_BEHIND and s and _trust(k, s, rec) and p and r.call and _gap(p, r.call) > tol
            r.roles[k] = "used" if k in sup else "elsewhere" if k in ag or far else "not"
        r.verdict = verdict(rec, r, _clade(rec, s))
        r.placed = placed_by(rec, r, tags)
    return out


def present(rows: list[Row], summary: dict | None, has: dict) -> list[str]:
    """The signals that get a close-up track, in TRACK_ORDER (the Hi-C trace first): the satellite arrays always; the
    LTR family with tes/ (has['te']); with Hi-C (has['hic']) the contact with the other centromeres, else the hub;
    CHG methylation with its input (has['meth']); every other signal where it passed its genome test or is used or
    points elsewhere on any chromosome, a secondary Hi-C readout or the repeat-family cluster only where cenfuse
    consulted it in this genome."""
    s = summary or {}
    acted = {k for r in rows for k, v in r.roles.items() if v != "not"}
    genome = consulted(rows)
    keep = {"sat"} | {k for k in ("te", "meth") if has.get(k)}
    if has.get("hic"):
        keep.add("hic" if "hic" in acted or any(r.hic_scored or r.places.get("hic") for r in rows) else "hub")
    keep |= {k for k in ("hub", "flank", "fold", "align", "cpg", "at", "rep", "fam")
             if (k not in ACTED or k in genome) and (k in acted or _status(k, s) in ("pass", "support", "lead"))}
    out = [k for k in TRACK_ORDER if k in keep]
    if "hub" in out and "hic" not in out:                   # the hub is the Hi-C trace
        out.insert(0, out.pop(out.index("hub")))
    return out


# ------------------------------------------------------------------------------------------------ verdict
def placed_by(rec: dict, row: Row, tags: set | None = None) -> frozenset:
    """The signals that placed the call, by the verdict's branches: the satellite array or the LTR family's cluster a
    call kept on the DNA sequence stands on (with Hi-C or the AT-rich trough where it chose between two satellite
    clusters; none for a long tandem repeat, which is no satellite array), the Hi-C readouts used where Hi-C placed or
    moved it, the signal that placed it from the DNA alone. Agreeing signals are used but did not place it. Empty
    without a call."""
    if row.call is None:
        return frozenset()
    tags = _sides(rec)[0] if tags is None else tags
    ev, src = str(rec.get("seq_evidence", ".")), str(rec.get("call_source", "."))
    hic = frozenset(k for k in HIC_KEYS if row.roles.get(k) == "used") or frozenset({"hic"})
    if src == "sequence":
        if "hic_chose_alt" in tags:
            return frozenset({"sat"}) | hic
        if "at_chose_alt" in tags:
            return frozenset({"sat", "at"})
        return frozenset({"te"} if ev == "retrotransposon" else set() if ev == "long_unit_tandem" else {"sat"})
    if src == "composition":
        return frozenset({"rep"} if "rep_only" in tags else {"fam"} if "fam_only" in tags else
                         {"at", "hub"} if "at_denovo_pair" in tags else {"at"})
    return hic


def distance_text(bp: int) -> str:
    """A distance in words: '350 kb' below 1 Mb, else '2.3 Mb'."""
    kb = round(bp / 1e3)
    return f"{max(kb, 1) if bp > 0 else 0} kb" if kb < 1000 else f"{bp / 1e6:.1f} Mb"


def _join(xs: list) -> str:
    xs = list(dict.fromkeys(xs))
    return xs[0] if len(xs) == 1 else ", ".join(xs[:-1]) + " and " + xs[-1]


def _cap(x: str) -> str:
    return x[:1].upper() + x[1:]


def verdict(rec: dict, row: Row, clade: str = "LTR") -> str:
    """How the call was made, in one plain sentence: the head by the branches of interpret.basis (what placed it, and
    the signals that agree), then each signal that points elsewhere, a second satellite cluster not chosen, and a CpG
    dip on another array, each with its distance from the call."""
    if row.call is None:
        return "No centromere called"
    tags, _, _ = _sides(rec)
    ev, src = str(rec.get("seq_evidence", ".")), str(rec.get("call_source", "."))
    used = [k for k in SIGNALS if row.roles.get(k) == "used"]

    def others(skip) -> list:                            # the agreeing signals besides the one that placed the call
        xs = [SIGNALS[k].short.format(clade=clade) for k in used if k not in skip]
        return list(dict.fromkeys(xs + ([f"younger {clade} copies"] if "age" in _toks(rec.get("support")) else [])))

    if src == "sequence":
        if "hic_chose_alt" in tags:
            head = "Second satellite cluster, chosen by Hi-C"
        elif "at_chose_alt" in tags:
            head = "Second satellite cluster, chosen by the AT-rich trough"
        elif ev == "retrotransposon":
            x = others({"te"})
            head = f"{clade} cluster" + (f", backed by {_join(x)}" if x else " alone")
        else:
            kind, x = KIND.get(ev, KIND["satellite"]), others({"sat"})
            head = (f"{_join(['satellite'] + x)} agree" if "sat" in used and x else
                    f"{_join(x)} {'agrees' if len(x) == 1 else 'agree'}, on a {kind}" if x else f"{kind} alone")
    elif src == "composition":
        head = ("Placed by repeat density alone (a lead)" if "rep_only" in tags else
                "Placed by a repeat-family cluster alone (a lead)" if "fam_only" in tags else
                "AT-rich trough and Hi-C hub agree" if "at_denovo_pair" in tags else
                "Placed by the AT-rich trough alone")
    else:                                                # Hi-C placed or moved the call
        x = others(HIC_KEYS)
        two = bool(tags & TWO_HIC)
        head = ("Placed by two Hi-C measures" if two else "Placed by Hi-C") + \
            (f", backed by {_join(x)}" if x else "" if two else " alone")
    tails = []
    for k, sig in SIGNALS.items():
        if k == "cpg" or row.roles.get(k) != "elsewhere":
            continue
        p = row.places.get(k)
        verb = "overruled" if k == "sat" and OVERRULED & set(_toks(rec.get("seq_flags"))) else "disagrees"
        tails.append(f"{sig.where.format(clade=clade)} "
                     f"{distance_text(_gap(p, row.call)) + ' away' if p else 'elsewhere'} {verb}")
    if row.alt and not tags & {"hic_chose_alt", "at_chose_alt"}:
        tails.append(f"second satellite cluster {distance_text(_gap(row.alt, row.call))} away")
    cpg_alt = _iv(rec.get("cpg_alt_start"), rec.get("cpg_alt_end"))
    if "cpg_elsewhere" in tags and cpg_alt:
        tails.append(f"CpG dip on another array {distance_text(_gap(cpg_alt, row.call))} away")
    return _cap(head) + "".join(f"; {t}" for t in tails)
