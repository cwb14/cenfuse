"""fuse._decide, split into ordered rule functions (future.md #63), decides every chromosome exactly as the single cascade of
0.14.0 did (ref_decide_014.py): the same call, relation, confidence, source, support, against, flags, rule, second lines and
methylation item, on random evidence drawn so that the lines agree, sit near each other or point elsewhere often."""
import math
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from ref_decide_014 import decide_014  # noqa: E402

from cenfuse import fuse  # noqa: E402

L = 50_000_000
BIN = 250_000
SITES = [0, 400_000, 5_000_000, 9_000_000, 10_000_000, 10_250_000, 10_500_000, 11_000_000, 12_000_000, 20_000_000,
         30_000_000, 44_000_000, 49_500_000]


def _iv(r, p=0.15):
    if r.random() < p:
        return -1, -1
    s = r.choice(SITES)
    return s, min(L, s + r.choice([50_000, 250_000, 400_000, 1_000_000, 3_000_000]))


def _evidence(r: random.Random) -> fuse.Evidence:
    s0, s1 = _iv(r, 0.25)
    a0, a1 = _iv(r, 0.6) if s0 >= 0 else (-1, -1)
    seq = dict(start=s0, end=s1, evidence="satellite" if s0 >= 0 else "none",
               confidence=r.choice(["high", "medium", "low"]) if s0 >= 0 else "none", flags=".", main_family="f",
               alt_start=a0, alt_end=a1)
    te = None
    if r.random() < 0.4:
        t0, t1 = _iv(r, 0)
        te = dict(clade="CRM", start=t0, end=t1, decisive=r.random() < 0.6, strong=r.random() < 0.7)
    h = None
    if r.random() < 0.75:
        p0, _ = _iv(r, 0)
        c0, c1 = (p0, p0 + BIN) if r.random() < 0.5 else (max(0, p0 - BIN), p0 + 2 * BIN)   # the engine's cores start >= 0
        d0, _ = _iv(r, 0.3)
        h = dict(bin_size=BIN, gate_pass=r.random() < 0.5, gate_support=r.random() < 0.3, anch_strong=r.random() < 0.6,
                 anch_level=r.choice(["high", "medium", "low", "none"]), anch_peak_start=p0, anch_peak_end=p0 + BIN,
                 anch_core_start=c0 if r.random() < 0.9 else -1, anch_core_end=c1,
                 anch_flags=r.choice(["", "", "hub", "telomeric", "rDNA"]), gate_refs=r.choice([None, "high", "low", "mixed"]),
                 denovo_peak_start=d0, denovo_peak_end=d0 + BIN if d0 >= 0 else -1,
                 denovo_core_start=max(0, d0 - BIN) if d0 >= 0 else -1, denovo_core_end=d0 + 2 * BIN if d0 >= 0 else -1,
                 denovo_confidence=r.choice(["high", "medium", "low", "none"]),
                 denovo_flags=r.choice(["", "", "hub", "not_converged", "rDNA", "cycle_unstable"]),
                 gate_genome_untested=r.random() < 0.3, end_pass=r.random() < 0.3, end_higher=r.random() < 0.5,
                 masked=[] if r.random() < 0.6 else [(9_500_000, 9_900_000), (10_600_000, 10_900_000)])
        if r.random() < 0.3:
            k0, k1 = _iv(r, 0)
            h["te_pair"] = dict(core=(k0, k1), status=r.choice(["pass", "support"]))
        if r.random() < 0.3:
            q0, _ = _iv(r, 0)
            h["cis"] = dict(peak=(q0, q0 + BIN), status=r.choice(["pass", "support"]), interior=r.random() < 0.5)
    fold = align = None
    if h is not None and r.random() < 0.35:
        f0, f1 = _iv(r, 0)
        fold = dict(status=r.choice(["pass", "fail"]), core_start=f0, core_end=f1, agrees=r.random() < 0.6, point=f0)
    if h is not None and r.random() < 0.3:
        g0, g1 = _iv(r, 0)
        align = dict(status=r.choice(["pass", "fail"]), core_start=g0, core_end=g1, agrees=r.random() < 0.6, point=g0,
                     fold_status=r.choice(["pass", "fail"]), fold_p=r.choice([0.001, 0.03, 0.5]))
    at = None
    if r.random() < 0.35:
        u0, u1 = _iv(r, 0)
        at = dict(start=u0, end=u1, gc=0.3, status=r.choice(["pass", "fail", "untested"]), validated=r.random() < 0.3,
                  support=r.random() < 0.3, agrees=r.random() < 0.5)
    meth = None
    if r.random() < 0.3:
        lab = lambda: dict(label=r.choice(["", "centromere_like", "dip", "other"]), score=r.random())
        meth = {k: lab() for k in ("sequence_call", "hic_peak", "provisional") if r.random() < 0.8}
    age = dict(p=r.choice([1e-4, 0.05]), premise=r.choice([1e-4, 0.3])) if r.random() < 0.2 else None
    bp = r.choice([0, 1_000, 10_000_000])
    return fuse.Evidence("c", L, seq, te, h, meth, (lambda a, b: bp), fold, align, at, age)


def _same(a: dict, b: dict) -> bool:
    for k in set(a) | set(b):
        x, y = a.get(k), b.get(k)
        if k == "meth_call":
            if x.get("label") != y.get("label"):
                return False
            sx, sy = x.get("score"), y.get("score")
            if not (sx == sy or (isinstance(sx, float) and isinstance(sy, float) and math.isnan(sx) and math.isnan(sy))):
                return False
        elif x != y:
            return False
    return True


def test_split_rules_decide_as_the_single_cascade():
    r = random.Random(14)
    n, fired = 20_000, {}
    for i in range(n):
        ev = _evidence(r)
        new, old = fuse._decide(ev), decide_014(ev)
        assert _same(new, old), (i, ev, new, old)
        fired[new["rule"]] = fired.get(new["rule"], 0) + 1
    assert len(fired) >= 18, fired                          # the draw reaches most rules


def test_the_cascade_order_is_declared():
    names = [f.__name__ for f in fuse.CASCADE]
    assert names[0] == "_axis_pair" and names[-1] == "_kept"
    assert [f.__name__ for f in fuse.PHASES] == ["_te_hic_pair", "_cis", "_at_places"]
