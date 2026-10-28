"""Chromosome morphology and assembly completeness from the integrated call (docs/design.md, 'Morphology').

The integrated call fixes where a chromosome's arms meet. Each called chromosome gets its centromeric index (short arm
/ length), arm ratio (long / short) and Levan class, measured from the CpG kinetochore's midpoint where the CpG line
placed one (the finest position of the centromere), else from the call's midpoint. Classes (Levan, Fredga & Sandberg
1964, by arm ratio): metacentric <= 1.7, submetacentric <= 3.0, subtelocentric <= 7.0, acrocentric beyond (telocentric
included: a sequence cannot show a missing short arm). When the call (or kinetochore) spans two classes, both are given.

short_arm: the sequence end nearer the centromere; '.' when the call (or kinetochore) spans the chromosome's middle,
where the arms cannot be ordered (CHM13 chr1: its p arm is 50 kb the longer by sequence).

Completeness, per chromosome: telomere arrays at its ends (cenfuse.telomere), assembly gaps (runs of >= GAP_MIN_N N)
on it and inside the call, and t2t = telomere arrays at both ends and no gap."""
from __future__ import annotations

import re

import numpy as np
import pandas as pd
import pysam

LEVAN = ((1.7, "metacentric"), (3.0, "submetacentric"), (7.0, "subtelocentric"))   # arm ratio (long / short)
                           # upper bounds of the classes, beyond them acrocentric: Levan, Fredga & Sandberg 1964
GAP_MIN_N = 10             # an assembly gap is a run of >= 10 N: the satellite engine's min_gap_n (assemblers write
                           # gaps as runs of >= 10 N; shorter runs are ambiguity codes)
SCAN_CHUNK = 10_000_000    # bp read at once in the gap scan (memory)
SHORT = {"metacentric": "m", "submetacentric": "sm", "subtelocentric": "st", "acrocentric": "t"}
ORDER = list(SHORT)
COLS = ["cen_point", "short_arm", "ci", "arm_ratio", "morphology", "telomeres", "gaps", "call_gaps", "t2t"]


def scan_gaps(fasta, chroms: dict, min_n: int = GAP_MIN_N, chunk: int = SCAN_CHUNK) -> dict:
    """{chrom: [(start, end)]} 0-based half-open runs of >= min_n N (either case), read chunk by chunk; a run that
    crosses a chunk edge is one gap."""
    pat, out = re.compile("[Nn]+"), {}
    with pysam.FastaFile(str(fasta)) as fa:
        for c, L in chroms.items():
            runs = []
            for s in range(0, L, chunk):
                for m in pat.finditer(fa.fetch(c, s, min(L, s + chunk))):
                    a, b = s + m.start(), s + m.end()
                    if runs and runs[-1][1] == a:
                        runs[-1] = (runs[-1][0], b)
                    else:
                        runs.append((a, b))
            out[c] = [(a, b) for a, b in runs if b - a >= min_n]
    return out


def levan(r: float) -> str:
    """Levan class of an arm ratio (long / short)."""
    return next((name for bound, name in LEVAN if r <= bound), "acrocentric")


def _ratio(x: float, L: int) -> float:
    short = min(x, L - x)
    return np.inf if short <= 0 else (L - short) / short


def classes(s: int, e: int, L: int) -> list:
    """Levan classes met over [s, e): from the point nearest the chromosome's middle to the point nearest an end."""
    far = max(_ratio(s, L), _ratio(e, L))
    mid = 1.0 if s <= L / 2 <= e else min(_ratio(s, L), _ratio(e, L))
    return ORDER[ORDER.index(levan(mid)):ORDER.index(levan(far)) + 1]


def _num(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def columns(t: pd.DataFrame, chroms: dict, telo: dict | None, gaps: dict) -> pd.DataFrame:
    """t with COLS appended (calls 1-based inclusive, as in the calls table; '.' without a call)."""
    sides = {}
    for c, side, *_ in (telo or {}).get("ends", []):
        sides.setdefault(c, set()).add(side)
    rows = []
    for r in t.itertuples(index=False):
        L = int(chroms.get(r.chrom, r.chrom_len))
        g = gaps.get(r.chrom, [])
        sd = sides.get(r.chrom, set())
        tel = "both" if len(sd) == 2 else next(iter(sd)) if sd else "none"
        d = dict(cen_point=".", short_arm=".", ci=".", arm_ratio=".", morphology=".", telomeres=tel, gaps=len(g),
                 call_gaps=".", t2t="yes" if tel == "both" and not g else "no")
        s, e = _num(r.start), _num(r.end)
        if s is not None and e is not None:
            s -= 1
            ks, ke = _num(getattr(r, "cpg_kt_start", None)), _num(getattr(r, "cpg_kt_end", None))
            a, b = (ks - 1, ke) if ks is not None and ke is not None else (s, e)
            p = int(round((a + b) / 2))
            ratio = _ratio(p, L)
            side = "." if a <= L / 2 <= b else "start" if p < L - p else "end"   # spans the middle: arms unordered
            d.update(cen_point=p, short_arm=side, ci=round(min(p, L - p) / L, 3),
                     arm_ratio="inf" if np.isinf(ratio) else round(float(ratio), 2),
                     morphology="/".join(classes(a, b, L)), call_gaps=sum(1 for x, y in g if x < e and y > s))
        rows.append(d)
    return pd.concat([t.reset_index(drop=True), pd.DataFrame(rows, columns=COLS)], axis=1)


def summary(t: pd.DataFrame) -> dict:
    """Genome-level fields: the karyotype (Levan class of each high or medium call at its point: counts; low
    calls and chromosomes without a call counted apart: their morphology is a guess), complete chromosomes,
    calls holding a gap, and high or medium calls on chromosomes that are not metacentric with the short arm at the
    sequence end (the cytogenetic convention puts the short arm, p, first; metacentric arms are too alike to tell)."""
    called = t[t.cen_point != "."]
    k = {x: 0 for x in SHORT.values()}
    sure = called[called.confidence.isin(["high", "medium"])]
    for r in sure.itertuples(index=False):
        k[SHORT[levan(float(r.arm_ratio))]] += 1
    parts = [f"{x}:{n}" for x, n in k.items() if n] + \
        [f"{x}:{n}" for x, n in (("low", len(called) - len(sure)), ("no_call", len(t) - len(called))) if n]
    return dict(karyotype=",".join(parts) or ".", t2t_chromosomes=f"{int((t.t2t == 'yes').sum())}/{len(t)}",
                calls_with_gaps=int((called.call_gaps.astype(int) > 0).sum()),
                short_arm_at_end=int(sum(r.short_arm == "end" and levan(float(r.arm_ratio)) != "metacentric"
                                         for r in sure.itertuples(index=False))))
