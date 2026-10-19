"""Sequence classes. Every cytosine gets one class, first match wins: satellite family > transposon class (LTRquest
clade, innermost element by nesting depth; EDTA superfamily) > gene body > other. Overlaps are resolved once per
chromosome by painting non-overlapping segments, so classifying a cytosine is one binary search."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

OTHER = "other"
TIERS = (("sat:",), ("LTR:", "TE:"), ("gene",))   # priority order of class-name prefixes


def tier(name: str) -> int:
    for i, prefixes in enumerate(TIERS):
        if name.startswith(prefixes):
            return i
    return len(TIERS)


@dataclass
class ClassMap:
    names: list[str]                                         # class index -> name; the last is 'other'
    segs: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]  # chrom -> (starts, ends, class index), sorted

    @property
    def other(self) -> int:
        return len(self.names) - 1

    def classify(self, chrom: str, pos: np.ndarray) -> np.ndarray:
        return classify(self.segs.get(chrom), pos, self.other)


def classify(segs, pos: np.ndarray, other: int) -> np.ndarray:
    """Class index of each 0-based position (segments from paint(); positions outside every segment -> other)."""
    out = np.full(len(pos), other, dtype=np.int16)
    if segs is None or len(segs[0]) == 0 or len(pos) == 0:
        return out
    s, e, lab = segs
    j = np.searchsorted(s, pos, side="right") - 1
    ok = (j >= 0) & (pos < e[np.clip(j, 0, None)])
    out[ok] = lab[j[ok]]
    return out


def paint(starts: np.ndarray, ends: np.ndarray, rank: np.ndarray, labels: np.ndarray):
    """Non-overlapping segments (starts, ends, label); each piece of sequence takes the label of the covering
    interval with the lowest rank (= highest priority). Adjacent pieces with the same label are merged."""
    if len(starts) == 0:
        z = np.zeros(0, np.int64)
        return z, z, np.zeros(0, np.int16)
    bp = np.unique(np.concatenate([starts, ends]))
    lab = np.full(len(bp) - 1, -1, dtype=np.int32)
    lo, hi = np.searchsorted(bp, starts), np.searchsorted(bp, ends)
    for i in np.argsort(-rank, kind="stable"):     # worst first: better ranks overwrite
        lab[lo[i]:hi[i]] = labels[i]
    s, e = bp[:-1], bp[1:]
    keep = lab >= 0
    s, e, lab = s[keep], e[keep], lab[keep]
    if len(lab) == 0:
        z = np.zeros(0, np.int64)
        return z, z, np.zeros(0, np.int16)
    new = np.ones(len(lab), bool)
    new[1:] = (lab[1:] != lab[:-1]) | (s[1:] != e[:-1])
    first = np.flatnonzero(new)
    return s[first], np.maximum.reduceat(e, first), lab[first].astype(np.int16)


def build(features: list[pd.DataFrame], chroms) -> ClassMap:
    """features: DataFrames[chrom, start, end, cls, depth] (io.read_satellites / read_ltrquest / read_edta /
    read_genes). Rank: tier, then deeper nesting first, then shorter interval first."""
    parts = [x[["chrom", "start", "end", "cls", "depth"]] for x in features if len(x)]
    if not parts:
        return ClassMap([OTHER], {})
    f = pd.concat(parts, ignore_index=True)
    f = f[f.chrom.isin(set(chroms))]
    names = sorted(set(f.cls), key=lambda n: (tier(n), n)) + [OTHER]
    if len(names) > np.iinfo(np.int16).max:
        raise ValueError(f"{len(names)} classes exceed the int16 class index")
    idx = {n: i for i, n in enumerate(names)}
    f = f.assign(lab=f.cls.map(idx), t=f.cls.map(tier), length=f.end - f.start)
    f = f.sort_values(["t", "depth", "length"], ascending=[True, False, True], kind="stable")
    f["rank"] = np.arange(len(f))
    segs = {}
    for c, g in f.groupby("chrom", sort=False):
        segs[c] = paint(g.start.to_numpy(np.int64), g.end.to_numpy(np.int64), g["rank"].to_numpy(),
                        g.lab.to_numpy(np.int32))
    return ClassMap(names, segs)
