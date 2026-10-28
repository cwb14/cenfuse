"""The CpG line: the kinetochore inside each centromere call, from the centromere dip region (CDR) of CpG methylation.

In primates the CENP-A domain sits in a hypomethylated dip of the otherwise methylated active satellite array (the
CDR: Altemose 2022, Gershman 2022, Logsdon 2024). The satellite and Hi-C lines locate the array (0.8-5.6 Mb in T2T
mammals); the dip locates the kinetochore inside it (tens to hundreds of kb).

- **Profile**: mean CpG methylation per BIN (CDR-Finder's 5-kb window) on the bins that centromeric-role arrays
  cover for at least half, inside the final call. Bins outside the arrays are left out: the unmethylated DNA between
  arrays (CpG islands, genes) would be the deepest valleys.
- **Dip**: the most prominent valley of the profile (methylation higher on both sides; topographic prominence, so an
  array edge is never a dip), extended to its half-depth run inside its run of array bins: the bins below the midpoint
  between its bottom and the array median.
- **Dip test** (per call): the dip is longer than in the call's own bins shuffled. N_PERM shuffles, P = (1 + #shuffles
  at least as long) / (1 + N_PERM). Noise makes narrow valleys; a coherent hypomethylated run is not a chance
  arrangement of the array's values.
- **Positional test** (per call): the dip is longer than when the chromosome's CpG values are rotated by N_ROT random
  offsets (methylation from elsewhere on the chromosome, put on the call's own CpGs and array bins). Same size, same
  bins, same CpG density: only the place differs. A dip that methylation from anywhere would make is not a centromere
  mark.
- **Gate** (per chromosome, leave-one-chromosome-out as the Hi-C gate): references are the other chromosomes' high
  calls (high + medium when fewer than GATE_MIN_REF + 1 are high) that can be tested. k of n hold a positional dip (P
  <= DIP_P); the binomial tail at DIP_P (each reference's chance under the null) must be <= GATE_P. No majority is
  required: the positional test is calibrated and conservative (passing primates: 28-95% of tested calls; rotated
  tracks 4-9%, docs/benchmark.md).
- **Role**: where the gate passes, the kinetochore is the dip of every call whose dip test P <= DIP_P. The line never
  moves a call or changes its confidence, and it cannot tell a wrong array (an inactive array has dips too). Tested and
  not adopted: choosing among satellite clusters by their dips (docs/design.md).

Inputs: the CG files of the deciding methylation samples (provenance, docs/plants.md): modkit bedMethyl (5mC + 5hmC
counted methylated), ENCODE bedMethyl (coverage and percent methylated: WGBS pipelines, Megalodon) or a 4-column
bedGraph of the methylated fraction (0-1 or 0-100: public tracks without coverage), tabix-indexed. A bin's methylation pools its sites' calls (methylated over valid calls, every site with a call), the
binomial estimate at any coverage; a bedGraph site counts as one call of its fraction. Several samples are pooled site
by site. Malformed rows are skipped and counted."""
from __future__ import annotations

import io

import numpy as np
import pandas as pd
import pysam
from scipy.signal import find_peaks
from scipy.stats import binom

from .meth.io import ENCODE_COLS, MODKIT_COLS
from .params import GATE_MIN_REF, GATE_P
from .util import log

BIN = 5_000               # CDR-Finder's window (Mastrorosa et al. 2024): >= 4 bins across the smallest CDRs (~20 kb)
MIN_SITES = 10            # CpG rows per bin with data (a 5-kb alpha-satellite bin holds ~30-60 CpGs; heuristic)
ARRAY_COVER = 0.5         # a bin is array sequence when arrays cover most of it (majority rule)
MIN_BINS = 10             # bins with data for a call to be tested (50 kb: a dip needs shoulders; heuristic)
N_PERM = 200              # shuffles and ...
N_ROT = 200               # ... rotations per call: P resolution 0.005, ten times finer than DIP_P
DIP_P = 0.05              # a significant dip (conventional level); each reference's chance under the gate's null
SEED = 7                  # fixed seed: reruns give identical P values
CEN_ROLES = ("centromeric", "cen_associated")


class Track:
    """One CG file: modkit bedMethyl (>= 18 columns), ENCODE bedMethyl (11-17: coverage and percent methylated in
    columns 10-11, modkit's own first 11) or a bedGraph of methylated fraction (4 columns)."""

    def __init__(self, path):
        self.path = str(path)
        self.tb = pysam.TabixFile(self.path)
        self.contigs = set(self.tb.contigs)
        self.kind, self.scale, self.bad, self.read = "bedgraph", 1.0, 0, 0
        first = []
        for c in self.tb.contigs:
            for row in self.tb.fetch(c):
                first.append(row.split())
                if len(first) >= 1000:
                    break
            if len(first) >= 1000:
                break
        self.ncol = max(4, pd.Series([len(p) for p in first]).mode().iat[0]) if first else 4
        if first and self.ncol >= ENCODE_COLS:
            self.kind = "bedmethyl" if self.ncol >= MODKIT_COLS else "encode"
        elif first:
            vals = pd.to_numeric(pd.Series([p[3] if len(p) >= 4 else "" for p in first]), errors="coerce")
            self.scale = 100.0 if vals.max() > 1.0 else 1.0

    def sites(self, chrom: str, s: int | None = None, e: int | None = None):
        """(positions, methylated calls, valid calls) of the CpG rows of chrom (in [s, e) when given), sorted by
        position. modkit bedMethyl: 5mC + 5hmC calls of each 'm' row with a call; ENCODE bedMethyl: coverage x percent
        / 100 of coverage; bedGraph: the fraction as one call."""
        empty = np.zeros(0, np.int64), np.zeros(0), np.zeros(0)
        if chrom not in self.contigs:
            return empty
        lines = list(self.tb.fetch(chrom, s, e)) if s is not None else list(self.tb.fetch(chrom))
        if not lines:
            return empty
        cols = {"bedmethyl": [1, 3, 9, 11, 13], "encode": [1, 9, 10]}.get(self.kind, [1, 3])
        c, parsed = self._columns(lines, cols)
        p = c[1]
        if self.kind == "bedmethyl":                     # 5hmC ('h') rows: their counts are the 'm' row's N_other_mod
            n, m, h = c[9], c[11], c[13]
            ok = np.isfinite(p) & np.isfinite(n) & np.isfinite(m) & np.isfinite(h)
            self.bad += len(lines) - parsed + int((~ok).sum())
            ok &= (c[3] != "h") & (n >= 1)
            pos, meth, cov = p[ok].astype(np.int64), m[ok] + h[ok], n[ok]
        elif self.kind == "encode":
            n, pct = c[9], c[10]
            ok = np.isfinite(p) & np.isfinite(n) & np.isfinite(pct) & (pct >= 0) & (pct <= 100)
            self.bad += len(lines) - parsed + int((~ok).sum())
            ok &= n >= 1
            pos, meth, cov = p[ok].astype(np.int64), n[ok] * pct[ok] / 100, n[ok]
        else:
            v = c[3]
            ok = np.isfinite(p) & np.isfinite(v)
            self.bad += len(lines) - parsed + int((~ok).sum())
            pos, meth = p[ok].astype(np.int64), v[ok] / self.scale
            cov = np.ones(len(pos))
        self.read += len(pos)
        o = np.argsort(pos, kind="stable")
        return pos[o], meth[o], cov[o]

    def _columns(self, lines: list, cols: list):
        """({column: values}, rows parsed): numeric columns as float (NaN where not a number), column 3 as text for
        bedMethyl. pandas' C parser with the file's column count (so a malformed first row cannot reshape the table);
        rows it cannot place are skipped. A plain split per row when the C parser refuses the block."""
        num = lambda x: pd.to_numeric(x, errors="coerce").to_numpy(float)
        try:
            d = pd.read_csv(io.StringIO("\n".join(lines)), sep=r"\s+", header=None, names=range(self.ncol),
                            usecols=cols, dtype={3: str} if self.kind == "bedmethyl" else None, on_bad_lines="skip")
            out = {k: (d[k].to_numpy(object) if k == 3 and self.kind == "bedmethyl" else num(d[k])) for k in cols}
            return out, len(d)
        except (ValueError, pd.errors.ParserError):
            need = max(cols) + 1
            rows = [x for x in (r.split() for r in lines) if len(x) >= need]
            out = {k: (np.array([x[k] for x in rows], dtype=object) if k == 3 and self.kind == "bedmethyl"
                       else num(pd.Series([x[k] for x in rows], dtype=object))) for k in cols}
            return out, len(rows)

    def close(self):
        self.tb.close()


def track_kind(path) -> str:
    """'bedmethyl' or 'bedgraph'."""
    t = Track(path)
    t.close()
    return t.kind


def bedgraph_on_cpg(path, fasta, windows: list, span: int = 200_000, max_rows: int = 20_000) -> tuple[float, int]:
    """(share of a bedGraph's rows that sit on reference CpGs, rows checked) in the first `span` of each window: the C
    at the row start, or one base later or earlier (public tracks differ by one base), whichever offset fits most rows.
    The bedMethyl check (provenance.check_assembly) reads strands and contexts; a bedGraph has neither."""
    hits, n = np.zeros(3), 0
    with pysam.FastaFile(str(fasta)) as fa, pysam.TabixFile(str(path)) as tb:
        contigs = set(tb.contigs)
        for c, s, e in windows:
            if c not in contigs or n >= max_rows:
                continue
            e = min(e, s + span)
            o = max(0, s - 2)
            seq = np.frombuffer(fa.fetch(c, o, e + 3).upper().encode(), np.uint8)
            pos = []
            for r in tb.fetch(c, s, e):
                try:
                    pos.append(int(r.split("\t", 2)[1]))
                except (IndexError, ValueError):
                    continue
            pos = np.array(pos[: max_rows - n], np.int64) - o
            pos = pos[(pos >= 1) & (pos < len(seq) - 2)]
            for k, d in enumerate((0, 1, -1)):
                hits[k] += int(((seq[pos + d] == ord("C")) & (seq[pos + d + 1] == ord("G"))).sum())
            n += len(pos)
    return (float(hits.max()) / n if n else float("nan")), n


def array_cover(arrays: pd.DataFrame, s: int, e: int) -> np.ndarray:
    """Share of each BIN of [s, e) covered by the given arrays (overlaps counted once)."""
    nb = (e - s + BIN - 1) // BIN
    mask = np.zeros(nb * BIN, bool)
    for a0, a1 in zip(arrays.start, arrays.end):
        a0, a1 = max(int(a0), s), min(int(a1), e)
        if a1 > a0:
            mask[a0 - s:a1 - s] = True
    return mask.reshape(nb, BIN).mean(axis=1)


class Binned:
    """The bins of one call: which site falls in which bin, and which bins count (MIN_SITES rows, array sequence)."""

    def __init__(self, pos: np.ndarray, s: int, cover: np.ndarray):
        self.j = (pos - s) // BIN
        self.nb = len(cover)
        self.cnt = np.bincount(self.j, minlength=self.nb)
        self.ok = (self.cnt >= MIN_SITES) & (cover >= ARRAY_COVER)

    def profile(self, meth: np.ndarray, cov: np.ndarray) -> np.ndarray:
        """Pooled methylated fraction per bin (NaN where the bin does not count)."""
        sm = np.bincount(self.j, weights=meth, minlength=self.nb)
        sn = np.bincount(self.j, weights=cov, minlength=self.nb)
        return np.where(self.ok, sm / np.maximum(sn, 1e-12), np.nan)


def best_valley(y: np.ndarray):
    """(prominence, bottom index, first and last index of its run) of the most prominent local minimum of y, searched
    in runs of bins with data (one missing bin is interpolated, two or more break a run); None without one."""
    idx = np.flatnonzero(np.isfinite(y))
    if len(idx) < 3:
        return None
    best = None
    for run in np.split(idx, np.flatnonzero(np.diff(idx) > 2) + 1):
        if len(run) < 3:
            continue
        a, b = run[0], run[-1]
        seg = pd.Series(y[a:b + 1]).interpolate().to_numpy()
        pk, pr = find_peaks(-seg, prominence=0)
        if len(pk):
            i = int(np.argmax(pr["prominences"]))
            if best is None or pr["prominences"][i] > best[0]:
                best = (float(pr["prominences"][i]), a + int(pk[i]), int(a), int(b))
    return best


def dip(y: np.ndarray):
    """(length in bins, first bin, last bin) of the half-depth run around the most prominent valley, inside the valley's
    run of array bins (it never extends across DNA between arrays, which the shuffles do not see); (0, -1, -1)
    without a valley."""
    v = best_valley(y)
    if v is None:
        return 0, -1, -1
    _, j, a, b = v
    z = pd.Series(y[a:b + 1]).interpolate().to_numpy()
    thr = (np.nanmedian(y) + z[j - a]) / 2
    lo = hi = j - a
    while lo > 0 and z[lo - 1] < thr:
        lo -= 1
    while hi < len(z) - 1 and z[hi + 1] < thr:
        hi += 1
    return hi - lo + 1, a + lo, a + hi


def dip_test(y: np.ndarray, rng) -> dict:
    """The dip of y and its shuffle P (NaN when y has fewer than MIN_BINS bins with data)."""
    fin = y[np.isfinite(y)]
    if len(fin) < MIN_BINS:
        return dict(p=np.nan, run=0, lo=-1, hi=-1, bins=len(fin))
    run, lo, hi = dip(y)
    if run == 0:
        return dict(p=1.0, run=0, lo=-1, hi=-1, bins=len(fin))
    ge = sum(dip(rng.permutation(fin))[0] >= run for _ in range(N_PERM))
    return dict(p=(1 + ge) / (1 + N_PERM), run=run, lo=lo, hi=hi, bins=len(fin))


def positional_p(binned: Binned, meth: np.ndarray, cov: np.ndarray, i0: int, i1: int, run: int, rng) -> float:
    """P that methylation from elsewhere on the chromosome makes a dip at least `run` bins long in this call: the
    chromosome's site calls rotated by N_ROT random offsets onto the call's sites i0..i1 (offsets that bring back the
    call's own values are excluded when the chromosome is long enough)."""
    n, m = len(meth), i1 - i0
    if run == 0 or n < 2:
        return 1.0
    offs = rng.integers(m, n - m + 1, N_ROT) if n > 3 * m else rng.integers(1, n, N_ROT)
    idx = np.arange(i0, i1)
    ge = sum(dip(binned.profile(meth[(idx + o) % n], cov[(idx + o) % n]))[0] >= run for o in offs)
    return (1 + ge) / (1 + N_ROT)


def gate(t: pd.DataFrame) -> pd.DataFrame:
    """Leave-one-chromosome-out gate (module doc) on the per-chromosome table t (chrom, confidence, pos_p)."""
    t = t.copy()
    tested = t.pos_p.notna()
    ref = t[tested & t.confidence.eq("high")]
    if len(ref) < GATE_MIN_REF + 1:
        ref = t[tested & t.confidence.isin(["high", "medium"])]
    sig = ref.pos_p <= DIP_P
    ks, ns, ps, ok = [], [], [], []
    for c in t.chrom:
        r = sig[ref.chrom != c]
        k, n = int(r.sum()), len(r)
        p = float(binom.sf(k - 1, n, DIP_P)) if n else np.nan
        ks.append(k)
        ns.append(n)
        ps.append(p)
        ok.append(bool(n >= GATE_MIN_REF and p <= GATE_P))
    t["gate_k"], t["gate_n"], t["gate_p"], t["gate_pass"] = ks, ns, ps, ok
    return t


def run(files: list, calls: dict, arrays: pd.DataFrame) -> dict:
    """files: CG files of the deciding samples; calls: {chrom: (start, end, confidence)} final calls (0-based
    half-open; start -1 = none); arrays: the sequence line's arrays (chrom, start, end, role). Returns dict(table: one
    row per chromosome, bad_rows)."""
    rng = np.random.default_rng(SEED)
    tracks = [Track(f) for f in files]
    cen = arrays[arrays.role.isin(CEN_ROLES)]
    rows = []
    try:
        for c, (s, e, conf) in calls.items():
            r = dict(chrom=c, confidence=conf, dip_p=np.nan, pos_p=np.nan, dip_bins=0, tested_bins=0, kt_start=-1,
                     kt_end=-1, array_meth=np.nan, kt_meth=np.nan)
            rows.append(r)
            if s < 0:
                continue
            got = [t.sites(c) for t in tracks]
            pos = np.concatenate([g[0] for g in got])
            meth = np.concatenate([g[1] for g in got])
            cov = np.concatenate([g[2] for g in got])
            o = np.argsort(pos, kind="stable")
            pos, meth, cov = pos[o], meth[o], cov[o]
            i0, i1 = np.searchsorted(pos, s), np.searchsorted(pos, e)
            b = Binned(pos[i0:i1], s, array_cover(cen[(cen.chrom == c) & (cen.end > s) & (cen.start < e)], s, e))
            y = b.profile(meth[i0:i1], cov[i0:i1])
            d = dip_test(y, rng)
            r.update(dip_p=d["p"], dip_bins=d["run"], tested_bins=d["bins"])
            if np.isfinite(d["p"]):
                r["pos_p"] = positional_p(b, meth, cov, i0, i1, d["run"], rng)
            if d["run"]:
                r.update(kt_start=s + d["lo"] * BIN, kt_end=min(e, s + (d["hi"] + 1) * BIN),
                         array_meth=float(np.nanmedian(y)), kt_meth=float(np.nanmean(y[d["lo"]:d["hi"] + 1])))
    finally:
        bad = sum(t.bad for t in tracks)
        unread = [t for t in tracks if t.bad and not t.read]
        for t in tracks:
            t.close()
    for t in unread:                                     # a format problem, not a few malformed rows
        log.warning(f"CpG: {t.path}: no row could be read ({t.ncol} columns); CpG is read from modkit bedMethyl (18 "
                    f"columns), ENCODE bedMethyl (11: coverage, percent methylated) or a 4-column bedGraph")
    if bad:
        log.warning(f"CpG: {bad:,} malformed rows skipped")
    t = gate(pd.DataFrame(rows))
    shown = t.gate_pass & (t.dip_p <= DIP_P)
    t.loc[~shown, ["kt_start", "kt_end"]] = -1
    log.info(f"CpG: dips in {int((t.dip_p <= DIP_P).sum())}/{int(t.dip_p.notna().sum())} tested calls, positional in "
             f"{int((t.pos_p <= DIP_P).sum())}; gate passes on {int(t.gate_pass.sum())}/{len(t)} chromosomes; "
             f"kinetochores on {int(shown.sum())}")
    return dict(table=t, bad_rows=bad)


COLS = ["cpg_gate", "cpg_gate_k", "cpg_gate_n", "cpg_gate_p", "cpg_dip_p", "cpg_pos_p", "cpg_kt_start", "cpg_kt_end",
        "cpg_array_meth", "cpg_kt_meth"]


def columns(t: pd.DataFrame, res: dict | None) -> pd.DataFrame:
    """The calls table with the CpG columns appended ('.' without the CpG line; coordinates 1-based inclusive as the
    other columns). cpg_gate: pass / fail / untested (fewer than GATE_MIN_REF testable references)."""
    t = t.copy()
    if res is None:
        for c in COLS:
            t[c] = "."
        return t
    x = {r.chrom: r for r in res["table"].itertuples(index=False)}
    na = "."

    def one(c, f):
        r = x.get(c)
        return f(r) if r is not None else na
    status = lambda r: "pass" if r.gate_pass else ("untested" if r.gate_n < GATE_MIN_REF else "fail")
    num = lambda v: v if np.isfinite(v) else na
    t["cpg_gate"] = [one(c, status) for c in t.chrom]
    t["cpg_gate_k"] = [one(c, lambda r: r.gate_k) for c in t.chrom]
    t["cpg_gate_n"] = [one(c, lambda r: r.gate_n) for c in t.chrom]
    t["cpg_gate_p"] = [one(c, lambda r: num(r.gate_p)) for c in t.chrom]
    t["cpg_dip_p"] = [one(c, lambda r: num(r.dip_p)) for c in t.chrom]
    t["cpg_pos_p"] = [one(c, lambda r: num(r.pos_p)) for c in t.chrom]
    t["cpg_kt_start"] = [one(c, lambda r: r.kt_start + 1 if r.kt_start >= 0 else na) for c in t.chrom]
    t["cpg_kt_end"] = [one(c, lambda r: r.kt_end if r.kt_start >= 0 else na) for c in t.chrom]
    t["cpg_array_meth"] = [one(c, lambda r: num(r.array_meth)) for c in t.chrom]
    t["cpg_kt_meth"] = [one(c, lambda r: num(r.kt_meth) if r.kt_start >= 0 else na) for c in t.chrom]
    return t
