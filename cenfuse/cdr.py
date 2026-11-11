"""The CpG line: the kinetochore inside each centromere call, from the centromere dip region (CDR) of CpG methylation,
and a contest where a call's arrays hold no CDR while another array of its chromosome does.

In primates the CENP-A domain sits in a hypomethylated dip of the otherwise methylated active satellite array (the
CDR: Altemose 2022, Gershman 2022, Logsdon 2024); stickleback, sparrows, duck and some reptiles and amphibians show it
too, other vertebrates (chicken, bats, many fish) do not (docs/benchmark.md). The satellite and Hi-C lines locate the
array (0.8-5.6 Mb in T2T mammals); the dip locates the kinetochore inside it (tens to hundreds of kb).

- **Profile** (per chromosome): CpG methylation per BIN (CDR-Finder's 5-kb window), each bin pooling the calls of the
  SMOOTH bins around it (methylated over valid); the background is the running median of the profile over +-SURR.
- **Dips**: runs of at least MIN_RUN bins below FRAC x the background (one missing bin bridged), on array bins only
  (bins that arrays cover for at least half: the unmethylated DNA between arrays, CpG islands and genes, would be the
  deepest valleys). A dip's depth is its mean methylation over its mean background; the deepest dip in a call's
  centromeric-role arrays is its kinetochore candidate.
- **Positional test** (per call): the dip is deeper than the deepest dip that the chromosome's methylation, rotated by
  N_ROT random offsets (any offset, so the observed window is one of the rotations: an exact test up to Monte Carlo
  error; excluding offsets near the call made P <= 0.05 three times as common as 5% when calls span a sixth of their
  chromosome), makes on the same array bins. Same bins, same window: only where the methylation comes from differs.
  Depth, not width or deficit: the partially methylated domains of cultured cells (the T2T chimpanzee track) and
  subtelomeric caps make wide shallow dips anywhere; a CDR falls to 10-35% of its array (survey of the CENP-A truth
  genomes, docs/benchmark.md). A call spanning much of its chromosome cannot reach a small P: its test is weak, never
  anti-conservative.
- **Gate** (per chromosome, leave-one-chromosome-out as the Hi-C gate): references are the other chromosomes' high
  calls (high + medium when fewer than GATE_MIN_REF + 1 are high) that can be tested; their positional P values,
  combined by Fisher's method (independent chromosomes, valid P values), must give P <= GATE_P. Every reference counts
  by its P, so many calls with moderate dips (wide ape calls: few can reach P <= 0.05 on their own) pass as well as a
  few with strong ones; gate_k counts the references with P <= DIP_P.
- **Kinetochore**: where the gate passes, the deepest dip of every call with one.
- **Contest**: where the gate passes, a call whose centromeric arrays have data (MIN_BINS) and no dip, while the
  deepest dip on the chromosome's other satellite arrays passes the positional test on its array, has its CDR elsewhere:
  the fusion contests the call (fuse._cpg: one level down, as a contesting Hi-C line) and never moves it (in the
  truth genomes such a dip was the centromere once and elsewhere twice, docs/benchmark.md). The dip is reported.
- **Coverage holes**: arrays with data in fewer than HOLE x the genome's share of bins are reported: tracks piled up
  at MAPQ 60 lose most reads in arrays (they map at MAPQ 20-59).

Inputs: the CG files of the deciding methylation samples (provenance, docs/plants.md): modkit bedMethyl (5mC + 5hmC
counted methylated), ENCODE bedMethyl (coverage and percent methylated: WGBS pipelines, Megalodon), pb-CpG-tools
(coverage and estimated modified count) or a 4-column bedGraph of the methylated fraction (0-1 or 0-100: public tracks
without coverage), tabix-indexed. A bin's methylation pools its sites' calls (methylated over valid calls, every site
with a call), the binomial estimate at any coverage; a bedGraph site counts as one call of its fraction. Several
samples are pooled site by site. Malformed rows are skipped and counted."""
from __future__ import annotations

import io

import numpy as np
import pandas as pd
import pysam
from scipy.stats import chi2

from .meth.io import ENCODE_COLS, MODKIT_COLS, open_tabix
from .params import GATE_MIN_REF, GATE_P
from .util import log

BIN = 5_000               # CDR-Finder's window (Mastrorosa et al. 2024): >= 4 bins across the smallest CDRs (~20 kb)
SMOOTH = 3                # bins pooled into each bin's methylation (15 kb): a lone 5-kb bin of a sparse track is noise
MIN_SITES = 5             # CpG rows in a SMOOTH window for its methylation to count (heuristic; alpha satellite holds ~30-60
                          # per 5 kb, so only bins with almost no data are dropped)
SURR = 500_000            # the background: running median over +-SURR, wider than CDRs (20-200 kb), narrower than arrays
FRAC = 0.5                # a dip: methylation below half its background (CDRs fall to 10-35% of their arrays; survey)
MIN_RUN = 3               # bins in a dip (15 kb): the narrowest CDRs are ~20 kb
ARRAY_COVER = 0.5         # a bin is array sequence when arrays cover most of it (majority rule)
MIN_BINS = 10             # array bins with data for a call to be tested (50 kb: a dip needs shoulders; heuristic)
N_ROT = 200               # rotations per call: P resolution 0.005, ten times finer than DIP_P (every offset when fewer)
DIP_P = 0.05              # a significant dip (conventional level); each reference's chance under the gate's null
BRIDGE = 2                # bins of non-array sequence an array run spans (an alternative dip's array; heuristic)
HOLE = 0.8                # arrays with data in < HOLE x the genome's share of bins: a coverage hole (MAPQ-60 tracks: 0.36-0.68,
                          # MAPQ >= 1: >= 0.97, four truth genomes)
SEED = 7                  # fixed seed: reruns give identical P values
CEN_ROLES = ("centromeric", "cen_associated")
PBCPG_COLS = 9            # pb-CpG-tools: score, type (Total/hap1/hap2), coverage (6), estimated modified count (7)
PBCPG_TYPES = ("Total", "hap1", "hap2")


class Track:
    """One CG file: modkit bedMethyl (>= 18 columns), ENCODE bedMethyl (11-17: coverage and percent methylated in
    columns 10-11, modkit's own first 11), pb-CpG-tools (9: coverage and estimated modified count in columns 6-7) or
    a bedGraph of methylated fraction (4 columns)."""

    def __init__(self, path):
        self.path = str(path)
        self.tb = open_tabix(self.path)
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
        elif first and self.ncol == PBCPG_COLS and np.mean([len(p) > 4 and p[4] in PBCPG_TYPES for p in first]) > 0.5:
            self.kind = "pbcpg"
        elif first:
            vals = pd.to_numeric(pd.Series([p[3] if len(p) >= 4 else "" for p in first]), errors="coerce")
            self.scale = 100.0 if vals.max() > 1.0 else 1.0

    def sites(self, chrom: str, s: int | None = None, e: int | None = None):
        """(positions, methylated calls, valid calls) of the CpG rows of chrom (in [s, e) when given), sorted by
        position. modkit bedMethyl: 5mC + 5hmC calls of each 'm' row with a call; ENCODE bedMethyl: coverage x percent
        / 100 of coverage; pb-CpG-tools: the estimated modified count of the coverage; bedGraph: the fraction as one
        call."""
        empty = np.zeros(0, np.int64), np.zeros(0), np.zeros(0)
        if chrom not in self.contigs:
            return empty
        lines = list(self.tb.fetch(chrom, s, e)) if s is not None else list(self.tb.fetch(chrom))
        if not lines:
            return empty
        cols = {"bedmethyl": [1, 3, 9, 11, 13], "encode": [1, 9, 10], "pbcpg": [1, 5, 6]}.get(self.kind, [1, 3])
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
        elif self.kind == "pbcpg":
            n, m = c[5], c[6]
            ok = np.isfinite(p) & np.isfinite(n) & np.isfinite(m) & (m >= 0) & (m <= n)
            self.bad += len(lines) - parsed + int((~ok).sum())
            ok &= n >= 1
            pos, meth, cov = p[ok].astype(np.int64), m[ok], n[ok]
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
    """'bedmethyl', 'encode', 'pbcpg' or 'bedgraph'."""
    t = Track(path)
    t.close()
    return t.kind


def bedgraph_on_cpg(path, fasta, windows: list, span: int = 200_000, max_rows: int = 20_000) -> tuple[float, int]:
    """(share of a bedGraph's rows that sit on reference CpGs, rows checked) in the first `span` of each window: the C
    at the row start, or one base later or earlier (public tracks differ by one base), whichever offset fits most rows.
    The bedMethyl check (provenance.check_assembly) reads strands and contexts; a bedGraph has neither."""
    hits, n = np.zeros(3), 0
    with pysam.FastaFile(str(fasta)) as fa, open_tabix(path) as tb:
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


# ---- the profile and its dips
def bins(pos: np.ndarray, meth: np.ndarray, cov: np.ndarray, nb: int) -> np.ndarray:
    """(3, nb): methylated calls, valid calls and CpG rows per BIN of a sequence of nb bins (rows beyond it dropped)."""
    keep = (pos >= 0) & (pos < nb * BIN)
    j = pos[keep] // BIN
    return np.vstack([np.bincount(j, weights=meth[keep], minlength=nb), np.bincount(j, weights=cov[keep], minlength=nb),
                      np.bincount(j, minlength=nb).astype(float)])


def profile(b: np.ndarray) -> np.ndarray:
    """Methylated fraction per bin, pooling the calls of the SMOOTH bins around it; NaN where those hold fewer than
    MIN_SITES rows."""
    k, h = np.ones(SMOOTH), SMOOTH // 2
    m, c, n = (np.convolve(x, k)[h:h + len(x)] for x in b)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where((n >= MIN_SITES) & (c > 0), m / c, np.nan)


def background(y: np.ndarray) -> np.ndarray:
    """Running median of y over +-SURR (NaN where fewer than a quarter of the window's half has data)."""
    half = SURR // BIN
    return pd.Series(y).rolling(2 * half + 1, center=True, min_periods=max(5, half // 4)).median().to_numpy()


def dips(y: np.ndarray, bg: np.ndarray, mask: np.ndarray) -> list:
    """[(depth, first bin, end bin)] of the runs of >= MIN_RUN mask bins below FRAC x bg, deepest first; one missing bin
    inside a run is bridged. depth = the run's mean methylation over its mean background."""
    with np.errstate(invalid="ignore"):
        low = (y < FRAC * bg) & mask
    gap = ~np.isfinite(y) & mask & np.r_[False, low[:-1]] & np.r_[low[1:], False]
    low = low | gap
    d = np.diff(np.r_[0, low.astype(np.int8), 0])
    out = []
    for i, j in zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)):
        if j - i >= MIN_RUN:
            out.append((float(np.nanmean(y[i:j]) / np.nanmean(bg[i:j])), int(i), int(j)))
    return sorted(out)


def array_bins(arrays: pd.DataFrame, chrom: str, nb: int) -> np.ndarray:
    """The bins of chrom (nb of them) that the arrays cover for at least ARRAY_COVER (overlaps counted once)."""
    a = arrays[arrays.chrom == chrom]
    cover = np.zeros(nb)
    end = -1
    for s, e in sorted(zip(a.start.astype(int), a.end.astype(int))):
        s, e = max(s, end, 0), min(e, nb * BIN)
        if e <= s:
            continue
        end = e
        i, j = s // BIN, (e - 1) // BIN
        if i == j:
            cover[i] += e - s
        else:
            cover[i] += (i + 1) * BIN - s
            cover[j] += e - j * BIN
            cover[i + 1:j] += BIN
    return cover / BIN >= ARRAY_COVER


def _run(mask: np.ndarray, a: int, b: int) -> tuple:
    """The run of mask bins around [a, b), spanning up to BRIDGE bins outside the mask."""
    lo, hi = a, b
    while lo > 0 and mask[max(0, lo - 1 - BRIDGE):lo].any():
        lo -= 1
    while hi < len(mask) and mask[hi:hi + 1 + BRIDGE].any():
        hi += 1
    while lo < a and not mask[lo]:
        lo += 1
    while hi > b and not mask[hi - 1]:
        hi -= 1
    return lo, hi


class Chrom:
    """One chromosome's CpG profile: methylation y and its background bg per BIN, the centromeric-role and all
    satellite-array bins. Its verdicts on an interval are memoized (the fusion asks for a call more than once)."""

    def __init__(self, b: np.ndarray, cen: np.ndarray, sat: np.ndarray, length: int | None = None):
        self.y = profile(b)
        self.bg = background(self.y)
        self.cen, self.sat, self.nb = cen, sat, len(self.y)
        self.length = length or self.nb * BIN
        self._memo = {}

    def window(self, s: int, e: int) -> tuple:
        return max(0, s // BIN), min(self.nb, -(-e // BIN))

    def deepest(self, i: int, j: int, mask: np.ndarray, o: int = 0):
        """The deepest dip of bins [i, j) on mask, the methylation moved by o bins (rotation); None without one."""
        idx = (np.arange(i, j) + o) % self.nb if o else slice(i, j)
        d = dips(self.y[idx], self.bg[idx], mask)
        return d[0] if d else None

    def positional_p(self, i: int, j: int, mask: np.ndarray, depth: float) -> float:
        """P that methylation from elsewhere on the chromosome makes a dip at least `depth` deep on the bins [i, j)
        of mask: N_ROT rotations by offsets drawn from all others (every one when there are fewer)."""
        if self.nb < 2:
            return 1.0
        offs = np.arange(1, self.nb) if self.nb - 1 <= N_ROT else \
            np.random.default_rng([SEED, i, j]).integers(1, self.nb, N_ROT)
        ge = 0
        for o in offs:
            d = self.deepest(i, j, mask, int(o))
            ge += d is not None and d[0] <= depth
        return (1 + ge) / (1 + len(offs))

    def assess(self, s: int, e: int) -> dict:
        """The call [s, e) on its centromeric arrays: bins with data, the share of array bins with data, the deepest
        dip (depth 1 without one) and its positional P, the kinetochore candidate and the methylation of both."""
        key = ("assess", s, e)
        if key in self._memo:
            return self._memo[key]
        i, j = self.window(s, e)
        m = self.cen[i:j]
        y = self.y[i:j]
        has = np.isfinite(y) & m
        r = dict(tested_bins=int(has.sum()), array_bins=int(m.sum()), array_data=float(has.sum() / m.sum()) if m.any()
                 else np.nan, dip_depth=np.nan, pos_p=np.nan, kt_start=-1, kt_end=-1,
                 array_meth=float(np.median(y[has])) if has.any() else np.nan, kt_meth=np.nan)
        if r["tested_bins"] >= MIN_BINS:
            d = self.deepest(i, j, m)
            if d is None:
                r.update(dip_depth=1.0, pos_p=1.0)
            else:
                depth, a, b = d
                r.update(dip_depth=depth, pos_p=self.positional_p(i, j, m, depth), kt_start=max(s, (i + a) * BIN),
                         kt_end=min(e, (i + b) * BIN), kt_meth=float(np.nanmean(y[a:b])))
        self._memo[key] = r
        return r

    def elsewhere(self, s: int, e: int):
        """(start, end) of the CDR-like dip on another satellite array of the chromosome when the call [s, e) has data
        and no dip on its centromeric arrays (the contest, module doc); None otherwise. The caller checks the gate."""
        key = ("elsewhere", s, e)
        if key in self._memo:
            return self._memo[key]
        out = None
        r = self.assess(s, e)
        if r["tested_bins"] >= MIN_BINS and r["kt_start"] < 0:
            i, j = self.window(s, e)
            other = self.sat.copy()
            other[i:j] = False
            d = dips(self.y, self.bg, other)
            if d:
                depth, a, b = d[0]
                lo, hi = _run(other, a, b)
                if self.positional_p(lo, hi, other[lo:hi], depth) <= DIP_P:
                    out = (a * BIN, min(b * BIN, self.length))
        self._memo[key] = out
        return out


def prepare(files: list, arrays: pd.DataFrame, lengths: dict | None = None, chroms=None) -> dict:
    """files: CG files of the deciding samples; arrays: the sequence line's arrays (chrom, start, end, role); lengths:
    {chrom: length} (else the last row, call or array end). Returns dict(chroms: {chrom: Chrom}, bad_rows, data: the
    genome's share of bins with data)."""
    tracks = [Track(f) for f in files]
    cen = arrays[arrays.role.isin(CEN_ROLES)]
    names = list(lengths) if lengths else list(chroms or [])
    out, with_data, total = {}, 0, 0
    try:
        for c in names:
            got = [t.sites(c) for t in tracks]
            if not any(len(g[0]) for g in got):
                continue
            n = lengths[c] if lengths else max([int(g[0][-1]) + 1 for g in got if len(g[0])] +
                                               [int(x) for x in arrays.end[arrays.chrom == c]])
            nb = -(-n // BIN)
            ch = Chrom(sum(bins(*g, nb) for g in got), array_bins(cen, c, nb), array_bins(arrays, c, nb), n)
            out[c] = ch
            with_data += int(np.isfinite(ch.y).sum())
            total += nb
    finally:
        bad = sum(t.bad for t in tracks)
        unread = [t for t in tracks if t.bad and not t.read]
        for t in tracks:
            t.close()
    for t in unread:                                     # a format problem, not a few malformed rows
        log.warning(f"CpG: {t.path}: no row could be read ({t.ncol} columns); CpG is read from modkit bedMethyl (18 "
                    f"columns), ENCODE bedMethyl (11: coverage, percent methylated), pb-CpG-tools (9) or a 4-column "
                    f"bedGraph")
    if bad:
        log.warning(f"CpG: {bad:,} malformed rows skipped")
    return dict(chroms=out, bad_rows=bad, data=with_data / total if total else np.nan)


def gate(t: pd.DataFrame) -> pd.DataFrame:
    """Leave-one-chromosome-out gate (module doc) on the per-chromosome table t (chrom, confidence, pos_p)."""
    t = t.copy()
    tested = t.pos_p.notna()
    ref = t[tested & t.confidence.eq("high")]
    if len(ref) < GATE_MIN_REF + 1:
        ref = t[tested & t.confidence.isin(["high", "medium"])]
    ks, ns, ps, ok = [], [], [], []
    for c in t.chrom:
        r = ref.pos_p[ref.chrom != c]
        k, n = int((r <= DIP_P).sum()), len(r)
        p = float(chi2.sf(-2 * np.log(r.clip(lower=1e-300)).sum(), 2 * n)) if n else np.nan
        ks.append(k)
        ns.append(n)
        ps.append(p)
        ok.append(bool(n >= GATE_MIN_REF and p <= GATE_P))
    t["gate_k"], t["gate_n"], t["gate_p"], t["gate_pass"] = ks, ns, ps, ok
    return t


def assess(prep: dict, calls: dict, gates: pd.DataFrame | None = None, contested=None) -> pd.DataFrame:
    """One row per chromosome of calls {chrom: (start, end, confidence)} (0-based half-open; start -1 = none): Chrom.assess,
    the gate (computed here, or taken from `gates`: the provisional calls' gate, so the final calls are judged by the
    references that decided them), the kinetochore where the gate passes and the contest's dip (alt_start, alt_end):
    for the chromosomes in `contested` (the calls the fusion contested), or for every call where the gate passes."""
    rows = []
    for c, (s, e, conf) in calls.items():
        r = dict(chrom=c, confidence=conf, tested_bins=0, array_bins=0, array_data=np.nan, dip_depth=np.nan,
                 pos_p=np.nan, kt_start=-1, kt_end=-1, array_meth=np.nan, kt_meth=np.nan)
        ch = prep["chroms"].get(c)
        if s >= 0 and ch is not None:
            r.update(ch.assess(s, e))
        rows.append(r)
    t = pd.DataFrame(rows)
    if gates is None:
        t = gate(t)
    else:
        g = gates.set_index("chrom")
        t["gate_k"] = t.chrom.map(g.gate_k.to_dict()).fillna(0).astype(int)
        t["gate_n"] = t.chrom.map(g.gate_n.to_dict()).fillna(0).astype(int)
        t["gate_p"] = t.chrom.map(g.gate_p.to_dict()).astype(float)
        t["gate_pass"] = t.chrom.map(g.gate_pass.to_dict()).eq(True)
    t.loc[~t.gate_pass, ["kt_start", "kt_end"]] = -1
    alt = [prep["chroms"][c].elsewhere(s, e) if ok and s >= 0 and c in prep["chroms"] and
           (contested is None or c in contested) else None for (c, (s, e, _)), ok in zip(calls.items(), t.gate_pass)]
    t["alt_start"] = [a[0] if a else -1 for a in alt]
    t["alt_end"] = [a[1] if a else -1 for a in alt]
    return t


def array_data(t: pd.DataFrame, prep: dict) -> float:
    """The median share of array bins with data over the calls with MIN_BINS array bins, over the genome's share."""
    x = t.array_data[t.array_bins >= MIN_BINS]
    return float(x.median() / prep["data"]) if len(x) and prep["data"] > 0 else np.nan


def summary(t: pd.DataFrame, prep: dict) -> float:
    """Log the line's counts and a coverage hole; returns the array data ratio (array_data)."""
    ad = array_data(t, prep)
    shown = t.kt_start >= 0
    log.info(f"CpG: dips in {int((t.dip_depth < 1).sum())}/{int(t.dip_depth.notna().sum())} tested calls, positional "
             f"in {int((t.pos_p <= DIP_P).sum())}; gate passes on {int(t.gate_pass.sum())}/{len(t)} chromosomes; "
             f"kinetochores on {int(shown.sum())}; another array's dip on {int((t.alt_start >= 0).sum())}")
    if np.isfinite(ad) and ad < HOLE:
        log.warning(f"CpG: the arrays of the calls have data in {ad:.0%} of the genome's share of bins: a coverage hole "
                    f"(reads in arrays map at MAPQ 20-59; tracks piled up at MAPQ 60 lose most of them). Pile up at "
                    f"MAPQ >= 1 (pb-CpG-tools --min-mapq 1, modkit pileup without a MAPQ filter)")
    return ad


def run(files: list, calls: dict, arrays: pd.DataFrame, lengths: dict | None = None) -> dict:
    """prepare + assess for one set of calls (tests, and callers without a provisional gate). Returns dict(table,
    bad_rows, array_data, chroms)."""
    prep = prepare(files, arrays, lengths, chroms=list(calls))
    t = assess(prep, calls)
    return dict(table=t, bad_rows=prep["bad_rows"], array_data=summary(t, prep), chroms=prep["chroms"])


COLS = ["cpg_gate", "cpg_gate_k", "cpg_gate_n", "cpg_gate_p", "cpg_dip_depth", "cpg_pos_p", "cpg_kt_start",
        "cpg_kt_end", "cpg_array_meth", "cpg_kt_meth", "cpg_array_data", "cpg_alt_start", "cpg_alt_end"]


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
    t["cpg_dip_depth"] = [one(c, lambda r: num(r.dip_depth)) for c in t.chrom]
    t["cpg_pos_p"] = [one(c, lambda r: num(r.pos_p)) for c in t.chrom]
    t["cpg_kt_start"] = [one(c, lambda r: r.kt_start + 1 if r.kt_start >= 0 else na) for c in t.chrom]
    t["cpg_kt_end"] = [one(c, lambda r: r.kt_end if r.kt_start >= 0 else na) for c in t.chrom]
    t["cpg_array_meth"] = [one(c, lambda r: num(r.array_meth)) for c in t.chrom]
    t["cpg_kt_meth"] = [one(c, lambda r: num(r.kt_meth) if r.kt_start >= 0 else na) for c in t.chrom]
    t["cpg_array_data"] = [one(c, lambda r: num(r.array_data)) for c in t.chrom]
    t["cpg_alt_start"] = [one(c, lambda r: r.alt_start + 1 if r.alt_start >= 0 else na) for c in t.chrom]
    t["cpg_alt_end"] = [one(c, lambda r: r.alt_end if r.alt_start >= 0 else na) for c in t.chrom]
    return t
