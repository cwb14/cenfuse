"""Provenance: does each input come from the reference plant (the assembled genotype), from another plant, or can't
cenfuse tell (docs/plants.md)?

- On this assembly: methylation rows must sit on reference cytosines of their context (fail fast otherwise).
- Genotype (methylation samples). modkit writes no row for a cytosine whose reads all show another base. At a
  symmetric CHG site (CAG/CTG) a fixed genetic difference removes both strands' rows; read errors, failed calls and
  heterozygous sites remove one strand at a time (ONT reads one strand per read). Assembly errors are shared by every
  plant's reads, so a difference that one sample has and another lacks is a real genotype difference of the first:
  it is not the reference plant.
- Karyotype (every input with coverage). A segment that an input lacks while another input covers it shows that the
  input's plant lacks DNA the reference has. The absence of such segments shows nothing.
"""
from __future__ import annotations

import itertools
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import numpy as np
import pandas as pd
import pysam
from scipy.stats import poisson

from .meth.params import DEPTH_HI
from .params import (GENO_ALPHA, GENO_FLANK, GENO_MIN_READS, GENO_POWER, GENO_WINDOW, KARYO_ABSENT,
                     KARYO_MIN_WINDOWS, KARYO_PRESENT)
from .util import CenfuseError, log

MOTIF_MIN = 0.5         # rows on reference motifs below this share: the file belongs to another assembly (error)
MOTIF_WARN = 0.99       # ... below this: a warning (a few rows off-motif are possible at assembly edits)
END_MARGIN = 0.05       # genotype windows stay out of the outer 5% of each chromosome (subtelomeric repeats)
MAX_TRIES = 12          # windows tried per target position before giving up (gap-rich chromosomes)
GAP_SHARE = 0.1         # a karyotype window with more N than this is not judged
SAT_SHARE = 0.5         # ... nor one mostly satellite array (mapping, not DNA, sets its coverage); heuristic
KEY = 1 << 34           # site key = chromosome index x KEY + position (chromosomes < 17 Gb)
VERDICTS = ("reference", "other", "unknown")


@dataclass
class Genotype:
    sample: str
    rows: int               # methylation rows read in the genotype windows (both strands)
    on_motif: float         # share of those rows on reference CHG motifs
    tested: np.ndarray      # sorted keys of the tested symmetric CHG sites
    fixed: np.ndarray       # sorted keys of tested sites with both strands' rows missing (fixed differences)
    plus_only: int          # tested sites with only the + strand's row missing
    minus_only: int         # ... only the - strand's row missing

    @property
    def noise(self) -> float:
        """Per-site rate of both rows missing by chance: q+ x q-, add-one smoothed so it is never zero."""
        n = len(self.tested)
        return ((self.plus_only + 1) / (n + 2)) * ((self.minus_only + 1) / (n + 2))

    def per_kb(self) -> float:
        """Fixed differences per kb of tested bases (a site tests its C, its A/T and its G)."""
        return 1000 * len(self.fixed) / (3 * len(self.tested)) if len(self.tested) else float("nan")


def symmetric_chg(seq: str, offset: int) -> np.ndarray:
    """Positions (offset-based) of the + strand C of every CAG/CTG in seq (upper case): a CHG cytosine on each
    strand, the - strand one at +2."""
    s = np.frombuffer(seq.encode(), dtype=np.uint8)
    if len(s) < 3:
        return np.zeros(0, np.int64)
    ok = (s[:-2] == ord("C")) & ((s[1:-1] == ord("A")) | (s[1:-1] == ord("T"))) & (s[2:] == ord("G"))
    return np.flatnonzero(ok).astype(np.int64) + offset


def on_motif(seq: str, offset: int, pos: np.ndarray, strand: np.ndarray, context: str) -> np.ndarray:
    """Is each row's cytosine a reference cytosine of its context? seq covers [offset, ...) with 2 bp spare on both
    sides; strand '+', '-' or '.' (combined CpG rows sit on the + strand C)."""
    b = np.frombuffer(seq.encode(), dtype=np.uint8)
    i = pos - offset
    ok = (i >= 2) & (i < len(b) - 2)
    i = np.where(ok, i, 2)
    C, G = ord("C"), ord("G")
    plus = strand != "-"
    at = lambda d: b[i + d]
    if context == "CG":
        hit = np.where(plus, (at(0) == C) & (at(1) == G), (at(0) == G) & (at(-1) == C))
    elif context == "CHG":
        hit = np.where(plus, (at(0) == C) & (at(1) != G) & (at(2) == G), (at(0) == G) & (at(-1) != C) & (at(-2) == C))
    else:                                                  # CHH
        hit = np.where(plus, (at(0) == C) & (at(1) != G) & (at(2) != G), (at(0) == G) & (at(-1) != C) & (at(-2) != C))
    return ok & hit


def genotype_windows(chroms: dict, avoid: dict, fasta) -> list:
    """[(chrom, start, end)]: per chromosome, up to two GENO_WINDOW windows nearest 1/4 and 3/4 of its length, clear
    of the `avoid` intervals ({chrom: [(start, end)]}: sequence calls +- 2 D, satellite clusters), of N and of the
    outer END_MARGIN."""
    out = []
    with pysam.FastaFile(str(fasta)) as fa:
        for c, L in chroms.items():
            lo, hi = int(END_MARGIN * L), int((1 - END_MARGIN) * L)
            starts = list(range(lo, hi - GENO_WINDOW + 1, GENO_WINDOW))
            iv = avoid.get(c, [])
            starts = [s for s in starts if not any(a < s + GENO_WINDOW and b > s for a, b in iv)]
            chosen = []
            for target in (L / 4, 3 * L / 4):
                cand = sorted((s for s in starts if s not in chosen), key=lambda s: abs(s + GENO_WINDOW / 2 - target))
                for s in cand[:MAX_TRIES]:
                    if "N" not in fa.fetch(c, s, s + GENO_WINDOW).upper():
                        chosen.append(s)
                        break
            out += [(c, s, s + GENO_WINDOW) for s in sorted(chosen)]
    return out


def _strand_state(pos: np.ndarray, cov: np.ndarray, q: np.ndarray):
    """(present, spanned) at query positions q on one strand: a row at q; rows within GENO_FLANK on both sides with a
    mean of at least GENO_MIN_READS reads."""
    if not len(pos):
        z = np.zeros(len(q), bool)
        return z, z
    idx = np.searchsorted(pos, q)
    j = np.minimum(idx, len(pos) - 1)
    present = (idx < len(pos)) & (pos[j] == q)
    cs = np.concatenate([[0.0], np.cumsum(cov)])
    lo = np.searchsorted(pos, q - GENO_FLANK)
    hi = np.searchsorted(pos, q + GENO_FLANK, side="right")
    left = idx - lo
    right = hi - idx - present
    n = left + right
    tot = cs[hi] - cs[lo] - np.where(present, cov[j], 0)
    spanned = (left > 0) & (right > 0) & (tot >= GENO_MIN_READS * np.maximum(n, 1))
    return present, spanned


def _rows(tb, chrom: str, s: int, e: int):
    """{strand: (positions, reads)} of 5mC rows in [s, e); reads = N_valid + N_delete + N_fail + N_diff."""
    acc = {"+": ([], []), "-": ([], []), ".": ([], [])}
    bad = 0
    for line in tb.fetch(chrom, s, e):
        p = line.split()
        if len(p) < 17:
            bad += 1
            continue
        if p[3] == "h" or p[5] not in acc:               # 5hmC rows repeat the 5mC row's reads
            continue
        try:
            acc[p[5]][0].append(int(p[1]))
            acc[p[5]][1].append(int(p[9]) + int(p[14]) + int(p[15]) + int(p[16]))
        except ValueError:
            bad += 1
    return {k: (np.array(v[0], np.int64), np.array(v[1], float)) for k, v in acc.items()}, bad


def sample_genotype(args) -> Genotype:
    """Fixed differences of one sample's CHG file in the genotype windows (worker; picklable arguments)."""
    name, path, fasta, windows, index = args
    tested, fixed, rows, hits, po, mo = [], [], 0, 0, 0, 0
    with pysam.FastaFile(str(fasta)) as fa, pysam.TabixFile(str(path)) as tb:
        contigs = set(tb.contigs)
        for c, s, e in windows:
            if c not in contigs:
                continue
            o = max(0, s - 2)
            seq = fa.fetch(c, o, e + 2).upper()
            r, _ = _rows(tb, c, s, e)
            for strand in ("+", "-"):
                p = r[strand][0]
                rows += len(p)
                hits += int(on_motif(seq, o, p, np.full(len(p), strand), "CHG").sum())
            sites = symmetric_chg(seq[s - o:], s)
            sites = sites[sites + 2 < e]
            pp, pc = _strand_state(*r["+"], sites)
            mp, mc = _strand_state(*r["-"], sites + 2)
            ok = pc & mc
            key = index[c] * KEY + sites
            tested.append(key[ok])
            fixed.append(key[ok & ~pp & ~mp])
            po += int((ok & ~pp & mp).sum())
            mo += int((ok & pp & ~mp).sum())
    cat = lambda x: np.unique(np.concatenate(x)) if x else np.zeros(0, np.int64)
    return Genotype(name, rows, hits / rows if rows else float("nan"), cat(tested), cat(fixed), po, mo)


MOTIF_SPAN = 200_000     # bp of each window read for the assembly check (CHH files hold ~10x the rows of CHG)


def motif_share(path, context: str, fasta, windows: list, max_rows: int = 20_000) -> tuple[float, int]:
    """(share of rows on reference motifs of `context`, rows checked) in the first MOTIF_SPAN of the first windows of a
    methylation file."""
    hits = n = 0
    with pysam.FastaFile(str(fasta)) as fa, pysam.TabixFile(str(path)) as tb:
        contigs = set(tb.contigs)
        for c, s, e in windows:
            if c not in contigs or n >= max_rows:
                continue
            e = min(e, s + MOTIF_SPAN)
            o = max(0, s - 2)
            seq = fa.fetch(c, o, e + 2).upper()
            r, _ = _rows(tb, c, s, e)
            for strand, (p, _) in r.items():
                p = p[: max_rows - n]
                hits += int(on_motif(seq, o, p, np.full(len(p), strand), context).sum())
                n += len(p)
    return (hits / n if n else float("nan")), n


def check_assembly(samples, fasta, windows) -> list:
    """Every methylation file must sit on this assembly's cytosines of its context. Raises CenfuseError when fewer
    than MOTIF_MIN of its rows do; returns warnings for files below MOTIF_WARN."""
    warn = []
    for smp in samples:
        for ctx, path in smp.files.items():
            share, n = motif_share(path, "CG" if ctx.startswith("CG") else ctx, fasta, windows)
            if not n:
                continue
            if share < MOTIF_MIN:
                raise CenfuseError(f"methylation {smp.name} {ctx} ({path}): only {share:.0%} of {n} rows sit on "
                                   f"reference {ctx} cytosines; the file belongs to another assembly")
            if share < MOTIF_WARN:
                warn.append(f"methylation {smp.name} {ctx}: {share:.1%} of {n} rows on reference {ctx} cytosines")
    return warn


def genotypes(samples, fasta, chroms: dict, windows: list, threads: int) -> list:
    """Genotype per methylation sample (samples: cenmeth Sample objects with a CHG file)."""
    index = {c: i for i, c in enumerate(chroms)}
    jobs = [(s.name, s.files["CHG"], str(fasta), windows, index) for s in samples]
    if len(jobs) > 1 and threads > 1:
        with ProcessPoolExecutor(min(threads, len(jobs))) as ex:
            return list(ex.map(sample_genotype, jobs))
    return [sample_genotype(j) for j in jobs]


def private_tests(gts: list) -> pd.DataFrame:
    """For every ordered pair (sample, vs): fixed differences of `sample` that `vs` lacks, at sites tested in both,
    against the noise of `sample` (Poisson upper tail)."""
    rows = []
    for a, b in itertools.permutations(gts, 2):
        T = np.intersect1d(a.tested, b.tested, assume_unique=True)
        priv = len(np.setdiff1d(np.intersect1d(a.fixed, T, assume_unique=True), b.fixed, assume_unique=True))
        exp = a.noise * len(T)
        p = float(poisson.sf(priv - 1, exp)) if len(T) and priv else (1.0 if len(T) else float("nan"))
        rows.append(dict(sample=a.sample, vs=b.sample, tested=len(T), private=priv, expected=exp, p=p))
    return pd.DataFrame(rows, columns=["sample", "vs", "tested", "private", "expected", "p"])


def genotype_verdicts(gts: list, tests: pd.DataFrame) -> dict:
    """{sample: dict(verdict, reason, vs, p)} from the private-difference tests (docs/plants.md)."""
    out = {g.sample: dict(verdict="unknown", reason="", vs=".", p=float("nan")) for g in gts}
    if len(gts) < 2:
        for g in gts:
            out[g.sample]["reason"] = ("one methylation sample: assembly errors and another plant cannot be told "
                                       "apart")
        return out
    alpha = GENO_ALPHA / (len(gts) * (len(gts) - 1))
    sig = tests[tests.p <= alpha]
    rates = (sig.private - sig.expected) / sig.tested
    r_alt = float(rates.min()) if len(sig) else float("nan")
    for g in gts:
        mine = tests[tests["sample"] == g.sample]
        hit = mine[mine.p <= alpha]
        if len(hit):
            b = hit.loc[hit.p.idxmin()]
            out[g.sample] = dict(verdict="other", vs=b.vs, p=float(b.p),
                                 reason=f"{int(b.private)} fixed differences that {b.vs} lacks, at "
                                        f"{int(b.tested):,} shared sites ({b.expected:.1f} expected by chance)")
            continue
        against = tests[(tests.vs == g.sample) & (tests.p <= alpha)]
        if not len(against) or not np.isfinite(r_alt):
            out[g.sample]["reason"] = "no private differences, and no other sample differs from it"
            continue
        best, power = None, 0.0
        for k in against["sample"]:
            t = mine[mine.vs == k].iloc[0]
            crit = poisson.isf(alpha, t.expected) + 1
            pw = float(poisson.sf(crit - 1, t.expected + r_alt * t.tested))
            if pw > power:
                best, power = k, pw
        if power >= GENO_POWER:
            out[g.sample] = dict(verdict="reference", vs=best, p=float("nan"),
                                 reason=f"no fixed differences beyond chance against any sample (power {power:.2f} "
                                        f"against {best}, which differs from it)")
        else:
            out[g.sample]["reason"] = f"no private differences, but too few shared sites to tell (power {power:.2f})"
    return out


def window_profile(bin_cov: np.ndarray, bin_size: int, w: int, L: int) -> np.ndarray:
    """Coverage summed into windows of w bp (bins of bin_size bp, bin i starting at i x bin_size)."""
    nw = -(-L // w)
    idx = np.minimum((np.arange(len(bin_cov)) * bin_size) // w, nw - 1)
    return np.bincount(idx, weights=np.nan_to_num(bin_cov), minlength=nw)


def satellite_share(arrays: pd.DataFrame, chrom: str, L: int, w: int) -> np.ndarray:
    """Share of each w-bp window of chrom covered by satellite arrays (union)."""
    nw = -(-L // w)
    cov = np.zeros(nw)
    a = arrays[arrays.chrom == chrom].sort_values("start")
    cs, ce = -1, -1
    spans = []
    for s, e in zip(a.start.to_numpy(), a.end.to_numpy()):
        if s > ce:
            if ce > cs:
                spans.append((cs, ce))
            cs, ce = int(s), int(e)
        else:
            ce = max(ce, int(e))
    if ce > cs:
        spans.append((cs, ce))
    for s, e in spans:
        for i in range(max(0, s // w), min(nw, -(-e // w))):
            cov[i] += max(0, min(e, (i + 1) * w) - max(s, i * w))
    return cov / w


def skip_windows(chroms: dict, w: int, avoid: dict, fasta, arrays: pd.DataFrame | None = None) -> dict:
    """{chrom: bool per window}: windows not judged for missing segments: an `avoid` interval (the sequence calls),
    more than GAP_SHARE N, more than SAT_SHARE satellite array, or a chromosome end."""
    out = {}
    with pysam.FastaFile(str(fasta)) as fa:
        for c, L in chroms.items():
            nw = -(-L // w)
            sk = np.zeros(nw, bool)
            sk[0] = sk[-1] = True
            for a, b in avoid.get(c, []):
                sk[max(0, a // w):min(nw, -(-b // w))] = True
            if arrays is not None and len(arrays):
                sk |= satellite_share(arrays, c, L, w) > SAT_SHARE
            for i in range(nw):
                if not sk[i]:
                    seq = fa.fetch(c, i * w, min(L, (i + 1) * w)).upper()
                    sk[i] = seq.count("N") > GAP_SHARE * max(len(seq), 1)
            out[c] = sk
    return out


def medians(profiles: dict, skip: dict) -> dict:
    """{input: median coverage over its judged windows}."""
    med = {}
    for n, prof in profiles.items():
        v = [prof[c][~skip[c]] for c in prof if c in skip]
        v = np.concatenate(v) if v else np.zeros(0)
        med[n] = float(np.median(v)) if len(v) else float("nan")
    return med


def missing_segments(profiles: dict, skip: dict, w: int) -> dict:
    """{input: [dict(chrom, start, end, windows, covered_by)]}: runs of at least KARYO_MIN_WINDOWS windows where the
    input has < KARYO_ABSENT of its median while another input of the same data type has >= KARYO_PRESENT of its own
    (read length and mapping filters differ between types, so only same-type coverage shows that the DNA is
    mappable). Windows that cannot be judged neither count nor break a run: skipped windows, windows no same-type
    input covers, and pile-ups where the input has more than DEPTH_HI x its median (collapsed repeats draw reads
    from elsewhere; cenmeth's depth rule). A window where the input has coverage ends the run. profiles: {(kind,
    name): {chrom: coverage per window}}."""
    med = medians(profiles, skip)
    out = {n: [] for n in profiles}
    for n, prof in profiles.items():
        if not med[n] > 0:
            continue
        kind = n[0] if isinstance(n, tuple) else None
        others = [m for m in profiles if m != n and med[m] > 0 and (m[0] if isinstance(m, tuple) else None) == kind]
        for c, y in prof.items():
            if c not in skip:
                continue
            rel = y / med[n]
            absent = rel < KARYO_ABSENT
            pile = rel > DEPTH_HI
            cover = {m: profiles[m][c] / med[m] >= KARYO_PRESENT for m in others if c in profiles[m]}
            covered = np.any(np.vstack(list(cover.values())), axis=0) if cover else np.zeros(len(y), bool)
            run = []
            for i in range(len(y) + 1):
                if i < len(y) and (skip[c][i] or pile[i] or (absent[i] and not covered[i])):
                    continue
                if i < len(y) and absent[i]:
                    run.append(i)
                    continue
                if len(run) >= KARYO_MIN_WINDOWS:
                    by = sorted(":".join(m) if isinstance(m, tuple) else str(m) for m, cv in cover.items()
                                if cv[run].mean() >= 0.5)
                    out[n].append(dict(chrom=c, start=run[0] * w, end=(run[-1] + 1) * w, windows=len(run),
                                       covered_by=",".join(by)))
                run = []
    return out


def call_absent(profiles: dict, skip: dict, segments: dict, calls: dict, w: int) -> dict:
    """{(plant, chrom): True} where the integrated call lies inside a missing segment of one of the plant's inputs
    (the call's own windows cannot be judged, so the segment's judged windows flank it) and that input has
    < KARYO_ABSENT of its median over every window of the call."""
    med = medians(profiles, skip)
    out = {}
    for (kind, n), segs in segments.items():
        for c, call in calls.items():
            prof = profiles[(kind, n)].get(c)
            if call is None or prof is None or not med[(kind, n)] > 0:
                continue
            y = prof[call[0] // w:(call[1] - 1) // w + 1] / med[(kind, n)]
            if len(y) and (y < KARYO_ABSENT).all() and any(
                    s["chrom"] == c and s["start"] <= call[0] and call[1] <= s["end"] for s in segs):
                out[(n, c)] = True
    return out


def table(verdicts: dict, geno: dict, segments: dict, deciders: list, dec_lib, read_ends: dict) -> pd.DataFrame:
    """One row per input: its provenance verdict and the evidence behind it."""
    gts = {g.sample: g for g in geno.get("genotypes", [])}
    tests = geno.get("tests")
    rows = []
    for (kind, n), v in verdicts.items():
        g = gts.get(n) if kind == "meth" else None
        best = {}
        if g is not None and tests is not None and len(tests):
            mine = tests[tests["sample"] == n]
            if len(mine) and mine.p.notna().any():
                b = mine.loc[mine.p.idxmin()]
                best = dict(private_vs=b.vs, private=int(b.private), private_expected=round(float(b.expected), 2),
                            private_p=float(b.p), shared_sites=int(b.tested))
        seg = segments.get((kind, n), [])
        rows.append(dict(kind=kind, input=n, plant=n, verdict=v["verdict"], reason=v["reason"],
                         decides=(n in deciders) if kind == "meth" else (n == dec_lib),
                         on_motif=round(g.on_motif, 4) if g is not None and np.isfinite(g.on_motif) else np.nan,
                         sites_tested=len(g.tested) if g is not None else np.nan,
                         fixed_per_kb=round(g.per_kb(), 4) if g is not None and len(g.tested) else np.nan,
                         noise_per_site=g.noise if g is not None else np.nan,
                         private_vs=best.get("private_vs", "."), private=best.get("private", np.nan),
                         private_expected=best.get("private_expected", np.nan),
                         private_p=best.get("private_p", np.nan), shared_sites=best.get("shared_sites", np.nan),
                         missing_segments=";".join(f"{s['chrom']}:{s['start'] + 1}-{s['end']}" for s in seg) or ".",
                         read_ends=read_ends.get(n, np.nan) if kind == "hic" else np.nan))
    cols = ["kind", "input", "plant", "verdict", "reason", "decides", "on_motif", "sites_tested", "fixed_per_kb",
            "noise_per_site", "private_vs", "private", "private_expected", "private_p", "shared_sites",
            "missing_segments", "read_ends"]
    return pd.DataFrame(rows, columns=cols)


def combine(meth_names: list, hic_names: list, geno: dict, segments: dict) -> dict:
    """{(kind, name): dict(verdict, reason)} for every input. Missing segments or genotype evidence make an input
    'other'; genotype evidence alone makes it 'reference'; otherwise it takes a same-named input's verdict."""
    out = {}
    for kind, names in (("meth", meth_names), ("hic", hic_names)):
        for n in names:
            seg = segments.get((kind, n), [])
            g = geno.get(n) if kind == "meth" else None
            if seg:
                txt = "; ".join(f"{s['chrom']}:{s['start'] / 1e6:.0f}-{s['end'] / 1e6:.0f} Mb" for s in seg[:4])
                more = f" (+{len(seg) - 4})" if len(seg) > 4 else ""
                out[(kind, n)] = dict(verdict="other", reason=f"lacks DNA another input covers: {txt}{more}")
                if g and g["verdict"] == "reference":
                    out[(kind, n)]["reason"] += "; genotype matched the reference: conflict"
            elif g and g["verdict"] != "unknown":
                out[(kind, n)] = dict(verdict=g["verdict"], reason=g["reason"])
            else:
                out[(kind, n)] = dict(verdict="unknown",
                                      reason=(g["reason"] if g else "no genotype test for this input type"))
    for (kind, n), v in list(out.items()):
        if v["verdict"] != "unknown":
            continue
        twin = [(k2, n2) for (k2, n2), v2 in out.items() if n2 == n and k2 != kind and v2["verdict"] != "unknown"]
        if twin:
            k2 = twin[0][0]
            out[(kind, n)] = dict(verdict=out[twin[0]]["verdict"],
                                  reason=f"by name: the {'methylation sample' if k2 == 'meth' else 'Hi-C library'} "
                                         f"{n} is {out[twin[0]]['verdict']}")
    return out


def log_summary(genome: str, verdicts: dict) -> None:
    n = {v: sorted(f"{k}:{name}" for (k, name), d in verdicts.items() if d["verdict"] == v) for v in VERDICTS}
    log.info(f"{genome}: provenance: " + "; ".join(f"{v} {len(x)}" + (f" ({', '.join(x)})" if x else "")
                                                  for v, x in n.items()))
