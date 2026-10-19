"""Command line: `cenhic call` (centromeres from Hi-C pairs) and `cenhic simulate` (mis-mapping null)."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from . import __version__, params, plot, qc
from .balance import balanced, ice
from .bins import Bins, choose_bin_size, coarsen, matrix_from_pairs
from .confidence import bootstrap, leave_one_out, level, null_test
from .io import (InputError, PairStats, check_readable, cool_resolutions, is_cool, read_cool, read_fai,
                       read_intervals)
from .locate import locate, peak_core, rabl_test, start_positions, telomere_score
from .output import write_bed, write_bedgraph, write_qc, write_tsv
from .score import block_means

log = logging.getLogger("cenhic")


def _logging(verbose):
    h = logging.StreamHandler(sys.stderr)
    h.setFormatter(logging.Formatter("[cenhic %(asctime)s] %(message)s", "%H:%M:%S"))
    log.handlers[:] = [h]
    log.setLevel(logging.DEBUG if verbose else logging.INFO)
    log.propagate = True


def fast_checksum(*paths) -> str:
    """Cache key: a fast fingerprint (version, name, size, mtime, first and last MiB of each input), not a checksum
    of the whole file."""
    n = params.CACHE_FP_BYTES
    h = hashlib.md5(f"{__version__}:{params.BIN_LADDER}:{params.MAPQ}".encode())
    for p in paths:
        p = Path(p)
        st = p.stat()
        h.update(f"{p.name}:{st.st_size}:{st.st_mtime_ns}".encode())
        with open(p, "rb") as f:
            h.update(f.read(n))
            if st.st_size > 2 * n:
                f.seek(-n, 2)
                h.update(f.read())
    return h.hexdigest()[:16]


def _pair_warnings(stats):
    """Pair-level warnings; also re-emitted from the cached stats so a rerun warns like the first run."""
    if stats.get("source") != "pairs":
        return
    if not stats.get("mapq_checked"):
        log.warning("pairs carry no mapq1/mapq2 columns: uniqueness (MAPQ >= 30) was not checked")
    if stats.get("out_of_range"):
        log.warning(f"{stats['out_of_range']} pair ends lie beyond their sequence's .fai length (other assembly?); dropped")
    if stats.get("malformed"):
        log.warning(f"{stats['malformed']} malformed pair lines skipped")


def _warn_short(names_lengths, why):
    """One warning for sequences dropped as too short: count, total Mb, up to 3 names."""
    if names_lengths:
        ex = ", ".join(n for n, _ in names_lengths[:3])
        log.warning(f"{len(names_lengths)} sequences {why} are not called ({sum(L for _, L in names_lengths) / 1e6:.1f} Mb"
                    f" in total; e.g. {ex})")


def _base_matrix(args, seqs):
    """Base (finest) matrix and its Bins, from cache when the inputs are unchanged."""
    cool = is_cool(args.pairs)
    if cool:      # only resolutions inside the ladder's range: finer ones would not fit in memory genome-wide
        sizes = [r for r in cool_resolutions(args.pairs) if params.BIN_LADDER[0] <= r <= params.BIN_LADDER[-1]]
        if not sizes:
            raise InputError(f"{args.pairs}: no resolution between {params.BIN_LADDER[0]} and {params.BIN_LADDER[-1]} bp")
    else:
        sizes = list(params.BIN_LADDER)
    cache = Path(str(args.o) + ".cache")
    key = fast_checksum(args.pairs, args.fai)
    npz, meta = cache / f"{key}.npz", cache / f"{key}.json"
    base = Bins(seqs, sizes[0])
    if npz.exists() and meta.exists():
        try:
            M, stats = sp.load_npz(npz).tocsr(), json.loads(meta.read_text())
        except Exception as e:      # any unreadable cache (truncated zip, bad json, ...) is rebuilt, never fatal
            log.warning(f"cache {npz} is unreadable ({type(e).__name__}: {e}); rebuilding it")
        else:
            log.info(f"cache hit: {npz}")
            _pair_warnings(stats)
            return M, base, sizes, stats
    t0 = time.time()
    if cool:
        M, stats = read_cool(args.pairs, sizes[0], seqs), dict(source="cool")
    else:
        st = PairStats()
        M, n_out = matrix_from_pairs(args.pairs, seqs, base, st)
        stats = dict(source="pairs", read=st.read, kept=st.kept, cis=st.cis, dropped_mapq=st.dropped_mapq,
                     dropped_chrom=st.dropped_chrom, malformed=st.malformed, mapq_checked=st.mapq_checked,
                     out_of_range=n_out, unknown_chroms=st.unknown_chroms)
        if st.kept == 0:
            raise InputError(f"no pairs matched the .fai sequences; pairs use e.g. {st.unknown_chroms[:3]}, "
                             f".fai has e.g. {list(seqs)[:3]}")
        _pair_warnings(stats)
    if M.nnz == 0:
        raise InputError(f"no contacts on the .fai sequences (name mismatch?); .fai has e.g. {list(seqs)[:3]}")
    cache.mkdir(parents=True, exist_ok=True)
    sp.save_npz(npz, M)
    meta.write_text(json.dumps(stats))
    log.info(f"contact matrix built in {time.time() - t0:.0f} s")
    return M, base, sizes, stats


def _over_median(S, bins, members):
    """S divided by its chromosome median, per member chromosome; NaN where a chromosome has no finite score or a
    median <= 0 (no warnings: such chromosomes are simply left out)."""
    out = np.full(bins.N, np.nan)
    for k in np.where(members)[0]:
        sl = bins.sl(int(k))
        v = S[sl]
        if np.isfinite(v).any():
            med = np.nanmedian(v)
            if med > 0:
                out[sl] = v / med
    return out


def _boot_interval(peaks):
    """95% interval of the bootstrap peak bins, rounded outwards to whole bins (floor of the 2.5th percentile,
    ceil of the 97.5th), so the interval never excludes a bin the percentile falls inside."""
    lo, hi = np.percentile(peaks, [2.5, 97.5])
    return int(np.floor(lo)), int(np.ceil(hi))


def _annotation(paths, opt, attrs=True):
    """Intervals of every BED/GFF3 file given to `opt`; a file that yields none is warned about (wrong file?)."""
    dfs = []
    for p in paths:
        d = read_intervals(p, attrs=attrs)
        if d.empty:
            log.warning(f"{opt} {p}: no BED/GFF3 intervals read (wrong file?)")
        else:
            dfs.append(d)
    return pd.concat(dfs, ignore_index=True) if dfs else read_intervals(paths[0], attrs=attrs)


def run_call(args) -> dict:
    # every input is checked and the annotations parsed before the long steps, so a bad file fails in seconds
    check_readable([args.pairs, args.fai, *(args.satellites or []), *(args.repeats or [])])
    cool = is_cool(args.pairs)
    fai = read_fai(args.fai)
    sdf = _annotation(args.satellites, "--satellites") if args.satellites else None
    rep_iv = None
    if args.repeats:
        if cool:
            log.warning("--repeats is ignored with .cool/.mcool input: the read-class split needs the pairs")
        else:       # only the merged coordinates are kept (a repeat annotation can have millions of lines)
            rep_iv = qc.merge_intervals(_annotation(args.repeats, "--repeats", attrs=False))
    seqs = {c: L for c, L in fai.items() if L >= params.MIN_BINS * params.BIN_LADDER[0]}
    if len(seqs) < params.MIN_CALLABLE:
        raise InputError(f"fewer than {params.MIN_CALLABLE} sequences >= {params.MIN_BINS * params.BIN_LADDER[0] / 1e6:g} Mb in {args.fai}")
    _warn_short([(c, L) for c, L in fai.items() if c not in seqs],
                f"shorter than {params.MIN_BINS * params.BIN_LADDER[0] / 1e6:g} Mb")
    log.info(f"cenhic {__version__} call: {len(seqs)} sequences considered")
    M0, base, sizes, stats = _base_matrix(args, seqs)

    def get(size):
        b = Bins(seqs, size)
        if size == base.size:
            return M0, b
        if cool:
            return read_cool(args.pairs, size, seqs), b
        return coarsen(M0, base, b), b

    size, bin_table = choose_bin_size(get, sizes)
    if size is None:
        raise InputError(f"fewer than {params.MIN_CALLABLE} callable sequences at every bin size")
    if not bin_table.loc[bin_table.chosen, "meets_rule"].iloc[0]:
        log.warning(f"shallow library: < {params.MIN_EXPECTED_CEN} expected centromere contacts per bin even at "
                    f"{size} bp; running at {size} bp")
    M, bins = get(size)
    if size != base.size:
        del get, M0          # the finer base matrix is not needed once the chosen size is built (memory)
    members = bins.nb >= params.MIN_BINS
    _warn_short([(n, int(L)) for n, L, m in zip(bins.names, bins.lengths, members) if not m],
                f"with fewer than {params.MIN_BINS} bins at {size} bp")
    log.info(f"bin size {size} bp; {int(members.sum())} callable chromosomes")
    bias, good, n_it, conv = ice(M)
    if not conv:
        log.warning(f"ICE did not converge in {n_it} rounds")
    W = balanced(M, bias, good)
    mean = block_means(W, bins, good, members)
    # a chromosome is scored only through balanced contacts with other callable chromosomes
    no_contacts = members & ~(np.nan_to_num(mean, nan=0.0) > 0).any(axis=1)
    n_with = int((members & ~no_contacts).sum())
    if n_with < params.MIN_CALLABLE:
        raise InputError(f"only {n_with} callable chromosomes have inter-chromosomal contacts (need >= "
                         f"{params.MIN_CALLABLE}); do the pairs and the .fai come from the same assembly?")
    if no_contacts.any():
        nm = [bins.names[k] for k in np.where(no_contacts)[0]]
        log.warning(f"{len(nm)} callable chromosomes have no usable inter-chromosomal contacts and are not called "
                    f"(flag no_contacts): {', '.join(nm[:10])}{' ...' if len(nm) > 10 else ''}")
    T = telomere_score(W, bins, good, mean, members)
    rabl = rabl_test(T, bins, members)
    if rabl["p_floor"] > params.RABL_P:
        log.warning(f"Rabl test cannot reach P <= {params.RABL_P} with {rabl['n']} chromosomes "
                    f"(floor {rabl['p_floor']:.3g}); calls are flagged no_rabl")
    elif not rabl["passed"]:
        log.warning("no Rabl polarity detected: the method's premise fails; calls flagged no_rabl")
    loc = locate(W, bins, good, mean, members, start_positions(T, bins, members))
    unstable = {k for k, s in loc.spread.items() if s > params.CYCLE_SPREAD}
    if not loc.converged:
        log.warning(f"located in {loc.rounds} rounds (not converged)")
    else:
        log.info(f"located in {loc.rounds} rounds" +
                 (f" (cycle of {loc.cycle} states; kept the highest-scoring)" if loc.cycle else "") +
                 (f"; estimate moves > {params.CYCLE_SPREAD} bin within the cycle (flag cycle_unstable): "
                  f"{', '.join(bins.names[k] for k in sorted(unstable))}" if unstable else ""))
    pcs = {k: peak_core(loc.S, bins, k) for k in loc.est}
    pcs = {k: v for k, v in pcs.items() if v is not None}
    est = {k: loc.est[k] for k in pcs}
    log.debug("bootstrap")
    boot = bootstrap(M, bins, members, est, threads=args.t)
    lopo = leave_one_out(loc.oe, bins, good, members, est)
    log.debug("null")
    pnull = null_test(W, bins, good, mean, members, est, {k: v["peak"] for k, v in pcs.items()},
                      {k: v["ratio"] for k, v in pcs.items()})
    hub = qc.hub_bins(M, bins, good)
    sat, rdna = {}, {}
    if sdf is not None:
        sat, rdna = qc.satellite_centres(sdf, bins), qc.rdna_intervals(sdf)
    rows = []
    for k, name in enumerate(bins.names):
        r = dict.fromkeys(["peak_start", "peak_end", "core_start", "core_end", "peak_over_median", "boot_lo", "boot_hi",
                           "boot_in_core", "lopo_min", "lopo_max", "lopo_stable", "null_p", "sat_mid",
                           "peak_minus_sat_Mb", "rounds", "converged"], "")
        r.update(chrom=name, length=int(bins.lengths[k]), callable=bool(members[k]), bin_size=size, confidence="none",
                 flags="no_contacts" if no_contacts[k] else "")
        if k in pcs:
            pc = pcs[k]
            inside = float(np.mean((boot[k] >= pc["core_lo"]) & (boot[k] <= pc["core_hi"])))
            stable = bool(((lopo[k] >= pc["core_lo"]) & (lopo[k] <= pc["core_hi"])).all())
            conf = level(inside, stable, pnull[k])
            g = lambda j: int(bins.off[k] + j)
            bk, lk = boot[k][boot[k] >= 0], lopo[k][lopo[k] >= 0]     # -1 = no finite score (kept as outside the core)
            r.update(peak_over_median=round(pc["ratio"], 3), boot_in_core=round(inside, 3), lopo_stable=stable,
                     null_p=pnull[k], confidence=conf, rounds=loc.rounds, converged=loc.converged,
                     flags=",".join(qc.call_flags(k, pc, bins, hub, rdna, rabl["passed"], pnull[k], loc.converged,
                                                  cycle_unstable=k in unstable)))
            if len(bk):
                bl, bh = _boot_interval(bk)
                r.update(boot_lo=int(bins.start[g(bl)]), boot_hi=int(bins.end[g(bh)]))
            if len(lk):
                r.update(lopo_min=int(bins.start[g(lk.min())]), lopo_max=int(bins.end[g(lk.max())]))
            if conf != "none":
                r.update(peak_start=int(bins.start[g(pc["peak"])]), peak_end=int(bins.end[g(pc["peak"])]),
                         core_start=int(bins.start[g(pc["core_lo"])]), core_end=int(bins.end[g(pc["core_hi"])]))
            if k in sat:
                pm = (bins.start[g(pc["peak"])] + bins.end[g(pc["peak"])]) / 2
                r.update(sat_mid=int(sat[k]), peak_minus_sat_Mb=round((pm - sat[k]) / 1e6, 3))
        rows.append(r)
    out = Path(args.o)
    out.parent.mkdir(parents=True, exist_ok=True)
    ratio, tratio = _over_median(loc.S, bins, members), _over_median(T, bins, members)
    write_bed(rows, f"{out}.centromeres.bed", size)
    write_tsv(rows, f"{out}.centromeres.tsv")
    write_bedgraph(ratio, bins, members, f"{out}.score.bedgraph", "cenhic centromere-cluster score")
    write_bedgraph(tratio, bins, members, f"{out}.telomere.bedgraph", "cenhic telomere-pole score")
    sections = {"input": pd.DataFrame([stats]), "bin_size": bin_table,
                "rabl": pd.DataFrame([{k: v for k, v in rabl.items() if k != "per_chrom"}]),
                "rabl_per_chromosome": rabl["per_chrom"]}
    if rep_iv is not None:
        sections["read_class_split"] = qc.read_class_split(args.pairs, seqs, bins, members, pcs, est, rep_iv)
    write_qc(sections, f"{out}.qc.tsv")
    plot.overview(ratio, bins, members, rows, sat, rabl, f"{out}.overview.pdf")
    n = {lv: sum(r["confidence"] == lv for r in rows) for lv in ("high", "medium", "low", "none")}
    log.info(f"done: {n}")
    # [cenfuse local change; see VENDORED.md] the in-memory state, so cenfuse scores anchored partner windows on
    # the same matrix, balancing and bins without rebuilding them
    state = dict(M=M, W=W, bins=bins, good=good, mean=mean, members=members, loc=loc, pcs=pcs, hub=hub, rdna=rdna,
                 sat=sat, size=size, T=T)
    return dict(rows=rows, rabl=rabl, bin_table=bin_table, stats=stats, state=state)


def _at_least_one(x):
    try:
        n = int(x)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not an integer: {x!r}") from None
    if n < 1:
        raise argparse.ArgumentTypeError("must be >= 1")
    return n


def main(argv=None):
    ap = argparse.ArgumentParser(prog="cenhic", description="Functional centromeres from Hi-C centromere clustering.")
    ap.add_argument("--version", action="version", version=f"cenhic {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("call", help="call one centromere per chromosome from mapped Hi-C pairs")
    c.add_argument("--pairs", required=True, help="4DN .pairs[.gz], simple 6-column pairs, or .cool/.mcool")
    c.add_argument("--fai", required=True, help="genome .fai (sequence names and lengths)")
    c.add_argument("-o", required=True, help="output prefix")
    c.add_argument("--repeats", nargs="+", help="repeat annotation GFF3/BED (QC: read-class split; pairs input only)")
    c.add_argument("--satellites", nargs="+", help="cenrep GFF or BED of satellite centromeres (comparison, rDNA flag)")
    c.add_argument("-t", type=_at_least_one, default=4,
                   help="threads for the bootstrap [4]; each holds a few copies of the contact matrix, so memory "
                        "grows with -t (P. annua, 500 kb, with --repeats: ~1.2 GB at -t 1, ~2.4 GB at -t 8)")
    c.add_argument("-v", "--verbose", action="store_true", help="per-step progress")
    s = sub.add_parser("simulate", help="mis-mapping null: simulated cis-only pairs mapped like the library")
    s.add_argument("--fasta", required=True, help="genome FASTA (indexed); its minimap2 index needs ~4.3 GB RAM per Gb")
    s.add_argument("--pairs", required=True, help="the real pairs, not .cool (cis count for scaling)")
    s.add_argument("-o", required=True, help="prefix of a previous `cenhic call` (reads <o>.centromeres.tsv)")
    s.add_argument("--fastq", help="one real FASTQ for the read-length distribution [150 bp]")
    s.add_argument("-n", type=_at_least_one, default=30_000_000, help="simulated pairs, >= 1 [30000000]")
    s.add_argument("-t", type=_at_least_one, default=12, help="minimap2 threads (index and mapping) [12]")
    s.add_argument("-v", "--verbose", action="store_true", help="per-step progress")
    a = ap.parse_args(argv)
    _logging(a.verbose)
    try:
        if a.cmd == "call":
            run_call(a)
        else:
            from .simulate import run_simulate
            run_simulate(a)
    except (InputError, OSError, EOFError, gzip.BadGzipFile) as e:
        print(f"cenhic: error: {e}", file=sys.stderr)
        sys.exit(2)
