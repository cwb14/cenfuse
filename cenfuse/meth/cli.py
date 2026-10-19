"""cenmeth run: inputs -> class bins (cached in <output folder>/cenmeth_cache/, shared by runs in that folder) ->
depth filter -> calibration at the anchors (if given) -> candidates -> outputs. Errors in the input stop the run with
'cenmeth: error: ...' (exit 2)."""
from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

import numpy as np

from . import __version__, bins, calibrate, candidates, classes, io, output, params, plot, track
from .scan import scan

log = logging.getLogger("cenmeth")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="cenmeth", description="Centromere-like regions from DNA methylation "
                                "(exploratory evidence, never bona fide centromere calls).")
    p.add_argument("--version", action="version", version=f"cenmeth {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="scan samples for centromere-like methylation; calibrate on anchors if given")
    r.add_argument("--samples", required=True, help="TSV: sample, CHG [, CG, CHH] = modkit bedMethyl (bgzip+tabix)")
    r.add_argument("--fai", required=True, help="genome .fai (same assembly as every other input)")
    te = r.add_mutually_exclusive_group(required=True)
    te.add_argument("--ltrquest", help="LTRquest *_LTRs_all_depth_LTR_cleaned.gff3 (preferred: clades, nesting)")
    te.add_argument("--edta", help="EDTA *.TEanno.gff3 (superfamily classes; coarser)")
    r.add_argument("--satellites", required=True,
                   help="cenrep <genome>_satellite_arrays.gff, or any BED/GFF3 with the family as name")
    r.add_argument("--genes", help="GFF3 with gene features (gene-body class)")
    r.add_argument("--anchors", help="known centromeres: cenhic BED, cenrep GFF3, or BED (e.g. CENH3 ChIP)")
    r.add_argument("--exclude", help="BED of regions to drop; optional 4th column = sample (blank = all)")
    r.add_argument("-o", "--out", required=True, help="output prefix, e.g. out/x")
    r.add_argument("-t", "--threads", type=int, default=8, help="parallel processes [8]")
    r.add_argument("-v", "--verbose", action="store_true", help="per-step progress and checks")
    return p


def _check_names(what: str, names, chroms: list[str], fai: dict) -> None:
    """Fail when no feature lies on an analysed chromosome; warn about features on other sequences."""
    names = set(names)
    if not names & set(chroms):
        raise io.InputError(f"{what}: no feature on an analysed chromosome (e.g. {', '.join(sorted(names)[:3])} vs "
                            f".fai {', '.join(chroms[:3])}); same assembly?")
    unknown = names - set(fai)
    if unknown:
        log.warning("%s: features on %d sequences absent from the .fai ignored (e.g. %s)", what, len(unknown),
                    ", ".join(sorted(unknown)[:3]))


def run(a) -> None:
    t0 = time.time()
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fai = io.read_fai(a.fai)
    names, G, D = track.chromosome_set(fai)
    if not names:
        raise io.InputError(f"no sequence in {a.fai} is long enough to analyse ({params.MIN_LEN_D} x the expected "
                            f"centromere size); chromosome-scale assembly needed")
    chroms = {c: fai[c] for c in names}
    geo = track.geometry(D)
    log.info("version %s: %d chromosomes, genome %.3g Gb, D = %.2f Mb (window %d bins, ring %d-%d bins)",
             __version__, len(names), G / 1e9, D / 1e6, geo.w, geo.lo, geo.hi)
    if D < 2 * params.BIN:
        raise io.InputError(f"expected centromere size D = {D / 1e3:.0f} kb is under two {params.BIN // 1000}-kb bins "
                            f"({len(names)} sequences counted as chromosomes): the assembly is too fragmented or the "
                            f"genome too small for this tool")
    samples = io.read_samples(a.samples)
    for s in samples:
        for ctx, f in s.files.items():
            ctg = io.tabix_contigs(f)
            if not set(ctg) & set(names):
                raise io.InputError(f"sample {s.name} {ctx}: no analysed chromosome in {f} (e.g. "
                                    f"{', '.join(ctg[:3])} vs .fai {', '.join(names[:3])}); same assembly?")
    feats = []
    if a.ltrquest:
        feats.append(io.read_ltrquest(a.ltrquest))
    else:
        log.warning("EDTA transposon classes (superfamilies) are coarser than LTRquest clades")
        feats.append(io.read_edta(a.edta))
    _check_names("transposons", feats[-1].chrom, names, fai)
    feats.append(io.read_satellites(a.satellites))
    _check_names("satellites", feats[-1].chrom, names, fai)
    if a.genes:
        feats.append(io.read_genes(a.genes))
        _check_names("genes", feats[-1].chrom, names, fai)
    anchors, src = {}, "none"
    if a.anchors:
        adf, src = io.read_anchors(a.anchors)
        bad = sorted(set(adf.chrom) - set(fai))
        if bad:
            raise io.InputError(f"anchors on sequences absent from the .fai: {', '.join(bad[:5])}")
        anchors = calibrate.anchor_table(adf, chroms)
        log.info("anchors: %d (%s)", len(anchors), src)
        if not anchors:
            log.warning("no anchor on an analysed chromosome: uncalibrated run (direction down, every candidate "
                        "labelled 'dip')")
    else:
        log.warning("no anchors: uncalibrated run (direction down, every candidate labelled 'dip')")
    excl = io.read_exclusions(a.exclude) if a.exclude else None
    if excl is not None:
        nosample = sorted(set(excl["sample"]) - {""} - {s.name for s in samples})
        if nosample:
            raise io.InputError(f"exclusions {a.exclude}: column 4 names no sample ({', '.join(nosample[:5])}); "
                                f"leave column 4 blank (or '.') to exclude a region in every sample")
        unknown = sorted(set(excl.chrom) - set(fai))
        if unknown:
            log.warning("exclusions on sequences absent from the .fai ignored: %s", ", ".join(unknown[:5]))
    cmap = classes.build(feats, names)
    log.info("sequence classes: %d", len(cmap.names))
    ann = [x for x in (a.ltrquest or a.edta, a.satellites, a.genes) if x]
    key = bins.fingerprint(ann, "|".join(cmap.names))     # same files, other analysed chromosomes -> other classes
    raw = bins.collect(samples, chroms, cmap, key, out.parent / "cenmeth_cache" / key, max(1, a.threads))
    filt, dstats = bins.depth_filter(raw, samples, chroms, excl)
    for d in dstats:
        if not np.isfinite(d.median_depth):
            raise io.InputError(f"sample {d.sample}: no CHG cytosine with >= {params.MIN_CALLS_SITE} calls")
        if d.median_depth < params.LOW_DEPTH:
            log.warning("sample %s: median depth %.1fx < %gx; few repeat cytosines observed, candidates rest on "
                        "unique DNA", d.sample, d.median_depth, params.LOW_DEPTH)
        log.debug("sample %s: median depth %.1f, dropped %d (depth) + %d (excluded) of %d bins", d.sample,
                  d.median_depth, d.dropped_depth, d.dropped_excluded, d.bins_with_data)
    K = len(cmap.names)
    abin = {c: x["bin"] for c, x in anchors.items()}
    cal, direction, dsrc = None, "down", "default"
    if anchors:
        s1 = scan(filt, samples, chroms, K, geo, abin, None)
        cal = calibrate.calibrate(s1.pooled, s1.totals, cmap.names)
        direction, dsrc = calibrate.learned_direction(cal)
    sig = calibrate.signature(cal)
    if anchors and sig.empty:
        if cal.p.isna().all():
            log.warning("too few anchors with data to calibrate (need %d): every candidate labelled 'dip'",
                        params.MIN_ANCHORS)
        else:
            log.warning("no consistent methylation signature at the anchors: every candidate labelled 'dip'")
    sc = scan(filt, samples, chroms, K, geo, abin, direction)
    for s in samples:
        if not any(np.isfinite(y).any() for y in sc.adj[s.name].values()):
            log.warning("sample %s: no usable window (>= %d calls) on any analysed chromosome; no track or candidates",
                        s.name, params.MIN_CALLS)
    cand = candidates.table(sc.picks, sig, cmap.names, anchors, geo, D, chroms)
    if not cand.empty:
        cand.insert(cand.columns.get_loc("D") + 1, "anchor_source", src)
    rep = candidates.anchored(sc, anchors, cand, direction, sig, cmap.names) if anchors else None
    if rep is not None and not rep.empty:
        empty = sorted(set(rep.anchor) - set(rep[rep.calls >= params.MIN_CALLS].anchor))
        if empty:
            log.warning("%d anchors without data (< %d calls in every sample): %s", len(empty), params.MIN_CALLS,
                        ", ".join(empty[:5]))
    sig_txt = ",".join(f"{r.context}:{r.cls}:{r.direction}" for r in sig.itertuples()) or "none"
    header = (f"#cenmeth {__version__} D={round(D)} direction={direction}({dsrc}) anchors={src} "
              f"signature={sig_txt}")
    output.write_bed(cand, header, f"{out}.candidates.bed")
    output.write_tsv(cand, f"{out}.candidates.tsv")
    if rep is not None:
        output.write_tsv(rep, f"{out}.anchors.tsv")
    if cal is not None:
        output.write_tsv(cal, f"{out}.signature.tsv")
    for s in samples:
        output.write_bedgraph(sc.adj[s.name], chroms, f"{out}.{s.name}.adjusted.bedgraph")
    tiers = {}
    for s in samples:
        t = tiers.setdefault(s.name, {})
        for cb in filt[(s.name, "CHG")].values():
            for k, n in zip(cb.cls, cb.n):
                tier_name = ("satellite", "transposon", "gene", "other")[classes.tier(cmap.names[k])]
                t[tier_name] = t.get(tier_name, 0) + int(n)
    n_like = int((cand.label == "centromere_like").sum()) if not cand.empty else 0
    runinfo = {"version": __version__, "chromosomes": len(names), "genome_size": int(G), "D": round(D),
               "window_bp": geo.w * params.BIN, "ring_bp": f"{geo.lo * params.BIN}-{geo.hi * params.BIN}",
               "direction": direction, "direction_source": dsrc, "anchors_source": src, "anchors": len(anchors),
               "classes": K, "signature": sig_txt, "candidates": len(cand), "centromere_like": n_like}
    output.write_tsv(output.qc_rows(runinfo, dstats, tiers), f"{out}.qc.tsv")
    footer = (f"Line = composition-adjusted CHG per sample (one row each; 0 = class reference); boxes = candidates "
              f"(blue = centromere_like, grey = dip); band/triangle = anchor. D = {D / 1e6:.2f} Mb; "
              f"signature: {sig_txt}.")
    plot.overview(sc.adj, chroms, cand, anchors, footer, f"{out}.overview.pdf")
    log.info("done in %.0f s: %d candidates (%d centromere_like); outputs %s.*", time.time() - t0, len(cand),
             n_like, out)


class _Fmt(logging.Formatter):
    def format(self, r):
        return f"cenmeth: {'warning: ' if r.levelno >= logging.WARNING else ''}{r.getMessage()}"


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    h = logging.StreamHandler(sys.stderr)
    h.setFormatter(_Fmt())
    log.handlers[:] = [h]
    log.setLevel(logging.DEBUG if a.verbose else logging.INFO)
    log.propagate = False
    try:
        run(a)
    except io.InputError as e:
        print(f"cenmeth: error: {e}", file=sys.stderr)
        return 2
    return 0
