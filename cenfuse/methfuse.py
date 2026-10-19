"""The methylation line (cenmeth's machinery, calibrated on the integrated calls).

Anchors are the provisional calls of the other chromosomes. The genome calibration (all anchors) gives the scan
direction and cenmeth's own outputs. Then for each anchored chromosome x, a second calibration leaves x's anchor
out; its signature, sig_x, is the only one ever used to judge places on x. So x's own call never shapes the test of
x. Methylation support for an interval: the most extreme window of the pooled adjusted CHG track inside it,
scored against sig_x; `centromere_like` at a score >= METH_SCORE. Methylation never moves a call alone; a `dip` is
not evidence against (cenmeth).

Plants (docs/plants.md): the support of the reference call pools only the deciding samples (the reference plant's,
or samples of unknown origin). Every sample is also scored at every site of its chromosome (the reference lines'
places, the other libraries' Hi-C peaks, every plant's centromere_like candidates); plants.py compares the plants'
signature scores at each site."""
from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .meth import bins as mbins
from .meth import calibrate, candidates, classes, io, output, params as mp, plot as mplot, track
from .meth.candidates import sign, signature_match
from .meth.scan import scan
from .params import METH_SCORE
from .util import CenfuseError, log

SUPPORT_COLS = ["chrom", "item", "start", "end", "sample", "bin_start", "value", "with_data", "matched", "score",
                "label"]
SITE_COLS = ["chrom", "site", "start", "end", "mid", "origins", "sample", "label", "score", "value", "with_data"]


@dataclass
class MethPrep:
    samples: list
    chroms: dict              # analyzed chromosomes (cenfuse's set, long enough for cenmeth)
    cmap: object
    filt: dict
    dstats: list
    geo: object
    D: float
    G: float
    notes: list = field(default_factory=list)
    coverage: dict = field(default_factory=dict)   # sample -> chrom -> summed CHG N_valid per 50-kb bin, before the
                                                   # depth filter (karyotype test of provenance)


def prepare(g, arrays_gff, cache_dir, threads: int) -> MethPrep | None:
    """cenmeth's inputs, classes and depth-filtered class bins for genome g (GenomeInputs). None (with a warning)
    when the genome cannot be analyzed: methylation is optional evidence, so it never stops the run."""
    try:
        fai = io.read_fai(f"{g.fasta}.fai")                 # the whole assembly: genome size G for D
        names, G, D = track.chromosome_set(fai)
        chroms = {c: fai[c] for c in names if c in g.chroms}
        if not chroms:
            log.warning(f"methylation {g.name}: no chromosome long enough for cenmeth (>= {mp.MIN_LEN_D} D); skipped")
            return None
        if D < 2 * mp.BIN:
            log.warning(f"methylation {g.name}: expected centromere size D = {D / 1e3:.0f} kb is under two "
                        f"{mp.BIN // 1000}-kb bins; skipped")
            return None
        samples = io.read_samples(g.meth)
        for s in samples:
            for ctx, f in s.files.items():
                if not set(io.tabix_contigs(f)) & set(chroms):
                    raise CenfuseError(f"methylation {g.name}: sample {s.name} {ctx} has no analyzed chromosome "
                                       f"({f}); same assembly?")
        feats = [read_tes(g.tes), io.read_satellites(arrays_gff)]
        if g.genes:
            feats.append(io.read_genes(g.genes))
        cmap = classes.build(feats, list(chroms))
        ann = [x for x in (g.tes, arrays_gff, g.genes) if x]
        key = mbins.fingerprint(ann, "|".join(cmap.names))
        raw = mbins.collect(samples, chroms, cmap, key, cache_dir / key, max(1, threads))
        cov = {s.name: {c: np.bincount(raw[(s.name, "CHG")][c].bin, weights=raw[(s.name, "CHG")][c].v,
                                       minlength=mbins.n_bins(L)).astype(float) for c, L in chroms.items()}
               for s in samples}
        excl = io.read_exclusions(g.exclude) if g.exclude else None
        filt, dstats = mbins.depth_filter(raw, samples, chroms, excl)
    except io.InputError as e:
        raise CenfuseError(f"methylation {g.name}: {e}") from None
    for d in dstats:
        if not np.isfinite(d.median_depth):
            raise CenfuseError(f"methylation {g.name}: sample {d.sample} has no CHG cytosine with >= "
                               f"{mp.MIN_CALLS_SITE} calls")
        if d.median_depth < mp.LOW_DEPTH:
            log.warning(f"methylation {g.name}: sample {d.sample} median depth {d.median_depth:.1f}x < "
                        f"{mp.LOW_DEPTH:g}x; few repeat cytosines observed")
    return MethPrep(samples, chroms, cmap, filt, dstats, track.geometry(D), D, G, coverage=cov)


def read_tes(path) -> pd.DataFrame:
    """cenmeth's LTRquest classes (clade, else superfamily_unclassified; nesting depth) from tes/, with the TE line's
    type rule: any feature type containing 'LTR_retrotransposon' (cenrep), so EDTA + TEsorter files work as well."""
    d = io.read_intervals(path)
    d = io._need(d[d.type.str.contains("LTR_retrotransposon", na=False)], path, "LTR_retrotransposon features")
    clade, sf = io.attr(d["attrs"], "clade"), io.attr(d["attrs"], "superfamily")
    un, sf_un = clade.str.lower().isin(io.UNCLASSIFIED), sf.str.lower().isin(io.UNCLASSIFIED)
    cls = np.where(~un, "LTR:" + clade, np.where(~sf_un, "LTR:" + sf + "_unclassified", "LTR:unclassified"))
    depth = pd.to_numeric(io.attr(d["attrs"], "depth"), errors="coerce").fillna(0).astype(int)
    return pd.DataFrame(dict(chrom=d.chrom, start=d.start, end=d.end, cls=cls, depth=depth)).reset_index(drop=True)


def _anchor_frame(anchors: dict, chroms: dict) -> pd.DataFrame:
    rows = [(c, s, e, (s + e) // 2, f"{c}_call") for c, (s, e) in anchors.items() if c in chroms and s >= 0]
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "mid", "name"])


def _class_arrays(prep: MethPrep, c: str, raw: dict | None = None):
    """{sample: {ctx: (num, vok)}} and the pooled {ctx: (num, vok)} (classes x bins) on chromosome c. raw, if given,
    receives the pooled unadjusted CHG sums per bin: raw[c] = (methylated calls, valid calls)."""
    K, nb = len(prep.cmap.names), mbins.n_bins(prep.chroms[c])
    per, pooled = {}, {}
    for smp in prep.samples:
        for ctx in smp.files:
            _, v, m = prep.filt[(smp.name, ctx)][c].dense(K, nb)
            if raw is not None and ctx == "CHG":
                acc = raw.setdefault(c, [np.zeros(nb), np.zeros(nb)])
                acc[0] += m.sum(axis=0)
                acc[1] += v.sum(axis=0)
            num, vok = track.class_sums(v, m, prep.geo)
            per.setdefault(smp.name, {})[ctx] = (num, vok)
            acc = pooled.setdefault(ctx, [np.zeros_like(num), np.zeros_like(vok)])
            acc[0] += num
            acc[1] += vok
    return per, {k: tuple(v) for k, v in pooled.items()}


def _calibrate(prep: MethPrep, samples: list, anc: dict, abin: dict, names, K: int):
    """(calibration, direction, its source, genome signature, {chrom: leave-one-out signature}) from the samples'
    data at the anchors."""
    cal, direction, dsrc, sig, loco = None, "down", "default", calibrate.signature(None), {}
    if len(anc) < mp.MIN_ANCHORS or not samples:
        return cal, direction, dsrc, sig, loco
    s1 = scan(prep.filt, samples, prep.chroms, K, prep.geo, abin, None)
    cal = calibrate.calibrate(s1.pooled, s1.totals, names)
    direction, dsrc = calibrate.learned_direction(cal)
    sig = calibrate.signature(cal)
    mlog = logging.getLogger("cenmeth")
    level = mlog.level
    mlog.setLevel(logging.ERROR)                          # its FDR warning would repeat once per chromosome
    try:
        for x in anc:
            pooled_x = {ctx: {c: v for c, v in per.items() if c != x} for ctx, per in s1.pooled.items()}
            loco[x] = calibrate.signature(calibrate.calibrate(pooled_x, s1.totals, names))
    finally:
        mlog.setLevel(level)
    return cal, direction, dsrc, sig, loco


def _pool(per: dict, names: list) -> dict:
    """{ctx: (num, vok)} summed over the named samples ({} when none has data)."""
    out = {}
    for n in names:
        for ctx, (num, vok) in per.get(n, {}).items():
            acc = out.setdefault(ctx, [np.zeros_like(num), np.zeros_like(vok)])
            acc[0] += num
            acc[1] += vok
    return {k: tuple(v) for k, v in out.items()}


def _signature_bins(arr: dict, sig: pd.DataFrame, names, j0: int, j1: int) -> np.ndarray:
    """Bins j0..j1 where at least one signature entry has MIN_CALLS calls (the score can be computed there)."""
    idx = {n: i for i, n in enumerate(names)}
    ok = np.zeros(j1 - j0 + 1, bool)
    for r in sig.itertuples():
        if r.context in arr and r.cls in idx:
            ok |= arr[r.context][1][idx[r.cls], j0:j1 + 1] >= mp.MIN_CALLS
    return ok


def _best_window(arr: dict, s: int, e: int, nb: int, w: int, direction: str, sig=None, names=None):
    """(bin, adjusted track, signature data there) of the most extreme adjusted CHG window inside [s, e) widened by
    half a window, among the bins where the signature has data when there are any; None without data. Choosing the
    extreme first and asking for signature data afterwards lands low-depth samples on bins without it."""
    if "CHG" not in arr:
        return None
    y = track.adjusted(*arr["CHG"])
    j0, j1 = max(0, s // mp.BIN - w // 2), min(nb - 1, (max(e, s + 1) - 1) // mp.BIN + w // 2)
    seg = sign(direction) * y[j0:j1 + 1]
    fin = np.isfinite(seg)
    if not fin.any():
        return None
    ok = fin & _signature_bins(arr, sig, names, j0, j1) if sig is not None else fin
    pool = ok if ok.any() else fin
    return j0 + int(np.argmax(np.where(pool, seg, -np.inf))), y, bool(ok.any())


def support(arr: dict, sig: pd.DataFrame, names, s: int, e: int, nb: int, prep: MethPrep, direction: str) -> dict:
    hit = _best_window(arr, s, e, nb, prep.geo.w, direction, sig, names)
    if hit is None:
        return dict(bin_start=-1, value=np.nan, with_data=0, matched=0, score=np.nan, label="no_data")
    j, y, _ = hit
    wd, mt, score = signature_match({k: v[0][:, j] for k, v in arr.items()}, {k: v[1][:, j] for k, v in arr.items()},
                                    sig, names)
    lab = "centromere_like" if wd > 0 and score >= METH_SCORE else ("dip" if wd > 0 else "no_signature_data")
    return dict(bin_start=j * mp.BIN, value=float(y[j]), with_data=wd, matched=mt, score=score, label=lab)


def merge_sites(raw: list, W: int, L: int) -> list:
    """Sites of one chromosome from [(origin, start, end)] in priority order: each is judged in a window W centered on
    its midpoint, and a site within W/2 of an earlier one joins it. Returns [dict(mid, start, end, origins)]."""
    groups = []
    for origin, s, e in raw:
        if s < 0 or e <= s:
            continue
        mid = (int(s) + int(e)) // 2
        for gr in groups:
            if abs(mid - gr["mid"]) < W // 2:
                if origin not in gr["origins"]:
                    gr["origins"].append(origin)
                break
        else:
            groups.append(dict(mid=mid, origins=[origin]))
    for gr in groups:
        gr["start"], gr["end"] = max(0, gr["mid"] - W // 2), min(L, gr["mid"] + W // 2)
    return groups


def run(prep: MethPrep, anchors: dict, items: dict, prefix, width_min: int = 0, decide=None, sites=None) -> dict:
    """anchors: {chrom: (start, end)} (provisional calls >= medium); items: {chrom: [(item, start, end)]} the places
    to judge (sequence call, Hi-C peak, TE core, provisional call). Every place is judged in a window of one width
    centered on it, W = max(2 D, width_min) (the Hi-C bin), so a wide interval gets no more chances than a narrow
    one. decide: the samples whose pooled data support the reference call (None: all). sites: {chrom: [(origin,
    start, end)]} more places at which every sample is scored (plants). Writes cenmeth's outputs at prefix.* and
    returns dict(support, candidates, direction, signature, loco, adj, n_anchors, sites)."""
    W = int(max(2 * prep.D, width_min))
    names, K = prep.cmap.names, len(prep.cmap.names)
    anc = calibrate.anchor_table(_anchor_frame(anchors, prep.chroms), prep.chroms)
    abin = {c: a["bin"] for c, a in anc.items()}
    cal, direction, dsrc, sig, loco = _calibrate(prep, prep.samples, anc, abin, names, K)
    if len(anc) >= mp.MIN_ANCHORS:
        log.info(f"methylation: {len(anc)} anchors; genome signature "
                 f"{','.join(f'{r.context}:{r.cls}:{r.direction}' for r in sig.itertuples()) or 'none'}; "
                 f"leave-one-out signatures non-empty on {sum(not v.empty for v in loco.values())}/{len(loco)}")
    else:
        log.warning(f"methylation: {len(anc)} anchors with data (need {mp.MIN_ANCHORS}); uncalibrated, support not "
                    f"assessed")
    sc = scan(prep.filt, prep.samples, prep.chroms, K, prep.geo, abin, direction)
    cand = candidates.table(sc.picks, sig, names, anc, prep.geo, prep.D, prep.chroms)
    if not cand.empty:
        lo = [signature_match(p.num, p.vok, loco.get(p.chrom, sig), names) for p in sc.picks]
        cand.insert(cand.columns.get_loc("label") + 1, "loco_label",
                    ["centromere_like" if wd > 0 and s >= METH_SCORE else "dip" for wd, _, s in lo])
        cand.insert(cand.columns.get_loc("loco_label") + 1, "loco_score", [s for _, _, s in lo])
    deciders = [s.name for s in prep.samples] if decide is None else [n for n in decide]
    if set(deciders) == {s.name for s in prep.samples}:        # the reference call is judged on the deciders only:
        d_direction, d_sig, d_loco = direction, sig, loco      # other plants never shape its signature either
    elif deciders:
        _, d_direction, _, d_sig, d_loco = _calibrate(prep, [s for s in prep.samples if s.name in deciders], anc,
                                                      abin, names, K)
        log.info(f"methylation: signature of the deciding samples ({', '.join(deciders)}): "
                 f"{','.join(f'{r.context}:{r.cls}:{r.direction}' for r in d_sig.itertuples()) or 'none'}")
    else:
        d_direction, d_sig, d_loco = direction, calibrate.signature(None), {}
    cand_sites = {}
    if not cand.empty:
        for r in cand[cand.loco_label == "centromere_like"].itertuples(index=False):
            cand_sites.setdefault(r.chrom, []).append((f"meth:{r.sample}", int(r.window_start), int(r.window_end)))
    rows, site_rows, pooled_adj, raw = [], [], {}, {}
    for c in prep.chroms:
        per, pooled_all = _class_arrays(prep, c, raw)
        pooled = _pool(per, deciders)
        nb = mbins.n_bins(prep.chroms[c])
        pooled_adj[c] = track.adjusted(*pooled_all["CHG"]) if "CHG" in pooled_all else np.full(nb, np.nan)
        sg = loco.get(c, sig)
        dsg = d_loco.get(c, d_sig)
        if sg.empty:
            continue
        for item, s0, e0 in items.get(c, []):
            if s0 < 0:
                continue
            mid = (s0 + e0) // 2
            s, e = max(0, mid - W // 2), min(prep.chroms[c], mid + W // 2)
            if pooled and not dsg.empty:
                rows.append(dict(chrom=c, item=item, start=s, end=e, sample="pooled",
                                 **support(pooled, dsg, names, s, e, nb, prep, d_direction)))
            for smp, arr in per.items():
                rows.append(dict(chrom=c, item=item, start=s, end=e, sample=smp,
                                 **support(arr, sg, names, s, e, nb, prep, direction)))
        raw_sites = list(items.get(c, [])) + list((sites or {}).get(c, [])) + cand_sites.get(c, [])
        for i, gr in enumerate(merge_sites(raw_sites, W, prep.chroms[c]), 1):
            for smp, arr in per.items():
                sp = support(arr, sg, names, gr["start"], gr["end"], nb, prep, direction)
                site_rows.append(dict(chrom=c, site=f"{c}:S{i}", start=gr["start"], end=gr["end"], mid=gr["mid"],
                                      origins=",".join(gr["origins"]), sample=smp, label=sp["label"],
                                      score=sp["score"], value=sp["value"], with_data=sp["with_data"]))
    sup = pd.DataFrame(rows, columns=SUPPORT_COLS)
    _write(prep, prefix, cand, sc, anc, direction, dsrc, sig, cal, loco, names)
    return dict(support=sup, candidates=cand, direction=direction, direction_source=dsrc, signature=sig, loco=loco,
                adj=sc.adj, pooled_adj=pooled_adj, raw_chg={c: tuple(v) for c, v in raw.items()},
                n_anchors=len(anc), D=prep.D, chroms=prep.chroms, samples=[s.name for s in prep.samples],
                deciders=deciders, decide_signature=d_sig, W=W, sites=pd.DataFrame(site_rows, columns=SITE_COLS))


def _write(prep, prefix, cand, sc, anc, direction, dsrc, sig, cal, loco, names) -> None:
    sig_txt = ",".join(f"{r.context}:{r.cls}:{r.direction}" for r in sig.itertuples()) or "none"
    header = (f"#cenmeth (cenfuse) D={round(prep.D)} direction={direction}({dsrc}) anchors=cenfuse "
              f"signature={sig_txt}")
    output.write_bed(cand, header, f"{prefix}.candidates.bed")
    output.write_tsv(cand, f"{prefix}.candidates.tsv")
    if anc:
        output.write_tsv(candidates.anchored(sc, anc, cand, direction, sig, names), f"{prefix}.anchors.tsv")
    if cal is not None:
        output.write_tsv(cal, f"{prefix}.signature.tsv")
        rows = [dict(left_out=x, context=r.context, cls=r.cls, mean_dev=r.mean_dev, direction=r.direction)
                for x, s in loco.items() for r in s.itertuples()]
        output.write_tsv(pd.DataFrame(rows, columns=["left_out", "context", "cls", "mean_dev", "direction"]),
                         f"{prefix}.loco_signatures.tsv")
    for s in prep.samples:
        output.write_bedgraph(sc.adj[s.name], prep.chroms, f"{prefix}.{s.name}.adjusted.bedgraph")
    run = {"D": round(prep.D), "genome_size": int(prep.G), "chromosomes": len(prep.chroms), "direction": direction,
           "direction_source": dsrc, "anchors": len(anc), "classes": len(names), "signature": sig_txt,
           "candidates": len(cand)}
    tiers = {}
    for s in prep.samples:
        t = tiers.setdefault(s.name, {})
        for cb in prep.filt[(s.name, "CHG")].values():
            for k, n in zip(cb.cls, cb.n):
                nm = ("satellite", "transposon", "gene", "other")[classes.tier(names[k])]
                t[nm] = t.get(nm, 0) + int(n)
    output.write_tsv(output.qc_rows(run, prep.dstats, tiers), f"{prefix}.qc.tsv")
    footer = (f"cenfuse anchors (provisional calls). Line = composition-adjusted CHG per sample; boxes = candidates "
              f"(blue = centromere_like, gray = dip). D = {prep.D / 1e6:.2f} Mb; signature: {sig_txt}.")
    mplot.overview(sc.adj, prep.chroms, cand, anc, footer, f"{prefix}.overview.pdf")
