"""Co-localization: every candidate feature of every line, merged into loci, with the value of each line at each
locus and a kind (docs/design.md, 'Loci').

Members are:
- satellite clusters of every family (>= LOCI_MIN_SAT_BP), with their role and rDNA tag;
- the core of every concentrated LTR clade;
- the anchored and de novo Hi-C peaks, and runs of Hi-C hub bins;
- per-plant methylation candidates that are centromere_like (leave-one-out signature) or plant-specific;
- the sequence call and the integrated call.

- telomere arrays at the chromosome ends, and interstitial arrays of the telomere unit (cenfuse.telomere);
- 45S arrays (NORs) and 5S arrays (cenfuse.rdna).

Satellite families that are dispersed or minor in this genome are background, not features, and are not members
(their 3D association is still tested per family). Members that overlap form one locus; Hi-C members are padded by
one bin (their placement error) for this, and nothing else is, so dense satellite content cannot chain a locus
along a whole chromosome. Its kind is the first that matches:
centromere, sequence_centromere, telomere, interstitial_telomeric, rDNA, subtelomeric_block, interstitial_block,
alt_centromere_candidate, hic_hub, hic_peak, te_core, satellite_block, meth_candidate, other. feature_signatures() then summarizes every line's signal per kind, so
each kind's joint signature (satellite, LTR, Hi-C, methylation) can be read and compared."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pysam

from .params import HIC_TOL_BINS
from .seqcalls import bp_by_role
from .teline import composition

FEATURE_ROLES = ("centromeric", "cen_associated", "interstitial", "subtelomeric", "rDNA")
HIC_TYPES = ("hic_anchored_core", "hic_denovo_core", "hic_hub")
LOCUS_COLS = ["locus", "chrom", "start", "end", "length", "kind", "lines", "n_lines", "members",
              "sat_bp_centromeric", "sat_bp_cen_associated", "sat_bp_interstitial", "sat_bp_subtelomeric",
              "sat_bp_dispersed", "sat_bp_other", "rdna", "telomere_bp", "te_clades", "te_marker_elements",
              "te_elements", "te_top_clades", "hic_anch_max", "hic_denovo_max", "hic_telomere_max", "hic_hub",
              "meth_pooled_extreme", "meth_raw_chg", "meth_samples_centromere_like", "meth_samples_specific",
              "gene_ratio", "gap_bp", "dist_to_call"]
KIND_ORDER = ["centromere", "sequence_centromere", "telomere", "interstitial_telomeric", "rDNA", "subtelomeric_block",
              "interstitial_block", "alt_centromere_candidate", "hic_hub", "hic_peak", "te_core", "satellite_block",
              "meth_candidate", "other"]


def _runs(mask: np.ndarray):
    """(first, last) index of each run of True."""
    d = np.diff(np.concatenate([[0], mask.astype(int), [0]]))
    return list(zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1) - 1))


def members(chrom: str, L: int, fused: dict, clusters: pd.DataFrame, tel: dict, hic: dict | None,
            meth: dict | None, telo: dict | None = None, nors=(), five_s=()) -> list:
    """[(start, end, type, line, label, extra)] candidate features on one chromosome."""
    out = []
    for c, s, e, n, unit in nors or ():
        if c == chrom:
            out.append((s, e, "rdna_array", "satellite", f"45S x{n}", dict(units=n, unit_bp=unit)))
    for c, s, e, n, unit in five_s or ():
        if c == chrom:
            out.append((s, e, "rdna5s_array", "satellite", f"5S x{n}", dict(units=n, unit_bp=unit)))
    for c, side, s, e, bp, gap_end in (telo or {}).get("ends", []):
        if c == chrom:
            out.append((s, e, "telomere_array", "telomere", telo["unit"], dict(side=side, bp=bp, gap_to_end=gap_end)))
    for c, s, e, bp in (telo or {}).get("interstitial", []):
        if c == chrom:
            out.append((s, e, "telomere_interstitial", "telomere", telo["unit"], dict(bp=bp)))
    if fused.get("call"):
        out.append((*fused["call"], "integrated_call", "call", fused["relation"], {}))
    if fused.get("seq"):
        out.append((*fused["seq"], "sequence_call", "satellite", fused["seq_evidence"], {}))
    for r in clusters[clusters.chrom == chrom].itertuples(index=False):
        if r.role in FEATURE_ROLES or r.rdna not in (".", ""):
            out.append((r.start, r.end, "satellite_cluster", "satellite", r.family,
                        dict(role=r.role, rdna=r.rdna, bp=r.bp)))
    marker = tel["marker"][0] if tel.get("marker") else None
    for clade, cores in tel.get("clades", []):
        if chrom in cores:
            c = cores[chrom]
            out.append((int(c["start"]), int(c["end"]), "te_core", "te", clade,
                        dict(marker=clade == marker, k=c["k"], n=c["n"])))
    if hic is not None and chrom in hic["rows"]:
        h = hic["rows"][chrom]
        if h["anch_peak_start"] >= 0:
            out.append((h["anch_core_start"], h["anch_core_end"], "hic_anchored_core", "hic", h["anch_level"],
                        dict(strong=bool(h["anch_strong"]), gate=bool(h["gate_pass"]))))
        if h["denovo_core_start"] >= 0:
            out.append((h["denovo_core_start"], h["denovo_core_end"], "hic_denovo_core", "hic",
                        h["denovo_confidence"], {}))
        k = hic["index"].get(chrom)
        if k is not None:
            b = hic["bins"]
            hub = hic["hub"][b.sl(k)]
            for i, j in _runs(hub):
                out.append((int(b.start[b.off[k] + i]), int(b.end[b.off[k] + j]), "hic_hub", "hic", "hub", {}))
    if meth is not None and not meth["candidates"].empty:
        cd = meth["candidates"]
        cd = cd[(cd.chrom == chrom) & ((cd.loco_label == "centromere_like") | cd.sample_specific)]
        for r in cd.itertuples(index=False):
            out.append((int(r.window_start), int(r.window_end), "meth_candidate", "meth", r.sample,
                        dict(label=r.loco_label, specific=bool(r.sample_specific), p=r.plant_specific_p,
                             value=r.value, rank=r.rank)))
    return [m for m in out if m[1] > m[0]]


def merge(mem: list, tol: int) -> list:
    """Groups of members whose intervals overlap, Hi-C members padded by tol bp (one bin) on each side."""
    pad = lambda m: (m[0] - tol, m[1] + tol) if m[2] in HIC_TYPES else (m[0], m[1])
    mem = sorted(mem, key=lambda m: pad(m))
    groups, cur, cur_end = [], [], None
    for m in mem:
        s, e = pad(m)
        if cur and s > cur_end:
            groups.append(cur)
            cur, cur_end = [], None
        cur.append(m)
        cur_end = e if cur_end is None else max(cur_end, e)
    if cur:
        groups.append(cur)
    return groups


def kind(group: list) -> str:
    types = {m[2] for m in group}
    roles = {m[5].get("role") for m in group if m[2] == "satellite_cluster"}
    rdna = any(m[5].get("rdna", ".") not in (".", "", None) for m in group if m[2] == "satellite_cluster")
    if "integrated_call" in types:
        return "centromere"
    if "sequence_call" in types:
        return "sequence_centromere"
    if "telomere_array" in types:
        return "telomere"
    if "telomere_interstitial" in types:
        return "interstitial_telomeric"
    if rdna or "rDNA" in roles or types & {"rdna_array", "rdna5s_array"}:
        return "rDNA"
    if "subtelomeric" in roles:
        return "subtelomeric_block"
    if "interstitial" in roles:
        return "interstitial_block"
    alt = any(m[2] == "hic_anchored_core" and m[5].get("strong") and m[5].get("gate") for m in group) or \
        any(m[2] == "te_core" and m[5].get("marker") for m in group) or \
        any(m[2] == "meth_candidate" and m[5].get("label") == "centromere_like" for m in group)
    if alt:
        return "alt_centromere_candidate"
    if "hic_hub" in types:
        return "hic_hub"
    if types & {"hic_anchored_core", "hic_denovo_core"}:
        return "hic_peak"
    if "te_core" in types:
        return "te_core"
    if "satellite_cluster" in types:
        return "satellite_block"
    if "meth_candidate" in types:
        return "meth_candidate"
    return "other"


def _describe(m) -> str:
    s, e, typ, line, label, extra = m
    tag = {"satellite_cluster": lambda: f"{label}[{extra.get('role')}]",
           "te_core": lambda: f"{label}{'*' if extra.get('marker') else ''}",
           "meth_candidate": lambda: f"{label}:{extra.get('label')}{':specific' if extra.get('specific') else ''}",
           "hic_anchored_core": lambda: f"anchored:{label}{':strong' if extra.get('strong') else ''}",
           "hic_denovo_core": lambda: f"denovo:{label}"}.get(typ, lambda: label)()
    return f"{typ}:{tag}:{s + 1}-{e}"


def locus_values(chrom, s, e, group, seq, tel, hic, meth, genes, chrom_len, fa) -> dict:
    roles = bp_by_role(seq.arrays, chrom, s, e)
    v = {f"sat_bp_{r}": roles.get(r, 0) for r in ("centromeric", "cen_associated", "interstitial", "subtelomeric",
                                                  "dispersed")}
    v["sat_bp_other"] = sum(b for r, b in roles.items() if r not in ("centromeric", "cen_associated",
                                                                       "interstitial", "subtelomeric", "dispersed"))
    v["rdna"] = ",".join(sorted({m[5]["rdna"] for m in group if m[2] == "satellite_cluster"
                                 and m[5].get("rdna", ".") not in (".", "")}
                                | {f"{'45S' if m[2] == 'rdna_array' else '5S'}_array:{m[5]['units']}x{m[5]['unit_bp']}"
                                   for m in group if m[2] in ("rdna_array", "rdna5s_array")})) or "."
    v["te_clades"] = ",".join(sorted({m[4] for m in group if m[2] == "te_core"})) or "."
    v["telomere_bp"] = sum(m[5]["bp"] for m in group if m[2] in ("telomere_array", "telomere_interstitial"))
    dens = tel.get("density", {}).get(chrom)
    v["te_marker_elements"] = int(dens[1][s // dens[0]: (e - 1) // dens[0] + 1].sum()) if dens else np.nan
    if tel.get("mids"):
        v["te_elements"], v["te_top_clades"] = composition(tel["mids"], chrom, s, e)
    else:
        v["te_elements"], v["te_top_clades"] = np.nan, "."
    for key, track in (("hic_anch_max", "S"), ("hic_denovo_max", "denovo_S"), ("hic_telomere_max", "telomere_S")):
        v[key] = np.nan
        if hic is not None and chrom in hic["index"]:
            k, b = hic["index"][chrom], hic["bins"]
            y = hic[track][b.sl(k)]
            med = np.nanmedian(y) if np.isfinite(y).any() else np.nan
            j0, j1 = s // b.size, (e - 1) // b.size
            seg = y[j0:j1 + 1]
            if med and med > 0 and np.isfinite(seg).any():
                v[key] = round(float(np.nanmax(seg) / med), 3)
    v["hic_hub"] = any(m[2] == "hic_hub" for m in group)
    v["meth_pooled_extreme"] = v["meth_raw_chg"] = np.nan
    v["meth_samples_centromere_like"] = v["meth_samples_specific"] = 0
    if meth is not None and chrom in meth.get("raw_chg", {}):
        from .meth import params as mp
        mm, vv = meth["raw_chg"][chrom]
        j0, j1 = s // mp.BIN, (e - 1) // mp.BIN + 1
        if vv[j0:j1].sum() > 0:
            v["meth_raw_chg"] = round(float(mm[j0:j1].sum() / vv[j0:j1].sum()), 4)
    if meth is not None and chrom in meth["pooled_adj"]:
        from .meth import params as mp
        from .meth.candidates import sign
        y = meth["pooled_adj"][chrom]
        seg = y[s // mp.BIN:(e - 1) // mp.BIN + 1]
        if np.isfinite(seg).any():
            sg = sign(meth["direction"])
            v["meth_pooled_extreme"] = round(float(sg * np.nanmax(sg * seg)), 4)
        mc = [m for m in group if m[2] == "meth_candidate"]
        v["meth_samples_centromere_like"] = len({m[4] for m in mc if m[5].get("label") == "centromere_like"})
        v["meth_samples_specific"] = len({m[4] for m in mc if m[5].get("specific")})
    v["gene_ratio"] = np.nan
    if genes is not None and chrom in genes and len(genes[chrom]) and e > s:
        n = int(np.searchsorted(genes[chrom], e) - np.searchsorted(genes[chrom], s))
        v["gene_ratio"] = round((n / (e - s)) / (len(genes[chrom]) / chrom_len), 3)
    v["gap_bp"] = fa.fetch(chrom, s, e).upper().count("N") if fa is not None else np.nan
    return v


def table(genome: str, chroms: dict, fused: dict, seq, clusters, tel, hic, meth, genes, fasta, telo=None,
          nors=(), five_s=()) -> tuple:
    """(loci DataFrame, members list) for one genome. fused: {chrom: dict(call, seq, relation, seq_evidence)}."""
    rows, feats = [], []
    tol = HIC_TOL_BINS * hic["bins"].size if hic is not None else 0
    with pysam.FastaFile(str(fasta)) as fa:
        for c, L in chroms.items():
            mem = members(c, L, fused.get(c, {}), clusters, tel, hic, meth, telo, nors, five_s)
            call = fused.get(c, {}).get("call")
            for i, grp in enumerate(merge(mem, tol), 1):
                s, e = min(m[0] for m in grp), max(m[1] for m in grp)
                lid = f"{genome}:{c}:L{i}"
                lines = sorted({m[3] for m in grp if m[3] != "call"})
                dist = "." if not call else max(0, max(s, call[0]) - min(e, call[1]))
                rows.append(dict(locus=lid, chrom=c, start=s + 1, end=e, length=e - s, kind=kind(grp),
                                 lines=",".join(lines) or ".", n_lines=len(lines),
                                 members=";".join(_describe(m) for m in grp), dist_to_call=dist,
                                 **locus_values(c, s, e, grp, seq, tel, hic, meth, genes, L, fa)))
                feats += [(lid, m) for m in grp]
    return pd.DataFrame(rows, columns=LOCUS_COLS), feats


SIG_COLS = ["kind", "loci", "chromosomes", "median_length", "sat_bp_centromeric", "sat_bp_interstitial",
            "sat_bp_subtelomeric", "telomere_bp", "te_elements_per_mb", "te_marker_per_mb", "common_te_clade",
            "hic_anch_max", "hic_telomere_max", "hic_hub_share", "meth_pooled_extreme", "meth_raw_chg", "gene_ratio",
            "gap_share", "hic_3d_class", "hic_3d_fold", "hic_3d_p"]
KIND_CLASS = {"centromere": "centromeres (sequence calls)", "sequence_centromere": "centromeres (sequence calls)",
              "telomere": "telomere arrays", "rDNA": "rDNA arrays (all families)",
              "subtelomeric_block": "subtelomeric blocks (all families)",
              "interstitial_block": "interstitial blocks (all families)"}


def feature_signatures(lt: pd.DataFrame, contacts: pd.DataFrame | None) -> pd.DataFrame:
    """One row per locus kind: the median signal of every line (per-Mb densities for TEs), so each kind's joint
    signature can be read, and the 3D self-association of the matching feature class (Hi-C)."""
    rows = []
    fc = {} if contacts is None or contacts.empty else {r.feature_class: r for r in contacts.itertuples(index=False)}
    for kind in KIND_ORDER:
        d = lt[lt.kind == kind]
        if d.empty:
            continue
        mb = d.length / 1e6
        med = lambda col: round(float(np.nanmedian(pd.to_numeric(d[col], errors="coerce"))), 4) \
            if pd.to_numeric(d[col], errors="coerce").notna().any() else np.nan
        top = d.te_top_clades[d.te_top_clades != "."].str.split(",").str[0].str.split(":").str[0]
        cls = fc.get(KIND_CLASS.get(kind, ""))
        rows.append(dict(kind=kind, loci=len(d), chromosomes=d.chrom.nunique(), median_length=int(d.length.median()),
                         sat_bp_centromeric=med("sat_bp_centromeric"), sat_bp_interstitial=med("sat_bp_interstitial"),
                         sat_bp_subtelomeric=med("sat_bp_subtelomeric"), telomere_bp=med("telomere_bp"),
                         te_elements_per_mb=round(float(np.nanmedian(pd.to_numeric(d.te_elements, errors="coerce")
                                                                     / mb)), 2) if d.te_elements.notna().any() else np.nan,
                         te_marker_per_mb=round(float(np.nanmedian(pd.to_numeric(d.te_marker_elements, errors="coerce")
                                                                   / mb)), 2) if d.te_marker_elements.notna().any() else np.nan,
                         common_te_clade=top.mode().iat[0] if len(top) else ".",
                         hic_anch_max=med("hic_anch_max"), hic_telomere_max=med("hic_telomere_max"),
                         hic_hub_share=round(float(d.hic_hub.astype(bool).mean()), 3),
                         meth_pooled_extreme=med("meth_pooled_extreme"), meth_raw_chg=med("meth_raw_chg"),
                         gene_ratio=med("gene_ratio"),
                         gap_share=round(float((pd.to_numeric(d.gap_bp, errors="coerce") / d.length).median()), 4),
                         hic_3d_class=cls.feature_class if cls is not None else ".",
                         hic_3d_fold=cls.fold if cls is not None else np.nan,
                         hic_3d_p=cls.p_enriched if cls is not None else np.nan))
    return pd.DataFrame(rows, columns=SIG_COLS)
