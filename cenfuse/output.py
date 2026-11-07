"""Output files (README, Outputs): integrated calls (TSV, BED, GFF3), loci (TSV, GFF3 of every member), anchored Hi-C
bedGraph, and the run summary. The provenance and plants tables are written by the pipeline (provenance.table,
plants.score)."""
from __future__ import annotations

from datetime import date
from pathlib import Path
from urllib.parse import quote

import numpy as np
import pandas as pd

from . import __version__, morphology
from .cdr import DIP_P
from .params import AXIS_LEAN_P, GATE_MIN_REF
from .sat.karyotype import sentence

RGB = {"high": "0,114,178", "medium": "86,180,233", "low": "153,153,153"}
# GFF3 types are Sequence Ontology terms (validators and annotation tools check them); cenfuse's own kind and member
# type are kept as the kind= and member= attributes, and a member names its locus with locus= (GFF3's Parent= means
# "part of" in SO, and a satellite block is no part of a centromere). The CpG kinetochore is SO's
# regional_centromere_central_core, "the region of a modular centromere where the kinetochore is formed".
SO_KIND = {"centromere": "centromere", "sequence_centromere": "biological_region", "long_unit_array": "tandem_repeat",
           "telomere": "telomere",
           "interstitial_telomeric": "telomeric_repeat", "rDNA": "tandem_repeat", "subtelomeric_block": "satellite_DNA",
           "interstitial_block": "satellite_DNA", "alt_centromere_candidate": "biological_region",
           "hic_hub": "biological_region", "hic_peak": "biological_region", "te_core": "repeat_region",
           "satellite_block": "satellite_DNA", "meth_candidate": "epigenetically_modified_region",
           "other": "biological_region"}
SO_MEMBER = {"rdna_array": "tandem_repeat", "rdna5s_array": "tandem_repeat", "long_unit_array": "tandem_repeat",
             "telomere_array": "telomeric_repeat",
             "telomere_interstitial": "telomeric_repeat", "integrated_call": "centromere",
             "sequence_call": "biological_region", "satellite_cluster": "satellite_DNA", "te_core": "repeat_region",
             "hic_anchored_core": "biological_region", "hic_denovo_core": "biological_region",
             "hic_hub": "biological_region", "meth_candidate": "epigenetically_modified_region"}


def write_tsv(df: pd.DataFrame, path) -> None:
    df.to_csv(path, sep="\t", index=False, float_format="%.6g", na_rep=".")


def _attrs(d: dict) -> str:
    return ";".join(f"{k}={quote(str(v), safe=' ,:+*[]._-')}" for k, v in d.items() if v not in (None, "", "."))


def write_calls(t: pd.DataFrame, genome: str, chroms: dict, out: Path) -> None:
    write_tsv(t, out / f"{genome}.centromeres.tsv")
    c = t[t.start != "."]
    with open(out / f"{genome}.centromeres.bed", "w") as f:
        f.write(f"#cenfuse {__version__} {genome}: integrated centromere calls (score: high 1000, medium 600, "
                f"low 300)\n")
        for r in c.itertuples(index=False):
            sc = {"high": 1000, "medium": 600, "low": 300}.get(r.confidence, 0)
            f.write(f"{r.chrom}\t{int(r.start) - 1}\t{r.end}\t{r.relation}:{r.confidence}\t{sc}\t.\t"
                    f"{int(r.start) - 1}\t{r.end}\t{RGB.get(r.confidence, '0,0,0')}\n")
    with open(out / f"{genome}.centromeres.gff3", "w") as f:
        f.write("##gff-version 3\n")
        for ch, L in chroms.items():
            f.write(f"##sequence-region {ch} 1 {L}\n")
        f.write(f"# cenfuse {__version__} {date.today()} genome={genome}\n")
        for r in c.itertuples(index=False):
            a = dict(ID=f"cen_{r.chrom}", relation=r.relation, confidence=r.confidence, source=r.call_source,
                     support=r.support, against=r.against, flags=r.flags, sequence_call=(
                         f"{r.seq_start}-{r.seq_end}" if r.seq_start != "." else "."),
                     sequence_evidence=r.seq_evidence, main_family=r.main_family, offset_bp=r.offset_bp)
            f.write(f"{r.chrom}\tcenfuse\tcentromere\t{r.start}\t{r.end}\t.\t.\t.\t{_attrs(a)}\n")
            if getattr(r, "cpg_kt_start", ".") != ".":           # the CpG line's kinetochore inside the call
                k = dict(ID=f"kt_{r.chrom}", Parent=f"cen_{r.chrom}", Name="kinetochore", dip_p=r.cpg_dip_p,
                         kt_meth=r.cpg_kt_meth, array_meth=r.cpg_array_meth)
                f.write(f"{r.chrom}\tcenfuse\tregional_centromere_central_core\t{r.cpg_kt_start}\t{r.cpg_kt_end}"
                        f"\t.\t.\t.\t{_attrs(k)}\n")


def write_loci(loci: pd.DataFrame, feats: list, genome: str, chroms: dict, out: Path) -> None:
    write_tsv(loci, out / f"{genome}.loci.tsv")
    kinds, chrom_of = dict(zip(loci.locus, loci.kind)), dict(zip(loci.locus, loci.chrom))
    with open(out / f"{genome}.features.gff3", "w") as f:
        f.write("##gff-version 3\n")
        for ch, L in chroms.items():
            f.write(f"##sequence-region {ch} 1 {L}\n")
        f.write(f"# cenfuse {__version__} {date.today()} genome={genome}: loci (type = kind) and their members\n")
        for r in loci.itertuples(index=False):
            f.write(f"{r.chrom}\tcenfuse\t{SO_KIND.get(r.kind, 'biological_region')}\t{r.start}\t{r.end}\t.\t.\t.\t"
                    f"{_attrs(dict(ID=r.locus, Name=r.kind, kind=r.kind, lines=r.lines, n_lines=r.n_lines))}\n")
        for k, (lid, (s, e, typ, line, label, extra)) in enumerate(feats, 1):
            ch = chrom_of[lid]                               # not parsed from the ID: names may hold ':'
            a = dict(ID=f"{lid}.m{k}", locus=lid, Name=f"{typ}:{label}", member=typ, line=line, label=label,
                     locus_kind=kinds.get(lid))      # locus=, not Parent=: a member lies in the locus, it is no SO part
            a.update({(x if x not in a else f"{line}_{x}"): v for x, v in extra.items()
                      if not isinstance(v, float) or np.isfinite(v)})
            so = SO_MEMBER.get(typ, "biological_region")
            f.write(f"{ch}\tcenfuse\t{so}\t{s + 1}\t{e}\t.\t.\t.\t{_attrs(a)}\n")


def write_bedgraph(S, bins, members, path, title: str) -> None:
    """Score per bin over its chromosome median (callable chromosomes; NaN bins left out)."""
    with open(path, "w") as f:
        f.write(f'track type=bedGraph name="{title}"\n')
        for k in np.where(members)[0]:
            sl = bins.sl(int(k))
            y = S[sl]
            if not np.isfinite(y).any():
                continue
            med = np.nanmedian(y)
            if not med > 0:
                continue
            for j in np.flatnonzero(np.isfinite(y)):
                g = int(bins.off[k] + j)
                f.write(f"{bins.names[k]}\t{bins.start[g]}\t{bins.end[g]}\t{y[j] / med:.4g}\n")


def _p(x) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return 1.0


def architecture(d: dict) -> str:
    """What the deciding Hi-C library showed about the nucleus (0.10.1; reporting only, from tests already run):
    rabl_axis = the arm fold, or the trans arm alignment with the fold leaning the same way (fuse.axis_line),
    validated the prior-free peaks (centromeres clustered at one pole, arms aligned: Type I, Hoencamp et al. 2021);
    centromere_clustering = the anchored gate passed (the centromeres meet in trans on most chromosomes);
    weak_centromere_clustering = only the support tier passed (beyond chance, not on most: territorial nuclei);
    untested = no gate references and no axis line; not_detected = none of these (a territorial or holocentric
    nucleus, or sequence calls wrong genome-wide: Hi-C alone cannot tell which)."""
    n = int(d.get("hic_chromosomes") or 0)
    alone = bool(n) and int(d.get("hic_gate_untested") or 0) == n
    if "hic_axis" in d:                                    # the decision fuse.genome_axis made (not re-derived from
        if d["hic_axis"] in ("fold", "align"):            # rounded P values)
            return "rabl_axis"
    elif d.get("fold_gate") == "pass" or (d.get("align_gate") == "pass" and _p(d.get("fold_p")) <= AXIS_LEAN_P):
        return "rabl_axis"
    if int(d.get("hic_gate_pass") or 0):
        return "centromere_clustering"
    if int(d.get("hic_gate_support") or 0):
        return "weak_centromere_clustering"
    return "untested" if alone else "not_detected"


def summary_row(genome: str, g, t: pd.DataFrame, hic: dict | None, meth: dict | None, tel: dict,
                telo: dict | None = None, verdicts: dict | None = None, ptab: pd.DataFrame | None = None,
                nors=(), cpg: dict | None = None, five_s=(), at: dict | None = None, sat_test: dict | None = None) -> dict:
    """One genome: what each line found, whether it was validated, calls by relation and confidence, the inputs'
    provenance and the plants' verdicts."""
    called = t[t.start != "."]
    rel = t.relation.value_counts().to_dict()
    conf = called.confidence.value_counts().to_dict()
    d = dict(genome=genome, chromosomes=len(t), lines=",".join(g.lines), called=len(called),
             high=conf.get("high", 0), medium=conf.get("medium", 0), low=conf.get("low", 0),
             concordant=rel.get("concordant", 0), offset=rel.get("offset", 0), contested=rel.get("contested", 0),
             functional_only=rel.get("functional_only", 0), sequence_only=rel.get("sequence_only", 0),
             none=rel.get("none", 0),
             seq_called=int(((t.seq_start != ".") & ~t["flags"].astype(str).str.contains("long_unit_unsupported")).sum()),
             partial_chain_calls=int(t.seq_flags.astype(str).str.contains("partial_chain").sum()),
             sat_one_per_chromosome=sentence(sat_test) if sat_test else ".",
             te_marker=tel["marker"][0] if tel.get("marker") else ".",
             te_decisive=bool(tel["marker"][2]) if tel.get("marker") else ".",
             te_seed=tel["signature"].seed if tel.get("signature") is not None else ".",
             te_specificity=round(tel["signature"].specificity, 1) if tel.get("signature") is not None else ".",
             te_verdict=tel["signature"].why if tel.get("signature") is not None else ".",
             at_gate=(at["status"] if at else "."), at_window=(at.get("window") or "." if at else "."),
             at_p=(f"{at['p']:.3g}" if at else "."),
             at_hic=(f"{at['hic']['k']}/{at['hic']['n']} {at['hic']['status']}" if at and at.get("hic") else "."),
             telomere_unit=(telo or {}).get("unit") or ".", telomere_ends=len((telo or {}).get("ends", [])),
             telomere_interstitial=len((telo or {}).get("interstitial", [])),
             rdna_arrays=len(nors or ()), rdna_chromosomes=",".join(sorted({x[0] for x in nors or ()})) or ".",
             rdna5s_arrays=len(five_s or ()), rdna5s_chromosomes=",".join(sorted({x[0] for x in five_s or ()})) or ".")
    if "cen_point" in t:
        d.update(morphology.summary(t))
    if hic is not None:
        h = hic["table"]
        d.update(hic_bin=hic["bin_size"], hic_rabl=hic["rabl"]["passed"],
                 hic_gate_pass=int(h.gate_pass.sum()),
                 hic_gate_support=int(h.gate_support.sum()) if "gate_support" in h else 0,
                 hic_gate_untested=int((~h.gate_pass & (h.gate_n < GATE_MIN_REF)).sum()),
                 hic_gate_refs=h.gate_refs.iat[0] if "gate_refs" in h and len(h) else ".", hic_chromosomes=len(h),
                 hic_anch_on_seq=int(h.anch_agrees_seq.sum()), hic_denovo_on_seq=int(h.denovo_agrees_seq.sum()),
                 hic_denovo_on_anch=int(h.denovo_agrees_anch.sum()), hic_seq_called=int(h.seq_start.ge(0).sum()),
                 hic_sequence_conflict=hic["conflict"]["conflict"] if "conflict" in hic else ".",
                 te_hic=(f"{hic['te_anchored']['k']}/{hic['te_anchored']['n']} chromosomes (P "
                         f"{hic['te_anchored']['p']:.2g}): {hic['te_anchored']['status']}") if hic.get("te_anchored") else ".",
                 cis_line=(f"{'peak' if hic['cis']['sign'] > 0 else 'dip'} at the core on {hic['cis']['k']}/{hic['cis']['n']} "
                           f"chromosomes (P {hic['cis']['p']:.2g}): {hic['cis']['status']}") if hic.get("cis") else ".")
        fo = hic.get("fold")
        d.update(fold_gate=fo["status"] if fo else ".", fold_on_denovo=f"{fo['k']}/{fo['n']}" if fo else ".",
                 fold_p=f"{fo['p']:.3g}" if fo else ".")
        al = (fo or {}).get("align")
        d.update(align_gate=al["status"] if al else ".", align_on_denovo=f"{al['k']}/{al['n']}" if al else ".",
                 align_p=f"{al['p']:.3g}" if al else ".", hic_axis=hic.get("axis") or ".")
        fc = hic.get("classes")
        if fc is not None and len(fc):
            for name, key in (("centromeres (sequence calls)", "hic_3d_centromeres"), ("chromosome ends", "hic_3d_ends")):
                r = fc[fc.feature_class == name]
                if len(r):
                    d[key] = f"{r.fold.iat[0]:.2f}x P={r.p_enriched.iat[0]:.3g}"
        d["hic_architecture"] = architecture(d)
    if meth is not None:
        sig = meth["signature"]
        d.update(meth_samples=len(meth["samples"]), meth_anchors=meth["n_anchors"], meth_direction=meth["direction"],
                 meth_signature=",".join(f"{r.context}:{r.cls}:{r.direction}" for r in sig.itertuples()) or "none",
                 meth_call_support=int((called.meth_call == "centromere_like").sum()),
                 meth_deciding=",".join(meth.get("deciders", [])) or ".")
    if cpg is not None:
        c = cpg["table"]
        d.update(cpg_gate_pass=int(c.gate_pass.sum()),
                 cpg_dips=f"{int((c.dip_p <= DIP_P).sum())}/{int(c.dip_p.notna().sum())}",
                 cpg_positional=f"{int((c.pos_p <= DIP_P).sum())}/{int(c.pos_p.notna().sum())}",
                 cpg_kinetochores=int((c.kt_start >= 0).sum()))
    if verdicts:
        for v in ("reference", "other", "unknown"):
            d[f"inputs_{v}"] = ",".join(sorted(f"{k}:{n}" for (k, n), x in verdicts.items() if x["verdict"] == v)) or "."
        d["hic_library"] = t.hic_library.iat[0] if "hic_library" in t and len(t) else "."
    if ptab is not None and len(ptab):
        d["plants"] = ptab.plant.nunique()
        for v in ("shifted", "possible_shift", "call_absent", "conflict"):
            x = ptab[ptab.verdict == v]
            d[f"plants_{v}"] = ",".join(f"{r.plant}:{r.chrom}" for r in x.itertuples()) or "."
    return d
