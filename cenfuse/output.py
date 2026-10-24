"""Output files (README, Outputs): integrated calls (TSV, BED, GFF3), loci (TSV, GFF3 of every member), anchored Hi-C
bedGraph, and the run summary. The provenance and plants tables are written by the pipeline (provenance.table,
plants.score)."""
from __future__ import annotations

from datetime import date
from pathlib import Path
from urllib.parse import quote

import numpy as np
import pandas as pd

from . import __version__
from .params import GATE_MIN_REF

RGB = {"high": "0,114,178", "medium": "86,180,233", "low": "153,153,153"}


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


def write_loci(loci: pd.DataFrame, feats: list, genome: str, chroms: dict, out: Path) -> None:
    write_tsv(loci, out / f"{genome}.loci.tsv")
    kinds = dict(zip(loci.locus, loci.kind))
    with open(out / f"{genome}.features.gff3", "w") as f:
        f.write("##gff-version 3\n")
        for ch, L in chroms.items():
            f.write(f"##sequence-region {ch} 1 {L}\n")
        f.write(f"# cenfuse {__version__} {date.today()} genome={genome}: loci (type = kind) and their members\n")
        for r in loci.itertuples(index=False):
            f.write(f"{r.chrom}\tcenfuse\t{r.kind}\t{r.start}\t{r.end}\t.\t.\t.\t"
                    f"{_attrs(dict(ID=r.locus, lines=r.lines, n_lines=r.n_lines))}\n")
        for k, (lid, (s, e, typ, line, label, extra)) in enumerate(feats, 1):
            ch = lid.split(":")[1]
            a = dict(ID=f"{lid}.m{k}", Parent=lid, line=line, label=label, locus_kind=kinds.get(lid))
            a.update({(x if x not in a else f"{line}_{x}"): v for x, v in extra.items()
                      if not isinstance(v, float) or np.isfinite(v)})
            f.write(f"{ch}\tcenfuse\t{typ}\t{s + 1}\t{e}\t.\t.\t.\t{_attrs(a)}\n")


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


def summary_row(genome: str, g, t: pd.DataFrame, hic: dict | None, meth: dict | None, tel: dict,
                telo: dict | None = None, verdicts: dict | None = None, ptab: pd.DataFrame | None = None,
                nors=()) -> dict:
    """One genome: what each line found, whether it was validated, calls by relation and confidence, the inputs'
    provenance and the plants' verdicts."""
    called = t[t.start != "."]
    rel = t.relation.value_counts().to_dict()
    conf = called.confidence.value_counts().to_dict()
    d = dict(genome=genome, chromosomes=len(t), lines=",".join(g.lines), called=len(called),
             high=conf.get("high", 0), medium=conf.get("medium", 0), low=conf.get("low", 0),
             concordant=rel.get("concordant", 0), offset=rel.get("offset", 0), contested=rel.get("contested", 0),
             functional_only=rel.get("functional_only", 0), sequence_only=rel.get("sequence_only", 0),
             none=rel.get("none", 0), seq_called=int((t.seq_start != ".").sum()),
             te_marker=tel["marker"][0] if tel.get("marker") else ".",
             te_decisive=bool(tel["marker"][2]) if tel.get("marker") else ".",
             telomere_unit=(telo or {}).get("unit") or ".", telomere_ends=len((telo or {}).get("ends", [])),
             rdna_arrays=len(nors or ()), rdna_chromosomes=",".join(sorted({x[0] for x in nors or ()})) or ".")
    if hic is not None:
        h = hic["table"]
        d.update(hic_bin=hic["bin_size"], hic_rabl=hic["rabl"]["passed"],
                 hic_gate_pass=int(h.gate_pass.sum()),
                 hic_gate_untested=int((~h.gate_pass & (h.gate_n < GATE_MIN_REF)).sum()), hic_chromosomes=len(h),
                 hic_anch_on_seq=int(h.anch_agrees_seq.sum()), hic_denovo_on_seq=int(h.denovo_agrees_seq.sum()),
                 hic_denovo_on_anch=int(h.denovo_agrees_anch.sum()), hic_seq_called=int(h.seq_start.ge(0).sum()),
                 hic_sequence_conflict=hic["conflict"]["conflict"] if "conflict" in hic else ".")
        fo = hic.get("fold")
        d.update(fold_gate=fo["status"] if fo else ".", fold_on_denovo=f"{fo['k']}/{fo['n']}" if fo else ".",
                 fold_p=f"{fo['p']:.3g}" if fo else ".")
        fc = hic.get("classes")
        if fc is not None and len(fc):
            for name, key in (("centromeres (sequence calls)", "hic_3d_centromeres"), ("chromosome ends", "hic_3d_ends")):
                r = fc[fc.feature_class == name]
                if len(r):
                    d[key] = f"{r.fold.iat[0]:.2f}x P={r.p_enriched.iat[0]:.3g}"
    if meth is not None:
        sig = meth["signature"]
        d.update(meth_samples=len(meth["samples"]), meth_anchors=meth["n_anchors"], meth_direction=meth["direction"],
                 meth_signature=",".join(f"{r.context}:{r.cls}:{r.direction}" for r in sig.itertuples()) or "none",
                 meth_call_support=int((called.meth_call == "centromere_like").sum()),
                 meth_deciding=",".join(meth.get("deciders", [])) or ".")
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
