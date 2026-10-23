"""The sequence line: the satellite engine's outputs read back. They come from this run (OUT/satellite) or
from an earlier run (--satellite-dir). Coordinates are converted to 0-based half-open."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

import numpy as np
import pandas as pd

from .params import LOCI_CLUSTER_GAP, LOCI_MIN_SAT_BP
from .sat.classify import make_clusters
from .util import CenfuseError, log

CALL_COLS = ["chrom", "chrom_len", "start", "end", "evidence", "confidence", "flags", "main_family", "sat_bp",
             "te_clade", "te_start", "te_end", "alt_start", "alt_end"]


@dataclass
class SeqLine:
    genome: str
    calls: pd.DataFrame        # CALL_COLS; start/end -1 without a call, te_start/te_end -1 without a core
    arrays: pd.DataFrame       # chrom, start, end, family, role, rdna, monomer_len, n_monomers, in_centromere
    families: pd.DataFrame     # satellite_families.tsv rows of this genome (empty when absent)
    source: Path

    def call(self, chrom: str):
        r = self.calls.loc[self.calls.chrom == chrom]
        return None if r.empty or r.iloc[0].start < 0 else r.iloc[0]


def gff_attrs(col: str) -> dict:
    out = {}
    for kv in str(col).split(";"):
        k, sep, v = kv.partition("=")
        if sep:
            out[k.strip()] = unquote(v.strip())
    return out


def _int(v, minus: int = 0) -> int:
    try:
        return int(v) - minus
    except (TypeError, ValueError):
        return -1


def read_calls(path: Path, chroms: dict) -> pd.DataFrame:
    d = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    need = {"chrom", "cen_start", "cen_end", "evidence", "confidence"}
    if not need <= set(d.columns):
        raise CenfuseError(f"{path}: not a satellite-engine centromere summary (needs {', '.join(sorted(need))})")
    rows = []
    by = {r["chrom"]: r for r in d.to_dict("records")}
    for c, L in chroms.items():
        r = by.get(c)
        if r is None:
            rows.append(dict(chrom=c, chrom_len=L, start=-1, end=-1, evidence="none", confidence="none", flags=".",
                             main_family=".", sat_bp=0, te_clade=".", te_start=-1, te_end=-1, alt_start=-1,
                             alt_end=-1))
            continue
        call = r["evidence"] != "none"
        rows.append(dict(chrom=c, chrom_len=L, start=_int(r["cen_start"], 1) if call else -1,
                         end=_int(r["cen_end"]) if call else -1, evidence=r["evidence"],
                         confidence=r["confidence"] if call else "none", flags=r.get("flags", ".") or ".",
                         main_family=r.get("main_family", "."), sat_bp=max(0, _int(r.get("cen_sat_bp"))),
                         te_clade=r.get("te_clade", ".") or ".", te_start=_int(r.get("te_core_start"), 1),
                         te_end=_int(r.get("te_core_end")), alt_start=_int(r.get("alt_start"), 1) if call else -1,
                         alt_end=_int(r.get("alt_end")) if call else -1))
    missing = sorted(set(chroms) - set(by))
    if missing:
        log.warning(f"{path.name}: no row for {len(missing)} chromosomes (e.g. {', '.join(missing[:3])}): no "
                    f"sequence call there")
    return pd.DataFrame(rows, columns=CALL_COLS)


def read_arrays(path: Path, chroms: dict) -> pd.DataFrame:
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "satellite_DNA" or p[0] not in chroms:
                continue
            a = gff_attrs(p[8])
            rows.append((p[0], int(p[3]) - 1, int(p[4]), a.get("family", "."), a.get("role", "."),
                         a.get("rdna", "."), _int(a.get("monomer_len")), _int(a.get("n_monomers")),
                         a.get("in_centromere", "false") == "true"))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "family", "role", "rdna", "monomer_len",
                                       "n_monomers", "in_centromere"])


def read(sat_dir, genome: str, chroms: dict) -> SeqLine:
    sat_dir = Path(sat_dir)
    summ, arr = sat_dir / f"{genome}_centromere_summary.tsv", sat_dir / f"{genome}_satellite_arrays.gff"
    for p in (summ, arr):
        if not p.exists():
            raise CenfuseError(f"missing satellite-engine output {p}")
    fam_path = sat_dir / "satellite_families.tsv"
    fams = pd.read_csv(fam_path, sep="\t", dtype=str, keep_default_na=False) if fam_path.exists() else pd.DataFrame()
    if "genome" in fams:
        fams = fams[fams.genome == genome].reset_index(drop=True)
    return SeqLine(genome, read_calls(summ, chroms), read_arrays(arr, chroms), fams, sat_dir)


def satellite_clusters(arrays: pd.DataFrame) -> pd.DataFrame:
    """Clusters of each family's arrays per chromosome (the satellite engine's cluster rule), at least LOCI_MIN_SAT_BP of
    satellite: chrom, start, end, bp, family, role, rdna."""
    out = []
    for (c, fam), d in arrays.groupby(["chrom", "family"], sort=False):
        role = d.role.mode().iat[0] if len(d) else "."
        rd = next((x for x in d.rdna if x not in (".", "")), ".")
        for cl in make_clusters(list(zip(d.start, d.end)), LOCI_CLUSTER_GAP):
            if cl.bp >= LOCI_MIN_SAT_BP:
                out.append((c, int(cl.start), int(cl.end), int(cl.bp), fam, role, rd))
    return pd.DataFrame(out, columns=["chrom", "start", "end", "bp", "family", "role", "rdna"])


def bp_by_role(arrays: pd.DataFrame, chrom: str, s: int, e: int) -> dict:
    """Satellite bp inside [s, e) per family role (arrays clipped to the interval)."""
    a = arrays[(arrays.chrom == chrom) & (arrays.end > s) & (arrays.start < e)]
    if a.empty:
        return {}
    clip = (np.minimum(a.end, e) - np.maximum(a.start, s)).clip(lower=0)
    return {k: int(v) for k, v in clip.groupby(a.role).sum().items() if v > 0}
