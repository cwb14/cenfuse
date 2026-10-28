"""Input readers: .fai, samples table, bedMethyl regions, GFF3/BED intervals and the typed annotation readers.
All intervals are returned 0-based half-open."""
from __future__ import annotations

import gzip
import io as _io
import logging
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

import numpy as np
import pandas as pd
import pysam

log = logging.getLogger("cenfuse-meth")
CONTEXTS = ("CHG", "CG", "CHH")
UNCLASSIFIED = {"", ".", "na", "unknown", "mixture", "none"}
EDTA_SKIP = {"repeat_region", "long_terminal_repeat", "target_site_duplication"}
MODKIT_COLS = 18          # modkit pileup bedMethyl: N_valid (column 10), N_mod (12), N_other_mod (14)
ENCODE_COLS = 11          # ENCODE bedMethyl (WGBS pipelines, Megalodon): coverage (10), percent methylated (11); modkit's
                          # first 11 columns are the same layout
FORMATS = "modkit bedMethyl (18 columns) or ENCODE bedMethyl (11: coverage, percent methylated)"


class InputError(Exception):
    """Bad input: the CLI prints 'cenfuse-meth: error: <message>' and exits 2."""


def read_fai(path) -> dict[str, int]:
    """Sequence name -> length, in file order."""
    try:
        d = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1], names=["c", "L"], dtype={"c": str})
    except Exception as e:
        raise InputError(f"cannot read .fai {path}: {e}") from e
    if d.empty or not pd.api.types.is_integer_dtype(d.L):
        raise InputError(f".fai {path} is empty or has non-integer lengths")
    return dict(zip(d.c, d.L.astype(int)))


@dataclass
class Sample:
    name: str
    files: dict  # context ("CHG", "CG", "CHH") -> Path


def read_samples(path, require_chg: bool = True) -> list[Sample]:
    """TSV with columns sample, CHG and optional CG, CHH (paths relative to the table's folder allowed). With
    require_chg False (cenfuse: its CpG line needs only CG), every sample needs CHG or CG."""
    path = Path(path)
    try:
        d = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    except Exception as e:
        raise InputError(f"cannot read samples table {path}: {e}") from e
    need = ("CHG",) if require_chg else ("CHG", "CG")
    if "sample" not in d.columns or not set(need) & set(d.columns):
        raise InputError(f"samples table {path} needs columns 'sample' and {' or '.join(repr(x) for x in need)} "
                         f"(optional 'CG', 'CHH'); found {list(d.columns)}")
    if d.empty:
        raise InputError(f"samples table {path} has no rows")
    bad = [n for n in d["sample"] if not n or "/" in n or any(ch.isspace() for ch in n)]
    if bad:
        raise InputError(f"samples table {path}: bad sample name(s) {bad[:5]}; a sample name must be non-empty, "
                         f"without '/' or whitespace (it names output files)")
    if d["sample"].duplicated().any():
        raise InputError(f"duplicate sample names in {path}: {sorted(set(d['sample'][d['sample'].duplicated()]))}")
    out = []
    for r in d.to_dict("records"):
        files = {}
        for ctx in CONTEXTS:
            v = r.get(ctx, "").strip()
            if not v:
                if ctx == "CHG" and require_chg:
                    raise InputError(f"sample {r['sample']}: no CHG file (CHG is required)")
                continue
            p = Path(v) if Path(v).is_absolute() else path.parent / v
            if not p.exists():
                raise InputError(f"sample {r['sample']}: {ctx} file not found: {p}")
            idx = [x for x in (Path(f"{p}.tbi"), Path(f"{p}.csi")) if x.exists()]
            if not idx:
                raise InputError(f"sample {r['sample']}: {ctx} file has no tabix index (.tbi/.csi): {p} "
                                 f"(make one with: tabix -p bed {p})")
            if idx[0].stat().st_mtime < p.stat().st_mtime:
                log.warning("sample %s: %s index is older than the file (re-index if the file changed): %s",
                            r["sample"], ctx, idx[0])
            files[ctx] = p
        if not files:
            raise InputError(f"sample {r['sample']}: no methylation file (give CHG or CG)")
        out.append(Sample(r["sample"], files))
    return out


def tabix_contigs(path) -> list[str]:
    with pysam.TabixFile(str(path)) as tb:
        return list(tb.contigs)


def bedmethyl_layout(ncol: int) -> str:
    """'modkit' (>= MODKIT_COLS columns), 'encode' (ENCODE_COLS up to modkit's), else 'other'."""
    return "modkit" if ncol >= MODKIT_COLS else "encode" if ncol >= ENCODE_COLS else "other"


def usual_ncol(lines: list) -> int:
    """The field count of most of the first 100 rows (whitespace-separated), so one malformed row cannot decide it."""
    return int(pd.Series([len(x.split()) for x in lines[:100]]).mode().iat[0]) if lines else 0


def table_layout(path) -> tuple[str, int]:
    """(bedmethyl_layout, columns) of a tabixed methylation table, from the first rows that hold data."""
    with pysam.TabixFile(str(path)) as tb:
        for c in tb.contigs:
            lines = [r for _, r in zip(range(100), tb.fetch(c))]
            if lines:
                n = usual_ncol(lines)
                return bedmethyl_layout(n), n
    return "other", 0


def read_bedmethyl(tb: pysam.TabixFile, chrom: str, start: int, end: int, min_calls: int):
    """Cytosines of chrom in [start, end) with N_valid >= min_calls from an open tabix file ->
    (pos, n_valid, n_methylated, n_malformed). modkit: methylated = N_mod + N_other_mod (5mC + 5hmC); rows with
    modification code 'h' are skipped (their counts are already in the 'm' row's N_other_mod); tab- and space-separated
    columns both parse. ENCODE: valid = coverage, methylated = coverage x percent / 100 (rounded). Other tables stop
    with the formats read (InputError)."""
    empty = (np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0, np.int64), 0)
    if chrom not in tb.contigs:
        return empty
    lines = list(tb.fetch(chrom, start, end))
    if not lines:
        return empty
    ncol = usual_ncol(lines)
    layout = bedmethyl_layout(ncol)
    if layout == "other":
        name = tb.filename.decode() if isinstance(tb.filename, bytes) else str(tb.filename)
        raise InputError(f"{name}: rows of {ncol} columns; methylation is read from {FORMATS}")
    if layout == "encode":
        d = pd.read_csv(_io.StringIO("\n".join(lines)), sep=r"\s+", header=None, names=range(ncol),
                        usecols=[1, 9, 10], dtype=float, on_bad_lines="skip")
        ok = (d.notna().all(axis=1) & d[10].between(0, 100)).to_numpy()
        bad = int((~ok).sum()) + (len(lines) - len(d))
        d = d[ok & (d[9] >= min_calls).to_numpy()]
        n = d[9].to_numpy(np.int64)
        return d[1].to_numpy(np.int64), n, np.rint(n * d[10].to_numpy() / 100).astype(np.int64), bad
    d = pd.read_csv(_io.StringIO("\n".join(lines)), sep=r"\s+", header=None, usecols=[1, 3, 9, 11, 13],
                    dtype={1: float, 3: str, 9: float, 11: float, 13: float}, on_bad_lines="skip")
    ok = d.notna().all(axis=1).to_numpy()
    bad = int((~ok).sum()) + (len(lines) - len(d))
    d = d[ok & (d[3] != "h").to_numpy()]
    d = d[d[9] >= min_calls]
    return (d[1].to_numpy(np.int64), d[9].to_numpy(np.int64), (d[11] + d[13]).to_numpy(np.int64), bad)


def _is_gff(path) -> bool:
    name = Path(path).name.lower().removesuffix(".gz")
    return name.endswith((".gff", ".gff3"))


def attr(col: pd.Series, key: str) -> pd.Series:
    """Value of one GFF3 attribute (URL-decoded), '' when absent."""
    v = col.str.extract(rf"(?:^|;)\s*{key}=([^;]*)", expand=False).fillna("")
    return v.map(unquote) if v.str.contains("%").any() else v


def read_intervals(path) -> pd.DataFrame:
    """GFF3 or BED -> DataFrame[chrom, start, end, source, type, name, attrs, thick_start, thick_end] (0-based
    half-open). GFF3 by extension (.gff/.gff3, optionally .gz); everything else is BED. BED name = column 4,
    thick = columns 7-8 when numeric (else NaN: BED6+ files may hold text there)."""
    path = Path(path)
    gff = _is_gff(path)
    try:  # no comment='#': GFF3 attributes may contain '#' (LTRquest family_rep=chr:1-2#LTR/Gypsy/Reina)
        raw = pd.read_csv(path, sep="\t", header=None, dtype=str, names=range(9 if gff else 12),
                          skip_blank_lines=True)
    except Exception as e:
        raise InputError(f"cannot read {path}: {e}") from e
    if gff and (raw[0] == "##FASTA").any():
        raw = raw.iloc[:int(np.argmax((raw[0] == "##FASTA").to_numpy()))]
    raw = raw[~raw[0].fillna("#").str.startswith(("#", "track", "browser"))]
    cols = ["chrom", "start", "end", "source", "type", "name", "attrs", "thick_start", "thick_end"]
    if raw.empty:
        return pd.DataFrame(columns=cols)
    try:
        if gff:
            if raw[8].isna().any():
                raise InputError(f"{path}: GFF3 lines need 9 tab-separated columns")
            a = raw[8]
            name = attr(a, "Name")
            d = pd.DataFrame(dict(chrom=raw[0], start=raw[3].astype(np.int64) - 1, end=raw[4].astype(np.int64),
                                  source=raw[1], type=raw[2], name=name.where(name != "", attr(a, "ID")), attrs=a,
                                  thick_start=np.nan, thick_end=np.nan))
        else:
            if raw[2].isna().any():
                raise InputError(f"{path}: BED lines need at least 3 tab-separated columns")
            d = pd.DataFrame(dict(chrom=raw[0], start=raw[1].astype(np.int64), end=raw[2].astype(np.int64),
                                  source="", type="", name=raw[3].fillna(""), attrs="",
                                  thick_start=pd.to_numeric(raw[6], errors="coerce"),
                                  thick_end=pd.to_numeric(raw[7], errors="coerce")))
    except ValueError as e:
        raise InputError(f"{path}: non-integer coordinates ({e})") from e
    bad = d.end <= d.start
    if bad.any():
        log.warning("%s: %d intervals with end <= start skipped", path, int(bad.sum()))
        d = d[~bad]
    return d.reset_index(drop=True)


def _need(d: pd.DataFrame, path, what: str) -> pd.DataFrame:
    if d.empty:
        raise InputError(f"{path}: no {what} found")
    return d


def read_ltrquest(path) -> pd.DataFrame:
    """LTRquest GFF3 -> DataFrame[chrom, start, end, cls, depth]: one row per LTR_retrotransposon, class 'LTR:<clade>'
    (unclassified clades -> 'LTR:<superfamily>_unclassified' or 'LTR:unclassified'), depth = nesting depth."""
    d = read_intervals(path)
    d = _need(d[d.type == "LTR_retrotransposon"], path, "LTR_retrotransposon features (is this an LTRquest GFF3?)")
    clade, sf = attr(d["attrs"], "clade"), attr(d["attrs"], "superfamily")
    un, sf_un = clade.str.lower().isin(UNCLASSIFIED), sf.str.lower().isin(UNCLASSIFIED)
    cls = np.where(~un, "LTR:" + clade, np.where(~sf_un, "LTR:" + sf + "_unclassified", "LTR:unclassified"))
    depth = pd.to_numeric(attr(d["attrs"], "depth"), errors="coerce").fillna(0).astype(int)
    return pd.DataFrame(dict(chrom=d.chrom, start=d.start, end=d.end, cls=cls, depth=depth)).reset_index(drop=True)


def read_edta(path) -> pd.DataFrame:
    """EDTA *.TEanno.gff3 -> DataFrame[chrom, start, end, cls, depth]: class 'TE:<Classification>' (superfamily
    level, e.g. TE:LTR/Gypsy); container and part features (repeat_region, LTRs, TSDs) skipped; depth 0."""
    d = read_intervals(path)
    d = d[~d.type.isin(EDTA_SKIP)]
    c = attr(d["attrs"], "Classification")
    d, c = d[c != ""], c[c != ""]
    _need(d, path, "features with a Classification= attribute (is this an EDTA TEanno GFF3?)")
    return pd.DataFrame(dict(chrom=d.chrom, start=d.start, end=d.end, cls="TE:" + c, depth=0)).reset_index(drop=True)


def read_satellites(path) -> pd.DataFrame:
    """Tandem-repeat arrays -> DataFrame[chrom, start, end, cls]: class 'sat:<family>'. Satellite-engine
    <genome>_satellite_arrays.gff: satellite_DNA features, family= attribute. Other GFF3: every feature, family from
    family= or Name=. BED: family = name column."""
    d = read_intervals(path)
    if _is_gff(path):
        if (d.type == "satellite_DNA").any():
            d = d[d.type == "satellite_DNA"]
        fam = attr(d["attrs"], "family")
        fam = fam.where(fam != "", d.name.fillna(""))
    else:
        fam = d.name.fillna("")
    d, fam = d[fam != ""], fam[fam != ""]
    _need(d, path, "tandem-repeat arrays with a family name (GFF3 family=/Name= attribute, or BED name column)")
    return pd.DataFrame(dict(chrom=d.chrom, start=d.start, end=d.end, cls="sat:" + fam, depth=0)).reset_index(drop=True)


def read_genes(path) -> pd.DataFrame:
    """GFF3 'gene' features (or every BED interval) -> DataFrame[chrom, start, end, cls='gene']."""
    d = read_intervals(path)
    if _is_gff(path):
        d = d[d.type == "gene"]
    _need(d, path, "gene features")
    return pd.DataFrame(dict(chrom=d.chrom, start=d.start, end=d.end, cls="gene", depth=0)).reset_index(drop=True)


def anchor_source(path) -> str:
    """'hic' for a Hi-C engine BED (first line '#cenfuse-hic'), 'satellite' for a GFF3 with centromere features,
    otherwise 'user' (e.g. CENH3 ChIP domains)."""
    try:
        opener = gzip.open if str(path).lower().endswith(".gz") else open
        with opener(path, "rt") as f:
            first = f.readline()
    except (OSError, UnicodeDecodeError, gzip.BadGzipFile) as e:
        raise InputError(f"cannot read anchors {path}: {e}") from e
    if first.startswith("#cenfuse-hic"):
        return "hic"
    if _is_gff(path):
        d = read_intervals(path)
        if (d.type == "centromere").any():
            return "satellite"
    return "user"


def read_anchors(path) -> tuple[pd.DataFrame, str]:
    """Known centromeres -> (DataFrame[chrom, start, end, mid, name], source). Midpoint = middle of the BED thick
    span when it lies inside the interval (Hi-C engine: the peak bin), else the interval's middle. GFF3 with centromere features: only
    'centromere' features."""
    src = anchor_source(path)
    d = read_intervals(path)
    if src == "satellite":
        d = d[d.type == "centromere"]
    _need(d, path, "anchor intervals")
    ts, te = d.thick_start.to_numpy(float), d.thick_end.to_numpy(float)
    ok = np.isfinite(ts) & np.isfinite(te) & (te > ts) & (ts >= d.start.to_numpy()) & (te <= d.end.to_numpy())
    mid = np.where(ok, (np.nan_to_num(ts) + np.nan_to_num(te)) // 2, (d.start.to_numpy() + d.end.to_numpy()) // 2)
    name = d.name.fillna("").where(d.name.fillna("") != "", d.chrom + "_anchor")
    out = pd.DataFrame(dict(chrom=d.chrom, start=d.start, end=d.end, mid=mid.astype(np.int64), name=name))
    return out.reset_index(drop=True), src


def read_exclusions(path) -> pd.DataFrame:
    """BED of regions to drop -> DataFrame[chrom, start, end, sample]; sample '' (blank, '.', absent) = all."""
    d = read_intervals(path)
    s = d.name.fillna("").astype(str).str.strip().replace(".", "")
    return pd.DataFrame(dict(chrom=d.chrom, start=d.start, end=d.end, sample=s)).reset_index(drop=True)
