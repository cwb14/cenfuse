"""Per-plant scoring (docs/plants.md): does each plant's centromere sit at the reference call, or elsewhere?

Methylation: every plant is scored at the same sites of each chromosome (methfuse.run). A site is at the call when
its window overlaps the call. A plant is at the call if it is centromere_like at such a site; a shift needs a dip
there, centromere_like at a site away from the call, and between-plant support: the plant is less centromeric than
the others at the call, or more at that site (Benjamini-Hochberg over every test of the genome).
The between-plant statistic d is the plant's signature score minus the median of the other plants' scores at the
same site (the score defines centromere_like, and comparing plants on the same sequence cancels composition). Its
null on chromosome x is its spread at the calls of the other chromosomes, where plants mostly share the centromere
(leave one chromosome out): the plant's own values, or every plant's when the plant has fewer than MIN_NULL. A new
value is judged against those n values with Student's t prediction interval, t = (d - mean) / (sd sqrt(1 + 1/n)) on
n - 1 degrees of freedom, which carries the uncertainty of a spread estimated from about a dozen chromosomes. (The
synthetic genomes showed why each choice matters: d from the extreme adjusted-CHG bin, a null from the arms, or a
normal tail on a median-absolute-deviation scale from nine values all made shared centromeres look lost.) A gain is tested only at sites that do not come from the plant's own candidates
alone, which were selected on its data.
Hi-C: a library that does not decide is tested with partner windows at the reference calls and its own gate. A far
strong peak is a shift only where the deciding library's peak is elsewhere: where both libraries peak at the same
place, the plant shares the reference plant's functional centromere (even if the reference call is contested).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

from .meth.calibrate import bh
from .params import HIC_FAR_BINS, HIC_TOL_BINS, PLANT_Q
from .util import gap

PLANT_COLS = ["plant", "provenance", "chrom", "verdict", "meth_verdict", "call_site", "call_label", "call_score",
              "loss_q", "alt_site", "alt_start", "alt_end", "alt_origins", "alt_label", "alt_score", "gain_q",
              "shift_bp", "hic_library", "hic_verdict", "hic_peak_start", "hic_peak_end"]
NO_DATA = ("no_data", "no_signature_data")
MIN_NULL = 5             # call-site values needed for the null spread (4 degrees of freedom); heuristic


def _own_only(origins: str, plant: str) -> bool:
    return set(origins.split(",")) == {f"meth:{plant}"}


def site_tests(sites: pd.DataFrame, calls: dict) -> pd.DataFrame:
    """The site table with, per row: inside (the site's window overlaps the reference call), d and its z against the
    call-site null, the test it carries ('loss' at the plant's call site, 'gain' at sites away from the call that do
    not come only from the plant's own candidates) and its Benjamini-Hochberg q over every test."""
    s = sites.reset_index(drop=True)
    s["inside"] = [c in calls and calls[c] is not None and a < calls[c][1] and b > calls[c][0]
                   for c, a, b in zip(s.chrom, s.start, s.end)]
    sc = pd.to_numeric(s.score, errors="coerce").where(~s.label.isin(NO_DATA))
    s["plants_with_data"] = sc.notna().groupby([s.chrom, s.site]).transform("sum").astype(int)
    scv = sc.to_numpy(float)
    d = np.full(len(s), np.nan)
    for idx in s.groupby(["chrom", "site"]).indices.values():          # positions of one site's rows
        v = scv[idx]
        fin = np.isfinite(v)
        if fin.sum() < 3:
            continue
        for a in np.flatnonzero(fin):
            o = np.delete(v, a)
            d[idx[a]] = v[a] - float(np.median(o[np.isfinite(o)]))
    s["d"] = d
    s["call_site"] = False
    for (c, smp), d in s[s.inside].groupby(["chrom", "sample"]):
        best = sc.loc[d.index]
        i = best.idxmax() if best.notna().any() else (d.mid - sum(calls[c]) / 2).abs().idxmin()
        s.loc[i, "call_site"] = True
    s["t"], s["df"] = np.nan, np.nan
    ref = s[s.call_site & s.d.notna()]
    for (c, smp), idx in s.groupby(["chrom", "sample"]).groups.items():
        v = ref[(ref["sample"] == smp) & (ref.chrom != c)].d.to_numpy(float)
        if len(v) < MIN_NULL:
            v = ref[ref.chrom != c].d.to_numpy(float)
        if len(v) < MIN_NULL:
            continue
        sd = float(np.std(v, ddof=1))
        if sd > 0:
            s.loc[idx, "t"] = (s.loc[idx, "d"] - float(np.mean(v))) / (sd * np.sqrt(1 + 1 / len(v)))
            s.loc[idx, "df"] = len(v) - 1
    s["p_gain"] = student_t.sf(s.t, s.df)
    s["p_loss"] = student_t.cdf(s.t, s.df)
    s["test"] = ""
    ok = (s.plants_with_data >= 3) & np.array([calls.get(c) is not None for c in s.chrom], bool)
    s.loc[s.call_site & ok & s.p_loss.notna(), "test"] = "loss"
    gain = (~s.inside) & ok & s.p_gain.notna() & np.array([not _own_only(o, p) for o, p in
                                                           zip(s.origins, s["sample"])], bool)
    s.loc[gain, "test"] = "gain"
    p = np.where(s.test == "loss", s.p_loss, np.where(s.test == "gain", s.p_gain, np.nan)).astype(float)
    s["q"] = bh(p)
    return s


def _meth_empty(verdict: str = "no_data") -> dict:
    return dict(meth_verdict=verdict, call_site=".", call_label=".", call_score=np.nan, loss_q=np.nan, alt_site=".",
                alt_start=-1, alt_end=-1, alt_origins=".", alt_label=".", alt_score=np.nan, gain_q=np.nan)


def meth_verdict(d: pd.DataFrame, absent: bool) -> dict:
    """One plant x chromosome from its rows of the tested site table."""
    out = _meth_empty()
    if absent:
        out["meth_verdict"] = "call_absent"
        return out
    c0 = d[d.call_site]
    if c0.empty:
        return out
    r0 = c0.iloc[0]
    out.update(call_site=r0.site, call_label=r0.label, call_score=r0.score,
               loss_q=r0.q if r0.test == "loss" else np.nan)
    alt = d[(~d.inside) & (d.label == "centromere_like")]
    if len(alt):
        sig = alt[(alt.test == "gain") & (alt.q <= PLANT_Q)]
        a = (sig if len(sig) else alt).sort_values("score", ascending=False).iloc[0]
        out.update(alt_site=a.site, alt_start=int(a.start), alt_end=int(a.end), alt_origins=a.origins,
                   alt_label=a.label, alt_score=a.score, gain_q=a.q if a.test == "gain" else np.nan)
    loss = np.isfinite(out["loss_q"]) and out["loss_q"] <= PLANT_Q
    gainq = np.isfinite(out["gain_q"]) and out["gain_q"] <= PLANT_Q
    if r0.label in NO_DATA:
        out["meth_verdict"] = "no_data"
    elif r0.label == "centromere_like" and not loss:
        out["meth_verdict"] = "shared"
    elif r0.label == "dip" and len(alt):
        out["meth_verdict"] = "shifted" if (loss or gainq) else "possible_shift"
    else:
        out["meth_verdict"] = "unclear"
    return out


def hic_verdict(r: dict | None, call, bin_size: int, ref: dict | None = None, ref_bin: int = 0) -> dict:
    """One non-deciding library x chromosome (a row of its hicfuse table); ref: the deciding library's row."""
    out = dict(hic_verdict="no_data", hic_peak_start=-1, hic_peak_end=-1)
    if r is None or not r.get("callable") or call is None or int(r.get("anch_peak_start", -1)) < 0:
        return out
    pk = (int(r["anch_peak_start"]), int(r["anch_peak_end"]))
    core = (int(r["anch_core_start"]), int(r["anch_core_end"]))
    out.update(hic_peak_start=pk[0], hic_peak_end=pk[1])
    if not r.get("gate_pass"):
        out["hic_verdict"] = "uninformative"
    elif r.get("anch_level", "none") != "none" and gap(*pk, *call) <= HIC_TOL_BINS * bin_size:
        out["hic_verdict"] = "shared"
    elif r.get("anch_strong") and gap(*pk, *call) > HIC_FAR_BINS * bin_size and \
            gap(*core, *call) > HIC_TOL_BINS * bin_size:
        same = ref is not None and int(ref.get("anch_peak_start", -1)) >= 0 and ref.get("anch_level", "none") != "none" \
            and gap(*pk, int(ref["anch_peak_start"]), int(ref["anch_peak_end"])) <= HIC_TOL_BINS * max(bin_size, ref_bin)
        out["hic_verdict"] = "shared" if same else "shifted"
    else:
        out["hic_verdict"] = "unclear"
    return out


def combine(vs: list) -> str:
    """Plant verdict from its lines' verdicts (methylation first)."""
    v = [x for x in vs if x]
    if "call_absent" in v:
        return "call_absent"
    if "shifted" in v:
        return "conflict" if "shared" in v else "shifted"
    if "shared" in v:
        return "shared"
    if "possible_shift" in v:
        return "possible_shift"
    if all(x in ("no_data", "uninformative") for x in v):
        return "no_data"
    return "unclear"


def score(chroms: list, calls: dict, provenance: dict, meth_sites: pd.DataFrame | None, absent: dict,
          libraries: dict, reference: dict | None = None) -> tuple:
    """(plants table, tested site table). calls: {chrom: (start, end) | None} the integrated calls; provenance:
    {(kind, name): dict(verdict, ...)}; absent: {(plant, chrom): True} where the call lies in a missing segment of
    the plant; libraries: {name: hicfuse result} of the libraries that do not decide; reference: the deciding
    library's hicfuse result."""
    st = site_tests(meth_sites, calls) if meth_sites is not None and len(meth_sites) else None
    plants = sorted({n for _, n in provenance if (("meth", n) in provenance and st is not None) or n in libraries})
    refrows = {r["chrom"]: r for r in reference["table"].to_dict("records")} if reference is not None else {}
    ref_bin = int(reference["bin_size"]) if reference is not None else 0
    rows = []
    for p in plants:
        prov = (provenance.get(("meth", p)) or provenance.get(("hic", p)) or {}).get("verdict", "unknown")
        lib = libraries.get(p)
        hrows = {r["chrom"]: r for r in lib["table"].to_dict("records")} if lib is not None else {}
        for c in chroms:
            call = calls.get(c)
            r = dict(plant=p, provenance=prov, chrom=c)
            mv = None
            if st is not None and ("meth", p) in provenance:
                r.update(meth_verdict(st[(st.chrom == c) & (st["sample"] == p)], bool(absent.get((p, c)))))
                mv = r["meth_verdict"]
            else:
                r.update(_meth_empty("."))
            hv = None
            if lib is not None:
                r.update(hic_library=p, **hic_verdict(hrows.get(c), call, int(lib["bin_size"]), refrows.get(c),
                                                      ref_bin))
                hv = r["hic_verdict"]
            else:
                r.update(hic_library=".", hic_verdict=".", hic_peak_start=-1, hic_peak_end=-1)
            if call is None:
                r["verdict"] = "no_call"
            elif absent.get((p, c)):
                r["verdict"] = "call_absent"
            else:
                r["verdict"] = combine([mv, hv])
            shift = "."
            if r["verdict"] in ("shifted", "possible_shift", "conflict") and call is not None:
                if r["alt_start"] >= 0 and mv in ("shifted", "possible_shift"):
                    shift = gap(r["alt_start"], r["alt_end"], *call)
                elif r["hic_peak_start"] >= 0 and hv == "shifted":
                    shift = gap(r["hic_peak_start"], r["hic_peak_end"], *call)
            r["shift_bp"] = shift
            rows.append(r)
    return pd.DataFrame(rows, columns=PLANT_COLS), st


def one_based(df: pd.DataFrame, cols) -> pd.DataFrame:
    """A copy for output with 1-based starts (cenfuse's TSVs are 1-based inclusive; -1 = absent stays -1)."""
    out = df.copy()
    for c in cols:
        if c in out:
            v = pd.to_numeric(out[c], errors="coerce")
            out[c] = np.where(v >= 0, v + 1, v).astype("int64") if v.notna().all() else np.where(v >= 0, v + 1, v)
    return out


def call_counts(tab: pd.DataFrame, chroms: list) -> pd.DataFrame:
    """Per chromosome: plants whose centromere is at the call, shifted, absent, or unresolved."""
    rows = []
    for c in chroms:
        d = tab[tab.chrom == c]
        rows.append(dict(chrom=c, plants_shared=int((d.verdict == "shared").sum()),
                         plants_shifted=int(d.verdict.isin(["shifted"]).sum()),
                         plants_absent=int((d.verdict == "call_absent").sum()),
                         plants_unresolved=int(d.verdict.isin(["possible_shift", "conflict", "unclear",
                                                                "no_data"]).sum())))
    return pd.DataFrame(rows)
