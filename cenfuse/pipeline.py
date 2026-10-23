"""Stage orchestration: satellite -> genotype -> hic -> karyotype -> provisional -> meth -> fuse -> plants -> loci ->
output -> plot.

The satellite engine keeps its own checkpoints (OUT/satellite/work). The genotype, hic, karyotype and meth
stages are checkpointed per genome (OUT/work/<stage>/<name>.pkl + .done), keyed by cenfuse's version and fingerprints
of every input they read (including the calls they are anchored on), so a changed input recomputes them and an
unchanged rerun resumes. The later stages are cheap and always rerun.

Provenance (docs/plants.md): genotype and karyotype tell, per input, whether it comes from the reference plant,
another plant, or cannot be told. Only reference and unknown inputs decide the reference calls; every plant is
scored on its own (plants.py)."""
from __future__ import annotations

import json
import pickle
import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__, fuse, hicfuse, inputs, loci, methfuse, output, plants, provenance, rdna, seqcalls, \
    telomere, teline
from .meth import io as mio
from .meth import params as mp
from .meth import track as mtrack
from .params import ANCHOR_LEVELS, HIC_MAX_THREADS, KARYO_WINDOW
from .sat.genomes import load_gene_gff
from .sat.params import Params
from .sat.util import SatelliteError
from .util import CenfuseError, file_fingerprint, log, obj_fingerprint

STAGES = ["satellite", "rdna", "genotype", "hic", "karyotype", "meth"]
CEN_ROLES = ("centromeric", "cen_associated")
_PKG = Path(__file__).parent


def code_key(*modules, packages=()) -> str:
    """Hash of the source a stage runs: its cenfuse modules, every module of the engine packages it calls, and
    cenfuse's params and util. Any code change recomputes the stage."""
    files = [_PKG / m for m in (*modules, "params.py", "util.py")]
    for pk in packages:
        files += sorted((_PKG / pk).glob("*.py"))
    return obj_fingerprint([f.read_text() for f in files])


def feature_classes(seq: seqcalls.SeqLine, tel: dict, telo: dict, nors=()) -> dict:
    """Feature classes whose 3D self-association the Hi-C line tests: the sequence centromeres (positive control in
    Rabl nuclei), the chromosome ends, telomere arrays, 45S arrays (NORs: the nucleolus), interstitial /
    subtelomeric / rDNA satellite blocks pooled and per family, and each concentrated LTR clade's cores."""
    out = {"centromeres (sequence calls)": dict(kind="centromere", label="sequence calls",
                                                iv=[(r.chrom, r.start, r.end) for r in seq.calls.itertuples(index=False)
                                                    if r.start >= 0]),
           "chromosome ends": dict(kind="ends", label="outermost bins", iv=[])}
    if telo.get("ends"):
        out["telomere arrays"] = dict(kind="telomere", label=telo["unit"],
                                      iv=[(c, s, e) for c, _, s, e, *_ in telo["ends"]])
    if nors:
        out["45S rDNA arrays (NORs)"] = dict(kind="nor", label="18S motif", iv=[(c, s, e) for c, s, e, *_ in nors])
    cl = seqcalls.satellite_clusters(seq.arrays)
    pooled = {"interstitial blocks (all families)": cl.role == "interstitial",
              "subtelomeric blocks (all families)": cl.role == "subtelomeric",
              "rDNA arrays (all families)": (cl.role == "rDNA") | (cl.rdna != ".")}
    for name, sel in pooled.items():
        d = cl[sel]
        if len(d):
            out[name] = dict(kind=name.split()[0], label="pooled", iv=list(zip(d.chrom, d.start, d.end)))
    for fam, d in cl.groupby("family"):
        out[f"satellite {fam}"] = dict(kind=f"satellite:{d.role.iat[0]}", label=fam, iv=list(zip(d.chrom, d.start, d.end)))
    for clade, cores in tel.get("clades", []):
        out[f"LTR {clade} cores"] = dict(kind="te_core", label=clade,
                                         iv=[(c, int(x["start"]), int(x["end"])) for c, x in cores.items()])
    return out


def avoid_intervals(seq: seqcalls.SeqLine, D: float) -> dict:
    """{chrom: [(start, end)]} where reads mismap and genotypes or coverage are not judged: sequence calls +- 2 D and
    satellite clusters."""
    out = {}
    for r in seq.calls.itertuples(index=False):
        if r.start >= 0:
            out.setdefault(r.chrom, []).append((max(0, int(r.start - 2 * D)), int(r.end + 2 * D)))
    for r in seqcalls.satellite_clusters(seq.arrays).itertuples(index=False):
        out.setdefault(r.chrom, []).append((int(r.start), int(r.end)))
    return out


def plant_sites(seq: seqcalls.SeqLine, hic, others: dict) -> dict:
    """{chrom: [(origin, start, end)]} where every plant is scored besides the methylation items: centromeric-family
    satellite clusters, the deciding library's prior-free peak, and the anchored peaks of the other libraries."""
    out = {}
    cl = seqcalls.satellite_clusters(seq.arrays)
    for r in cl[cl.role.isin(CEN_ROLES)].itertuples(index=False):
        out.setdefault(r.chrom, []).append((f"satellite:{r.family}", int(r.start), int(r.end)))
    if hic is not None:
        for r in hic["table"].itertuples(index=False):
            if r.denovo_peak_start >= 0:
                out.setdefault(r.chrom, []).append(("denovo_hic", int(r.denovo_peak_start), int(r.denovo_peak_end)))
    for n, h in others.items():
        for r in h["table"].itertuples(index=False):
            if r.callable and r.anch_peak_start >= 0:
                out.setdefault(r.chrom, []).append((f"hic:{n}", int(r.anch_peak_start), int(r.anch_peak_end)))
    return out


def deciding_library(hics: dict, verdicts: dict):
    """The Hi-C library that decides the reference calls: a reference library, else one of unknown origin. Among
    them, the one whose gate passes on the most chromosomes (its peaks land on the reference's sequence calls of the
    other chromosomes, the reason an unknown library may decide at all), then the most reads. None when every
    library comes from another plant."""
    rank = {"reference": 0, "unknown": 1}
    cand = sorted((rank[verdicts[("hic", n)]["verdict"]], -int(h["table"].gate_pass.sum()),
                   -float(np.sum(h["coverage"])), n) for n, h in hics.items()
                  if verdicts[("hic", n)]["verdict"] in rank)
    return cand[0][-1] if cand else None


class Run:
    def __init__(self, genome_dir, out=None, threads: int = 8, satellite_dir=None, seed_engine: str = "fastan",
                 seed_cache=None, classify_from: str = "", force: str | None = None, sat_params: dict | None = None):
        self.gdir = Path(genome_dir).resolve()
        self.out = Path(out).resolve() if out else self.gdir / "cenfuse"
        self.work = self.out / "work"
        self.threads = max(1, int(threads))
        self.sat_import = Path(satellite_dir).resolve() if satellite_dir else None
        self.seed_engine, self.seed_cache, self.classify_from = seed_engine, seed_cache, classify_from
        self.force = force
        self.sat_params = sat_params or {}           # tests only: satellite-engine Params overrides (scaled-down constants)
        self.genomes = inputs.discover(self.gdir)
        self._prep = {}
        self._nors = {}
        if self.sat_import and not self.sat_import.is_dir():
            raise CenfuseError(f"--satellite-dir {self.sat_import}: not a directory")

    # ------------------------------------------------------------------ checkpoints
    def _cached(self, stage: str, name: str, fp: dict):
        done, pk = self.work / stage / f"{name}.done", self.work / stage / f"{name}.pkl"
        if self.force in (stage, "all") or not (done.exists() and pk.exists()):
            return None
        try:
            if json.loads(done.read_text()) != fp:
                return None
            with open(pk, "rb") as f:
                return pickle.load(f)
        except (OSError, ValueError, pickle.UnpicklingError, EOFError):
            return None

    def _store(self, stage: str, name: str, fp: dict, obj) -> None:
        d = self.work / stage
        d.mkdir(parents=True, exist_ok=True)
        tmp = d / f"{name}.pkl.tmp"
        with open(tmp, "wb") as f:
            pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
        tmp.replace(d / f"{name}.pkl")
        (d / f"{name}.done").write_text(json.dumps(fp))

    def _meth_prep(self, g, sat_dir: Path):
        """The methylation engine's inputs for genome g, built once per run (the class bins are also cached on disk)."""
        if g.name not in self._prep:
            arrays = sat_dir / f"{g.name}_satellite_arrays.gff"
            self._prep[g.name] = methfuse.prepare(g, arrays, self.work / "meth" / "meth_cache", self.threads)
        return self._prep[g.name]

    def _meth_inputs_fp(self, g, sat_dir: Path) -> dict:
        arrays = sat_dir / f"{g.name}_satellite_arrays.gff"
        samples = pd.read_csv(g.meth, sep="\t", dtype=str).fillna("")
        files = sorted(str(Path(v) if Path(v).is_absolute() else g.meth.parent / v)
                       for col in ("CHG", "CG", "CHH") if col in samples for v in samples[col] if v)
        return dict(samples=file_fingerprint(g.meth), files=[file_fingerprint(f) for f in files if Path(f).exists()],
                    ann=[file_fingerprint(x) for x in (g.tes, arrays, g.genes, g.exclude) if x])

    # ------------------------------------------------------------------ stages
    def stage_satellite(self) -> Path:
        if self.sat_import:
            log.info(f"satellite: reading earlier satellite outputs from {self.sat_import} (--satellite-dir)")
            return self.sat_import
        from .sat.pipeline import STAGES as SAT_STAGES, Pipeline
        view = inputs.make_view(self.genomes, self.gdir, self.work / "genomes" / self.gdir.name)
        force = self.force if self.force in SAT_STAGES else ("seeds" if self.force in ("satellite", "all") else None)
        try:
            Pipeline(view, self.threads, self.seed_cache, Params(seed_engine=self.seed_engine,
                                                                 classify_from=self.classify_from, **self.sat_params),
                     force=force, out_dir=self.out / "satellite").run()
        except SatelliteError as e:
            raise CenfuseError(f"satellite engine: {e}") from None
        return self.out / "satellite"

    def stage_rdna(self, g) -> list:
        """45S arrays (NORs) from the universal 18S motif (rdna.py), checkpointed by the assembly index and the code."""
        fp = dict(version=__version__, code=code_key("rdna.py"), fai=file_fingerprint(f"{g.fasta}.fai"),
                  chroms=obj_fingerprint(g.chroms))
        hit = self._cached("rdna", g.name, fp)
        if hit is None:
            hit = rdna.detect(g.fasta, g.chroms)
            self._store("rdna", g.name, fp, hit)
        self._nors[g.name] = hit
        return hit

    def stage_genotype(self, g, samples: list, windows: list) -> dict:
        """Fixed differences of every methylation sample and the private-difference tests (provenance)."""
        files = sorted(str(f) for s in samples for f in s.files.values())
        fp = dict(version=__version__, code=code_key("provenance.py"), fai=file_fingerprint(f"{g.fasta}.fai"),
                  files=[file_fingerprint(f) for f in files], names=[s.name for s in samples],
                  windows=obj_fingerprint(windows))
        hit = self._cached("genotype", g.name, fp)
        if hit is not None:
            log.info(f"genotype: {g.name} checkpoint found")
            return hit
        t0 = time.time()
        warn = provenance.check_assembly(samples, g.fasta, windows)          # stops the run on another assembly
        gts = provenance.genotypes(samples, g.fasta, g.chroms, windows, self.threads)
        tests = provenance.private_tests(gts)
        res = dict(genotypes=gts, tests=tests, verdicts=provenance.genotype_verdicts(gts, tests), warnings=warn,
                   windows=windows)
        self._store("genotype", g.name, fp, res)
        log.info(f"genotype: {g.name}: {len(windows)} windows, sites tested "
                 f"{', '.join(f'{x.sample} {len(x.tested):,}' for x in gts)} ({time.time() - t0:.0f} s)")
        return res

    def stage_hic(self, g, lib, seq: seqcalls.SeqLine, sat_dir: Path, classes: dict | None):
        name = g.name if lib.name == "hic" else f"{g.name}.{lib.name}"
        calls = {r.chrom: (int(r.start), int(r.end), r.confidence) for r in seq.calls.itertuples(index=False)}
        gff = self.work / "hic" / f"{g.name}.sequence_calls.gff3"
        gff.parent.mkdir(parents=True, exist_ok=True)
        with open(gff, "w") as f:                     # Hi-C engine --satellites: centromere features for its sat_mid
            f.write("##gff-version 3\n")
            for c, (s, e, _) in calls.items():
                if s >= 0:
                    f.write(f"{c}\tcenfuse\tcentromere\t{s + 1}\t{e}\t.\t.\t.\tID=cen_{c}\n")
        arrays = sat_dir / f"{g.name}_satellite_arrays.gff"
        fai = self.work / "genomes" / f"{g.name}.chromosomes.fai"
        fai.parent.mkdir(parents=True, exist_ok=True)
        txt = "".join(f"{c}\t{L}\t0\t0\t0\n" for c, L in g.chroms.items())
        if not fai.exists() or fai.read_text() != txt:     # the Hi-C engine's matrix cache is keyed by the file's mtime
            fai.write_text(txt)
        fp = dict(version=__version__, code=code_key("hicfuse.py", packages=("hic",)),
                  pairs=file_fingerprint(lib.path), chroms=obj_fingerprint(g.chroms),
                  calls=obj_fingerprint(calls), arrays=file_fingerprint(arrays), classes=obj_fingerprint(classes),
                  nors=obj_fingerprint(self._nors.get(g.name, [])))
        hit = self._cached("hic", name, fp)
        if hit is not None:
            log.info(f"hic: {name} checkpoint found")
            return hit
        t0 = time.time()
        (self.out / "hic").mkdir(parents=True, exist_ok=True)
        res = hicfuse.run(lib.path, fai, self.out / "hic" / name, [gff, arrays], calls,
                          min(self.threads, HIC_MAX_THREADS), classes, nors=self._nors.get(g.name, []))
        self._store("hic", name, fp, res)
        log.info(f"hic: {name} done ({time.time() - t0:.0f} s)")
        return res

    def stage_karyotype(self, g, sat_dir: Path, w: int, seq: seqcalls.SeqLine) -> dict:
        """Windows not judged for missing segments (the sequence calls, gaps, satellite-dominated windows, ends), and
        each methylation sample's coverage before the depth filter."""
        calls = {}
        for r in seq.calls.itertuples(index=False):
            if r.start >= 0:
                calls.setdefault(r.chrom, []).append((int(r.start), int(r.end)))
        fp = dict(version=__version__, code=code_key("provenance.py", "methfuse.py", packages=("meth",)),
                  fai=file_fingerprint(f"{g.fasta}.fai"), w=w, calls=obj_fingerprint(calls),
                  arrays=obj_fingerprint(seq.arrays[["chrom", "start", "end"]].values.tolist()),
                  meth=self._meth_inputs_fp(g, sat_dir) if g.meth else None)
        hit = self._cached("karyotype", g.name, fp)
        if hit is not None:
            log.info(f"karyotype: {g.name} checkpoint found")
            return hit
        t0 = time.time()
        cov = {}
        if g.meth:
            try:
                prep = self._meth_prep(g, sat_dir)
                cov = prep.coverage if prep is not None else {}
            except CenfuseError as e:
                log.warning(f"{g.name}: methylation coverage unavailable: {e}")
        res = dict(skip=provenance.skip_windows(g.chroms, w, calls, g.fasta, seq.arrays), meth_cov=cov, w=w)
        self._store("karyotype", g.name, fp, res)
        log.info(f"karyotype: {g.name} done ({time.time() - t0:.0f} s)")
        return res

    def stage_meth(self, g, sat_dir: Path, anchors: dict, items: dict, width_min: int, deciders: list,
                   sites: dict):
        fp = dict(version=__version__, code=code_key("methfuse.py", packages=("meth",)),
                  assembly=file_fingerprint(f"{g.fasta}.fai"), width_min=width_min, chroms=obj_fingerprint(g.chroms),
                  anchors=obj_fingerprint(anchors), items=obj_fingerprint(items), deciders=list(deciders),
                  sites=obj_fingerprint(sites), **self._meth_inputs_fp(g, sat_dir))
        hit = self._cached("meth", g.name, fp)
        if hit is not None:
            log.info(f"meth: {g.name} checkpoint found")
            return hit
        t0 = time.time()
        prep = self._meth_prep(g, sat_dir)
        if prep is None:
            return None
        (self.out / "methylation").mkdir(parents=True, exist_ok=True)
        res = methfuse.run(prep, anchors, items, self.out / "methylation" / g.name, width_min, decide=deciders,
                           sites=sites)
        self._store("meth", g.name, fp, res)
        log.info(f"meth: {g.name} done ({time.time() - t0:.0f} s)")
        return res

    # ------------------------------------------------------------------ assembly
    @staticmethod
    def evidences(g, seq: seqcalls.SeqLine, tel: dict, hic, meth_items=None) -> list:
        marker = tel.get("marker")
        hrows = {} if hic is None else {r["chrom"]: r for r in hic["table"].to_dict("records")}
        cen = seq.arrays[seq.arrays.role.isin(CEN_ROLES)]
        by = {c: d for c, d in cen.groupby("chrom")}

        def sat_bp(c):
            d = by.get(c)

            def f(s, e):
                if d is None:
                    return 0
                return int((np.minimum(d.end, e) - np.maximum(d.start, s)).clip(lower=0).sum())
            return f
        out = []
        for r in seq.calls.itertuples(index=False):
            te = None
            if marker and r.chrom in marker[1]:
                core = marker[1][r.chrom]
                te = dict(clade=marker[0], start=int(core["start"]), end=int(core["end"]), decisive=bool(marker[2]))
            elif r.te_start >= 0:                              # imported satellite calls made with tes/
                te = dict(clade=r.te_clade, start=int(r.te_start), end=int(r.te_end),
                          decisive=r.te_clade.lower() == "crm")
            h = hrows.get(r.chrom)
            if h is not None and not h["callable"]:
                h = None
            meth = None if meth_items is None else meth_items.get(r.chrom, {})
            out.append(fuse.Evidence(r.chrom, int(r.chrom_len), dict(start=int(r.start), end=int(r.end),
                                                                    evidence=r.evidence, confidence=r.confidence,
                                                                    flags=r.flags, main_family=r.main_family,
                                                                    alt_start=int(r.alt_start),
                                                                    alt_end=int(r.alt_end)),
                                     te, h, meth, sat_bp(r.chrom)))
        return out

    def provenance_of(self, g, seq, sat_dir: Path, tel: dict, telo: dict):
        """Provenance of every input, then the Hi-C libraries: (samples, geno, hics, kar, profiles, segments,
        verdicts)."""
        samples = []
        if g.meth:
            try:
                samples = mio.read_samples(g.meth)
            except mio.InputError as e:
                log.warning(f"{g.name}: methylation line skipped: {e}")
                g.notes.append(f"methylation skipped: {e}")
                g.meth = None
        D = mtrack.chromosome_set(mio.read_fai(f"{g.fasta}.fai"))[2]
        avoid = avoid_intervals(seq, D)
        geno = {}
        if samples:
            windows = provenance.genotype_windows(g.chroms, avoid, g.fasta)
            geno = self.stage_genotype(g, samples, windows)
            for x in geno["warnings"]:
                log.warning(f"{g.name}: {x}")
        known = {n: v["verdict"] for n, v in geno.get("verdicts", {}).items()}
        classes = feature_classes(seq, tel, telo, self._nors.get(g.name, []))
        hics = {}
        for lib in g.hic:
            try:
                hics[lib.name] = self.stage_hic(g, lib, seq, sat_dir, None if known.get(lib.name) == "other"
                                                else classes)
            except CenfuseError as e:
                log.warning(f"{g.name}: Hi-C library {lib.name} skipped: {e}")
                g.notes.append(f"Hi-C {lib.name} skipped: {e}")
        w = max([KARYO_WINDOW] + [int(h["bin_size"]) for h in hics.values()])
        if len(samples) >= 2 or len(hics) >= 2:            # a missing segment needs another input of the same type
            kar = self.stage_karyotype(g, sat_dir, w, seq)
        else:
            kar = dict(skip={}, meth_cov={}, w=w)
        profiles = {}
        for n, cov in kar["meth_cov"].items():
            profiles[("meth", n)] = {c: provenance.window_profile(cov[c], mp.BIN, w, L)
                                     for c, L in g.chroms.items() if c in cov}
        for n, h in hics.items():
            b = h["bins"]
            idx = {x: i for i, x in enumerate(b.names)}
            profiles[("hic", n)] = {c: provenance.window_profile(h["coverage"][b.sl(idx[c])], b.size, w, L)
                                    for c, L in g.chroms.items() if c in idx}
        segments = provenance.missing_segments(profiles, kar["skip"], w) if len(profiles) >= 2 else {}
        verdicts = provenance.combine([s.name for s in samples], list(hics), geno.get("verdicts", {}), segments)
        provenance.log_summary(g.name, verdicts)
        return samples, geno, hics, kar, profiles, segments, verdicts

    # ------------------------------------------------------------------ driver
    def run(self) -> None:
        t0 = time.time()
        self.out.mkdir(parents=True, exist_ok=True)
        log.info(f"cenfuse {__version__}: {self.gdir.name}: " + "; ".join(
            f"{g.name} ({len(g.chroms)} chromosomes; {', '.join(g.lines)})" for g in self.genomes))
        for g in self.genomes:
            for n in g.notes:
                log.info(f"{g.name}: {n}")
        sat_dir = self.stage_satellite()
        summary, figs = [], []
        for g in self.genomes:
            seq = seqcalls.read(sat_dir, g.name, g.chroms)
            tel = teline.te_line(g.tes, g.chroms)
            telo = telomere.detect(g.fasta, g.chroms)
            nors = self.stage_rdna(g)
            samples, geno, hics, kar, profiles, segments, verdicts = self.provenance_of(g, seq, sat_dir, tel, telo)
            deciders = [s.name for s in samples if verdicts[("meth", s.name)]["verdict"] != "other"]
            dec_lib = deciding_library(hics, verdicts)
            hic = hics.get(dec_lib) if dec_lib else None
            others = {n: h for n, h in hics.items() if n != dec_lib}
            if hics and hic is None:
                log.warning(f"{g.name}: every Hi-C library comes from another plant; Hi-C does not decide the "
                            f"reference calls (plants are still scored)")
            if samples and not deciders:
                log.info(f"{g.name}: every methylation sample comes from another plant; methylation scores the "
                         f"plants and does not decide the reference calls")
            if hic is not None:
                hic["conflict"] = c = fuse.sequence_conflict(hic["table"])
                if c["conflict"]:
                    log.warning(f"{g.name}: Hi-C's prior-free calls are confident on {c['denovo_confident']}/{c['n']} "
                                f"chromosomes but agree with the sequence calls on only {c['denovo_on_seq']}, and the "
                                f"gate fails: the two lines disagree genome-wide, so at least one is wrong. The "
                                f"sequence calls are kept; give tes/ (CRM cores) or methylation to decide "
                                f"(cenfuse_summary.tsv: hic_sequence_conflict)")
            evs = self.evidences(g, seq, tel, hic)
            prov = [fuse.decide(e) for e in evs]
            meth, items = None, {}
            if g.meth:
                anchors = {e.chrom: d["call"] for e, d in zip(evs, prov) if d["call"] and d["confidence"] in
                           ANCHOR_LEVELS}
                for e, d in zip(evs, prov):
                    it = [("sequence_call", e.seq["start"], e.seq["end"])]
                    peak = fuse.hic_peak(e)
                    if peak:
                        it.append(("hic_peak", *peak))
                    if e.te:
                        it.append(("te_core", e.te["start"], e.te["end"]))
                    if d["call"]:
                        it.append(("provisional", *d["call"]))
                    items[e.chrom] = [x for x in it if x[1] >= 0]
                try:
                    meth = self.stage_meth(g, sat_dir, anchors, items, hic["bin_size"] if hic else 0, deciders,
                                           plant_sites(seq, hic, others))
                except CenfuseError as e:
                    log.warning(f"{g.name}: methylation line skipped: {e}")
                    g.notes.append(f"methylation skipped: {e}")
            mi = None
            if meth is not None:
                pooled = meth["support"][meth["support"]["sample"] == "pooled"]
                mi = {c: {r.item: dict(label=r.label, score=r.score) for r in d.itertuples(index=False)}
                      for c, d in pooled.groupby("chrom")}
                evs = self.evidences(g, seq, tel, hic, mi)
            final = [fuse.decide(e) for e in evs]
            t = fuse.table(evs, final)
            calls = {e.chrom: d["call"] for e, d in zip(evs, final)}
            absent = provenance.call_absent(profiles, kar["skip"], segments, calls, kar["w"])
            ptab, psites = plants.score(list(g.chroms), calls, verdicts, meth["sites"] if meth is not None else None,
                                        absent, others, hic)
            t["meth_samples"] = ",".join(deciders) if meth is not None and deciders else "."
            t["hic_library"] = dec_lib or "."
            t["hic_provenance"] = verdicts[("hic", dec_lib)]["verdict"] if dec_lib else "."
            if len(ptab):
                t = t.merge(plants.call_counts(ptab, list(t.chrom)), on="chrom", how="left")
            else:
                for col in ("plants_shared", "plants_shifted", "plants_absent", "plants_unresolved"):
                    t[col] = "."
            output.write_calls(t, g.name, g.chroms, self.out)
            ends = {n: int(np.sum(h["coverage"])) for n, h in hics.items()}
            ptable = provenance.table(verdicts, geno, segments, deciders if meth is not None else [], dec_lib, ends)
            output.write_tsv(ptable, self.out / f"{g.name}.provenance.tsv")
            if len(ptab):
                output.write_tsv(plants.one_based(ptab, ("alt_start", "hic_peak_start")),
                                 self.out / f"{g.name}.plants.tsv")
            if psites is not None:
                output.write_tsv(plants.one_based(psites, ("start", "mid")), self.out / f"{g.name}.plant_sites.tsv")
            fused = {e.chrom: dict(call=d["call"], relation=d["relation"], seq_evidence=e.seq["evidence"],
                                   seq=(e.seq["start"], e.seq["end"]) if e.seq["start"] >= 0 else None)
                     for e, d in zip(evs, final)}
            hic_l = None
            if hic is not None:
                hic_l = dict(rows={r["chrom"]: r for r in hic["table"].to_dict("records")}, bins=hic["bins"],
                             index={n: i for i, n in enumerate(hic["bins"].names)}, hub=hic["hub"], S=hic["S"],
                             denovo_S=hic["denovo_S"], telomere_S=hic["telomere_S"])
                output.write_bedgraph(hic["S"], hic["bins"], hic["members"],
                                      self.out / f"{g.name}.hic_anchored.bedgraph",
                                      f"cenfuse {g.name} anchored Hi-C score (x chromosome median)")
            genes = load_gene_gff(g.genes, g.chroms)[0] if g.genes else None
            lt, feats = loci.table(g.name, g.chroms, fused, seq, seqcalls.satellite_clusters(seq.arrays), tel, hic_l,
                                   meth, genes, g.fasta, telo, nors)
            output.write_loci(lt, feats, g.name, g.chroms, self.out)
            fc = hic["classes"] if hic is not None else None
            if fc is not None:
                output.write_tsv(fc, self.out / f"{g.name}.feature_contacts.tsv")
            output.write_tsv(loci.feature_signatures(lt, fc), self.out / f"{g.name}.feature_signatures.tsv")
            summary.append(output.summary_row(g.name, g, t, hic, meth, tel, telo, verdicts, ptab, nors))
            figs.append(dict(genome=g, seq=seq, tel=tel, hic=hic, meth=meth, calls=t, loci=lt, telo=telo,
                             plants=ptab, psites=psites, verdicts=verdicts, libraries=others))
            log.info(f"{g.name}: {len(t[t.start != '.'])}/{len(t)} chromosomes called; relations "
                     f"{t.relation.value_counts().to_dict()}; confidence "
                     f"{t[t.start != '.'].confidence.value_counts().to_dict()}")
            if len(ptab):
                sh = ptab[ptab.verdict.isin(["shifted", "call_absent"])]
                log.info(f"{g.name}: plants: {ptab.plant.nunique()} scored; "
                         + ("; ".join(f"{r.plant} {r.chrom} {r.verdict}" for r in sh.itertuples()) or "no shift"))
        output.write_tsv(pd.DataFrame(summary), self.out / "cenfuse_summary.tsv")
        try:
            from . import plot
            plot.overview(figs, self.out)
        except Exception as e:                       # a figure never costs the results
            log.warning(f"overview figure failed ({type(e).__name__}: {e}); results are complete")
        log.info(f"done ({(time.time() - t0) / 60:.1f} min); outputs in {self.out}")
