"""Stage orchestration: satellite -> telomere -> rdna -> rdna5s -> gaps -> attrough -> genotype -> hic -> karyotype -> fold
-> provisional -> meth -> fuse -> cpg -> plants -> morphology -> interpret -> output -> loci -> report -> plot.

The satellite engine keeps its own checkpoints (OUT/satellite/work), and next to them the chromosome set they were
made on (chromosome_set.json): a resume on another set stops (Run._check_chromosome_set). The telomere, rdna, rdna5s,
gaps, attrough, genotype, hic, fold, karyotype, meth and cpg stages are checkpointed per genome (OUT/work/<stage>/<name>.pkl +
.done), keyed by the code they run and fingerprints of every input they read (including the calls they are anchored on; most also
by cenfuse's version), so a changed input recomputes them and an unchanged rerun resumes. The later stages are cheap and
always rerun.

Provenance (docs/plants.md): genotype and karyotype tell, per input, whether it comes from the reference plant,
another plant, or cannot be told. Only reference and unknown inputs decide the reference calls; every plant is
scored on its own (plants.py)."""
from __future__ import annotations

import ast
import functools
import json
import os
import pickle
import socket
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from . import __version__, ageline, atline, cdr, fold, fuse, hicfuse, inputs, interpret, linetests, loci, manifest, \
    methfuse, morphology, output, plants, provenance, rdna, seqcalls, telomere, teline
from .hic.balance import balanced, ice
from .hic.bins import Bins, coarsen
from .hic.io import is_cool, read_cool
from .hic.score import block_means
from .meth import io as mio
from .meth import params as mp
from .meth import track as mtrack
from .params import ANCHOR_LEVELS, GATE_MIN_REF, HIC_FAR_BINS, HIC_MAX_THREADS, HIC_TOL_BINS, KARYO_WINDOW
from .sat import signature
from .sat.karyotype import MAX_SEQUENCES
from .sat.genomes import load_gene_gff
from .sat.params import Params, parse_borrow
from .sat.util import SatelliteError
from .runtime import Heartbeat, Usage, atomic_open, atomic_write_text, fasta_md5, file_fingerprint, release_lock, try_lock
from .util import CenfuseError, log, obj_fingerprint

STAGES = ["satellite", "telomere", "rdna", "rdna5s", "gaps", "attrough", "genotype", "hic", "fold", "karyotype", "meth", "cpg"]
CEN_ROLES = ("centromeric", "cen_associated")
CHROM_RECORD = "chromosome_set.json"      # in satellite/work: {genome: sorted names} the checkpoints stand on
LOCK = ".cenfuse.lock"                    # in OUT: 'host pid start command' of the run using the folder (runtime.try_lock)
_PKG = Path(__file__).parent


def _module_file(name: str) -> Path | None:
    """The file of cenfuse module 'a.b' (a/b.py, or a/b/__init__.py for a package), or None outside the package."""
    base = _PKG.joinpath(*name.split(".")) if name else _PKG
    for f in (base.with_suffix(".py"), base / "__init__.py"):
        if f.is_file():
            return f
    return None


def _imports(f: Path) -> set:
    """The cenfuse modules a file imports, at any depth of its code (lazy imports in functions included), resolved to files:
    'from . import x', 'from .sat import te', 'from ..util import gap', 'import cenfuse.x'."""
    pkg = ".".join(f.relative_to(_PKG).with_suffix("").parts[:-1]) if f.name != "__init__.py" else \
        ".".join(f.relative_to(_PKG).parent.parts)
    out = set()
    for node in ast.walk(ast.parse(f.read_text())):
        names = []
        if isinstance(node, ast.ImportFrom):
            if node.level:
                parts = pkg.split(".") if pkg else []
                parts = parts[:len(parts) - (node.level - 1)] if node.level > 1 else parts
                base = ".".join(parts + ([node.module] if node.module else []))
            elif node.module and (node.module == "cenfuse" or node.module.startswith("cenfuse.")):
                base = node.module[len("cenfuse."):] if node.module != "cenfuse" else ""
            else:
                continue
            names = [base] + [f"{base}.{a.name}" if base else a.name for a in node.names]
        elif isinstance(node, ast.Import):
            names = [a.name[len("cenfuse."):] for a in node.names if a.name.startswith("cenfuse.")]
        for n in names:
            m = _module_file(n)
            if m is not None:
                out.add(m)
    return out


def _normalized(f: Path) -> str:
    """What a module does, without comments, docstrings or layout: the dump of its syntax tree, docstrings removed. Editing
    a comment or a docstring recomputes nothing."""
    tree = ast.parse(f.read_text())
    for node in ast.walk(tree):
        body = getattr(node, "body", None)
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)) and body and \
                isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) and isinstance(body[0].value.value, str):
            node.body = body[1:] or [ast.Pass()]
    return ast.dump(tree, include_attributes=False)


@functools.lru_cache(maxsize=None)
def _stage_method(stage: str) -> str:
    """The syntax tree of Run.stage_<stage> in this module, docstring stripped: the stage's own orchestration (what it reads,
    in which order it calls the engines) is code that can change its result."""
    tree = ast.parse((_PKG / "pipeline.py").read_text())
    run = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == "Run")
    fn = next(n for n in run.body if isinstance(n, ast.FunctionDef) and n.name == f"stage_{stage}")
    if fn.body and isinstance(fn.body[0], ast.Expr) and isinstance(fn.body[0].value, ast.Constant):
        fn.body = fn.body[1:] or [ast.Pass()]
    return ast.dump(fn, include_attributes=False)


# modules that change no result, kept out of every code key: how files are fingerprinted, written and locked, and what a run
# records (runtime.py: an atomic write or a lock rule edited there recomputes nothing; manifest.py: run_info.json)
KEYLESS = {"runtime.py", "manifest.py"}


@functools.lru_cache(maxsize=None)
def code_key(*modules, packages=(), stage: str | None = None, data=()) -> str:
    """Key of the code a stage runs: its cenfuse modules, every module of the engine packages it calls, cenfuse's params and
    util, and everything those import inside cenfuse (the import closure), each as its syntax tree without comments or
    docstrings; the stage's own method of Run (stage); and the content of the package data files it reads (data, relative to
    the package). A change to what that code does recomputes the stage; a comment or a docstring does not."""
    todo = [_PKG / m for m in (*modules, "params.py", "util.py")]
    for pk in packages:
        todo += sorted((_PKG / pk).glob("*.py"))
    seen = set()
    while todo:
        f = todo.pop()
        if f in seen or f.relative_to(_PKG).as_posix() in KEYLESS:
            continue
        seen.add(f)
        todo += [x for x in _imports(f) if x not in seen]
    return obj_fingerprint([[str(f.relative_to(_PKG)), _normalized(f)] for f in sorted(seen)] +
                           [["stage", _stage_method(stage) if stage else None],
                            ["data", [file_fingerprint(_PKG / d) for d in data]]])


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


def hic_rows(table: pd.DataFrame, nors=None) -> dict:
    """{chrom: hicfuse row} with gate_genome_untested (no chromosome of the genome could test Hi-C: fewer than
    GATE_MIN_REF references; fuse's Hi-C alone) and the NOR test of the prior-free core (denovo_nor_flag)."""
    untested = not bool((table.gate_n >= GATE_MIN_REF).any()) if len(table) else True
    return {r["chrom"]: denovo_nor_flag(dict(r, gate_genome_untested=untested), nors)
            for r in table.to_dict("records")}


def denovo_nor_flag(h: dict, nors) -> dict:
    """h with 'rDNA' added to its prior-free flags when the prior-free core lies within HIC_TOL_BINS mappable bins of
    a 45S array (rdna.detect): NORs meet in the nucleolus, so a prior-free peak beside one is no clean centromere
    candidate (the anchored peaks get the same test in hicfuse.nor_peaks; the engine's own flags only know
    rDNA-tagged satellite arrays)."""
    if not nors or h.get("denovo_core_start", -1) < 0:
        return h
    s, e, b = int(h["denovo_core_start"]), int(h["denovo_core_end"]), int(h["bin_size"])
    near = any(c == h["chrom"] and hicfuse.visible_gap(s, e, int(a), int(z), h.get("masked") or (),
                                                       HIC_FAR_BINS * b) <= HIC_TOL_BINS * b
               for c, a, z, *_ in nors)
    fl = [f for f in str(h.get("denovo_flags") or "").split(",") if f]
    if near and "rDNA" not in fl:
        h = dict(h, denovo_flags=",".join(fl + ["rDNA"]))
    return h


def cis_entry(cis: dict | None, chrom: str) -> dict | None:
    """The cis line's entry for one chromosome as fuse reads it (Evidence.hic['cis']), or None where it does not act: the extreme,
    the genome status, what tested it ('ltr': hicfuse.cis_line; 'trans': hicfuse.cis_trans_line; a 0.12 result names none and is
    the LTR one), whether it was read outside the terminal eighths (interior: fuse leaves a call there alone) and whether this
    chromosome's own trans peak lies at its extreme (pair)."""
    if not hicfuse.cis_acts(cis) or chrom not in cis.get("peaks", {}):
        return None
    return dict(peak=cis["peaks"][chrom], status=cis["status"], tested_by=cis.get("tested_by", "ltr"),
                interior=bool(cis.get("interior")), pair=bool((cis.get("pairs") or {}).get(chrom)))


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
        try:
            legacy = frozenset(parse_borrow(classify_from))      # --classify-from targets keep the 0.8 chromosome rule
        except SatelliteError as e:
            raise CenfuseError(str(e)) from None
        self.genomes = inputs.discover(self.gdir, legacy=legacy)
        self._prep = {}
        self._nors = {}
        self._chg_ok = {}                            # genome -> the CHG methylation line can run (tes/, CHG files)
        if self.sat_import and not self.sat_import.is_dir():
            raise CenfuseError(f"--satellite-dir {self.sat_import}: not a directory")

    # ------------------------------------------------------------------ checkpoints
    @staticmethod
    def _assembly(g) -> str:
        """Key of the assembly a stage reads: its index and the whole FASTA's md5 (runtime.fasta_md5). A same-length edit
        (hard-masking, polishing) keeps the .fai byte for byte, and must still recompute what reads the sequence."""
        return obj_fingerprint([file_fingerprint(f"{g.fasta}.fai"), fasta_md5(g.fasta)])

    def _cached(self, stage: str, name: str, fp: dict):
        """The checkpoint of stage for name, or None when it must be computed: forced (--force STAGE or all), absent, made
        by other code or on other inputs (its key differs), or unreadable (a torn or foreign pickle: said, then recomputed).
        A stage that reads another's result has that result in its key, so it reruns when the result changed."""
        done, pk = self.work / stage / f"{name}.done", self.work / stage / f"{name}.pkl"
        if self.force in (stage, "all") or not (done.exists() and pk.exists()):
            return None
        try:
            if json.loads(done.read_text()) != fp:
                log.debug(f"{stage}: {name}: checkpoint made by other code or on other inputs; recomputed")
                return None
            with open(pk, "rb") as f:
                return pickle.load(f)
        except Exception as e:                   # OSError, a torn pickle, a library upgrade (ModuleNotFoundError, ...)
            log.warning(f"{stage}: {name}: checkpoint {pk.name} unreadable ({type(e).__name__}: {e}); recomputed")
            return None

    def _store(self, stage: str, name: str, fp: dict, obj) -> None:
        d = self.work / stage
        d.mkdir(parents=True, exist_ok=True)
        tmp = d / f"{name}.pkl.tmp"
        with open(tmp, "wb") as f:
            pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
        tmp.replace(d / f"{name}.pkl")
        atomic_write_text(d / f"{name}.done", json.dumps(fp))

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
                    ann=[fasta_md5(x) for x in (g.tes, arrays, g.genes, g.exclude) if x])      # whole files: annotations are edited

    # ------------------------------------------------------------------ satellite chromosome set
    def _check_chromosome_set(self, sat_work: Path) -> None:
        """Stop when the satellite checkpoints in sat_work were made on another chromosome set than this run analyzes.
        Why: the engine's resume check covers its code and parameters, not the chromosomes, and 0.9 changes the set
        (gradual scaffold tails, <genome>.chromosomes.txt): a resume would reuse seeds, families and calls made on other
        sequences. The earlier set is the record's (CHROM_RECORD) or, in a folder made before 0.9 (no record), the 0.8
        rule's on the .fai. No seeds.done: nothing is checkpointed. A force that redoes the seeds, the first stage the
        set shapes, allows the change."""
        if not (sat_work / "seeds.done").exists() or self.force in ("satellite", "all", "seeds"):
            return
        record, now = sat_work / CHROM_RECORD, {g.name: set(g.chroms) for g in self.genomes}
        if record.exists():
            try:
                then = {name: set(names) for name, names in json.loads(record.read_text()).items()}
            except (OSError, ValueError, AttributeError, TypeError):
                raise CenfuseError(f"{record}: cannot be read, so the chromosomes the satellite checkpoints were made "
                                   "on are unknown; rerun with --force satellite") from None
            how = ""
        else:
            then = {g.name: set(inputs.chromosome_set(g.fai, legacy=True)) for g in self.genomes}
            how = " (cenfuse 0.8's rule)"
        source, changed = {g.name: g.chrom_source for g in self.genomes}, []
        for name in [*then, *(n for n in now if n not in then)]:           # the record's genomes, then the new ones
            old, new = then.get(name), now.get(name)
            if old == new:
                continue
            if old is None:
                changed.append(f"{name}: not among the genomes of the checkpoints, {len(new)} chromosomes now")
            elif new is None:
                changed.append(f"{name}: {len(old)} chromosomes when they were made, no longer in {self.gdir}")
            else:
                changed.append(f"{name}: {len(old)} chromosomes when they were made{how}, {len(new)} now "
                               f"({source[name]}): {len(old - new)} dropped, {len(new - old)} added")
        if changed:
            raise CenfuseError(f"{sat_work}: the satellite checkpoints (seeds, families, calls) were made on another "
                               f"chromosome set: {'; '.join(changed)}; rerun with --force satellite")
        log.debug(f"satellite: the chromosome set is the one the checkpoints in {sat_work} were made on")

    def _record_chromosome_set(self, sat_work: Path) -> None:
        """Write {genome: sorted chromosome names} beside the satellite checkpoints (CHROM_RECORD), the set
        _check_chromosome_set compares on the next run. Called after the engine is constructed and before it runs:
        what its stages checkpoint stands on this set whether the run ends or not, while the earlier record would
        describe checkpoints already wiped (a failed forced run, then a return to the earlier set, would resume on
        them). Written to a temporary file and moved over the record, so it is never torn."""
        sat_work.mkdir(parents=True, exist_ok=True)
        path = sat_work / CHROM_RECORD
        tmp = path.with_name(f"{path.name}.tmp")
        tmp.write_text(json.dumps({g.name: sorted(g.chroms) for g in self.genomes}, indent=1))
        tmp.replace(path)

    def _warn_imported_set(self, sat_dir: Path) -> None:
        """Warn when an imported satellite summary (--satellite-dir) has rows for other sequences than the chromosome
        set analyzed now. Why: its families, prevalence tests and calls were made on those sequences (cenfuse 0.8's
        set, or a list), and an import is taken as it is. A summary that is absent or not one is left to seqcalls.read,
        which stops the run with the reason."""
        for g in self.genomes:
            path = sat_dir / f"{g.name}_centromere_summary.tsv"
            try:
                rows = set(pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, usecols=["chrom"]).chrom)
            except (OSError, ValueError):
                continue
            outside, without = rows - set(g.chroms), set(g.chroms) - rows
            if outside or without:
                log.warning(f"{g.name}: the imported {path.name} has rows for {len(rows)} sequences, but cenfuse "
                            f"analyzes {len(g.chroms)} chromosomes ({len(outside)} rows on sequences outside the set, "
                            f"{len(without)} chromosomes without a row): the satellite calls were made on another "
                            f"chromosome set; rerun without --satellite-dir to make them on this one")

    # ------------------------------------------------------------------ stages
    def stage_satellite(self) -> Path:
        if self.sat_import:
            log.info(f"satellite: reading earlier satellite outputs from {self.sat_import} (--satellite-dir)")
            imported = manifest.read(self.sat_import) or {}
            made = imported.get("cenfuse")
            if imported and imported.get("status") != "done":
                log.warning(f"satellite: the run that made {self.sat_import} did not finish ({manifest.NAME}: status "
                            f"{imported.get('status')}): its outputs may be incomplete")
            if made and made != __version__:
                log.warning(f"satellite: the outputs in {self.sat_import} were made by cenfuse {made} ({manifest.NAME}); "
                            f"this is {__version__}: they keep {made}'s satellite rules")
            elif not made:
                log.info(f"satellite: no {manifest.NAME} beside the imported outputs: the cenfuse version that made them "
                         f"is unknown, and they keep its satellite rules")
            self._warn_imported_set(self.sat_import)
            return self.sat_import
        from .sat.pipeline import STAGES as SAT_STAGES, Pipeline
        sat_work = self.out / "satellite" / "work"
        self._check_chromosome_set(sat_work)          # first, before anything is written: a stopped run changes nothing
        view = inputs.make_view(self.genomes, self.gdir, self.work / "genomes" / self.gdir.name)
        force = self.force if self.force in SAT_STAGES else ("seeds" if self.force in ("satellite", "all") else None)
        try:
            pipe = Pipeline(view, self.threads, self.seed_cache,
                            Params(seed_engine=self.seed_engine, classify_from=self.classify_from, **self.sat_params),
                            force=force, out_dir=self.out / "satellite")
            # after construction (force wipe, stale check done) every checkpoint, present or future, stands on this set
            self._record_chromosome_set(sat_work)
            pipe.run()
        except SatelliteError as e:
            raise CenfuseError(f"satellite engine: {e}") from None
        return self.out / "satellite"

    def stage_telomere(self, g) -> dict:
        """Telomere unit, end arrays and interstitial arrays (telomere.py), checkpointed by the assembly index and the
        code (the interstitial scan reads the whole genome: about a minute per 3 Gb)."""
        fp = dict(code=code_key("telomere.py", stage="telomere"), fai=self._assembly(g),
                  chroms=obj_fingerprint(g.chroms))
        hit = self._cached("telomere", g.name, fp)
        if hit is None:
            hit = telomere.detect(g.fasta, g.chroms)
            self._store("telomere", g.name, fp, hit)
        return hit

    def stage_rdna(self, g) -> list:
        """45S arrays (NORs) from the universal 18S motif (rdna.py), checkpointed by the assembly index and the code."""
        fp = dict(code=code_key("rdna.py", stage="rdna"), fai=self._assembly(g),
                  chroms=obj_fingerprint(g.chroms))
        hit = self._cached("rdna", g.name, fp)
        if hit is None:
            hit = rdna.detect(g.fasta, g.chroms)
            self._store("rdna", g.name, fp, hit)
        self._nors[g.name] = hit
        return hit

    def stage_rdna5s(self, g) -> list:
        """5S rDNA arrays from the Rfam 5S model (rdna.detect_5s), checkpointed by the assembly index and the code. They
        are features only: unlike NORs, 5S arrays do not gather in the nucleolus, so Hi-C does not use them."""
        fp = dict(code=code_key("rdna.py", stage="rdna5s", data=("sat/data/rfam_rrna.cm",)), fai=self._assembly(g),
                  chroms=obj_fingerprint(g.chroms))
        hit = self._cached("rdna5s", g.name, fp)
        if hit is None:
            hit = rdna.detect_5s(g.fasta, g.chroms, self.threads)
            self._store("rdna5s", g.name, fp, hit)
        return hit

    def stage_gaps(self, g) -> dict:
        """Assembly gaps (runs of >= GAP_MIN_N N) per chromosome (morphology.scan_gaps), checkpointed by the assembly
        index and the code: one pass over the chromosomes, seconds per Gb."""
        fp = dict(code=code_key("morphology.py", stage="gaps"), fai=self._assembly(g),
                  chroms=obj_fingerprint(g.chroms))
        hit = self._cached("gaps", g.name, fp)
        if hit is None:
            hit = morphology.scan_gaps(g.fasta, g.chroms)
            self._store("gaps", g.name, fp, hit)
            n = sum(len(v) for v in hit.values())
            log.info(f"gaps: {g.name}: {n} assembly gaps on {sum(1 for v in hit.values() if v)} chromosomes")
        return hit

    def stage_attrough(self, g, seq: seqcalls.SeqLine) -> dict | None:
        """The AT-trough line (atline.py), only in a genome where most chromosomes have no sequence call (a long-unit
        array alone is none: fuse.decide) or a low one with a competing cluster the line could choose (fuse._at_places): elsewhere
        the line would never act, and it reads the whole genome. Checkpointed by the assembly index and the code."""
        missing = seq.calls.start.lt(0) | seq.calls.evidence.eq("long_unit_tandem") | \
            (seq.calls.confidence.eq("low") & seq.calls.alt_start.ge(0))
        if 2 * int(missing.sum()) <= len(seq.calls) or len(g.chroms) > MAX_SEQUENCES:
            return None
        fp = dict(version=__version__, code=code_key("atline.py", "sat/karyotype.py", stage="attrough"), fai=self._assembly(g),
                  chroms=obj_fingerprint(g.chroms))
        hit = self._cached("attrough", g.name, fp)
        if hit is None:
            t0 = time.time()
            hit = atline.scan(atline.gc_bins(g.fasta, g.chroms), g.chroms)
            self._store("attrough", g.name, fp, hit)
            log.info(f"attrough: {g.name} done ({time.time() - t0:.0f} s)")
        return hit

    def stage_genotype(self, g, samples: list, windows: list) -> dict:
        """Fixed differences of every methylation sample and the private-difference tests (provenance)."""
        files = sorted(str(f) for s in samples for f in s.files.values())
        fp = dict(version=__version__, code=code_key("provenance.py", stage="genotype"), fai=self._assembly(g),
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

    def stage_hic(self, g, lib, seq: seqcalls.SeqLine, sat_dir: Path, classes: dict | None, te_cores: dict | None = None):
        name = g.name if lib.name == "hic" else f"{g.name}.{lib.name}"
        calls = {r.chrom: (int(r.start), int(r.end), r.confidence) for r in seq.calls.itertuples(index=False)}
        gff = self.work / "hic" / f"{g.name}.sequence_calls.gff3"
        gff.parent.mkdir(parents=True, exist_ok=True)
        with atomic_open(gff) as f:                   # Hi-C engine --satellites: centromere features for its sat_mid
            f.write("##gff-version 3\n")
            for c, (s, e, _) in calls.items():
                if s >= 0:
                    f.write(f"{c}\tcenfuse\tcentromere\t{s + 1}\t{e}\t.\t.\t.\tID=cen_{c}\n")
        arrays = sat_dir / f"{g.name}_satellite_arrays.gff"
        fai = self.work / "genomes" / f"{g.name}.chromosomes.fai"
        fai.parent.mkdir(parents=True, exist_ok=True)
        txt = "".join(f"{c}\t{L}\t0\t0\t0\n" for c, L in g.chroms.items())
        if not fai.exists() or fai.read_text() != txt:     # the Hi-C engine's matrix cache is keyed by the file's mtime
            atomic_write_text(fai, txt)
        fp = dict(version=__version__, code=code_key("hicfuse.py", packages=("hic",), stage="hic"),
                  pairs=file_fingerprint(lib.path) if not lib.more else obj_fingerprint(
                      [file_fingerprint(f) for f in lib.files]), chroms=obj_fingerprint(g.chroms),
                  calls=obj_fingerprint(calls), arrays=file_fingerprint(arrays), classes=obj_fingerprint(classes),
                  nors=obj_fingerprint(self._nors.get(g.name, [])), te=obj_fingerprint(te_cores))
        hit = self._cached("hic", name, fp)
        if hit is not None:
            mc = hit.get("matrix_cache")
            if mc and not (Path(f"{self.out / 'hic' / name}.cache") / Path(mc).name).exists():
                # the arm fold and the alignment read the engine's contact matrix (stage_fold): without it they would be
                # dropped from a resumed run, so its results would differ from the first one's
                log.warning(f"hic: {name}: checkpoint found but its contact-matrix cache is gone; Hi-C is recomputed so the "
                            f"Rabl-axis lines can read it")
            else:
                log.info(f"hic: {name} checkpoint found")
                return hit
        t0 = time.time()
        (self.out / "hic").mkdir(parents=True, exist_ok=True)
        res = hicfuse.run(lib.path if not lib.more else lib.files, fai, self.out / "hic" / name, [gff, arrays], calls,
                          min(self.threads, HIC_MAX_THREADS), classes, nors=self._nors.get(g.name, []), te_cores=te_cores)
        self._store("hic", name, fp, res)
        log.info(f"hic: {name} done ({time.time() - t0:.0f} s)")
        return res

    def stage_fold(self, g, lib, hic: dict) -> dict | None:
        """Arm-fold calls (fold.py) at the Hi-C run's bin size, from the Hi-C engine's matrix cache (or the .cool),
        validated against the de novo Hi-C peaks. Checkpointed apart from hic: one ICE and a dense cis block per
        chromosome, seconds to minutes per genome."""
        name = g.name if lib.name == "hic" else f"{g.name}.{lib.name}"
        b, t = hic["bins"], hic["table"]
        seqs = {c: int(L) for c, L in zip(b.names, b.lengths)}
        cool = is_cool(lib.path)       # the cache file under this run's folder (the folder may have moved since)
        npz = None if cool else Path(f"{self.out / 'hic' / name}.cache") / Path(hic["matrix_cache"]).name
        if npz is not None and not npz.exists():
            log.warning(f"{g.name}: Hi-C matrix cache {npz} missing; arm fold skipped")
            return None
        den = t[["chrom", "callable", "denovo_peak_start", "denovo_peak_end", "denovo_core_start", "denovo_core_end",
                 "denovo_confidence"]].values.tolist()
        fp = dict(version=__version__, code=code_key("fold.py", packages=("hic",), stage="fold"),
                  matrix=(obj_fingerprint([file_fingerprint(f) for f in lib.files]) if lib.more else
                          file_fingerprint(lib.path)) if cool else file_fingerprint(npz),
                  size=int(b.size), chroms=obj_fingerprint(seqs), members=obj_fingerprint(hic["members"].tolist()),
                  denovo=obj_fingerprint(den))
        hit = self._cached("fold", name, fp)
        if hit is not None:
            log.info(f"fold: {name} checkpoint found")
            return hit
        if cool:
            M = sum(read_cool(f, b.size, seqs) for f in lib.files)
        else:
            M = sp.load_npz(npz).tocsr()
            base = Bins(seqs, int(hic["base_size"]))
            if b.size != base.size:
                M = coarsen(M, base, b)
        bias, good, _, _ = ice(M)
        W = balanced(M, bias, good)
        calls = fold.scan(W, b, good, hic["members"])
        res = dict(calls=calls, **fold.validate(calls, t))
        cen, cores = fold.align_centers(t, b)              # partners of the alignment: the prior-free peaks
        acalls = fold.align_scan(W, b, good, block_means(W, b, good, hic["members"]), hic["members"], cen, cores)
        res["align"] = dict(calls=acalls, **fold.validate(acalls, t, label="trans arm alignment"))
        self._store("fold", name, fp, res)
        return res

    def stage_karyotype(self, g, sat_dir: Path, w: int, seq: seqcalls.SeqLine) -> dict:
        """Windows not judged for missing segments (the sequence calls, gaps, satellite-dominated windows, ends), and
        each methylation sample's coverage before the depth filter."""
        calls = {}
        for r in seq.calls.itertuples(index=False):
            if r.start >= 0:
                calls.setdefault(r.chrom, []).append((int(r.start), int(r.end)))
        fp = dict(version=__version__, code=code_key("provenance.py", "methfuse.py", packages=("meth",), stage="karyotype"),
                  fai=self._assembly(g), w=w, calls=obj_fingerprint(calls),
                  arrays=obj_fingerprint(seq.arrays[["chrom", "start", "end"]].values.tolist()),
                  meth=self._meth_inputs_fp(g, sat_dir) if g.meth else None)
        hit = self._cached("karyotype", g.name, fp)
        if hit is not None:
            log.info(f"karyotype: {g.name} checkpoint found")
            return hit
        t0 = time.time()
        cov = {}
        if self._chg_ok.get(g.name):
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
        fp = dict(version=__version__, code=code_key("methfuse.py", packages=("meth",), stage="meth"),
                  assembly=self._assembly(g), width_min=width_min, chroms=obj_fingerprint(g.chroms),
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

    def stage_cpg(self, g, files: list, calls: dict, seq: seqcalls.SeqLine) -> dict:
        """The CpG line (cdr.py): a kinetochore inside each final call where CpG dips are validated genome-wide.
        Checkpointed by its CG files, the final calls and the arrays."""
        arr = seq.arrays[["chrom", "start", "end", "role"]]
        fp = dict(version=__version__, code=code_key("cdr.py", "meth/params.py", "meth/io.py", stage="cpg"),
                  files=[file_fingerprint(f) for f in files], calls=obj_fingerprint(calls),
                  arrays=obj_fingerprint(arr.values.tolist()))
        hit = self._cached("cpg", g.name, fp)
        if hit is not None:
            log.info(f"cpg: {g.name} checkpoint found")
            return hit
        t0 = time.time()
        res = cdr.run(files, calls, arr)
        self._store("cpg", g.name, fp, res)
        log.info(f"cpg: {g.name} done ({time.time() - t0:.0f} s)")
        return res

    # ------------------------------------------------------------------ assembly
    @staticmethod
    def evidences(g, seq: seqcalls.SeqLine, tel: dict, hic, meth_items=None, nors=None, at=None, ages=None) -> list:
        marker = tel.get("marker")
        sig = tel.get("signature")
        hrows = {} if hic is None else hic_rows(hic["table"], nors)
        fo = (hic or {}).get("fold")
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
        pairs = (hic or {}).get("te_anchored") or {}             # hicfuse.anchored_on_cores: the LTR cores tested by Hi-C
        for r in seq.calls.itertuples(index=False):
            te = None
            if marker and r.chrom in marker[1]:
                core = marker[1][r.chrom]
                te = dict(clade=core.get("clade", marker[0]), start=int(core["start"]), end=int(core["end"]),
                          decisive=bool(marker[2]), strong=core.get("strong", True), share=core.get("share"),
                          n=core.get("n"), seed=sig.seed if sig is not None else None)
            elif r.te_start >= 0 and sig is not None:          # a core too thin to count, where the signature is known (weak_core)
                te = dict(clade=r.te_clade, start=int(r.te_start), end=int(r.te_end), decisive=False, strong=False)
            elif r.te_start >= 0:                              # imported satellite calls made with tes/
                te = dict(clade=r.te_clade, start=int(r.te_start), end=int(r.te_end),
                          decisive=r.te_clade.lower() == "crm")
            h = hrows.get(r.chrom)
            if h is not None and not h["callable"]:
                h = None
            if h is not None and pairs.get("status") in ("pass", "support") and r.chrom in pairs["peaks"]:
                x = sig.cores[r.chrom]
                h = dict(h, te_pair=dict(core=(int(x["start"]), int(x["end"])), status=pairs["status"], clade=marker[0]))
                if te is None:                                 # a core that acts only through Hi-C: the tables report it
                    te = dict(clade=x.get("clade", marker[0]), start=int(x["start"]), end=int(x["end"]), decisive=False, strong=False,
                              share=x.get("share"), n=x.get("n"), seed=sig.seed)
            c_entry = cis_entry((hic or {}).get("cis"), r.chrom)
            if h is not None and c_entry is not None:
                h = dict(h, cis=c_entry)
            meth = None if meth_items is None else meth_items.get(r.chrom, {})
            f = fo["calls"].get(r.chrom) if fo else None
            f = dict(f, status=fo["status"], agrees=fo["agrees"].get(r.chrom)) if f else None
            al = (fo or {}).get("align")
            a = al["calls"].get(r.chrom) if al else None
            a = dict(a, status=al["status"], agrees=al["agrees"].get(r.chrom), fold_status=fo["status"],
                     fold_p=fo["p"]) if a else None
            out.append(fuse.Evidence(r.chrom, int(r.chrom_len), dict(start=int(r.start), end=int(r.end),
                                                                    evidence=r.evidence, confidence=r.confidence,
                                                                    flags=r.flags, main_family=r.main_family,
                                                                    alt_start=int(r.alt_start),
                                                                    alt_end=int(r.alt_end),
                                                                    overruled_start=int(r.overruled_start),
                                                                    overruled_end=int(r.overruled_end),
                                                                    spacing=int(r.spacing),
                                                                    te_core_elements=int(r.te_core_elements),
                                                                    te_call_elements=int(r.te_call_elements),
                                                                    rejected_family=r.rejected_family,
                                                                    rejected_start=int(r.rejected_start),
                                                                    rejected_end=int(r.rejected_end)),
                                     te, h, meth, sat_bp(r.chrom), f, a, atline.on_chromosome(at, r.chrom),
                                     age=(ages or {}).get(r.chrom)))
        return out

    def provenance_of(self, g, seq, sat_dir: Path, tel: dict, telo: dict):
        """Provenance of every input, then the Hi-C libraries: (samples, geno, hics, kar, profiles, segments,
        verdicts)."""
        samples = []
        if g.meth:
            try:
                samples = mio.read_samples(g.meth, require_chg=False)
            except mio.InputError as e:
                log.warning(f"{g.name}: methylation line skipped: {e}")
                g.notes.append(f"methylation skipped: {e}")
                g.meth = None
        samples, dropped = provenance.drop_non_modkit_chg(samples)
        for name, f, what in dropped:
            log.warning(f"{g.name}: methylation {name} CHG ({f}) is {what}: the CHG line needs modkit bedMethyl (its "
                        f"provenance test reads modkit's counts of failed and differing calls); CHG not used")
            g.notes.append(f"methylation {name}: CHG not modkit bedMethyl, not used")
        chg = [s for s in samples if "CHG" in s.files]
        self._chg_ok[g.name] = bool(g.tes and samples and len(chg) == len(samples))
        if samples and g.tes and not self._chg_ok[g.name]:       # without tes/ the input note says it already
            log.info(f"{g.name}: the CHG methylation line needs a CHG file for every sample: skipped" +
                     ("; CG files locate kinetochores (CpG line)" if any("CG" in s.files for s in samples) else ""))
        # the same D as the methylation line's (methfuse.prepare): it divides by the chromosomes cenfuse analyzes
        D = mtrack.chromosome_set(mio.read_fai(f"{g.fasta}.fai"), n=len(g.chroms))[2]
        avoid = avoid_intervals(seq, D)
        geno = {}
        if samples:
            windows = provenance.genotype_windows(g.chroms, avoid, g.fasta)
            if chg:                                     # genotypes are read at CHG sites
                geno = self.stage_genotype(g, chg, windows)
            cg_only = [s for s in samples if "CHG" not in s.files]
            warn = geno.get("warnings", []) + (provenance.check_assembly(cg_only, g.fasta, windows) if cg_only else [])
            for smp in samples:                         # bedGraph CG files (no strand or context): CpG positions
                f = smp.files.get("CG")
                if f is None or cdr.track_kind(f) != "bedgraph":
                    continue
                share, n = cdr.bedgraph_on_cpg(f, g.fasta, windows)
                if n and share < provenance.MOTIF_MIN:
                    raise CenfuseError(f"methylation {smp.name} CG ({f}): only {share:.0%} of {n} rows sit on "
                                       f"reference CpGs; the file belongs to another assembly")
                if n and share < provenance.MOTIF_WARN:
                    warn.append(f"methylation {smp.name} CG: {share:.1%} of {n} rows on reference CpGs")
            for x in warn:
                log.warning(f"{g.name}: {x}")
        known = {n: v["verdict"] for n, v in geno.get("verdicts", {}).items()}
        classes = feature_classes(seq, tel, telo, self._nors.get(g.name, []))
        sig = tel.get("signature")
        te_cores = {c: (int(x["start"]), int(x["end"])) for c, x in sig.cores.items()} if sig is not None else None
        hics = {}
        for lib in g.hic:
            try:
                hics[lib.name] = self.stage_hic(g, lib, seq, sat_dir, None if known.get(lib.name) == "other"
                                                else classes, te_cores)
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
        """Lock the output folder (a second run on it stops; a stale lock is taken over), write the run manifest
        (run_info.json), run, and record how the run ended in the manifest, also when it fails."""
        t0 = time.time()
        self.out.mkdir(parents=True, exist_ok=True)
        lock = self.out / LOCK
        if not try_lock(lock, f"cenfuse run {self.gdir}"):
            try:
                owner = lock.read_text().strip()
            except OSError:
                owner = "unknown"
            raise CenfuseError(f"{self.out} is in use by another cenfuse run ({owner}); wait for it to end, or delete "
                               f"{lock} if no run is active")
        try:
            info = manifest.start(self.out, self)
            use = Usage()
            try:
                with Heartbeat(lock), use:
                    incomplete = self._run(t0)
            except BaseException as e:
                manifest.finish(self.out, info, "failed", t0, error=f"{type(e).__name__}: {e}", cost=use.cost())
                raise
            cost = use.cost()
            manifest.finish(self.out, info, "done", t0, incomplete=incomplete, cost=cost)
            log.info(f"cost: {cost['cpu_hours']:.2f} CPU h"
                     + (f", peak memory {cost['peak_memory_gb']:.1f} GB ({cost['memory']})"
                        if cost["peak_memory_gb"] is not None else "") + f" ({manifest.NAME})")
        finally:
            release_lock(lock)

    def _run(self, t0: float) -> list:
        """The stages and outputs of every genome; returns the outputs that failed ('<genome>: report', 'figure'), which
        are also in cenfuse_summary.tsv (incomplete) and the report: a failed report or figure never costs the tables, and
        is never silent."""
        log.info(f"cenfuse {__version__}: {self.gdir.name}: " + "; ".join(
            f"{g.name} ({len(g.chroms)} chromosomes; {', '.join(g.lines)})" for g in self.genomes))
        for g in self.genomes:
            for n in g.notes:
                log.info(f"{g.name}: {n}")
        sat_dir = self.stage_satellite()
        summary, figs, sections, failed = [], [], [], {}
        for g in self.genomes:
            seq = seqcalls.read(sat_dir, g.name, g.chroms)
            tel = teline.te_line(g.tes, g.chroms, sat_dir, g.name)
            telo = self.stage_telomere(g)
            spare = inputs.left_out(g.fasta, g.fai, g.chroms, telo.get("unit"))       # future.md #8: flagged, never added
            if spare:
                names = ", ".join(f"{c} ({L / 1e6:.2f} Mb)" for c, L in spare[:10]) + (f" and {len(spare) - 10} more" if len(spare) > 10 else "")
                g.notes.append(f"{len(spare)} sequence{'s' if len(spare) > 1 else ''} left out of the chromosome set "
                               f"hold{'' if len(spare) > 1 else 's'} telomere arrays at both ends ({names}): complete molecules, "
                               f"such as a small chromosome or a B chromosome; list them in {g.name}.chromosomes.txt to analyze them")
                log.warning(f"{g.name}: {g.notes[-1]}")
            telo = dict(telo, left_out=spare)
            nors = self.stage_rdna(g)
            five_s = self.stage_rdna5s(g)
            gaps = self.stage_gaps(g)
            at = self.stage_attrough(g, seq)
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
                hic["fold"] = self.stage_fold(g, next(x for x in g.hic if x.name == dec_lib), hic)
                hic["axis"] = fuse.genome_axis(hic["fold"])
                hic["conflict"] = c = fuse.sequence_conflict(hic["table"])
                folded = hic["axis"] is not None
                if c["conflict"]:
                    log.warning(f"{g.name}: Hi-C's prior-free calls are confident on {c['denovo_confident']}/{c['n']} "
                                f"chromosomes but agree with the sequence calls on only {c['denovo_on_seq']}, and the "
                                f"gate fails: the two lines disagree genome-wide, so at least one is wrong. " +
                                (f"The {'arm fold (cis Hi-C)' if hic['axis'] == 'fold' else 'trans arm alignment'}"
                                 " sides with the prior-free calls: sequence calls both contradict are moved or "
                                 "lowered" if folded else
                                 "The sequence calls are kept; give tes/ (LTR cores, against which Hi-C is then tested) to decide") +
                                " (cenfuse_summary.tsv: hic_sequence_conflict, hic_axis)")
            at = atline.with_hic(at, hic, g.chroms)
            ages = ageline.for_genome(g.tes, g.chroms, tel.get("signature"), hic["table"] if hic is not None else None,
                                      seq.calls) if g.tes else {}
            evs = self.evidences(g, seq, tel, hic, nors=self._nors.get(g.name, []), at=at, ages=ages)
            prov = [fuse.decide(e) for e in evs]
            meth, items = None, {}
            if self._chg_ok.get(g.name):
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
                evs = self.evidences(g, seq, tel, hic, mi, nors=self._nors.get(g.name, []), at=at, ages=ages)
            final = [fuse.decide(e) for e in evs]
            t = fuse.table(evs, final)
            calls = {e.chrom: d["call"] for e, d in zip(evs, final)}
            cg = [s.files["CG"] for s in samples if "CG" in s.files and s.name in deciders]
            cpg = None
            if cg:
                fin = {e.chrom: (*d["call"], d["confidence"]) if d["call"] else (-1, -1, "none")
                       for e, d in zip(evs, final)}
                try:
                    cpg = self.stage_cpg(g, cg, fin, seq)
                except (OSError, ValueError) as e:
                    log.warning(f"{g.name}: CpG line skipped: {e}")
                    g.notes.append(f"CpG line skipped: {e}")
            elif samples:
                log.info(f"{g.name}: no CG file from a deciding methylation sample; CpG line not run")
            t = cdr.columns(t, cpg)
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
            t = morphology.columns(t, g.chroms, telo, gaps)
            ictx = interpret.context(seq, tel, bool(g.hic), bool(g.tes), telo)
            t = interpret.columns(t, ictx)
            output.write_calls(t, g.name, g.chroms, self.out)
            lines = linetests.table(t, int(hic["bin_size"]) if hic is not None else None)
            output.write_tsv(lines, self.out / f"{g.name}.line_tests.tsv")
            sig = tel.get("signature")
            if sig is not None:
                signature.write_table(sig, self.out / f"{g.name}.te_signature.tsv")
            ends = {n: int(np.sum(h["coverage"])) for n, h in hics.items()}
            ptable = provenance.table(verdicts, geno, segments, deciders if meth is not None else [], dec_lib, ends)
            output.write_tsv(ptable, self.out / f"{g.name}.provenance.tsv")
            if len(ptab):
                output.write_tsv(plants.one_based(ptab, ("alt_start", "hic_peak_start")),
                                 self.out / f"{g.name}.plants.tsv")
            if psites is not None:
                output.write_tsv(plants.one_based(psites, ("start", "mid")), self.out / f"{g.name}.plant_sites.tsv")
            fused = {e.chrom: dict(call=d["call"], relation=d["relation"], seq_evidence=e.seq["evidence"],
                                   seq=(e.seq["start"], e.seq["end"]) if e.seq["start"] >= 0 else None,
                                   long_unit_dropped="long_unit_unsupported" in d["flags"])
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
                                   meth, genes, g.fasta, telo, nors, five_s)
            output.write_loci(lt, feats, g.name, g.chroms, self.out)
            fc = hic["classes"] if hic is not None else None
            if fc is not None:
                output.write_tsv(fc, self.out / f"{g.name}.feature_contacts.tsv")
            output.write_tsv(loci.feature_signatures(lt, fc), self.out / f"{g.name}.feature_signatures.tsv")
            summary.append(output.summary_row(g.name, g, t, hic, meth, tel, telo, verdicts, ptab, nors, cpg, five_s, at, seq.test))
            try:
                sections.append(interpret.genome_section(g.name, g, t, summary[-1], lt, fc, ictx, lines))
            except Exception as e:                   # the report never costs the results, and says what is missing
                log.warning(f"{g.name}: report section failed ({type(e).__name__}: {e}); the tables are complete")
                sections.append(f"## {g.name}\n\nThis section could not be written ({type(e).__name__}: {e}). The tables "
                                f"(`{g.name}.centromeres.tsv`, `{g.name}.loci.tsv`) are complete.\n")
                failed.setdefault(g.name, []).append("report")
            figs.append(dict(genome=g, seq=seq, tel=tel, hic=hic, meth=meth, calls=t, loci=lt, telo=telo,
                             plants=ptab, psites=psites, verdicts=verdicts, libraries=others,
                             rdna=dict(nors=list(nors or ()), five_s=list(five_s or ()))))
            log.info(f"{g.name}: {len(t[t.start != '.'])}/{len(t)} chromosomes called; relations "
                     f"{t.relation.value_counts().to_dict()}; confidence "
                     f"{t[t.start != '.'].confidence.value_counts().to_dict()}")
            if len(ptab):
                sh = ptab[ptab.verdict.isin(["shifted", "call_absent"])]
                log.info(f"{g.name}: plants: {ptab.plant.nunique()} scored; "
                         + ("; ".join(f"{r.plant} {r.chrom} {r.verdict}" for r in sh.itertuples()) or "no shift"))
        if sections:
            try:
                interpret.write_report(sections, self.out, self.gdir)
            except Exception as e:                   # the report never costs the results
                log.warning(f"report failed ({type(e).__name__}: {e}); the tables are complete")
                for g in self.genomes:
                    failed.setdefault(g.name, []).append("report")
        try:
            from . import plot
            plot.overview(figs, self.out)
        except Exception as e:                       # a figure never costs the results
            log.warning(f"overview figure failed ({type(e).__name__}: {e}); the tables are complete")
            for g in self.genomes:
                failed.setdefault(g.name, []).append("figure")
        for row in summary:
            row["incomplete"] = ",".join(dict.fromkeys(failed.get(row["genome"], []))) or "."
        output.write_tsv(pd.DataFrame(summary), self.out / "cenfuse_summary.tsv")
        lost = [f"{n}: {x}" for n, xs in failed.items() for x in dict.fromkeys(xs)]
        if lost:
            log.warning(f"done ({(time.time() - t0) / 60:.1f} min) with incomplete outputs ({'; '.join(lost)}); the tables "
                        f"are complete; outputs in {self.out}")
        else:
            log.info(f"done ({(time.time() - t0) / 60:.1f} min); outputs in {self.out}")
        return lost
