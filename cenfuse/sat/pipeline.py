"""Stage orchestration. Checkpoints: <complex>/centromere/work/<stage>.done; a rerun resumes by default."""
import json
import os
import time
from dataclasses import replace
from pathlib import Path

import pandas as pd
import pysam
import pyhmmer

from . import __version__, annotate, centromere, classify, families, hor, karyotype, output, plot, seeds, signature, te
from ..runtime import atomic_copy, atomic_csv, atomic_write_text, fasta_md5, file_fingerprint, unshare_tree
from .genomes import discover_genomes, gene_gff, load_genes
from .params import MAX_CHROMOSOMES, TE_DERIVED_FRAC, Params, format_borrow, parse_borrow, to_dict
from .util import SatelliteError, log, require_tools

STAGES = ["seeds", "families", "annotate", "classify", "centromeres", "hor", "output", "plot"]
# results version per stage: bump a stage when its code can change outputs; a checkpoint's fingerprint carries the
# highest version among it and earlier stages, so a bump invalidates that stage and later ones only.
# 0.3.0 = C1/C2/I1 review fixes; 0.3.2 = classify: strong-family precedence per chromosome;
# 0.3.3 = annotate: shortest-unit rule (drop_multimer_hits); 0.3.4 = annotate: coverage-based overlap resolution;
# 0.3.5 = families: shared monomer frame for related families + near-identical family merge (0.97);
#         centromeres: diffuse flag; plot: stacked bars capped at one bin;
# 0.3.6 = families: trim run-through (terminally redundant) consensus + chain groups; annotate: chain per group.
# 0.3.7 = classify: precedence only for dominant families (largest cluster on >= strong_prevalence of chromosomes).
# 0.4.0 = any number of genomes; annotate: evalue_ref_len; classify: classify_from (borrowed families for contig
#         assemblies) + genes/<genome>.gff3 gene loci; centromeres: borrowed genomes; plot: N panels.
# 0.4.1 = plot: gene density line grey, under the satellite bars.
# 0.4.2 = classify: 45S rDNA spacer repeats (arrays flanked by 18S/5.8S/25S genes) are never centromeric.
# 0.4.3 = inputs: bgzip *.fa.gz genomes (results unchanged: no stage version bump).
# 0.4.4 = families: models founded on a seed whose consensus spans its whole TRF period (base_unit).
# 0.6.0 = classify: dominant column; centromeres: retrotransposon evidence (tes/), honesty cap without it.
# 0.7.0 = centromeres: with decisive (CRM) evidence a satellite call away from the core is replaced by the core.
# 0.4.0 (seeds) = stage trf renamed seeds; FasTAN engine (default) beside TRF.
# 0.10.0 = classify: dispersed repeats (tandem test), dominance per satellite kind, acrocentric test;
#          centromeres: competing clusters reported (alt_start/alt_end) and low (T2T mammal benchmark).
# 0.10.1 = classify: the acrocentric test's interior rival is judged per satellite kind; centromeres: a second
#          acrocentric satellite unrelated to the call competes with it (T2T mouse).
# 0.10.2 = centromeres: a flags-only (non-CRM) LTR marker lifts the no-dominant cap only for calls it supports.
# 0.10.3 = centromeres: a decisive core replaces a weak satellite call beside it (weak_satellite_beside_core).
# 0.10.4 = centromeres: a significant CRM core replaces a flags-only marker's core on its chromosome.
# seeds 0.11.0: FasTAN runs with address-space randomization off (fastan.fixed_layout), so reruns are identical.
# families 0.12.0: clusters founded on two-copy seeds only are not families (families.build_families).
# annotate 0.12.0: follows families (no change of its own; a fingerprint takes the highest version up to its stage).
# classify 0.12.0: TE floor counts independent insertion sites (te.core); a decisive core replaces a confident
# satellite call only if the clade marks it significantly more (te.contradicts). centromeres 0.12.1: the same test for
# a call beside the core (satellite_beside_core). (0.10.5-0.11.0: v8 test variants.)
# annotate 0.13.0: long models search with their low-complexity tracts neutral (annotate.neutralize).
# classify 0.13.1: the period test of `related` is symmetric; dominant_satellite; te_rejected; ltr_clade counts every
# overlapping element. centromeres 0.13.1: calls on chains of partial hits (cap_chain_calls); rejected_family, element
# counts. output 0.13.1: copy spacing; new columns.
# classify 0.14.0: the LTR signature (families recruited leave one chromosome out; a specific lineage other than CRM may
# decide). centromeres 0.14.0: the one-per-chromosome test of the centromeric satellites caps calls; a core must be specific on
# its chromosome to overrule a confident call; chain calls narrow to a decisive core. output 0.14.0: te_signature, satellite_test.
# classify 0.15.0: the LTR signature's core is the youngest of several comparable peaks (age_choice), a family concentrated over the
# chromosomes can seed it (pooled_candidates), a signature whose cores are not specific on their chromosomes does not decide; centromeres
# 0.15.0: a significant minority core replaces a call that is low on its own (minority_core); output 0.15.0: core_age_from, core_alt.
# centromeres 0.16.0: where the satellites failed their test (or none dominates), the decisive core an array agrees with is the call
# (satellite_in_core).
STAGE_VERSION = {"seeds": "0.11.0", "families": "0.12.0", "annotate": "0.13.0", "classify": "0.15.0",
                 "centromeres": "0.16.0", "hor": "0.3.0", "output": "0.15.0", "plot": "0.4.1"}
# earliest stage each parameter affects; a stage's checkpoint fingerprint covers its own and all earlier stages'
STAGE_PARAMS = {
    "seeds": ["seed_engine", "seed_chunk", "seed_overlap"],
    "families": ["min_period", "max_period", "min_seed_span", "fam_ident", "multi_ident", "superfam_ident",
                 "absorb_ident", "dedup_ident", "chain_ident", "min_fam_span", "max_families", "hmm_sample", "hmm_sample_bp", "min_hmm_seqs",
                 "nhmmer_evalue", "array_gap_factor", "min_array_monomers"],
    "annotate": ["evalue_ref_len"],
    "classify": ["chrom_min_bp", "cluster_gap", "min_concentration", "max_subtel_frac",
                 "subtel_frac_len", "max_gene_ratio", "strong_prevalence", "strong_cluster_bp", "classify_from",
                 "rdna_flank", "rdna_arrays", "rdna_array_frac"],
    "centromeres": ["low_sat_bp", "competing_ratio", "fallback_min_bp", "cen_assoc_frac", "min_gap_n",
                    "diffuse_min_len", "diffuse_max_frac", "lu_flank",
                    "lu_window", "lu_window_overlap", "lu_min_unit", "lu_max_unit", "lu_min_ident",
                    "contig_end_bp"],
    "hor": ["hor_ident", "hor_min_cov", "hor_max_len"],
    "output": ["tandem_min_len"],
    "plot": [],
}
TOOLS = ["makeblastdb", "blastn", "dustmasker", "mafft", "cd-hit-est", "minimap2", "cmsearch"]        # + the seed engine's
# the first stage that reads the optional inputs genes/ (the gene veto) and tes/ (the LTR signature): from it on, a checkpoint is
# keyed also by those files' content, and the stages from it on are cheap next to the search, so a stale one is recomputed
# (adding tes/ to a finished run is the commonest reason) where a stale search stage stops the run
INPUT_STAGE = "classify"


class Pipeline:
    def __init__(self, complex_dir, threads: int, seed_cache, p: Params, force: str = None, out_dir=None):
        if p.seed_engine not in seeds.ENGINES:
            raise SatelliteError(f"seed engine {p.seed_engine!r} is not one of: {', '.join(seeds.ENGINES)}")
        self.cdir = Path(complex_dir).resolve()
        self.name = self.cdir.name
        self.out = Path(out_dir).resolve() if out_dir else self.cdir / "centromere"
        self.work = self.out / "work"
        self.threads = threads
        legacy = self.work / "trf_cache"                 # satellite engine < 0.9: its TRF results stay usable
        self.cache = Path(seed_cache).resolve() if seed_cache else legacy if legacy.is_dir() else self.work / "seed_cache"
        genomes = discover_genomes(self.cdir)
        self.borrow = self._check_borrow(parse_borrow(p.classify_from), {g.name for g in genomes})
        self.p = p = replace(p, classify_from=format_borrow(self.borrow))
        # sources are classified before the genomes that borrow from them
        self.genomes = [g for g in genomes if g.name not in self.borrow] + \
                       [g for g in genomes if g.name in self.borrow]
        for g in self.genomes:
            if g.name not in self.borrow and len(g.chroms) > MAX_CHROMOSOMES:
                log.warning(f"{g.name}: {len(g.chroms)} sequences; classify treats each as a chromosome "
                            f"(prevalence = share of sequences). Use a chromosome-level FASTA, or --classify-from "
                            f"for a contig assembly")
        self.work.mkdir(parents=True, exist_ok=True)
        self._inputs = None
        stale = [st for st in STAGES if self._done(st) and not self._same(st)]
        if stale and (force is None or STAGES.index(force) > STAGES.index(stale[0])):
            why = "; ".join(f"{st}: {self._why(st)}" for st in stale)
            if STAGES.index(stale[0]) < STAGES.index(INPUT_STAGE):
                raise SatelliteError(f"changed since checkpointed in {self.work} ({why}); the search stages are not redone "
                                     f"unasked: rerun with --force {stale[0]} (or an earlier stage)")
            log.warning(f"changed since checkpointed in {self.work} ({why}): recomputing from {stale[0]}")
            force = stale[0]
        self.force = force
        if force:
            for st in STAGES[STAGES.index(force):]:
                (self.work / f"{st}.done").unlink(missing_ok=True)
        for st in STAGES:                  # an accepted stamp of before 0.14 gets the full key: a later FASTA edit is then seen
            if self._done(st) and "fasta" not in self._stored(st) and self._same(st):
                atomic_write_text(self.work / f"{st}.done", json.dumps(self.fingerprint(st)))
        atomic_write_text(self.work / "params.json", json.dumps(to_dict(p), indent=1))
        self._genes, self._te, self._sig = {}, {}, {}

    @staticmethod
    def _check_borrow(borrow: dict, names: set) -> dict:
        for tgt, srcs in borrow.items():
            for n in [tgt] + srcs:
                if n not in names:
                    raise SatelliteError(f"--classify-from: unknown genome {n!r} (genomes: {', '.join(sorted(names))})")
            if tgt in srcs:
                raise SatelliteError(f"--classify-from: {tgt} borrows from itself")
            for s in srcs:
                if s in borrow:
                    raise SatelliteError(f"--classify-from: source {s} of {tgt} is itself borrowed")
        return borrow

    def fingerprint(self, stage: str) -> dict:
        """A stage's checkpoint key: the highest results version up to it, the parameters of it and earlier stages, and from
        INPUT_STAGE on the content fingerprints of every genome's genes/ and tes/ files (None when absent)."""
        cur = to_dict(self.p)
        keys = [k for st in STAGES[:STAGES.index(stage) + 1] for k in STAGE_PARAMS[st]]
        upto = STAGES[:STAGES.index(stage) + 1]
        version = max((STAGE_VERSION[st] for st in upto), key=lambda v: tuple(int(x) for x in v.split(".")))
        fp = {"version": version, "params": {k: cur[k] for k in keys},
              "fasta": {g.name: fasta_md5(g.fasta) for g in self.genomes}}
        if STAGES.index(stage) >= STAGES.index(INPUT_STAGE):
            fp["inputs"] = self.inputs_fp()
        return fp

    def _same(self, stage: str) -> bool:
        """The stage's stamp matches its key. A stamp made before 0.14 holds no FASTA md5: it is judged on the rest, and only
        while every genome's FASTA is older than the stamp (a FASTA written since may have been edited, and its content cannot
        be checked against the stamp); __init__ then rewrites such a stamp with the full key, so a later edit is caught by its
        md5. A stamp made since 0.14 must match the md5 too: a same-length edit (hard-masking, polishing) is never resumed on."""
        old, new = self._stored(stage), self.fingerprint(stage)
        if "fasta" not in old:
            new = {k: v for k, v in new.items() if k != "fasta"}
            if old == new and not self._fasta_older(stage):
                return False
        return old == new

    def _fasta_older(self, stage: str) -> bool:
        """Every genome's FASTA (the file a link points to) was last written before the stage's stamp."""
        try:
            t = (self.work / f"{stage}.done").stat().st_mtime
            return all(Path(os.path.realpath(g.fasta)).stat().st_mtime <= t for g in self.genomes)
        except OSError:
            return False

    def inputs_fp(self) -> dict:
        """{genome: {'genes': fingerprint or None, 'tes': fingerprint or None}} of the optional inputs classify onward read."""
        if self._inputs is None:
            out = {}
            for g in self.genomes:
                gg, tt = gene_gff(self.cdir, g.name), te.te_gff(self.cdir, g.name)
                out[g.name] = dict(genes=fasta_md5(gg) if gg else None, tes=fasta_md5(tt) if tt else None)   # whole files: edited by hand
            self._inputs = out
        return self._inputs

    def _why(self, stage: str) -> str:
        """What differs between a stage's checkpoint key and the current one, in words."""
        old, new = self._stored(stage), self.fingerprint(stage)
        out = []
        if old.get("version") != new.get("version"):
            out.append(f"results version {old.get('version')} -> {new.get('version')}")
        op, np_ = old.get("params") or {}, new.get("params") or {}
        ch = sorted(k for k in set(op) | set(np_) if op.get(k) != np_.get(k))
        if ch:
            out.append("parameters " + ", ".join(ch))
        if "fasta" not in old and not self._fasta_older(stage):
            out.append("the FASTA is newer than this checkpoint, made before cenfuse 0.14 keyed the FASTA's content")
        if "fasta" in old and old["fasta"] != new.get("fasta"):
            out.append("FASTA content (" + ", ".join(sorted(g for g in set(old["fasta"]) | set(new.get("fasta") or {})
                                                            if old["fasta"].get(g) != (new.get("fasta") or {}).get(g))) + ")")
        if old.get("inputs") != new.get("inputs"):
            oi = old.get("inputs")
            if oi is None:
                out.append("made before inputs were keyed")
            else:
                diff = sorted(f"{g} {k}" for g in set(oi) | set(new.get("inputs") or {})
                              for k in ("genes", "tes") if (oi.get(g) or {}).get(k) != ((new.get("inputs") or {}).get(g) or {}).get(k))
                out.append("inputs " + ", ".join(diff))
        return "; ".join(out) or "key differs"

    def _done(self, stage: str) -> bool:
        return (self.work / f"{stage}.done").exists()

    def _stored(self, stage: str) -> dict:
        try:
            return json.loads((self.work / f"{stage}.done").read_text() or "{}")
        except json.JSONDecodeError:
            return {}

    # ---------------------------------------------------------------- helpers
    def _p(self, stage: str, name: str) -> Path:
        return self.work / stage / name

    def genes(self, g):
        if g.name not in self._genes:
            self._genes[g.name] = load_genes(self.cdir, g.name, g.chroms)
        return self._genes[g.name]

    def seed_table(self, g):
        return seeds.load_seeds(self._p("seeds", f"{g.name}.seeds.tsv.gz"))

    def fam_table(self):
        return pd.read_csv(self._p("families", "families.tsv"), sep="\t",
                           dtype={"superfamily": str, "rdna": str, "chain_group": str},
                           keep_default_na=False)

    def fam_info(self) -> dict:
        f = self.fam_table()
        return {r.fam_id: dict(superfamily=r.superfamily, length=int(r.consensus_len), rdna=r.rdna,
                               chain_group=getattr(r, "chain_group", r.fam_id))
                for r in f.itertuples(index=False)}

    def monomers(self, g):
        return pd.read_csv(self._p("annotate", f"{g.name}.monomers.tsv.gz"), sep="\t", dtype={"chrom": str})

    def arrays(self, g):
        return pd.read_csv(self._p("annotate", f"{g.name}.arrays.tsv"), sep="\t", dtype={"chrom": str})

    def fams(self, g, stage: str = "centromeres"):
        name = g if isinstance(g, str) else g.name
        return pd.read_csv(self._p(stage, f"{name}.families.tsv"), sep="\t")

    def regions(self, g):
        return centromere.load_regions(self._p("centromeres", f"{g.name}.regions.tsv"))

    # ---------------------------------------------------------------- driver
    def run(self) -> None:
        require_tools(TOOLS + seeds.ENGINES[self.p.seed_engine].TOOLS)
        if not families.RRNA_CM.exists():
            raise SatelliteError(f"rRNA models missing: {families.RRNA_CM} (part of cenfuse; reinstall it)")
        t0 = time.time()
        log.info(f"satellite engine {__version__}: {self.name} [{', '.join(g.name for g in self.genomes)}], "
                 f"{self.threads} threads")
        rerun = False                                   # once a stage runs, every later stage runs too
        for s in STAGES:
            if not rerun and self._done(s):
                log.info(f"{s}: checkpoint found, skipping")
                continue
            rerun = True
            t = time.time()
            log.info(f"{s}: start")
            (self.work / s).mkdir(exist_ok=True)
            n = unshare_tree(self.work / s)               # a hard-linked copy's files are not rewritten
            if n:
                log.debug(f"{s}: {n} files shared with another folder (hard links) unlinked before recomputing")
            getattr(self, f"stage_{s}")()
            atomic_write_text(self.work / f"{s}.done", json.dumps(self.fingerprint(s)))
            log.info(f"{s}: done ({time.time() - t:.0f} s)")
        log.info(f"all stages done ({(time.time() - t0) / 60:.1f} min); outputs in {self.out}")

    # ---------------------------------------------------------------- stages
    def stage_seeds(self):
        paths = seeds.seed_genomes(self.genomes, self.cache, self.threads, self.p)
        for g in self.genomes:
            link = self._p("seeds", f"{g.name}.seeds.tsv.gz")
            link.unlink(missing_ok=True)
            link.symlink_to(paths[g.name])
            crashes = seeds.crash_path(paths[g.name])           # exists only if FasTAN crashed on this genome
            link = self._p("seeds", f"{g.name}.fastan_crashes.tsv")
            link.unlink(missing_ok=True)                        # a link left by an earlier run of this stage
            if crashes.exists():
                link.symlink_to(crashes)

    def stage_families(self):
        families.build_families({g.name: self.seed_table(g) for g in self.genomes}, self.genomes,
                                self.work / "families", self.threads, self.p)

    def stage_annotate(self):
        fam = self.fam_table()
        hmms = []
        if len(fam):
            with pyhmmer.plan7.HMMFile(self._p("families", "families.hmm")) as f:
                hmms = list(f)
        fam_len = dict(zip(fam.fam_id, fam.consensus_len))
        groups = dict(zip(fam.fam_id, fam.chain_group)) if "chain_group" in fam else None
        for g in self.genomes:
            annotate.annotate_genome(g, hmms, fam_len, self.threads, self.p,
                                     self._p("annotate", f"{g.name}.monomers.tsv.gz"),
                                     self._p("annotate", f"{g.name}.arrays.tsv"), groups,
                                     out_cost=self._p("annotate", f"{g.name}.search_cost.tsv"))

    def stage_classify(self):
        """Family character is judged over the complex (one family = one sequence): the subtelomeric vote counts the
        chromosome-level genomes where a family is terminal (classify.family_metrics: an interior genome is never
        overruled; terminal blocks on most chromosomes of a genome decide there); an rDNA-spacer tag in one genome
        holds in all (rRNA genes are often missing from an assembly, so their absence next to the arrays is weak
        evidence). With tes/, each family's share of array bp inside LTR retrotransposons is recorded (in_ltr)."""
        info = self.fam_info()
        vote = {}                                         # subtelomeric votes over the chromosome-level genomes
        for g in self.genomes:
            if g.name not in self.borrow:
                for f, (s, b) in classify.terminal_share(self.arrays(g), g.chroms, self.genes(g), self.p).items():
                    if b >= self.p.chrom_min_bp:
                        k, n = vote.get(f, (0, 0))
                        vote[f] = (k + (s >= self.p.max_subtel_frac * b), n + 1)
        metrics, rds, elsewhere = {}, {}, {}
        for g in self.genomes:
            mono = pd.read_csv(self._p("annotate", f"{g.name}.monomers.tsv.gz"), sep="\t",
                               usecols=["family", "start", "end", "hmm_from", "hmm_to", "model_len"])
            tandem = classify.tandem_support(self.arrays(g), mono, self.seed_table(g), self.p)
            del mono
            metrics[g.name] = classify.family_metrics(self.arrays(g), g.chroms, self.genes(g), self.p, info, vote,
                                                      tandem)
            if g.name not in self.borrow:
                fams = metrics[g.name][0]
                elig = list(fams.family[fams.eligible.astype(bool)])
                rds[g.name] = rd = families.rdna_neighborhood(g.fasta, g.chroms, self.arrays(g), elig,
                                                               self._p("classify", f"rdna_{g.name}"), self.threads,
                                                               self.p.rdna_flank, self.p.rdna_arrays)
                for f, (n, k) in rd.items():
                    if n and k / n >= self.p.rdna_array_frac and f not in elsewhere:
                        elsewhere[f] = (g.name, k, n)
        for g in self.genomes:
            fams, mains = metrics[g.name]
            if g.name in self.borrow:
                src = {s: self.fams(s, "classify") for s in self.borrow[g.name]}
                fams, per_chrom = classify.borrow_centromeric(fams, mains, src, self.p)
                how = f" (borrowed from {'+'.join(src)})"
            else:
                rd = rds[g.name]
                other = {f: v for f, v in elsewhere.items() if v[0] != g.name}
                fams = classify.apply_rdna(fams, rd, self.p, other)
                tagged = [f"{f} ({k}/{n})" for f, (n, k) in rd.items() if n and k / n >= self.p.rdna_array_frac]
                if tagged:
                    log.info(f"classify: {g.name}: rDNA spacer repeats (not centromeric): {', '.join(tagged)}")
                tev = self._te_verdict(g, fams, mains, info)
                fams, per_chrom, rejected = classify.choose_with_te(fams, mains, g.chroms, self.p, info, tev,
                                                                    self.arrays(g))
                if rejected:
                    log.info(f"classify: {g.name}: away from the {tev[0]} retrotransposon centromere (not "
                             f"centromeric): {', '.join(f'{f} ({k}/{n})' for f, (k, n) in rejected.items())}")
                how = ""
            gff = te.te_gff(self.cdir, g.name)
            if gff is not None:                               # share of each family inside LTR retrotransposons
                iv = te.load_te_intervals(gff, g.chroms)
                longest = max((int((e - s).max()) for s, e, _ in iv.values()), default=0)
                ov = te.ltr_overlap(self.arrays(g), iv, longest)
                fams["in_ltr"] = [round(ov.get(f, (float("nan"), "."))[0], 3) for f in fams.family]
                fams["ltr_clade"] = [ov.get(f, (0, "."))[1] for f in fams.family]
            atomic_csv(fams, self._p("classify", f"{g.name}.families.tsv"), sep="\t", index=False)
            prim = fams[fams.cls == "centromeric"]
            log.info(f"classify: {g.name}: centromeric{how} = {', '.join(prim.family) or 'none'}; "
                     f"{len(per_chrom)}/{len(g.chroms)} sequences have one")

    def te_evidence(self, g):
        """(clade, {chrom: significant core}, decisive) of the centromere-marking LTR clade (te.marker), or None
        (cached)."""
        if g.name not in self._te:
            self._te[g.name] = self._te_evidence(g)
        return self._te[g.name]

    def _te_verdict(self, g, fams, mains, info):
        """The genome's TE evidence, with the signature's verdict (signature.decisive) made against what the satellite line
        says without the TE line, and written beside the classify checkpoint so a resume reads the same verdict."""
        table = self._p("classify", f"{g.name}.te_signature.tsv")
        table.unlink(missing_ok=True)
        self._te.pop(g.name, None)
        tev = self.te_evidence(g)
        s = self._sig.get(g.name)
        if s is None:
            return tev
        base, _ = classify.choose_centromeric(fams, mains, g.chroms, self.p, info, self.arrays(g))
        test = karyotype.satellite_test(self.arrays(g), set(base.family[base.cls == "centromeric"]), g.chroms, self.p)
        signature.decisive(s, bool(base.dominant_satellite.astype(bool).any()), test["status"])
        signature.write_table(s, table)
        log.info(f"{g.name}: LTR signature: {s.why}")
        self._te[g.name] = s.marker()
        return self._te[g.name]

    def _te_evidence(self, g):
        gff = te.te_gff(self.cdir, g.name)
        if gff is None:
            log.info(f"{g.name}: no tes/{g.name}.gff3 or .gff [.gz]; satellite evidence only")
            return None
        s = signature.build(gff, g.chroms)
        if s is not None:
            self._sig[g.name] = s
            verdict, why = signature.read_verdict(self._p("classify", f"{g.name}.te_signature.tsv"))
            if verdict is not None:                           # made in the classify stage of this folder
                s.decisive, s.why = verdict, why
            return s.marker()
        tes = te.load_te_gff(gff, g.chroms)
        if not tes:
            n = te.te_census(gff, g.chroms)
            log.warning(f"{g.name}: tes/{gff.name}: no LTR retrotransposon with a lineage on these sequences "
                        f"({n['ltr']:,} without one, {n['foreign']:,} on other sequence names); satellite evidence only")
            return None
        cands = te.candidate_clades(tes, g.chroms)
        m = te.choose_marker(tes, g.chroms, cands)
        if m is None:
            log.info(f"{g.name}: no LTR clade concentrated at one locus on most chromosomes; "
                     f"satellite evidence only")
            return None
        crm = sum(1 for x in m[1].values() if "clade" in x)
        log.info(f"{g.name}: LTR clades concentrated at one locus per sequence: "
                 f"{', '.join(f'{c} ({len(x)}/{len(g.chroms)})' for c, x in cands) or 'none'}; marker {m[0]}"
                 f"{' (decisive)' if m[2] else ' (flags only)'}" + (f"; CRM cores on {crm} sequences" if crm else ""))
        return m

    def stage_centromeres(self):
        for g in self.genomes:
            arrays = self.arrays(g)
            fams = self.fams(g, "classify")
            primary = set(fams.family[fams.cls == "centromeric"])
            borrowed = g.name in self.borrow
            tev = None if borrowed else self.te_evidence(g)
            dom = "dominant_satellite" if "dominant_satellite" in fams else "dominant"    # a chain is no satellite
            dominant = fams[dom].astype(str).eq("True").any() if dom in fams else True
            capped = centromere.honesty_cap(tev, borrowed, bool(primary), dominant)

            acro = frozenset(fams.family[fams["acrocentric"].astype(str).eq("True")]) if "acrocentric" in fams \
                else frozenset()
            rivals = classify.acro_rivals(acro, primary, self.fam_info())
            chosen = {}                                      # classify's per-chromosome choice (defines)
            for r in fams[fams.cls == "centromeric"].itertuples(index=False):
                for c in str(r.defines).split(","):
                    if c not in ("", "."):
                        chosen.setdefault(c, set()).add(r.family)

            test = dict(n=len(g.chroms), features=0, covered=0, null_median=float("nan"), p=1.0, best_possible=1.0, status="none") \
                if borrowed else karyotype.satellite_test(arrays, primary, g.chroms, self.p)     # a contig is no chromosome
            atomic_write_text(self._p("centromeres", f"{g.name}.satellite_test.json"), json.dumps(test))
            log.info(f"centromeres: {g.name}: centromeric satellites: {karyotype.sentence(test)}")

            def call_with(prim):
                return centromere.call_regions(g, arrays, prim, self.genes(g), self.seed_table(g),
                                               self.work / "centromeres" / f"tmp_{g.name}", self.threads, self.p,
                                               borrowed=borrowed, te=tev, capped=capped, acro=acro, rivals=rivals,
                                               per_chrom=chosen, not_one=not borrowed and test["status"] == "fail",
                                               weak_cores=None if g.name not in self._sig else
                                               {c: x for c, x in self._sig[g.name].cores.items() if not x["significant"]})
            rows, lu = call_with(primary)
            te_fams = set(fams.family[pd.to_numeric(fams.get("in_ltr"), errors="coerce") >= TE_DERIVED_FRAC]) \
                if "in_ltr" in fams else set()
            for r in rows.values():                          # main family mostly inside LTR retrotransposons
                if r.get("main_family") in te_fams:
                    r["flags"] = "te_derived" if r["flags"] in (".", "") else r["flags"] + ",te_derived"
            if not borrowed:                                 # (a contig's own family metrics mean nothing)
                with pysam.FastaFile(str(g.fasta)) as fa:
                    n = centromere.cap_chain_calls(rows, arrays, primary, classify.chain_families(fams[fams.family.isin(primary)]),
                                                   centromere.derived_families(fams), bool(tev and tev[2]), self.p,
                                                   cores=tev[1] if tev else None, fa=fa, genes=self.genes(g))
                if n:
                    log.info(f"centromeres: {g.name}: {n} calls stand on chains of partial hits, not tandem satellites: "
                             f"at most medium")
            if "te_rejected" in fams:                        # the array each call would otherwise have stood on
                rejected = set(fams.family[fams.te_rejected.astype(str).eq("True")])
                for c, r in rows.items():
                    x = None
                    if rejected and r["evidence"] != "none":
                        x = centromere.rejected_cluster(arrays[arrays.chrom == c], rejected, r["cen_start"], r["cen_end"], self.p)
                    if x:
                        r.update(rejected_family=x[0], rejected_start=x[1], rejected_end=x[2])
            if capped:
                log.info(f"centromeres: {g.name}: no dominant satellite and no decisive retrotransposon evidence; "
                         f"calls capped at low confidence" + (f" unless the {tev[0]} cores support them" if tev else ""))
            fams = classify.assign_cen_associated(fams, arrays, rows, self.p)
            fams = classify.assign_roles(fams, self.p, borrowed)
            n = fams.role.value_counts()
            log.info(f"centromeres: {g.name}: family roles: "
                     + (", ".join(f"{r} {n[r]}" for r in classify.ROLES + ["unassessed"] if r in n) or "none (no families)"))
            centromere.write_regions(rows, self._p("centromeres", f"{g.name}.regions.tsv"))
            atomic_csv(lu, self._p("centromeres", f"{g.name}.longunit.tsv"), sep="\t", index=False)
            atomic_csv(fams, self._p("centromeres", f"{g.name}.families.tsv"), sep="\t", index=False)
            conf = pd.Series([r["confidence"] for r in rows.values()]).value_counts().to_dict()
            log.info(f"centromeres: {g.name}: {conf}; {len(lu)} long-unit arrays")

    def stage_hor(self):
        for g in self.genomes:
            fams = self.fams(g)
            primary = set(fams.family[fams.cls == "centromeric"])
            mh, occ, units = hor.run_hor(g, self.monomers(g), self.arrays(g), self.regions(g), primary,
                                         self.work / "hor" / f"tmp_{g.name}", self.threads, self.p)
            atomic_csv(mh, self._p("hor", f"{g.name}.monomer_hor.tsv.gz"), sep="\t", index=False)
            atomic_csv(occ, self._p("hor", f"{g.name}.hor_occ.tsv"), sep="\t", index=False)
            atomic_csv(units, self._p("hor", f"{g.name}.hor_units.tsv"), sep="\t", index=False)
            log.info(f"hor: {g.name}: {len(units)} HOR units, {len(occ)} occurrences, "
                     f"{(mh.hor_unit != '.').sum()}/{len(mh)} classed monomers in HORs")

    def genome_data(self, g) -> dict:
        rd = lambda n: pd.read_csv(self._p("hor", n), sep="\t", dtype={"chrom": str})
        return dict(genome=g, arrays=self.arrays(g), monomers=self.monomers(g), fams=self.fams(g),
                    regions=self.regions(g),
                    longunit=pd.read_csv(self._p("centromeres", f"{g.name}.longunit.tsv"), sep="\t",
                                         dtype={"chrom": str}),
                    mon_hor=rd(f"{g.name}.monomer_hor.tsv.gz"), hor_occ=rd(f"{g.name}.hor_occ.tsv"),
                    hor_units=rd(f"{g.name}.hor_units.tsv"), seeds=self.seed_table(g), genes=self.genes(g),
                    fam_info=self.fam_info())

    def stage_output(self):
        fam = self.fam_table()
        data = {g.name: self.genome_data(g) for g in self.genomes}
        names = output.final_names(fam, {k: d["arrays"] for k, d in data.items()}, self.name)
        atomic_csv(pd.DataFrame(sorted(names.items()), columns=["fam_id", "name"]), self._p("output", "names.tsv"), sep="\t",
                   index=False)
        for g in self.genomes:
            output.write_genome(data[g.name], names, self.out, self.p)
        output.write_families(fam, data, names, self._p("families", "families.hmm"),
                              self._p("families", "families.fa"), self.out)
        for g in self.genomes:
            for src, name in ((self._p("classify", f"{g.name}.te_signature.tsv"), f"{g.name}_te_signature.tsv"),
                              (self._p("centromeres", f"{g.name}.satellite_test.json"), f"{g.name}_satellite_test.json")):
                if src.exists():
                    atomic_copy(src, self.out / name)
        dropped = self._p("families", "dropped_clusters.tsv")       # absent in a folder of an earlier version
        if dropped.exists():
            atomic_copy(dropped, self.out / "dropped_clusters.tsv")
        for g in self.genomes:
            self._log_costliest(g, data[g.name]["fams"], names)

    def _log_costliest(self, g, fams, names: dict) -> None:
        """Log the searches that cost most CPU (annotate's work table) by final family name and role: a dispersed
        family that cost much decides nothing, which only the roles, known after annotate, can show. A family without
        arrays in this genome is 'absent', as in satellite_families.tsv. A folder made before the table existed logs
        nothing; a table or roles that cannot be read log a warning instead."""
        path = self._p("annotate", f"{g.name}.search_cost.tsv")
        try:
            top = annotate.costliest(pd.read_csv(path, sep="\t", dtype={"family": str})) if path.exists() else []
            role = dict(zip(fams.family, fams.role)) if top else {}
            shown = [f"{names.get(f, f)} ({role.get(f, 'absent')}, {100 * share:.0f}%)" for f, share in top]
        except Exception as e:          # a diagnostic line: it must never stop the output stage before output.done
            log.warning(f"output: {g.name}: costliest searches not logged, {path} or the family roles cannot be read "
                        f"({type(e).__name__}: {e})")
            return
        if shown:
            log.info(f"output: {g.name}: costliest searches: " + ", ".join(shown))

    def stage_plot(self):
        names = dict(pd.read_csv(self._p("output", "names.tsv"), sep="\t").itertuples(index=False))
        panels = []
        for g in self.genomes:
            fams = self.fams(g)
            a = self.arrays(g)
            a = a.assign(cls=a.family.map(dict(zip(fams.family, fams.cls))).fillna("other"),
                         family=a.family.map(names))
            panels.append(dict(genome=g.name, chroms=g.chroms, arrays=a, regions=self.regions(g),
                               genes=self.genes(g), borrowed=g.name in self.borrow))
        plot.plot_overview(self.out / "centromere_overview.pdf", panels)
