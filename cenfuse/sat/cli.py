"""Command line: python -m cenfuse.sat run GENOME_DIR [-o OUTDIR] [-t THREADS] [--seed-engine {fastan,trf}] [--seed-cache DIR]
[--classify-from T=S] [--force STAGE] [-v]. Every threshold is an internal constant (params.py); none is an option."""
import argparse
from pathlib import Path

from ..cli import _threads
from ..inputs import chromosome_set
from . import __version__
from .params import Params
from .util import SatelliteError, log, setup_logging


def _warn_short(g) -> None:
    """Name the sequences cenfuse's chromosome rule leaves out (organelles, unplaced scaffolds): the one-per-chromosome test
    and family prevalence count every sequence as a chromosome. `cenfuse run` hands the engine its chromosome set."""
    keep = chromosome_set(g.chroms)
    short = [c for c in g.chroms if c not in keep]
    if short:
        names = ", ".join(repr(c) for c in short[:5]) + (f" and {len(short) - 5} more" if len(short) > 5 else "")
        log.warning(f"{g.name}: {len(short)} of {len(g.chroms)} sequences look too short to be chromosomes ({names}); "
                    f"the one-per-chromosome test counts them as chromosomes, which can cap every satellite call at low: "
                    f"give a FASTA of the chromosomes only")


def main(argv=None) -> int:
    from .pipeline import STAGES, Pipeline
    from .seeds import ENGINES
    ap = argparse.ArgumentParser(prog="python -m cenfuse.sat", description=(
        "Centromeric satellite, HOR and tandem-repeat annotation. Finds satellite families de novo, shared across "
        "all genomes in GENOME_DIR, and calls one centromere per chromosome with its evidence and confidence."))
    ap.add_argument("--version", action="version", version=f"cenfuse-sat {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run all stages; resumes from checkpoints in OUTDIR/work/",
                       description="Run all stages; a rerun resumes from checkpoints in OUTDIR/work/.")
    r.add_argument("genome_dir", nargs="?", type=Path, metavar="GENOME_DIR",
                   help="dir with one FASTA per genome (*.fa/*.fasta/*.fna, plain or bgzip, indexed with "
                        "samtools faidx); optional GENOME_DIR/genes/<genome>.gff3|.gff|.gtf[.gz] (gene loci) and "
                        "GENOME_DIR/tes/<genome>.gff3|.gff[.gz] (LTR retrotransposons with clade=, or DANTE_LTR)")
    r.add_argument("--complex-dir", dest="complex_dir", type=Path, help=argparse.SUPPRESS)   # alias of GENOME_DIR
    r.add_argument("-o", "--outdir", type=Path, help="output dir [GENOME_DIR/centromere]")
    r.add_argument("-t", "--threads", type=_threads, default=8, help="threads [8]")
    r.add_argument("--seed-engine", choices=list(ENGINES), default=Params.seed_engine,
                   help="tandem-repeat finder for the family seeds [fastan]")
    r.add_argument("--seed-cache", type=Path, metavar="DIR",
                   help="seed result cache, shareable between runs [OUTDIR/work/seed_cache, or trf_cache if present]")
    r.add_argument("--trf-cache", dest="seed_cache", type=Path, help=argparse.SUPPRESS)    # alias of --seed-cache
    r.add_argument("--classify-from", default="", metavar="T=S[+S2],...",
                   help="genome T (e.g. a contig assembly) takes the centromeric families of genome(s) S")
    r.add_argument("--force", choices=STAGES, help="recompute this stage and all later ones")
    r.add_argument("-v", "--verbose", action="store_true", help="per-step progress and sanity checks")
    a = ap.parse_args(argv)
    gdir = a.genome_dir or a.complex_dir
    if gdir is None or (a.genome_dir and a.complex_dir and a.genome_dir != a.complex_dir):
        ap.error("give one GENOME_DIR")
    setup_logging(a.verbose)
    try:
        pipe = Pipeline(gdir, a.threads, a.seed_cache, Params(seed_engine=a.seed_engine, classify_from=a.classify_from),
                        force=a.force, out_dir=a.outdir)
        for g in pipe.genomes:
            if g.name not in pipe.borrow:
                _warn_short(g)
        pipe.run()
    except SatelliteError as e:
        log.error(str(e))
        return 1
    return 0
