"""Command line: cenfuse run GENOME_DIR [-o OUT] [-t THREADS] [--satellite-dir DIR] [--seed-engine {fastan,trf}]
[--seed-cache DIR] [--classify-from T=S] [--force STAGE] [-v]. Every threshold is an internal constant
(cenfuse/params.py and the engines' params); none is an option."""
from __future__ import annotations

import argparse
import os
import signal
import sys
from pathlib import Path

from . import __version__
from .util import CenfuseError, setup_logging


def _threads(x: str) -> int:
    try:
        n = int(x)
    except ValueError:
        raise argparse.ArgumentTypeError(f"must be a whole number, not {x!r}") from None
    if n < 1:
        raise argparse.ArgumentTypeError(f"must be at least 1, not {n}")
    return n


_MAIN = os.getpid()


def _terminate(signum, frame):
    """SIGTERM (a scheduler's kill, `kill PID`) ends the run as an exception, so the output folder's lock is released and
    run_info.json records that the run failed. A worker process forked from the run inherits this handler; it dies as SIGTERM
    would, at once, so a pool does not finish its queued work first."""
    if os.getpid() != _MAIN:
        signal.signal(signum, signal.SIG_DFL)
        os.kill(os.getpid(), signum)
        return
    raise SystemExit(128 + signum)


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv[:1] == ["control"]:                     # negative controls for one input (cenfuse.controls)
        from .controls import main as control
        return control(argv[1:])
    from .pipeline import STAGES
    from .sat.pipeline import STAGES as SAT_STAGES
    ap = argparse.ArgumentParser(prog="cenfuse", description=(
        "Centromeres from satellite arrays, LTR retrotransposons, Hi-C and methylation: each line's own call, "
        "their agreement, one integrated call per chromosome, and the other features they find, on one map."))
    ap.add_argument("--version", action="version", version=f"cenfuse {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("control", add_help=False, help="write a negative control of one input (randomized Hi-C pairs, shuffled GC "
                                                     "blocks or LTR positions): cenfuse control -h")
    r = sub.add_parser("run", help="run every evidence line present and integrate them (resumes)",
                       description="Run every evidence line present in GENOME_DIR and integrate them. A rerun "
                                   "resumes from checkpoints.")
    r.add_argument("genome_dir", type=Path, metavar="GENOME_DIR",
                   help="one FASTA per genome (+ .fai; bgzip also .gzi) [+ <g>.chromosomes.txt: the chromosomes to "
                        "analyze, one name per line; default: a length rule]; optional genes/<g>.gff3|.gff|.gtf, "
                        "tes/<g>.gff3|.gff (LTR_retrotransposon with clade=, or DANTE_LTR; make it with the LTRquest "
                        "command in the README), "
                        "hic/<g>.pairs[.gz]|.cool|.mcool (or HiC-Pro .allValidPairs[.gz], Juicer "
                        ".merged_nodups.txt[.gz]) or hic/<g>.tsv (one "
                        "library per plant: sample, file), meth/<g>.tsv (samples table: sample, CHG and/or CG, "
                        "optional CHH: modkit bedMethyl; CG also ENCODE bedMethyl or a bedGraph) "
                        "[+ meth/<g>.exclude.bed]. A Hi-C library and a methylation sample with the same name are "
                        "the same plant")
    r.add_argument("-o", "--outdir", type=Path, help="output dir [GENOME_DIR/cenfuse]")
    r.add_argument("-t", "--threads", type=_threads, default=8, help="threads [8]")
    r.add_argument("--satellite-dir", type=Path, metavar="DIR",
                   help="reuse the satellite outputs of an earlier run (<g>_centromere_summary.tsv, <g>_satellite_arrays.gff) "
                        "instead of running the satellite engine")
    r.add_argument("--seed-engine", choices=["fastan", "trf"], default="fastan",
                   help="tandem-repeat finder for satellite family seeds [fastan]")
    r.add_argument("--seed-cache", type=Path, metavar="DIR", help="seed cache shareable between runs")
    r.add_argument("--classify-from", default="", metavar="T=S[+S2],...",
                   help="genome T (e.g. a contig assembly) takes the centromeric satellite families of genome(s) S")
    r.add_argument("--force", choices=["all"] + STAGES + SAT_STAGES, metavar="STAGE",
                   help=f"recompute STAGE, or every stage (all): {', '.join(STAGES)}; a satellite-engine stage "
                        f"({', '.join(SAT_STAGES)}) recomputes the engine from there. A later stage reruns when what it reads "
                        f"changed, as on any rerun")
    r.add_argument("-v", "--verbose", action="store_true", help="per-step progress and sanity checks")
    a = ap.parse_args(argv)
    setup_logging(a.verbose)
    from .pipeline import Run
    from .sat.util import SatelliteError
    global _MAIN
    _MAIN = os.getpid()
    signal.signal(signal.SIGTERM, _terminate)
    try:
        Run(a.genome_dir, a.outdir, a.threads, a.satellite_dir, a.seed_engine, a.seed_cache, a.classify_from,
            a.force).run()
    except (CenfuseError, SatelliteError) as e:      # bad input or a failed precondition, also from an engine's reader
        print(f"cenfuse: error: {e}", file=sys.stderr)
        return 2
    return 0
