"""Command line: cenfuse run GENOME_DIR [-o OUT] [-t THREADS] [--satellite-dir DIR] [--seed-engine {fastan,trf}]
[--seed-cache DIR] [--classify-from T=S] [--force STAGE] [-v]. Every threshold is an internal constant
(cenfuse/params.py and the engines' params); none is an option."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .util import CenfuseError, log, setup_logging


def main(argv=None) -> int:
    from .pipeline import STAGES
    from .sat.pipeline import STAGES as SAT_STAGES
    ap = argparse.ArgumentParser(prog="cenfuse", description=(
        "Centromeres from satellite arrays, LTR retrotransposons, Hi-C and methylation: each line's own call, "
        "their agreement, one integrated call per chromosome, and the other features they find, on one map."))
    ap.add_argument("--version", action="version", version=f"cenfuse {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run every evidence line present and integrate them (resumes)",
                       description="Run every evidence line present in GENOME_DIR and integrate them. A rerun "
                                   "resumes from checkpoints.")
    r.add_argument("genome_dir", type=Path, metavar="GENOME_DIR",
                   help="one FASTA per genome (+ .fai; bgzip also .gzi); optional genes/<g>.gff3, tes/<g>.gff3 "
                        "(LTR_retrotransposon with clade=), hic/<g>.pairs[.gz]|.cool|.mcool or hic/<g>.tsv (one "
                        "library per plant: sample, file), meth/<g>.tsv (samples table: sample, CHG[, CG, CHH] "
                        "modkit bedMethyl) [+ meth/<g>.exclude.bed]. A Hi-C library and a methylation sample with "
                        "the same name are the same plant")
    r.add_argument("-o", "--outdir", type=Path, help="output dir [GENOME_DIR/cenfuse]")
    r.add_argument("-t", "--threads", type=int, default=8, help="threads [8]")
    r.add_argument("--satellite-dir", type=Path, metavar="DIR",
                   help="reuse the satellite outputs of an earlier run (<g>_centromere_summary.tsv, <g>_satellite_arrays.gff) "
                        "instead of running the satellite engine")
    r.add_argument("--seed-engine", choices=["fastan", "trf"], default="fastan",
                   help="tandem-repeat finder for satellite family seeds [fastan]")
    r.add_argument("--seed-cache", type=Path, metavar="DIR", help="seed cache shareable between runs")
    r.add_argument("--classify-from", default="", metavar="T=S[+S2],...",
                   help="genome T (e.g. a contig assembly) takes the centromeric satellite families of genome(s) S")
    r.add_argument("--force", choices=STAGES + SAT_STAGES, metavar="STAGE",
                   help=f"recompute STAGE ({', '.join(STAGES)}; or a satellite-engine stage: "
                        f"{', '.join(SAT_STAGES)})")
    r.add_argument("-v", "--verbose", action="store_true", help="per-step progress and sanity checks")
    a = ap.parse_args(argv)
    setup_logging(a.verbose)
    from .pipeline import Run
    try:
        Run(a.genome_dir, a.outdir, a.threads, a.satellite_dir, a.seed_engine, a.seed_cache, a.classify_from,
            a.force).run()
    except CenfuseError as e:
        print(f"cenfuse: error: {e}", file=sys.stderr)
        return 2
    return 0
