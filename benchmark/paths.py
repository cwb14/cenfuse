"""Every root the benchmark reads, in one place.

The truth files, the stored runs and the earlier releases' tables live outside the repository (too large, and some truth
is unpublished). DATA is the root they live under: $CENFUSE_BENCH_DATA, else the folder they were made in (DEFAULT_DATA).
The names below default to the 0.16.1 release, whose runs and tables made docs/benchmark.md; release.sh takes any of
them from the environment instead.

  python -m benchmark.paths [--shell]     # print the roots (--shell: as bash defaults for release.sh)"""
import argparse
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
MANIFEST = HERE / "truth.tsv"            # one row per genome: kind, rule, tolerance, truth files (relative to DATA)
CANONICAL = HERE / "canonical.tsv"       # genome labels, lineages and the arm each doc table shows
DOCS_SRC = HERE / "docs_src"             # the template of docs/benchmark.md

DEFAULT_DATA = Path("/data2/chris/centromere_benchmark")
DATA = Path(os.environ.get("CENFUSE_BENCH_DATA") or DEFAULT_DATA)

# development folders, nested: v2 = 0.1.0, v11 = 0.12, v12 = 0.13, v13 = 0.14, v14 = 0.15, v16 = 0.16
V2 = DATA / "v2"
V5 = V2 / "v3/v4/v5"
V11 = V5 / "v6/v7/v8/v9/v10/v10.1/v10.2/v10.3/v10.4/v11"
V12 = V11 / "v12"
V13 = V12 / "v13"
V14 = V13 / "v14"
V16 = V14 / "v15/v16"

# the 0.16.1 release: the 0.16.0 runs (runs/rel16), the 51 that 0.16.1 can change rerun (runs/rel161); runs/rel161m
# links each arm to its newest run
RUNS = V16 / "runs/rel161m"                      # release root: RUNS/<arm>/<id>
JOBS = V13 / "bench/jobs_rel14.tsv"              # arm, id, inputs, threads, satellite dir, reference run (260 arms)
RESULTS = V16 / "results/docs161_rel161m"        # the chain's stored outputs (read-only reference)
JOBS_CMP = RESULTS / "jobs_cmp16.tsv"            # the same arms; column 6 = the 0.15.0 run of each
JOBS_RABL = V12 / "bench/jobs_rel12_run.tsv"     # the 220 arms of 0.12's own run (the Rabl-axis table)
SCORES = RESULTS / "release16"                   # this release's scores (score.py prefix)
TABLES = RESULTS / "md16"                        # its computed doc tables (fill_docs.py TABLES)
PREV = V14 / "results/docs15_rel15c/release15"   # 0.15.0's scores (score.py prefix)
PREV_RUNS = V14 / "runs/rel15c"                  # 0.15.0's runs (score_dev.py reference)
PREV_REDO = None                                 # reference runs redone after the release (none for 0.15.0)
EXTRA_REF = None                                 # a second reference for score_dev.py (none)
CIS_SURVEY = V11 / "results/cis_survey4.tsv"     # the cis-signal survey of every Hi-C genome (0.12)
ASM_REPORTS = V2 / "data/truth/src"              # NCBI assembly reports: accession -> chromosome name (K2 check)
MOUSE_NAMES = V5 / "data/mouse/truth/name_map.tsv"   # mouse accession -> chromosome name (K2 check)

SHELL = {"RUNS": RUNS, "JOBS": JOBS, "JOBS_CMP": JOBS_CMP, "JOBS_RABL": JOBS_RABL, "PREV": PREV,
         "PREV_RUNS": PREV_RUNS, "PREV_REDO": PREV_REDO, "EXTRA_REF": EXTRA_REF, "CIS_SURVEY": CIS_SURVEY}


def local(p) -> Path:
    """A path a job table wrote under the default root, moved under DATA (job tables hold absolute paths)."""
    p = Path(p)
    try:
        return DATA / p.relative_to(DEFAULT_DATA)
    except ValueError:
        return p


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shell", action="store_true", help="print bash defaults (VAR=${VAR-path}) for release.sh")
    a = ap.parse_args(argv)
    if a.shell:
        for k, v in SHELL.items():
            print(f'{k}="${{{k}-{v or ""}}}"')   # an unset name takes the default; a set one, even empty, is kept
        return
    print(f"DATA\t{DATA}")
    for k, v in {"MANIFEST": MANIFEST, "CANONICAL": CANONICAL, "DOCS_SRC": DOCS_SRC, **SHELL, "SCORES": SCORES,
                 "TABLES": TABLES, "ASM_REPORTS": ASM_REPORTS, "MOUSE_NAMES": MOUSE_NAMES}.items():
        print(f"{k}\t{v or ''}")


if __name__ == "__main__":
    main()
