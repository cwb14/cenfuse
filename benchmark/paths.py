"""Every root the benchmark reads, in one place.

The truth files, the stored runs and the earlier releases' tables live outside the repository (some truth is
unpublished). DATA is the root they live under: /data2/chris/centromere_benchmark, or $CENFUSE_BENCH_DATA. The names
below default to the 0.13.0 release, whose runs and tables made README.md and docs/benchmark.md.

  python -m benchmark.paths [--shell]     # print the roots (--shell: as bash defaults for release.sh)"""
import argparse
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
MANIFEST = HERE / "truth.tsv"            # one row per genome: kind, rule, tolerance, truth files (relative to DATA)
CANONICAL = HERE / "canonical.tsv"       # genome labels, lineages and the arm each doc table shows (from v11/results)
DOCS_SRC = HERE / "docs_src"             # templates of README.md and docs/benchmark.md, and their hand-written texts

DEFAULT_DATA = Path("/data2/chris/centromere_benchmark")
DATA = Path(os.environ.get("CENFUSE_BENCH_DATA") or DEFAULT_DATA)

# development folders (future.md, "Where things live"): v2 = 0.1.0, v11 = 0.12, v12 = 0.13
V2 = DATA / "v2"
V5 = V2 / "v3/v4/v5"
V11 = V5 / "v6/v7/v8/v9/v10/v10.1/v10.2/v10.3/v10.4/v11"
V12 = V11 / "v12"

# the 0.13.0 release
RUNS = V12 / "runs/rel13f"                       # release root: RUNS/<arm>/<id>
JOBS = V12 / "bench/jobs_rel13.tsv"              # arm, id, inputs, threads, satellite dir, reference run (260 arms)
JOBS_CMP = V12 / "bench/jobs_rel13_cmp.tsv"      # the same arms; column 6 = the 0.12 run of each
JOBS_RABL = V12 / "bench/jobs_rel12_run.tsv"     # the 220 arms of 0.12's own run (the Rabl-axis table)
RESULTS = V12 / "results"                        # the chain's stored outputs (read-only reference)
SCORES = RESULTS / "release13"                   # this release's scores (score.py prefix)
PREV = RESULTS / "release12all"                  # 0.12's scores, held-out arms included (score.py prefix)
PREV_RUNS = V11 / "runs/rel12"                   # 0.12's runs (score_dev.py reference)
PREV_REDO = V11 / "runs/rel12b"                  # the 20 runs 0.12 redid after its held-out change
EXTRA_REF = V12 / "runs/rel13"                   # 0.13's first release pass (score_dev.py, a second reference)
CIS_SURVEY = V11 / "results/cis_survey4.tsv"     # the cis-signal survey of every Hi-C genome (0.12 bench/cis_survey4.py)
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
            print(f'{k}="${{{k}-{v}}}"')       # an unset name takes the default; a set one, even empty, is kept
        return
    print(f"DATA\t{DATA}")
    for k, v in {"MANIFEST": MANIFEST, "CANONICAL": CANONICAL, "DOCS_SRC": DOCS_SRC, **SHELL, "SCORES": SCORES,
                 "ASM_REPORTS": ASM_REPORTS, "MOUSE_NAMES": MOUSE_NAMES}.items():
        print(f"{k}\t{v}")


if __name__ == "__main__":
    main()
