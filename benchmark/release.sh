#!/usr/bin/env bash
# release.sh: score a release's stored runs, write every table docs/benchmark.md uses, and fill its template
# (benchmark/README.md).
#
#   bash benchmark/release.sh OUT [TAG]
#
# OUT (created) receives, TAG (default 16) in the names:
#   release$TAG.{chrom,runs}.tsv .txt .err   the scores (score.py)    cmp$TAG.{diffs,summary}.tsv .txt  vs reference runs
#   rows_release$TAG.txt   rows beside PREV's   md$TAG/ md$TAG.log   the doc tables (+ one_per_chromosome.tsv, cis_trans.tsv)
#   lines_alone$TAG.*  k2_release$TAG.*  randomized_release$TAG.*  headroom$TAG.*  precision$TAG.tsv (+ md$TAG/TABLE_PRECISION.md)
#   dev$TAG.vs_prev.*  dev$TAG.vs_extra.*    calls changed vs PREV_RUNS and vs EXTRA_REF (score_dev.py)
#   docs/docs/benchmark.md   the filled template (DOCS=<repo root> rewrites the repo's own)
# Inputs default to the 0.16.0 release (python -m benchmark.paths); set any to override, even to empty: RUNS (release
# root, RUNS/<arm>/<id>), JOBS, JOBS_CMP, JOBS_RABL, PREV (previous release's scores, a prefix), PREV_RUNS, PREV_REDO
# (empty: none), EXTRA_REF (empty: skipped), CIS_SURVEY, DOCS,
# CENFUSE_BENCH_DATA (data root), PY (a python with numpy and pandas; default python). Stops if that python imports
# cenfuse from anywhere but this repository.
set -uo pipefail
[ $# -ge 1 ] || { awk 'NR > 1 && !/^#/ {exit} NR > 1' "$0" >&2; exit 2; }     # usage: the comment above
OUT=$(mkdir -p "$1" && cd "$1" && pwd -P) || exit 2
TAG=${2:-16}
REPO=$(cd "$(dirname "$0")/.." && pwd -P)
PY=${PY:-python}
export PYTHONPATH="$REPO${PYTHONPATH:+:$PYTHONPATH}" OMP_NUM_THREADS=1 NUMPY_MADVISE_HUGEPAGE=0
got=$("$PY" -c 'import os, cenfuse; print(os.path.realpath(os.path.dirname(os.path.dirname(cenfuse.__file__))))') || exit 2
[ "$got" = "$REPO" ] || { echo "release.sh: $PY imports cenfuse from $got, not $REPO" >&2; exit 2; }
eval "$("$PY" -m benchmark.paths --shell)" || exit 2
DOCS=${DOCS-$OUT/docs}
fail=()
run() { "$PY" -m "benchmark.$1" "${@:2}"; }
step() {   # step MODULE ARGS...: one table; a failure is reported and the chain goes on
    run "$@" || { echo "release.sh: $1 failed" >&2; fail+=("$1"); }
}
echo "release.sh: $RUNS -> $OUT (tag $TAG)"
run score "$JOBS" "$RUNS" "$OUT/release$TAG" > "$OUT/release$TAG.txt" 2> "$OUT/release$TAG.err" \
    || { echo "release.sh: scoring failed (see $OUT/release$TAG.err)" >&2; exit 1; }
S=$OUT/release$TAG M=$OUT/md$TAG
step compare_rel "$JOBS_CMP" "$RUNS" "$OUT/cmp$TAG" > "$OUT/cmp$TAG.txt" 2>&1
step rows "$S" "$PREV" > "$OUT/rows_release$TAG.txt" 2>&1
step md_tables "$S" "$M" --runs "$RUNS" > "$M.log" 2>&1
step md_013 "$M" --scores "$S" --runs "$RUNS" --jobs "$JOBS" --cis-survey "$CIS_SURVEY" >> "$M.log" 2>&1
step rabl_table "$M/TABLE_RABL.md" --runs "$RUNS" --jobs "$JOBS_RABL" >> "$M.log" 2>&1
step precision "$S.chrom.tsv" "$M/TABLE_PRECISION.md" --tsv "$OUT/precision$TAG.tsv" >> "$M.log" 2>&1
step lines_alone "$OUT/lines_alone$TAG.tsv" --runs "$RUNS" > "$OUT/lines_alone$TAG.txt" 2>&1
step score_k2_rel "$OUT/k2_release$TAG.tsv" --runs "$RUNS" > "$OUT/k2_release$TAG.txt" 2>&1
step randomized "$JOBS" "$RUNS" "$OUT/randomized_release$TAG.tsv" > "$OUT/randomized_release$TAG.txt" 2>&1
step score_dev "$RUNS" "$JOBS" "$OUT/dev$TAG.vs_prev" --ref "$PREV_RUNS" ${PREV_REDO:+--redo "$PREV_REDO"} \
    > "$OUT/dev$TAG.vs_prev.txt" 2>&1
[ -z "$EXTRA_REF" ] || step score_dev "$RUNS" "$JOBS" "$OUT/dev$TAG.vs_extra" --ref "$EXTRA_REF" \
    > "$OUT/dev$TAG.vs_extra.txt" 2>&1
step headroom "$RUNS" "$JOBS" "$OUT/headroom$TAG.tsv" > "$OUT/headroom$TAG.txt" 2>&1
step fill_docs "$M" "$DOCS"
if [ ${#fail[@]} -gt 0 ]; then
    echo "release.sh: FAILED: ${fail[*]} (logs in $OUT)" >&2
    exit 1
fi
echo "FINISHED $(date -Is)"
