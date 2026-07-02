#!/usr/bin/env bash
# Paired KLS/CKTSO benchmark: run both solvers back-to-back per matrix so
# slow machine-state drift cancels inside each pair.
set -u
KLS_BENCH=${1:?kls_bench path}
CKTSO_COMPARE=${2:?cktso_compare path}
MATRIX_DIR=${3:?matrix dir}
MANIFEST=${4:?manifest}
OUT_KLS=${5:?kls jsonl}
OUT_CK=${6:?cktso jsonl}
THREADS=${THREADS:-4}
TIMEOUT=${TIMEOUT:-120}
: > "$OUT_KLS"
: > "$OUT_CK"
while IFS= read -r name; do
  case "$name" in ''|'#'*) continue;; esac
  matrix=$(find "$MATRIX_DIR" -name "${name}.mtx" | head -1)
  [ -z "$matrix" ] && { echo "skip $name (not found)" >&2; continue; }
  for side in kls ck; do
    if [ "$side" = kls ]; then
      out=$(timeout "$TIMEOUT" "$KLS_BENCH" "$matrix" --orientation auto \
        --repeat 1 --refactor-repeat 3 --threads "$THREADS" --json 2>/dev/null)
      dst=$OUT_KLS
    else
      out=$(timeout "$TIMEOUT" "$CKTSO_COMPARE" "$matrix" "$THREADS" 1 3 \
        2>/dev/null)
      dst=$OUT_CK
    fi
    if [ -n "$out" ]; then
      printf '%s' "$out" | python3 -c "
import json,sys
d=json.load(sys.stdin)
d['matrix']='$matrix'
print(json.dumps(d))" >> "$dst" 2>/dev/null || echo "{\"matrix\":\"$matrix\",\"status\":\"parse_error\"}" >> "$dst"
    else
      echo "{\"matrix\":\"$matrix\",\"status\":\"timeout_or_fail\"}" >> "$dst"
    fi
  done
  echo "done $name" >&2
done < "$MANIFEST"
