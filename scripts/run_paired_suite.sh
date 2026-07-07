#!/usr/bin/env bash
# Paired multi-solver benchmark: run every solver back-to-back per matrix so
# slow machine-state drift cancels inside each group.
set -u
KLS_BENCH=${1:?kls_bench path}
CKTSO_COMPARE=${2:?cktso_compare path}
MATRIX_DIR=${3:?matrix dir}
MANIFEST=${4:?manifest}
OUT_KLS=${5:?kls jsonl}
OUT_CK=${6:?cktso jsonl}
SUBTREELU_COMPARE=${7:-}
OUT_ST=${8:-}
THREADS=${THREADS:-4}
TIMEOUT=${TIMEOUT:-120}
PASSES=${PASSES:-1}
: > "$OUT_KLS"
: > "$OUT_CK"
[ -n "$OUT_ST" ] && : > "$OUT_ST"
sides="kls ck"
[ -n "$SUBTREELU_COMPARE" ] && sides="kls ck st"
passes_sides=""
p=0
while [ "$p" -lt "$PASSES" ]; do
  passes_sides="$passes_sides $sides"
  p=$((p + 1))
done
while IFS= read -r name; do
  case "$name" in ''|'#'*) continue;; esac
  matrix=$(find "$MATRIX_DIR" -iname "${name}.mtx" | head -1)
  [ -z "$matrix" ] && { echo "skip $name (not found)" >&2; continue; }
  for side in $passes_sides; do
    if [ "$side" = kls ]; then
      out=$(timeout "$TIMEOUT" "$KLS_BENCH" "$matrix" --orientation auto \
        --repeat 1 --refactor-repeat 10 --threads "$THREADS" --json 2>/dev/null)
      dst=$OUT_KLS
    elif [ "$side" = ck ]; then
      out=$(timeout "$TIMEOUT" "$CKTSO_COMPARE" "$matrix" "$THREADS" 1 10 \
        2>/dev/null)
      dst=$OUT_CK
    else
      out=$(timeout "$TIMEOUT" "$SUBTREELU_COMPARE" "$matrix" "$THREADS" 1 10 \
        2>/dev/null)
      dst=$OUT_ST
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
