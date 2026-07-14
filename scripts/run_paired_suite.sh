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
KLU_COMPARE=${9:-}
OUT_KLU=${10:-}
THREADS=${THREADS:-4}
KLS_BACKEND=${KLS_BACKEND:-auto}
TIMEOUT=${TIMEOUT:-120}
PASSES=${PASSES:-1}
ROTATE_SIDES=${ROTATE_SIDES:-0}
# 5 refactors suffice: the harnesses report the first refactor (which
# carries KLS's deferred engine preps) separately from the steady-state
# average, so the 99-iteration SPICE cycle is reconstructed exactly as
# shot + (first + solve) + 98*(steady + solve).
REFACTOR_REPEAT=${REFACTOR_REPEAT:-5}
# The measured cycle starts at the initial factorization.  Extra full-factor
# timing passes would mutate numeric state before the first refactor and are
# therefore disabled for this paired experiment.
FACTOR_REPEAT=${FACTOR_REPEAT:-0}
: > "$OUT_KLS"
: > "$OUT_CK"
[ -n "$OUT_ST" ] && : > "$OUT_ST"
[ -n "$OUT_KLU" ] && : > "$OUT_KLU"
sides="kls ck"
[ -n "$SUBTREELU_COMPARE" ] && sides="kls ck st"
[ -n "$KLU_COMPARE" ] && sides="$sides klu"
passes_sides=""
p=0
while [ "$p" -lt "$PASSES" ]; do
  pass_sides=$sides
  if [ "$ROTATE_SIDES" -ne 0 ]; then
    shift_count=$((p % 4))
    while [ "$shift_count" -gt 0 ]; do
      first=${pass_sides%% *}
      rest=${pass_sides#* }
      [ "$rest" = "$pass_sides" ] && break
      pass_sides="$rest $first"
      shift_count=$((shift_count - 1))
    done
  fi
  passes_sides="$passes_sides $pass_sides"
  p=$((p + 1))
done
while IFS= read -r name; do
  case "$name" in ''|'#'*) continue;; esac
  matrix=$(find "$MATRIX_DIR" -iname "${name}.mtx" | head -1)
  [ -z "$matrix" ] && { echo "skip $name (not found)" >&2; continue; }
  failed_sides=""
  for side in $passes_sides; do
    case " $failed_sides " in
      *" $side "*) continue;;
    esac
    if [ "$side" = kls ]; then
      out=$(timeout "$TIMEOUT" "$KLS_BENCH" "$matrix" --orientation auto \
        --repeat 1 --factor-repeat "$FACTOR_REPEAT" \
        --refactor-repeat "$REFACTOR_REPEAT" --threads "$THREADS" \
        --backend "$KLS_BACKEND" --json 2>/dev/null)
      dst=$OUT_KLS
    elif [ "$side" = ck ]; then
      out=$(timeout "$TIMEOUT" "$CKTSO_COMPARE" "$matrix" "$THREADS" 1 \
        "$REFACTOR_REPEAT" "$FACTOR_REPEAT" \
        2>/dev/null)
      dst=$OUT_CK
    elif [ "$side" = st ]; then
      out=$(timeout "$TIMEOUT" "$SUBTREELU_COMPARE" "$matrix" "$THREADS" 1 \
        "$REFACTOR_REPEAT" "$FACTOR_REPEAT" \
        2>/dev/null)
      dst=$OUT_ST
    else
      out=$(timeout "$TIMEOUT" "$KLU_COMPARE" "$matrix" --repeat 1 \
        --factor-repeat "$FACTOR_REPEAT" --refactor-repeat "$REFACTOR_REPEAT" \
        --json 2>/dev/null)
      dst=$OUT_KLU
    fi
    if [ -n "$out" ]; then
      printf '%s' "$out" | python3 -c "
import json,sys
d=json.load(sys.stdin)
d['matrix']='$matrix'
print(json.dumps(d))" >> "$dst" 2>/dev/null || echo "{\"matrix\":\"$matrix\",\"status\":\"parse_error\"}" >> "$dst"
    else
      echo "{\"matrix\":\"$matrix\",\"status\":\"timeout_or_fail\"}" >> "$dst"
      # A timeout or solver failure is deterministic for this fixed matrix
      # and configuration; do not spend the remaining passes repeating it.
      failed_sides="$failed_sides $side"
    fi
  done
  echo "done $name" >&2
done < "$MANIFEST"
