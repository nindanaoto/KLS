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
# Optional explicit subset/order for pairwise audits.  The default preserves
# the historical behavior inferred from the supplied harness paths.
PAIRED_SIDES=${PAIRED_SIDES:-}
SKIP=${SKIP:-0}
SOLVE_REPEAT=${SOLVE_REPEAT:-1}
INPUT_INDEX=${INPUT_INDEX:-64}
# The harnesses report the first refactor separately from the remaining
# average.  Five refactors are useful for a quick screen; use a larger value
# when adaptive engine trials must be diluted into the long-run steady rate.
REFACTOR_REPEAT=${REFACTOR_REPEAT:-5}
# Keep the released-demo-compatible unchanged loop as the default.  Set
# REFACTOR_VALUES=rank-preserving for the separable nonlinear-SPICE check, or
# entrywise for a deterministic certificate-rejection counter-workload, or
# localized-entrywise for a fixed one-in-1024 cyclic column window.
REFACTOR_VALUES=${REFACTOR_VALUES:-unchanged}
REFACTOR_VALUE_AMPLITUDE=${REFACTOR_VALUE_AMPLITUDE:-0.001}
case "$REFACTOR_VALUES" in
  unchanged|rank-preserving|entrywise|localized-entrywise) ;;
  *) echo "invalid REFACTOR_VALUES: $REFACTOR_VALUES" >&2; exit 2;;
esac
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
if [ -n "$PAIRED_SIDES" ]; then
  selected_sides=""
  for side in $PAIRED_SIDES; do
    case "$side" in
      kls) ;;
      ck) ;;
      st)
        [ -n "$SUBTREELU_COMPARE" ] && [ -n "$OUT_ST" ] || {
          echo "PAIRED_SIDES requests st without a SubtreeLU harness/output" >&2
          exit 2
        }
        ;;
      klu)
        [ -n "$KLU_COMPARE" ] && [ -n "$OUT_KLU" ] || {
          echo "PAIRED_SIDES requests klu without a KLU harness/output" >&2
          exit 2
        }
        ;;
      *) echo "invalid PAIRED_SIDES member: $side" >&2; exit 2;;
    esac
    case " $selected_sides " in
      *" $side "*) echo "duplicate PAIRED_SIDES member: $side" >&2; exit 2;;
    esac
    selected_sides="${selected_sides:+$selected_sides }$side"
  done
  case " $selected_sides " in
    *" kls "*) ;;
    *) echo "PAIRED_SIDES must include kls" >&2; exit 2;;
  esac
  [ "${selected_sides#* }" != "$selected_sides" ] || {
    echo "PAIRED_SIDES must select at least two solvers" >&2
    exit 2
  }
  sides=$selected_sides
fi
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
selected_index=0
while IFS= read -r name; do
  case "$name" in ''|'#'*) continue;; esac
  if [ "$selected_index" -lt "$SKIP" ]; then
    selected_index=$((selected_index + 1))
    continue
  fi
  selected_index=$((selected_index + 1))
  matrix=$(find "$MATRIX_DIR" -iname "${name}.mtx" | head -1)
  [ -z "$matrix" ] && { echo "skip $name (not found)" >&2; continue; }
  failed_sides=""
  for side in $passes_sides; do
    case " $failed_sides " in
      *" $side "*) continue;;
    esac
    if [ "$side" = kls ]; then
      out=$(timeout "$TIMEOUT" "$KLS_BENCH" "$matrix" --orientation auto \
        --repeat "$SOLVE_REPEAT" --factor-repeat "$FACTOR_REPEAT" \
        --refactor-repeat "$REFACTOR_REPEAT" --threads "$THREADS" \
        --no-transpose-solve \
        --input-index "$INPUT_INDEX" \
        --refactor-values "$REFACTOR_VALUES" \
        --refactor-value-amplitude "$REFACTOR_VALUE_AMPLITUDE" \
        --backend "$KLS_BACKEND" --json 2>/dev/null)
      dst=$OUT_KLS
    elif [ "$side" = ck ]; then
      out=$(timeout "$TIMEOUT" "$CKTSO_COMPARE" "$matrix" "$THREADS" \
        "$SOLVE_REPEAT" \
        "$REFACTOR_REPEAT" "$FACTOR_REPEAT" \
        "$REFACTOR_VALUES" "$REFACTOR_VALUE_AMPLITUDE" \
        2>/dev/null)
      dst=$OUT_CK
    elif [ "$side" = st ]; then
      out=$(timeout "$TIMEOUT" "$SUBTREELU_COMPARE" "$matrix" "$THREADS" \
        "$SOLVE_REPEAT" \
        "$REFACTOR_REPEAT" "$FACTOR_REPEAT" \
        "$REFACTOR_VALUES" "$REFACTOR_VALUE_AMPLITUDE" \
        2>/dev/null)
      dst=$OUT_ST
    else
      out=$(timeout "$TIMEOUT" "$KLU_COMPARE" "$matrix" \
        --repeat "$SOLVE_REPEAT" \
        --factor-repeat "$FACTOR_REPEAT" --refactor-repeat "$REFACTOR_REPEAT" \
        --refactor-values "$REFACTOR_VALUES" \
        --refactor-value-amplitude "$REFACTOR_VALUE_AMPLITUDE" \
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
