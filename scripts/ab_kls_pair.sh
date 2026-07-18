#!/bin/sh
# Interleaved same-moment A/B of two kls_bench binaries over a manifest.
# Usage: ab_kls_pair.sh BIN_A BIN_B MATRIX_DIR MANIFEST OUT_A OUT_B
set -u
BIN_A=${1:?bin A}; BIN_B=${2:?bin B}; DIR=${3:?matrix dir}; MAN=${4:?manifest}
OUT_A=${5:?out A}; OUT_B=${6:?out B}
THREADS=${THREADS:-8}
INPUT_INDEX=${INPUT_INDEX:-64}
TIMEOUT=${TIMEOUT:-120}
PASSES=${PASSES:-2}
RR=${REFACTOR_REPEAT:-20}
FR=${FACTOR_REPEAT:-1}
SR=${SOLVE_REPEAT:-1}
REFACTOR_VALUES=${REFACTOR_VALUES:-unchanged}
REFACTOR_VALUE_AMPLITUDE=${REFACTOR_VALUE_AMPLITUDE:-0.001}
SKIP=${SKIP:-}
ENV_A=${ENV_A:-}
ENV_B=${ENV_B:-}
: > "$OUT_A"; : > "$OUT_B"
pass=1
while [ "$pass" -le "$PASSES" ]; do
  grep -v '^#' "$MAN" | while read -r name; do
    [ -z "$name" ] && continue
    case " $SKIP " in *" $name "*) continue;; esac
    matrix=$(find "$DIR" -iname "${name}.mtx" | head -1)
    [ -z "$matrix" ] && { echo "skip $name (not found)" >&2; continue; }
    if [ $((pass % 2)) -eq 1 ]; then sides="A B"; else sides="B A"; fi
    for side in $sides; do
      if [ "$side" = A ]; then
        bin=$BIN_A; out=$OUT_A; side_env=$ENV_A
      else
        bin=$BIN_B; out=$OUT_B; side_env=$ENV_B
      fi
      res=$(timeout "$TIMEOUT" env $side_env "$bin" "$matrix" --orientation auto \
        --repeat "$SR" --factor-repeat "$FR" --refactor-repeat "$RR" \
        --input-index "$INPUT_INDEX" \
        --refactor-values "$REFACTOR_VALUES" \
        --refactor-value-amplitude "$REFACTOR_VALUE_AMPLITUDE" \
        --threads "$THREADS" --json 2>/dev/null)
      if [ -n "$res" ]; then
        printf '%s\n' "$res" >> "$out"
      else
        printf '{"matrix":"%s","timeout_or_fail":true}\n' "$matrix" >> "$out"
      fi
    done
    echo "pass $pass done $name"
  done
  pass=$((pass + 1))
done
