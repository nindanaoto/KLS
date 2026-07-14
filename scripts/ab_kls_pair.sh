#!/bin/sh
# Interleaved same-moment A/B of two kls_bench binaries over a manifest.
# Usage: ab_kls_pair.sh BIN_A BIN_B MATRIX_DIR MANIFEST OUT_A OUT_B
set -u
BIN_A=${1:?bin A}; BIN_B=${2:?bin B}; DIR=${3:?matrix dir}; MAN=${4:?manifest}
OUT_A=${5:?out A}; OUT_B=${6:?out B}
THREADS=${THREADS:-8}
TIMEOUT=${TIMEOUT:-120}
PASSES=${PASSES:-2}
RR=${REFACTOR_REPEAT:-20}
SKIP=${SKIP:-}
: > "$OUT_A"; : > "$OUT_B"
pass=1
while [ "$pass" -le "$PASSES" ]; do
  grep -v '^#' "$MAN" | while read -r name; do
    [ -z "$name" ] && continue
    case " $SKIP " in *" $name "*) continue;; esac
    matrix=$(find "$DIR" -iname "${name}.mtx" | head -1)
    [ -z "$matrix" ] && { echo "skip $name (not found)" >&2; continue; }
    for side in A B; do
      if [ "$side" = A ]; then bin=$BIN_A; out=$OUT_A; else bin=$BIN_B; out=$OUT_B; fi
      res=$(timeout "$TIMEOUT" "$bin" "$matrix" --orientation auto \
        --repeat 1 --factor-repeat 1 --refactor-repeat "$RR" \
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
