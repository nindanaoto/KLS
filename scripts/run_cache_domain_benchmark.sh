#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
IMAGE="$ROOT/build/kls-benchmark-private.sif"
VENDOR_DIR=${KLS_BENCH_VENDOR_DIR:-"$ROOT/build/benchmark-vendor"}
MATRIX_DIR="$ROOT/data/suitesparse"
MANIFEST="$ROOT/bench/suitesparse_circuit_manifest.txt"
OUTPUT_DIR="$ROOT/build/cache-domain-smallest"
DOMAIN=smallest
THREADS=8
PASSES=8
TIMEOUT=120

usage() {
  echo "usage: $0 [--image FILE] [--matrix-dir DIR] [--manifest FILE]" >&2
  echo "          [--output-dir DIR] [--domain smallest|largest]" >&2
  echo "          [--threads N] [--passes N] [--timeout SECONDS]" >&2
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --image) IMAGE=$2; shift 2 ;;
    --matrix-dir) MATRIX_DIR=$2; shift 2 ;;
    --manifest) MANIFEST=$2; shift 2 ;;
    --output-dir) OUTPUT_DIR=$2; shift 2 ;;
    --domain) DOMAIN=$2; shift 2 ;;
    --threads) THREADS=$2; shift 2 ;;
    --passes) PASSES=$2; shift 2 ;;
    --timeout) TIMEOUT=$2; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) usage; echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done
case "$DOMAIN" in smallest|largest) ;; *) echo "invalid domain: $DOMAIN" >&2; exit 2;; esac
case "$THREADS:$PASSES:$TIMEOUT" in
  *[!0-9:]*|0:*|*:0:*|*:0) echo "threads, passes, and timeout must be positive integers" >&2; exit 2;;
esac
if [ "$THREADS" -ne 8 ]; then
  echo "the paper headline configurations are pinned to --threads 8" >&2
  exit 2
fi
if [ $((PASSES % 4)) -ne 0 ]; then
  echo "passes must be a multiple of four for balanced launch order" >&2
  exit 2
fi

ENGINE=${APPTAINER:-}
if [ -z "$ENGINE" ]; then
  ENGINE=$(command -v apptainer || command -v singularity || true)
fi
if [ -z "$ENGINE" ]; then echo "Apptainer or Singularity is required" >&2; exit 1; fi
if [ ! -f "$IMAGE" ]; then
  echo "missing image: $IMAGE" >&2
  echo "run scripts/build_benchmark_container.sh first" >&2
  exit 1
fi
if [ ! -f "$VENDOR_DIR/cktso/rocky8_x64_gcc850/libcktso.so" ] ||
   [ ! -f "$VENDOR_DIR/cktso/rocky8_x64_gcc850/cktso.lic" ] ||
   [ ! -f "$VENDOR_DIR/subtreelu/lib/libsubtree_lu.so" ]; then
  echo "missing pinned solver checkouts under: $VENDOR_DIR" >&2
  echo "run scripts/build_benchmark_container.sh first" >&2
  exit 1
fi
if [ ! -d "$MATRIX_DIR" ]; then echo "missing matrix directory: $MATRIX_DIR" >&2; exit 1; fi
if [ ! -f "$MANIFEST" ]; then echo "missing manifest: $MANIFEST" >&2; exit 1; fi
MATRIX_DIR=$(cd "$MATRIX_DIR" && pwd)
MANIFEST=$(cd "$(dirname "$MANIFEST")" && pwd)/$(basename "$MANIFEST")
OUTPUT_DIR=$(mkdir -p "$OUTPUT_DIR" && cd "$OUTPUT_DIR" && pwd)
mkdir -p "$OUTPUT_DIR/tmp"

CPUS=$(python3 "$ROOT/scripts/cache_domain_cpus.py" \
  --select "$DOMAIN" --threads "$THREADS")
echo "Host LLC domains:" >&2
python3 "$ROOT/scripts/cache_domain_cpus.py" >&2
echo "Benchmarking on $DOMAIN LLC physical-core CPUs: $CPUS" >&2
case "$DOMAIN" in
  smallest) CONFIG_ID=headline_llc32_t8; EXPECTED_CPUS=8,9,10,11,12,13,14,15 ;;
  largest) CONFIG_ID=headline_llc96_t8; EXPECTED_CPUS=0,1,2,3,4,5,6,7 ;;
esac
if [ "$CPUS" != "$EXPECTED_CPUS" ]; then
  echo "paper configuration $CONFIG_ID requires CPUs $EXPECTED_CPUS, but topology selection returned $CPUS" >&2
  exit 1
fi

BUILD_DIR=/work/KLS/build-container-benchmark
BIND_PATHS="$ROOT:/work/KLS,$MATRIX_DIR:/work/matrices:ro,$MANIFEST:/work/manifest.txt:ro,$OUTPUT_DIR:/work/output,$OUTPUT_DIR/tmp:/tmp,$VENDOR_DIR/cktso:/opt/cktso:ro,$VENDOR_DIR/subtreelu:/opt/subtreelu:ro"
EXEC_FLAGS=(--no-home --no-mount home,tmp,hostfs,cwd --bind "$BIND_PATHS" --pwd /work/KLS)
IMAGE_ID=$(sha256sum "$IMAGE" | awk '{print $1}')
IMAGE_STAMP="$ROOT/build-container-benchmark/.kls-benchmark-image-id"
CMAKE_FRESH=()
if [ ! -f "$IMAGE_STAMP" ] || [ "$(cat "$IMAGE_STAMP")" != "$IMAGE_ID" ]; then
  CMAKE_FRESH=(--fresh)
fi
"$ENGINE" exec "${EXEC_FLAGS[@]}" "$IMAGE" \
  cmake "${CMAKE_FRESH[@]}" -S /work/KLS -B "$BUILD_DIR" \
    -DCMAKE_BUILD_TYPE=Release \
    -DKLS_BUILD_CKTSO_COMPARE=ON -DCKTSO_ROOT=/opt/cktso \
    -DKLS_BUILD_SUBTREELU_COMPARE=ON -DSUBTREELU_ROOT=/opt/subtreelu
printf '%s\n' "$IMAGE_ID" > "$IMAGE_STAMP"
"$ENGINE" exec "${EXEC_FLAGS[@]}" "$IMAGE" \
  cmake --build "$BUILD_DIR" -j "$THREADS" \
    --target kls_bench cktso_compare subtreelu_compare klu_compare

"$ENGINE" exec "${EXEC_FLAGS[@]}" "$IMAGE" \
  python3 /work/KLS/scripts/run_paper_campaign.py \
    --root /work/KLS \
    --manifest /work/manifest.txt \
    --matrix-dir /work/matrices \
    --configs /work/KLS/bench/paper_campaign_configs.json \
    --output-dir /work/output/campaign \
    --kls-bench "$BUILD_DIR/kls_bench" \
    --ck-bench "$BUILD_DIR/cktso_compare" \
    --st-bench "$BUILD_DIR/subtreelu_compare" \
    --klu-bench "$BUILD_DIR/klu_compare" \
    --passes "$PASSES" --timeout "$TIMEOUT" --only "$CONFIG_ID"

"$ENGINE" exec "${EXEC_FLAGS[@]}" "$IMAGE" \
  python3 /work/KLS/scripts/reduce_paper_campaign.py \
    --results /work/output/campaign \
    --json /work/output/summary.json \
    --markdown /work/output/summary.md \
    --csv /work/output/matrix-ratios.csv \
    --latex /work/output/summary-rows.tex \
    --only "$CONFIG_ID"

{
  echo "image_sha256=$IMAGE_ID"
  echo "kls_commit=$(git -C "$ROOT" rev-parse HEAD)"
  echo "kls_tracked_changes=$(git -C "$ROOT" status --short --untracked-files=no | wc -l)"
  echo "kls_diff_sha256=$(git -C "$ROOT" diff --binary HEAD | sha256sum | awk '{print $1}')"
  echo "cache_domain=$DOMAIN"
  echo "cpus=$CPUS"
  echo "threads=$THREADS"
  echo "passes=$PASSES"
  echo "manifest_sha256=$(sha256sum "$MANIFEST" | awk '{print $1}')"
  echo "cktso_commit=b9b2dae065fa776365ba009cba6d50888fdede85"
  echo "subtreelu_commit=9930fd2cadff043f44bd0c9b48a63403f6ed218e"
} > "$OUTPUT_DIR/run-metadata.txt"

echo "Strict paired summary: $OUTPUT_DIR/summary.md" >&2
