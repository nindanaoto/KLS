#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
IMAGE=${1:-"$ROOT/build/kls-benchmark-private.sif"}
VENDOR_DIR=${KLS_BENCH_VENDOR_DIR:-"$ROOT/build/benchmark-vendor"}
ENGINE=${APPTAINER:-}
if [ -z "$ENGINE" ]; then
  ENGINE=$(command -v apptainer || command -v singularity || true)
fi
if [ -z "$ENGINE" ]; then
  echo "Apptainer or Singularity is required" >&2
  exit 1
fi

mkdir -p "$(dirname "$IMAGE")" "$VENDOR_DIR"
echo "Building the Ubuntu 26.04 benchmark toolchain image." >&2
if ! "$ENGINE" build --fakeroot --force "$IMAGE" "$ROOT/containers/benchmark.def"; then
  # Some container hosts permit execution with a privileged Apptainer engine
  # but prohibit the unprivileged mount namespace that --fakeroot requires.
  # Keep the normal rootless route first, then use explicitly available sudo
  # rather than making the reproducible build depend on manual retry steps.
  if ! sudo -n "$ENGINE" build --force "$IMAGE" "$ROOT/containers/benchmark.def"; then
    echo "Container build failed with --fakeroot and passwordless sudo." >&2
    exit 1
  fi
fi

if [ ! -d "$VENDOR_DIR/cktso/.git" ]; then
  git clone https://github.com/chenxm1986/cktso.git "$VENDOR_DIR/cktso"
fi
git -C "$VENDOR_DIR/cktso" fetch origin
git -C "$VENDOR_DIR/cktso" checkout --detach b9b2dae065fa776365ba009cba6d50888fdede85
cp "$VENDOR_DIR/cktso/license/cktso.lic" \
  "$VENDOR_DIR/cktso/rocky8_x64_gcc850/cktso.lic"

if [ ! -d "$VENDOR_DIR/subtreelu/.git" ]; then
  git clone https://github.com/THU-numbda/SubtreeLU.git "$VENDOR_DIR/subtreelu"
fi
git -C "$VENDOR_DIR/subtreelu" fetch origin
git -C "$VENDOR_DIR/subtreelu" checkout --detach 9930fd2cadff043f44bd0c9b48a63403f6ed218e

echo "Third-party solvers remain separate, ignored host checkouts and are" >&2
echo "bind-mounted read-only at runtime. Follow their license terms." >&2
echo "$IMAGE"
