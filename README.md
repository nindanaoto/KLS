# KLS

KLS is a standalone sparse direct solver project targeting SPICE-style
workloads: one symbolic analysis followed by many numeric factorizations,
refactorizations, and solves with a fixed sparse structure.

This repository currently contains the first working baseline:

- A stable C API in `include/kls/kls.h`
- A vendored SuiteSparse-derived 64-bit symbolic/numeric engine
- CSC and CSR input paths with 32-bit or 64-bit index arrays
- AMD-first automatic symbolic ordering with explicit AMD, COLAMD, and natural
  ordering controls
- Factor, refactor, solve, transpose-solve, and statistics APIs
- A MatrixMarket benchmark tool
- A small correctness smoke test
- A SuiteSparse Matrix Collection downloader script for public benchmark cases

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The main build is self-contained. KLS vendors the SuiteSparse-derived KLU,
AMD, COLAMD, and BTF C sources from Trilinos under `third_party/suitesparse`
as the current in-tree serial engine.

## Benchmark

```sh
./build/kls_bench ../cktso/demo/add20.mtx --repeat 20 --refactor-repeat 20 --json
```

The benchmark reports analysis, factorization, refactorization, solve,
transpose-solve, residual, fill, flop, and memory statistics.

To fetch public SuiteSparse Matrix Collection matrices listed in the manifest:

```sh
python3 scripts/fetch_suitesparse.py --manifest bench/suitesparse_circuit_manifest.txt --out data/suitesparse
```

To run every downloaded matrix and compute the SPICE-cycle geometric mean:

```sh
python3 scripts/run_bench_suite.py --kls-bench build/kls_bench --matrix-dir data/suitesparse --jsonl build/kls_suite.jsonl
```

The suite metric is:

```text
analysis + factor + solve + 99 * (refactor + solve)
```

An optional CKTSO comparison tool can be built when you provide a local CKTSO
distribution:

```sh
cmake -S . -B build-cktso -DKLS_BUILD_CKTSO_COMPARE=ON -DCKTSO_ROOT=/path/to/cktso
cmake --build build-cktso -j --target cktso_compare
LD_LIBRARY_PATH=/path/to/cktso/rocky8_x64_gcc850 ./build-cktso/cktso_compare matrix.mtx 16 10
```

CKTSO must be licensed correctly according to its own distribution
requirements, usually by colocating the license file with the selected shared
library.

## Status

This is a functional baseline, not yet a CKTSO-beating implementation. The next
algorithmic work is to evolve the vendored numeric kernel into KLS-owned
parallel factor/refactor/solve kernels while keeping the public API and
benchmark harness stable.

The `threads` option is accepted for API stability but is not used by the
current serial KLS numeric engine.

## License

KLS is licensed under LGPL-2.1-or-later. The current in-tree solver engine
includes SuiteSparse-derived KLU, AMD, COLAMD, BTF, and UFconfig sources from
Trilinos; see `THIRD_PARTY_NOTICES.md` for attribution.
