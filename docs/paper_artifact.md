# KLS paper artifact

This repository is the release artifact for the KLS paper.  It contains the
solver source, the exact deterministic value sequence, all four benchmark
harnesses, the confirmatory/diagnostic configuration matrix, and the strict
paired reducer.  The paper's primary H100 quantity is measured directly; the
legacy `spice_cycle_seconds` projection remains in JSON only for comparison
with older runs.

## Build and test

```sh
git submodule update --init --recursive
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DKLS_BUILD_CKTSO_COMPARE=ON -DCKTSO_ROOT=/path/to/cktso \
  -DKLS_BUILD_SUBTREELU_COMPARE=ON -DSUBTREELU_ROOT=/path/to/subtreelu
cmake --build build -j --target \
  kls_bench klu_compare cktso_compare subtreelu_compare kls_smoke
ctest --test-dir build --output-on-failure
```

KLU is built from the pinned vendored source.  CKTSO and SubtreeLU are not
redistributable in this repository; provide the versions and checksums stated
in the campaign metadata.  A KLS/KLU-only reproduction needs neither package.

## Direct lifecycle contract

For `--lifecycle-systems 100`, every harness executes

```text
analyze(pattern)
factor(values[0])
solve(values[0])
for generation = 1..99:
    values[generation] = deterministic_update(values[0], generation)
    refactor(values[generation])
    solve(values[generation])
```

The JSON field `measured_lifecycle_seconds` is the sum of those measured
operations.  `bench/bench_value_sequence.h` is the normative update
definition.  Every solve uses
`x_true[i] = 1 + 0.01*(i mod 17)` and a right-hand side recomputed from the
current values.  With verification enabled, a record is valid only if the
initial/final and every post-refactor relative L2 residual are at most
`1e-8`.

## Campaign

`bench/paper_campaign_configs.json` defines:

- the 110-matrix H100 comparison in each 8-core LLC domain;
- 1/2/4/8-core scaling in each domain and 16-core full-socket scaling;
- rank-preserving, 10% entrywise, and localized 10% update sensitivity;
- expected-horizon hints of 0, 10, 100, and 1000 systems;
- AUTO, serial, row-off, EGraph-off, panel-off, fast-factor-off, BTF-off,
  forced-AMD, and forced-METIS diagnostic arms.

The scaling, sensitivity, and ablation experiments use the deterministic
30-matrix manifest `bench/paper_scaling_subset_manifest.txt`.  CPU lists in the
configuration file describe the evaluation host and must be regenerated, not
blindly reused, on a host with different cache topology.

Run and reduce the campaign with:

```sh
python3 scripts/run_paper_campaign.py \
  --root . \
  --manifest bench/suitesparse_paper_manifest.txt \
  --matrix-dir data/suitesparse \
  --configs bench/paper_campaign_configs.json \
  --output-dir results/paper-v2 \
  --kls-bench build/kls_bench \
  --ck-bench build/cktso_compare \
  --st-bench build/subtreelu_compare \
  --klu-bench build/klu_compare \
  --passes 8

python3 scripts/reduce_paper_campaign.py \
  --results results/paper-v2 \
  --json results/paper-v2/summary.json \
  --markdown results/paper-v2/summary.md \
  --csv results/paper-v2/matrix-ratios.csv \
  --latex results/paper-v2/summary-rows.tex
```

The runner counterbalances launch order in four-pass blocks, records explicit
pass and launch identities, and rotates configurations sharing `pair_group`
inside each matrix/pass so scaling, sensitivity, and ablation arms remain
adjacent.  It retains all failures and resumes only missing keys.  It writes
binary/config/manifest hashes, the repository commit and diff hash, CPU
topology, commands, stderr tails, return codes, and wall times.  The reducer
accepts a matrix only when all configured passes are present and valid.  It
forms ratios within matched passes, takes the median log ratio for each matrix,
reports wins/ties/losses with a 2% tie band, and reports a matrix-bootstrap 95%
confidence interval for the corpus geometric mean.  Cross-configuration output
reports both arm/baseline time and baseline/arm speedup.

## One-use external validation

Do not reveal the v2 external corpus until the complete scientific stack is
committed.  The selector rejects a dirty or untracked solver, harness, or
campaign-script tree.  It excludes every exact matrix already named in the
development manifests and raw paper data, admits at most one unseen matrix per
SuiteSparse group, and freezes the metadata snapshot and selection seed.

```sh
python3 scripts/build_external_validation_manifest.py \
  --protocol kls-external-validation-v2 \
  --index build/ssstats.csv --root . --count 24 \
  --seed kls-external-validation-2026-09-v2 \
  --exclude-glob 'bench/*manifest*.txt' \
  --exclude-glob 'build/**/*manifest*' \
  --exclude-glob '../KLS-paper/evaluation/*manifest*' \
  --exclude-glob '../KLS-paper/evaluation/raw/*.jsonl'
```

After reveal, run the resulting manifest once with the same direct H100
configuration and make no solver, harness, threshold, or reducer changes.  A
failed or incomplete pass remains a disclosed failure; it is never replaced by
a favorable rerun.
