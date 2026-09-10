# Frozen external validation

This campaign tests whether the KLS result survives matrices and value
sequences that were not used to develop the current solver.  It is designed
to answer the overfitting question on the current machine; repeating the same
protocol on another machine is a separate, later test.

## Confirmatory boundary

The one-use manifest freezes KLS commit
`ce959450b1b2a1affe60ab2813863b3c1124dad7` and hashes the complete tracked
solver tree (`CMakeLists.txt`, `include`, `src`, and `third_party`).  It was
selected from the pinned 31-Oct-2023 SuiteSparse index after resolving every
entry in every pre-existing `bench/*manifest*.txt` file and excluding all 101
SuiteSparse groups represented there.  The 24 selected matrices use 24 of the
36 remaining eligible groups, at most one matrix per group.  All three unused
circuit-, semiconductor-, or electromagnetics-like groups are included; the
rest are balanced across size, density, and structural symmetry.

Do not change solver source, thresholds, default environment controls, or the
manifest after looking at results.  A failure remains a failure.  Fixes and
new tuning belong to a later development revision and require a new unseen
suite.

The group label is an imperfect proxy for independence, and most unused
groups are not circuit problems.  Report circuit-like and complete-suite
results separately when interpreting domain performance; do not discard the
out-of-domain rows from coverage.

## Prepare matrices

Use the same pinned index that generated the manifest, not the live collection
index:

```sh
python3 scripts/build_external_validation_manifest.py \
  --root . --index build/ssstats.csv --check

python3 scripts/fetch_suitesparse.py \
  --index build/ssstats.csv \
  --manifest bench/suitesparse_external_validation_v1_manifest.txt \
  --out data/suitesparse-external-validation-v1
```

The builder refuses a tracked solver-source change and verifies that the
committed selection still follows the recorded metadata, exclusions, seed,
and source freeze.

## Run the confirmatory campaign

Build the same four harnesses used by the paper, then run:

```sh
python3 scripts/run_external_validation.py \
  --root . \
  --manifest bench/suitesparse_external_validation_v1_manifest.txt \
  --matrix-dir data/suitesparse-external-validation-v1 \
  --kls-bench build/kls_bench \
  --cktso-compare build/cktso_compare \
  --subtreelu-compare build/subtreelu_compare \
  --klu-compare build/klu_compare \
  --output-dir build/external-validation-v1
```

The fixed campaign is:

- thread counts 1, 4, and 8;
- rank-preserving changes at amplitude 0.001;
- independent entrywise changes at amplitudes 0.001 and 0.1;
- localized entrywise changes at amplitude 0.1;
- 20 changed-value refactors, 100 solves, and four counterbalanced solver
  orders per matrix/configuration;
- residual verification after every refactor, with a `1e-8` limit;
- a 300-second process limit recorded as a coverage failure, never imputed as
  a successful timing.

This is 4,608 solver processes.  It is safe to run one configuration first,
then resume into the same output directory:

```sh
python3 scripts/run_external_validation.py [the same arguments] \
  --only t8_entry_p001
```

Completed configurations have marker files and are not overwritten.  The
runner records source, manifest, binary, host, and configuration hashes in
`campaign.json`; it refuses to mix a changed campaign into the directory.  It
also removes inherited KLS/competitor tuning variables and supplies only the
documented verification and single-threaded-BLAS environment.

## Analyze without survivorship bias

```sh
python3 scripts/analyze_external_validation.py \
  --results build/external-validation-v1 \
  --json build/external-validation-v1/summary.json \
  --markdown build/external-validation-v1/report.md
```

The reducer requires all four passes for a solver/matrix to count as valid,
reports missing, invalid, and timed-out rows by name, and computes pairwise
performance only on mutually valid rows.  It reports the predeclared H10,
H100, and H1000 projections from the same measured components.  A partial
campaign exits nonzero unless `--allow-incomplete` is explicitly used for an
interim report.

The strongest evidence against overfitting would be stable aggregate wins,
coverage, and route choices across entrywise/localized workloads and thread
counts—not another all-wins count on the separable workload.  Thin margins
within 2% are reported as ties.

## Later cross-machine run

Copy the repository at the frozen commit, the manifest, pinned metadata, and
matrices.  Rebuild all harnesses normally, use a new output directory, and run
the identical command.  Different binary hashes are expected and recorded;
the frozen source-tree digest and campaign definition must remain identical.
