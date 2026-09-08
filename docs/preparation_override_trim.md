# Deferred-preparation override removal

Relative to `69f5d37`, removed seven development overrides:

- `KLS_DIRECT_FORCED_ROW_PREP`
- `KLS_DIRECT_FORCED_ROW_SKIP_COLUMN_PREPS`
- `KLS_DIRECT_FORCED_ROW_PATTERN_ONLY`
- `KLS_SKIP_DEFERRED_ROW_SOLVE_SEED_FOR_ROW`
- `KLS_ENABLE_DEFERRED_SNODE_OVERLAP`
- `KLS_ENABLE_DEFERRED_PREP_OVERLAP`
- `KLS_DISABLE_LARGE_SPARSE_PREP_OVERLAP`

The default preparation sequence remains intact, including the automatic
overlap predicate, sort-before-reader ordering, worker joins, thread-creation
fallback, solve seeding, row-model selection, and SNB preparation. The earlier
SNB binary-search fix and worker-affinity repair are unchanged. Explicit row
execution itself remains supported; only its experimental preparation bypasses
were removed.

Removed four constant `&& 1` predicates. Retained phase-timing diagnostics and
detached the `snode_start` and `i32` markers from obsolete `else` clauses so
they report their actual phases when tracing is enabled. Replaced unsupported
U-supernode and BTF grouped-state README instructions with a retirement notice;
the separate supported BTF scalar-run executor documentation remains.

Net reduction: **60 source lines** and **254 README lines**, excluding this
report. Public API/ABI and reserved statistics are unchanged.

## Validation

Release and ASan/UBSan CTest each pass **5/5**, with no compiler warnings from
the edited sources. An additional Release smoke run on CPUs 8–15 passes with
all seven retired switches set to 1. Removed-switch searches in production
source and `git diff --check` are clean.

Artifacts are in `build/prep-override-trim-oLNb7b/`. The frozen before binary is
the committed SNB repair. Compiler settings, dependency revisions, linked
targets, external objects, and runtime libraries match exactly across builds.
The harness records the source patch, binary/matrix hashes, and runner hash.

Benchmarks run sequentially with alternating before/after order and tracing
disabled: H100 entrywise changes of amplitude 0.001, every refactor verified
against `1e-8`, an 80GiB address-space limit, and a 90-second launch timeout.
TSOPF_FS_b9_c6 receives 20 pairs at one thread and on each eight-thread cache
domain. Fifteen other matrices receive four pairs on each eight-thread domain.
This bounded matrix-lifecycle check is not a full paper-corpus or Xyce rerun.

## Results

All **360/360 main-trial launches pass**, with maximum refactor residual
`5.45371652e-9`. TSOPF_FS_b9_c6 median paired lifecycle changes are **+0.19%**
at one thread, **-0.79%** on the 32MiB eight-thread domain, and **+0.41%** on
the 96MiB domain. Its one-thread first-refactor median is 24.345ms before and
24.350ms after: the SNB setup repair is retained.

The initial four-pair controls have two slowdowns above 2%:
TSOPF_RS_b9_c6 at 32MiB (+3.75%) and ASIC_680ks at 96MiB (+2.25%). Both receive
16 fresh alternating pairs. Those rechecks show **-0.41%** and **+0.36%**,
respectively; neither initial slowdown repeats. Other initial control changes
range from -3.99% to +1.42%. Raw per-matrix and per-phase results remain in
`validation-summary.json` and `recheck-summary.json` under the artifact directory.

All **64/64 recheck launches pass**, for **424/424 total**. These bounded
results show no repeatable material regression; they do not establish universal
performance parity. The final Release binary matches the frozen after binary,
SHA256 `3172dbf6881d78773005b1ec95d7c1559cd91ff91138516093b4b991fb8dddb2`.
