# Remove vendor factor-wrapper diagnostics

Baseline: `ca42199`. Remove `KLS_TRACE_FILL` and `KLS_TRACE_BLOCK_TIME`
diagnostics from the vendored KLU factor wrapper, reducing source by 73 net
lines. This removes logging, associated timespec variables, unconditional
clock reads, and now-unneeded POSIX/time declarations. The scale and kernel
factor calls, numeric statistics, allocations, error handling, and pivot
permutations remain unchanged. Both integer-width variants include this
source. No tests or benchmark scripts reference the removed trace modes.

The kernel's construction and pipeline profiler remains for a separate
review; this patch does not remove its counters or numerical controls.

## Validation

Release CTest passes 6/6 (4.11 s); ASan CTest passes 6/6 (17.21 s).
Both builds are warning-free; `git diff --check` passes.

All 579 focused H100 launches pass residual validation. Median paired
lifecycle changes range from -1.594% to +2.334% versus accepted and from
-3.273% to +0.604% versus cumulative `372f9d1`. TSOPF FS changes +0.370%
and +0.466% on the two eight-thread CPU sets. Onetone1 changes -0.408%
and -0.199% on the two single-thread configurations.

The near-threshold LeGresley_2508 result (+2.334% incremental) receives an
independent 24-triplet confirmation on CPUs 8–15. All 72 launches pass;
the repeated lifecycle change is -0.176% incremental and -0.004% cumulative.
The initial slowdown does not reproduce. No alignment fix is needed.

All 651 launches have audited jobs, pinned commands, source diffs, executable
and matrix hashes, exits, residuals and recomputed medians. The build proof
explicitly records the reviewed dependency patch SHA-256:
`b2acdf56770ea9addd4a7f6cab4b4827c67a531573736e2fea998adf8d8eecc6`.
It does not silently treat changed dependency sources as identical. Source
and build files stayed fixed throughout timing. This focused validation is
not a full-paper or cross-machine performance guarantee.

Artifacts: `build/prep-trim-repair-dBcoha/factor-wrapper-trace-trim-*` and
`factor-wrapper-trace-confirm-repeat*`. Runner: `validate-factor-trace-trim.py`.
Main audit: `audit-factor-trace-trim.py`. Confirmation uses
`audit-eligibility-repeat.py` with one 24-triplet LeGresley job, total 72,
and the same explicit dependency-patch hash.
