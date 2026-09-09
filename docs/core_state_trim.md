# Remove unused core state and failure details

Baseline: `7c9364c`. This batch removes 75 net source lines across six files:

- Four write-only first-factor failure details and their pure diagnostic
  helper. Keep the failure reason, failed-row tracking, and all retry decisions.
- Two write-only private transpose-plan counts. Keep the local tail-entry
  count used to construct slices; request no unused thread-bound maximum.
- The unused supernode parallel-work share and its exclusive accumulator.
  Keep the cooperative-work share used by scheduling.
- The unused dense-producer-run identifier. Keep its group and length.

Public headers, numerical policies, arithmetic, and dependency sources are
unchanged. The remaining forced-pivot debug hook is used by the separator
pivot/restart correctness test and is not dead diagnostic code.

## Validation

Release and ASan CTest each pass 6/6 (4.08 s and 17.14 s respectively).
Both build logs contain no warnings or errors; `git diff --check` passes.

All 579 focused H100 launches pass residual validation. Median paired
lifecycle changes range from -2.153% to +1.476% versus the accepted baseline
and from -3.709% to +0.976% versus cumulative baseline `372f9d1`.
No case reaches the 2.5% regression-review threshold. TSOPF_FS_b9_c6 changes
-0.365% on CPUs 8–15 and -0.122% on CPUs 0–7 versus accepted. No alignment
adjustment is needed.

The completed audit verifies exact jobs, pinned commands, source revision
and diff, build provenance, binary and matrix hashes, successful exits,
residuals, and recomputed medians. Source and binaries remained fixed during
timing. This focused campaign does not establish full-paper or cross-machine
performance equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/core-state-trim-*`.
Runner: `validate-core-state-trim.py`. Audit: `audit-focused-trim.py`.
