# Cumulative medium paper-corpus check after resumed slimming

Compare committed `bc859b4` (three cleanup chunks, 270 fewer C-source lines)
with the resumed-pass baseline `372f9d1`. Both retain the SNB alignment fix.
No additional source deletion proceeds while this checkpoint is unresolved.

The authoritative medium manifest contains 93 locally available matrices.
The screen covers eight threads on CPUs 8-15 and 0-7: 186 planned cases.
It uses H100 direct lifecycles with entrywise changes, amplitude 0.001,
automatic solver policies and 64-bit input indices. Every changed refactor
is checked at a 1e-8 residual limit. Three alternating baseline/current pairs
screen each case. Suspects at +2.5% or more get 12 additional alternating
pairs before review; that trigger is deliberately below the approximately
3% notable-regression stopping rule and is not itself a regression verdict.

Runs are sequential with clean environments, a 90-second per-launch timeout,
and a 64 GiB address-space cap. Baseline limitations are recorded separately
from candidate-only failures. Matrices are resolved case-insensitively and
uniquely; none was missing. This is the medium paper corpus, not the large
overnight supplement, all paper configurations, or Xyce.

Artifacts are under `build/resumed-corpus-chhiOV/`: `scan.py`, frozen
executables, exact commands and results in JSONL, per-case summaries, and
metadata containing matrix/binary/manifest/runner hashes and build provenance.
The baseline hash is checked against the accepted alignment-fix campaign;
the candidate hash matches the preceding rejection-cleanup validation.
The captured source diff is empty. No benchmark overlaps builds or tests.

## Completed result

The runner exited successfully after all 186 cases. The summary matches the
metadata job list exactly, and both frozen binary hashes still match their
recorded hashes. Of 1,264 launches, 1,254 were valid, producing 181 valid
case comparisons and five inconclusive cases:

- `ss1`, CPUs 8-15: both versions exceeded the 90-second launch limit.
- `OPF_10000` and `OPF_3754`, both CPU domains: both versions exited with
  `KLS benchmark failed: solve failed (-6)` and no JSON result.

These limitations are not candidate-only regressions or successful checks.
`ss1` on CPUs 0-7 completed with a -0.22% paired lifecycle change.

Seven initial slowdown flags received 12-pair repeats. Their median paired
lifecycle changes were `circuit_1` -0.76%, `hcircuit` +0.51%, `rajat17`
-2.28%, `rajat18` -1.45%, `rajat21` -0.18%, `gemat11` +1.42%, and
`qh1484` -1.07%. Only `hcircuit` used CPUs 8-15; the others used CPUs 0-7.
None persisted above the 2.5% review trigger. The `gemat11` repeat was slower
in all 12 pairs, so this is a small observed slowdown, not a zero-regression
claim; its median remains below the approximately 3% stopping rule.

Among valid cases, using the repeat when present, the largest median
lifecycle increase was +2.46% on `rajat25`, CPUs 8-15 (three screen pairs).
The historical `TSOPF_FS_b9_c6` regression case measured +0.34% on CPUs
8-15 and -3.19% on CPUs 0-7. No notable stopping regression was confirmed.
The medium-corpus checkpoint permits the next cleanup trial; it does not
prove the five inconclusive cases or untested configurations regression-free.
