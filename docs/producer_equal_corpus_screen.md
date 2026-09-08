# Shared-path cumulative medium-corpus screen

Compare numeric source `78a4c89` against resumed-pass baseline `372f9d1`:
seven cleanup chunks, 467 fewer net C-source lines. Both retain the SNB
alignment fix. No source edits, builds or correctness tests overlap timing.

This corrects the executable-path confound documented in
`add32_corpus_launch_check.md`. Before every sequential launch, the chosen
frozen binary is copied to a shared temporary executable path and its
SHA-256 verified. Both sides launch through that same path with identical
workload arguments. Copying and hashing are outside the lifecycle timer.
Frozen inputs remain unchanged.

## Scope and method

All 93 matrices in `bench/suitesparse_paper_medium_manifest.txt` are present
and uniquely resolved. Test eight threads on CPUs 8-15 and CPUs 0-7: 186
cases. Three alternating baseline/current pairs screen each case; a +2.5%
median lifecycle flag receives 12 additional pairs. The stopping rule is
approximately 3% repeatable lifecycle slowdown or a correctness regression.

Use automatic solver policies, 64-bit input indices, H100 direct lifecycles
with entrywise-changed refactors at amplitude 0.001, and verification of
every changed refactor at 1e-8. Each launch has a 90-second timeout and a
64 GiB address-space cap. The environment is cleaned and CPU sets pinned.
This is not the large supplement, all paper configurations, Xyce, or a
cross-machine campaign.

## Completed result

All 186 cases were processed; the runner exited successfully. The summary
matches the metadata job list exactly. Of 1,384 launches, 1,374 are valid,
yielding 181 valid case comparisons and five inconclusive cases:

- `ss1`, CPUs 8-15: both versions timed out after 90 seconds.
- `OPF_10000` and `OPF_3754`, both CPU domains: both versions exited with
  `KLS benchmark failed: solve failed (-6)` and no JSON result.

These are not candidate-only regressions or successful checks. `ss1` on
CPUs 0-7 completed at +0.09% paired lifecycle change.

All 12 initial slowdown flags fell below the review threshold on repeat.
Repeat medians: ASIC_680ks +0.64%, ASIC_100k -1.18%, circuit_1 +0.27%,
circuit_4 +1.55%, adder_dcop_01 +0.14%, adder_trans_01 +0.59%, bips98_606
-2.48%, circuit204 +0.76%, LeGresley_2508 +2.20%, mult_dcop_03 -0.65%,
TSOPF_FS_b39_c19 +0.55%, ww_vref_6405 +0.29%. Exact configurations and
individual samples are recorded in the summaries and JSONL.

Using the repeat when present, the largest valid median lifecycle increase
was +2.46% on `mimo46x46_system`, CPUs 0-7 (three screen pairs).
The earlier confounded cases passed: `add32` -0.26% / -1.06%, and
`hcircuit` -0.75% / +0.33%, on CPUs 8-15 / CPUs 0-7 respectively.
`TSOPF_FS_b9_c6` measured -1.33% / +0.19% on those same CPU sets.

No notable stopping regression was confirmed. The medium-corpus gate
permits further cleanup, but does not prove the inconclusive cases or
untested workloads regression-free. Future timed comparisons should use
shared execution paths to avoid the demonstrated launch-path confound.

## Evidence

`build/resumed-corpus-chhiOV/scan-producer-equal.py` and `producer-equal0-*`
contain the runner, complete progress log, exact commands/results, summaries,
and metadata with build provenance, matrix/binary/manifest/runner hashes,
limits and execution-copy policy. The recorded source diff is empty.
Final verification confirms unchanged frozen binary hashes and that every
recorded command uses the shared execution path.
