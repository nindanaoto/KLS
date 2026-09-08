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

The screen is in progress. Completion or a verified stopping regression
must be established before making an acceptance claim for the wider corpus.

At the commit-and-push checkpoint, cases 0-42 (43 of 186 planned cases)
had completed with valid results and no confirmed stopping regression.
The initial slowdowns on `circuit_1` (CPUs 0-7) and `hcircuit` (CPUs 8-15)
did not persist above the review threshold in their 12-pair repeats.
The campaign remains running; this snapshot is not a final benchmark result.
