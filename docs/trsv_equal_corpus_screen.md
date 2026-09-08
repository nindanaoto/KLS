# Cumulative corpus screen after optional-output and gate cleanup

Candidate: `6f16d4a`; resumed-pass baseline: `372f9d1`.
The four chunks since the preceding cumulative checkpoint removed 63 net
C-source lines: redundant native-panel eligibility plumbing, discarded output
locals, unused reserve arguments/prefix counts, and the dead compact-TRSV AUTO
gate. No source changes or builds overlap this screen.

All 186 jobs in the 93-matrix medium paper manifest completed, with eight
threads on CPUs 8–15 and 0–7. Each screen uses three alternating pairs, followed
by twelve additional pairs for initial lifecycle slowdowns of at least 2.5%.
Both frozen versions execute through one shared path with per-launch copy/hash
verification. The H100 workload uses 99 entrywise-changed refactors and 100
solves, automatic policies, and residual checks at 1e-8. Limits are 90 seconds
and 64 GiB of address space per launch.

Artifacts: `build/resumed-corpus-chhiOV/trsv-equal0*`; runner:
`scan-trsv-equal.py` in that directory. Audit verified all ordered jobs, the
empty source diff, manifest/runner/frozen-binary/matrix hashes, shared launch
paths, launch counts, exit status and residual checks for valid records, and
paired lifecycle medians against raw records.

There are 1,312 launches, 1,302 valid, covering 181 valid cases. Five cases are
inconclusive: both versions time out for `ss1` on CPUs 8–15, and both fail with
`solve failed (-6)` for `OPF_10000` and `OPF_3754` in each configuration.
Matching baseline limitations are neither passes nor regressions.

All nine initial flags fall below the review threshold in twelve-pair repeats:

| Matrix | CPUs | Repeat lifecycle change |
| --- | --- | ---: |
| ASIC_320ks | 8–15 | +0.065% |
| bcircuit | 8–15 | +0.745% |
| rajat17 | 0–7 | -0.484% |
| bips98_1142 | 8–15 | -1.068% |
| bips98_606 | 8–15 | +2.320% |
| case9 | 8–15 | +1.189% |
| hvdc2 | 0–7 | -1.351% |
| mimo28x28_system | 8–15 | +0.321% |
| TSOPF_RS_b9_c6 | 8–15 | +0.672% |

Using repeats where available, observed valid changes range from -12.648%
(`mult_dcop_03`, CPUs 0–7, three pairs) to +2.320% (`bips98_606`, CPUs 8–15).
This screening range is not a claim of confirmed speedups for every negative
result. TSOPF_FS_b9_c6 measures +1.402% / -1.542% on CPUs 8–15 / 0–7.
TSOPF_RS_b9_c6 initially measures +3.192% on CPUs 8–15, but its independent
repeat measures +0.672%; that slowdown again does not persist. No persistent
review flag or approximately 3% stopping-rule regression was found.

Retain the accepted cleanup and resume the low-risk audit. This gate does not
prove zero loss, cross-machine equivalence, or exhaustion of cleanup candidates;
the broader goal remains active.
