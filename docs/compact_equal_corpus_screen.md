# Cumulative corpus screen after compact and tracing cleanup

Candidate: `168a415`; resumed-pass baseline: `372f9d1`.
This gate follows the four cleanup chunks documented in `panel_trace_trim.md`,
`scalar_progress_trim.md`, `compact_publish_trim.md`, and `compact_shape_trim.md`:
175 further net C-source lines removed since the preceding corrected corpus
screen. No additional source changes or builds overlapped timing.

All 186 cases in the 93-matrix medium paper manifest completed, using eight
threads on CPUs 8–15 and 0–7. Each screen uses three alternating baseline/current
pairs; an initial lifecycle slowdown of at least 2.5% triggers twelve additional
pairs. Both frozen binaries execute through the same temporary path, copied and
SHA-256 verified before each launch. The H100 workload uses 99 entrywise-changed
refactors and 100 solves, automatic policies, and residual validation at 1e-8.
Each launch has a 90-second timeout and a 64-GiB address-space limit.

Artifacts are under `build/resumed-corpus-chhiOV/compact-equal0*`; the runner is
`scan-compact-equal.py` in that directory. Audit confirmed all 186 ordered jobs,
the empty source diff, manifest/runner/frozen-binary/matrix hashes, shared
execution paths, successful exit status and residual validation for every valid
launch, launch counts, and paired lifecycle medians against the raw records.

There are 1,264 launches, of which 1,254 are valid, covering 181 valid cases.
Five cases are inconclusive: both versions time out for `ss1` on CPUs 8–15,
and both report `solve failed (-6)` for `OPF_10000` and `OPF_3754` in each
configuration. These matching baseline limitations are not passes or regressions.

Seven initial flags all fall below the review threshold in their twelve-pair
repeats:

| Matrix | CPUs | Repeat lifecycle change |
| --- | --- | ---: |
| circuit_3 | 0–7 | +1.752% |
| rajat16 | 0–7 | -1.412% |
| rajat18 | 8–15 | -0.552% |
| rajat21 | 0–7 | +0.703% |
| mult_dcop_01 | 8–15 | +0.128% |
| mult_dcop_02 | 8–15 | +0.672% |
| zeros_nopss_13k | 8–15 | +0.240% |

Using repeat medians where available, valid case changes range from -5.945%
(`circuit_1`, CPUs 8–15) to +2.417% (`mult_dcop_02`, CPUs 0–7, three pairs).
`TSOPF_FS_b9_c6` measures +0.437% / -0.105%, and `onetone1` measures
+0.476% / -0.439%, for CPUs 8–15 / 0–7 respectively. No persistent review
flag or approximately 3% stopping-rule regression was found.

Retain the accepted cleanup and resume the low-risk source audit. This gate
does not prove zero performance loss, cross-machine equivalence, or exhaustion
of cleanup candidates; the broader cleanup goal remains active.
