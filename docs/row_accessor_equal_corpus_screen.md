# Cumulative corpus checkpoint after row-helper cleanup

Candidate: `eaec5d7`; resumed-pass baseline: `372f9d1`.
Four source chunks since the preceding TRSV checkpoint remove 68 net C-source
lines. This cumulative check covers all 93 medium-paper matrices on both
eight-thread cache domains, using H100 entrywise-changed refactors and solves.

All 186 matrix/domain cases were processed: 1288 launches, 1278 valid.
There are 181 valid comparisons and five inconclusive cases. Both versions
hit the 90-second limit for ss1 on CPUs 8–15. Both versions report
`solve failed (-6)` for OPF_10000 and OPF_3754 on both domains (process exit
1, solver status -6). These shared failures are not candidate regressions
and do not count as valid performance evidence.

Eight initial positive screening flags received twelve-pair repeats:

| Matrix / domain | Initial change | Repeat change |
| --- | ---: | ---: |
| circuit_1 / CPUs 0–7 | +2.936% | -2.066% |
| circuit_3 / CPUs 0–7 | +2.973% | +1.201% |
| circuit_4 / CPUs 0–7 | +2.654% | -0.277% |
| rajat17 / CPUs 8–15 | +2.780% | -1.065% |
| scircuit / CPUs 8–15 | +4.963% | +0.199% |
| circuit204 / CPUs 0–7 | +4.168% | -1.519% |
| Hamrle2 / CPUs 8–15 | +2.635% | +1.650% |
| LeGresley_2508 / CPUs 8–15 | +2.525% | -0.092% |

None persists at the 2.5% review trigger. Among valid final comparisons
(repeat where available, otherwise three-pair screen), changes range from
-7.969% for fpga_trans_01 on CPUs 0–7 to +2.423% for fpga_dcop_01 on CPUs
0–7. These three-pair extremes are not independently confirmed speed claims.
TSOPF_FS_b9_c6 measures +2.114% / +1.051% on CPUs 8–15 / CPUs 0–7;
TSOPF_RS_b9_c6 measures -3.101% / -3.674%. The results establish no persistent
approximately 3% stopping-rule regression, not universal zero loss.

Both versions use the same sequentially copied, hash-verified execution path.
Limits remain 90 seconds and 64 GiB address space per launch. No builds or
source edits overlap timing. Artifacts are under
`build/resumed-corpus-chhiOV/row-accessor-equal0*`; runner:
`scan-row-accessor-equal.py`. The candidate is frozen at
`build/prep-trim-repair-dBcoha/row-accessor-trim-auto-kls_bench`.

Audit verifies all ordered manifest jobs, runner/manifest/matrix/binary
hashes, build provenance, empty source diff, complete commands, rotated
pair ordering, launch counts, residual validation at 1e-8, and all paired
and per-version medians recomputed from raw records. The final error audit
distinguishes the OPF solver status from its process exit code and verifies
that each inconclusive outcome occurs on both sides.

This checkpoint permits the next low-risk source audit. The broader cleanup
goal remains active.
