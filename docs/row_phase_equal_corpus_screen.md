# Row-phase cleanup cumulative corpus screen

Candidate source: `161a14d` (report revision `53c6bcc`); resumed-pass
baseline: `372f9d1`. This checkpoint covers the four diagnostic cleanups
since `0f8578e`, removing 292 net C-source lines.

The complete 93-matrix paper medium manifest was screened on CPUs 8–15
and CPUs 0–7: 186 cases, 1,432 launches, 1,422 valid launches, 181 valid
comparisons, and five inconclusive comparisons. H100 entrywise-changed
refactors use AUTO policies and residual validation at 1e-8. Each case
starts with three alternating baseline/current pairs; positive lifecycle
changes of at least 2.5% receive an independent 12-pair confirmation.
Both binaries run sequentially through the same hash-verified executable
path. No builds or source edits overlap the screen.

All 14 initial flags clear on confirmation:

| Matrix | CPUs | Initial change | Confirmation |
| --- | --- | ---: | ---: |
| hcircuit | 0–7 | +2.511% | -2.208% |
| ckt11752_tr_0 | 8–15 | +3.030% | -0.753% |
| rajat17 | 0–7 | +2.834% | +0.001% |
| rajat18 | 0–7 | +3.511% | -0.263% |
| trans5 | 0–7 | +3.685% | -0.203% |
| fpga_dcop_01 | 0–7 | +3.438% | -0.181% |
| fpga_trans_01 | 8–15 | +3.920% | +0.325% |
| hvdc2 | 8–15 | +2.829% | +0.329% |
| hvdc2 | 0–7 | +4.994% | -0.177% |
| mult_dcop_02 | 0–7 | +9.678% | -0.553% |
| mult_dcop_03 | 0–7 | +3.580% | +0.265% |
| power197k | 8–15 | +3.368% | +0.544% |
| rajat12 | 0–7 | +5.001% | -0.343% |
| zeros_nopss_13k | 0–7 | +3.148% | -0.031% |

Final valid paired-median changes range from -9.497% (LeGresley_2508,
CPUs 8–15) to +2.231% (rajat17, CPUs 8–15). Unflagged three-pair screens
are not independently confirmed performance bounds. TSOPF_FS_b9_c6
changes are +0.159% / +1.301%, and TSOPF_RS_b9_c6 changes are +0.449% /
-1.258%, on CPUs 8–15 / CPUs 0–7 respectively.

Inconclusive outcomes match both versions and the preceding corpus screen:
ss1 on CPUs 8–15 times out at 90 seconds in both; OPF_10000 and OPF_3754
on both domains exit with process status 1 and `solve failed (-6)`.
These cases do not establish performance equivalence. ss1 on CPUs 0–7
is valid and measures +0.230%.

Artifacts: `build/resumed-corpus-chhiOV/row-phase-equal0*`; runner:
`scan-row-phase-equal.py`; frozen candidate:
`build/prep-trim-repair-dBcoha/row-phase-trim-auto-kls_bench`.
The final audit verifies the complete ordered manifest, revision and empty
source diff, build provenance, runner/binary/matrix hashes, every command,
all valid residuals and exits, matching invalid outcomes, record counts,
and all summary paired values and medians recomputed from raw records.

No repeatable approximately 3% lifecycle regression or new correctness
failure was found in this screen. This permits the next low-risk cleanup;
it does not prove zero loss or cross-machine generality. The broader
cleanup goal remains active.
