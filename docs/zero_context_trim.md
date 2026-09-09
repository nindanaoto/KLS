# Remove always-zero ordering contexts

Baseline: `65e0541`. Remove five private thread-local flags that are
zero-initialized, assigned only zero, and have no address escapes or macro
setters: `kls_lowdeg_ndp_class`, `kls_spiked_ndp_class`,
`kls_medium_partial_static_metis_ctx`,
`kls_low_work_symmetric_partial_diagonal_pts_ctx`, and
`kls_large_weak_diagonal_static_metis_ctx`.

Remove their reset stores and obsolete comments; simplify matching and METIS
conditions and unreachable leaf-count choices. Keep live environment controls,
generic tuning contexts, and numerical behavior. Net reduction: 58 lines.

## Validation

Release CTest passes 6/6 in 4.13 s; ASan passes 6/6 in 17.21 s. Builds succeed
without reported warnings, and `git diff --check` passes. Code and TLS layout
change, so this batch receives performance checks rather than a binary-
equivalence waiver.

All 579 initial focused H100 launches pass correctness and provenance audits.
Initial median paired lifecycle flags: first-factor circuit_4 +4.796% against
accepted and +8.909% against cumulative baseline `372f9d1`; TSOPF_RS_b9_c6
+2.899% cumulative. LeGresley_2508 is near the review threshold at +1.958%
cumulative. These results triggered longer repeats without source changes.

The 216-launch confirmation uses 24 rotating-order triplets per case:

| Case | Change vs accepted | Change vs cumulative |
| --- | ---: | ---: |
| circuit_4, first factor | +0.129% | +0.055% |
| TSOPF_RS_b9_c6 | -0.194% | -0.433% |
| LeGresley_2508 | +0.397% | -1.368% |

The initial slowdowns do not persist in the longer repeat. No alignment fix
is needed. All 795 launches pass residual validation; audits verify jobs,
commands, source revision/diff, provenance, executable/matrix hashes, exits,
and recomputed medians. Source and binaries remain fixed during timing.
This focused evidence does not establish full-paper or cross-machine
performance equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/zero-context-trim-*` and
`zero-context-confirm-*`. Runner: `validate-zero-context-trim.py`.
Initial audit: `audit-focused-trim.py`; confirmation audit uses
`audit-eligibility-repeat.py` with explicit first/repeat jobs and total 216.
