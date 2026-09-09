# Remove unreachable ordering policy remnants

Source baseline: `f8c7cf6`. Remove 75 net lines: two ordering paths gated by
private thread-local selectors that are never assigned, their write-only
selection flags, associated METIS tuning branches, and orphaned enum constants.
The selectors are zero-initialized, have no address escapes or macro setters,
and cannot activate these branches. Preserve explicit values of surviving
enum constants and all live ordering decisions.

## Validation

Release CTest passes 6/6 in 4.10 s; rebuilt ASan CTest passes 6/6 in 17.33 s.
Both builds succeed without reported warnings. `git diff --check` passes.
The release text section shrinks by 64 bytes and code layout changes, so
binary equivalence cannot replace performance validation for this batch.

All 579 focused H100 launches pass the provenance/result audit. Median paired
lifecycle changes range from -3.974% to +0.886% against the accepted executable
and from -3.005% to +0.091% against cumulative baseline `372f9d1`.
No case reaches the 2.5% review threshold; no alignment fix is needed.
TSOPF_FS_b9_c6 changes +0.321% and +0.783% on the two eight-thread CPU sets.

The accepted executable is `core-state-trim-auto-kls_bench`, whose allocated
code/data match the subsequent `f8c7cf6` build except for its build-ID note
(see `unused_inline_helpers_trim.md`). The campaign kept source and binaries
fixed, rotated execution order, and checked each changed-value refactor.
The audit verifies exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals, and recomputed medians. This focused
campaign does not establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/unreachable-policy-trim-*`.
Runner: `validate-unreachable-policy-trim.py`. Audit: `audit-focused-trim.py`.

Read-only follow-up inspection identified five further always-zero private
context flags: `kls_lowdeg_ndp_class`, `kls_spiked_ndp_class`,
`kls_medium_partial_static_metis_ctx`,
`kls_low_work_symmetric_partial_diagonal_pts_ctx`, and
`kls_large_weak_diagonal_static_metis_ctx`. Their cleanup is not part of this
validated patch.
