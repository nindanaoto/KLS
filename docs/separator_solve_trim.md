# Remove opt-in separator solve and partition overrides

Baseline: `1c7994c`. Remove the separator-based solve path end to end.
The only constructor of its private plan was gated by
`KLS_ROW_SOLVE_SEPARATOR`; no default path could create the plan. Remove
that builder and its declaration, the private worker and mode-three
dispatch, its pool runner, the diagonal executor, the dependent hybrid-BTF
capability/executor, and solve/publication branches requiring that plan.
Remove three private plan fields and cleanup. Keep the shared first-factor
queue-plan type, normal serial BTF solve, ordinary parallel row solves,
transpose solves, and their numerical and allocation fallbacks.

Also remove `KLS_ROW_REFACTOR_REBALANCE_PRIVATE` and
`KLS_ROW_REFACTOR_PROMOTE_CROSS_THREAD_PRIVATE`. Retain default dense-spiked
rebalancing, but remove the opt-in cross-thread promotion pass. None of the
three controls is referenced by tests, benchmark drivers, scripts or existing
documentation. They are no longer supported. The public C API is unchanged.
Correct the stale publication comment about the removed hybrid solver.
Net reduction: 831 source lines across four files, including six helpers.

## Validation

Final-source builds succeed without reported warnings. Release CTest passes
6/6 in 4.26 s; ASan CTest passes 6/6 in 17.51 s. Existing tests cover normal
BTF, parallel row and transpose solves. `git diff --check` passes. Searches
confirm no references to the removed private plan or helpers remain.
Executable `size` text falls by 9,056 bytes; data/BSS are unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -1.719% to +2.760% versus accepted and from
-3.450% to +1.683% versus cumulative baseline `372f9d1`. First-factor
1138_bus on eight threads/32 MB placement crosses the 2.5% review threshold
in that initial sample. A further 24 rotated triplets (72 launches) measure
+0.503% versus accepted and -0.346% versus cumulative. All 651 launches
are valid. The notable slowdown does not persist; no alignment fix is used.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries stay fixed during timing. These focused results do not establish
full-paper or cross-machine equivalence. The retained dense-spiked rebalance
is justified by source-path equivalence, not exhaustive corpus coverage.

Artifacts: `build/prep-trim-repair-dBcoha/separator-solve-trim-*` and
`separator-solve-confirm-first-*`. Runner: `validate-separator-solve-trim.py`.
Audits: `audit-focused-trim.py` and `audit-eligibility-repeat.py` adapted to
the first-factor 1138_bus confirmation job.
