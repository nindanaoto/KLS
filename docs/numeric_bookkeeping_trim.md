# Remove write-only numeric bookkeeping

Baseline: `9708fc0`. Remove three private solver fields and their assignments:
`numeric_from_pipe`, `padded_run_count`, and `refactor_l_sorted_entries`.
Whole-source inspection finds no reads. Assignment expressions only read
local counts or the thread configuration; they have no side effects.
Keep the actual padded-run arrays, sorted-value storage, numeric preparation,
and thread configuration unchanged.

Also remove an obsolete comment referring to the deleted pipe flag and a
no-op `#undef KLS_FREE_GREEDY_MATCH` preceding its only definition. The net
reduction is 20 source lines. Public headers and dependency sources are
unchanged.

## Validation

Release CTest passes 6/6 (4.15 s); ASan CTest passes 6/6 (17.30 s).
Both builds are warning-free; `git diff --check` passes.

All 579 focused H100 launches pass residual validation. Median paired
lifecycle changes range from -2.508% to +0.668% versus accepted and from
-3.188% to +0.409% versus cumulative `372f9d1`. No case crosses the 2.5%
review threshold. TSOPF_FS_b9_c6 improves 1.171% on CPUs 8–15 and 2.508%
on CPUs 0–7 versus accepted. Onetone1 changes +0.068% and -0.034% on the
two single-thread CPU configurations. No additional alignment fix is needed.

The audit verifies exact jobs, pinned commands, source diff, build provenance,
matrix and executable hashes, exits, residuals, and recomputed medians.
Source and binaries remained fixed during timing. This focused comparison
does not establish full-paper or cross-machine performance equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/numeric-bookkeeping-trim-*`.
Runner: `validate-numeric-bookkeeping-trim.py`.
Audit: `audit-focused-trim.py`.
