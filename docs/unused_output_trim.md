# Remove unused private helper outputs

Source baseline: `e984bd1`. Remove four output parameters that every caller
passes as NULL:

- `kls_pivot_tail_plan_local_range`: contiguity output (two callers).
- `kls_row_refactor_group_trailing_slice`: start output (six callers).
- `kls_compact_singleton_run_solve_profile`: singleton count and maximum-run
  outputs (two callers).

Update declarations and callers, and remove only the associated output
initialization/stores. Keep live counts, validation, eligibility decisions,
slice construction, and public statistics. Net reduction: 26 source lines.

## Validation

Release CTest passes 6/6 in 4.15 s; ASan passes 6/6 in 17.24 s. Builds succeed
without reported warnings; `git diff --check` passes. Generated code differs,
so this cleanup receives a timing campaign rather than an equivalence waiver.

All 579 focused H100 launches pass the provenance/result audit. Median paired
lifecycle changes range from -1.850% to +0.502% versus accepted, and from
-3.251% to +1.668% versus cumulative baseline `372f9d1`. No case reaches the
2.5% review threshold. No alignment adjustment is needed.

The accepted executable is `detndp-state-trim-auto-kls_bench`, corresponding
to `06b6f9d`; intervening `e984bd1` changes comments only and preserves the
non-comment token sequence. The audit checks exact jobs, commands, source
revision/diff, build provenance, matrix/binary hashes, exits, residuals, and
recomputed medians. Source and binaries remain fixed during timing.
Focused results do not establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/unused-output-trim-*`.
Runner: `validate-unused-output-trim.py`. Audit: `audit-focused-trim.py`.
