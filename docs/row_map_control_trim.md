# Remove unused row and compact-map controls

Baseline: `bcf7f9e`. Remove eleven environment controls with no references
in tests, scripts, benchmark drivers or existing documentation:

- `KLS_DISABLE_ROW_TASK_FLOOR`
- `KLS_DISABLE_ROW_LIGHT_SNODE`
- `KLS_DISABLE_GENERIC_HOISTED_SNODE_WORKER`
- `KLS_DISABLE_LEAN_PACKED_INPUT`
- `KLS_ENABLE_DIRTY_ROW_TRANSPOSE_PLAN_SOLVE`
- `KLS_DISABLE_DIRTY_ROW_TRANSPOSE_PLAN_SOLVE`
- `KLS_ENABLE_SETTLED_COMPACT_MAP32_REFACTOR`
- `KLS_DISABLE_LOW_WORK_SINGLE_BLOCK_COMPACT_MAP32_POLICY`
- `KLS_DISABLE_COMPACT_MAP32_TOURNAMENT`
- `KLS_DISABLE_COMPACT_MAP32_VENDOR_PREFLIGHT`
- `KLS_DISABLE_COMPACT_MAP32_EARLY_VERDICT`

Keep default row task-floor admission, supernode construction and worker
eligibility, packed-input profitability and width checks, transpose-plan
admission and fallback, and compact-map selection and numerical checks.
The removed overrides are no longer supported; the public C API is unchanged.

Removing the overrides also makes the vendor-first tournament unreachable:
on reaching the trial, the repeated-update contract must hold and the
choice must be zero. A true low-work predicate would already have set the
choice positive and returned through the settled compact branch. Nothing
between these checks changes the predicate's inputs. Therefore all remaining
trials use compact-first ordering. Remove the dead vendor preflight and
early-verdict blocks, their thresholds, and ordering alternatives. Retain
the live six-sample comparison, state-2 residual certificate, median-of-three
verdict and expected-horizon payback checks.

Net reduction: 104 source lines across two files.

## Validation

Final Release and ASan builds pass without reported warnings. Release CTest
passes 6/6 in 4.09 s; ASan passes 6/6 in 17.18 s. `git diff --check` passes.
Executable text decreases by 1,512 bytes; data and BSS are unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -1.479% to +0.594% versus accepted and from
-3.242% to +0.458% versus cumulative baseline `372f9d1`. No case reaches
the 2.5% review threshold. No alignment fix is needed.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed throughout timing. These focused results do not
establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/row-map-control-trim-*`.
Runner: `validate-row-map-control-trim.py`. Audit: `audit-focused-trim.py`.
