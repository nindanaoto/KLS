# Remove fixed private mode arguments

Baseline: `108b702`. Remove six arguments whose complete caller sets,
including macro-generated wrappers, pass one:

- `balanced_target_sum` in the lean supernode helper: retain balanced
  arithmetic and remove the unreachable sequential accumulation branch.
- `scaled_input_mode` in the scaled hoisted worker: retain scaling and remove
  the unreachable unscaled input branch.
- `wait_for_dependencies` in cached pair/multi consumers: retain checks.
- `check_pivots` in parallel tail refactor: always request pivot checks.
- `allow_pivot` in partial first-factor finish: retain pivot handling.

PTS `forward_only` was a scan false positive: the second caller passes one
and contains a cast that the preliminary scanner filtered out. Retain both
PTS modes. No public interface changes. Net reduction: 80 source lines.

## Validation

An initial build caught an incomplete removal of the unreachable unscaled
branch; it was corrected before validation. The final release build and ASan
build succeed without reported warnings. Release CTest passes 6/6 in 4.06 s;
ASan passes 6/6 in 17.25 s. `git diff --check` passes.

Generated code differs, so the final source receives the full focused H100
comparison. All 579 launches pass the provenance/result audit. Median paired
lifecycle changes range from -1.618% to +1.071% versus accepted, and from
-3.193% to +1.135% versus cumulative baseline `372f9d1`. No case reaches the
2.5% review threshold; no alignment fix is needed.

The audit verifies exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals, and recomputed medians. Source and
binaries remain fixed during timing. These focused results do not establish
full-paper or cross-machine performance equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/fixed-mode-trim-*`.
Runner: `validate-fixed-mode-trim.py`. Audit: `audit-focused-trim.py`.
