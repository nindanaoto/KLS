# Remove unused EGraph preparation and admission controls

Baseline: `63d1519`. Remove twelve controls unused by tests, scripts,
benchmark drivers and existing documentation:

- `KLS_DISABLE_PARALLEL_PERMUTED_VALUE_PREP`
- `KLS_DISABLE_MODERATE_FRAGMENTED_PARALLEL_VALUE_PREP`
- `KLS_ENABLE_EGRAPH_PARALLEL_SCALE_ROWS`
- `KLS_ENABLE_EGRAPH_PARALLEL_SCALE_PERMUTE`
- `KLS_DISABLE_PLAIN_SINGLE_UNSCALED_EGRAPH`
- `KLS_ENABLE_LEAN_SINGLE_BLOCK_MAP32_REFACTOR`
- `KLS_ENABLE_LEAN_BTF_MAP32_REFACTOR`
- `KLS_ENABLE_EGRAPH_SEPARATOR_PRIVATE`
- `KLS_DISABLE_COMPLETE_SEPARATOR_PRIVATE_SETTLE`
- `KLS_DISABLE_EGRAPH_PIPELINE_FUSION_TRIAL`
- `KLS_DISABLE_RELAXED_EGRAPH_INDEPENDENCE`
- `KLS_EGRAPH_SIZE_FLOOR`

Keep the default preparation/scaling thresholds, persistent-pool reuse,
allocation and serial fallbacks, compact-map eligibility, separator
capability checks and measured choices. Keep default fusion and independence
admission. Replace the parsed size-floor helper with one shared 5000-row
constant at all three call sites. Correct a stale reference to the removed
width trial. The paper-used `KLS_DISABLE_EGRAPH_REFACTOR` remains supported.

The removed overrides, including the historical preparation alias, are no
longer supported. Public C API unchanged. Net reduction: 43 source lines
across three files.

## Validation

The first build exposed two additional callers of the removed size-floor
helper. All three callers now use the shared constant; final Release and
ASan builds pass without reported warnings. Release CTest passes 6/6 in
4.14 s; ASan passes 6/6 in 17.08 s. `git diff --check` passes.
Executable text decreases by 664 bytes; data and BSS are unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -1.833% to +0.628% versus accepted and from
-3.086% to +1.418% versus cumulative baseline `372f9d1`. No case reaches
the 2.5% review threshold. No alignment fix is needed.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed throughout timing. These focused results do not
establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/egraph-admission-control-trim-*`.
Runner: `validate-egraph-admission-control-trim.py`. Audit: `audit-focused-trim.py`.
