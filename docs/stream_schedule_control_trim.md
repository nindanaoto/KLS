# Remove unused stream and schedule controls

Baseline: `9d8c5ea`. Remove eight controls unused by tests, scripts and
benchmark drivers:

- `KLS_DISABLE_PARALLEL_SORT`
- `KLS_CLUSTER_WIDTH_ALPHA`
- `KLS_DISABLE_DENSE_TAIL_REFACTOR`
- `KLS_DISABLE_EGRAPH_CLUSTER_CUT_TRIAL`
- `KLS_DISABLE_EGRAPH_STREAM_KERNEL_TRIAL`
- `KLS_EXPERIMENT_SNODE_TAIL_CHUNK128`
- `KLS_EXPERIMENT_SNODE_TAIL_CHUNK144`
- `KLS_EXPERIMENT_SNODE_TAIL_MASKED_REMAINDER`

Keep default sorting admission/fallback, the 2.5-worker incumbent cut and
alpha-3/alpha-4 trial plans, dense-tail refresh/fallback and measured stream
selection. After removing forced tail experiments, all three private tail
flags are always reset together and assigned the same applicability value.
Replace them with one `snode_wide_tail` flag; retain all tail kernels and
their original length checks. Net reduction: 37 source lines and two private
fields across three files. Public API unchanged; removed overrides are no
longer supported. Update the historical cluster-width note accordingly.

Test-used SNB controls, the EGraph floor hook and the paper configuration's
panel-refactor disable control remain intact.

## Validation

Release and ASan builds pass without reported warnings. Release CTest passes
6/6 in 4.13 s; ASan passes 6/6 in 17.16 s. `git diff --check` passes.
Executable text decreases by 832 bytes; data and BSS are unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -0.779% to +2.302% versus accepted and from
-5.144% to +0.681% versus cumulative baseline `372f9d1`.

Repeat both eight-thread TSOPF_FS configurations because private-state
layout changed, and forced-first-factor 1138_bus because its +2.302% result
was near the 2.5% review threshold. All 216 additional launches pass audits:

| Case | Configuration | Versus accepted | Versus cumulative |
| --- | --- | ---: | ---: |
| TSOPF_FS_b9_c6 | t8_32 | +1.210% | -2.206% |
| TSOPF_FS_b9_c6 | t8_96 | +1.168% | -1.477% |
| 1138_bus, forced first factor | t8_32 | +1.180% | -0.402% |

The repeats show a small approximately 1.2% cost, not proof of identical
performance. It remains below the 2.5% review threshold; no alignment change
is made. All repeated cases remain faster than the cumulative baseline.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed throughout timing; documentation changes follow audits.
These focused results do not establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/stream-schedule-control-trim-*`
and `stream-schedule-control-confirm-*`. Runner:
`validate-stream-schedule-control-trim.py`. Audits: `audit-focused-trim.py`
and `audit-eligibility-repeat.py`, adapting the latter's first-factor job
and count to 1138_bus / t8_32 / 24 triplets / 72 launches.
