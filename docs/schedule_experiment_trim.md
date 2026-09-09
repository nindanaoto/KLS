# Remove opt-in scheduling experiments

Baseline: `ed8edf7`. Remove three off-by-default scheduling controls and
the private implementations reachable only through them:

- `KLS_ROW_REFACTOR_ORDER_PRIVATE_WITH_EXTERNAL_WAITS`: remove the mask
  builder, mode-two worker branches, mask storage and cleanup.
- `KLS_ROW_REFACTOR_STEAL_PRIVATE`: remove the opt-in private-group
  redistribution and its thread-partition reset.
- `KLS_ROW_REFACTOR_CALIBRATE_COMPONENTS`: remove per-component timing,
  calibration installation, scale arrays and calibrated work estimates.

No tests, benchmark drivers, scripts or existing documentation reference
these overrides. They are no longer supported. Keep the default private
ordering validation, normal work estimates, dynamic queues, shared live
dependency-wait helper, and numerical/error fallbacks. The public C API
is unchanged. Net reduction: 371 source lines across four files, including
two helpers and six private fields. Alternative work-model overrides remain
for a separate review; the overhead model also has a live default selector.

## Validation

Both builds succeed without reported warnings. Release CTest passes 6/6
in 4.07 s and ASan CTest passes 6/6 in 17.17 s. `git diff --check` passes.
Executable `size` text decreases by 3,512 bytes; data/BSS are unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -1.526% to +1.887% versus accepted and from
-3.054% to +1.921% versus cumulative baseline `372f9d1`. Four modest
slowdowns receive confirmation runs with the same frozen source/binaries:

| Case | Triplets | Versus accepted | Versus cumulative |
| --- | ---: | ---: | ---: |
| TSOPF_FS_b9_c6, eight threads/96 MB | 24 | +0.994% | -1.883% |
| TSOPF_RS_b9_c6, eight threads/32 MB | 24 | -0.596% | -0.183% |
| mimo46x46_system, eight threads/96 MB | 24 | -0.411% | -0.104% |
| twotone, eight threads/32 MB | 9 | +0.997% | +0.157% |

All 243 confirmation launches pass the audit, for 822 valid launches total.
No notable regression is observed; no alignment fix is needed. Small timing
differences are reported rather than treated as exact equivalence.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Focused
results do not establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/schedule-experiment-trim-*` and
`schedule-experiment-confirm-repeat-*`.
Runner: `validate-schedule-experiment-trim.py`. Audits: `audit-focused-trim.py`
and `audit-eligibility-repeat.py` adapted to the four confirmation jobs.
