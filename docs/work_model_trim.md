# Remove override-only separator work models and leaf policies

Baseline: `a86a7e3`. Remove six environment controls not referenced by
tests, benchmark drivers, scripts or existing documentation:

- `KLS_ROW_REFACTOR_SEPARATOR_WORK_MODEL`
- `KLS_ROW_REFACTOR_GROUP_OVERHEAD`
- `KLS_ROW_REFACTOR_SEPARATOR_BETA`
- `KLS_ROW_REFACTOR_SPLIT_PAST_LEAF`
- `KLS_ROW_REFACTOR_PROMOTE_OVERWEIGHT_LEAF`
- `KLS_ROW_REFACTOR_MAX_LEAF_PROMOTIONS`

Remove the rows/square-root/cube-root work models, optional leaf-skipping and
promotion branches, their counters, and override parsing. Keep both default
work models: ordinary symbolic work and the structurally selected
dense-spiked work plus 2000 estimate. Keep the existing balance constant,
default leaf stopping behavior, topology validation and allocation fallback.
The private model selector becomes a boolean. These overrides are no longer
supported; the public C API is unchanged. Net reduction: 119 source lines.

## Validation

Both builds succeed without reported warnings. Release CTest passes 6/6
in 4.11 s; ASan CTest passes 6/6 in 17.23 s. `git diff --check` passes.
Executable `size` text decreases by 2,282 bytes and data by 16 bytes;
BSS is unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -2.797% to +0.686% versus accepted and from
-2.636% to +0.842% versus cumulative baseline `372f9d1`. No case reaches
the 2.5% review threshold. No alignment fix is needed. The twotone gain is
from three triplets and is not presented as a demonstrated optimization.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed during timing. Focused results do not establish
full-paper or cross-machine equivalence; in particular, retaining the
dense-spiked default is supported by source-path equivalence, not a claim
that this focused corpus exhaustively exercises that structural class.

Artifacts: `build/prep-trim-repair-dBcoha/work-model-trim-*`.
Runner: `validate-work-model-trim.py`. Audit: `audit-focused-trim.py`.
