# Remove unused worker diagnostic alternatives

Baseline: `988b64a`. Remove eight environment overrides not referenced by
tests, benchmark drivers, scripts, or existing documentation:

- `KLS_DISABLE_SCALED_ROW_VALUES_SPECIALIZATION`
- `KLS_DISABLE_SCALED_GROUPED_SPECIALIZATION`
- `KLS_DISABLE_SCALED_FMA_UPDATE`
- `KLS_DISABLE_SCALED_SCALE_SPECIALIZATION`
- `KLS_DISABLE_SCALED_BRANCHLESS_UPDATE`
- `KLS_DISABLE_ROW_SCALE_HOIST`
- `KLS_ROW_SCALAR_SIMPLE`
- `KLS_ROW_TELEMETRY_SHARED`

Remove four diagnostic-only worker variants, the now-constant fused-update
argument and nonfused arithmetic branch, three private fields, and alternate
scalar/scaling/telemetry routing. Retain the default numerical workers,
representation checks, sum/max scaling, scalar fallback, worker-local counters
and their public statistics. The general worker still needs its branchless
argument because surviving callers use both modes. These environment
overrides are no longer supported; the public C API is unchanged.

Net reduction: 92 source lines. The release executable's `size` text column
falls from 3,085,183 to 3,076,255 bytes (8,928 bytes); data/BSS are unchanged.

## Validation

Both builds succeed without reported warnings. After build completion,
release CTest passes 6/6 in 4.14 s and ASan CTest passes 6/6 in 17.27 s.
`git diff --check` passes. Existing smoke tests include scaled serial and
parallel row refactoring and assertions on retained telemetry.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -0.567% to +2.032% versus accepted and from
-3.018% to +1.586% versus cumulative baseline `372f9d1`.
The largest incremental slowdown is first-factor bcircuit on eight threads
with the 32 MB placement. A further 24 rotated triplets (72 launches)
measure -0.269% versus accepted and -0.371% versus cumulative. All 651
launches are valid; no persistent notable regression is observed. No
alignment fix is needed.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals, and recomputed medians. Source and
binaries stay fixed during timing. This focused campaign does not establish
full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/worker-control-trim-*` and
`worker-control-confirm-first-*`. Runner: `validate-worker-control-trim.py`.
Audits: `audit-focused-trim.py` and `audit-eligibility-repeat.py` adapted to
the single first-factor bcircuit job.
