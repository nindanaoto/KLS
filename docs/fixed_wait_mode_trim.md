# Remove constant worker-wait fields

Baseline: `6f539b9`. Both refactor-pool `busy_wait` fields are initialized
to one and never changed after the preceding override removal. Delete both
fields, their initializations, two worker-loop checks and the caller-side
conditional around its existing bounded completion spin.

Keep the bounded spin loops, adaptive EGraph spin counts, generation and
shutdown checks, mutexes, condition-variable sleep fallback, and completion
wait. No synchronization ordering, default threshold or numerical kernel
changes. Net reduction: six source lines and two private fields. Public C
API unchanged.

## Validation

Release and ASan builds pass without reported warnings. Release CTest passes
6/6 in 4.06 s; ASan passes 6/6 in 17.20 s. `git diff --check` passes.
Executable text decreases by 64 bytes; data and BSS are unchanged.

All 579 focused H100 launches pass the provenance/result audit. Paired
lifecycle changes range from -1.345% to +1.531% versus accepted and from
-2.969% to +0.409% versus cumulative baseline `372f9d1`. No case crosses
the 2.5% regression review threshold. No alignment change is needed.

Audits verify exact jobs, commands, source revision/diff, build provenance,
binary/matrix hashes, exits, residuals and recomputed medians. Source and
binaries remain fixed throughout timing. These focused results do not
establish full-paper or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/fixed-wait-mode-trim-*`.
Runner: `validate-fixed-wait-mode-trim.py`. Audit: `audit-focused-trim.py`.
