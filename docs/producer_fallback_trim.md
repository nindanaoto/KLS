# Remove unreachable scalar producer fallback

Baseline: `e6b9840` (numeric source `bc859b4`), following the completed
medium-corpus checkpoint. Remove seven net C-source lines in
`kls_first_factor.inc`: the saved `active_target_count`, its always-true
comparison with unchanged `target_count`, and the unreachable in-lock
fallback update. The existing detached-update block remains scoped.

Between the removed assignment and comparison, nothing modifies or takes
the address of `target_count`. Trace updates affect separate counters;
admission failures return before the comparison. Thus the detached path
was always selected. All allocation checks, copies, ownership gates,
mutex operations, update calls and cleanup on that path are unchanged.

Release and ASan/UBSan CTest each pass 6/6. The release benchmark's `.text`
is byte-identical to the frozen executable used by the completed corpus
screen. More strongly, the entire executables compare equal after
`objcopy --strip-debug --remove-section .note.gnu.build-id` on separate
output copies. Their common normalized SHA-256 is
`ae66f1177ea3b9006892f8d04752829e9f88dafa223f5e9a8911137e6748f01f`.
The common `.text` SHA-256 is
`8049bf6c2f54d4468e212b6155104fd4347bc74ef40ee3d7bd1deec3fd55ccd9`.

Artifacts are under `build/producer-fallback-iDiTCf/`; the reference is
`build/resumed-corpus-chhiOV/screen0-kls_bench`. This establishes unchanged
release runtime code/data/layout for this build, so no new timing campaign
is needed for this source-only removal. It is not a cross-compiler claim.
Retain the previous corpus result and its five inconclusive cases; the
broader cleanup goal remains active.
