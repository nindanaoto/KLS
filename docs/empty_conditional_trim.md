# Remove empty conditional remnants

Baseline: `d559d5f`. Net reduction: 10 C-source lines across three files.

Remove two empty first-factor conditionals and an inert ordering threshold
check. Replace the empty schedule-build conditional with an explicit discarded
return value; the schedule builder still executes and elapsed accounting is
unchanged. Live predicates, numeric policies and alignment hints remain intact.

Warning-free Release and ASan builds pass CTest 6/6 each (3.68 and 15.45
seconds). Complete Release executables compare byte-identical after applying
objcopy --strip-debug --remove-section=.note.gnu.build-id to both versions.
Runtime code, data and layout are unchanged in this build, so no separate
timing rerun is needed for this chunk. This does not establish equivalence
for arbitrary other compilers or configurations.

Original SHA-256: dc5f0246478e2cd5da333f592cf557c3117a9c0e6c15ec3d26b32a92cc233357.
Rebuilt SHA-256: b0f708cc34b4bf1c0b89b954dbd5759d0253fbde49fe57ef1ba64d803fecd87d.
Artifacts: `build/prep-trim-repair-dBcoha/empty-block-*` (frozen baseline and
before/after stripped executables).

The broader cleanup goal remains active.
