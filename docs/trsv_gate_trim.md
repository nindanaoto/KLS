# Collapse dead compact-TRSV AUTO gate

Preceding baseline: `c79d758`; resumed-pass baseline: `372f9d1`.

The private `kls_compact_supernode_trsv_auto_allows` helper discarded both
arguments and always returned zero. Its only caller wrapped that result with
null/negative/positive state tests; that wrapper also had only one caller.
Replace both helpers with the equivalent null-safe positive-state predicate.
The explicit TRSV opt-in remains available; AUTO still declines. Workspace
limits, arithmetic, and numerical checks are unchanged. This removes 20 net
C-source lines.

Release and ASan/UBSan builds are clean, and CTest passes 6/6 in each. The smoke
test explicitly enables compact TRSV and checks that its counters record use.
Stripped executables differ after removing symbols and build ID; comparison
copies are under `build/trsv-gate-djrYK9/`. The SNB kernel remains 64-byte
aligned at 0x64540, size 0x1931.

All 651 launches are valid: 36 onetone, 495 automatic-policy, 48 forced-first,
and 72 confirmation launches. Artifacts:
`build/prep-trim-repair-dBcoha/trsv-gate-trim-*` and `trsv-gate-repeat*`;
runner: `validate-trsv-gate-trim.py`; preceding frozen executable:
`reserve-output-accepted-kls_bench`. Three versions rotate through one shared
hash-verified execution path with pinned H100 changed-refactor workloads and
residual validation at 1e-8. No builds or edits overlap timing.

TSOPF_RS_b9_c6 initially flags at +2.961% incrementally and +2.090% cumulatively
over twelve repetitions. An independent twenty-four-repetition confirmation
measures +1.100% / +1.357%, with no failures. The initial near-3% slowdown does
not persist. No other initial control reaches the 2.5% review trigger.
TSOPF_FS_b9_c6 cumulative results are -0.380%, +0.003%, and +0.046% on
CPUs 8–15, CPUs 0–7, and CPU 8 respectively.

Audit verifies expected jobs/counts, runner/binary/matrix hashes, exact pending
source diff, shared path, forced-first arguments, successful exits, residual
validity, and paired medians from raw records. Retain the chunk: the repeat
does not establish an approximately 3% stopping-rule regression. This is not
proof of zero loss or cross-machine equivalence.

Four chunks have now removed 63 net C-source lines since the last cumulative
corpus checkpoint. Repeat that shared-path gate before further deletions.
The broader low-risk cleanup goal remains active.
