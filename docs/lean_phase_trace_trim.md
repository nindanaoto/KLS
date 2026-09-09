# Remove output-only lean-phase diagnostics

Preceding source baseline: `161a14d`; cumulative checkpoint: `f162cf5`;
resumed-pass baseline: `372f9d1`.

Remove `KLS_TRACE_LEAN_PHASES`, its phase clocks, and its two output blocks.
The two value-map setup results were used only in logging; replace those
locals with conditional calls preserving their predicates, evaluation order,
and side effects. Value-map construction, numeric work, and policy decisions
remain intact. No paper runner or test consumes this diagnostic output.
The cleanup removes 27 net C-source lines.

Release and ASan builds succeed without the unused-variable warnings from
the intermediate edit. Final CTest passes 6/6 in both builds (3.74 and
15.72 seconds). The SNB kernel remains 64-byte aligned at 0x641c0, size
0x1931; no placement workaround was added.

All 579 focused launches pass: 36 onetone, 495 AUTO, and 48 forced-first
controls. Three frozen versions rotate through one hash-verified execution
path, with pinned H100 entrywise-changed refactors and residual validation
at 1e-8. No source edits or builds overlap timing.

No case reaches the positive 2.5% review trigger. Incremental paired-median
lifecycle changes range from -1.589% to +1.013%; cumulative changes range
from -2.694% to +1.804%. TSOPF_FS_b9_c6 incremental changes are -0.114%,
-0.011%, and -0.072% on CPUs 8–15, CPUs 0–7, and CPU 8; cumulative
changes are -0.821%, +0.194%, and -0.143%. These results do not establish
zero loss or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/lean-phase-trim-*`; runner:
`validate-lean-phase-trim.py`; preceding frozen executable:
`row-phase-accepted-kls_bench`. The final audit verifies all expected
jobs/counts, exact pending source diff, revision, build provenance,
runner/binary/matrix hashes, complete commands, successful exits, residual
validity, and every reported median recomputed from raw records.

This is the first source chunk after the row-phase cumulative checkpoint.
The broader low-risk cleanup goal remains active.
