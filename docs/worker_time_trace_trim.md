# Remove output-only row-worker timing diagnostics

Preceding baseline: `b2a95e5`; resumed-pass baseline: `372f9d1`.

Remove `KLS_TRACE_ROW_WORKER_TIMES` output, diagnostic-only group/task/
component clocks, the local plan census, and eight private worker timing
fields plus their shared enable flag. This removes 142 C-source lines.
Live component-calibration sampling and installation, scheduling, numerical
work, and public numerical-work counters remain intact. The component
lookup helper remains because it has non-diagnostic callers. No paper
runner or test consumes the removed output.

Release and ASan builds succeed; CTest passes 6/6 in each (3.75 and
15.59 seconds). The SNB kernel remains 64-byte aligned at 0x641c0, size
0x1931. The private struct layout changes are included in performance
validation; no padding or placement workaround was added.

All 579 focused launches pass: 36 onetone, 495 AUTO, and 48 forced-first
controls. Three frozen versions rotate through one hash-verified execution
path, with pinned H100 entrywise-changed refactors and residual validation
at 1e-8. No source edits or builds overlap timing.

No case reaches the positive 2.5% review trigger. Incremental paired-median
lifecycle changes range from -2.948% to +1.152%; cumulative changes range
from -3.046% to +1.361%. TSOPF_FS_b9_c6 incremental changes are +0.530%,
-0.441%, and -0.114% on CPUs 8–15, CPUs 0–7, and CPU 8; cumulative
changes are -0.622%, -0.831%, and -0.194%. These measurements do not
establish zero loss or cross-machine equivalence.

Artifacts: `build/prep-trim-repair-dBcoha/worker-time-trim-*`; runner:
`validate-worker-time-trim.py`; preceding frozen executable:
`lean-phase-accepted-kls_bench`. The final audit verifies expected
jobs/counts, exact pending source diff, revision, build provenance,
runner/binary/matrix hashes, complete commands, successful exits, residual
validity, and every reported median recomputed from raw records.

This is the second source chunk after cumulative checkpoint `f162cf5`,
bringing the reduction since that checkpoint to 169 C-source lines.
The broader low-risk cleanup goal remains active.
