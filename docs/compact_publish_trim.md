# Remove unreachable compact-panel publication

Preceding baseline: `3ca9d25`; resumed-pass baseline: `372f9d1`.
The private compact supernode update had two callers: the producer passed
`publish_panel = 0`, and the ordinary wrapper passed a local initialized to
zero and never modified. No other source reference or indirect use exists.
Remove the unreachable publication branch, flag, wrapper local and now-unused
statistics parameter. This removes 11 net C-source lines. The separate
live panel-cache publication path and all compact arithmetic remain unchanged.

Release and ASan/UBSan builds are clean; CTest passes 6/6 in each.
The generated release code differs even after removing debug information,
symbol tables and build ID, so performance neutrality was not assumed from
the unreachable branch alone. Comparison artifacts are in
`build/compact-publish-Vo8i2Y/`. The SNB entry remains 64-byte aligned
(0x64c80, size 0x1931).

All 579 timed launches pass: 36 onetone controls, 495 automatic-policy
controls and 48 forced-first controls. Runner and artifacts are
`build/prep-trim-repair-dBcoha/validate-compact-publish-trim.py` and
`compact-publish-trim-*`; preceding binary:
`scalar-progress-accepted-kls_bench`. All three binaries execute through
one shared temporary path after per-launch copy/hash verification. Pinned,
rotating H100 lifecycles check every changed refactor at 1e-8. Provenance,
commands, hashes and exact source diffs are checked or recorded; no builds
or correctness tests overlap timing.

Across tested cases, median paired lifecycle changes range from -2.23% to
+2.03% incrementally, and -1.09% to +1.30% versus the resumed-pass baseline.
The largest incremental increase is `twotone`, CPUs 8-15, at +2.03%, with
-1.09% cumulative change. `TSOPF_FS_b9_c6` on CPUs 8-15 is +0.55% / +1.23%;
onetone1 is +0.27% / +0.20% on CPU 8 and -0.11% / +0.03% on CPU 0.
No case reaches the 2.5% review trigger or approximately 3% stopping rule.
Retain this chunk; the cleanup goal remains active. This focused check is
not a fresh full-corpus or cross-machine validation.
