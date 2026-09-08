# Remove constant compact-update shape flag

Preceding baseline: `7103c5b`; resumed-pass baseline: `372f9d1`.
Both callers of the private supernode-run wrapper passed `shape_known = 1`;
the wrapper forwarded it unchanged to the compact helper. The other compact
caller also passed 1. There are no indirect references. Remove the flag
from declarations, definitions and calls, and remove the never-entered
`!shape_known` validation branch. Add a comment documenting the caller
contract. This removes 20 net C-source lines.

Active bounds, pivot, row-length, workspace and numerical checks remain.
Callers still identify runs through the existing supernode machinery; no
previously executed shape check is bypassed by this change. Cached/compact/
scalar dispatch and arithmetic are unchanged.

Release and ASan/UBSan builds are clean; CTest passes 6/6 in each. The
release executable differs after stripping debug information, symbol tables
and build ID, so timing validation was required. Binary-comparison copies
are under `build/compact-shape-4EMQik/`.

All 579 timing launches pass: 36 onetone controls, 495 automatic-policy
controls and 48 forced-first controls. Artifacts:
`build/prep-trim-repair-dBcoha/compact-shape-trim-*`; runner:
`validate-compact-shape-trim.py`; preceding frozen binary:
`compact-publish-accepted-kls_bench`. Three versions rotate through pinned
H100 lifecycles using a shared temporary executable path with per-launch
copy/hash verification. Every changed refactor is checked at 1e-8. Commands,
build provenance, hashes and exact source diffs are checked or recorded;
no builds/tests overlap timing.

Median paired lifecycle changes span -1.18% to +0.77% incrementally and
-2.84% to +0.93% cumulatively. Onetone1 measures +0.50% / +0.11% on CPU 8,
and +0.77% / +0.83% on CPU 0. TSOPF_FS_b9_c6 measures -0.19% / +0.62% on
CPUs 8-15 and +0.34% / -0.48% on CPUs 0-7. No case reaches the 2.5% review
trigger or approximately 3% stopping rule. Retain this chunk.

The last four chunks removed 175 net C-source lines since the corrected
medium-corpus checkpoint. Repeat that shared-path cumulative screen before
further deletions; focused validation is not corpus-wide or cross-machine
proof. The broader cleanup goal remains active.
