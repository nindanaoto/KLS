# Remove unused panel-gate plumbing and RCM alias

Preceding baseline: `168a415` (source unchanged by report commit `428c475`);
resumed-pass baseline: `372f9d1`.

The private native-panel eligibility wrapper forwarded `trailing_len` to a
helper that discarded it. That helper also repeated the negative environment
check already enforced by `kls_native_row_panel_enabled_for_run`. Remove the
unused argument and helper, retaining its width check in the wrapper. Panel
selection, automatic-disable behavior, pivot checks, and arithmetic are unchanged.
Also remove the unused `rcm_queue` alias and its void cast; the actual RCM queue
continues to use the permutation scratch. This removes 13 net C-source lines.

Release and ASan/UBSan builds complete cleanly and CTest passes 6/6 in each.
Binary copies under `build/panel-gate-AtDu7N/` differ even after removing debug
information, symbol tables, and build ID, so timing validation was required.
`kls_snb_refactor` remains 64-byte aligned at 0x64540, with size 0x1931.

All 579 focused launches are valid: 36 onetone controls, 495 automatic-policy
controls, and 48 forced-first controls. Artifacts are under
`build/prep-trim-repair-dBcoha/panel-gate-trim-*`; runner:
`validate-panel-gate-trim.py`; preceding frozen executable:
`compact-shape-accepted-kls_bench`. Three versions rotate through a shared
execution path, copied and hash-verified for each launch, using the established
pinned H100 changed-refactor workload with residual validation at 1e-8.
No builds or source changes overlap timing.

Audit verifies the expected jobs and launch counts, runner/binary/matrix hashes,
exact pending source diff, shared execution path, forced-first arguments,
successful exits, residual validity, and paired medians against raw records.
Lifecycle changes range from -4.001% to +1.209% incrementally and -2.448% to
+2.203% cumulatively. The largest cumulative slowdown is LeGresley_2508;
it remains below the 2.5% review trigger. TSOPF_FS_b9_c6 measures +0.210%,
+0.104%, and +0.258% cumulatively on CPUs 8–15, CPUs 0–7, and CPU 8.
Neither the review trigger nor the approximately 3% stopping rule is reached.

Retain this chunk. The preceding cumulative medium-corpus screen is documented
in `compact_equal_corpus_screen.md`; this focused screen does not replace a
future cumulative gate or prove cross-machine performance equivalence. The
broader low-risk cleanup audit remains active.
