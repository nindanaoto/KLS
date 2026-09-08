# Remove discarded optional-output locals

Preceding baseline: `40307fa`; resumed-pass baseline: `372f9d1`.

Remove nine locals whose only use was receiving helper outputs subsequently
discarded with void casts: four suffix offsets in ragged-batch eligibility,
one phase row count, and four pipeline counters in the parallel block worker.
Pass NULL through the helpers' existing optional-output contracts. Every write
to these outputs is null-guarded; none controls execution. Retain helper outputs
and statistics consumed by other callers, including live suffix offsets used
to construct panels. This removes 20 net C-source lines without changing solver
selection, numerical checks, or arithmetic.

Release and ASan/UBSan builds are clean; CTest passes 6/6 in each. Stripped
executables still differ after removing symbol tables and build ID, so timings
were required. Comparison copies: `build/discarded-output-kEJ4lq/`.
`kls_snb_refactor` remains 64-byte aligned at 0x64540, size 0x1931.

All 579 launches are valid: 36 onetone controls, 495 automatic-policy controls,
and 48 forced-first controls. Artifacts:
`build/prep-trim-repair-dBcoha/discarded-output-trim-*`; runner:
`validate-discarded-output-trim.py`; preceding frozen binary:
`panel-gate-accepted-kls_bench`. The three versions rotate through the same
hash-verified execution path using pinned H100 changed-refactor workloads with
residual validation at 1e-8. No builds or source edits overlap timing.

Audit verifies job coverage, launch counts, runner/binary/matrix hashes, exact
pending source diff, shared path, forced-first arguments, exit status, residual
validity, and paired lifecycle medians against raw records. Incremental changes
range from -0.510% to +1.344%; cumulative changes range from -3.562% to +1.276%.
No control reaches the 2.5% review trigger or approximately 3% stopping rule.
TSOPF_FS_b9_c6 cumulative changes are +0.287%, -0.431%, and +0.584% on
CPUs 8–15, CPUs 0–7, and CPU 8 respectively. Retain this chunk.

The broader low-risk audit remains active. Focused validation does not prove
cross-machine equivalence or replace periodic cumulative corpus validation.
