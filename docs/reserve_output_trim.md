# Remove unused reserve arguments and discarded prefix counts

Preceding baseline: `48d3940`; resumed-pass baseline: `372f9d1`.

Remove the unused `used` parameters from the private integer/double panel-cache
reserve helpers and their three calls. Reallocation still uses the same growth,
capacity, and overflow checks. The live cache counts remain in use elsewhere.
Also replace three discarded prefix-count outputs in the prefix-validation
caller with NULL; their helper already guards all output writes. The repair
caller still receives and uses those counts. This removes 10 net C-source lines
without changing numerical checks, allocation sizes, or solver choices.

Release and ASan/UBSan builds are clean; CTest passes 6/6 in each. Executables
differ after stripping all symbols and build ID, requiring timing validation.
Comparison copies: `build/reserve-output-3P28et/`. `kls_snb_refactor` remains
64-byte aligned at 0x64540, size 0x1931.

All 579 focused launches are valid: 36 onetone controls, 495 automatic-policy
controls, and 48 forced-first controls. Artifacts:
`build/prep-trim-repair-dBcoha/reserve-output-trim-*`; runner:
`validate-reserve-output-trim.py`; preceding frozen binary:
`discarded-output-accepted-kls_bench`. Three versions rotate through a shared
hash-verified execution path using pinned H100 changed-refactor workloads,
with residual validation at 1e-8. No builds or edits overlap timing.

Audit confirms expected jobs/counts, runner/binary/matrix hashes, exact pending
source diff, shared execution path, forced-first arguments, successful exits,
residual validity, and paired lifecycle medians against raw records. Changes
range from -2.195% to +1.079% incrementally and -3.712% to +1.106% cumulatively.
TSOPF_FS_b9_c6 cumulative changes are -0.586%, -0.187%, and -0.128% on
CPUs 8–15, CPUs 0–7, and CPU 8 respectively. No control reaches the 2.5%
review trigger or approximately 3% stopping rule. Retain this chunk.

The broader cleanup audit remains active; this focused gate does not establish
cross-machine performance equivalence or exhaustion of low-risk candidates.
