# Remove numeric-free tracking

Preceding baseline: `93d2751`; cumulative baseline: `372f9d1`.

Remove 108 C-source lines across KLS, the vendored KLU free routine, and
the comparison stub: pointer snapshots, freed-pointer ring, backtraces,
log/check helpers, callers and callback shim. Actual freeing, numerical
work and alignment hints are unchanged. No tests or scripts use the tracker.

This intentionally changes dependency sources: the sole third_party edit
removes the 20-line KLU logging shim. The runner and audit pin the reviewed
dependency-patch SHA-256
`1bcb48f46d4ed4124bae3b4b6a9e078d5300e5cde897d564f947c83050bcb478`.
Provenance explicitly records this patch alongside the original dependency
tree. All other compiler/link/flags, gitlink, external-object and runtime
library evidence must match the baseline. Default paired-build checks remain
fail-closed; these runs use the explicit reviewed-patch opt-in.

Warning-free Release and ASan builds pass CTest 6/6 each (3.70 and 15.43
seconds), including rebuilding KLU comparison executables without the stub.
All 579 focused launches pass residual validation at 1e-8. Incremental
paired-median lifecycle changes range from -2.570% to +3.244%; cumulative
changes range from -2.977% to +1.568%. TSOPF_FS incremental changes are
-2.570%, -1.728%, and +0.022% on CPUs 8–15, CPUs 0–7, and CPU 8.

1138_bus with forced first factor initially flags +3.244% incremental.
An independent 24-triplet repeat passes all 72 launches and measures
-1.223% incremental and -2.058% cumulative. The slowdown does not reproduce;
no alignment intervention is justified. Total valid launches: 651.

Sequential pinned H100 comparisons rotate frozen binaries through one
hash-verified execution path without overlapping builds or source edits.
Both audits verify exact diff, revision, reviewed dependency patch, other
build provenance, hashes, jobs/counts, commands, exits, residuals and medians.
No new full corpus run was made; this is not zero-loss or cross-machine proof.

Artifacts: `build/prep-trim-repair-dBcoha/nfree-trim-*`, including repeats;
runner: `validate-nfree-trim.py`; audits: `audit-nfree-trim.py` and
`audit-nfree-repeat.py`; baseline executable: `nfree-accepted-kls_bench`.

The broader cleanup goal remains active.
