# Remove compact-solve diagnostic scans

Baseline: `a610e94`. Remove the 187-line `KLS_TRACE_RUNS` block from compact
solve preparation, including stored-zero/tiny-value counts, contiguous-run
counts, the nested `KLS_TRACE_TAIL` density/flop profile and border-density
scan, and hypothetical gap-padding simulations. All counts and temporary
storage were local and used only for printed diagnostics.

The retained i32 solve mirror, RHS mapping, singleton-run preparation,
numeric kernels, factor values, and dispatch policies are unchanged.
The two retired environment variables no longer affect execution; no
source/test/script references remain.

Release and ASan/UBSan CTest each pass 6/6. The timed comparison is sequential
and does not overlap builds or tests. Artifacts live under
`build/prep-trim-repair-dBcoha/compact-census*`, recording source diff,
commands, binary/matrix hashes, checked build provenance and phase timings.
H100 runs use automatic policy and entrywise changed values, verifying
every changed-system residual at 1e-8.

Three executables rotate order: `original` is `a610e94` (the immediately
preceding chunk), `start` is pushed `88a177e` (the start of this cleanup
sequence), and `fix` is this candidate. Both incremental and cumulative
comparisons matter for the stopping condition. A repeatable lifecycle loss
of roughly 3% or more, or any correctness failure, triggers investigation
before further slimming. The focused comparison is not a full paper or
Xyce rerun.

All **282/282 launches pass** across seven matrices and ten configurations.
No notable slowdown was observed. TSOPF eight-thread lifecycle changes are
+0.14% / -0.36% versus the preceding commit and +0.21% / -0.78% versus the
pushed starting point (CPUs 8-15 / 0-7). The largest incremental slowdown is
onetone1 at +1.00%, which is only +0.36% cumulatively. Other cumulative
changes range from -2.68% to +0.62%. Small control samples are not a proof
of identical performance for all inputs.

The source-wide follow-up audit found another bounded candidate: the
`kls_snode_trace_*` batch counters have only definitions and updates, with
no consumers in any source file. They will be handled in a separate chunk,
leaving the numeric batching algorithms and general trace output intact.
