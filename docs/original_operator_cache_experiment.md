# Original-operator structural cache experiment

Status: unqualified experiment, compared against `adcc9f0`.

The current qualified-competitor profiles identify certificate checking as a
large solver cost and original-operator reconstruction as a prominent BJT stack
sample. This experiment targets the latter without changing certificate math,
numerical kernels, correction budgets, or the accuracy policy.

The cache retains the validated internal pattern/permutations, orientation,
and a caller-value gather map. Every hit requires matching dimensions, map
presence, orientation, and exact array contents; pointer identity is not enough.
Changed layouts take the original fully validated reconstruction path. Fresh
values are gathered into a new allocation and published only after successful
preparation. An invalid update or allocation failure cannot partially overwrite
the existing snapshot. Optional cache allocation failure falls back to the
uncached implementation. Clearing the operator releases the cache.

The extra retained storage is up to `(3*n + 2*nnz + 1)*sizeof(UF_long)` plus
`nnz*sizeof(int64_t)`. Each cache hit still compares structural data and allocates
a fresh value array. Whether this tradeoff improves performance must be measured.

`tests/kls_original_operator_cache.c` compares cache updates against full rebuilds
across orientations, optional permutations, in-place structural/map edits,
value updates, nonfinite application-managed inputs, and allocation failures.
The standalone ASan/UBSan test passed. Full solver tests and performance gates
are orchestrated by `/tmp/kls-refinement-first/qualify_operator_cache_v2.py`, with
evidence under `/tmp/kls-refinement-first/operator-cache-v2*`. The first attempt
stopped at the ownership test's assumption that every update made six required
allocations. The revised test separately injects the cached value-allocation
failure and all six forced-rebuild failures, retaining the same invalid-epoch
and unchanged-snapshot assertions. The failed attempt is preserved under
`operator-cache/`; it did not reach performance testing.

The revised candidate passed all 37 CTests, including the expanded ownership
fault coverage. Nine fresh BJT-8 pairs showed only about 0.9% median total-time
improvement, below the 2% incremental-benefit criterion by itself. The full
seven-circuit performance comparison is running; no retention or speed-leader
claim is supported yet.

The three-block screening gate passed: MOS13 improved 2.68%/3.19% at
one/eight threads, with no configuration exceeding the 2% regression limit.
Nine independent confirmation blocks are now running. Exact qualification is
queued by `/tmp/kls-refinement-first/queue_operator_accuracy.py`, gated on a
passing confirmation and completion of the serialized worker-prototype tests.
The audit uses the frozen cache-only candidate, not the newer working tree;
archive reuse does not replace fresh exact audits. No gate is considered passed
merely because its execution has been queued.

Promotion requires the full seven-circuit, one/eight-thread no-regression and
incremental-benefit gates, followed by fresh exact accuracy qualification.
This experiment alone does not establish that KLS is fastest. Worker reuse for
certificate checking remains a separate measured optimization target.

The completed campaign contains 343 runs (336 paired-comparison timings and
seven references). Both screening and nine-block confirmation passed the
configured no-regression and benefit gates. Confirmation median time reductions
at one/eight threads were: BJT 0.32%/0.77%, MOS13 2.14%/3.21%, mem_plus
1.54%/1.69%, mux8 1.23%/1.92%, Rallpack3 1.84%/2.91%, and RAM2k
1.66%/1.65%. DAC improved 1.29% at one thread but was 0.33% slower at eight
threads, with substantial paired variation. These are incremental comparisons
against the qualified KLS baseline, not competitor wins. Authoritative results:
`/tmp/kls-refinement-first/operator-cache-v2-performance/confirmation-summary.json`.
Exact accuracy qualification is still pending; passing timing gates alone does
not authorize promotion. The separate worker-reuse prototype has entered its
serialized build/test stage.

The first combined Release test run exposed a test-harness defect: assertions
were disabled in the standalone cache fixture, removing its setup calls and
causing a null-pointer access. Keeping assertions enabled in that fixture fixes
the test without changing production code. The revised combined tree passed
all 40 CTests, the worker-plan ASan/UBSan check, and the GCC/OpenMP exact-bound
checks (`/tmp/kls-refinement-first/worker-reuse-tests-v2/status.json`). The first
failed test run remains preserved. The cache-only accuracy audit was then
started directly with `qualify_operator_accuracy.py`; the original queue had
correctly stopped on the test failure and was not reused.

Final cache-only accuracy qualification passed all six configurations with
complete coverage: 65,674 certified successes and two explicitly reported
MOS13 rejections (one per thread configuration). Both MOS13 simulations
recovered with valid reference waveforms, maximum absolute difference
7.99999622103087e-10. BJT and mem_plus waveforms matched exactly, with no
rejections. Authoritative final evidence is
`/tmp/kls-refinement-first/operator-cache-v2-accuracy/status.json` and its six
per-configuration audit reports. Together with the completed performance gate,
this qualifies the cache-only candidate as an incremental improvement, not as
the fastest solver across the current corpus.
