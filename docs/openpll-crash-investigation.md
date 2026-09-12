# OpenPLL extracted-DCO crash investigation (2026-09-12)

Status: one uninitialized-metadata defect is reproduced and fixed; the previously
observed segmentation fault is not yet attributed conclusively to that defect.
Do not describe the crash as fixed solely from passing reruns.

## Matched KLU control

The same Xyce executable and extracted DCO deck were run with KLU instead of KLS:
V3 candidate, SS, 1.62 V, 125 C, coarse 1 / fine 8, one rank, 30 ns,
RELTOL=1e-5, ABSTOL=1e-9, maximum step 10 ps. KLU completed with 14,127 linear
solves and zero failed linear solves. Its five-period mean over 10–30 ns was
276.647176 MHz; a completed KLS debugger run gave 276.648277 MHz (about 4 ppm
apart). These short traces are diagnostic, not numerical-convergence signoff.
Concurrent diagnostic execution makes their elapsed times unsuitable as a
performance comparison.

## Reproduced metadata defect

Valgrind on the isolated OpenPLL matrix/RHS replay reports an uninitialized read
on the first factorization in `kls_compute_nicslu_task_flow_model`, originating
from `Numeric->Ulen` allocated by `trilinos_klu_l_factor`.

The bundled native KLU factorizer does not populate `Lip`, `Uip`, `Llen`, `Ulen`
and `LUsize` for singleton BTF blocks: those blocks have no stored strict L/U
entries, and KLU's block kernels skip these fields. KLS's global metadata
consumers inspect them. Heap contents therefore leak into KLS's diagnostics and
representation decisions.

The bundled factorizer now explicitly sets that empty metadata to zero. The
`kls_singleton_metadata` regression uses a poison allocator and a matrix mixing
a two-by-two block with two singleton blocks. It deterministically fails before
the change and passes afterward, while checking solves and value refactoring.
The KLS smoke and transformed-residual tests also pass. A 500-call OpenPLL replay
under Valgrind after the fix reports zero errors and no leaks.

## Segmentation fault still under investigation

The original coarse-1/fine-8 crash reads an invalid address at the scalar
`x[li[p]]` load in `kls_parallel_refactor_block`. The faulting instruction is
`vmulsd`, not an aligned SIMD load. A native factor row index or its source
storage is invalid at that point; this alone does not identify who corrupted it.

Disabling strict aliasing did not prevent the original failure. Subsequent
30 ns debugger reruns (including runs preserving address randomization) have
completed. Public factor/refactor/solve boundary guards check native LU offsets,
triangular row indices and cached index consistency. The captured 2,048-call
prefix also replays successfully. These negative reproduction results do not
prove that the crash is repaired.

Full-process memory checking and a native-run stop-on-fault diagnostic are being
used to distinguish internal KLS corruption from corruption originating in the
Xyce integration or device evaluation. Fresh numerical factorization remains a
previously tested diagnostic workaround, not an established root-cause fix.

## Artifacts

Workspace artifacts: `../build/openpll-verification/kls-root-cause/` relative to
this repository. Important files include:

- `manifest.json`, `c1-klu.log`, `matched-frequency.json`;
- `replay-valgrind.txt` (original uninitialized read);
- `test-singleton-before.log`, `test-singleton-after.log`, `tests.log`;
- `replay-valgrind-fix1.txt` (500 calls, zero memory errors);
- `c1-calls.gz`, `trace_wrappers.inc`, `replay.c` (local diagnostic format);
- guard, trap and device-sanitizer build scripts and logs.

The trace format is a local, native-ABI diagnostic format, not a supported public
file format. It retains a capped prefix of calls, not the full 30 ns run.


## Follow-up reproduction assessment

The crash is not currently reproducible on demand. The saved failing noalias
executable completed three 30 ns native runs with a fault hook and another
30 ns run traced by its parent without any injected library. The saved standard
executable (`Xyce-before-fresh-factor`, SHA256
`0f10c51b8e8a130598588ff6ea3c15aa13090650368e223a90ba7703ca3beb98`)
also completed the historical V3 coarse-0/fine-247 and V2 coarse-0/fine-8 decks.
These are passing runs of unfixed executables, so they cannot establish that
the singleton initialization repaired the original crash.

The saved noalias link command includes the noalias KLS and vendor archives;
the suspected mixed strict-aliasing linkage is not present. Of 2,429 saved
input hashes, the circuit, PDK and executable still match. Only two Python
analysis scripts differ, neither of which executes during direct saved-deck
reruns. The original output filename and stdout pipe were not reproduced
exactly by these runs; their possible effect on heap layout remains untested.

Additional diagnostic results:

- A 2,048-call original-vendor replay passes guards for disjoint packed L/U
  storage, native/cached row indices and input maps checked against the current
  permutation. Maximum relative solution difference is approximately 1.07e-13.
- Initial internal state has supernode batching and compact 32-bit solve
  indices; no moderate-BTF tournament or row-engine adoption was observed.
  With one thread, model-based row preparation exits before consuming the
  task-flow recommendation affected by the metadata defect.
- The BSIM4 B4/B4p61 AddressSanitizer run completed 30 ns. This instruments
  those device objects, not the entire application.
- A diagnostic protecting interior pages of numerical-factor allocations
  between KLS API calls completed 30 ns. No external write into those protected
  pages was observed; this does not protect writes occurring inside KLS calls.
- The original 500-call Valgrind replay with only the known task-flow metadata
  warnings suppressed had zero additional errors (659 suppressed reports in
  two contexts).

`run_ptrace.py` captures the fault before MPI's signal handler, retains raw
registers and siginfo, and generates a core using the target's parent PID for
debugger attachment. Its intentional-fault self-test successfully captured a
core; that self-test is not a reproduced KLS failure. All completed native
capture runs above finished without a fault. Full-process Valgrind was deliberately stopped at approximately 15.30 ns
when PLL development resumed. No additional error had been reported; this is
an interrupted diagnostic, not a passing 30 ns test. The crash investigation
is parked with its capture tools and evidence retained.
