# SubtreeLU callstream diagnostics

These optional tools capture public SubtreeLU calls and replay one analysis
generation against the original library. They are not part of KLS or competitive
timing. Never preload the recorder during a performance comparison.

`bench/subtreelu_callstream_capture.c` interposes the four public C++ symbols
for analyze, factorize, refactorize and solve. It delegates original arguments
and returns unchanged. Set `KLS_SUBTREELU_CALLSTREAM_DIR` to an existing directory
and preload the compiled shared object into an isolated diagnostic process.
Each analysis writes an exclusive generation directory containing `initial.csr`
and `calls.bin`. Current coefficients and RHS are flushed before each call;
return codes and initialized successful outputs are flushed after it. A stalled
call therefore retains its complete entry, without a fabricated return.
Recording errors emit `SUBTREELU_STREAM_RECORDING_ERROR` and invalidate capture.
They do not change the delegated solver call. Missing original symbols exit125.

`bench/subtreelu_callstream_replay.cpp` accepts:

```
replay GENERATION_DIR --validate-only
replay GENERATION_DIR --execute-xyce
replay GENERATION_DIR --execute-xyce --include-unreturned
```

Preflight validates complete file extents, CSR, flags, statuses and call order
before constructing a solver. The explicit Xyce preset sets pivot tolerance to
1000 millionths and memory growth to150 percent, matching the audited adapter.
Every recorded return, successful permutation and solution byte must match;
a mismatch stops with exit3. A trailing entry is retained and not invoked
(exit4), unless explicitly requested. If that fresh invocation returns, exit5
identifies a new result without a captured return to match. Use an external
supervisor with a diagnostic time cap and kill/wait the complete process group.

Build with the existing Apptainer dependency image. Inside the image, given the
SubtreeLU checkout at `/work/SubtreeLU`:

```sh
gcc -O2 -std=c11 -fPIC -shared bench/subtreelu_callstream_capture.c \
  -o /tmp/subtreelu_capture.so -ldl -pthread
g++ -O2 -std=c++17 -I/work/SubtreeLU/include \
  bench/subtreelu_callstream_replay.cpp -L/work/SubtreeLU/lib -lsubtree_lu \
  -Wl,-rpath,/work/SubtreeLU/lib -fopenmp -o /tmp/subtreelu_replay
```

The binary format uses native int32, binary64 and uint64. It requires frozen
completed capture files and nonoverlapping public calls on each instance.
It preserves solve aliasing but does not serialize constructors, destructors,
prior analysis generations, arbitrary parameter changes, cross-call input
addresses or internal scheduling. A successful replay cannot disprove a
state-dependent or scheduling-dependent stall. The original peer timeouts
remain unqualified; diagnostic results never replace competitive observations.

The fixture covers two coefficient generations, separate/in-place solves,
singular refactorization, factorization recovery and original API rejections.
Original-library API teams1/8 × OpenMP teams1/8 × null/provided permutations
passed independent review: hooked/unhooked statuses, output bytes, immediate
floating-point control/exception state and errno match; every captured payload
was checked; all eight strict replays match.26 malformed streams are rejected
before any public call, and trailing entries exercise exit4/5. A supplemental
check covers empty patterns with null/present index pointers. Evidence is in
`build/subtreelu-callstream-fixture-qualification-20261018`, its independent
review, and `build/subtreelu-callstream-empty-pattern-20261018` and review.
The installed recorder only splits a misleadingly indented return statement;
replay changes only its introductory comment relative to those qualified sources.
Full-circuit capture/replay and stall attribution remain separate work.
