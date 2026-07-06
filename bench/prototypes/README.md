# Refactor-kernel prototypes (task #12 and its refuted predecessors)

Standalone benchmarks against the vendored KLU (`libkls_suitesparse.a`).
Build (from repo root, after building `build/`):

    gcc -O3 -march=native -o sb bench/prototypes/supernodal_bench.c \
      -Ithird_party/suitesparse/KLU/Include -Ithird_party/suitesparse/BTF/Include \
      -Ithird_party/suitesparse/AMD/Include -Ithird_party/suitesparse/COLAMD/Include \
      -Ithird_party/suitesparse/UFconfig build/libkls_suitesparse.a -lm

    ./sb <matrix.mtx> [reps] [zeta] [wnarrow]

- `supernodal_bench.c` — task #12 round 3 (in progress): amalgamated
  supernodes, COLUMN-major dense L panels as primary storage, per-edge
  TRSM + GEMM-into-temp + column scatter, hybrid scalar path (w<wnarrow).
  Measured t1 vs the landed batch kernel: mc2depi 1.25x WIN, ASIC_320ks
  0.77x, rajat25 0.44x (edge-waste diagnosis pending).  Verify mc2depi /
  rajat25 by residual, not Udiag (chaotic for all kernels there).
- `supernodal_bench_r1.c.bak` — round 1 (row-major panels), for reference.
- `snode_bench.c` — the landed per-column chunked batch kernel's prototype
  (`mine-panel` mode = production shape) + refuted dense-mirror mode 2.
- `consumer_gemm.c` — refuted joint consumer-run designs (scattered-W and
  compact-row variants).
- `consumer_share.c` — measures the GEMM-able flop share per consumer run.
- `match_tool.c` — standalone SPRAL Hungarian matcher (the block-ordering
  prototype pipeline; superseded by the in-tree raw-value matcher).

Full history and measurements: Plan.md and the session memory
(dense-endgame-design).
