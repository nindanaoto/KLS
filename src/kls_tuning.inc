/* Compile-time performance policy for KLS.
 *
 * This is the single catalog for empirical crossover points, admission
 * thresholds, batching sizes, and resource budgets. These values are not
 * correctness requirements: they may be retuned as the implementation,
 * compiler, or target machines change. Changes should be validated with the
 * correctness suite and a representative benchmark corpus, not one matrix.
 *
 * Mathematical constants, sentinel and state values, integer/ABI limits,
 * storage-layout masks, SIMD widths, and ordinary loop/index literals do not
 * belong here. Keeping that boundary prevents this file from becoming a
 * misleading list of every number used by the implementation.
 */
#ifndef KLS_TUNING_INC
#define KLS_TUNING_INC

/* Workload lifecycle and repeated-factorization admission. */
#define KLS_LEAN_LIFECYCLE_REAUDIT_SAMPLES 4
#define KLS_LEAN_SNODE_AVX512_MIN_TARGETS 32
#define KLS_REPEATED_UPDATE_MIN_REFACTORIZATIONS 16
#define KLS_COMPACT_CYCLE_MIN_RELATIVE_SAVING 0.02
#define KLS_COMPACT_CYCLE_MIN_SECONDS_SAVED 1.0e-6
#define KLS_COMPACT_CYCLE_PREPARATION_PAYBACK 2.0

/* Row-oriented refactor batching and dense-panel crossovers. */
#define KLS_ROW_REFACTOR_BATCH_MIN_ROWS 8u
#define KLS_ROW_REFACTOR_BATCH_MAX_ROWS 16u
#define KLS_ROW_REFACTOR_SMALL_SORT_MAX 32u
#define KLS_ROW_REFACTOR_DENSE_MIN_WORK 1024.0
#define KLS_ROW_REFACTOR_COMPACT_PANEL_MIN_WORK 32768.0
#define KLS_ROW_REFACTOR_COMPACT_PANEL_MIN_WORK_PER_ENTRY 8.0
/* Native panels retain their separate floor because their setup cost differs
 * from the generic compact-panel path. Arithmetic intensity is gated below. */
#define KLS_ROW_REFACTOR_NATIVE_PANEL_AUTO_MIN_WORK 32768.0
#define KLS_ROW_REFACTOR_DENSE_BLOCK_ROWS 32u
#define KLS_ROW_REFACTOR_BLOCKED_TRAILING_MIN_ROWS 8u
#define KLS_ROW_REFACTOR_BLOCKED_TRAILING_MIN_COLS 16u
#define KLS_ROW_REFACTOR_BATCH_SUPERNODE_MIN_WORK 32768.0
#define KLS_ROW_REFACTOR_BATCH_SUPERNODE_MIN_WORK_PER_ENTRY 8.0
#define KLS_ROW_REFACTOR_SEPARATOR_BALANCE_BETA 1.2
#define KLS_DENSE_HELP_TAIL_SLICE 128u

/* Row-first pipeline construction and execution. */
#define KLS_BTF_SCALAR_RUN_EXEC_MIN_ROWS 16u
#define KLS_BTF_SCALAR_RUN_EXEC_MAX_ROWS 1024u
#define KLS_ROW_FIRST_PIPELINE_PREFIX_CACHE_REBUILD_MAX_ROWS 32768u
#define KLS_ROW_FIRST_MISSED_PANEL_CACHE_MAX_STORED_ENTRIES 4194304u
#define KLS_ROW_FIRST_PRODUCER_BATCH_MIN_SAVED_STREAM 1024u

/* Triangular solve work sharing. */
#define KLS_ROW_SOLVE_DENSE_TAIL_MIN_NNZ 300000u
#define KLS_ROW_SOLVE_DENSE_TAIL_MIN_FRACTION 0.70
#define KLS_ROW_SOLVE_TRAPEZOID_SLICES 8u
#define KLS_ROW_SOLVE_PARALLEL_RECT_MIN_NNZ KLS_ROW_SOLVE_DENSE_TAIL_MIN_NNZ
#define KLS_ROW_SOLVE_PARALLEL_MIN_ENTRIES_PER_SYNC \
  KLS_ROW_SOLVE_DENSE_TAIL_MIN_NNZ
#define KLS_ROW_SOLVE_CLUSTER_ALPHA_NUMERATOR 2u
#define KLS_PTS_MAX_THREADS 16

/* Supernode batching and EGraph execution. */
#define KLS_SNODE_DEFAULT_MIN_BATCH 3u
#define KLS_SNODE_DEFAULT_MIN_BATCH_WORK 192u
#define KLS_SNODE_MAX_BATCH 64
#define KLS_SNODE_TAIL_CHUNK 32
#define KLS_SNODE_MAX_FUSE 4
#define KLS_EGRAPH_CACHED_SUPERNODE_MIN_ROWS 16u
#define KLS_EGRAPH_CACHED_SUPERNODE_MIN_WORK 512.0
#define KLS_EGRAPH_CACHED_SUPERNODE_MIN_WORK_PER_ENTRY \
  KLS_ROW_REFACTOR_BATCH_SUPERNODE_MIN_WORK_PER_ENTRY
/* Minimum work needed to amortize EGraph worker dispatch. */
#define KLS_EGRAPH_REFACTOR_MIN_FLOPS_PER_THREAD 2.5e5
#define KLS_EGRAPH_REFACTOR_MIN_SIZE 5000u
#define KLS_EGRAPH_POOL_SPIN_ITERS 200000u
/* The caller participates as worker zero. Background workers get only a short
 * grace period between SPICE refactors, avoiding contention in serial phases. */
#define KLS_EGRAPH_WORKER_SPIN_ITERS 4096u
#define KLS_SNB_MAX_ENTRIES 150000000
/* A lone task above this size is divided across the whole worker crew. */
#define KLS_SNB_COOP_MIN_DOUBLES 32768

/* Ordering and nested-dissection policy. */
#define KLS_NICSLU_PARALLEL_R1_THRESHOLD 2.0
#define KLS_NICSLU_PARALLEL_R2_THRESHOLD 50.0
#define KLS_NICSLU_TASK_FLOW_SYNC_COST 1.0
#define KLS_METIS_NDP_MIN_ROWS 30000u
#define KLS_METIS_NDP_MIN_LEAF_ROWS 200u
#define KLS_METIS_NDP_TARGET_DIVISOR 1000u
#define KLS_METIS_NDP_SUBTREE_THRESHOLD_MIN_ROWS 200000u
#define KLS_METIS_NDP_MAX_OVERPARTITION_LEAVES 1024u

/* Fast factor/refactor pipeline admission. */
#define KLS_FAST_FACTOR_PIPELINE_REFACTOR_MIN_WORK 100000000.0
#define KLS_FAST_FACTOR_PIPELINE_REFACTOR_MIN_SHARE 0.95

#endif /* KLS_TUNING_INC */
