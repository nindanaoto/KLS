# Consumer-block experiment removal

Relative to `d5033a3`, this batch removes the opt-in serial consumer-block
supernodal refactor engine from `src/kls_egraph_refactor.inc`:

- `kls_supernodal_refactor_block` and `kls_supernodal_mapped_refactor`;
- `KLS_SN_CONSUMER_MAX_WIDTH` and `KLS_ENABLE_SN_CONSUMER_BLOCKS`;
- the exclusive dispatch and moderate-BTF incumbent bookkeeping;
- the now-unnecessary comparator forward declaration. The comparator's
  definition and production first-factor caller remain.

The paper configurations did not enable this experiment. Normal mapped
refactorization, moderate-BTF comparisons, production supernode/panel kernels,
and dense-tail code remain. No public API, statistics layout, production numeric
threshold, or affinity policy changed. The TSOPF affinity repair and all
unrelated slimming are retained.

The source diff removes **296 lines** with no replacement implementation.

## Validation

The paired harness uses frozen before/after binaries with matching linked-target
build provenance, including optimized SPRAL. Before/after launch order
alternates, one solver runs at a time, and both eight-thread cache domains are
tested. Each launch uses 100 systems with entrywise changes of amplitude 0.001,
verification of every refactor at residual limit `1e-8`, and one-thread auxiliary
library pools. TSOPF_FS_b9_c6 receives 40 pairs per eight-thread domain plus 20
one-thread pairs; eight controls receive four pairs per domain.

Artifacts are in `build/consumer-trim-3sLeWG/`. The previous full campaign and
affinity-repair evidence remain unchanged.

Final Release and ASan/UBSan CTest each pass **5/5**, including the repeated
parallel transpose fixture and caller-affinity checks. Builds have no new
warnings; removed-symbol searches and `git diff --check` are clean.

TSOPF_FS_b9_c6's median paired lifecycle changes are -0.22% at one thread,
+0.32% on the eight-thread 32MiB domain, and -0.50% on the eight-thread 96MiB
domain. No eight-thread launch in either version exceeded 0.6ms average steady
refactor time (the descriptive cutoff used in the affinity diagnosis). Fixed
after medians were 0.410ms and 0.415ms. The affinity repair remains effective
in this sample.

All **328 initial launches passed**. Initial control lifecycle changes
(positive means slower) were:

| Matrix | 32MiB, 8 threads | 96MiB, 8 threads |
| --- | ---: | ---: |
| TSOPF_FS_b9_c1 | -0.39% | -0.03% |
| TSOPF_RS_b9_c6 | -1.08% | -0.17% |
| onetone1 | +0.98% | -1.38% |
| twotone | +0.84% | -0.01% |
| rajat25 | +0.41% | -1.03% |
| adder_dcop_01 | +0.92% | +4.39% |
| 1138_bus | +0.67% | -7.56% |
| bips98_1142 | -1.96% | +1.83% |

The two largest absolute control differences received independent 20-pair
rechecks on 96MiB. All 80 additional launches passed; paired changes were
+0.65% for adder_dcop_01 and -2.12% for 1138_bus. Their steady-refactor medians
were essentially unchanged (16.725us to 16.709us and 9.204us to 9.182us).
The initial large differences did not reproduce at the same magnitude.
Recheck evidence is in `build/consumer-recheck-56jrY2/`.

In total **408/408 launches passed**. No large repeatable regression was
established in these checks; they do not establish full-corpus parity or
guarantee exact performance equality. The final Release binary is identical
to the frozen after binary used for both measurement batches.
