# Repair timer and scalar trace cleanup after alignment recovery

This resumed cleanup pass starts at pushed `372f9d1`. Retain the general
SNB alignment fix and all preceding slimming. The earlier pass stopped at
the cumulative onetone1 regression; its stopping report is historical, not
a claim that no candidates remained. Historical reference: `88a177e`.

Remove optional ETree repair phase timing (`KLS_TRACE_FAST_REJECT_REPAIR_TIMING`)
and scalar internal/output touch accounting from the first-factor pipeline.
Also remove `scalar_dep_l_entries`: it was incremented exactly alongside
`scalar_dep_calls`, with identical saturation and aggregation. The remaining
dependency count reports that information without duplicate state.

There are no changes to arithmetic, dependency discovery, repair acceptance,
allocation policy, numerical safeguards, scheduling decisions, or public
statistics. Private trace structures shrink. General pipeline progress and
producer trace accounting remain. The summary helper retains archived trace
support and explicitly documents the removed scalar fields. Net C-source
reduction: 83 lines.

Release and ASan/UBSan CTest each pass 6/6. A trace-enabled release smoke
run passes and emits 24 general pipeline trace records without the removed
scalar fields or ETree timer output, even with the retired timing control
set. The summary helper parses the current log successfully.

The SNB kernel stays at 0x65000 with size 0x1931. Its disassembly differs
from `372f9d1` only in the relocation of its call to the prelude (which
moves from 0x64bf0 to 0x64be0); the numerical body and internal instruction
addresses remain unchanged. This verifies the hot body's alignment recovery
is retained, rather than assuming diagnostic removal cannot affect emitted
code. The prelude's new placement still requires performance validation.

Focused H100 artifacts are under `build/prep-trim-repair-dBcoha/scalar-trace-*`:
preceding aligned baseline (`original`), historical `88a177e` (`start`), and
candidate (`fix`) rotate in pinned, sequential runs with clean environments,
verified build provenance, exact commands, hashes, source diff, and numerical
validation of every entrywise-changed refactor at a 1e-8 residual limit.
No builds/tests overlap timing. This is not a full paper or Xyce rerun.

The initial campaign passes 375/375 launches: 36 onetone1 comparisons
(six triplets per cache domain), 303 automatic-policy controls across six
other matrices and nine configurations, and 36 forced-KLS-first controls.
Onetone1 changes by +0.14% and +0.23% versus `372f9d1`, and -0.15% and
-0.38% versus `88a177e`. The alignment recovery remains intact on both domains.

The eight-triplet serial TSOPF_FS_b9_c6 screen flags +2.65% lifecycle time,
concentrated in steady refactoring (marginal medians 1.032 vs 0.990 ms).
It is repeated before acceptance: 32 rotating triplets (96/96 valid launches)
measure +0.30% lifecycle versus `372f9d1` and +0.24% versus `88a177e`.
The small-sample slowdown does not reproduce. No further source change
was made during that investigation. Evidence: `scalar-trace-tspof-recheck*`.

The chunk is retained after 471/471 numerically valid launches and no
repeatable approximately 3% lifecycle regression. This is not proof of zero
slowdown on every matrix. The candidate source diff matches the recorded
campaign diff. Forced-first controls exercise KLS-first but do not establish
timed parallel-pipeline coverage; the trace-enabled smoke supplies correctness
coverage of the pipeline.

The broader low-risk goal remains active. Next candidate: producer
internal/output touch counters and their output parameters. Their callers
use the totals only for trace aggregation; actual producer update work,
dependency handling, and the NOT_ROOT catch-up decision must remain intact.
That candidate has not been removed or performance-validated here.
