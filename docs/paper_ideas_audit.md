# Paper Ideas Audit

This note records which ideas from the reference papers are present in KLS and
which remain open. The intent is to keep KLS development focused on general
solver algorithms instead of tuning individual benchmark matrices.

## Current Conclusion

The first-factor separator-private pivot path is now closer to SubtreeLU's
private-mode rule: private pivot search is restricted to columns owned by the
same private worker. This fixes the earlier `pre2` diagnostic where a private
row selected a pipeline-owned pivot candidate and then failed only at the
post-selection owner check. The current trace
`build/kls_pre2_pivot_reject_no_partial_trace35.stderr` instead fails the
private phase at `failed_row=130192` with `reason=pivot-reject`, after
`38134` rows in that worker and `480390` total completed private rows across
workers. That is the expected paper-level boundary: private mode cannot safely
continue when the row needs a cross-owner dynamic pivot.

Partial publication of those completed private rows was tested and rejected for
pivot-related failures. With all completed private rows published, the active
tail shrank to `149238` rows, but relaxing the tail to allow cross-component
active pivots later failed at `reason=pivot-exchange`: the selected active
column conflicted with already-published private U rows. With a first-failed-row
cutoff, the active tail grew to `499436` rows and still did not provide a
CKTSO-scale path. KLS therefore only allows partial separator-private
publication for non-pivot private failures; `pivot-reject`, `pivot-owner`,
`pivot-exchange`, `pivot-verify`, and append failures discard private work and
fall back to the ordinary full row pipeline. On `pre2`, the guarded fallback
records `private_rows=0`, `partial_active=0`, and the full row pipeline reaches
`589824/629628` rows in the 35s trace. Closing this gap still requires the
larger paper-level mechanism: a pivot-compatible grouped-current/deferred-swap
numeric owner, not just committing private prefixes after a dynamic-pivot
failure.

A follow-up repair of the initial partitioned separator plan was prototyped and
rejected before commit. The idea promoted private rows with invalid left
dependencies into the pipeline, then rebuilt the private queues instead of
dropping the SubtreeLU Algorithm 6-style partitioned plan. On `pre2`, the
repaired plan validated and kept the partitioned shape, but it over-expanded
the private phase: `build/kls_pre2_partition_repair_trace45.stderr` reported
`final_partitioned=1`, `final_private_rows=623574`, and only `6054` pipeline
rows. The private phase then failed late at `failed_row=471797` with
`reason=pivot-reject` after `188557` completed rows in the failing worker,
`29266` dynamic pivots, and `16666624` worker-local U entries, before falling
back to the full row pipeline. The resulting full-pipeline trace reached only
`131072/629628` rows under the 45s cap, far worse than the guarded legacy split
trace. This rejects ownership-repair promotion as a default; preserving the
partitioned plan needs a pivot-compatible private/pipeline owner, not simply
more private rows.

The separator-tree row refactor now has a validation-gated ordered-private
executor for SubtreeLU-style private subdomains. The separator FLOP queue
already promotes any private group reached from pipeline work back into the
pipeline side. KLS now additionally validates that every private-to-private
successor stays in the same worker's ordered private list and appears later in
that list. When the proof passes, private groups execute without spinning on
`row_pipeline_remaining_preds` and without atomically releasing private
successors; they still release pipeline successors through the guarded ready
queue. This makes the retained separator-private phase closer to the paper
algorithm instead of running thread-local subdomains through the generic
pipeline dependency scaffold. The new diagnostics
`row_refactor_last_separator_flop_ordered_private` and
`row_refactor_separator_flop_ordered_private_run_count` show whether this
path was selected.

Correctness passed `git diff --check`, `cmake --build build -j2`,
`./build/kls_smoke`, and `ctest --test-dir build --output-on-failure`. On the
forced row-refactor top-five CKTSO-gap slice, the ordered-private run
`build/kls_ordered_private_forced_row_gap5_t4_r1_ref3_timeout120.jsonl`
measured `3.8786s` SPICE-cycle geomean versus `4.2265s` in the stored forced
row current artifact
`build/kls_forced_row_refactor_current_gap5_t4_r1_ref3_timeout120.jsonl`
(`1.09x` geomean speedup, wins on all five rows). The ordered-private stat was
active for all three refactor repeats on the large ASIC rows
(`ASIC_320ks`, `ASIC_320k`, `ASIC_100ks`) and inactive on the tiny private
queues in `gemat12` and `rajat03`. This is retained as a solid paper-aligned
row-refactor cleanup, but it does not close the main gap: the same focused
default artifact remains `1.4664s`, so forced row refactor is still about
`2.65x` slower than the normal KLS path on this slice.

The row-refactor group builder now reports scalar batching blocker diagnostics:
`row_refactor_group_scalar_candidate_count`,
`row_refactor_group_scalar_candidate_rows`,
`row_refactor_group_scalar_short_count`,
`row_refactor_group_scalar_short_rows`, and the first stop reason split across
level mismatch, internal dependency, next segment, max width, and matrix end.
The focused forced-row probe
`build/kls_row_group_blockers_gap5_t4_r1_ref3_timeout120.jsonl` measured
`4.0429s` SPICE-cycle geomean and showed that internal row dependencies are not
the scalar batching blocker. On `ASIC_320ks`, `204,113` of `219,053` groups
were single-row groups; scalar candidates covered `509,279` rows, `405,067`
candidate rows were discarded by the minimum-width rule, and `213,854` of
`214,348` scalar candidates stopped first at the same-level boundary. `ASIC_320k`
was similar: `230,012` of `242,821` groups were singles and `237,505` of
`238,002` scalar candidates stopped first at level mismatch.

A level-relaxed scalar batching prototype was tested and rejected before
commit. The prototype allowed scalar batches to cross row levels while still
requiring no internal row dependencies, directly testing whether the SubtreeLU
private-refactor gap was mostly scheduler fragmentation. It reduced group counts
on the ASIC rows (`ASIC_320k` `242,821 -> 218,620`, singles
`230,012 -> 203,963`; `ASIC_100ks` `54,297 -> 40,871`, singles
`43,358 -> 28,692`), but the top-five forced-row geomean still worsened
slightly (`4.0429s -> 4.0664s`) in
`build/kls_relaxed_scalar_batch_gap5_t4_r1_ref3_timeout120.jsonl`. This rejects
"coarsen scalar private groups across row levels" as the clear slow-case gap
closer. The missing paper-aligned piece is therefore not just fewer ready-queue
groups; it remains a production row-major supernode/current-state numeric
executor that fuses or avoids the scalar row/current-state work itself.

The compact producer-consumer row update then gained a portable packed scaled
accumulator for its retained trailing panel update. A callgrind run with
collection toggled around `kls_refactor` on one-thread forced-row `ASIC_100ks`
showed `kls_row_refactor_try_compact_supernode_update` at about `34%` of
isolated refactor instructions, with the inner
`trailing_workspace += lij * dep_panel` loop accounting for about `16%`. KLS now
routes that contiguous row-panel operation through one always-inline helper
used by both the normal compact update and the prefactor compact-supernode
update. This is not a BLAS threshold change; it keeps the LGPL KLS-owned
SubtreeLU-style trailing accumulation while making the packed row update a
first-class numeric kernel primitive.

Against the previous committed control
`build/kls_control_c76_pass3_gap5_t4_r1_ref3_timeout120.jsonl`, the retained
helper run
`build/kls_accum_scaled_pass3_gap5_t4_r1_ref3_timeout120.jsonl` improved the
forced row-refactor top-five CKTSO-gap pass-3 median geomean from `4.0007s` to
`3.7969s` (`0.949x`). The three matrices that execute compact supernode
updates improved together by `0.896x` geomean: `ASIC_320ks` `26.77s -> 23.94s`,
`ASIC_320k` `31.46s -> 27.16s`, and `ASIC_100ks` `15.54s -> 14.48s`. `gemat12`
and `rajat03` do not execute compact-supernode updates in this run and moved
slightly backward within the usual small-matrix noise. This is retained as a
real row numeric-kernel improvement, but it does not change the broader
conclusion: KLS still needs the larger production row/segment panel executor to
make forced row refactor competitive with the default EGraph path and CKTSO.

The same refactor-only callgrind pass exposed one retained scheduler-build
cleanup that is worth keeping but not enough to change direction. Row-refactor
successor ordering sorted every successor list through a malloc/qsort/free path,
even though the row task graph is dominated by tiny fanout lists. KLS now sorts
lists of at most `32` groups in place with the same work-descending and
group-id tie order, while preserving qsort for large root lists. On one-thread
forced-row `ASIC_100ks`, the callgrind instruction total moved from `3.8385B`
after the packed accumulator to `3.7862B`, and the visible libc qsort-family
entries fell from tens of millions of instructions to only the remaining
large-list sort. The same-machine top-five forced-row pass-3 comparison against
the previous committed control `0ed2402` measured
`build/kls_control_0ed2402_forcedrow_gap5_t4_r1_ref3_p3_timeout120.jsonl` at
`3.7495s` geomean and
`build/kls_small_sort_forcedrow_gap5_t4_r1_ref3_p3_timeout120.jsonl` at
`3.7331s` (`0.996x`). This is retained as a scheduler-overhead reduction, not
as the missing CKTSO/SubtreeLU row numeric executor.

A direct current-panel scalar scatter prototype was then tested and rejected
before commit. The trial targeted the remaining scalar fallback inside compact
dense row groups: after direct input initialized the current dense/trailing row
panel, scalar external dependencies wrote current-group and trailing targets
straight into that panel while leaving truly external targets in `x`. This
looked paper-aligned because it tried to avoid the row-panel round trip through
the sparse workspace, but the target checks and mixed scatter pattern lost badly.
The focused top-five forced-row run
`build/kls_current_panel_scalar_forcedrow_gap5_t4_r1_ref3_p3_timeout120.jsonl`
measured `4.4153s` geomean versus `3.7331s` for the retained small-sort build;
`ASIC_100ks` regressed from about `14.01s` to `25.77s`. This rejects a
per-scalar current-panel scatter shortcut. The remaining useful path is still a
larger grouped producer/current executor that batches the target mapping once
per panel or selected run, not once per scalar dependency.

An active producer-bucket grouped advance executor was prototyped and rejected
before commit. The trial kept the retained position-coded Algorithm 5
pre-prefix state plan, but replaced the runtime scan over active current
columns with per-producer buckets that requeued a current only after its next
U dependency had been applied. This preserves each current column's triangular
order and directly tests whether the missing paper-level owner is mainly the
active-current scan around grouped producer work. It is not: on the same-source
five-matrix CKTSO-gap focus manifest
`build/kls_group_advance_bucket_default_gap5_t4_r1_ref3_timeout120.jsonl`
measured `7.8971s` SPICE-cycle geomean, while the active-bucket prototype
`build/kls_group_advance_bucket_gap5_t4_r1_ref3_timeout120.jsonl` measured
`22.0281s`. Refactor time roughly doubled on the hard ASIC rows
(`ASIC_320k` `0.0989s -> 0.1679s`, `ASIC_320ks` `0.0802s -> 0.1643s`,
`ASIC_100ks` `0.0444s -> 0.0923s`) and `G2_circuit` hit the 120s ceiling.
The active path did run (`ASIC_320k` seeded `3,916` current states and
`3,718,316` sparse rows), so the loss is not lack of coverage. This rejects
"bucket the already retained sparse current states" as the clear CKTSO-gap
closer; the remaining refactor work still has to avoid or fuse the huge
unbounded advance/current-state stage itself.

The latest refactor probe adds default-safe bounded owner diagnostics for the
exact-shape group-L/refactor surface. The new
`refactor_supernode_consumer_plan_shape_bounded_advance_*` stats count exact
shape groups where every run reaches the retained producer panel after at most
`128` scalar advance dependencies, then separately count the subset whose cached
panel update entries pay for the advance work. Correctness passed
`git diff --check`, `cmake --build build -j2`, `./build/kls_smoke`,
`KLS_REFACTOR_SUPERNODE_CONSUMER_PLAN_OUTPUT_STATS=1 ./build/kls_smoke`, and
`ctest --test-dir build --output-on-failure`.

The result rejects a small bounded-prefix refactor owner as the clear CKTSO-gap
closer. With group-L cache stats enabled, the focused top-five gap run
`build/kls_bounded_refactor_owner_groupcache_gap5_t4_r1_ref3_timeout120.jsonl`
completed all rows at `1.5700s` SPICE-cycle geomean. The exact-shape surface is
large on the slow ASIC rows, but the bounded-and-payoff subset is tiny:
`ASIC_320ks` has `746,229` exact-shape run rows and `2.145e9` advance work, but
only `23,262` bounded rows and `348` bounded-payoff rows; `ASIC_320k` has
`755,745` exact-shape rows and `2.475e9` advance work, but only `21,697`
bounded rows and `96` payoff rows; `ASIC_100ks` has `1,047,047` exact-shape rows
and `5.691e9` advance work, but only `35,262` bounded rows and `4,416` payoff
rows. `rajat03` has more bounded coverage (`27,202` of `35,706` rows), but only
`345` payoff rows.

This sharpens the refactor target. KLS does not just need a guard around small
advance prefixes; the missing paper-level mechanism must own or avoid the huge
unbounded advance stage. A paper-shaped producer-panel/multi-current executor
still needs to batch producer work across currents, fuse advancement with prefix
application, and avoid publishing per-current prepared states that the scalar
consumer later replays.

The bounded-owner diagnostic limit is now runtime-configurable through
`KLS_REFACTOR_PLAN_GROUP_L_BOUNDED_ADVANCE_MAX_DEPS`, with the existing `128`
dependency bound kept as the default and reported in
`refactor_supernode_consumer_plan_shape_bounded_advance_dep_limit`. This is a
diagnostic knob only: it changes retained-plan accounting, not numeric
execution. Validation passed `git diff --check`, `cmake --build build -j2`,
`./build/kls_smoke`, `ctest --test-dir build --output-on-failure`, and a
malformed-value smoke check with
`KLS_REFACTOR_PLAN_GROUP_L_BOUNDED_ADVANCE_MAX_DEPS=not-a-number`.

Rerunning the same top-five exact-shape/group-cache probe at wider bounds
confirms that merely increasing the bounded prefix window is not enough. The
new artifacts
`build/kls_bounded_refactor_owner_512_groupcache_gap5_t4_r1_ref3_timeout120.jsonl`
and
`build/kls_bounded_refactor_owner_2048_groupcache_gap5_t4_r1_ref3_timeout120.jsonl`
completed cleanly. At the default `128` bound, the focused set had `107,423`
bounded exact-shape run rows and only `5,205` payoff rows. At `512`, bounded
rows grew to `259,289`, but payoff rows only reached `7,605`. At `2048`, nearly
the whole exact-shape surface became bounded (`2,438,635` of `2,584,727` run
rows), yet payoff rows stayed at `7,605` and payoff update entries stayed at
`376,831`. On the hard ASIC rows this is especially clear: `ASIC_320k` went
from `21,697` bounded rows at `128` to `754,674` at `2048`, but remained at
only `96` payoff rows. The next refactor implementation should therefore not
be a larger bounded-prefix owner; it needs the broader grouped live-workspace or
row-major supernode owner that reduces the advance stage itself.

The previous compact retained-state probe implemented an opt-in compact plan for
Algorithm 5 grouped pre-prefix advance:
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_ADVANCE_COMPACT_STATE=1`.
For the compatible direct-prefix handoff path, the current-state row plan now
retains only the advance closure plus prefix output rows instead of cloning the
full scatter/current/suffix/target row set. It falls back to the full retained
state when final-state or suffix-completion modes are requested. Correctness
passed `git diff --check`, `cmake --build build -j2`, `./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_ADVANCE_COMPACT_STATE=1 ./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_ADVANCE_COMPACT_STATE=1 KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_ADVANCE_POS=1 ./build/kls_smoke`,
and `ctest --test-dir build --output-on-failure`.

The focused top-five CKTSO-gap result rejects compact state as a default. The
current default control measured `1.3839s` SPICE-cycle geomean in
`build/kls_compact_state_default_gap5_t4_r1_ref3_timeout120.jsonl`.
Compact state alone measured `3.3963s` in
`build/kls_compact_state_gap5_t4_r1_ref3_timeout120.jsonl`; compact state plus
the position map measured `2.6173s` in
`build/kls_compact_state_pos_gap5_t4_r1_ref3_timeout120.jsonl`, essentially tied
with the full-state position control at `2.6128s`
(`build/kls_pos_control_gap5_t4_r1_ref3_timeout120.jsonl`). The row-count
reduction is real: on `ASIC_320ks`, retained current-state rows dropped from
`2,568,461` in the full-state position control to `910,070`, and seeded rows
dropped from `3,564,401` to `1,326,362`; `ASIC_320k` dropped from `2,925,515`
to `1,067,018` retained rows; `ASIC_100ks` dropped from `783,675` to `195,090`.
Refactor time did not move enough (`ASIC_320ks` `0.1827s` full-state position
control vs `0.1809s` compact+position), so retained sparse-state size is not the
first-order CKTSO gap.

A current-source grouped-final-state probe rejects a tempting shortcut for that
executor. The trial changed the guarded direct-prefix final-state path in
`kls_egraph_refactor_apply_algorithm5_payoff_group_prefix_items` from
current-major replay to a producer-row/current grouped loop. The first version
used sorted retained-state lookups inside the dense update loop; the focused
top-five run measured `4.9728s` geomean in
`build/kls_finalstate_grouped_groupcomplete_gap5_t4_r1_ref3_timeout120.jsonl`
versus `1.4722s` for the same-source default control
`build/kls_finalstate_grouped_default_gap5_t4_r1_ref3_timeout120.jsonl`.
Precomputing prefix and dense-target state positions once per batch did not fix
the issue: `build/kls_finalstate_grouped_prepos_groupcomplete_gap5_t4_r1_ref3_timeout120.jsonl`
still measured `4.9402s`. The source change was removed. This closes the
"just flip the final-state loop order" shortcut; the required paper-level owner
must carry efficient per-current state-position metadata with the grouped
numeric task, not rebuild or binary-search state maps during completion.

A batched-hash version of the same final-state shortcut was also tested and
rejected before commit. The trial added an opt-in
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_FINAL_STATE_HASH=1`
path that built one retained-state hash table for all currents in a grouped
completion batch and then applied the final-state prefix in producer-row/current
order. Correctness passed `cmake --build build -j2`, `./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_COMPLETE=1 ./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_COMPLETE=1 KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_FINAL_STATE_HASH=1 ./build/kls_smoke`,
and `ctest --test-dir build --output-on-failure`. The focused top-five
CKTSO-gap run measured `3.1472s` geomean in
`build/kls_finalstate_hash_groupcomplete_hash_gap5_t4_r1_ref3_timeout120.jsonl`
versus `2.5542s` for the same-source group-complete control
`build/kls_finalstate_hash_groupcomplete_gap5_t4_r1_ref3_timeout120.jsonl` and
`1.4463s` for the default control
`build/kls_finalstate_hash_default_gap5_t4_r1_ref3_timeout120.jsonl`. This
confirms that a lookup-only owner does not fill the paper gap; the missing
piece is still reducing the grouped prefix/suffix numeric work itself.

A direct scalar BTF dispatch split was tested and rejected before commit. The
prototype added a plain scalar dependency helper for the default BTF refactor
loop, passed `git diff --check`, `cmake --build build -j2`, `./build/kls_smoke`,
and `ctest --test-dir build --output-on-failure`, then compared against a fresh
same-source top-five CKTSO-gap baseline. The baseline
`build/kls_scalar_plain_baseline_gap5_t4_r1_ref3_timeout120.jsonl` measured
`1.4524s` SPICE-cycle geomean. The split helper measured `1.4547s` in
`build/kls_scalar_plain_candidate_gap5_t4_r1_ref3_timeout120.jsonl`, and a
three-repeat check measured `1.4609s` in
`build/kls_scalar_plain_candidate_gap5_t4_r3_ref3_timeout120.jsonl`. It slightly
helped `ASIC_320ks` but regressed `ASIC_100ks`, so the source change was
reverted. This is useful negative evidence: shaving format dispatch inside
`kls_egraph_refactor_apply_btf_scalar_dep` is not the clear missing paper
mechanism. The large gap remains reducing the number of scalar dependency
applications through a grouped producer-panel/multi-current executor.

A compact final-state variant was tested and rejected before commit. Extending
the compact row plan to direct-prefix final-state completion has to retain the
pivot row, current-column L rows, prefix target rows, and every remaining suffix
dependency closure. On the hard ASIC rows that closure was effectively the same
size as the full retained state (`ASIC_320ks` and `ASIC_320k` still reported
`2,568,461` and `2,925,515` retained current-state rows). Worse, enabling the
compact knob together with group completion also enabled grouped advance prep
and increased seed work: the top-five run
`build/kls_group_complete_compact_state_gap5_t4_r1_ref3_timeout120.jsonl`
measured `3.6382s` geomean versus `2.5457s` for
`build/kls_group_complete_control_gap5_t4_r1_ref3_timeout120.jsonl`. This
confirms that the final-state row set cannot be shrunk enough to matter; the
remaining useful direction is a real grouped suffix/current-state executor, not
another retained-row pruning mode.

A SubtreeLU-inspired suffix `2P` fanout guard was also tested and rejected
before commit. The trial let
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_SUFFIX_GROUP_ADVANCE=1` build
a producer batch only when the retained fanout for that trigger was at least
twice the worker count, matching the paper's "small unfinished supernodes should
wait" rule. Correctness passed default smoke and suffix-group smoke, but the
focused top-five run
`build/kls_suffix_group_2p_gap5_t4_r1_ref3_timeout120.jsonl` measured `2.8081s`
geomean versus `1.4498s` for the same-binary default control
`build/kls_suffix_2p_default_gap5_t4_r1_ref3_timeout120.jsonl`, and it was
worse than the earlier suffix position-cache run. The rejected guard confirms
that the current suffix-group path is dominated by retained-state ownership and
cursor setup, not by too many tiny fanout batches.

The previous refactor probe turns the grouped pre-prefix advance result
into a position-coded executor. With
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_ADVANCE_POS=1`, KLS now
retains `advance occurrence -> current-state row position` spans for each
Algorithm 5 pre-prefix dependency and streams the grouped producer L column
without binary-searching or hashing sparse state rows in the inner update loop.
Correctness passed `cmake --build build -j2`, `git diff --check`,
`./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_ADVANCE_POS=1 ./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_ADVANCE_POS=1 KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_FINAL_STATE=1 ./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_ADVANCE_POS=1 KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_COMPLETE=1 ./build/kls_smoke`,
and `ctest --test-dir build --output-on-failure`.

The focused same-binary result is still rejected as a default CKTSO-gap closer.
The position-coded path measured `2.5584s` top-five SPICE-cycle geomean in
`build/kls_advance_pos_gap5_t4_r1_ref3_timeout120.jsonl`, better than the
same-source hash-only grouped advance control at `2.8508s`
(`build/kls_advance_hash_control_gap5_t4_r1_ref3_timeout120.jsonl`) but still
well behind the default control at `1.5000s`
(`build/kls_advance_pos_default_gap5_t4_r1_ref3_timeout120.jsonl`). Adding the
hash table for the remaining final row copies did not help (`2.5906s` in
`build/kls_advance_pos_hash_gap5_t4_r1_ref3_timeout120.jsonl`). On the hard
ASIC rows the opt-in path definitely ran (`ASIC_320ks` seeded `4,010` current
states / `3,564,401` sparse rows, `ASIC_320k` `3,916` / `3,718,316`,
`ASIC_100ks` `638` / `586,074`), but refactor time still roughly doubled
against default (`0.1869s` vs `0.0832s` on `ASIC_320ks`, `0.1983s` vs
`0.1017s` on `ASIC_320k`, `0.0858s` vs `0.0440s` on `ASIC_100ks`). This narrows
the missing paper-level piece: row-position lookup was a real overhead, but not
the first-order CKTSO gap. That result motivated the compact retained-state test
above.

The refactor focus now has a clearer missing paper-level surface before the
suffix probes: Algorithm 5 payoff groups have substantial duplicate producer
work in the pre-prefix advance section. KLS records this with the new
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_ADVANCE_MAP=1` descriptor,
which builds `group -> advance producer column -> (run, U position)` and reports
unique/duplicate/shared fanout counters. The focused artifact
`build/kls_alg5_advance_map_gap5_t4_r1_ref3_timeout120.jsonl` completed all five
rows with `1.7971s` geomean, which is not a speed claim because the map is still
descriptor-only. The important result is the hard-ASIC sharing: `ASIC_320ks` has
`190,994` grouped advance deps but only `53,758` unique producer keys
(`137,236` duplicate occurrences, max fanout `151`, `26,825,778` duplicate
update-entry proxy work), `ASIC_320k` has `194,219` / `61,682` / `132,537` /
`165` / `25,161,014`, and `ASIC_100ks` has `32,849` / `14,443` / `18,406` /
`53` / `2,166,521`. This makes the next CKTSO/SubtreeLU-aligned fix a true
grouped advance owner feeding the prefix executor, not BLAS thresholding,
suffix-only replay, or per-current row-state micro tuning.

The guarded suffix-advance follow-up filled the next direct Algorithm 5
surface, but it is rejected as a working CKTSO-gap closer. The implementation
adds a retained `current_up` cursor, a suffix-trigger map from completed
dependencies to prepared current slots, and a live current-state suffix executor
that advances READY current states as later dependencies finish. It also extends
the current-state row plan to include suffix dependency rows and their L-target
rows, which is required for the paper-style live state update. Because repeated
four-thread runs showed residual drift on `ASIC_320k`, the executor is now gated
behind the explicit unsafe opt-in pair
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_SUFFIX_ADVANCE=1` and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_SUFFIX_ADVANCE_UNSAFE=1`.
The normal suffix flag alone no longer enables the multicore live-state
executor.

The rejected benchmark artifact
`build/kls_alg5_suffix_advance_reset_gap5_t4_r1_ref3_timeout120.jsonl`
completed all five focused rows with no process failures but measured
`3.1956s` SPICE-cycle geomean, worse than the final-trigger-only surface at
`2.5581s` and the default runs around `1.4s`. The three active ASIC rows
allocated `6,277,651` retained current-state rows and advanced `1,108,677`
suffix dependencies, but that replay performed `409,264,540` suffix update
entries. A rerun,
`build/kls_alg5_suffix_advance_reset_gap5_t4_r1_ref3_timeout120_rerun.jsonl`,
was worse at `3.3224s`. Isolated repeated `ASIC_320k` probes under the unsafe
suffix executor produced residuals around `1e-8` to `4e-8`, while the default
and final-trigger controls stayed at `2.04e-15`. This confirms the next missing
paper-level piece is not BLAS thresholding and not per-current suffix replay; it
is a correct grouped owner for live multi-current state, or a different
supernodal task decomposition that avoids mutating partially owned current
states from multiple suffix triggers.

The producer-triggered suffix follow-up extends the retained suffix map into a
reverse `producer column -> unique suffix-map entry` trigger table and reuses
the scalar in-order current-state suffix helper when a producer column
finishes. This is closer to the paper's producer-side scheduling surface than
the plain suffix trigger, but it is still not the grouped owner described by
the Algorithm 5/SubtreeLU-style algorithms: every current slot is claimed and
advanced independently through the scalar current-state path. It is therefore
quarantined behind the explicit unsafe opt-in pair
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_SUFFIX_PRODUCER_ADVANCE=1` and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_SUFFIX_PRODUCER_ADVANCE_UNSAFE=1`.
The single producer flag alone is intentionally inert.

The focused producer-trigger artifact
`build/kls_alg5_suffix_producer_advance_unsafe_gap5_t4_r1_ref3_timeout120.jsonl`
regressed to `3.1927s` SPICE-cycle geomean and produced relative-residual drift
on `ASIC_320ks` (`3.08376272E-8`). The active rows advanced `1,109,202` runtime
suffix dependencies and an estimated `409,499,854` suffix update entries, which
matches the expected large work surface but does not close the CKTSO gap. This
rejects "drive more scalar current-state suffix advancement from producer
triggers" as the next lead cause. The missing paper-level piece remains a true
grouped live-workspace executor with clear completion ownership, not BLAS
thresholding or small local replay changes.

After the current-state row plan was fixed to include suffix targets for
producer-triggered advancement, the producer path was rerun. The default
completion-state branch still drifted on `ASIC_320ks`
(`build/kls_alg5_suffix_producer_advance_unsafe_rowplanfix_gap5_t4_r1_ref3_timeout120.jsonl`,
`2.18e-8` relative residual, `3.1941s` geomean). Forcing the live-state
completion branch made the five-row run numerically clean while still advancing
about `1.109M` suffix dependencies and `409.6M` update entries
(`build/kls_alg5_suffix_producer_live_state_gap5_t4_r1_ref3_timeout120.jsonl`),
but it measured `3.1440s` geomean. That is still slower than the same-build
default at `1.4427s` and much slower than the CKTSO medium artifact near
`0.5505s`. So the paper gap is not merely "allow more suffix advancement"; KLS
needs to amortize the duplicate producer work through a grouped numeric owner.

The grouped suffix-trigger follow-up adds an opt-in
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_SUFFIX_GROUP_ADVANCE=1`
executor that uses the retained suffix-map trigger table directly instead of
replaying the scalar producer helper. To preserve U-order and avoid the residual
drift seen in the unsafe scalar suffix probes, it only applies a dependency when
the retained current-state cursor is exactly at that producer (`current_up ==
target_up`). This is numerically clean on the focused gap run, but it is not a
CKTSO-gap closer: `build/kls_alg5_suffix_group_advance_gap5_t4_r1_ref3_timeout120.jsonl`
measured `2.6927s` geomean versus the same-build default
`build/kls_default_compare_gap5_t4_r1_ref3_timeout120.jsonl` at `1.4427s` and
the CKTSO medium artifact at about `0.5505s`. The strict in-order grouped path
advanced only `1,249` suffix dependencies and `492,359` update entries, all on
`ASIC_320k`, so it avoids the correctness failure by giving up most of the
paper's intended live overlap. This rejects BLAS thresholding and
single-producer in-order suffix triggering as the first cause; the remaining
paper-level gap is a supernode/window owner that can safely process a broad
finished-prefix frontier, closer to SubtreeLU's split wait/update/wait schedule,
instead of advancing one dependency per current-state cursor.

The prepared-state suffix catch-up extends that safe grouped path without
re-enabling the unsafe live-state replay: when prefix preparation still owns a
current slot in `PREPARING`, it now advances already-finished suffix dependencies
inside that private current-state window and immediately hands fully advanced
currents to the final-state completer. The focused opt-in probe
`build/kls_alg5_prepared_suffix_group_gap5_t4_r1_ref3_timeout120_rerun.jsonl`
completed cleanly and advanced far more suffix work (`98,106` / `134,593` /
`3,857` deps on the three ASIC rows, with `24` finished currents) than the
strict trigger-only path. It still regressed to `2.7651s` geomean, while the
same-source default
`build/kls_prepared_suffix_default_gap5_t4_r1_ref3_timeout120.jsonl` stayed at
`1.4287s`. This confirms the missing CKTSO/SubtreeLU step is not just "catch up
more suffix dependencies"; the current per-current row-state lookup and retained
state ownership model must be replaced by a grouped multi-current producer/window
executor that amortizes the duplicate suffix producers.

A bounded prepared-suffix batch owner was also prototyped and rejected before
commit. The trial added an opt-in
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_PREPARED_SUFFIX_BATCH=1` path
that, while a grouped prefix batch still owned its current states, advanced only
suffix producers shared by at least two owned currents and left singleton
producers for the normal later trigger path. The source built cleanly with
`cmake --build build -j2`, and the focused benchmark completed all five rows
with clean residuals, but it regressed sharply: the new path measured `4.0739s`
geomean in
`build/kls_prepared_suffix_batch_gap5_t4_r1_ref3_timeout120.jsonl` versus
`2.8150s` for the existing suffix-group path
`build/kls_prepared_suffix_batch_suffixgroup_gap5_t4_r1_ref3_timeout120.jsonl`
and `1.4142s` for the same-source default
`build/kls_prepared_suffix_batch_default_gap5_t4_r1_ref3_timeout120.jsonl`.
This rejects "batch only the already-owned prepared suffix frontier" as the
missing paper mechanism; it still spends too much effort on sparse retained
state setup and row-position validation before a coarse numeric owner exists.

The suffix-group executor now has a first producer-column batched implementation
behind the same opt-in flag
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_SUFFIX_GROUP_ADVANCE=1`.
The retained current-state row lists are sorted once during plan construction,
and the runtime suffix trigger now claims all eligible READY current states for
one completed producer, captures their U coefficient, then streams that
producer's L column across the batch. This is directly aligned with the paper
direction of amortizing one finished producer over a current window, and it
removes the earlier single-current trigger loop from the opt-in path.
Correctness passed `cmake --build build -j2`,
`ctest --test-dir build --output-on-failure`, `./build/kls_smoke`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_SUFFIX_GROUP_ADVANCE=1 ./build/kls_smoke`.
The final focused run
`build/kls_batched_suffix_group_gap5_t4_r1_ref3_timeout120_final.jsonl`
measured `2.7996s` top-five SPICE-cycle geomean. That is a real improvement
over the prior sorted single-current suffix-group probe
`build/kls_sorted_state_suffix_group_gap5_t4_r1_ref3_timeout120.jsonl`
at `4.1072s` geomean (`1.47x` speedup), with the largest improvements on
`ASIC_320k` and `ASIC_320ks`. It is still not a CKTSO-gap closer: the refreshed
same-binary default control
`build/kls_batched_suffix_default_gap5_t4_r1_ref3_timeout120_final.jsonl`
measured `1.4446s`, so the batched opt-in path is still `1.94x` slower than
default on the focused set. The result confirms batching the producer trigger
is necessary substrate, but the remaining paper gap is broader than this
one-producer sparse suffix update: KLS still needs a true grouped live-workspace
or supernodal window owner that avoids retained-state lookup/cursor work as the
dominant cost.

The suffix-group batch follow-up now retains the validated row positions for a
producer L column inside each claimed current-state batch, so the executor no
longer binary-searches the retained row list once for validation and again for
the update. This moves the opt-in path a little closer to a real grouped
numeric object: map rows once, stream numeric values once. Correctness passed
`cmake --build build -j2`, `ctest --test-dir build --output-on-failure`,
`./build/kls_smoke`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_SUFFIX_GROUP_ADVANCE=1 ./build/kls_smoke`.
The focused artifact
`build/kls_suffix_group_position_cache_gap5_t4_r1_ref3_timeout120.jsonl`
measured `2.7521s` top-five geomean, a small `1.017x` speedup over the prior
batched artifact (`2.7996s`). It remains far from the same-binary default
control at `1.4446s` (`1.91x` slower), so retained row-position caching is
useful substrate but not the missing CKTSO mechanism. The first-order gap is
still the owner/scheduler level: the path prepares and owns retained current
states per current cursor, while the papers' advantage comes from a broader
grouped live workspace or supernodal window that avoids this per-current
retained-state traffic.

The suffix-group window probe adds a separate opt-in
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_SUFFIX_GROUP_WINDOW=1` path
that keeps the claimed current-state batch alive after the first producer
column and advances later already-finished producers inside that owned window.
To keep the path aligned with the papers' shared-producer payoff rather than a
scalar tail replay, it only continues when a later producer has at least two
valid current-state consumers in the batch. Correctness passed
`cmake --build build -j2`, `git diff --check`,
`ctest --test-dir build --output-on-failure`, `./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_SUFFIX_GROUP_ADVANCE=1 ./build/kls_smoke`,
and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_SUFFIX_GROUP_WINDOW=1 ./build/kls_smoke`.
The focused artifact
`build/kls_suffix_group_window_fanout_gap5_t4_r1_ref3_timeout120.jsonl`
completed cleanly but measured `8.2147s` geomean, `2.98x` slower than the
position-cache suffix-group path and `5.69x` slower than the same-binary
default. The fanout guard reduced the unguarded window's suffix work, but the
hard rows still advanced too much retained-state traffic:
`ASIC_320ks` processed `228,582` suffix deps and `81,956,367` update entries,
`ASIC_320k` processed `337,750` / `141,881,301`, and `ASIC_100ks` processed
`11,048` / `3,980,511`. This rejects "keep a per-current retained cursor window
alive and push more suffix deps through it" as the CKTSO-gap closer. The next
refactor work should use the grouped producer/current-state surface without
materializing or cursor-advancing each current independently, e.g. the bounded
group-L state shape as an owned streaming executor.

The latest Algorithm 5 grouped-prefix work now includes a targetless
direct-prefix variant, an advance-seed probe, a retained-current-row state probe,
and a final-state probe, without changing BLAS thresholds. The
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_PREP=1` path
requests the selected payoff groups and runtime prefix workspace, but it does
not allocate runtime target slots. When a producer-prefix trigger fires, it
claims eligible current slots, prepares their prefix workspaces, applies only
the internal prefix triangular work, and later lets the ordinary ragged-L path
consume the ready prefix while streaming dense suffix and L-trailing updates
directly into the current column workspace. The newer
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_ADVANCE_SEED=1`
probe records post-advance row-workspace slots plus skipped U coefficients
from the same Algorithm 5 advance descriptors, so a consuming column can jump
over ready remaining advance dependencies before consuming the prepared prefix.
The newer
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_CURRENT_STATE=1`
probe stores and restores values for the compact per-current-slot row-state map,
plus skipped U coefficients. The newer
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_FINAL_STATE=1`
probe goes one step further: the grouped producer pass applies the prepared
prefix's dense suffix and L-trailing effects into the retained current-state row
map, so the consuming column restores that final state, writes the skipped
advance and prefix U coefficients, and jumps past the prefix. This fills another
direct CKTSO/SubtreeLU paper gap, but remains rejected as a default because the
retained-state restore volume is still too high in this executor shape.
The sparse-restore follow-up
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_SPARSE_RESTORE=1`
keeps the same final retained state but restores only rows needed by the
remaining scalar continuation, clearing unneeded input-scatter rows before
resuming the existing loop. A cheap structural upper bound rejects the sparse
row-set construction unless it can skip a material share of retained rows. This
is the first guarded attempt to remove the current full-row copy without yet
replacing the scalar continuation with a live grouped workspace.
The direct-complete follow-up
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_COMPLETE=1`
also uses the final retained state, but tries to complete the whole BTF-local
column from that state instead of returning to the scalar continuation. The
safe probe publishes terminal columns directly from the retained final state
when the prepared prefix covers all U dependencies; otherwise it computes on a
private retained-state copy and commits U, diagonal, and L values only after the
pivot is nonzero. The sparse-delta follow-up
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_SPARSE_DELTA=1`
uses the same direct-complete setup but replaces the private copy with a sparse
delta accumulator over the retained row map. The direct-complete remaining
stream now validates dependencies and retained-state rows as it computes into
private/local state, instead of pre-scanning the whole suffix and then scanning
it again. These probes prove the direct arithmetic is correct and make the copy
cost measurable, but they deliberately remain opt-in: the next lead cause is
still live grouped current workspaces rather than BLAS thresholds.
Correctness passed `cmake --build build -j2`,
`ctest --test-dir build --output-on-failure`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_PREP=1 ./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_ADVANCE_SEED=1 ./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_CURRENT_STATE=1 ./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_FINAL_STATE=1 ./build/kls_smoke`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_SPARSE_RESTORE=1 ./build/kls_smoke`;
the direct-complete follow-up also passed
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_COMPLETE=1 ./build/kls_smoke`
and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_SPARSE_DELTA=1 ./build/kls_smoke`.
The focused current-state value probe
`build/kls_direct_prefix_current_state_gap5_t4_r1_ref3_timeout120.jsonl`
completed all five rows with `2.2981s` SPICE-cycle geomean, worse than the
direct-prefix probe at `2.1306s`, the lighter advance-seed probe at `2.2868s`,
and the same-source default at `1.3810s`. The hard ASIC rows did exercise the
new path: `ASIC_320ks` restored `2,018` current states / `120,963` skipped deps
/ `1,792,538` state rows, `ASIC_320k` restored `1,972` / `113,755` /
`1,872,554`, and `ASIC_100ks` restored `346` / `16,109` / `320,274`. This
confirms the next paper-level executor should avoid per-current restore and keep
current workspaces live or batch producer-side advancement in place, rather than
revisiting BLAS thresholds first.

The final-state probe
`build/kls_direct_prefix_final_state_gap5_t4_r1_ref3_timeout120.jsonl` completed
all five focused rows with `2.2742s` SPICE-cycle geomean. It did skip the later
prefix consumer: `ASIC_320ks` recorded `2,018` final-state restores,
`226,588` skipped U deps including prefix rows, `1,792,538` restored state rows,
and `126,823` direct prefix rows with zero retained target slots. `ASIC_320k`
recorded `1,973` / `202,111` / `1,873,376` / `107,808`, and `ASIC_100ks`
recorded `350` / `27,565` / `323,524` / `11,938`. The timing stays in the same
bad range as the current-state restore probe, so the missing performance feature
is not another retained-state shortcut. It is the paper's live grouped
multi-current workspace execution, where producer-side prefix work is applied to
current workspaces that do not have to be copied back into the scalar column path.

The guarded sparse-restore run
`build/kls_direct_prefix_sparse_restore_guarded_gap5_t4_r1_ref3_timeout120.jsonl`
completed the same five focused rows with `2.2665s` SPICE-cycle geomean. The
guard mostly rejected sparse row-set construction on the ASIC cases:
`ASIC_320ks` recorded `2,018` restores / `226,588` skipped deps / `1,792,325`
restored rows, `ASIC_320k` recorded `1,973` / `202,111` / `1,873,119`, and
`ASIC_100ks` recorded `350` / `27,565` / `323,524`. The continuation set is
therefore nearly the whole retained state on the gap cases, so sparse restore is
not enough to close the CKTSO gap. The next direct paper-level target is avoiding
the scalar continuation after the retained prefix, for example by keeping grouped
current workspaces live through column completion rather than copying them into
the scalar scratch.

The private-state direct-complete run
`build/kls_direct_prefix_complete_private_state_gap5_t4_r1_ref3_timeout120.jsonl`
also completed the same five rows, but regressed to `3.3422s` SPICE-cycle
geomean. The ASIC rows solved correctly but paid the private retained-state
copy and remaining dependency stream: `ASIC_320ks` recorded `1,360` direct
completion runs / `702,618` skipped deps / `1,147,897` copied rows,
`ASIC_320k` recorded `1,357` / `742,805` / `1,224,196`, and `ASIC_100ks`
recorded `309` / `179,851` / `282,959`. This rejects per-column direct
completion with copied retained state as a CKTSO-gap closer. The missing paper
mechanism is still the grouped live workspace executor: Algorithm 5 current
states should stay owned by the grouped pipeline through completion, not be
copied into a private scalar/direct column path.

The terminal-live/private-copy follow-up
`build/kls_direct_prefix_complete_terminal_private_gap5_t4_r1_ref3_timeout120.jsonl`
completed all five rows with `3.1603s` SPICE-cycle geomean. Terminal no-copy
completion removed only a tiny retained-copy slice: `ASIC_320ks` still copied
`1,147,422` rows over `1,360` direct-complete runs, `ASIC_320k` copied
`1,223,719` rows over `1,357` runs, and `ASIC_100ks` copied `282,959` rows
over `309` runs. The sparse-delta A/B run
`build/kls_direct_prefix_sparse_delta_gap5_t4_r1_ref3_timeout120.jsonl`
also completed all five rows and reported zero copied retained rows for those
same runs, but measured `3.1847s` geomean. That rejects retained-row copy as
the sole large missing part: after copy removal, the hard rows still stream
`702,618`, `742,805`, and `179,851` remaining dependencies through a scalar
per-current path. The paper-level fix must therefore be a grouped owner that
advances multiple current workspaces and their remaining dependency stream
together, not just a lower-copy scalar direct-complete path.

The single-pass direct-complete follow-up
`build/kls_direct_prefix_complete_single_pass_gap5_t4_r1_ref3_timeout120.jsonl`
removed that duplicate suffix validation pass and completed the same five rows
with `2.3314s` SPICE-cycle geomean. `ASIC_320ks` improved to `16.5794s`,
`ASIC_320k` to `19.4593s`, and `ASIC_100ks` to `10.3019s` while preserving the
same direct-complete run/dependency surface. The sparse-delta single-pass A/B
`build/kls_direct_prefix_sparse_delta_single_pass_gap5_t4_r1_ref3_timeout120.jsonl`
also improved to `2.4101s` geomean with zero copied retained rows. This proves
the previous direct-complete loss was partly an executor artifact, not just
Algorithm 5 arithmetic. It still does not beat the current default top-five run
`build/kls_current_default_gap5_t4_r1_ref3_timeout120.jsonl` at `1.4629s`, so
the remaining paper gap is still broader grouped post-prefix execution rather
than promoting direct-complete as-is.

The guarded state-ragged suffix experiment
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_STATE_RAGGED_SUFFIX=1`
tries the next direct paper-level idea inside that direct-complete private state:
when the remaining U stream is contiguous and matches an existing ragged-L
supernode panel, it applies that grouped dense/trailing update against the
retained state map before falling back to the scalar dependency loop. Correctness
passed `KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_STATE_RAGGED_SUFFIX=1 ./build/kls_smoke`,
but focused ASIC probes reject this broad probing form. On `ASIC_320ks`,
`build/kls_direct_prefix_state_ragged_asic320ks_t4_r1_ref3_timeout120.jsonl`
measured `22.3606s` versus `16.6473s` for the same direct-complete path with the
suffix experiment off in
`build/kls_direct_prefix_complete_state_ragged_gated_fixed_gap5_t4_r1_ref3_timeout120.jsonl`.
The loss is structural: the opt-in suffix issued `187,757` ragged-L probes, but
`175,399` were panel misses and `3,848` were stream rejects. On `ASIC_100ks`,
the sequential diagnostic
`build/kls_direct_prefix_state_ragged_asic100ks_seq_t4_r1_ref3_timeout45.jsonl`
measured `12.0946s` and showed the same pattern (`95,890` probes, `92,865`
panel misses). This closes off another tempting micro-path: KLS should not
rediscover grouped suffix shapes one dependency at a time. The remaining paper
gap is still a retained multi-current Algorithm 5 executor that owns the grouped
workspaces and target slots up front.

Correctness passed `cmake --build build -j2`,
`ctest --test-dir build --output-on-failure`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_ADVANCE_SEED=1 ./build/kls_smoke`.
The focused probe `build/asic100ks_direct_prefix_prep_probe.json` reported an
EGraph refactor, runtime workspace rows `20,859`, runtime target slots `0`,
`350` consumed prefix-prep runs, `11,938` prefix rows, `3.65M` direct
suffix/trailing target entries, zero retained target slots, and `1.92e-15`
relative residual.
The matching advance-seed probe
`build/asic100ks_direct_prefix_advance_seed_probe.json` preserved the same
`1.92e-15` residual and reported `346` advance-seed runs, `16,109` skipped
advance deps, and `113,341` restored row slots, but refactor time rose to about
`0.099s` on the single probe.

The targetless direct-prefix executor is a partial paper-aligned improvement,
but still not the CKTSO closer. The top-five focused run
`build/kls_direct_prefix_prep_gap5_t4_r1_ref3_timeout120.jsonl` measured
`2.1306s` SPICE-cycle geomean. That improves the retained-target grouped probe
`build/kls_group_prefix_prep_gap5_t4_r1_ref3_timeout120.jsonl` at `2.3654s`
because the large ASIC rows no longer allocate or replay retained target slots,
but it is still slower than the same-source default
`build/kls_dense_direct_default_gap5_t4_r1_ref3_timeout120.jsonl` at `1.3810s`.
For example `ASIC_320k` recorded `1,973` targetless prefix-prep runs,
`107,808` prefix rows, `41.9M` direct suffix/trailing target entries, and zero
retained target slots. This matches the paper-level diagnosis: avoiding retained
targets helps, but KLS is still not executing CKTSO/SubtreeLU's full
row-workspace prefactor/postfactor pipeline with enough overlap and locality to
beat the default EGraph path.
The advance-seed top-five run
`build/kls_direct_prefix_advance_seed_gap5_t4_r1_ref3_timeout120.jsonl`
completed all five rows with `2.2868s` SPICE-cycle geomean. It did skip real
work (`ASIC_320ks`: `2,017` seed runs, `120,895` skipped deps, `664,375`
restored slots; `ASIC_320k`: `1,972`/`113,755`/`738,427`; `ASIC_100ks`:
`346`/`16,109`/`113,341`), but the row-slot restore traffic made it slower
than targetless direct-prefix prep. This rejects advance seeding as currently
implemented and sharpens the remaining paper gap: KLS needs a true shared
multi-current row-workspace pipeline, not a per-current restore of the
post-advance state.

The latest retained-target cleanup follows the Algorithm 5 descriptor more
directly without revisiting BLAS thresholds. In the scalar payoff slot-accum and
prefix-prep diagnostics, dense suffix targets now write by their retained slot
ordinal and only the irregular L-trailing targets build a row-stamp map. This
keeps `KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_PREP=1`
preparation-only rather than silently enabling the known-losing scalar replay.
Correctness passed `cmake --build build -j2`,
`ctest --test-dir build --output-on-failure`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_PREP=1 ./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_PREFIX_PREP=1 ./build/kls_smoke`,
and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_SLOT_ACCUM=1 ./build/kls_smoke`.
Focused top-five CKTSO-gap probes still reject scalar current replay:
same-source default
`build/kls_dense_direct_default_gap5_t4_r1_ref3_timeout120.jsonl` measured
`1.3810s` SPICE-cycle geomean, prefix prep measured `2.0942s` in
`build/kls_dense_direct_prefix_prep_gap5_t4_r1_ref3_timeout120.jsonl`, and slot
accumulation measured `2.1187s` in
`build/kls_dense_direct_slotaccum_gap5_t4_r1_ref3_timeout120.jsonl`. A
`GROUP_PREP` one-matrix check,
`build/asic100ks_group_prep_planonly_after_dense_direct.json`, reported zero
runtime payoff workspace and zero last-refactor prefix/slot-accum runs,
confirming it remains a descriptor path. The missing CKTSO-sized step is still
the compact grouped multi-current prefix-advance executor backed by the retained
advance slots, target slots, current slots, and current-run maps.

The latest Algorithm 5 prefactor change makes the guarded EGraph prefactor
slice automatic for retained EGraph pipelines with at least
`KLS_FAST_FACTOR_PIPELINE_REFACTOR_MIN_WORK` modeled dependency work, with
`KLS_ENABLE_EGRAPH_ALGORITHM5_PREF_UPDATE=1` and `=0` retained as explicit
force/disable A/B controls. The path now covers the unscaled BTF kernel, the
single-block unscaled/scaled kernels, and the generic scaled/fallback kernel.
This is a paper-semantics coverage change, not a BLAS or CPU-threshold change:
a blocked pipeline column can consume later already-finished U predecessors
when the structural safety scan proves that earlier unapplied predecessors
cannot still write that workspace entry. Correctness passed
`cmake --build build -j2`, `ctest --test-dir build --output-on-failure`,
`./build/kls_smoke`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_EXEC=1 ./build/kls_smoke`.
The focused top-ten CKTSO-gap A/B run now favors the modeled-work auto gate:
`build/kls_pref_auto_gap10_t4_r1_ref3_timeout120.jsonl` measured `2.2600s`
geomean, versus explicit prefactor-off
`build/kls_pref_default_off_gap10_t4_r1_ref3_timeout120.jsonl` at `2.4045s`
and forced-on `build/kls_pref_default_on_gap10_t4_r1_ref3_timeout120.jsonl` at
`2.4027s`. Compared with explicit off, the auto gate had seven wins and one
loss over 2%. Against CKTSO on the same ten rows, KLS is still `2.41x` slower
geomean (`2.2600s` versus `0.9376s`), so the main conclusion is unchanged: the
next first cause is the grouped multi-current Algorithm 5 executor rather than
BLAS thresholding, CPU-specific dispatch thresholds, or this scalar prefactor
slice.

The latest scalar payoff-exec rerun removes another ambiguity in the retained
paper path. A direct target-slot accumulator was tested and rejected before
commit: it reused the retained Algorithm 5 `group_target_cols` slot map inside
the existing scalar payoff ragged-L replay, but still processed one current
column at a time and paid an extra slot-map lookup on every update. The
top-five CKTSO-gap artifact
`build/kls_alg5_slotaccum_payoff_exec_gap5_t4_r1_ref3_timeout120.jsonl`
measured `2.1229s` geomean with
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_EXEC=1`, versus the same-source
default `build/kls_slotaccum_default_gap5_t4_r1_ref3_timeout120.jsonl` at
`1.4114s`. The ASIC rows did execute the retained payoff runs:
`ASIC_320ks`/`ASIC_320k`/`ASIC_100ks` reported `3,596`/`3,719`/`1,112`
last-refactor ragged-L runs and `55.7M`/`52.0M`/`5.1M` modeled update entries.
This confirms that the paper-level missing piece is not simply "let selected
scalar runs fire" or "publish scalar target slots"; it is the grouped
multi-current executor that advances many current workspaces and
applies/publishes a producer prefix as one batch.

The same conclusion holds after rerunning the existing producer-panel grouping
switches on the current auto-prefactor source. On the top-five CKTSO-gap focus,
the same-source default
`build/kls_slotaccum_default_gap5_t4_r1_ref3_timeout120.jsonl` measured
`1.4114s` geomean. The current opt-in grouped producer-panel executor measured
`1.6218s` in
`build/kls_current_group_exec_gap5_t4_r1_ref3_timeout120.jsonl`, grouped-batch
execution measured `1.6850s` in
`build/kls_current_group_batch_exec_gap5_t4_r1_ref3_timeout120.jsonl`, and the
shape-claim scheduler measured `1.6672s` in
`build/kls_current_shape_claims_gap5_t4_r1_ref3_timeout120.jsonl`. These runs
all stayed on the no-CBLAS EGraph path for the hard ASIC rows and still lost.
The useful paper target is therefore not another BLAS guard and not the current
single-current/producer-panel replay switches; it is persistent grouped
Algorithm-5 state with one batch advancing many current workspaces.

The latest retained-plan change moves the grouped-output direction from a
counter-only observation to retained data. With
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_SHAPE_TARGETS=1`, KLS now stores
per-producer-shape target-row maps for the compact publish surface that a
multi-current output accumulator would need. This is still off-by-default and
does not alter numeric execution, but it confirms the scale of the paper-level
opportunity: `ASIC_100ks` collapsed `133,615,959` raw shape-publish entries to
`159,004` retained target rows, and `ASIC_320k` collapsed `109,078,358` raw
entries to `75,161` target rows in focused one-refactor probes. The next
executor should consume these retained targets directly; BLAS thresholding is
not the first cause from this point.

The clean rebuild rerun after the rejected selector probe confirms that BLAS is
not the first-order explanation from this point forward. The same source rebuilt
from the clean tree measured `1.4454s` top-five CKTSO-gap geomean in
`build/kls_current_rebuilt_gap5_t4_r1_ref3_timeout120.jsonl`, while the existing
Algorithm-5 payoff-claims scalar replay path measured `1.8426s` in
`build/kls_alg5_claims_rebuilt_gap5_t4_r1_ref3_timeout120.jsonl`, with `rajat03`
regressing most sharply. Forced KLS-first row-up-looking factorization also lost
at `1.6118s` in `build/kls_klsfirst_forced_gap5_t4_r1_ref3_timeout120.jsonl`.
The component comparison against
`build/cktso_paper_medium93_t4_timeout120.jsonl` shows the SPICE-cycle gap on
the hard ASIC rows is dominated by repeated numeric refactorization: KLS solve
time is close on `ASIC_320ks`, while refactor time is roughly `3.6x` CKTSO and
initial factorization is even farther behind but amortized over the 99 repeated
refactors. KLS fill is only modestly worse on those rows, so the next paper-level
work should stay on row/producer-panel grouped numeric execution and
checked-tail scheduling, not BLAS thresholds or matrix-specific ordering tweaks.
KLS now also has an explicit
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_PREP=1` preparation flag
for that work: it builds the retained payoff-group plan, and now only builds the
Algorithm-5 U-supernode pattern when the payoff plan actually finds groups,
without enabling the known-losing scalar
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_EXEC=1` replay path. This keeps
future grouped-executor experiments separated from scalar current replay and
keeps no-group cases from paying retained-pattern setup cost under the prep flag.
Correctness passed `cmake --build build -j2`,
`ctest --test-dir build --output-on-failure`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_PREP=1 ./build/kls_smoke`.
The focused top-five diagnostic
`build/kls_alg5_group_prep_gap5_t4_r1_ref3_timeout120.jsonl` populated the hard
ASIC grouped surfaces (`129` groups / `3,596` currents on `ASIC_320ks` and `102`
/ `1,112` on `ASIC_100ks`) while reporting zero last-refactor ragged-L update
runs, so it is preparation only. Its `1.7807s` geomean is not a speed claim; the
extra retained-plan/schedule work must be consumed by a real grouped numeric
executor before this path can close the CKTSO gap.

The June 30, 2026 advance-slot diagnostic adds the missing row-surface
measurement for that executor. KLS now records, per retained Algorithm-5 payoff
group, how many distinct local workspace rows would be touched by advancing the
selected current columns to the producer prefix. This changes stats only, not
numeric execution. Correctness passed `cmake --build build -j2`,
`ctest --test-dir build --output-on-failure`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_PREP=1 ./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_PLAN=1 ./build/kls_smoke`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_EXEC=1 ./build/kls_smoke`.
The focused artifact
`build/kls_alg5_advance_slots_gap5_t4_r1_ref3_timeout120.jsonl` reports that the
hard ASIC advance surface is comparable to the target surface, not a full-block
explosion: `ASIC_320ks` has `909,325` advance slots versus `784,673` target
slots, `ASIC_320k` has `1,065,691` versus `971,718`, and `ASIC_100ks` has
`194,019` versus `188,230`. The largest single run needs `1,884` advance slots.
This makes the plausible paper-aligned next fix a compact grouped
multi-current prefix-advance executor backed by the retained group/current/run
maps. It also reinforces that BLAS thresholding is not the first cause from
here: the issue is still repeated scalar workspace advancement and publication,
not small dense kernel dispatch.

KLS now retains those compact advance rows instead of only counting them. The
payoff-plan post-pass materializes per-selected-run advance slot pointers and
row lists plus per-group advance slot counts, reusing the same row surface
validated by the diagnostic above. This is still a descriptor/substrate change,
not a promoted numeric executor: the existing scalar payoff replay remains
guarded and known-losing, but the next grouped executor no longer needs a
full-block workspace or hot-path row discovery to advance selected current
columns. Correctness passed `cmake --build build -j2`,
`ctest --test-dir build --output-on-failure`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_PREP=1 ./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_PLAN=1 ./build/kls_smoke`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_EXEC=1 ./build/kls_smoke`.
The focused prep artifact
`build/kls_alg5_advance_rows_gap5_t4_r1_ref3_timeout120.jsonl` completed all
five rows with `1.7659s` geomean and preserved the same hard-ASIC surfaces:
`909,325` retained advance slots on `ASIC_320ks`, `1,065,691` on `ASIC_320k`,
and `194,019` on `ASIC_100ks`.

KLS now also retains the producer-keyed advance surface behind
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_ADVANCE_MAP=1`. This mirrors the
suffix producer map but covers the pre-prefix dependencies that each selected
current must advance before the grouped prefix panel can run. The descriptor
records `refactor_supernode_algorithm5_payoff_group_advance_unique_deps`,
`...advance_duplicate_deps`, `...advance_shared_deps`,
`...advance_max_dep_fanout`, and
`...advance_duplicate_update_entries`, plus the opt-in map arrays
`group -> advance producer column -> (run, U position)`. Correctness passed
`cmake --build build -j2`, `git diff --check`, `./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_ADVANCE_MAP=1 ./build/kls_smoke`,
and `ctest --test-dir build --output-on-failure`. The focused descriptor run
`build/kls_alg5_advance_map_gap5_t4_r1_ref3_timeout120.jsonl` exposes a large
shared pre-prefix advance target: the two large ASIC rows have about `70%`
duplicate advance producer occurrences inside their retained payoff groups
(`137,236/190,994` on `ASIC_320ks` and `132,537/194,219` on `ASIC_320k`). This
fills the descriptor gap left by the row-slot plan; the remaining implementation
gap is the actual grouped advance executor that consumes this producer map.

The EGraph prefactor path now has the corresponding guarded supernode-shaped
prefactor slice when `KLS_ENABLE_REFACTOR_U_SUPERNODE_RAGGED_L=1` is explicitly
requested. A blocked BTF EGraph column may consume a contiguous finished
U-supernode producer run during Algorithm-5 prefactorization, but only after a
run-aware safety scan proves that no earlier unapplied predecessor can still
write any row in that run. This fills the paper's "use finished producer
supernodes while skipping unfinished predecessors" semantic without promoting
the losing ragged producer cache to the default path. Correctness passed
`cmake --build build -j2`, `ctest --test-dir build --output-on-failure`,
`KLS_ENABLE_REFACTOR_U_SUPERNODE_RAGGED_L=1 ./build/kls_smoke`, and
`KLS_ENABLE_REFACTOR_U_SUPERNODE_RAGGED_L=1 KLS_ENABLE_EGRAPH_ALGORITHM5_PREF_UPDATE=1 ./build/kls_smoke`.
The focused top-five rerun still rejects this as a CKTSO-gap closer:
`build/kls_pref_supernode_default_gap5_t4_r1_ref3_timeout120.jsonl` measured
`1.4631s` geomean, while the opt-in ragged-prefactor run
`build/kls_pref_supernode_ragged_gap5_t4_r1_ref3_timeout120.jsonl` measured
`1.6197s`. The opt-in path executed retained ragged producer work on
`ASIC_320ks`, `rajat03`, and `ASIC_100ks`, then auto-disabled after low useful
coverage. This confirms the missing large piece is still the multi-current
workspace/accumulator executor, not another scalar or single-workspace
producer-cache variant.

A structural selector probe was also rejected. Raising the exact-EGraph work
floor for the small compact dominant-BTF class from `5.0e5` to `5.0e6` moved
`rajat03` from EGraph to the mapped path, but worsened the focused top-five
CKTSO-gap geomean to `1.4493s` in
`build/kls_smallcompact_workfloor_gap5_t4_r1_ref3_timeout120.jsonl`. The
`rajat03` repeated refactor itself slowed from the recent EGraph baseline around
`0.00069s` to `0.00103s`. This rules out a broad "small compact EGraph is the
problem" selector fix for the observed gap.

The latest forced-row rerun adds explicit row-group shape/work diagnostics and
rejects the small-BLAS hypothesis for the current slow cases. The slow default
ASIC rows still report zero CBLAS update counters, and forced row refactor with
`KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=0` shows that the row work is not trapped
primarily in tiny independent batches: `ASIC_320ks` has about 270.5M of 285.2M
modeled row work in dense groups, `ASIC_100ks` about 308.0M of 319.1M, and
`onetone1` about 415.1M of 421.3M. Enabling the native/compact dense-panel path
does execute hundreds of retained panels, but only modestly improves the
forced-row top-five geomean from about 4.33s to about 4.03s and still leaves it
far behind the current default top-five around 1.40s. The missing paper-scale
piece is therefore a broader CKTSO/SubtreeLU-style dense producer-panel numeric
executor that reduces row-panel update traffic and synchronization, not another
guard around BLAS calls for small cases. The artifacts are
`build/kls_group_shape_diag_forcedrow_gap10_t4_r1_ref3_timeout120.jsonl`,
`build/kls_group_shape_diag_forcedrow_native_auto_gap5_t4_r1_ref3_timeout120.jsonl`,
and `build/kls_current_gap10_continuation_t4_r1_ref3_timeout120.jsonl`.

The June 29, 2026 Algorithm-5 split diagnostic makes the next SubtreeLU gap
more concrete. With `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN=1`, the
top-ten CKTSO-gap artifact
`build/kls_alg5_split_diag_gap10_t4_r1_ref3_timeout120.jsonl` reports 892
producer panels wide enough for the paper's `2P` split rule, 19,325 covered
panel rows, and 15,757 prefix rows before the `P`-row suffix. Their retained
consumer-plan runs cover 156,935 candidate runs and 4,484,827 run rows. The
prefix side alone covers 156,660 runs, 3,921,698 run rows, and about 1.397B
modeled update work, compared with only about 0.136B work in the earlier exact
identical-prefix batch diagnostic. That prefix work is about 79.7% of retained
plan `L` work and 71.4% of retained `L+internal` work across these ten rows.
This is a much broader paper-aligned target than another BLAS-size guard: KLS
needs a producer-panel prefix/suffix executor with readiness tracking and
multiple consumer workspaces, not a narrower exact-prefix replay.

The follow-up payoff diagnostic narrows the executor shape further. The final
top-ten run
`build/kls_alg5_payoff_diag_gap10_t4_r1_ref3_timeout120.jsonl` still reports
the same 156,660 Algorithm-5 prefix runs and about 1.397B modeled update work,
but advancing every current column to those producer offsets would cost about
11.29B modeled work. Whole producer panels are therefore the wrong unit:
only 3 full panels are payoff-positive. The useful target is the low-advance
subset inside each producer panel: 594 producer-panel subsets contain 13,980
payoff-positive runs, 670,921 prefix rows, about 251.4M update work, and about
132.9M advance work, with up to 322 selected runs in one panel. That is still
a materially broader target than exact-prefix batching, but it is not enough
to justify replaying complete panels or expensive scalar current tails. The
next executor should select low-advance current columns inside each retained
Algorithm-5 panel, advance those workspaces together, apply the shared producer
prefix, and leave high-advance currents on the scalar path.

KLS now retains that missing executor shape explicitly when
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_PLAN=1` or
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_EXEC=1` is set. The retained
payoff plan stores distinct current columns per selected producer group and a
compact workspace-row estimate for the future multi-current batch executor.
It now also keeps each selected run's current-slot index plus a compact
workspace pointer for that current slot, so the future executor has the direct
gather/publish layout it needs rather than only aggregate counters.
The focused top-five rerun
`build/kls_alg5_workspace_map_gap5_t4_r1_ref3_timeout120.jsonl` populated
those counters without changing default execution: `ASIC_320ks` has `129`
multi-current groups, `3,596` distinct current workspaces, max `196` currents
in one group, and `157,067` compact workspace rows; `ASIC_320k` has `121`
groups, `3,755` workspaces, max `233` currents, and `146,027` workspace rows.
This confirms the paper-level target is a real multi-current workspace
executor, not a BLAS size guard or another single-current ragged replay.
The next trigger diagnostic now wires producer-prefix completion to the retained
Algorithm-5 payoff descriptor with
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_CLAIMS=1`. It is deliberately
surface-guarded: a naive unguarded run serialized large Sandia future-current
groups behind one producer worker and timed out, which directly confirms that
the paper gap is the shared multi-current work queue/batch executor. With the
guard, large Sandia rows keep plan counters but record zero claims, while
`rajat03` records 12 payoff groups, 76 current workspaces, and 4 claimed
columns. The guarded top-five CKTSO-gap artifact
`build/kls_alg5_claims_plan_guard_gap5_t4_r1_ref3_timeout120.jsonl` completed
with no failures and a `1.397s` geomean versus `1.430s` for the same-session
default artifact, but this should be read as bounded trigger validation rather
than a general CKTSO gap closer.

KLS has **not** implemented every paper idea that is still worth trying. It has
implemented the ideas that can be layered around the current KLU-derived
Gilbert-Peierls kernel: BTF, AMD/COLAMD/METIS ordering policy, explicit SCOTCH
ordering, CAMD refinement, auto scaling policy, pivot-checked reuse, static
row-pivoting trials with dual-potential matching-derived equilibration, and
BTF-block parallel refactorization with a solver-owned worker pool. KLS now
also keeps precomputed refactor scatter metadata for unscaled serial repeated
refactors and for a narrow scaled many-fringe dominant-BTF subset, records an
exact no-pivot EGraph level schedule from the numeric U pattern, and consumes
that schedule in a guarded cluster/pipeline intra-block refactor path for large
single-block matrices, including KLU row-scaled cases whose scale factors can
be recomputed and permuted safely, and inside large dominant BTF blocks whose
off-block entries can be refreshed from the retained map. A narrow many-fringe
dominant-BTF shape can also consume the same exact EGraph in an all-pipeline
mode, and very high-work unscaled single-block factors can do the same,
avoiding cluster barriers and waiting only on actual U-pattern predecessors.
KLS also now schedules extremely fragmented unscaled BTF refactors with one
substantial but non-dominant block through the exact EGraph path when measured
work is high enough, so the large block is not left serial behind hundreds of
thousands of singleton blocks. Auto scaling now starts those full-rank
fragmented BTF shapes in no-scale/no-recheck mode, because KLU row scaling
prevents the EGraph refactor from running and roughly doubled repeated refactor
time on the ASIC 680k/680ks structural class.
The row-permuted solve wrapper now keeps a solver-owned dense permutation
workspace instead of allocating it on every forward or transpose solve, trimming
repeated SPICE-cycle solve overhead for static-pivoting cases while leaving the
KLU triangular kernels unchanged.
Solve calls also now update only solve timing, status, and memory counters
instead of refreshing the full numeric stats snapshot after every triangular
solve.
The barriered EGraph
cluster levels now use FLOP-estimated per-thread slices instead of equal column
slices, and the no-pivot pipeline tail now uses an atomic dynamic work cursor
instead of static per-thread strides, which are retained pieces of the
CKTSO/SubtreeLU load-balance idea.
The same schedule pass now records adjacent row-major U-pattern supernode
candidate counts, covered rows, maximum width, dense-block entries, and shared
trailing entries. The gated row refactor also reports executable segment groups,
segment-covered rows, segment maximum width, and dense/trailing work. These are
diagnostic bridges toward SubtreeLU-style row/segment storage and BLAS-friendly
updates. Dense-eligible gated row-refactor groups now also compute their
internal dense `L` block, dense `U` block, and shared trailing `U` block
directly in KLS-owned row-major mirrors before scattering back to KLU on
success. Row-refactor groups now also retain a reverse group graph and report
group dependency edges, root groups, leaf groups, and maximum group fanout,
giving the row/segment layer the explicit task graph needed by future
SubtreeLU-style private/pipeline queues and CKTSO-style tail schedulers. The
experimental row pipeline can now consume that graph through a bounded
successor-ready queue for the pipeline tail, and reports whether the last row
run used that queue plus cumulative queued group counts. Unchecked queued row
tails also skip the per-row completion bitmap, while checked row runs still keep
it for prefix-reject validation and report that usage separately. Row-pattern
analysis also records when input columns are already structurally covered by
`L`, the pivot, or `U`, allowing row kernels to skip redundant residual cleanup
loops. Dense row segments now keep deferred all-value scatter only for checked
pivot probes, while unchecked row refactors scatter completed dense rows
directly. Dense and sparse row segments now retain per-input-entry target maps
for rows whose raw values can be placed exactly into external work-vector
slots, row-major `L` mirrors, the pivot, or in-segment/trailing `U` mirrors;
the direct loaders consume those maps on repeated refactors instead of
rediscovering the segment slots. The row ready queue now also orders initially
ready tail groups and newly released successor groups by the retained
FLOP-style group work estimate, which moves the row scheduler closer to
SubtreeLU's workload-balanced queue generation. Checked row fast-factor rejects
now also report a
conservative group-tail restart scope from the retained row-group successor
graph, giving future CKTSO-style pivoting tail work an explicit row/segment
task-tail measurement. Checked KLS-owned row fast-factor/refactor passes also
test the guessed diagonal against the maximum absolute value already present in
the current U row before publishing it, matching CKTSO's row-wise pivot check
rather than relying only on later L-multiplier growth. The row ready queue also
keeps solver-owned workspace across repeated row refactors, avoiding
queue/bitmap/predecessor allocation churn in the experimental row scheduler.
KLU-compatible checked fast-refactor paths now enforce the same row-wise U
acceptance rule after the finalized U rows are available; a violation is
reported as an all-current fast reject and enters the existing KLS block/tail
repair path. This closes the pivot-acceptance mismatch for the column kernels,
but it is still delayed validation rather than the row-major CKTSO numeric
executor that can avoid doing rejected work in the first place.
Full-graph queued row runs now consume cached
root groups through a private-root cursor before falling back to the shared
ready queue for successor-released groups, trimming the first wave of shared
queue traffic without changing the retained row DAG. When a completed group
releases successors, the worker now keeps one ready successor as a local
continuation and enqueues the rest, moving the queued DAG scheduler another
step toward SubtreeLU-style private/pipeline execution; stats now expose how
often that local continuation path is used. Checked queued rejects now refresh
any missing prefix rows before accepting the prefix-tail repair classification,
so work ordering cannot turn an already-repairable prefix into a scheduler-race
miss. Experimental row refactors now cover single-block
factors and BTF
diagonal blocks; BTF off-block values are refreshed into KLU `Offx` from the
retained input map. Unchecked row refactors can also hand dirty KLS-owned
row-major `L`/`U` mirrors directly to guarded forward and transpose
solves, leaving KLU's column values stale until later factorization or other
non-row fallback needs them. For BTF factors the dirty row solve follows KLU's
block order and uses refreshed `Offx` coupling directly. This avoids making
eligible row-refactor solves pay an immediate KLU publish. KLU row-scaled row
refactors can now also recompute `Rs`, use unpermuted row scales for
fixed-position input loads, and restore `Rs` to pivot order after an accepted
pass. Dirty row solves consume both unscaled and KLU row-scaled normal
row mirrors directly, using KLU's pivot-order `Rs` semantics for the
right-hand-side load and transpose output. The guarded row-major solve also
runs under KLS's external static-matching row permutation and matching-derived
row/column scaling because `solve_impl` wraps the kernel solve with the same
right-hand-side/result transforms used by KLU. These pieces still do not change
the current default KLU-column numeric kernel.
KLS now also records `initial_factor_path`, `last_factor_path`, and largest
ordered diagonal-block factor ETree counters in `kls_stats`/`kls_bench`. These
diagnostics expose whether a solve still entered the KLU-derived pivoting
kernel (`klu_first` or `klu_fallback`) and the ETree upper-bound shape CKTSO
uses for pivoting-tail scheduling; they are evidence for the remaining
row/up-looking first-factor work, not a substitute for that kernel.
An env-gated `KLS_ENABLE_KLS_FIRST_FACTOR=1` path can now allocate and assemble
KLU-compatible numeric storage itself for no-scale and KLU row-scaled first
factorizations, handle singleton BTF blocks directly, and use the KLS-owned
pivoted block kernel for multi-column BTF blocks. It reports `kls_first` when it
succeeds, seeds KLS-owned row-major `L`/`U` value mirrors for guarded solves and
unchecked and checked repeated refactors, and falls back to the KLU first-factor
path otherwise. This is a first KLS-owned factorization scaffold, not the
default production row-major CKTSO-style factorization.
When this row up-looking first-factor path performs a dynamic column pivot and
the retained METIS `NodeNDP` map covers the full factor order, it now uses
SubtreeLU Algorithm 4's scoped `N'` maximum inside the current collapsed
component extent instead of the whole updated row; if no safe scoped candidate
exists for a required pivot, it rejects the row-up attempt instead of crossing
the separator pivot domain. When a separator-pipeline row has been safely
pre-updated and then needs a scoped pivot, the pipeline now publishes that row
under the ordered lock, updates both the committed and phase-local prefix U
entries for the column exchange, advances a column-order epoch, and makes
speculative suffix rows that began under an older epoch discard and restart
inside the same worker phase instead of tearing down and relaunching the whole
pipeline suffix.
Benchmark stats separate separator-domain dynamic pivots, would-have-crossed
fallback candidates, and strict rejects.
Its static-pivot
preprocessing has a cheap exact sparse maximum-log-product assignment path for
small candidates, derives MC64-style row/column scales from that exact
assignment's dual potentials when the scaled candidate is accepted, and can
improve medium row matchings with bounded alternating cycles beyond the
pair-swap pass. KLS now enables the pinned BSD-licensed SPRAL
scaling subset by default, so it can use Hungarian matching/scaling as a
pre-factor MC64-adjacent candidate for small and medium weak-diagonal matrices,
SPRAL auction matching/scaling for very large weak-diagonal dominant-block
matrices, and a post-factor value-gated Hungarian trial for dense
high-off-diagonal-pivot cases while still allowing
`KLS_ENABLE_SPRAL_SCALING=OFF` builds. Fast factorization can now repair an
unsafe unscaled BTF diagonal block by
restarting that block with pivoting and then retrying the checked no-pivot
factorization when later columns may not have been refreshed. If the failed
pass had refreshed all columns, or if a serial BTF pass rejected in the final
block, KLS validates the pivoted repair block tail directly and skips the
redundant checked refactor pass.
The same block-local repair now also covers KLS-owned scaled checked-refactor
failures; when the rejected pass left a valid prefix-current state, KLS
recomputes row scales and continues with the threaded checked BTF pool over only
later BTF blocks, falling back to the serial checked continuation when the pool
is not applicable. Scaled prefix-current rejects can also use the conservative
serial suffix restart inside the rejected block when the repaired block
preserves a validated non-empty live prefix. The preserved-column refresh used
by masked ETree-tail repair now consumes the same input-row `Rs` state, so
scaled exact-tail attempts can make omitted independent columns current instead
of immediately widening to block rebuild. Root rejected pivots can now use that
same scaled refresh to derive the shorter reject-only ETree closure instead of
being pinned to the prefix-current full suffix. Scaled all-refresh KLU-refactor
rejects now recompute the row scale vector back to input-row order before
attempting the same block-local repair and validated non-root serial suffix
restart. The newer row-first rejected-block rebuild consumes that same
input-row scale state and permutes `Rs` back to pivot order after accepting the
KLS-owned dynamic-column-pivot block, so scaled checked rejects are no longer
limited to the serial KLU-compatible block repair before the KLU fallback. If
that local repair ladder still fails, scaled rejects can now attempt the same
quality-checked KLS-first whole-numeric rebuild as unscaled rejects.
The threaded BTF worker pool now keeps a per-run completed-block bitmap and
reports prefix-current only when every diagonal block before a checked pivot
reject has finished. That lets the existing block repair and later-block
checked continuation consume more safe parallel fast-factor rejects without
claiming the full CKTSO tail algorithm.
Checked EGraph refactors now also retain the completed-column bitmap for
barriered cluster-only runs, not only for all-pipeline or pipeline-tail
schedules. A cluster-mode pivot reject can therefore be classified as
prefix-current when every earlier factor-order column is proven finished,
while same-level out-of-order cases remain conservative.
Strict tail-restart readiness now also validates that the repaired block
preserved the old prefix pivot order and that KLS can reconstruct the live KLU
prefix state a pivoting tail kernel would need, including finalized-L row
unfinalization, live `P`/`Pinv`, and symmetric-pruning `Lpend` boundaries.
KLS now uses that proof to run a conservative serial suffix restart with
pivoting for unscaled and scaled prefix-current/all-current
rejected blocks with a non-empty reusable prefix before falling back to full
block repair. Accepted serial suffix restarts now refresh only the off-diagonal
column suffix whose inverse row permutation may have changed, including the
scaled all-current subset where KLS first recomputes the row-scale vector back
to input-row order before the final pivot-order permutation.
Root-of-block rejects can now use the same KLS-owned pivoted block kernel, but
they are deliberately counted as full KLS block restarts rather than
serial-tail-ready events because no contiguous prefix can be reused. This is
still not CKTSO's full pipelined ETree-descendant tail factorization.
SPRAL static-pivot matches that are already close to complete can now be
finished with KLS's existing nonzero structural augmenting-path graph before
the row permutation is accepted. The SPRAL scaling vectors are retained only
when SPRAL itself found a complete weighted match; structurally augmented
matches fall back to KLS's matching-equilibration path so the scales remain
consistent with the final diagonal choice.

The remaining CKTSO gap is large enough that it should be treated as a missing
major algorithm, not an ordering-package tuning problem. On the selected large
`pre2` case, CKTSO completed the factor-only comparison inside the 120s cap
with about 28.4s cycle time, while KLS timed out under AMD/BTF, METIS/no-BTF,
SCOTCH/BTF, and auto policy variants. The CKTSO ordering supplement reports
`pre2` operation counts for CKTSO nested dissection and METIS in the same
range. After switching the very large pre-static MC64-adjacent path from exact
SPRAL Hungarian matching to SPRAL auction matching/scaling, `pre2` still times
out under the same 120s cap. Forced zero pivot tolerance fails as singular, and
`1e-5`/`1e-4` pivot tolerances still time out. This makes full MC64-quality
matching/scaling worth keeping, but not sufficient as the next expected gap
closer by itself. The current evidence points most strongly at CKTSO's
KLS-owned row/up-looking numeric factorization, EGraph fast factor with pivot
checks, and ETree-descendant pipelined tail restart machinery.
After rereading the local papers in `refs/`, KLS now also uses METIS's
`METIS_NodeNDP` entry point for threaded METIS analyses of matrices with at
least 30,000 rows. This asks vendored METIS for at least `log2(threads)` nested
dissection levels, matching SubtreeLU's minimum-depth separator-tree setup more
closely than plain `METIS_NodeND`, while retaining the existing CAMD refinement
and auto-ordering policy. A focused 12-row CKTSO-gap run showed this is neutral
as a standalone ordering change: 3.1772s geomean versus the saved 3.1671s
baseline, with no failures. KLS now retains the accepted `NodeNDP` top-level
component sizes as private leaf domains plus pipeline separator components, so
this is no longer only an ordering precursor. The experimental row-refactor
ready queue and KLS-first row-up-looking pivot path can consume those retained
components even when the separator map covers the dominant BTF block rather
than the full matrix, provided the accepted symbolic range can be matched
uniquely. This is still not evidence that ordering alone closes the gap because
the retained separator queues only became useful after KLS connected them to
the row numeric consumers; the remaining gap is no longer the queue metadata
itself, but CKTSO's complete pivoting ETree-tail scheduler and production
row/supernode numeric kernels.
The CKTSO paper in `refs/` is explicit that CKTSO's core factorization is a
row-major sparse up-looking factorization, and that the fast path combines
guessed EGraph pivot-checked refactorization with an ETree-scheduled pipelined
tail factorization when repivoting is needed. KLS no longer always hands large
eligible first factorizations to KLU: the conservative automatic `kls_first`
path can run a row-up-looking first factor and now seeds row-refactor mirrors
directly from those row entries. The remaining gap is that accepted factors are
still packed into a KLU-compatible numeric object for fallback coherence, and
KLS does not yet have CKTSO's fully row-major parallel first-factor/EGraph
executor. A debug trace on `pre2` confirmed the timeout occurs in
`trilinos_klu_l_kernel` from the fallback first-factor call. The same trace
showed SPRAL auction static pivoting finds 633565 weighted matches for 659033
rows; structurally augmenting that to a full
permutation made the static-pivoted AMD symbolic estimate worse than the
baseline, so it did not expose a credible low-risk static-pivot-only fix for
this slow case. A broader nested-dissection retry on that static-pivoted
dominant-BTF shape spent the probe timeout inside ordering analysis, so it was
not retained as a default policy.
The KLS-owned row-up first-factor scaffold now retains per-column entry counts
while generating each block, so packing no longer rescans all generated `L` and
`U` entries just to compute KLU column lengths. This is useful staging work for
native row/segment storage: on sampled KLS-first probes it moved `nxp1` initial
factor time to about 3.50s and `G2_circuit` to about 1.97s. The same scaffold
now also reserves its generated row-entry buffers from KLU's symbolic block
fill estimate before starting numeric updates, avoiding repeated large
realloc/copy waves on high-fill blocks. A focused rerun moved `nxp1` initial
factor time further to about 2.93-3.01s; `G2_circuit` stayed in the same rough
range. It did not close the CKTSO paper gap. `pre2` still times out under a
120s factor-only cap with `KLS_ENABLE_KLS_FIRST_FACTOR=1`, and `ASIC_320k`
repeated refactor remains on the column EGraph path because the current
row-group work model is still higher than the exact EGraph work. Therefore the
next direct CKTSO-aligned step is still a production row/segment numeric
engine, not simply enabling the current row-up scaffold by default.
The KLS-first pivoted fallback can now factor independent BTF diagonal blocks
concurrently when the row-up scaffold rejects. Each worker uses private
KLU-kernel scatter/workspace, a private precomputed `Offp` seed to avoid
boundary races between adjacent BTF blocks, and commits only block-local
`LUbx`/pivot/statistics into the KLS numeric object. Final off-diagonal entries
are still rebuilt from the accepted global pivot order. This closes the safe
BTF-level "factorization with pivoting over independent work" gap for the
current KLU-compatible numeric object, but it is still not CKTSO's
single-large-block ETree-descendant pipelined tail factorization.

The retained broader SPRAL post-factor trial is deliberately value-gated. A
plain broad gate improved several MC64-sensitive cases but regressed the medium
corpus because rejected or small accepted trials added setup cost. The retained
gate requires dense off-diagonal pivot evidence and accepts only when the SPRAL
candidate removes substantial pivoting pressure or materially reduces
factorization work/fill. On the 93-matrix medium paper corpus this moved KLS
geomean from 0.32198s to 0.31410s with the same three known failures. The new
exact-match wins were `hvdc1`, `OPF_10000`, `LeGresley_87936`, `rajat22`,
`rajat23`, `rajat24`, `mult_dcop_03`, and `TSOPF_FS_b39_c19`.
Benchmark artifacts now report whether METIS, SCOTCH, and SPRAL scaling were
compiled into the tested binary, and the paper-gap benchmark runner can require
SPRAL scaling explicitly. This matters for interpreting CKTSO/NICSLU-style
static-pivoting evidence: a no-SPRAL build on `rajat24` falls back to the
slow `metis`/KLU-scaling path with thousands of off-diagonal pivots, while the
SPRAL-enabled build selects the low-fill AMD unscaled static-match path with
20 off-diagonal pivots. The remaining top-row losses after that correction are
still dominated by repeated refactor throughput, not by the missing matching
hook.
Benchmark artifacts now also report `internal_index_bytes`, because KLS still
uses the SuiteSparse long-index numeric object internally even when callers
provide 32-bit input arrays. A same-machine KLU width diagnostic on current
slow rows showed a broad index-cache signal: system KLU32 refactor time was
about `0.73x` KLU64 on `ASIC_100ks` and about `0.74x` KLU64 on `G2_circuit`.
That does not close the full CKTSO gap by itself, but it is large enough to
keep a future dual-width KLS numeric backend on the major-work list alongside
the row/segment engine, rather than treating index width as a cosmetic API
detail.
A refreshed 12-row CKTSO-gap run against the saved 93-row CKTSO medium-paper
artifact makes the same point on the current code path. With `--kls-first-factor
on` and row solves enabled, the common-row KLS/CKTSO SPICE-cycle geomean ratio
is 2.4855, with all 12 rows still losing. The largest ratios are
`ASIC_320ks` at 3.267, `ASIC_320k` at 3.225, `G2_circuit` at 2.810, and
`ASIC_100ks` at 2.800. These rows report `last_factor_path=kls_fast_refactor`;
most use `initial_factor_path=kls_first`, while several use the pre-static
first-factor path. Their cycle time is dominated by repeated numeric refactor
time, not triangular solve time or first-factor fallback. On the sampled hard
rows, ordering sweeps left `auto` already choosing the best available AMD,
METIS, or SCOTCH candidate, and the row-refactor auto gate correctly stayed off
because the current row-group work estimate exceeded the exact EGraph work.
That combination makes the missing piece a numeric/storage algorithm, not an
untried ordering package or a one-matrix policy rule.
A newer top-10 CKTSO-gap run on the default production path
(`kls_current_gap10_t4_r3_timeout120.jsonl`) sharpens the same diagnosis after
the compact row-panel work. The largest ratios are `ASIC_320k` at 3.013,
`ASIC_320ks` at 2.997, `gemat12` at 2.491, `G2_circuit` at 2.477,
`ASIC_100ks` at 2.452, `transient` at 2.337, `rajat28` at 2.233,
`onetone2` at 2.223, `onetone1` at 2.169, and `rajat24` at 2.155. The updated
gap decomposition labels nine of those rows as
`column_egraph_refactor_missing_row_engine`: they use
`last_factor_path=kls_fast_refactor`, spend about 96.5%-100% of measured
EGraph work in the pipeline tail, and report zero row-refactor group work and
zero compact-panel work. `gemat12` is labelled
`klu_first_factor_missing_row_engine` because its loss is dominated by the
prestatic KLU first-factor path. Forcing the current experimental row refactor
does not convert this into a flag-selection issue: `G2_circuit` can build the
row groups but its forced row refactor is slower than the default EGraph path,
and after the singleton-BTF fix `ASIC_100ks` and `onetone2` also build row
groups but remain slower than the default path. This matches
the local CKTSO paper's stated distinction: CKTSO's fast path is a row-major
sparse up-looking numeric engine scheduled from a guessed EGraph, with row-wise
pivot checks and an ETree-descendant pipelined pivoting tail if the guess fails.
KLS's default fast path still updates KLU-compatible column storage through an
exact no-pivot EGraph schedule. The clear paper-backed missing part for these
slow rows is therefore a production row/segment-oriented up-looking numeric
engine, followed by the CKTSO/SubtreeLU pivoting-tail and separator/private-
pipeline machinery; compact row-panel kernels layered onto the current
experimental row mirror are only a precursor.
The row mirror scaffold now also handles singleton-heavy BTF partitions when
building forced row-refactor and row-solve patterns. The previous failure mode
was a KLU-storage boundary issue: singleton BTF blocks have no per-column
`L`/`U` pointer slices in the retained refactor cache, so row-pattern builders
must treat their in-block lengths as zero instead of reading raw numeric
length arrays. With that fixed, forced row refactor builds and runs on
the sampled slow rows. The row task builder now also coalesces long runs of
same-level independent scalar rows into bounded row-batch tasks, preserving
existing U-chain row segments and small-width parallel ready-queue cases. That
reduces forced-row task count on `ASIC_100ks` from 83090 groups to 54297, and
on `onetone2` from 27994 groups to 24409. It is still slower than the current
default EGraph/column path on both sampled rows, and the auto work gate
correctly leaves it off because row-group work exceeds exact EGraph dependency
work. This confirms the row/up-looking scaffold now covers and coarsens the
BTF structures that blocked it before, but it does not change the main paper
diagnosis: the gap needs a production row/segment numeric engine, not a flag
flip.

The fragmented-BTF scale policy improved the large recon artifact geomean over
the preceding EGraph build from 30.58s to 27.37s on the five completed common
rows. A refreshed eight-row selected-large reconstruction at `e610226`, after
the scaled single-block EGraph specialization, scored KLS about 1.28x slower
geomean than CKTSO with 120s timeouts penalized as 1000s, while scoring about
2.88x faster geomean than the saved KLU2 artifact under the same failure
penalty. KLS now wins this selected-large set on `TSOPF_FS_b39_c30` and
`rajat29`, ties the shared `Hamrle3` timeout, and still loses materially on
`pre2`, `nxp1`, `G3_circuit`, `rajat30`, and `ASIC_680k`. KLS still times out
on `pre2` and `Hamrle3`; CKTSO completes `pre2` and times out on `Hamrle3`.
This remains far too large to explain as a missing ordering package alone.

The remaining worthwhile ideas are not per-matrix tuning knobs. They require
new KLS-owned symbolic/numeric machinery:

- Production-scale MC64-equivalent maximum-weight matching with dual
  row/column scaling. KLS now derives dual row/column scaling from its
  KLS-owned exact sparse assignment path for small accepted candidates, builds
  BSD-licensed SPRAL Hungarian/auction matching from a pinned submodule by
  default, or can link to a system SPRAL install, and uses SPRAL for bounded
  pre-factor large weak-diagonal trials plus a value-gated post-factor trial
  for dense high-off-diagonal-pivot cases. This is still not a full production
  MC64-equivalent preprocessing stage. Existing
  MC64-style code can be reused only when its license is LGPL-compatible,
  permits source and binary redistribution with KLS, and allows preservation of
  upstream notices in KLS's third-party notice file. HSL MC64 and
  non-redistributable MC64 copies are out of scope for vendoring.
- A full intra-block parallel factor/refactor scheduler that consumes retained
  EGraph/ETree or separator-tree metadata with pivoting tail restart, beyond
  the current guarded no-pivot EGraph cluster/pipeline refactor.
- Full CKTSO-style fast factorization tail restart after a failed pivot check,
  beyond the current KLS-owned block-local restart and serial
  prefix-current/all-current tail subset.
- A production KLS-owned first-factor engine with row/segment storage, rather
  than the current env-gated KLU-compatible scaffold.
- SubtreeLU-style private/pipeline pivoting factorization from a retained
  separator tree. The current threaded METIS `NodeNDP` call now preserves the
  accepted top-level separator component sequence and private/pipeline row
  split, the experimental row-refactor ready queue can use that map for
  separator-private initial thread queues, and KLS-first row-up factorization
  now first tries the same Algorithm 6-style separator-tree split used by the
  no-pivot row-refactor path: dominant subtrees are collapsed into pipeline
  roots, child subtrees become private-thread candidates, and the candidate
  subtrees are greedily assigned by structural row-input work. The KLS-first
  consumer validates the resulting private ownership against the original
  block row dependencies before remapping; if validation fails it falls back to
  the older retained-component private/pipeline queue. When accepted, it remaps
  the block-local row/column order and `Pnum` to that queue, runs validated
  private rows concurrently, consumes pipeline rows with an Algorithm 3 atomic
  counter, and pre-updates pipeline rows from already-published private
  predecessors and earlier published pipeline-prefix predecessors before the
  ordered publish step. The pipeline row updater now records conservative
  row-supernode
  membership from published row-major `U` patterns and consumes consecutive
  ready predecessor rows, or the already-finished prefix of such a run, through
  one guarded supernode-run executor. That
  executor preserves the triangular in-run discovery rule: applying row `k`
  may create row `k+1` as the next ready dependency before the run continues,
  and it now validates the dense prefix/common-trailing shape so the trailing
  contribution is accumulated once before scatter. This covers Algorithm 4's
  ready-supernode update branch at the scalar compact-kernel level without
  claiming production BLAS-backed panel storage. When a
  safely pre-updated pipeline row needs a scoped dynamic pivot, KLS now
  publishes that pivot row under the ordered pipeline lock, updates both
  committed and phase-local prefix U entries for the column exchange, advances
  a column-order epoch, and makes speculative suffix rows that began under an
  older epoch discard and retry inside the same guarded pipeline phase.
  Separator-covered fast-factor block repairs now also try that Algorithm
  6-style split directly: the repaired block is remapped into private rows plus
  collapsed separator-pipeline roots, private ownership is validated after the
  remap, validated private rows run concurrently, and the separator roots finish
  through the restartable pivot-capable row pipeline. If the private proof is
  unsafe, KLS keeps the repaired block on the ordinary full row pipeline. KLS
  still lacks production coarse BLAS supernode storage and a complete
  checked-tail factor/refactor queue consumer.
- Broader supernodal row/segment updates in the sparse up-looking executor.
  KLS has exact-pattern, ragged single-producer, and default structural and
  work-gated fragmented multi-producer dense-panel updates, but these are
  still narrower than SubtreeLU's production supernodal coverage.

Several smaller dispatch experiments were tried and rejected because they helped
some benchmark cases while regressing others. Those are documented below so the
project does not drift toward benchmark-name-specific heuristics.

## Paper-by-Paper Coverage

| Reference | Implemented in KLS | Partial or open coverage |
| --- | --- | --- |
| Algorithm 907 / KLU | BTF preprocessing, fill-reducing ordering, row scaling modes, Gilbert-Peierls factorization with partial pivoting, no-pivot refactorization, and block back substitution are present through the vendored KLU-derived kernel. KLS adds automatic policy selection, serial refactor scatter metadata, and exact EGraph level metadata around these pieces. | KLS still inherits KLU's fundamentally sequential intra-block numeric kernel. |
| NICSLU | AMD-style ordering, optional static-pivoting preprocessing, optional SPRAL Hungarian/scaling trials, and the idea that parallel kernels should be selected by general structural/numeric evidence are represented in KLS policies. KLS now records exact no-pivot EGraph levels from the numeric U pattern and uses them in guarded large single-block, dominant-BTF-block, and fragmented non-dominant many-block refactor paths, including KLU row-scaled cases where scale handling is supported and work-estimated cluster-level thread slices. KLS also reports the NICSLU R1/R2 static parallel suitability model as `parallel_model_r1`, `parallel_model_r2`, and `parallel_model_recommends_parallel`, using the paper's 2.0 and 50.0 thresholds. After numeric factorization, KLS also evaluates NICSLU Algorithm 4's task-flow earliest-finish model on the actual U-dependency graph using `2*nnz(L(:,i))` update work, `nnz(L(:,k))` normalization work, and a unit dependency sync cost; benchmark artifacts expose the resulting work, finish time, speedup, dependency count, and recommendation. These models seed KLS-owned row/segment refactor preparation when an exact dependency schedule and the row-work gate agree. | Full production MC64 matching/scaling is not implemented. NICSLU's detailed ETree/EScheduler-guided intra-block factorization and full pivoting-aware ETree scheduling are not implemented. Earlier broader EGraph prototypes were rejected because they were not general wins on the current kernel/storage. |
| CKTSO | METIS nested-dissection ordering, guarded SCOTCH nested-dissection auto trials for large high-work symbolic candidates, constrained-minimum-degree-style CAMD refinement, combined ordering selection, pivot-checked fast factorization, default CKTSO-style row-wise guessed-diagonal checks in KLS-owned checked row fast/refactor passes, KLS-owned block-local restart after a failed fast-factor pivot check including root-of-block rejects, conservative serial prefix-current/all-current tail restart for validated non-root unscaled repaired blocks, range-aware row-first repair over retained topological pivoting-tail envelopes including non-contiguous active masks and compact retained ETree-tail row worklists, scaled block-local restart plus threaded checked continuation over later BTF blocks, scaled prefix-current/all-current block repair, scaled in-block serial tail restart after recomputing row scales to input-row order, threaded BTF worker-pool completed-block tracking for safe prefix-current rejects, guarded EGraph cluster/pipeline no-pivot refactors for single, dominant BTF, and selected fragmented many-block BTF shapes, work-balanced cluster-level refactor slices, cached row-permutation solve scratch, CKTSO Section V-style structure-adaptive triangular solve metadata/executor for normal and transpose single-RHS solves, and dual-potential plus optional SPRAL matching-derived equilibration trials are implemented in KLS at the KLU-wrapper layer. | CKTSO's maximum-weight matching with dual scaling is still only partially approximated because KLS does not have an always-on production MC64-equivalent weighted assignment stage. The retained ETree-tail executor exists inside KLS-owned block repair, but the full CKTSO clustered/pipelined fast-factor scheduler still is not implemented across the whole guessed-EGraph interruption path. Otherwise unsupported fast-factor failures still fall back to full pivoting factorization. |
| SubtreeLU | KLS vendors reproducible METIS/GKlib and SCOTCH submodules, uses METIS plus CAMD refinement, asks METIS `NodeNDP` for at least `log2(threads)` nested-dissection levels on larger threaded METIS analyses, and can keep SCOTCH from `auto` when its symbolic score is materially better on large high-work cases. KLS now retains accepted `NodeNDP` component sequences from METIS user-order callbacks, stitches them into a global BTF-aware separator forest with synthetic private components for blocks that did not run `NodeNDP`, reports the resulting global queue shape in stats/bench output, and uses the component map to build separator-private initial thread queues for the experimental row-refactor ready queue. No-pivot and checked row refactors can now consume the retained separator tree through a SubtreeLU Algorithm 6-style FLOP-balanced private/pipeline queue, with separator-crossing row groups forced into pipeline work. The KLS-first pivoting row-up-looking factorization now tries an Algorithm 6-style separator-tree split first: dominant subtrees become pipeline roots, child subtrees remain private candidates, candidates are assigned to private threads by block-local row-input work, and dependency validation falls back to the older retained-component queue if private ownership is unsafe. Accepted queues are remapped into the block-local row/column order and `Pnum`, execute private rows in worker-local entries, and consume pipeline rows through a guarded Algorithm 3 atomic counter. KLS also uses SubtreeLU Algorithm 4's scoped `N'` pivot maximum inside the current collapsed component extent, rejects unsafe cross-domain exchanges, applies ready predecessor row-supernode runs and ready prefixes during partial pipeline-row updates with in-run triangular discovery and a compact common-trailing accumulation, publishes compact-validated ready panels before the first consumer falls back to scalar/compact walking, and handles safely pre-updated weak-pivot pipeline rows with an ordered pivot publish plus column-order epoch retry for speculative suffix rows inside the same pipeline phase. This matches SubtreeLU's separator-domain pivot rule while preserving KLS's fallback to the pivoted block kernel for unsafe scoped rows. KLS also records row-major U-pattern supernode candidate diagnostics from the exact no-pivot refactor dependency pass, can retain the exact row-major U-supernode structural pattern behind `KLS_ENABLE_REFACTOR_U_SUPERNODE_PATTERN=1`, has scalar compact-panel producer and consumer updates, has KLS-owned scalar batched producer-to-consumer-row-subrange updates for exact multi-producer patterns, ragged single-producer suffix patterns in dense and independent row groups, and opt-in structural/work-gated contiguous and fragmented multi-producer row-panel updates for dense producer suffixes. Optional CBLAS experiments remain separate for completed-supernode row updates and unchecked blocked producer-panel `dtrsm`/`dgemm`. | KLS still does not use production SubtreeLU-style coarse supernodes/BLAS updates broadly enough for the paper slow cases, and the native row/segment numeric engine remains a scaffold around KLU-compatible packing. CKTSO's full checked-tail/pivoting executor also remains open. |

This means KLS has implemented or prototyped the ideas that can be layered
around the current KLU-derived data structures. It has **not** implemented all
paper ideas that are still worth trying. The remaining items are general solver
design work, not benchmark-specific tuning.

## Implemented

- KLU Algorithm 907 baseline: BTF preprocessing, AMD/COLAMD ordering,
  KLU-style row scaling modes, Gilbert-Peierls factorization with partial
  pivoting, no-pivot refactorization, block back substitution, diagnostics, and
  pivot tolerance controls are available through the vendored Trilinos
  SuiteSparse-derived sources.
- Nested-dissection ordering option: KLS vendors METIS/GKlib and SCOTCH as
  pinned submodules by default, and can optionally use compatible system METIS
  or SCOTCH installs.
- Combined ordering policy: KLS auto mode can choose AMD, COLAMD, or METIS,
  retry no-BTF symbolic analysis for structural cases where BTF is not useful,
  and promote expensive numeric factorizations to METIS when actual fill/flop
  evidence is better. The no-BTF retry covers both single-block BTF analyses
  and large dominant-block BTF analyses when the stripped fringe is small; it
  also covers many-block analyses whose large block leaves an inflated symbolic
  estimate and whose no-BTF retry cuts the symbolic score at least in half,
  matching KLU's warning that BTF can occasionally increase factor fill.
  Dominant/inflated METIS-start retries now require a known BTF symbolic score;
  otherwise KLS preserves BTF, because an unknown BTF estimate can make a bad
  no-BTF single-block symbolic look artificially preferable.
  Low-work dominant-BTF cases are excluded from this retry because keeping BTF
  is cheaper when the symbolic work estimate is already low. The METIS promotion
  gate also
  covers small BTF-dominant matrices when the first factorization shows many
  off-diagonal pivots and high actual fill/flop growth, and preserves the
  selected BTF/no-BTF mode when comparing a METIS promotion candidate. KLS also
  starts directly with METIS for narrow large-diagonal structural classes from
  the paper corpus: very low-degree full diagonals, sparse full diagonals with
  bounded but meaningful row/column degree, nearly full diagonals with a large
  dense degree spike, and large nearly diagonal sparse-spike or dense-spike
  patterns. The nearly diagonal spike class also starts without BTF so KLS does
  not first pay for an AMD symbolic analysis before settling on the same
  METIS/no-BTF numeric path. Medium spiked low-diagonal patterns with a large
  row/column degree spike also start directly with METIS and max scaling for
  TSOPF-style power-grid structures. The medium bounded-degree METIS start is
  limited to near-full diagonals, because Rommes/BIPS-style rows with a few
  percent missing diagonal entries measured faster on AMD despite fitting the
  same low-degree envelope. For large high-work no-BTF single-block analyses
  outside those direct-start classes, auto can also try METIS
  symbolically before the first numeric factorization and keep it when the
  symbolic fill score is clearly lower, avoiding a delayed post-factor METIS
  promotion. This preserves the papers' nested-dissection motivation without
  naming individual matrices.
- Constrained nested-dissection refinement: KLS can refine METIS rank groups
  with CAMD constraints, preserving nested-dissection rank shape while reducing
  local fill and flops. Large METIS orderings now use coarse rank-group CAMD
  refinement by default instead of limiting the CKTSO-style refinement to one
  medium structural class.
  Threaded METIS analyses for matrices with at least 30,000 rows now use
  `METIS_NodeNDP`, which enforces the top nested-dissection levels needed for a
  thread-count-sized separator-tree skeleton before the same CAMD refinement.
  Accepted METIS analyses also retain the `NodeNDP` top-level component sequence
  as private leaf domains plus pipeline separator components, including an
  ordering-position to component map. The row-refactor ready queue can consume
  that map to assign its first ready groups by separator component before the
  shared/local successor queue takes over.
- Scaling policy: KLS exposes KLU scale modes and has auto scale selection based
  on pattern and numeric evidence, including no-scale reuse where repeated
  SPICE refactorization benefits. Low-work dominant-BTF cases also start with
  no KLU row scaling, avoiding the scale recomputation cost while preserving
  the useful BTF decomposition. Large high-work METIS/no-BTF single-block paths
  that already selected max scaling skip redundant post-factor scale trial
  factorizations. Medium and large TSOPF-style spiked low-diagonal METIS starts
  begin with sum scaling and a lower `1e-4` pivot tolerance once their
  structural dominant-BTF shape is known, avoiding the earlier
  max-scale/default-tolerance discovery factorizations. Small spiked
  low-diagonal METIS starts use unscaled `0` mode with the same lower pivot
  tolerance and skip static-pivot trials that do not improve that structural
  class.
- Fast repeated factorization with pivot check: KLS reuses an existing numeric
  pattern, checks reused pivots against the selected threshold via L
  multipliers, and falls back to full pivoting factorization if the reused order
  is unsafe and cannot be repaired locally. For unscaled BTF patterns, KLS can
  restart the rejected diagonal block with pivoting, rebuild the global row
  permutation/off-diagonal entries around that block repair, and retry the
  checked fast factorization. KLS records the first rejected factor-order pivot,
  original matrix column, rejected block suffix, exact U-pattern descendant
  tail, ordered-block ETree successor path, sorted pivoting-tail worklist scope,
  and number of block restarts so future tail-restart work can distinguish
  late-tail failures from early failures and compare no-pivot versus
  pivoting-tail recomputation scopes.
- Static pivoting trial: KLS has value-aware greedy row matching, layered
  augmenting-path search for larger weak-diagonal candidates, and swap
  improvement for weak or high-off-diagonal-pivot medium matrices. Medium
  static-pivot candidates also run a bounded alternating-cycle pass that can
  apply profitable three- and four-row exchanges missed by pair swaps, while
  larger candidates keep the lower-fill layered/swap pattern. Small exact
  assignment candidates now also return assignment-dual row/column scaling, and
  KLS can still trial dual-potential matching-derived row/column equilibration
  for non-exact matches. It keeps the
  transformed candidate only when numeric quality and cost evidence justify it.
  With default SPRAL scaling enabled, large weak-diagonal candidates whose BTF
  analysis leaves one dominant block can run BSD-licensed SPRAL Hungarian
  matching/scaling before the first factorization; many-block BTF cases are
  excluded because the existing BTF path is already efficient there. After
  cheaper scale, ordering, and pivot-tolerance trials have run, KLS can also
  run a SPRAL Hungarian/scaling trial for dense high-off-diagonal-pivot cases.
  The accepted candidate must pass the normal numeric comparison plus a
  matching-specific value gate: it must remove substantial pivoting pressure or
  materially reduce factor work/fill, which rejects small cases where matching
  setup costs more than the repeated-refactor saving.
  The pre-factor static-pivoting gate also covers moderately sized matrices
  whose input values show a majority of weak or missing diagonal entries, plus
  medium-large mostly diagonal matrices with thousands of weak diagonal rows.
  These are general numeric-structure rules used by the frequency-domain and
  Rajat-family paper cases without naming individual benchmarks.
- Initial KLS-owned parallelism: repeated refactorization can run across
  independent BTF diagonal blocks for large high-flop cases where the block work
  is wide enough to offset thread overhead. The threaded path keeps a persistent
  worker pool on the solver instance so repeated SPICE refactors reuse workers
  and scratch storage instead of relaunching threads each cycle. Checked worker
  runs mark completed diagonal blocks and classify a pivot reject as
  prefix-current only when all earlier blocks completed.
- KLS-owned serial refactor metadata: for unscaled numeric patterns that are
  not handled by the threaded BTF worker pool, KLS precomputes the fixed
  scatter from factor-order columns to pivotal rows and input value positions.
  This removes repeated `Q`/`Pinv` structure lookups from SPICE refactor cycles
  while preserving the existing KLU-derived LU storage. It covers both
  single-block and serial BTF refactors. For serial BTF refactors, KLS also
  prepartitions each mapped column into off-block and diagonal-block entries, so
  repeated refactors do not reclassify the same BTF structure in the numeric
  scatter loop. A narrow scaled many-fringe dominant-BTF subset uses the same
  serial map after recomputing row scales and then permuting the scale vector
  back to pivot order.
- KLS-owned EGraph metadata and guarded refactor consumption: for threaded
  numeric runs, KLS can levelize the exact no-pivot refactor dependency graph
  from the actual U pattern after factorization. It retains level pointers and
  column lists and reports the level count, maximum level width, and dependency
  edge count in `kls_stats` and benchmark JSON. It also reports the
  CKTSO-style cluster/pipeline split point, pipeline-column count, approximate
  no-pivot update work, and tail work implied by the current level widths.
  Large single-block, high-flop dominant-BTF, and very fragmented unscaled
  many-block matrices with enough dependency work and level width can now
  consume this schedule through a work-estimated per-thread cluster-mode refactor,
  then switch to a no-pivot pipeline tail where each worker claims tail columns
  from an atomic cursor and waits only for actual U-pattern predecessors. The
  same block-aware EGraph kernel can run inside a large dominant BTF block and
  update Offx for entries above that block. The EGraph consumer now keeps a
  solver-owned worker pool and dense worker scratch across repeated refactors,
  plus generation-stamped solver-owned pipeline dependency markers, matching
  the CKTSO/NICSLU emphasis on retained scheduling state instead of relaunching
  threads or clearing a fresh done array on every SPICE step. Schedule
  construction is now
  limited to single-block, dominant-BTF, or high-work fragmented non-dominant
  many-block shapes with enough numeric or measured dependency work to consume
  it, while ordinary non-dominant many-block BTF and low-work dominant-BTF cases
  without enough dependency work skip the setup and stay on the BTF worker pool
  or mapped refactor paths.
  For the high-coverage many-fringe dominant-BTF class below the normal EGraph
  size floor, KLS can run the whole exact EGraph as an atomic topological
  pipeline with no cluster barriers. This is a retained SubtreeLU/CKTSO-aligned
  scheduler improvement for the current fixed-pivot LU storage.
  Very high-work single-block row-scaled factors can also use that full
  all-pipeline EGraph path once the existing scale recomputation/permutation
  support is available, so the scheduler does not leave a tiny barriered tail
  after hundreds of narrow cluster levels. Those scaled single-block EGraph
  refactors now use a dedicated column kernel that applies KLU row scales
  directly while loading the fixed input-position map, instead of paying the
  generic BTF-capable value-loader branch on every entry.
  A later fragmented many-block gate applies the same CKTSO-style lesson to
  unscaled BTF decompositions whose largest block covers less than half the
  matrix but still carries enough numeric work to justify intra-block
  scheduling. On the 93-matrix medium paper corpus this changed only
  `ASIC_680ks` schedule counters, cutting its cycle from 17.62s to 9.42s and
  moving KLS geomean from 0.28594s to 0.28392s with the same three known
  failures.
  A later retained low-work dominant-BTF gate lowers the EGraph schedule floor
  when measured dependency work is still material. On the 93-matrix medium
  paper corpus this moved KLS geomean from 0.32749s to 0.32198s with the same
  three known failures, mainly by cutting IBM `dc1/dc2/dc3/trans4/trans5`
  cycles by about 22-30% and `scircuit` by about 20%.
  This is still narrower than CKTSO's production pivoting machinery, but it is
  the first retained intra-block EGraph cluster/pipeline refactor path.
- SPICE-cycle orientation policy: KLS can analyze normal and transposed storage
  orientations and select the faster internal form for repeated solve cycles.
- LGPL project licensing and third-party notices: KLS itself is
  LGPL-2.1-or-later, with vendored-source attribution separated from the KLS
  license. MC64-equivalent preprocessing must stay inside that licensing
  boundary: HSL MC64 and solver-tree copies that retain HSL redistribution
  restrictions are not vendorable, while KLS can build the BSD-licensed SPRAL
  scaling subset from `third_party/spral`, use another LGPL-compatible
  redistributable MC64-style source with preserved notices, or use independent
  KLS code. The Rust `rwl/mc64` package is also BSD-licensed and useful as a
  reference, but it is a partial SPRAL translation and is not a better fit than
  the pinned SPRAL submodule for KLS's C/Fortran build.
- Paper-derived benchmark manifests: the full public SuiteSparse union from
  the local KLU, NICSLU, SubtreeLU, CKTSO papers and CKTSO ordering supplement
  resolves to 110 matrices. The routine medium subset contains 93 matrices, and
  the large supplement records the 17 excluded public paper cases for deliberate
  overnight tuning.

## Partially Implemented

- CKTSO-style static pivoting is only partial. KLS has a practical weighted row
  permutation, dual-potential matching-derived row/column equilibration, and
  optional SPRAL Hungarian/scaling trials, but not an always-on production
  MC64-style maximum-product matching algorithm with assignment dual scaling.
- NICSLU/CKTSO parallel scheduling is present at BTF-block granularity and, for
  large single-block cases, as a guarded exact-EGraph cluster/pipeline
  no-pivot refactor. This path now covers both unscaled factors and KLU
  row-scaled factors whose scale vector is recomputed before the EGraph refactor
  and permuted afterward. Failed checked passes now seed the retained pivoting
  tail from the interrupted guessed-EGraph unfinished set before prefix refresh,
  and KLS-owned block repair can consume that retained ETree tail as a compact
  topological row worklist. KLS still does not have CKTSO's full
  guessed-EGraph interruption scheduler around that executor.
- CKTSO fast factorization is present as pivot-checked reuse plus a KLS-owned
  BTF-block repair path. KLS validates the preserved-prefix live state and can
  enter a conservative serial pivoting-tail kernel for non-root rejects in that
  subset; threaded block repairs can consume retained ETree-tail rows through a
  compact topological row-pipeline worklist; root rejects can use the same
  KLS-owned pivoted kernel as full block restarts. KLS still does not implement
  the complete CKTSO fast-factor scheduler around that tail executor.
- SubtreeLU-style nested-dissection metadata is now retained from accepted
  METIS `NodeNDP` analyses as private/pipeline component queues. The
  experimental no-pivot row-refactor ready queue, checked row fast/refactor
  path, KLS-first pivoting row-up factorization, and separator-covered
  fast-factor block repair can all consume that map through validated Algorithm
  6-style private/pipeline splitters. The remaining separator-side gap is not
  queue retention anymore; it is using comparable scheduling inside CKTSO's
  complete ETree-descendant pivoting-tail executor and broader production
  row/supernode numeric storage.

## Not Implemented Yet

- Full production MC64-equivalent weighted matching and assignment-dual
  row/column scaling across the broad matrix set.
- NICSLU's ETree/EScheduler executor beyond the implemented R1/R2 and
  Algorithm 4 task-flow suitability counters.
- Intra-block parallel factorization with pivoting scheduled by an ETree.
- General ETree-based pivoting tail restart beyond the current guarded no-pivot
  EGraph cluster/pipeline path, KLS-owned retained-tail block repair, and serial
  prefix-current block-tail subset.
- CKTSO dual-mode cluster/pipeline fast factorization with pivot check.
- CKTSO's complete guessed-EGraph interruption scheduler around pipelined
  ETree-descendant tail factorization with pivoting after a pivot-check failure.
- SubtreeLU separator-tree collapse for checked-tail and refactor pivoting
  kernels outside the KLS-first row-up and separator-covered fast-factor
  block-repair paths.
- Broader/default SubtreeLU FLOP-balanced separator-tree partitioning for
  refactorization across BTF forests and checked-tail kernels.
- SubtreeLU constrained pivot search within nested-dissection subdomains.
- Production SubtreeLU supernodal updates with coarse BLAS kernels inside the
  sparse up-looking framework.

## Tried During This Audit

An intra-block levelized no-pivot refactor prototype was tested using the
existing U structure as the exact dependency graph, matching the NICSLU/CKTSO
EGraph refactorization idea. Broad versions were removed before commit because
they regressed dominant-block circuit cases where the current BTF-block
threaded refactor is faster. A later version replaced per-column mutex
scheduling with static per-thread level slices, closer to CKTSO cluster mode,
but still regressed the intended single-block case (`rajat15`) by about 16% on
the focused repeated-refactor sample. Keeping those broad dispatch rules would
have required case-specific tuning, which is not the desired direction for KLS.
After the EGraph path was narrowed and retained, the pipeline tail scheduler
was changed from fixed per-thread strides to an atomic dynamic work cursor,
matching CKTSO's on-the-fly tail assignment more closely while preserving KLS's
existing fixed-pivot column kernel. On the nine-row hard focus set, same-session
4-thread runs moved from about 9.07s geomean for the saved KLS artifact to
about 8.91-9.04s. A 90-common-row paper-medium run kept the same failure set as
the recent KLS medium baseline and improved geomean on common rows versus that
baseline, while KLS still trailed the CKTSO artifact by about 1.29x geomean on
the 90 common medium rows.

An EGraph-specific persistent worker-pool prototype was then tested to avoid
recreating pthreads and scratch arrays on every no-pivot EGraph refactor. The
idea matched CKTSO's emphasis on retained scheduling machinery, but
same-session A/B checks against commit `f36bd4f` did not show a general win:
`ASIC_680k` moved only from about 0.06887s to 0.06784s average refactor time,
while `rajat30` moved from about 0.30470s to 0.30623s and `nxp1` from about
0.64476s to 0.65076s. The prototype was removed because the current EGraph
runtime is dominated by numeric scatter/update work rather than pthread launch
overhead on the large CKTSO-gap rows.

Two later EGraph scheduler tweaks were also rejected on the same basis. A
dynamic work-claim path inside wide cluster levels was neutral on `ASIC_100ks`
but slower on `rajat24` and `transient`. Forcing clustered cases into a full
all-pipeline EGraph schedule was also slower on `ASIC_100ks`, `rajat24`, and
`transient`. These tests indicate that KLS is not mainly missing another
cluster/pipeline dispatch tweak in the current column-storage EGraph kernel.
The CKTSO and SubtreeLU papers point instead to row-major sparse up-looking
storage and row/supernode updates as the remaining large lever.

A broader KLS-owned serial no-pivot refactor path was also prototyped by
reusing the threaded BTF-block refactor kernel when thread-level parallelism was
not eligible. It passed correctness tests and helped some no-scale cases, but it
regressed other representative circuit cases such as `bcircuit` and `rajat03`.
The broad prototype was removed because it did not represent a general
improvement. A later narrower version kept only the precomputed unscaled serial
scatter metadata; same-session suite testing showed that subset as a modest
general win.

The refactor scatter map was also tried in scaled and threaded-worker forms.
The scaled extension improved a focused hard case but regressed the full suite,
and passing the map through the worker pool regressed representative threaded
BTF cases. KLS therefore keeps the map on the unscaled serial path and leaves
the threaded worker-pool scatter path separate.

The unscaled serial BTF scatter map was extended with a per-column boundary
between off-block entries and diagonal-block entries. This is a retained
general metadata improvement: a 25-matrix, 5-pass same-session extended-suite
A/B run improved the geomean from 0.03665s to 0.03651s and the median ratio to
0.993, with the largest win on `ckt11752_tr_0`. The single-block map path was
left in its original one-pass layout because single-block refactors have no
off-block entries and the partitioning setup did not help them.

A second layer of persistent refactor metadata was prototyped by precomputing
per-column L/U value and index pointers for mapped serial refactors. Focused
tests helped some repeated-refactor losses, but full-suite same-session testing
regressed the geometric mean because low-arithmetic and memory-sensitive cases
lost more than the hard cases gained. A structural gate based on flop count,
BTF shape, orientation, and arithmetic density reduced the regression but still
did not beat the committed baseline.

A narrower scaled single-block scatter map was also retried with a structural
gate for high off-diagonal-pivot cases. It repeatedly improved the targeted
scaled single-block case, but the added scaled path still perturbed unscaled
single-block hot cases enough to fail the no-regression bar. That prototype was
removed as well; scaled refactors still use the KLU-derived path unless they are
handled by the existing threaded BTF worker path.

A work-balanced BTF refactor scheduler was tested by sorting independent BTF
blocks by an LU-length work estimate before launching worker threads. It was
removed because the extra scheduling work did not improve the dominant-block
cases where BTF-block threading is eligible, and it regressed representative
threaded refactor timings on `ckt11752_tr_0` and `circuit_4`.

The reactive static-pivoting gate was also widened from medium matrices to
larger sparse matrices with high off-diagonal pivot counts. The existing
acceptance checks rejected `rajat22`, but only after paying a large matching
and refactorization trial cost; `rajat27` accepted a better numeric pattern but
lost on the repeated-SPICE metric because the one-time trial cost dominated.
The gate was restored to avoid converting matching into a broad overhead.

An auction-style weighted assignment pass was prototyped to move the current
greedy row matching closer to MC64's maximum-weight matching. It was removed
because the only current extended-suite matrices that select static pivoting
(`gemat11` and `gemat12`) already have good enough matched patterns; the extra
auction and dual-scaling work increased first-factor cost without improving
fill, refactor time, or residuals. A cheap weighted-gap guard avoided the
largest regression, but the guarded implementation still did not improve the
suite enough to justify the additional code.

The BTF-block threaded refactor scheduler was also changed from a mutex-protected
block counter to a C11 atomic work counter. This reduced scheduler overhead on
some BTF-threaded samples such as `ckt11752_tr_0`, but same-session median
suite testing showed no aggregate improvement and small regressions on other
cases. The mutex scheduler was kept until a broader scheduling change has a
clearer win.

The BTF-block parallel eligibility gate was relaxed to try more medium and
dominant-block structures. This was a general structural dispatch experiment,
not a matrix-name rule, but the added thread scheduling overhead regressed the
focused repeated-refactor samples that motivated the test. The conservative
high-flop, many-block eligibility gate was kept.

The no-pivot refactor update loops were also split into pivot-checking and
non-checking variants to reduce branch work in the fast repeated-factorization
path. The change passed correctness tests, but focused A/B timings were mixed
and did not show a general win; the simpler shared loop remains in place until a
larger KLS-owned numeric kernel makes this separation worthwhile.

The no-BTF symbolic retry was initially widened from single-block BTF cases to
dominant-block BTF cases, including METIS-started auto orderings. This is
consistent with KLU's observation that BTF can occasionally hurt, but the broad
version paid extra symbolic-analysis cost on unrelated multi-block matrices and
regressed the full 25-matrix same-session extended suite, with the candidate
geomean at 0.03816s versus the baseline at 0.03799s. A later narrower version
was retained: it only retries no-BTF for large BTF analyses whose largest block
covers at least 95% of the matrix with a small stripped fringe, and it only
accepts the no-BTF symbolic when the score improves by at least 20%. That
rescued the KLU-paper `Raj1` and `rajat24` cases from timeout/missing status
under the paper-medium run while keeping ASIC-style dominant-block cases on BTF.
The METIS promotion trial was also fixed to preserve the selected BTF/no-BTF
mode, which lets `Raj1` promote from no-BTF AMD to no-BTF METIS.

The retry was later extended only for many-block BTF analyses with a large but
not overwhelming dominant block. This targets cases like the IPSO HTC matrices,
where BTF created about 29k blocks but still left an 87% dominant block and a
much larger symbolic estimate than the no-BTF structure. The retained gate
requires at least 1024 BTF blocks, at least 100k rows, a largest block covering
80-95% of the matrix, no low symbolic-work estimate, and a no-BTF symbolic score
at most half of the BTF score. The METIS-started auto path now runs this strict
retry before returning; it still skips the older single-block retry to avoid
extra symbolic work on already-good METIS/BTF single-block cases. In same-session
checks, `HTC_336_4438` moved from the saved METIS/BTF scale-1 path
(`initial=13.14s`, `refactor=0.082s`) to METIS/no-BTF no-scale
(`initial=3.52s`, `refactor=0.103s`), while `HTC_336_9129` moved from
`initial=8.37s` to `1.93s`. Guard cases kept their earlier BTF decisions:
`G2_circuit` remained METIS/BTF, `transient` and `power197k` remained AMD/BTF,
and `ASIC_680ks` remained METIS/BTF. A parallel SCOTCH sweep did not justify an
auto SCOTCH policy: it was much slower than METIS on `rajat30` and `nxp1`, and
mixed or worse on `G2_circuit`, `transient`, and `power197k`.

The direct METIS retry was then tightened after guard checks showed a harmful
interaction with the ASIC 320k family. Those matrices have a very useful BTF
decomposition, but their METIS/BTF symbolic score is unknown while the no-BTF
symbolic has a huge explicit fill estimate. Accepting that no-BTF candidate sent
`ASIC_320k` to a 27M+27M single-block symbolic and caused a timeout, while the
retained METIS/BTF path has about 2.0M+2.0M numeric fill and refactors near
0.10s. KLS now requires a known current BTF score before accepting dominant or
inflated many-block no-BTF retries, and the broad large-low-degree structural
shortcut no longer starts no-BTF blindly. Clean checks restored `ASIC_320k` and
`ASIC_320ks` to METIS/BTF, while `HTC_336_4438` still selects METIS/no-BTF and
`power197k` keeps AMD/BTF.

The pre-factor static row-matching path was also tested without the augmenting
and swap-improvement pass, leaving only the initial greedy maximum-value
matching. This reduced setup cost on some static-pivot samples, but after
rebuilding a clean baseline the 25-matrix extended suite regressed from
0.03798s to 0.03869s. The full improvement pass remains enabled for pre-static
matching.

After adding the paper-derived benchmark corpus, the pre-factor static
row-matching gate was widened for moderately sized matrices with at least half
of rows having weak or missing diagonal entries. This retained the full
augmenting/swap improvement pass and is a general numeric-structure rule, not a
benchmark-name rule. On the AT&T `onetone1`/`onetone2` paper cases, the
focused SPICE-cycle geomean improved from about 99.5s to 9.3s and both cases
selected the static row permutation. The retained gate was later widened to
cover the AT&T `twotone` scale as well; `twotone` now completes the paper-medium
suite run and selects static pivoting, reducing off-diagonal pivots from about
9500 to about 1600. A still-broader attempt to cover `mac_econ_fwd500` scale
was not retained because it timed out with multi-GB memory use.

The same pre-factor static row-matching gate was then extended to
medium-large mostly diagonal matrices with thousands of weak diagonal rows and
very few missing diagonals. This retained a structural diagonal-completeness
guard, so earlier Rajat-family cases with fewer weak diagonals remain on the
normal dynamic-pivot path. On the paper-medium corpus, the retained gate
improved the KLS geomean from about 1.57s to 1.51s and reduced the affected
Rajat cases' off-diagonal pivots to 1-3. KLS is still much slower than CKTSO on
these cases, so this is a partial static-pivoting improvement rather than a
replacement for full MC64-style matching/scaling or a faster KLS-owned numeric
kernel.

The large-matrix static row-matching augment was then changed from one
independent breadth-first search per unmatched row to a Hopcroft-Karp-style
layered augment. This is a retained MC64-adjacent improvement to the existing
greedy matcher, not full weighted MC64 dual scaling. It improved the AT&T
`twotone` paper case in a same-session focused run from about 273s to about
124s on the 100-step SPICE-cycle estimate by finding a lower-fill static row
permutation, while leaving `onetone1`, `onetone2`, and the large Rajat static
cases essentially neutral. The 25-matrix extended suite still passed, with a
geomean of about 0.04409s versus about 0.04421s for the previous mainline run.
A prototype greedy-plus-layered perfect matching for `mac_econ_fwd500` reduced
matched-METIS refactor time from about 22s to about 16s, but that remained far
behind the CKTSO artifact at about 2.5s, so the `mac_econ_fwd500` gate was not
widened.

For medium-large static-row-matched matrices whose diagonal is both weak and
mostly missing, the matching-equilibration trial was narrowed to keep the row
permutation but prefer no numeric scaling. This is a structural rule for
80k-150k order, at most 1.5M nonzeros, at least half missing diagonal entries,
and at least half weak-or-missing diagonal rows. On AT&T `twotone`, it selected
`scale=-1`, reduced off-diagonal pivots from about 1514 to 219, skipped the
matching-equilibration setup, and improved the focused SPICE-cycle estimate
from about 124s to about 96-98s. It deliberately does not apply to the Rajat
static-pivot cases, whose diagonals are almost complete and were slower when
forced unscaled. This remains far behind the CKTSO artifact at about 26s on
the same SPICE-cycle formula.

After expanding the paper-medium corpus with the public CKTSO ordering
supplement cases, the post-factor METIS promotion gate was widened for small
BTF-dominant matrices whose initial AMD/COLAMD factorization produces both
many off-diagonal pivots and high fill/flops. This is a retained combined
ordering rule, not a matrix-name rule: the gate requires at most four BTF
blocks, a dominant block covering at least 90% of the matrix, at least 128
off-diagonal pivots, and enough actual factor work to amortize the METIS trial.
It selected METIS on the TSOPF/QY power-grid-style cases in the expanded paper
medium suite, improving the 87-common-row KLS geomean from about 0.333s to
0.322s and the common KLS/CKTSO ratio from about 1.48x to 1.43x. Guard cases
with no off-diagonal pivot pressure, BTF disabled, or many BTF blocks stayed on
the existing AMD path.

A medium spiked low-diagonal structural METIS start was retained for the large
TSOPF-style paper cases. The rule requires a 50k-125k order matrix, about
20-32 entries per column on average, a 20-35% diagonal fraction, no empty rows,
and one large row/column degree spike covering about 40-60% of the matrix. It
initially started with METIS while keeping BTF enabled, and left auto scaling
on the default max-scaling path rather than forcing no-scale. In same-session checks,
`TSOPF_FS_b39_c19` improved from the earlier about-568s SPICE-cycle estimate to
about 546s, and the widened bound rescued the large-supplement
`TSOPF_FS_b39_c30` from the 120s per-process timeout with a SPICE-cycle
estimate of about 920s. This is still much slower than CKTSO on `c30`, whose
same-session artifact is about 333s, but it turns a KLS timeout into a
completed structural paper case. Smaller TSOPF/QY guard cases retained their
previous policies.

Two additional low-degree structural METIS starts were retained to remove the
last current-auto 120s paper-medium timeouts. The very-low-degree full-diagonal
METIS start now allows row/column degree up to 8, covering `ss1` while still
leaving broader low-degree circuit families on their existing policies. A
separate sparse-diagonal low-degree class covers 150k-250k order matrices with
5-8 entries per column on average, 5-20% diagonal coverage, no empty rows or
columns, and row/column degree at most 64; it starts with METIS and sum scaling.
Both classes skip redundant post-factor auto scale trials. In same-session
120s-capped runs, `ss1` completed in about 104s wall time with a SPICE-cycle
estimate of about 3413s, and `mac_econ_fwd500` completed in about 107s wall
time with a SPICE-cycle estimate of about 2252s. This is a timeout/completeness
improvement, not a CKTSO win: the existing CKTSO artifacts are about 349s and
265s respectively on the same SPICE-cycle formula. The same scale-trial gate
also reduced `G2_circuit` initial factor time from about 2.8s to about 1.0s and
modestly improved `mc2depi`.

A selected large-supplement reconnaissance manifest was added for the smaller
large paper cases that are practical under a 120s per-process cap before an
overnight full-large run. In same-session KLS/CKTSO runs with one factor and
one refactor repeat, KLS beat CKTSO on `rajat29` by about 5%, but lost the
other common successful cases: about 3.1x on `ASIC_680k`, 4.0x on `rajat30`,
4.5x on `G3_circuit`, and 4.7x on `nxp1`. KLS also timed out on
`TSOPF_FS_b39_c30` and `pre2`, both of which CKTSO completed within the same
120s cap; both solvers timed out on `Hamrle3`. The later widened TSOPF
spiked-low-diagonal rule removes the `c30` timeout, but `pre2` remains an
unresolved large-case timeout. The result reinforces that the remaining
large-case gap is mostly very large single-block or near-single-block
numeric/refactor throughput, not the many-small-BTF-block class.

The remaining `pre2` selected-large timeout was then checked as a simple policy
question before attempting new numeric-kernel work. Eight 120s-capped variants
all timed out with one factor and one refactor repeat: AMD/max-scale,
AMD/no-scale, COLAMD/max-scale, METIS with max/sum/no-scale, AMD/no-BTF
max-scale, and METIS/no-BTF max-scale. Focused sweeps on the completed but
slow `nxp1` and `rajat30` cases also found no better simple dispatch: `nxp1`
needs the current auto no-BTF METIS path, while AMD/COLAMD time out and
BTF-enabled METIS is much worse; `rajat30` remains best under current auto
METIS/no-BTF max-scaling. These results make `pre2`, `nxp1`, and `rajat30`
poor candidates for another ordering/scale/BTF heuristic. They need the open
paper ideas around faster single-block numeric/refactor kernels, matching
quality, or EGraph/separator-tree scheduling.

SCOTCH was then added as a pinned, reproducible optional ordering package and
tested as an explicit `--ordering scotch` path. SCOTCH symbolic analysis
completed on `pre2`, but numeric factorization still timed out at 120s; a
`rajat30` SCOTCH numeric run was also not competitive before interruption near
the same cap. The result did not justify broad SCOTCH auto ordering or a new
matrix-specific rule. Later auto support therefore kept SCOTCH restricted to a
guarded large single-block symbolic trial that is accepted only on a clear
fill/work score win.
ParMETIS was deliberately not added in this pass because it is an
MPI/distributed-memory package, while the current KLS benchmark and solver path
is shared-memory and single-process.

The CKTSO paper's nested-dissection-plus-constrained-minimum-degree idea was
then strengthened for KLS METIS orderings by applying coarse rank-group CAMD
refinement to large METIS orderings, not only one medium structural class.
This helped the completed large cases: same-session `--ordering metis --no-btf`
runs finished `nxp1` with about 1.16s factor and 1.13s refactor averages, and
`rajat30` with about 0.70s factor and 0.70s refactor averages. Auto selected
the same METIS/no-BTF path for those cases. The change did not fix `pre2`:
both explicit METIS/no-BTF and auto still timed out at 120s. This reinforces
that `pre2` is a missing-major-algorithm case rather than a separator-package
case.

A follow-up coarse-grouping probe reduced the large METIS CAMD group size from
4096 to 1024 as a general CKTSO-style constrained-ordering experiment. It was
rejected: `pre2` no-BTF METIS symbolic fill worsened from about 121.5M to
124.7M nonzeros, so the retained 4096 grouping remains the better large
default.

The matching-derived equilibration pass was then moved closer to the
MC64/NICSLU/CKTSO preprocessing contract by first solving dual-potential
scaling constraints for the current greedy row match. When the constraints are
consistent, the matched diagonal is normalized to one and nonmatched entries
are bounded by the matched diagonal; when the greedy match leaves large
positive-cycle evidence, KLS falls back to the older heuristic balancing pass.
Same-machine checks on the retained static-pivot representatives (`gemat11`,
`gemat12`, `twotone`, and `rajat25`) produced the same fill, flops,
off-diagonal pivot counts, conditioning, and residuals as the previous
mainline. A temporary large-gate experiment also let `pre2` try this
preprocessing path with layered matching, but the factor-only run still timed
out at 120s. This makes the dual-potential pass a bounded MC64-adjacent
preprocessing cleanup, not the missing CKTSO-scale algorithm.

An LGPL-compatible medium-matrix SPRAL Hungarian-first static-pivot experiment
was also tested on the hard paper focus set. The structural gate targeted
50k-plus mostly complete diagonals with thousands of weak diagonal entries, so
majority-missing AT&T-style cases stayed on the retained greedy unscaled path.
It made `rajat28` select exact SPRAL matching, switch from KLU max scaling to
matching-derived no-scale values, and reduce off-diagonal pivots from 1 to 0,
but its refactor time worsened from about 0.162s to about 0.169s and the
9-row KLS/CKTSO focus ratio regressed from about 2.66x to about 2.72x. The
experiment was removed. This reinforces that BSD/LGPL-compatible MC64-style
matching is allowed and present through SPRAL, but widening it over medium
mostly diagonal cases is not enough to close the large CKTSO gap.

A bounded alternating-cycle improvement was then added after the existing
layered cardinality augment and pair-swap pass. This is not a full MC64
shortest-augmenting-path implementation: it only accepts profitable local
cycles of up to four rows, and it is intentionally limited to `n <= 50000`.
Same-session A/B against the previous commit showed useful medium effects:
`gemat12` reduced off-diagonal pivots from 11 to 7 and improved reciprocal
condition evidence, and `onetone2` reduced actual fill/flops. The same ungated
pass was rejected for larger static-pivot cases because `twotone` and
`rajat25` increased fill/flops despite comparable or better pivot counts, so
large static-pivot candidates stay on the lower-fill layered/swap pattern. A
large `pre2` analyze-only check still shows a dominant 629628-row block and
about `2.08e11` estimated flops under the current AMD/BTF symbolic path, so the
remaining CKTSO gap is still a missing numeric-kernel/scheduler issue.

A cheap exact sparse assignment stage was then added ahead of the greedy
static row matcher for small candidates. It solves a min-cost augmenting-path
problem on `col_max - log(abs(a_ij))`, which is equivalent to maximum-product
matching when a full assignment is found. The path is intentionally gated to
`n <= 4000`, `nnz <= 75000`, and `n * nnz <= 5e7`; a same-session HB/gemat
experiment with a wider gate improved `gemat12` off-diagonal pivots from 7 to
3, but raised initial factor/preprocessing time from about 0.034s to about
5.2s. With the retained gate, `gemat11` and `gemat12` stay on the previous
fast matcher, while synthetic 3000-row static-pivot smoke tests exercise exact
matching. This closes a small piece of MC64 functionality but confirms that
KLS still needs an optimized production MC64-equivalent matcher, not the
straightforward min-cost implementation, before the CKTSO-scale gap can close.

The exact assignment implementation was then changed from a generic
source/sink residual graph to a KLS-owned sparse row/column
shortest-augmenting-path matcher. This removes the extra source/sink edges and
keeps the code LGPL-compatible, but the policy gate remains conservative. A
same-machine wider-gate retest with the corrected source/sink semantics
selected exact matching for `gemat11` and `gemat12`. `gemat12` improved from 7
to 3 off-diagonal pivots, but initial factor/preprocessing time was still about
4.24s; `gemat11` stayed at 0 off-diagonal pivots and also paid about 4.24s.
Restoring the small gate returned `gemat11`/`gemat12` to the prior fast path
(`selected_exact_matching=false`, 0 and 7 off-diagonal pivots), while
`onetone2`, `twotone`, and `rajat25` also stayed off the exact path. The
lesson is that maximum-product matching alone is too expensive to use as a
medium/large default; KLS still needs MC64-quality scaling and acceptance plus
a more optimized assignment implementation before this can close the CKTSO
gap.

The exact sparse assignment path now keeps the assignment dual potentials when
the full match succeeds and converts them into row/column scaling for the
accepted candidate. For each original row matched to column `j`, KLS stores the
row scale at the permuted row `j` and the column scale at original column `j`;
the resulting scaled matched diagonal has unit magnitude, while all other
entries satisfy the dual reduced-cost bound. The smoke suite now requires the
3000-row scaled static-pivot fixture to report
`selected_exact_matching_scaling`, and benchmark JSON exports the same flag.
This fills the small exact-matching part of the NICSLU/CKTSO MC64 contract
without changing the conclusion above: larger paper cases still need an
optimized production MC64-equivalent implementation or a retained
LGPL-compatible SPRAL path, not a wider use of the straightforward exact
augmenting-path code.

An optional SPRAL hook was then added for BSD-licensed matching/scaling
support. It is deliberately not a default dependency and not a solver
replacement.
KLS first wired SPRAL auction matching as a large structural-deficit fallback,
then added SPRAL's MC64-like Hungarian unsymmetric matcher/scaler. The hook was
initially system-SPRAL only; KLS now also pins upstream SPRAL as a submodule and
builds just its scaling subset for reproducible LGPL-compatible MC64-adjacent
experiments. A wider same-cardinality Hungarian replacement was tested on the
AT&T `onetone2` and `twotone` static-pivot cases: it could reduce off-diagonal
pivots on `twotone`, but increased fill and refactor time, and it worsened
`onetone2`. The retained policy therefore uses optional SPRAL Hungarian when
it improves matching cardinality over the in-tree matcher, and can also factor
a same-cardinality SPRAL Hungarian candidate after the first numeric factor
only for expensive high-off-diagonal-pivot cases. The latter path is gated by
actual factor work/fill and still keeps the SPRAL candidate only when the
factored numeric evidence improves. This keeps a license-compatible
MC64-adjacent source available for hard structural-deficit and expensive
dynamic-pivot cases without letting a maximum-product match replace already
accepted KLS row matchings solely on weight. KLS now reports accepted SPRAL
Hungarian or auction row permutations separately as
`selected_spral_matching`, because the auction path is LGPL-compatible but not
an exact assignment path and should not be conflated with
`selected_exact_matching`.

The SPRAL Hungarian/scaling path was then promoted to a bounded pre-factor
large-matrix candidate for weak-diagonal matrices whose symbolic analysis
leaves one dominant BTF block. This is the closest retained KLS path to the
MC64 preprocessing described by NICSLU and CKTSO, while remaining
LGPL-compatible because it uses the BSD-licensed SPRAL scaling subset rather
than HSL MC64 or restricted solver-tree copies. Same-session checks show the
tradeoff clearly. With `KLS_ENABLE_SPRAL_SCALING=ON`, `mac_econ_fwd500`
completed a one-factor/one-refactor run inside the 120s cap, selecting static
pivoting and reducing repeated factor/refactor averages to about 5.22s; the
same command timed out at 120s in the no-SPRAL build. `rajat30` selected the
SPRAL static match, reduced off-diagonal pivots to one, and modestly improved
repeated factor/refactor times to about 0.74s/0.73s. A broad version also
matched `ASIC_680k`, but that regressed an already fast many-block BTF case
from about 0.15s repeated refactors to about 0.27s, so the retained gate was
narrowed toward dominant BTF structure. Later focused checks with the pinned
SPRAL submodule showed the same boundary on paper-medium cases: `power197k`
selected SPRAL matching, cut off-diagonal pivots from about 51k to about 1.8k,
and reduced repeated refactor from about 0.055s to about 0.008s, while
`HTC_336_4438` moved from the intended METIS/no-BTF path to AMD/BTF and
regressed from about 0.103s to about 0.190s. The large SPRAL pre-static gate
therefore now requires an existing BTF symbolic analysis with a dominant block,
preserving no-BTF ordering choices. `pre2` still timed out at 120s after trying
the SPRAL path; its analyze-only evidence remains a dominant 629628-row block
with about `2.08e11` estimated flops. This confirms that license-compatible
MC64-adjacent preprocessing is useful and worth keeping, but it does not
replace the missing CKTSO/SubtreeLU-style intra-block numeric/scheduling
machinery.

Because the retained SPRAL path is now guarded and materially improves a
paper-medium hard row, the pinned BSD scaling subset was promoted from an
opt-in component to the default build, while keeping
`KLS_ENABLE_SPRAL_SCALING=OFF` for C-only builds and
`KLS_USE_SYSTEM_SPRAL=ON` for system installations. Fresh default-vs-SPRAL
checks show why this is aligned with the solver goal: `power197k` moves from
the no-SPRAL default path with about 51k off-diagonal pivots and roughly 0.057s
refactors to the SPRAL static-match path with about 1.8k off-diagonal pivots
and about 0.008s refactors. The same guard keeps `HTC_336_4438` on its
METIS/no-BTF path, while `transient`, `onetone1`, `onetone2`, and `rajat28`
remain on their existing accepted paths.

The fast-factor pivot-check path was then made more diagnostic by recording the
first rejected factor-order pivot and original matrix column in `kls_stats` and
benchmark JSON. It now also records the rejected BTF block start/size, the
simple suffix length from the rejected pivot to the end of the block, and the
exact U-pattern descendant tail size/work inside that block. It now also
records an ordered-block ETree successor-path size/work estimate for the first
rejected pivot, and a sorted pivoting-tail worklist scope seeded from the
current refresh state, matching the CKTSO paper's pivoting-tail upper-bound
dependency idea more closely. This does not implement CKTSO's pipelined tail
factorization, but it is a required prerequisite: KLS can now measure whether
failed fast factorizations reject near the tail, where a pivoting-tail restart
could avoid recomputing the whole block, or near the front, where full fallback
is still expected, and can distinguish a true no-pivot dependency tail from a
broad suffix, a single ETree successor path, and the actual pivoting-tail plan.

The same fast-factor path was then extended for scaled serial refactors. When a
scaled pattern is using fast factorization with pivot checks, KLS now runs its
own checked refactor loop and interrupts at the first unsafe multiplier instead
of running a complete KLU refactor and scanning the completed factors
afterward. This is still not CKTSO tail restart because fallback remains a full
pivoting factorization, but it narrows the wasted work before fallback and uses
the same rejected-pivot coordinate needed by a future ETree-descendant restart.

The unscaled fast-factor path was then extended with a block-local restart
primitive. When a checked no-pivot fast factorization rejects a reused pivot in
a BTF diagonal block, KLS can refactor only that BTF block with pivoting,
splice the block's new row order into the global numeric permutation, rebuild
the unscaled off-diagonal entries from the updated inverse permutation, and
retry the checked fast factorization when later columns may not have been
refreshed. When the rejected fast pass had already refreshed all columns, or
when a serial BTF pass rejected in the final block, KLS now validates the
repaired block tail directly and returns without a second full checked
refactor. Smoke tests now cover both the unscaled repair and the scaled
fallback path, and representative static-pivot cases (`gemat12`, `onetone2`,
`twotone`, and `rajat25`) did not trigger unexpected block restarts. This is
useful CKTSO-aligned infrastructure, but it is still not CKTSO's production
tail restart: it does not retain an ETree/EGraph tail or restart only
descendant work inside a large single block.

The selected-large KLU2 comparison was also run with the same 120s cap and one
factor/refactor repeat. KLU2 completed only `rajat29`, `rajat30`, and
`ASIC_680k`; it timed out on `G3_circuit`, `pre2`, `nxp1`, `Hamrle3`, and
`TSOPF_FS_b39_c30`. On successful common rows, KLU2 was faster than KLS on
`rajat29` but slower on `rajat30` and `ASIC_680k`. After the widened TSOPF
spiked rule, KLS also completes `TSOPF_FS_b39_c30` where KLU2 timed out. The
remaining large-case primary competitor is therefore CKTSO rather than KLU2.

Extending the existing mapped refactor metadata to scaled refactors was
retested against current mainline using a clean `HEAD` worktree. The scaled
map was KLU-semantics-compatible after recomputing row scale factors before
the mapped scatter and permuting them afterward, but it did not improve the
focused scaled large/medium cases. On the seven-case focused set the candidate
geomean regressed by about 1%, with no wins over 2%, so the experiment was
removed again. This keeps the retained refactor map limited to unscaled serial
patterns until a broader EGraph/separator-tree numeric kernel exists.

A narrower EGraph refactor consumer was then retained for the high-work
unscaled single-block class. It reuses the exact U-pattern level schedule,
assigns each level to static per-thread slices, and uses private dense scratch
per worker with a barrier between levels. The gate requires a single BTF block,
no KLU row scaling, at least 100k rows, at least `1e9` estimated no-pivot
dependency work, and level width at least four times the requested thread
count; scaled cases, many-block BTF cases, and smaller matrices continue using
the existing mapped or BTF-worker paths. In same-session checks on the default
build, `G3_circuit` improved from the earlier about 34s/34s factor/refactor
averages to about 27.9s/28.1s with valid residuals. Guard cases stayed on their
previous paths: `nxp1` and `rajat30` remained scaled, `ASIC_680k` remained
many-block BTF, and small `bcircuit` remained below the large-case threshold.
This is still not the full CKTSO cluster/pipeline scheduler, but it is a
retained general implementation of one paper idea where the current metadata
shows enough work to amortize synchronization.

The retained EGraph refactor was then extended with a CKTSO-style no-pivot
pipeline tail. Levels before the recorded split still use barriered cluster
mode. Tail columns are assigned to workers without per-level barriers, and a
column waits only for the actual U-pattern predecessors it consumes. Completion
is published with C11 atomics, so LU writes from predecessor columns are visible
before dependent columns update. On `G3_circuit`, the same one-factor and
one-refactor focused check improved from about 27.9s/28.1s after the first
EGraph path to about 19.5s/19.4s with valid residuals. `G2_circuit` also
completed correctly on the same path. Scaled single-block and many-block guard
cases (`nxp1`, `rajat30`, `ASIC_680k`, and small `bcircuit`) stayed on their
existing paths. This closes the no-pivot EGraph cluster/pipeline piece for the
current KLS-owned storage, but not CKTSO's ETree-descendant restart with
pivoting after a failed pivot check.

The same EGraph cluster/pipeline refactor was then extended to large
single-block KLU row-scaled factors. The scaled path recomputes KLU's row
scale vector before the EGraph numeric update, divides each mapped matrix entry
by the unpermuted scale for its original row, and then permutes `Rs` back into
pivot order after the refactor. The normal mapped refactor path still rejects
scaled factors, so scaled maps are built lazily only when the EGraph path runs.
On same-session focused checks, `nxp1` improved from about 1.02s/1.01s
factor/refactor averages to about 0.73s/0.71s, and `rajat30` improved from
about 0.65s/0.64s to about 0.38s/0.39s, all with valid residuals. Guard cases
remained on their existing paths: `ASIC_680k` stayed a many-block BTF case,
`rajat29` stayed below the high-work EGraph gate, and small `bcircuit` remained
below the size threshold.

The EGraph pipeline tail then replaced mutex-based stop polling with an atomic
stop flag while keeping error details protected by the existing mutex. This is
a small synchronization reduction in the common no-error path: workers no
longer take a mutex while polling predecessor completion or checking level
boundaries. Same-session focused checks improved `G3_circuit` to about
19.0s/18.9s factor/refactor, `nxp1` to about 0.72s/0.71s, and `rajat30` to
about 0.37s/0.36s, with valid residuals.

A follow-up attempt to cache the per-worker EGraph dense scratch vectors was
rejected. A dedicated solver-owned scratch cache improved repeated `rajat30`
factor/refactor averages to about 0.34s/0.34s, but `nxp1` regressed slightly
to about 0.72s/0.72s and the long `G3_circuit` guard regressed to about
19.0s/19.2s. Since the effect was not a general win, the experiment was
removed instead of adding a size or matrix-shape gate.

The MC64 compatibility boundary was rechecked after allowing existing code if
it remains LGPL-compatible. The retained vendored route is still SPRAL's
BSD-3-Clause scaling subset: it is redistribution-compatible with KLS's
LGPL-2.1-or-later license, whereas HSL MC64 itself and restricted MC64 copies
from solver trees remain out of scope for vendoring. The policy does not
require every MC64-style implementation to originate in KLS; it requires any
copied or vendored implementation to be redistributable inside an LGPL KLS
distribution. CMake now enforces that distinction for the system-SPRAL path by
requiring `KLS_SYSTEM_SPRAL_LGPL_COMPATIBLE=ON`; the bundled path checks the
pinned SPRAL `LICENCE` file before building its scaling subset. A small BSD
Rust `mc64` crate exists as a partial SPRAL translation, but it does not
improve KLS's C integration story over the already pinned SPRAL Fortran/C
interface. Fresh SPRAL-enabled checks also confirm the policy should stay
guarded rather than become an unconditional default:
`rajat30` selected SPRAL matching, reduced off-diagonal pivots to one, and cut
initial factor time to about 5.7s, but its repeat-heavy factor/refactor
averages were about 0.67s/0.66s versus the faster current no-SPRAL EGraph path.
`nxp1` did not select SPRAL and stayed roughly neutral-to-slightly-worse. This
keeps license-compatible MC64-style code in KLS, but points the large remaining
CKTSO gap back to numeric scheduling and pivoting machinery rather than merely
importing another MC64 copy.

The EGraph refactor consumer was then generalized from hard-coded single-block
indices to block-local BTF indices. The same guarded cluster/pipeline schedule
can now run inside a large dominant BTF block, update the block-local LU
columns, and refresh Offx for entries above that BTF block. The dispatch remains
conservative: it requires a large block covering at least 75% of the matrix, so
ordinary many-block cases such as `ASIC_680k` stay on the existing BTF worker
pool. On same-session checks, the committed HEAD baseline timed out at 130s on
`mac_econ_fwd500` with one factor and one refactor, while the block-aware
EGraph path completed with valid residuals and about 14.9s/14.9s repeated
factor/refactor averages. Single-block guards remained in their previous
range: `nxp1` measured about 0.71s/0.71s and `rajat30` about 0.34s/0.33s in
focused repeated runs; `G3_circuit` stayed around 19.1s/19.1s. This is useful
dominant-block coverage, but not the missing major CKTSO algorithm: `pre2`
still timed out under a 125s cap.

A current large-recon run after the dominant-BTF change still shows the
remaining gap clearly. KLS completed six of the eight selected large cases and
timed out on `pre2` and `Hamrle3` under a 120s per-matrix cap. Against the
existing CKTSO four-thread run, KLS still lost all six common completed cases,
with the largest ratios on `nxp1`, `ASIC_680k`, `G3_circuit`, and `rajat30`.
The same run beat KLU2 on the larger common cases `ASIC_680k` and `rajat30`,
but still lost `rajat29`. A factor-only `pre2` probe also timed out under
125s, confirming that the unresolved `pre2` gap is first-factor numeric
machinery, not repeated-refactor scheduling.

The June 28, 2026 auto-input-width rerun keeps that conclusion. Artifact
`build/kls_input_auto_large_recon_t4_r1_ref1_timeout120.jsonl` completed six of
eight selected large rows with 32-bit input ingestion and a 38.19s completed-row
geomean; `pre2` and `Hamrle3` timed out at 120s. Against
`build/cktso_paper_large_recon_t4_timeout120.jsonl`, KLS still wins
`TSOPF_FS_b39_c30` decisively but loses the other completed common rows:
`rajat30` about 1.50x, `nxp1` about 1.40x, `ASIC_680k` about 1.30x,
`G3_circuit` about 1.08x, and `rajat29` about 1.02x. CKTSO completes `pre2`
in the saved artifact and times out only on `Hamrle3`. A same-session `pre2`
isolation sweep found that analyze-only normal AMD completes in about 5s with
a 629628-row dominant block, about 61.1M L and U entries, and about
`2.08e11` estimated flops, while factor-only default AMD, forced
`KLS_ENABLE_KLS_FIRST_FACTOR`, transpose AMD, METIS, and no-static-pivoting
probes all timed out at 120s. This makes `pre2` a cold first-factor kernel gap,
not an input-index, orientation, ordering, static-pivoting, or refactor-repeat
artifact.

A low-work dominant-BTF guard was then retained for Rajat-family large cases:
if BTF finds one block covering at least 95% of the matrix, the stripped fringe
is at most 5% but still nontrivial, and the symbolic flop estimate is below
`1e9`, auto keeps BTF and starts with no KLU row scaling. On `rajat29`, this
changed auto from AMD/no-BTF/max-scale to AMD/BTF/no-scale, reducing the
SPICE-cycle estimate from the old about 5.67s to about 3.74s. That now beats
the saved CKTSO artifact at about 5.28s and the saved KLU2 artifact at about
4.29s. Higher-work guards such as `rajat30`, `nxp1`, `Raj1`, and `rajat24`
still take the no-BTF retry when their symbolic evidence supports it.

A symbolic METIS retry was then retained before numeric factorization for
large high-work single-block no-BTF analyses. This is the same CKTSO/SubtreeLU
nested-dissection direction as the existing post-factor METIS promotion, but
it avoids first paying for an AMD numeric factorization when METIS already has
a clearly lower symbolic fill score. On same-session serial checks, `nxp1`
kept the same METIS/no-BTF numeric path but reduced the SPICE-cycle estimate
from the saved about 89.6s to about 78.9s, and `rajat30` moved from about
49.2s to about 42.4s. The low-work dominant-BTF guard stayed on `rajat29`, and
the many-block BTF guard stayed on `ASIC_680k`, so this is a general
symbolic-cost improvement rather than a benchmark-name policy.

The auto-scale gate was then tightened for the same large high-work
METIS/no-BTF single-block class. Earlier sweeps had already shown the completed
hard large cases should keep max row scaling, but auto mode still tried other
scale modes after the first factorization and rejected them. Skipping those
post-factor scale trials preserves the final factors and residuals while
cutting setup time: `rajat30` auto initial factor moved from about 4.06s to
about 1.89s, matching explicit METIS/no-BTF/max-scale, and `nxp1` moved to
about 2.77s while keeping the same METIS/no-BTF/max-scale EGraph path.
`rajat29` and `ASIC_680k` remained on their retained BTF policies.

The high-work METIS/no-BTF class was then moved one step earlier in analysis.
A structural direct-start gate now recognizes large nearly diagonal matrices
with a meaningful row/column spike in either a sparse-spike or dense-spike
density band. This lets `nxp1` and `rajat30` start directly with METIS/no-BTF
instead of first doing AMD symbolic analysis and then a METIS symbolic retry.
Same-session checks reduced `rajat30` analysis from about 3.9s to about 3.0s
while preserving the about 1.86s initial factor and about 0.31s refactor path,
and reduced `nxp1` analysis from about 2.6s to about 1.84s while preserving
the max-scale EGraph path. The low-work `rajat29` BTF guard and many-block
`ASIC_680k` guard stayed on their retained policies.

The same "start with the already accepted policy" idea was applied to
TSOPF-style spiked low-diagonal matrices. Auto mode had already learned that
`TSOPF_FS_b39_c30` wanted METIS/BTF, sum scaling, and `1e-4` pivot tolerance,
but it reached that state by first paying for max-scale/default-tolerance
factorizations. The structural dominant-BTF gate now starts this class directly
with sum scaling and `1e-4` tolerance and skips the redundant scale and pivot
trials. A focused `TSOPF_FS_b39_c30` check reduced initial factor time from
about 86s to about 12s while preserving the about 2.6s refactor path and valid
residual, moving the SPICE-cycle estimate from about 370s to about 290s versus
the saved CKTSO artifact at about 333s.

The same structural idea was then extended downward to small spiked
low-diagonal TSOPF/QY cases. Across the medium manifest, the small rule matches
the `TSOPF_FS_b9_c1`, `TSOPF_FS_b9_c6`, and `case9` shapes: 2k-20k rows,
about 8-14 nonzeros per row, 20-35% diagonal coverage, and row/column degree
spikes around 40-60% of the matrix order. Focused sweeps showed these cases
prefer METIS/BTF, KLU's unscaled `0` mode, and `1e-4` pivot tolerance. The
static-pivot trial was also skipped for this class because it added setup cost
without being selected. Current auto checks moved `TSOPF_FS_b9_c1` to about
0.13s on the SPICE-cycle estimate, `TSOPF_FS_b9_c6` to about 1.9s, and
`case9` to about 1.9s with valid residuals. The saved CKTSO artifact is still
faster on `TSOPF_FS_b9_c1` at about 0.061s, but KLS is ahead on the larger
`b9_c6` and `case9` rows.

EGraph schedule construction was then tightened to the same structural class
as the retained EGraph consumer. KLS had been building dependency-level
metadata for many-block and low-work BTF cases that could not use the
cluster/pipeline path, including `ASIC_680k` and the retained `rajat29`
low-work dominant-BTF policy. The new gate keeps schedules for high-work
single-block cases such as `nxp1` and `rajat30`, but skips them for the BTF
worker-pool and low-work mapped paths. Same-session checks showed schedule
metrics dropping to zero for `ASIC_680k` and `rajat29`, with `ASIC_680k`
initial factor/setup moving from about 1.59s to about 1.54s and `rajat29`
from about 0.26s to about 0.24s, while the high-work EGraph guards retained
their schedule metadata.

A CKTSO-style dynamic atomic assignment prototype for the EGraph pipeline tail
was tested and rejected. It replaced the static per-thread tail stride with a
shared atomic cursor, but `nxp1` and `rajat30` were neutral-to-slightly worse
and `G3_circuit` regressed to about 19.45s/19.41s factor/refactor. The
existing static tail assignment therefore remains the better fit for KLS's
current column storage and scratch model.

A CKTSO-style topological level-order all-pipeline cursor for huge single-block
EGraph refactors was also tested and rejected in favor of the current natural
column-order cursor. On same-session checks it regressed `nxp1` refactor time
to about 0.72s and `rajat30` to about 0.39s, so the missing gap is not simply
the absence of a level-ordered all-pipeline cursor. The retained improvement in
this area is narrower: scaled single-block EGraph refactors now dispatch to a
dedicated hot kernel that applies row scaling directly while loading the fixed
input-position map. Same-session probes improved `nxp1` repeated refactor from
about 0.301s to 0.287s and `rajat30` from about 0.264s to 0.218s with valid
residuals, while unscaled/BTF probes stayed valid.

Disabling the huge-single all-pipeline gate was then tested as a direct
cluster/pipeline split experiment. It regressed `nxp1` repeated refactor to
about 0.65s and `G3_circuit` to about 19.43s, so the current all-pipeline
natural cursor remains the better general choice for huge single-block factors.
A wider eight-way scalar scatter-subtract unroll was also tested and rejected:
it regressed `nxp1` repeated refactor to about 0.38s and did not provide a
general win on the quick large guards. The existing four-way generic scalar
scatter kernel is retained.

The BTF worker pool was then adjusted to fetch small ranges of diagonal blocks
per mutex acquisition when a matrix has many thousands of non-dominant BTF
blocks. This targets scheduler overhead on ASIC-style matrices with hundreds
of thousands of tiny blocks, while dominant-block and smaller-BTF cases still
fetch one block at a time for load balance. On same-session checks against a
clean baseline, `ASIC_680k` moved from about 0.173s/0.157s repeated
factor/refactor averages to about 0.159s/0.156s, while `ASIC_680ks` moved from
about 0.163s/0.170s to about 0.173s/0.163s. The net SPICE-cycle effect is
small but positive on the many-block ASIC guards, so this is a retained
scheduler-overhead cleanup, not a CKTSO-scale algorithmic fix.

The EGraph cluster scheduler was then adjusted from equal per-thread column
slices to contiguous slices balanced by the existing no-pivot column-work
estimate. This is a small retained part of the CKTSO/SubtreeLU load-balance
idea, scoped to the barriered cluster levels and leaving the existing static
stride pipeline tail unchanged. Same-session clean-baseline checks showed
modest but consistent wins on the EGraph path: `rajat30` moved from about
0.335s/0.330s repeated factor/refactor averages to about 0.317s/0.314s,
`nxp1` moved from about 0.710s/0.706s to about 0.691s/0.669s, and
`G3_circuit` moved from about 19.40s/19.33s to about 18.41s/18.45s. This
improves the retained no-pivot refactor consumer, but it does not address the
remaining first-factor timeout on `pre2` or implement CKTSO's pivoting tail
restart.

The dominant-BTF EGraph gate was then lowered for medium ASIC-style matrices
whose largest diagonal block is at least 90k rows, covers at least 95% of the
matrix, and has at least `5e8` actual factor flops. The EGraph consumer itself
still requires at least `2.5e8` no-pivot dependency work and enough level width,
so low-work dominant-BTF cases such as `rajat29` skip schedule construction.
This lets the retained cluster/pipeline refactor run inside the dominant block
of the `ASIC_100k` and `ASIC_320k` families while ordinary many-block
`ASIC_680k`/`ASIC_680ks` shapes remain on their existing paths. Focused checks
showed the largest wins on `ASIC_320k` and `ASIC_320ks`: repeated refactor time
moved from roughly 0.35s/0.29s to about 0.11s/0.08s with valid residuals.
`ASIC_100k` and `ASIC_100ks` improved more modestly, while the retained guard
kept `rajat29` schedule metrics at zero. This is a useful CKTSO-inspired
coverage extension, but CKTSO remains faster on these ASIC rows.

The same gate was then extended below the 90k dominant-block floor only for
high-work dominant-BTF shapes whose largest block covers at least 95% of the
matrix and whose measured factorization has at least `5e9` flops. This keeps
the rule tied to general work evidence rather than TSOPF names. On
`TSOPF_FS_b39_c19`, which has a 76215-row dominant block and about `1.54e10`
factor flops, the retained EGraph path reduced repeated factor/refactor
averages from about 4.9s/4.8s to about 1.5s/1.5s with valid residuals, cutting
the SPICE-cycle estimate to about 167s versus the saved CKTSO result near
603s. The previously rejected Rajat 80k-row class stays excluded because it
does not meet the high-work gate; focused `rajat28` checks still recorded zero
EGraph schedule metrics and stayed on the existing path.

The dominant-BTF EGraph gate was also probed down to 80k rows and 80k-row
largest blocks to see if the same policy should cover high-flop Rajat
dominant-BTF rows. It activated on `rajat20`, `rajat25`, and `rajat28`, but
slowed repeated refactors from the existing roughly 0.15s class to roughly
0.18-0.19s. The 90k retained floor remains the better general rule; the Rajat
rows need different numeric scheduling or pivoting work, not more EGraph
coverage with the current kernel.

The retained EGraph gate was later extended to a different dominant-BTF shape:
large-heavy blocks below the 95% coverage floor. The new branch requires the
largest block to cover at least 85% of the matrix, contain at least 100k rows,
have at most 20k total BTF blocks, and show at least `2e9` actual factor flops.
This activates on AT&T `twotone` but not on the previously rejected 80k-row
Rajat class. In clean checks, `twotone` built 2391 EGraph levels and moved
repeated refactor from the current worker-path class around 0.97-1.26s to about
0.34s, reducing the 100-cycle SPICE estimate from roughly 101s to roughly 38s.
Guards kept `rajat28` schedule metrics at zero, kept `ASIC_680ks` on the
many-small-block worker path, preserved `G2_circuit`, and restored
`ASIC_320k`/`ASIC_320ks` to their existing METIS/BTF EGraph path after the
no-BTF retry fix.

The same coverage-bounded idea was then extended to medium-heavy dominant BTF
blocks. The initial retained branch required 85-95% dominant-block coverage, at
least a 30k-row largest block, at most 5k BTF blocks, and at least `5e8` actual
factor flops; the upper coverage bound was deliberate so the previously
rejected 95%+ Rajat class remained excluded. In clean focused checks,
`onetone1` built 1705 EGraph levels and reduced repeated refactor from about
0.18s to about 0.06-0.065s, cutting the SPICE-cycle estimate from roughly 19s
to about 7s. `onetone2` stayed below that first work threshold with zero
schedule metrics, `rajat28` remained excluded, and `twotone` stayed on the
large-heavy EGraph path.

The high-coverage dominant-BTF schedule floor was then lowered separately from
the 85-95% medium-heavy branch. The retained rule still requires at least
`2e8` actual factor flops before building EGraph metadata and at least `1e8`
computed dependency work before the consumer runs. This activates a lower-work
95%+ dominant-block case without reopening the rejected Rajat class: `transient`
built 807 EGraph levels and reduced repeated refactor from about 0.057s to about
0.024s, cutting the focused SPICE-cycle estimate from about 7.0s to about
3.8s. At that point `onetone2` remained below the `2e8` high-coverage
factor-work floor with zero schedule metrics, `rajat28` remained excluded,
`onetone1` stayed on the medium-heavy EGraph branch, and `ASIC_320k` stayed on
the existing high-coverage path.

The EGraph worker launch path was also retried with a solver-owned persistent
worker/scratch pool, analogous to the retained BTF refactor pool. This was
removed before commit because it was not a general win after the cluster
work-balancing change: `rajat30` improved slightly from about 0.325s/0.313s
to about 0.316s/0.310s factor/refactor averages, `nxp1` was neutral to
slightly worse at about 0.685s/0.686s versus 0.684s/0.685s, and `G3_circuit`
regressed from about 18.40s/18.49s to about 18.73s/18.69s. The evidence points
back to numeric update structure, pivoting-tail machinery, or separator-tree
scheduling rather than thread-launch overhead.

The EGraph cluster/pipeline split threshold was also checked after the
work-balanced cluster slices. The retained policy switches to pipeline mode at
the first level narrower than `2 * threads`. A later handoff at `1 * threads`
cut the `rajat30` pipeline tail to 1156 columns but worsened factor/refactor
averages to about 0.346s/0.330s. An earlier handoff at `4 * threads` expanded
the tail to 3715 columns and measured about 0.328s/0.317s. The current
middle split remains the better general setting in this quick check.

The generic EGraph column kernel was also tested with a single-block fast path
that bypassed the per-column BTF block lookup and `R` checks when
`nblocks == 1`. This was removed before commit because it did not help the
scaled single-block guard: `rajat30` was only slightly positive at about
0.317s/0.313s versus 0.322s/0.314s, while `nxp1` regressed from about
0.671s/0.675s to repeated samples around 0.693s/0.675s and 0.696s/0.693s.
The block-lookup branch is therefore not the next useful source of the CKTSO
gap.

The EGraph pipeline completion array was also tested with retained per-column
level labels so cluster-phase columns could skip atomic completion stores and
pipeline waits for predecessors known to be before the split level. This was
removed before commit because it helped `rajat30` and repeated-factor `nxp1`
but regressed the largest EGraph guard: `rajat30` moved from about
0.316s/0.314s to 0.310s/0.305s, `nxp1` moved from about 0.687s/0.673s to
0.665s/0.674s, but `G3_circuit` regressed from about 18.39s/18.38s to
18.73s/18.75s. The retained pipeline publication remains the simpler
per-completed-column atomic store until a fuller scheduler changes the tail
execution model.

The remaining `pre2` timeout was rechecked after the MC64 licensing boundary
was clarified. A SPRAL-enabled build, using the retained BSD-licensed
Hungarian/scaling path, still timed out under a 180s one-factor cap. Lowering
the initial pivot tolerance also did not provide a usable CKTSO-style
pivot-reuse substitute: `--pivot-tol 0` reached a singular factor quickly,
while `1e-8` and `1e-4` still timed out under 120s. On completed large guards,
lower tolerances were not a general win: `rajat30` at `1e-4` slowed initial
factor and repeated refactor, `rajat30` at `1e-8` improved only refactor while
worsening initial factor and conditioning, and both `nxp1` lower-tolerance
checks regressed. This keeps `pre2` in the missing pivoting-tail/numeric-kernel
bucket rather than the tuning bucket.

A narrow many-block BTF worker-map retry was also prototyped for ASIC-style
structures after the CKTSO comparison showed a large `ASIC_680k` gap. The
prototype made the worker map path scale-aware and built the precomputed
input-position map only for non-dominant BTF patterns with at least 100k blocks
and at least half as many blocks as rows. It was removed because it was not a
general win: `ASIC_680k` had only a small refactor improvement but slower
repeated factor and setup, while `ASIC_680ks` regressed in both factor and
refactor. This confirms the previous broad worker-map rejection and points
ASIC-style gaps toward a different small-block numeric/storage kernel rather
than passing the existing map through the worker pool.

A separate worker-map policy was then retained for the opposite many-block
shape: dominant BTF matrices with a bounded block count and bounded input
size. KLS now builds the existing fixed-pivot input-position map for
pool-eligible dominant BTF cases when the largest block covers at least 75% of
the matrix, the block count is at most 20k, and the input has at most 3M
nonzeros. The worker-pool map path is scale-aware, so scaled Rajat cases can
reuse it without falling back to the original `Q`/`Pinv` scan. This is still
not the CKTSO numeric kernel, but it removes a general repeated-refactor
overhead from Rajat/AT&T-style dominant-block patterns while preserving the
previous ASIC tiny-block rejection. Focused one-pass checks improved
`twotone` from about 129s to about 101s on the SPICE-cycle estimate,
`onetone2` from about 5.3s to about 4.6s, and `rajat20`/`rajat25`/`rajat28`
from about 16.6s/18.2s/17.6s to about 14.7s/15.3s/15.5s. The large guards
remain on their intended paths: `rajat29` is excluded by the 3M-nnz gate and
stays near the 4.1s class, while `ASIC_680k` remains an extreme many-block
case. A full one-pass medium manifest with the retained policy completed 90 of
93 rows under the 120s cap, with the expected `ss1`/`mac_econ_fwd500`
timeouts and the known singular `bips07_1998`; common rows versus the previous
broad KLS artifact improved geomean by about 2%, though one-pass noise still
dominates many sub-millisecond rows.

The single-block EGraph schedule/consumer gate was then lowered from the
previous very-large-only floor to cover moderate single-block cases whose
actual factor work, LU fill, and measured dependency work are already high
enough to amortize schedule construction and worker scratch setup. This is a
paper-derived generalization of the CKTSO/NICSLU intra-block dependency-graph
idea, not a matrix-name rule: single-block scheduling now starts at about
`3e8` factor flops or 3M factor nonzeros, and the no-pivot EGraph consumer
requires about `1.5e8` dependency-work units. In same-session SPRAL-enabled
focused checks with four threads, `HTC_336_4438` changed from schedule-only
metadata to actual EGraph consumption and its repeated refactor average dropped
from about 0.105s to 0.041s; `rajat24` similarly dropped to about 0.096s
refactor average. A nine-row hard-focus JSONL improved geomean SPICE-cycle time
by about 1.17x versus the previous default-on artifact, with no loss over 2%.
The same comparison against CKTSO still leaves a large gap, about 2.91x
geomean on those focused rows, so this is retained as a useful KLS-owned
refactor threshold improvement rather than mistaken for CKTSO's full pivoting
tail scheduler.

A later moderate single-block refinement lowered the single-block EGraph floor
for unscaled one-block matrices between 30k and 100k rows when the factored
numeric object already has at least `5e7` measured flops, at least 1M LU
entries, enough level width, and at least `2e7` measured dependency-work units.
This is intentionally below the large single-block floor but still excludes
low-work single-block rows such as `bcircuit`, `ACTIVSg10K`, `ACTIVSg70K`, and
`OPF_10000`. In the current medium artifact the selector matches only
`rajat15`; it builds 633 levels and reduces repeated refactor from about
`0.0173s` to about `0.0101s`, moving the focused SPICE-cycle median from about
`2.16s` to about `1.47s`. A full one-pass medium run kept the same three
known failures and moved KLS geomean from the previous `0.3033s` artifact to
about `0.3023s`; the CKTSO ratio improved to about `1.145x` slower on the 90
common completed rows, and the KLU2 comparison improved to about a `2.02x`
geomean speedup on the 88 common completed rows. Large single-block guards
such as `G2_circuit` and `mc2depi` stay on the existing large EGraph path.

The medium-heavy dominant-BTF gate was then rechecked after the newer EGraph
cluster/pipeline scheduler and the single-block threshold work. The current
retained branch lowers the 85-95% coverage class to about `1.5e8` actual
factor flops and lets that class consume the EGraph path when measured
dependency work reaches about `8e7`. This brings `onetone2` into the same
structural policy as `onetone1`, without reopening the rejected 95%+ Rajat
class: `rajat28` still reports zero EGraph schedule metrics. In same-session
SPRAL-enabled focused checks with four threads, `onetone2` built 1010 EGraph
levels and repeated refactor dropped from about 0.040s to about 0.015s. The
nine-row hard-focus JSONL improved geomean SPICE-cycle time by about 1.09x
over the previous single-block-threshold artifact, and the KLS/CKTSO focused
geomean gap moved from about 2.91x to about 2.68x. CKTSO is still materially
faster, so the remaining gap still points to the larger pivoting-tail,
numeric-kernel, and solve-scheduler items rather than more ordering backends.

The single-block EGraph column fast path was then revisited with a narrower
scope than the earlier rejected scaled experiment. The retained version only
applies when the factor has one BTF block and no active KLU row scaling, so the
previous `nxp1`/`rajat30` scaled guards remain on the generic column kernel. In
that unscaled class, the hot loop bypasses repeated BTF-block lookup, BTF
off-block checks, and scaling branches. Same-session focused checks showed a
small net gain: the nine-row hard-focus JSONL improved by about 0.7% geomean
over the medium-dominant-block artifact, with `G2_circuit` moving from about
48.8s to about 47.8s and `HTC_336_4438` from about 8.64s to about 8.25s on the
SPICE-cycle estimate. A one-refactor `mc2depi` guard stayed in the same
refactor class, about 2.23s. This is retained as a minor EGraph kernel cleanup,
not as the missing CKTSO-scale scheduler.

Two follow-up EGraph hot-loop/scheduler probes were rejected after the retained
single-block fast path. First, the unscaled single-block fast path was changed
to trust the previously validated refactor map and U-pattern schedule, removing
per-entry bounds and dependency-order checks from the hot loop. This did not
help `G2_circuit`: refactor stayed in the same noisy 0.45s class, so the
validation branches are not the visible CKTSO gap. Second, the pipeline tail was
changed from the retained round-robin static assignment to contiguous
work-balanced ranges over the tail columns. This was a clear regression on
`G2_circuit`, moving refactor to about 0.94s because later contiguous ranges
wait behind earlier dependency ranges. The retained round-robin tail assignment
therefore remains the right fit for the current EGraph representation. A
separate LU pointer-cache prototype was also removed: caching L/U index and
value pointers with the schedule added memory and regressed the primary
unscaled EGraph guards (`G2_circuit` and `HTC_336_4438`) despite noise-driven
improvement on `rajat28`, which does not consume that cache. These results
make it unlikely that more small EGraph bookkeeping reductions will close the
remaining CKTSO gap.

Three later probes reached the same conclusion. Disabling the forced
all-pipeline path for huge unscaled single-block EGraph factors moved
`G2_circuit` and `mc2depi` backward, with only a small noisy `rajat30`
improvement, so the retained all-pipeline shape is still the better general
dispatch for that class. Retaining worker scratch buffers and stamped
pipeline-completion storage across refactors was numerically correct but
slower on the primary EGraph rows, which means per-refactor allocation is not
the main visible overhead. A focused `nxp1` unscaled-scale trial also did not
produce a safe general policy: no-scale candidates kept valid residuals and
sometimes lowered repeated refactor time, but the result was noisy, estimated
conditioning dropped by about five orders of magnitude versus max scaling, and
post-factor trial cost erased the possible cycle gain.

A follow-up `rajat28` policy sweep also confirmed that the remaining worst
focused-row gap is not a missing scale-mode or static-pivoting toggle. With the
retained static-pivoted AMD/BTF path, auto, no-scale, sum-scale, and max-scale
variants all stayed in the roughly 0.15-0.17s repeated-refactor class, while
CKTSO's saved four-thread artifact is about 0.011s. Disabling static pivoting
was worse: METIS/scale-auto paid about 10.37s initial factor and stayed around
0.21s repeated refactor, and METIS/no-scale paid about 13.93s initial factor
and about 1.54s repeated refactor. This keeps `rajat28` in the dominant-block
numeric-kernel bucket, not the MC64/preprocessing bucket. The MC64-compatible
boundary remains unchanged: use the pinned BSD-licensed SPRAL scaling subset,
a compatible system SPRAL, or independent KLS code; do not vendor HSL MC64 or
solver-tree copies that retain HSL redistribution restrictions.

The BTF worker-pool gate was then narrowed for low-work dominant decompositions
that have one 95%+ diagonal block, thousands of tiny fringe blocks, a largest
block below the current EGraph size floor, and less than `1e9` measured factor
flops. This is not matrix-name tuning: it is the structural case where the
dominant block still runs in one worker and the remaining fringe work is too
small to pay for threaded block scheduling. A broad first version also caught
`ckt11752_dc_1`, whose 172 BTF blocks still benefited from the pool, so the
retained gate requires at least 1024 BTF blocks. On the affected target set
(`ckt11752_dc_1`, `LeGresley_87936`, `rajat20`, `rajat25`, `rajat28`) the
narrow gate kept `ckt11752_dc_1` neutral, had no losses above 2%, and improved
geomean SPICE-cycle time by about 2.5% versus the prior dynamic-pipeline medium
artifact. The hard-focus run stayed essentially neutral against the second
dynamic-pipeline baseline (`9.07s` versus `9.04s` geomean) while improving
`rajat28` from about `17.36s` to about `15.84s` in that comparison. This is a
small policy cleanup. It does not change the main conclusion that the large
CKTSO gap on `rajat20/25/28`, `G2_circuit`, and similar rows requires a
CKTSO/SubtreeLU-style row-oriented numeric kernel, pivoting-tail restart, or
separator-tree/private-pipeline scheduler rather than another ordering or MC64
import.

To keep that diagnosis reproducible, `scripts/decompose_solver_gap.py` now
compares two benchmark JSONL files by phase contribution. On the current
nine-row hard-focus comparison against the saved CKTSO medium artifact, KLS is
still dominated by repeated refactorization: `rajat28` spends about 94% of its
SPICE-cycle estimate in repeated refactor work and that refactor component is
about 13.6x CKTSO's; `G2_circuit` spends about 93% there and is about 6.7x
CKTSO's; `onetone1`, `onetone2`, `rajat24`, and `transient` also have
refactor-component ratios around 2.7x to 3.4x. Solve ratios on those same rows
are only about 1.0x to 1.3x, and `twotone`/`power197k` solve is already faster
than CKTSO. This confirms that a CKTSO-style solve rewrite is secondary for the
current hard gap; the larger missing mechanism is the row-oriented
factor/refactor engine and its pivot-aware scheduler.

A narrow solve-workspace cache was tested and rejected after this diagnosis. The
prototype kept a solver-owned dense permutation buffer for the `row_perm` solve
path so static-pivot and exact-matching cases would not allocate/free an
`n`-entry buffer on every solve. A repeat-20 subset was tempting, with solve
geomean `0.945x` on `power197k`, `rajat20`, `onetone2`, `OPF_10000`,
`LeGresley_87936`, and a non-`row_perm` `G2_circuit` guard, but the broader
row-permutation paper set rejected it: `power197k`, `rajat20`, `rajat25`,
`rajat28`, `onetone1`, `onetone2`, `twotone`, `OPF_10000`,
`LeGresley_87936`, `rajat22`, `rajat23`, `rajat24`, `hvdc1`, and `hvdc2`
measured `1.007x` slower solve geomean against `a0583c5`, with the static
non-exact rows at `1.010x` slower. The current per-call buffer stays until the
larger KLS-owned LU storage makes a structure-adaptive solve worthwhile.

The scaled serial mapped BTF refactor was then enabled for the narrow
many-fringe dominant-BTF shape that the worker-pool narrowing intentionally
left serial: at least 1024 BTF blocks, a 95%+ largest block below the EGraph
floor, and at least `1e8` measured factor flops. The broad scale-aware mapped
prototype was rejected because it regressed unrelated low-work or small-fringe
scaled rows such as `dc1` and `ckt11752_dc_1`. The retained guard is structural
and, on the current medium paper artifact, selects only `rajat20` and
`rajat28`. On the targeted dominant-fringe set it improved geomean SPICE-cycle
time by about 2.4% with no losses over 2%; on the nine-row hard-focus set it was
about 1.2% faster than the previous many-fringe-serial artifact, mostly from
`rajat28`. The full 93-row medium manifest completed the same 90 rows as the
previous KLS artifact, with the known singular `bips07_1998` and the known
`ss1`/`mac_econ_fwd500` timeouts, and improved one-pass geomean from about
`0.3407s` to `0.3372s`. Against the saved CKTSO medium artifact, however, KLS
is still about `1.28x` slower geomean and `rajat28` remains about `9.3x` slower,
so this is another small fixed-pattern refactor cleanup rather than the missing
CKTSO-scale mechanism.

That same many-fringe dominant-BTF class was then revisited using the paper
idea that cluster barriers can dominate when exact dependencies are already
known. Instead of lowering the old barriered EGraph gate again, KLS now runs
this class through an all-pipeline exact-EGraph consumer: workers claim columns
in topological order and wait only for actual U-pattern predecessors. In the
current 93-row medium artifact this structural gate selects `rajat20`,
`rajat25`, and `rajat28`. On the repeat-heavy target set it improved geomean
SPICE-cycle time by about `1.71x` versus the previous guarded-scaled-map
artifact, with `rajat20/25/28` each moving to roughly 40% of their prior cycle
time. The full medium manifest kept the same 90 completed rows and the same
three failures, improving KLS one-pass geomean from about `0.3372s` to
`0.3275s`. Against CKTSO, the medium geomean gap moved from about `1.28x` to
about `1.24x`, and `rajat28` moved from about `9.3x` slower to about `3.8x`
slower. This is meaningful scheduler progress, but the remaining gap on
`G2_circuit`, ASIC rows, and the Rajat rows still points to the larger
row-oriented numeric kernel, pivot-aware scheduler, and separator/supernode
work described in the papers.

A small reactive static-pivoting payback gate was then retained. Earlier KLS
policy could launch the post-factor static row-matching trial on small,
low-work matrices where the accepted permutation improved pivoting and fill but
did not recover its setup cost over the 99-refactor SPICE-cycle estimate. The
new gate applies only after a first factorization has succeeded, only below
20k rows, and only when diagonal weakness is not severe; severe missing-diagonal
cases still use the pre-static path, and larger or higher-work cases keep the
existing matching policy. On the full medium paper manifest this kept the same
90 completed rows and the same known `bips07_1998`, `ss1`, and
`mac_econ_fwd500` failures, while improving KLS geomean from about `0.3141s`
to about `0.3096s`. Against the saved CKTSO artifact the geomean ratio moved
from about `1.190x` slower to about `1.173x` slower. The largest retained wins
were low-work `OPF_3754`, `bips98_*`, and `nopss_11k` cases; the phase
decomposition after this change still shows the remaining largest losses are
dominated by repeated refactor throughput (`G2_circuit`, ASIC, `mc2depi`, and
Rajat rows), not by another MC64-compatible matching import.

A follow-up low-work ordering/payback refinement tightened the medium
bounded-degree METIS-start class from 90% diagonal-present to 99%
diagonal-present. On the medium paper corpus this only affects the
`rajat03`/Rommes-BIPS/nopss structural class: `rajat03` remains a near-full
diagonal METIS start, while `bips98_606`, `bips98_1142`, `bips98_1450`, and
`nopss_11k` stay on AMD. The reactive static-match payback gate was also
extended to small many-block low-work cases with less than 5% weak and less
than 5% missing diagonal rows, preventing `bips98_1450` from replacing the
faster AMD factorization with a slower static-match candidate. A focused
five-row same-session check improved the geomean SPICE-cycle estimate from
about `0.0874s` under the previous auto policy to about `0.0564s`, while
keeping `rajat03` on METIS. The full 93-row medium paper run kept the same
three known failures (`bips07_1998`, `ss1`, and `mac_econ_fwd500`) and moved
KLS geomean from the previous saved `0.3096s` payback artifact to about
`0.3033s`; the CKTSO comparison ratio improved to about `1.149x` slower on the
90 common completed rows. Against the saved KLU2 artifact, KLS now wins 73 of
88 common completed rows with about a `2.01x` geomean speedup.

The EGraph floor was then lowered for a compact unscaled dominant-BTF shape:
95%+ largest-block coverage, 8-512 BTF blocks, a 10k-30k largest block, and at
least `2e7` measured factor flops. This is the smaller analogue of the retained
low-work dominant-BTF EGraph policy and intentionally excludes the TSOPF rows
with only two BTF blocks and the scaled `ckt11752_dc_1` shape. On the saved
medium artifact this selector matches only `coupled`. A focused three-pass run
with ten repeated refactors moved `coupled` from the prior `0.6345s` saved
SPICE-cycle estimate to a median `0.3783s`; the new row records 540 dependency
levels and about `1.19e7` dependency-work units. A one-pass full-medium run was
noise dominated overall (`0.3108s` geomean versus `0.3096s` in the previous
artifact), but the only row whose dependency schedule changed was `coupled`,
and that row improved from `0.6345s` to `0.3858s` in the full run. The retained
conclusion is narrow: compact dominant-BTF refactors can consume the exact
EGraph path, but this still does not address the much larger single-block and
ASIC refactor-kernel gap.

A scaled medium-dominant BTF EGraph gate was then retained for the adjacent
scaled shape that the unscaled compact gate deliberately skipped. The retained
selector requires BTF, 8-512 blocks, 95%+ largest-block coverage, a 30k-60k
largest block, active KLU scaling, at least `3e7` measured factor flops, at
least 1M LU entries, and at least `1.5e7` measured dependency-work units before
the consumer runs. In the current medium artifact this matches only
`ckt11752_dc_1`; a three-pass focused check built 1360 levels and reduced
repeated refactor from about `0.0100s` to about `0.0071s`, moving the focused
SPICE-cycle median from about `1.16s` to about `0.91s`. A full one-pass medium
run kept the same three known failures and moved geomean from about `0.3023s`
to about `0.3014s`; the CKTSO ratio improved to about `1.142x` slower and the
KLU2 comparison improved to about a `2.02x` geomean speedup. Nearby IBM
`dc*`/`trans*` scaled dominant-BTF guards keep their existing large-dominant
EGraph path.

The symbolic BTF retry policy was then extended to a large fragmented-BTF shape:
at least 100k rows, 8-512 BTF blocks, largest block below half the matrix, and
at least `5e7` estimated BTF flops. This is deliberately separate from the
dominant-BTF and inflated-many-block retries: it targets cases where BTF leaves
many off-block entries and the no-BTF symbolic score is at least 10% lower. In
the full medium artifact this changed only `hvdc2`, moving it from AMD/BTF to
AMD/no-BTF. The row's SPICE-cycle estimate dropped from about `4.88s` to about
`2.21s`, repeated refactor from about `0.0336s` to about `0.0154s`, and
off-diagonal pivots from 1019 to 8 with a valid residual. The full 93-row
medium run kept the same three known failures (`bips07_1998`, `ss1`, and
`mac_econ_fwd500`) and moved KLS geomean from about `0.3014s` to about
`0.2952s`; against saved CKTSO the ratio improved to about `1.118x` slower,
and against saved KLU2 the common-row geomean speedup improved to about
`2.06x`.

The BTF retry was then extended downward to a low-work medium many-block class:
4k-90k rows, at least 1024 BTF blocks, a 5-82% largest block, at most `2e7`
estimated BTF flops, and at least 95.5% structural diagonal coverage. The
no-BTF candidate is accepted only when symbolic fill stays within 1.75x and
estimated flops within 1.5x of the BTF symbolic. Focused checks showed forced
no-BTF wins for Bomhof `circuit_2/3/4`, Rommes BIPS/MIMO/NOPSS variants, and a
small `rajat22` improvement, while the diagonal-completeness gate kept the
noisy `rajat26/27` edge cases on the older BTF/static paths. The full 93-row
medium run again had the same three known failures and moved KLS geomean from
about `0.2952s` to about `0.2877s`. On the 90 rows common with saved CKTSO,
KLS is still about `1.09x` slower; on the 88 rows common with saved KLU2, KLS
is about `2.12x` faster. This is a retained general low-work BTF-overhead
reduction, not evidence that ordering alone closes the remaining CKTSO gap.

The largest-block ceiling was then widened from 80% to 82%, still under the
same low-work and diagonal-completeness guards. In the current medium corpus
this adds only `bips98_606`, another Rommes low-work matrix whose no-BTF
symbolic roughly halves the estimated flops without the tiny-block blow-up seen
on Sandia `mult_dcop_*` and `TSOPF_RS_b9_c6`. The full medium artifact again
kept the same three failures, moved KLS geomean to about `0.2853s`, narrowed
the CKTSO common-row ratio to about `1.081x`, and improved the KLU2 common-row
speedup to about `2.14x`.

The next retained policy pass addressed three low-work preprocessing costs
without changing the MC64 compatibility boundary. First, OPF-style bounded
degree matrices whose numeric diagonal is roughly half missing or weak now
start unscaled, and the score-gated single-block no-BTF retry is allowed down
to 12k rows for AMD/COLAMD auto candidates. This lets `OPF_3754` keep a
no-BTF/unscaled factor while still requiring symbolic evidence. Second, a
medium many-block, mostly diagonal, high-degree spike class whose largest BTF
block is 70-92% of the matrix now starts in KLU scale mode `0`; in the current
medium corpus the structural predicate matches only `rajat16`, `rajat17`,
`rajat18`, and `rajat26`, and it deliberately excludes the accepted-static
95%+ dominant Rajat rows. Third, the partial-weak pre-static matching gate now
requires at least five nonzeros per row on average, avoiding the expensive
rejected row-matching trial on sparse full-diagonal spike cases such as
`circuit_4` while preserving the denser Rajat static-match cases. The full
93-row medium run kept the same three failures, moved KLS geomean from about
`0.2853s` to about `0.2827s`, narrowed the CKTSO common-row ratio to about
`1.071x`, and improved the KLU2 common-row speedup to about `2.16x`. This is
still a preprocessing/payback refinement; the remaining large CKTSO gap is in
the KLS numeric refactor kernel and scheduling path.

The exact-EGraph all-pipeline mode was then extended from the retained
many-fringe dominant-BTF class to very high-work unscaled single-block
refactors. The selector is structural: one BTF block, no active KLU row
scaling, at least 100k rows, at least `1e9` measured factor flops, and at least
`1e9` measured no-pivot dependency-work before the EGraph consumer is allowed
to drop all cluster barriers. This targets the CKTSO/NICSLU observation that
pipeline mode can expose useful dependent-row overlap when the exact
dependency graph is large enough, while keeping lower-work unscaled HTC rows
and scaled `Raj1` on the existing barriered cluster/pipeline split. In focused
checks, `G2_circuit` switched from 501 cluster levels and 1617 pipeline-tail
columns to full-matrix all-pipeline execution, reducing repeated refactor from
about `0.447s` in the saved baseline artifact to about `0.392s`; `mc2depi`
also switched to full-matrix all-pipeline execution and moved from about
`2.16s` to about `2.10s` repeated refactor. A same-session full 93-row medium
comparison against a clean `HEAD` worktree kept the same known failures
(`bips07_1998`, `ss1`, and `mac_econ_fwd500`) and changed schedule metrics only
for `G2_circuit` and `mc2depi`; geomean moved from about `0.2869s` to about
`0.2859s`. The saved CKTSO comparison still shows KLS materially behind on
these rows, so this is a narrow scheduler improvement, not the missing
pivot-aware row-oriented numeric engine.

The same all-pipeline single-block selector was then widened to very high-work
row-scaled factors when KLU row scales and `Pnum` are already present, so KLS
can recompute `Rs`, execute the exact EGraph without cluster barriers, and
permute `Rs` back to pivot order after refactor. This keeps the same one-block,
100k-row, `1e9` factor-flop floor and does not affect lower-work scaled cases.
In the saved medium and large paper artifacts this structural gate matches only
`nxp1`. A same-session three-pass comparison against a clean baseline moved
`nxp1` median repeated refactor from about `0.677s` to about `0.649s`, median
fast-factor refactor from about `0.708s` to about `0.670s`, and the 100-step
SPICE-cycle estimate from about `75.7s` to about `73.4s`, while the solve time
stayed neutral.

Two follow-up scheduler probes were rejected after the scaled huge single-block
all-pipeline change. First, the all-pipeline work cursor was changed to claim
four-column dynamic chunks for modest average pipeline work. This reduced
atomic cursor traffic but delayed dependency publication inside each chunk and
was a major regression: on `rajat20`, `rajat25`, and `rajat28` repeated
refactor more than doubled, `rajat30` regressed by about 34%, and `nxp1` by
about 23%. Second, the fragmented non-dominant many-block EGraph shape used by
the ASIC 680k class was switched from the retained cluster/pipeline split to a
full all-pipeline run. That was neutral-to-worse: `ASIC_680k` was essentially
flat, `ASIC_680ks` regressed by about 4%, and the focused geomean regressed
slightly. The retained per-column dynamic cursor and barriered fragmented-BTF
split remain the better fit for the current fixed-pivot LU representation.

A third scheduler probe tried to skip `pipeline_done` publication for columns
with no successor in the retained U-pattern dependency graph. This was correct
but not useful on the current pipeline representation: a same-session focused
comparison against `23253eb` produced a `1.027x` repeated-refactor geomean
regression across `rajat20`, `rajat25`, `rajat28`, `rajat30`, `nxp1`,
`G2_circuit`, `mc2depi`, and `ASIC_680k`. The small wins on `rajat30`,
`G2_circuit`, and `mc2depi` did not offset `rajat25` at `1.088x`, `rajat28`
at `1.121x`, `nxp1` at `1.020x`, and `ASIC_680k` at `1.016x`. The extra
metadata load and changed publication pattern are therefore not a general
substitute for a larger CKTSO-style pivoting scheduler change.

A scaled refactor-map hot-loop probe was also rejected. The prototype stored the
original input row beside each retained refactor-map entry so scaled mapped
refactors could divide by `Rs[oldrow]` without chasing `row_idx[input_pos]`.
This looked like a cheap CKTSO-style row/segment metadata step, but the extra
map memory and load did not pay for itself on the current KLU storage. A
same-session focused comparison against `989360b` regressed repeated refactor by
`1.026x` geomean across `nxp1`, `rajat20`, `rajat28`, `Raj1`, `dc2`,
`G2_circuit`, `mc2depi`, `rajat25`, `rajat30`, and `ASIC_680k`; the scaled rows
alone regressed by `1.032x`, and the unscaled guards by `1.021x`. The prototype
was removed, leaving the retained refactor map limited to the map data that has
shown a general win.

The CKTSO paper was re-read after these scheduler probes because the remaining
gap is too large to explain by small EGraph bookkeeping. The missing mechanism
is larger and architectural: CKTSO's fast factorization is a row-oriented
up-looking kernel that first assumes the previous pivot order, schedules the
guessed EGraph in cluster/pipeline modes, checks the pivot against the current
row maximum, and, if the check fails, restarts only the ETree-descendant tail
with pivoting. KLS currently approximates only the no-pivot EGraph half on top
of KLU's column-oriented LU storage; the retained restart path invokes KLU's
serial pivoting kernel on the whole rejected BTF block. That is robust and
LGPL-compatible, but it cannot express CKTSO's partially recomputed
ETree-descendant tail or SubtreeLU's separator-tree private/pipeline queues.
The broad CKTSO gap on `G2_circuit`, `mc2depi`, `rajat20`, `rajat25`,
`rajat28`, `rajat30`, and `nxp1` should therefore be treated as evidence for a
new KLS-owned row/segment-oriented numeric layer, not another matrix-specific
ordering, scaling, MC64 import, or EGraph micro-optimization.

KLS now has a first KLS-owned scalar row-major no-pivot refactor scaffold behind
`KLS_ENABLE_ROW_REFACTOR=1` for unscaled single-block factors. It builds row
views of the fixed L pattern, U pattern, and factor-order input pattern, then
recomputes L row entries before writing each U row. Smoke checks on `add20` and
`G2_circuit` kept valid residuals, but the scalar row order is not a default
policy: on same-session one-thread checks, `add20` refactor was about
`0.00013s` versus the default `0.00006s`, and `G2_circuit` was about `1.19s`
versus `0.82s` over five repeated refactors. This narrows the missing paper
work: the row layer needs parallel scheduling and supernode/segment updates, not
just a scalar transposition of the KLU column kernel.

The row-major scaffold then gained an exact L-row dependency schedule and a
gated cluster-mode threaded execution path using the existing KLS worker pool.
This validates the CKTSO/SubtreeLU dependency direction without changing default
policy. Same-session four-thread checks kept valid residuals, but the
barriered scalar row tasks still lagged the retained column EGraph path:
`G2_circuit` measured about `0.86s` refactor with
`KLS_ENABLE_ROW_REFACTOR=1`, versus about `0.205s` for the default exact
U-pattern EGraph refactor. This points the next row-kernel work toward
supernode/segment updates and pipeline scheduling, not plain row-level cluster
barriers.

The row schedule was then upgraded to execute exact adjacent supernode-candidate
segments as single tasks, removing barriers between rows in the same dense-block
candidate. The grouped scheduler remained correct, but was still not a default
policy: on the same four-thread `G2_circuit` check, the grouped row path measured
about `0.88s` refactor, essentially no better than plain row-level barriers and
far behind the current column EGraph path. This confirms that merely coarsening
row tasks is not the missing paper mechanism; KLS needs a real segment kernel
that reuses shared trailing structure and eventually BLAS-style updates.

The grouped row path then gained a first true segment kernel: rows inside an
adjacent supernode-candidate segment use computed dense-block columns and one
shared trailing pattern instead of re-reading each row's full U pattern. This
kept valid residuals and moved `G2_circuit` from about `0.88s` to about `0.84s`
in the gated row path, but it remains far slower than the default column EGraph
path. The retained lesson is that the segment representation is now executable,
but the kernel still needs higher arithmetic intensity, such as dense triangular
mini-solves and batched trailing updates, before it can close the CKTSO gap.
The row-pattern builder now also retains the shared trailing `U` slice offset
for each executable multi-row group and the group dispatcher validates and
consumes that compact descriptor. This is deliberately a small SubtreeLU-style
row-supernode storage step, not a new tuning rule: the clear paper-level missing
piece for the slow cases remains a production row/up-looking factor/refactor
kernel with compact row-major supernode/segment updates and CKTSO's
ETree-descendant pivoting tail, rather than another ordering or scaling switch.

The row-segment path now exposes its executable workload separately from the
earlier adjacent-pattern candidate scan. `kls_stats`, `kls_bench` JSON/text, and
the gap decomposition script report row-refactor group count, group schedule
levels and maximum width, multirow segment count, covered segment rows, maximum
segment width, and dense/trailing segment work. These counters are only
populated when the gated row refactor pattern is built, so default benchmarks
continue to show zero. They let the next dense mini-solve or batched trailing
update experiment distinguish "no segment opportunity" from "segment
opportunity exists but the scalar kernel is still too weak." On `G2_circuit`
with four threads, the gated row path reported 113609 scheduled groups and 8035
multirow segments covering 44528 rows, max segment width 497, 606978 dense
entries, and 5033104 shared trailing entries. Its repeated refactor remained
about `0.906s` versus about `0.239s` for the default column EGraph path in the
same short run, so the counters expose real segment opportunity without
claiming the scalar row kernel is competitive yet.

The retained row-group pattern now also builds the reverse group dependency
graph, including unique group-to-group edges, root/leaf group counts, and
maximum group fanout. This is the row-segment analogue of the EGraph
root/leaf/fanout counters and is a direct bridge toward SubtreeLU-style
private/pipeline partitioning. It is intentionally metadata-only for now:
default KLS still uses the proven column/EGraph or KLU-backed refactor paths
unless the experimental row-refactor environment gates are enabled.

A first dense internal segment mini-solve was then added behind the same gated
row-refactor path. For exact multirow segments whose internal `L` pattern is a
dense lower block and whose dense/trailing work is large enough to amortize the
extra passes, KLS now gathers each row's external-update result into the
segment's own `L`/`U` value slots and applies the internal triangular updates
directly there. Smaller or non-dense-`L` segments keep the existing scalar
segment kernel. This is closer to the SubtreeLU supernode update shape than the
previous per-row scratch update, but it is still not default policy and not a
CKTSO-gap closer by itself: with four threads and five repeated refactors,
`G2_circuit` stayed valid but measured about `0.887s` for the gated row path
versus about `0.224s` for the default column EGraph path. A ten-refactor
`add20` smoke run also stayed valid at about `0.000236s`. The result indicates
that KLS now has an executable in-place dense segment solve scaffold, while the
larger missing piece is still compact row/segment numeric storage and batched
trailing updates with much higher arithmetic intensity.

The dense segment scaffold then stopped rediscovering segment shape on every
refactor. `kls_build_row_refactor_pattern` now precomputes dense-eligible group
flags, per-group trailing lengths, and the per-row internal-`L` offset used by
the gated mini-solve. `kls_stats`, `kls_bench`, and the gap decomposition script
also distinguish all executable row segments from the subset that can consume
the dense descriptor. This is a reusable row/segment descriptor step rather
than a CPU-specific hot-loop tweak. On a four-thread, five-refactor
`G2_circuit` check, KLS reported 8035 executable segments but only 662
dense-eligible segments covering 17070 rows; the gated row path stayed valid
and moved to about `0.846s`, while the default column EGraph path measured
about `0.201s`. On the `add20` ten-refactor smoke run, no segment crossed the
dense-work threshold and the path stayed valid at about `0.000223s`.

The dense segment mini-solve then split internal segment factorization from the
shared trailing-panel update. The dense path now first normalizes and applies
the internal lower/dense-upper block, checks pivots, and then runs the shared
trailing rows as a separate triangular update over the segment. This does not
change the flop count or make the gated row path competitive, but it matches
the paper direction more closely by isolating the panel operation that compact
row/segment storage or BLAS-style packing would later batch. Same-session
checks stayed valid: `G2_circuit` measured about `0.839s` repeated refactor for
the gated row path while the default column EGraph path measured about
`0.232s`, and `add20` stayed valid with no dense-eligible segments at about
`0.000274s`. The result is retained as row/segment scaffolding, not a
production dispatch candidate.

The separated trailing-panel step then gained worker-local compact storage.
Dense-eligible row segments now gather the shared trailing rows into a compact
row-major scratch panel, apply the triangular update there, and scatter the
finished panel back to the current KLU-owned value slots. This is still a
temporary bridge rather than persistent KLS-owned segment storage, but it
removes one layer of pointer chasing from the panel update and matches the
SubtreeLU supernode-storage direction more closely. Same-session checks stayed
valid: a clean `G2_circuit` gated row-refactor sample moved to about `0.798s`
repeated refactor with the same 662 dense-eligible segments covering 17070 rows,
while the default column EGraph path in a same-session check measured about
`0.223s`. The row path therefore remains gated, but this is the first retained
compact-panel step that measurably improves the experimental segment kernel.

The same compact scratch was then extended from the shared trailing panel to the
internal dense segment block. Dense-eligible row segments now gather the raw
internal lower block and dense upper block into one worker-local row-major
panel, perform the internal triangular update there, and scatter the normalized
`L`, dense `U`, and trailing `U` values back to the existing KLU slots only
after the segment is complete. This is still temporary scratch rather than a
persistent KLS numeric format, but it removes the main pointer-chasing loop from
the dense segment mini-solve. A clean four-thread `G2_circuit` sample stayed
valid and moved the gated row-refactor path to about `0.765s`; the default
column EGraph path in the same session measured about `0.205s`. This confirms
that compact segment storage is the right direction for the experimental row
kernel, while also confirming that the row kernel still needs much more work
before it can replace the current default.

The grouped row-refactor scheduler then gained work-balanced per-thread slices
inside each group level. The old grouped path assigned groups by static stride;
that can leave one thread with dense segment groups while others process mostly
singletons. KLS now builds persistent row-group level partitions from a
structural work estimate and reuses them across gated row refactors, mirroring
the work-balanced cluster slicing retained for the column EGraph path. This is
still a barriered level scheduler, not CKTSO's pivot-aware pipeline tail. It
stayed valid on the same focused checks: `G2_circuit` moved the gated row path
to about `0.739s` repeated refactor, while the same-session default column
EGraph path measured about `0.206s`. The small `add20` gated smoke remained
valid but measured about `0.000269s`, so this remains a large-row/segment
scaffolding step rather than a low-work dispatch policy.

The row refactor then gained a KLS-owned contiguous row-major `U` value mirror.
The retained KLU value slots are still updated so the existing triangular solve
path remains unchanged, but all gated row-refactor update loops now read
previous rows' `U` entries from KLS-owned row-major storage instead of following
one double pointer per value. This is a direct step toward the CKTSO/SubtreeLU
numeric format while preserving current semantics. Focused checks stayed valid:
the four-thread `G2_circuit` gated row path moved to about `0.532s` repeated
refactor, while the same-session default column EGraph path measured about
`0.234s`. The small `add20` gated smoke remained valid but moved from about
`0.000269s` to about `0.000291s`, confirming that this mirror is useful for
large row/update-heavy cases rather than as a low-work dispatch policy.

For patterns with dense-eligible row segments, the row-major `U` mirror then
became authoritative during the gated row refactor: KLS writes KLU's existing
`U` value slots only after the refactor succeeds. Sparse/no-dense row patterns
keep immediate KLU mirroring to avoid adding a final scatter pass where it is
not useful. This structurally gated policy preserves current solve semantics
while reducing pointer writes inside dense-segment hot loops. Focused checks
stayed valid: the four-thread `G2_circuit` gated row path measured about
`0.457s` repeated refactor, while longer no-dense checks remained stable
(`bcircuit` about `0.0050s`, `rajat22` about `0.00164s`, and `add20` about
`0.000209s` repeated refactor in same-session long-repeat samples).

The same dense-segment defer policy then gained a KLS-owned row-major `L` value
mirror. Sparse/no-dense patterns still write KLU `L` values immediately, but
dense-segment row refactors now write both `L` and `U` into KLS-owned mirrors
inside the hot loops and scatter both factor arrays back to KLU only after a
successful refactor. This keeps the current solve ABI while moving the
experimental row path another step toward native CKTSO/SubtreeLU-style numeric
storage. Focused checks stayed valid: the four-thread `G2_circuit` gated row
path measured about `0.441s` repeated refactor, and long no-dense samples
remained stable (`bcircuit` about `0.0046s`, `rajat22` about `0.00132s`, and
`add20` about `0.000175s` repeated refactor in same-session checks).

The dense-segment mini-solve then stopped using the worker-local dense panel as
the authoritative storage for dense-eligible row groups. For those groups, KLS
now writes raw internal `L` candidates, dense `U`, and trailing `U` values
directly into the persistent row-major mirrors, normalizes and updates those
mirrors in place, and leaves the older scratch-panel path as fallback for
non-deferred value storage. This removes one scratch-to-mirror copy layer while
preserving the final scatter to KLU for the current triangular solve ABI.
Focused checks stayed valid: four-thread `G2_circuit` measured about
`0.405s` repeated refactor in a two-repeat sample, down from the previous
retained `0.441s` row-mirror result. No-dense sanity checks still selected zero
dense segments and stayed numerically valid (`bcircuit`, `rajat22`, and
`add20`), while `mc2depi` exercised 1883 dense row segments with a valid
residual. This is a retained row/segment storage step, not the missing
CKTSO-style pivoting tail restart.

A selective dense-group scatter policy was then tested and rejected. The idea
was to let non-dense parallel row-refactor groups write KLU value slots
immediately while scattering only dense-group row mirrors after success. It was
structurally clean, but it moved the two intended dense-segment targets in the
wrong direction in same-session checks: `G2_circuit` regressed to about
`0.458s` repeated refactor from the retained native-mirror sample around
`0.405s`, and `mc2depi` regressed to about `2.39s` from about `2.19s`.
Keeping a uniform deferred scatter when dense row groups are present therefore
remains the better policy for the current hybrid KLS/KLU storage.

The native dense-group path was also tested with its repeated L/U shape checks
removed from the hot loop, trusting the precomputed dense-group descriptor in
the same way the older scratch-panel path does. This was also rejected:
same-session checks stayed numerically valid but moved `G2_circuit` to about
`0.445s` repeated refactor and `mc2depi` to about `2.28s`, both slower than
the retained native-mirror samples. The checks are therefore left in place until
KLS owns a fuller row/segment descriptor and numeric state that can make those
invariants cheaper to consume.

The experimental row/segment scheduler then gained a CKTSO-style dynamic
pipeline tail after barriered row-group cluster levels.
`kls_build_row_refactor_pattern` already stores levelized executable row
groups; the parallel row path now uses a `2 * threads` width split, processes
wide group levels with barriers, and lets workers claim the remaining
topological group sequence from an atomic cursor while waiting on predecessor
rows through the retained completion bitmap. Benchmark stats now report
row-refactor cluster levels, pipeline groups, pipeline rows, and pipeline work.
A constructed smoke matrix forces this path. A later instrumentation pass fixed
an important measurement ambiguity in this paragraph: enabling the checked row
fast-factor path does not by itself make the following unchecked `kls_refactor`
call use the row kernel. The row-pattern counters can therefore describe a
retained checked factor pass while the measured repeated refactor is still the
default column EGraph kernel. Benchmark JSON now reports
`row_refactor_last_run`, `row_refactor_last_checked`,
`row_refactor_last_parallel`, and row-refactor run counters, and `kls_bench`
accepts `--row-refactor env|off|refactor|checked|all` so future paper-suite
artifacts can separate checked-row factor probes from actual row-refactor
probes. The default KLS path is unchanged; this remains experimental
row-engine scaffolding until the row numeric factor/refactor path is validated
and faster across the broader paper matrix set.

`kls_bench` now has a deterministic diagonal-stress mode for measuring that
restart gap on real paper sparsity patterns without editing MatrixMarket files.
`--stress-diagonal-scale` and `--stress-diagonal-column` keep the first
factorization on the original values, then run repeated factor/refactor/solve
passes on values with selected diagonal entries scaled. With
`--repeat 1 --refactor-repeat 0`, the JSON `fast_rejected_*` fields can expose
the whole rejected BTF block, the suffix from the failed pivot, the exact
U-pattern descendant tail, the ordered-block ETree successor path, and the
refresh-state-seeded pivoting-tail scope that a CKTSO-style pivoting tail
restart would target. This is a diagnostic for architectural work, not a tuning
path for specific matrices.
The same diagnostic now also reports `fast_rejected_refresh_state`: unknown,
prefix-current, or all-current. This makes the benchmark artifact distinguish
KLS paths that can safely do block/tail continuation from parallel or KLU
all-refresh paths that must remain conservative until KLS owns the row/segment
numeric state needed for CKTSO-style ETree-descendant restart.

A follow-up removed redundant passes from that diagnostic path for unscaled
repairs whose surrounding numeric state is already current. After a failed
checked fast factorization, KLS tracks whether the failed pass refreshed all
columns, only a serial prefix, or stopped in an unknown parallel/EGraph state.
If all columns are current, or if the serial prefix reached the final BTF
block, KLS validates the repaired block and later block tail and returns
without immediately refactoring the same block again with the new fixed order;
unchanged earlier blocks are no longer rescanned on that shortcut path.
Same-session diagonal-stress checks kept valid residuals and reduced the
stressed `add20` factor path from the earlier `3.50s` to about `1.34s`, and
stressed `rajat03` from about `0.44s` to about `0.21s`. Other multi-block
repairs still retry the checked fast factorization because blocks after the
rejected block may not have been refreshed when the fast path stopped.

A block-level tail continuation was then retained for BTF cases. If the failed
checked pass stopped after a true prefix and the rejected block was not final,
KLS repairs the rejected block with pivoting and then refactors only later BTF
blocks with a checked continuation. The original unscaled serial continuation
was later opened to the threaded BTF worker pool, and the same input-row scale
invariant now lets KLU row-scaled repairs try that pool before falling back to
the serial checked scaled continuation. Same-session diagonal stress checks
against a clean `09d074a` baseline showed the intended non-final BTF benefit
with valid residuals: `coupled` moved from about `0.53-0.74s` to
`0.26-0.38s` factor time, and `circuit_1` moved from about
`0.012-0.014s` to `0.0096-0.0104s`. This is a useful CKTSO-aligned block-tail
step, but it still does not implement the larger missing single-block
ETree-descendant pivoting tail inside a KLS-owned row/segment numeric engine.

The experimental KLS-owned row refactor gained a serial checked fast
factorization pass. When `KLS_ENABLE_CHECKED_ROW_REFACTOR=1` explicitly enables
it, the pass checks the same no-pivot `L` multipliers used by the existing
column fast factorization, records the dependency pivot that fails, marks the
refreshed state as a serial prefix, and reuses the current block-restart
machinery. This is intentionally separate from `KLS_ENABLE_ROW_REFACTOR`, so
the existing unchecked repeated row-refactor path is unchanged. Focused checks
stayed valid: normal `add20` used the serial checked row path with no reject,
while stressed `add20`
(`--stress-diagonal-scale 1e-9`) reported a prefix-current fast reject at
the intended dependency pivot.

That checked row pass can now use the experimental parallel row/segment
scheduler when multiple threads are available. The threaded row tasks share the
same earliest-rejected-pivot recorder as the EGraph refactor path; dense and
non-dense row-segment kernels all check the same multiplier predicate. A
parallel checked-row rejection now reports
`KLS_FAST_REJECT_REFRESH_PREFIX` only when the retained completion bitmap proves
all rows before the rejected pivot have finished. If dense row segments deferred
their writes in KLS-owned row-major mirrors, KLS scatters only that proven
prefix into the KLU numeric object before reporting prefix-current. Otherwise
it remains
`KLS_FAST_REJECT_REFRESH_UNKNOWN`, because rows in the active level may finish
out of factor-order prefix before the stop flag is observed. Smoke coverage now
includes both a known weak-pivot 3-row prefix plus independent diagonal work to
force a two-thread row schedule, and a generated dense-row-segment case whose
prefix values change before a later rejected multiplier. The dense case checks
that the pivoting-tail restart consumes the published prefix and leaves a
bounded residual. This is a real CKTSO-style fast-factor-with-pivot-check
execution slice, but it still falls back to full block repair when the prefix
proof is unavailable rather than consuming the full ETree tail worklist.

The row-refactor pattern can now also retain the ordered-block ETree parent
array lazily for single-block row-major checked rejects. Rejected-pivot
diagnostics consume that cached parent after the first checked row-refactor
failure, instead of rebuilding the same CKTSO tail scope on later failures, while
the normal unchecked row-refactor path does not pay the setup cost. The smoke
test forces `KLS_ENABLE_CHECKED_ROW_REFACTOR` through a deterministic
pivot-reject repair so the retained row metadata, block restart, and ETree-tail
stats stay covered. This is a scheduler-metadata step toward the missing
pivoting tail restart, not a replacement for the row/segment-owned pivoting
factor kernel.

The checked row-refactor metadata now also retains reverse row dependencies
derived from the row-major `L` mirror. When a checked row refactor rejects a
pivot, KLS marks the exact row-successor tail, compacts it into increasing
row order, and reports the `fast_rejected_row_tail_*` diagnostics from that
retained topological list. KLS-owned pivot checks also record the factor row,
multiplier magnitude, accepted pivot magnitude, and candidate entry magnitude
that tripped the reject. For checked row-major rejects, KLS also scans the
retained row-tail list using the current prefix state and records the strongest
tail candidate row/value plus the number of tail rows that challenge the
current pivot tolerance. The diagnostic also records the candidate's retained
row-tail position and a `fast_rejected_tail_repair_ready` flag when the
checked row-major prefix is current and the best row-tail candidate satisfies
the same pivot-tolerance predicate that rejected the reused pivot. This gives a
row/segment pivoting tail kernel the local row candidate, tail location, and
value comparison that the older pivot-only diagnostics lacked.
This keeps the CKTSO restart target tied to KLS row storage rather than only
the KLU-column U-pattern or an ETree upper bound; the remaining missing piece
is the full CKTSO-style pipelined row/segment tail kernel.

The retained row-tail list remains an executable fallback seed for the
conservative pivoting-tail envelope. For prefix-current checked row-major
rejects with no saved interrupted worker seed, KLS can mark the retained
block-local row tail, close it through the ordered-block ETree, and record the
seed size as `fast_rejected_pivoting_tail_row_seed_columns`; otherwise the
CKTSO-style unfinished guessed-EGraph seed is preferred before row-tail and
suffix fallback plans. This is still not CKTSO's full concurrent pivoting-tail
executor, but it moves the checked row/segment metadata from diagnostics into
the restart planner.

The fallback pivoting block repair now also records the repaired row selected
at the rejected pivot, whether it matches the retained row-tail candidate, the
first pivot whose row changed in the repaired block, and how many changed
pivots occur before versus at/after the rejected pivot. These repair-delta
fields separate cases where the robust full-block KLU repair already preserves
the checked prefix and chooses the row-tail candidate from cases that would
need broader repivoting than a CKTSO-style local tail restart can safely
provide. KLS summarizes the compatible serial subset as
`fast_repaired_tail_restart_ready`, which now requires a prefix-current or
all-current non-root reject, zero changed pivots before the rejected pivot, a
changed suffix pivot, a topological pivoting-tail scope containing the reject,
and a validated live prefix replay. The checked row-tail candidate and
candidate match fields remain diagnostics; they are no longer required before
attempting the serial tail restart. When that flag is true, KLS also records
the fallback full-block repair work, the actual serial suffix column/work
count, and the saved-work estimate the local restart targets.
`scripts/summarize_tail_restart_opportunities.py` now turns those JSONL fields
into suite-level counts, blocker reasons, and largest saved-work opportunities
so the tail-kernel work can be prioritized from broad benchmark evidence. It
also totals and prints row-tail scope work when `fast_rejected_row_tail_*`
fields are present, and compares executed serial suffix-tail work against the
retained CKTSO-style pivoting-tail work to expose where the current fallback
overcomputes the planned restart set.
`scripts/run_bench_suite.py` now forwards the deterministic
`--stress-diagonal-scale` and `--stress-diagonal-column` controls to
`kls_bench`, so these tail-restart opportunity scans can be generated across
paper manifests instead of only from one-off piped benchmark rows.
The first focused stress probe shows why that distinction matters: stressed
`add20` chose the retained row-tail candidate at the rejected pivot, but the
full KLU block repair also changed one pivot before the rejected pivot, so a
future local tail kernel cannot treat tail-candidate availability alone as a
safe replacement for the fallback repair; its strict tail-restart-ready flag is
therefore false even though the full repaired-block work is about `1.66e6`.
The ordinary non-reject `G2_circuit` row-refactor path leaves these fields at
their sentinels.

EGraph checked-refactor rejects now distinguish a provably current prefix from
an unknown partial refresh when the retained pipeline-done generation is
available. KLS marks the reject as prefix-current only if every factor-order
column before the rejected pivot completed in the same EGraph generation; the
existing block-repair path can then continue with the serial BTF block tail
instead of discarding the whole fast pass. Focused diagonal-stress checks on
`coupled` show that the first observed EGraph reject is schedule-dependent:
when the dominant-block-start reject is observed after all earlier columns are
done, KLS records a prefix-current block repair and avoids the later full
checked-pass retry; otherwise it remains unknown and keeps the conservative
fallback. Rows such as stressed `onetone2`, where the completed-prefix proof
fails, remain classified as unknown.

The row-tail diagnostic was then extended beyond the single-block
row-refactor mirror. When that mirror is unavailable, KLS now scans the
rejected BTF block's numeric `L` columns, follows exact local row-successor
edges, and reports the affected tail work with the same per-column work
estimator used by the U-descendant and ETree-tail counters. This does not
change the numeric repair yet, but it exposes the row-tail scope on real
multi-block EGraph rejects. Repeated focused stress checks kept valid
residuals and reported nonzero row-tail scopes for `coupled` and `onetone2`,
where the previous JSON fields were zero.
The same diagnostic now preserves a checked-pivot row/value as a tail-candidate
row when that row lies in the recorded row-tail scope. The candidate is marked
repair-ready only for prefix-current rejects; unknown partial-refresh rejects
keep the conservative not-ready classification. Per-reject candidate and repair
fields are reset before each new reject so multi-restart stress rows do not mix
diagnostics from different attempts. The KLU-storage checked refactor kernels
now record the strongest violating entry in the rejected column instead of the
first entry visited by the sparse `L` pattern. This does not change the
accepted success-path factors or the fallback decision, but it makes the
CKTSO-tail diagnostic sharper: if the fallback pivot row still differs from
this strongest current-column candidate, then the missing paper mechanism is
broader pivoting-tail state, not merely an early-exit artifact of KLS's
diagnostic loop. The tail-opportunity summarizer now separates missing or
non-tolerance-valid candidates from unknown refresh state, all-current reject
state, and other not-ready cases so the paper-derived restart blocker is not
misreported as a pivot-search failure. A focused stressed probe on `coupled`,
`onetone2`, and `hvdc1` kept valid residuals; the two large repaired blockers
reported `tail_row == repair_row` and were instead classified as
`unknown_refresh_state`, confirming that the next missing CKTSO mechanism is
prefix-safe pivoting-tail execution rather than merely finding a local row
candidate.

The EGraph refactor worker scratch allocation was then narrowed for BTF paths:
single-block refactors still allocate one dense `n`-entry vector per worker,
but BTF EGraph workers only need block-local indices and now allocate
`maxblock` entries. This targets the fragmented-BTF/ASIC class where `n` can be
far larger than the largest diagonal block, reducing transient allocation and
zeroing without changing the numeric kernel or matrix-specific dispatch.
Focused post-change checks kept valid residuals; `ASIC_680k` refactor stayed
around `0.07-0.075s` and `rajat28` around `0.051-0.056s`, so this should be
treated as a memory-footprint cleanup rather than the missing CKTSO-scale timing
fix.

A small fragmented-BTF scheduler fix was then retained. The non-dominant
many-block EGraph path had hundreds of thousands of singleton BTF blocks at
level 0, but those singleton columns had zero scheduling weight because their
numeric work is outside the multi-column block LU estimates. The clustered
level partitioner could therefore hand a huge count of trivial-but-not-free
singleton columns to one worker before the substantial block work started. KLS
now assigns a unit balancing weight to singleton BTF columns only for the
fragmented non-dominant many-block shape. A same-session comparison against a
clean `d9f4c29` baseline moved `ASIC_680k` repeated refactor from
`0.0689-0.0695s` to `0.0569-0.0578s` with valid residuals. The gate is not
enabled for dominant-BTF all-pipeline schedules; `rajat28` kept the old
dependency-work and pipeline-work counters and stayed in the same
`0.0547-0.0566s` timing band. This is a general scheduling balance fix for the
fragmented BTF class, not another ordering, scaling, or MC64 import.

The fragmented-BTF EGraph path then stopped binary-searching the BTF boundary
array for every factor-order column. KLS already retains a fixed scatter map
for BTF refactors, so that map now also stores each column's owning BTF block
and the block-aware EGraph column kernel uses it directly. This targets the
same general ASIC-style shape with hundreds of thousands of BTF blocks, without
touching single-block EGraph rows. A five-row ASIC focus comparison against the
prior non-dominant many-block EGraph artifact improved geomean cycle time by
about 1.09x (`9.05s -> 8.31s`), with the largest win on `ASIC_680ks`
(`9.42s -> 6.79s`) and no >2% losses in that focus. Spot checks on `nxp1`,
`rajat28`, `coupled`, and `onetone2` stayed in their expected timing bands.
This is still a scheduler hot-loop cleanup, not the missing CKTSO row-oriented
pivoting factorization engine.

The initial-factor wrapper path then stopped recomputing KLU diagnostics
unconditionally after every optional policy probe. Auto-scale, METIS promotion,
and pivot-tolerance probes now report whether they actually replaced the
numeric object, and KLS recomputes `flops`/`rcond` only when the current numeric
state needs fresh diagnostics for the next decision or final stats. This removes
redundant O(U-pattern) flop scans and O(n) reciprocal-condition scans from the
first-factor path without changing ordering, scaling, or numeric acceptance.
Same-session focused checks were correspondingly modest and noisy:
`ASIC_680k` initial factor moved from about `1.72s` to `1.70s`, `rajat28` from
about `0.676s` to `0.668s`, while `rajat30` and `coupled` were neutral within
noise. This is retained as wrapper overhead cleanup; it does not replace the
larger row/segment numeric engine still needed for the CKTSO-scale gap.

The row-permuted solve wrapper then stopped allocating its dense permutation
scratch on every solve call. Static-pivoting cases reuse one solver-owned
workspace across repeated forward and transpose solves, and the smoke test now
calls the scaled pre-static forward and transpose solves twice to cover the
cached path. Same-session checks against prior artifacts showed forward-solve
improvements on the row-permuted focused rows (`onetone2` about
`0.00116s -> 0.00107s`, `hvdc1` about `0.000315s -> 0.000265s`,
`LeGresley_87936` about `0.00208s -> 0.00157s`, and `rajat23` about
`0.00172s -> 0.00157s`). TSOPF-style rows remained dominated by triangular
work and were neutral within noise. This is a retained repeated-solve overhead
cleanup, not CKTSO's missing structure-adaptive triangular solve.

The solve wrapper then stopped calling the full numeric-stats refresh after
every triangular solve. Solving does not change fill, FLOP, dependency, pivot,
or condition-estimate metadata, so KLS now refreshes only the last kernel
status and memory counters on the solve path. Focused checks kept residuals
valid and showed the expected small wrapper gain: `onetone2` forward solve
averaged about `0.00086s` in a ten-repeat run versus the previous
`0.00101-0.00107s` band, `rajat30` stayed in the retained KLU solve band at
about `0.0336s`, `nxp1` stayed around `0.0445s`, and `ASIC_680k` stayed around
`0.0098s`. This is a small solve-path cleanup, not the missing CKTSO
structure-adaptive solve.

Scaled fast-factor pivot rejection then gained the same block-local repair
path already used for unscaled checked refactors when the rejection came from a
KLS-owned checked path where `Rs` is still in input row order. KLS rebuilds
off-diagonal entries with the same row scaling convention used by KLU, runs the
KLU block factor kernel with the unpermuted row scale vector, then re-permutes
the scale vector after the new pivot order is installed. Later follow-up work
also lets scaled prefix-current repairs resume checked refactorization over
later BTF blocks through the threaded pool before serial fallback, and lets
validated non-root scaled prefix-current repairs use the serial suffix tail
restart inside the rejected block. Scaled KLU-refactor all-refresh rejects now
recompute row scales back to input-row order before trying the same repair path.
Default focused checks (`onetone2`, `rajat30`) stayed in the previous timing
band. This is a CKTSO-aligned coverage improvement for pivot-check restarts,
not a solution to the large `pre2`/`nxp1` gap.

KLS then exposed the fast-reject refresh-state classification in `kls_stats`
and benchmark JSON, and corrected the serial single-block mapped refactor to
record its pivot-check reject as prefix-current instead of unknown. Smoke
coverage now checks the prefix-current state for single-block, scaled
checked-refactor, and BTF block-restart repairs. This does not make KLS faster
by itself, but it turns the CKTSO tail-restart boundary into explicit benchmark
evidence instead of hidden control-flow knowledge.

A direct single-block row-major solve prototype was then tested and rejected.
The prototype built reusable row views of KLU's L and U factors with offsets
back into the existing LU storage, then used row dot-products for single-RHS
non-transpose solves. This matched the CKTSO paper's row-major triangular-solve
direction at a small scope, but it regressed the large single-block losses:
`rajat30` solve time rose from the retained KLU scatter-solve band of about
`0.034-0.036s` to about `0.099s` with three repeats and still about `0.071s`
with ten repeats, while `nxp1` rose from about `0.045s` to about `0.125s`.
The extra row metadata also added large memory traffic for 16-23M factor
entries. The result says the next triangular-solve attempt should not simply
transpose KLU columns into row-offset lists; it needs CKTSO's full
structure-adaptive partitioning with sparse-block/rectangular-slice decisions
and a storage layout designed for that access pattern.

The EGraph refactor path then stopped treating every repeated refactor as a
fresh thread launch. KLS now keeps a solver-owned EGraph worker pool plus dense
worker scratch vectors and dispatches each accepted EGraph refactor as a new
generation through that pool. Scratch is cleared only after early-abort paths
that can leave dense entries live, while successful no-pivot columns retain the
existing invariant that touched slots are zeroed as they go. This is still not
the missing CKTSO row-oriented pivoting-tail engine, but it removes repeated
thread launch and scratch allocation from the current exact-EGraph consumer.
A four-row focused probe with three passes and three refactors per sample
(`onetone2`, `rajat28`, `transient`, `power197k`) improved geomean cycle time
from the scratch-only artifact's `3.5166s` to `3.3196s`; compared with the
saved pre-pool blockmap artifact, `rajat28` moved from `6.3837s` to `5.3433s`
and `onetone2` from `1.8011s` to `1.7070s`, while the mostly initial-factor
`power197k` row stayed neutral. This is retained as a general repeated-refactor
overhead cleanup for EGraph-active SPICE rows.

The same retained EGraph state was then extended to the pipeline dependency
markers. Instead of allocating and zero-initializing an `n`-entry
`pipeline_done` array for every pipeline-capable refactor, KLS now keeps a
solver-owned atomic generation array and marks completed columns with the
current generation. Waiting workers compare against that generation, and the
array is cleared only if the unsigned generation counter wraps. This preserves
the existing dependency protocol while removing a repeated O(n) setup pass from
all-pipeline and pipeline-tail EGraph refactors. On the same four-row focused
probe, geomean cycle time improved from the worker-pool artifact's `3.3196s`
to `3.2480s`; `onetone2` moved from `1.7070s` to `1.5491s`, `transient`
nudged from `3.5853s` to `3.5510s`, and `power197k` stayed neutral because it
does not use the EGraph pipeline path in this run.

All-pipeline EGraph schedule setup was then narrowed by skipping the
per-level thread-slice table. Full all-pipeline runs set the cluster split to
zero and never enter the barriered cluster loop that consumes those slices, so
building the table only adds setup work and memory traffic. This is a small
setup cleanup, not a new numeric kernel: `G2_circuit` moved from the saved
generation-marker sample's `43.2928s` to `43.0747s`, while a three-pass
`rajat28` median was neutral within noise at `5.5145s` versus the prior
`5.4418s` artifact.

The huge single-block all-pipeline schedule then stopped storing and reading
the level-ordered column list. In that shape, factor order is already
topological for the no-pivot U-pattern dependencies, and the existing
generation-stamped wait markers still enforce readiness, so workers can claim
columns directly by factor-order index. The optimization is not used for
dominant-BTF all-pipeline schedules: a first broader trial made `rajat28`
noisier in the suite metric, and direct samples showed the BTF level-order path
should remain the structural default. With the narrowed rule, five direct
`rajat28` samples stayed in the same kernel band (`0.0446s` median refactor
versus the saved `0.0457s` no-thread-slices sample), while three `G2_circuit`
direct samples moved repeated refactor to about `0.218-0.241s` versus the saved
`0.398s` no-thread-slices sample.

The huge single-block unscaled EGraph column kernel then stopped repeating
structural checks that are already proven by the retained refactor map, exact
U-pattern schedule, and EGraph eligibility gate before worker launch. Numeric
singularity and pivot-threshold checks remain in the kernel. Three direct
`G2_circuit` samples moved median repeated refactor from the saved
natural-order sample's `0.2383s` to `0.2188s`, while direct `rajat28` guard
samples stayed in the same scaled dominant-BTF band at about `0.048s`. A
broader attempt to split the BTF value/scaling helper was rejected because it
pushed the scaled `rajat28` guard outside its timing band; that path stays on
the conservative helper logic.

A fresh four-thread medium paper-suite run was recorded after those EGraph
cleanups so future work uses current evidence instead of stale artifacts. With
`--repeat 1 --refactor-repeat 3 --timeout 120`, current KLS completed 90 of 93
medium rows, with `bips07_1998` singular and `ss1`/`mac_econ_fwd500` timing out
as before. Against the saved CKTSO medium artifact on the 90 common completed
rows, KLS is now `1.022x` slower geomean (`0.2698s` versus `0.2639s`), with 40
wins and 49 losses over 2%. Against the saved KLU2 artifact on 88 common rows,
KLS is `2.26x` faster geomean (`0.2385s` versus `0.5399s`) and wins 82 rows
over 2%. The remaining CKTSO losses are still dominated by repeated refactor
throughput: `rajat28`, `rajat25`, `ASIC_320k/320ks`, `rajat20`, `G2_circuit`,
`rajat03`, `ASIC_100ks`, `onetone2`, and `transient` lead the current gap.
Those rows, plus the timeout rows, are now tracked in
`bench/suitesparse_cktso_gap_manifest.txt` for focused regression checks before
rerunning the full medium suite.

The CKTSO gap review then found a concrete missing combination rather than a
generic ordering problem: the medium Rajat dominant-BTF rows benefited from
static row matching plus METIS nested-dissection, but the default fast-factor
path was accepting the static match with AMD. KLS now runs a narrowly gated
numeric METIS refinement for medium weak-diagonal static-match candidates whose
AMD symbolic leaves a many-block BTF with one dominant block. The METIS
candidate is kept only when actual factor fill/flops improve materially and
the reciprocal-condition estimate does not collapse. This changes only
`rajat20`, `rajat25`, and `rajat28` on the full medium run: they switch from
AMD/static to METIS/static, cut numeric fill by about 12-16%, and move their
SPICE-cycle estimates from `5.23s`, `5.10s`, and `5.33s` to `3.80s`,
`3.52s`, and `4.08s`. A 20-row CKTSO-gap guard improved by about 6% geomean
against the previous KLS artifact. The full medium suite remains essentially
flat within one-pass noise (`1.001x` versus the previous KLS artifact) and is
still `1.023x` slower than the saved CKTSO artifact on the 90 common completed
rows, while staying about `2.26x` faster than KLU2 on common rows. The remaining
large CKTSO gap is therefore still repeated refactor throughput on
ASIC/G2/Onetone-style rows, not a missing LGPL-compatible MC64 import or a
single broad ordering switch.

The CKTSO cluster/pipeline split was then rechecked because the CKTSO paper's
`alpha * thread_count` handoff is one of the few remaining explicit scheduling
knobs. KLS keeps CKTSO's `alpha=2` default, but now allows a structural
`alpha=4` handoff for compact large dominant-BTF schedules: unscaled matrices
with 128-512 BTF blocks, a 95%+ dominant block, max block at least 90k, and at
least `5e8` factor flops. This moves the pipeline tail earlier only for that
shape. A direct same-session five-pass A/B on the four touched medium rows
(`kls_asic_cluster_alpha2_ab_t4_p5_r10_timeout180.jsonl` versus
`kls_asic_cluster_alpha4_ab_t4_p5_r10_timeout180.jsonl`) improved geomean by
about `1.2%`: `ASIC_100k` improved `3.0%`, `ASIC_100ks` improved `1.7%`,
`ASIC_320k` improved `0.3%`, and `ASIC_320ks` regressed `0.3%`. Wider
`alpha=6`/`alpha=8` probes were mixed, and a full all-pipeline probe regressed
the same family, so this is retained only as a small structural scheduler
refinement. It does not explain the multi-x CKTSO refactor gap; that still
points to a row/segment-oriented numeric kernel, pivot-checked tail restart,
and eventually a SubtreeLU-style separator-tree/private-pipeline scheduler.

The EGraph refactor kernel then cached KLU's packed L/U column index/value
pointers at numeric-object scope. The previous worker hot loop recomputed the
same `LUbx + Lip/Uip + aligned-index-length` split for every U dependency and
for the output L column. The cache is rebuilt from current `LUbx` storage before
an EGraph refactor first runs, reused across repeated SPICE refactors, and freed
with numeric state or after a pivoting block restart that reallocates `LUbx`.
This is still KLU-storage based, not the row-major numeric engine CKTSO uses,
but it removes repeated wrapper-level pointer arithmetic from the exact-EGraph
no-pivot refactor path. A 10-row focused three-pass repeated-refactor probe
(`kls_lu_pointer_cache_probe_t4_p3_r10_timeout180.jsonl`) improved geomean by
about `2.5%` versus the saved hot-gap baseline and about `2.2%` versus the
current alpha=4 medium artifact on the same rows. The full 93-row medium run
(`kls_lu_pointer_cache_medium93_t4_r3_timeout120.jsonl`) kept the same three
known failures, improved the 27 EGraph-active rows by about `1.1%` geomean, and
was neutral-to-slightly-positive on all 90 completed rows (`0.998x` versus the
previous KLS artifact). The largest retained downside is `onetone2`, which was
about `8.5%` slower in the one-pass full run and about `2%` slower in the
focused run; an attempted structural gate for only unscaled high-work rows was
rejected because the added fallback branch path regressed the focused set
overall and lost the useful low-work `coupled` win.

The compact dominant-BTF EGraph gate was then extended downward for genuinely
small but still nontrivial refactor rows: unscaled matrices with 8-512 BTF
blocks, a 90%+ dominant block, max block in the 3k-10k range, and at least
`2e6` measured factor flops. The retained EGraph consumer still requires
measured dependency work before it runs, so this does not turn every tiny BTF
case into a threaded refactor. On the current medium artifact this structural
predicate matches only `rajat03` and `ACTIVSg2000`, both of which have zero
off-diagonal pivots and therefore are not MC64/matching failures. A five-pass
focused probe with ten refactors per sample
(`kls_small_compact_egraph_probe_t4_p5_r10.jsonl`) moved `rajat03` repeated
refactor from the saved `0.00127s` class to median `0.000696s`, and
`ACTIVSg2000` from `0.000821s` to `0.000379s`. Nearby guard rows such as
`gemat12` and `TSOPF_FS_b9_c1` did not build an EGraph schedule under this
gate. A one-pass CKTSO-gap focus run remained noisy on unchanged large rows,
so this is recorded as a narrow scheduler coverage improvement rather than a
claim that EGraph micro-gating closes the broad CKTSO gap. The MC64 boundary is
unchanged: existing MC64-style code may be reused only when it is
redistributable with LGPL KLS, such as the pinned BSD SPRAL scaling subset or a
verified compatible system library, not HSL MC64 or restricted solver-tree
copies.

An unscaled BTF-specific EGraph column kernel was then tested and rejected. The
prototype bypassed the generic scaling/value helper for unscaled BTF rows and
used the already validated refactor map directly. That looked aligned with
CKTSO's row/segment hot-loop emphasis, but the three-pass focused probe
(`kls_btf_unscaled_egraph_probe_t4_p3_r10_timeout180.jsonl`) regressed the
unscaled BTF EGraph set by about `2.5%` geomean. It helped some low-work IBM
and Rajat rows, but slowed the ASIC family by roughly 3-6% and hurt `coupled`
by about 20%, so the source change was removed. The result reinforces that
duplicating the current KLU-storage column kernel is not the missing
row-oriented CKTSO engine.

The single-block EGraph gate was then extended in the opposite direction:
low-work, unscaled, nearly no-pivot single-block rows. The retained rule is
structural: one BTF block, no KLU row scaling, 15k-250k rows, no more than
eight off-diagonal pivots, `4e6`-`5e7` measured factor flops, and at least
100k numeric LU entries. The EGraph consumer still requires measured
dependency work before it runs. On the current medium artifact this matches
`ACTIVSg10K`, `ACTIVSg70K`, `bcircuit`, `hvdc2`, and `OPF_10000`. A five-pass
focused probe with ten refactors per sample
(`kls_low_work_single_egraph_probe_t4_p5_r10_timeout180.jsonl`) improved every
touched row's repeated-refactor median versus the saved pointer-cache artifact:
`ACTIVSg10K` `0.00136s` to `0.00122s`, `ACTIVSg70K` `0.00338s` to `0.00219s`,
`bcircuit` `0.00414s` to `0.00304s`, `hvdc2` `0.01528s` to `0.01320s`, and
`OPF_10000` `0.00213s` to `0.00197s`. Guard rows that do not satisfy the new
predicate kept zero new EGraph dependency work. This is still a threshold
coverage improvement on the existing no-pivot EGraph refactor, not a
replacement for CKTSO's pivot-aware row-oriented factorization.

A narrow dense-fringe dominant-BTF scale policy was then retained for the IBM
`dc*` shape that still selected KLU max scaling despite having stable pivots
and valid unscaled residuals. The structural predicate is 100k-150k rows, 8-64
BTF blocks, a 99%+ dominant block, and 5-8 nonzeros per row; it only overrides
the scale when the existing value-based pattern policy would otherwise keep
max scaling (`2`). This last condition keeps the same-structure `trans4` and
`trans5` rows on their prior unscaled `-1` path. A five-pass scale-policy probe
with ten refactors per sample moved `dc1` from max scaling to scale `0` and
cut repeated refactor from `0.00810s` to `0.00629s`; `dc2` moved from
`0.00772s` to `0.00645s`; and `dc3` moved from `0.00836s` to `0.00648s`.
Relative residuals stayed around `1e-11` and off-diagonal pivots stayed zero
except for `dc1`, where the faster scale-0 candidate introduced five
off-diagonal pivots but remained numerically valid. Nearby scaled guards
(`ckt11752_dc_1`, `Raj1`, `rajat20`, and `rajat28`) did not change selected
scale. This is a scale-cost cleanup for a repeated-refactor class, not an MC64
import or a substitute for the missing CKTSO row-oriented kernel.

The low-work dominant-BTF EGraph gate was then probed below the older 90k
largest-block floor. A broad first version for 80%+ dominant blocks also
activated on `ckt11752_tr_0` and regressed that guard, so it was narrowed to
the high-coverage no-pivot subcase: unscaled BTF, 30k-120k rows, at most 5000
BTF blocks, a 95%+ dominant block in the 60k-90k range, zero off-diagonal
pivots, and `1e7`-`3e7` measured factor flops. In the current medium evidence
this matches only `LeGresley_87936`. A seven-pass focused probe with ten
refactors per sample
(`kls_low_work_high_coverage_btf_probe_t4_p7_r10_timeout180.jsonl`) moved
`LeGresley_87936` repeated refactor from the current focused `0.00647s` median
to `0.00464s`, while `ckt11752_tr_0` kept zero EGraph dependency work under
the narrowed rule. This is another measured-work scheduler coverage step for
the existing no-pivot EGraph refactor; it does not change the conclusion that
the remaining ASIC/G2/Rajat gap needs a different row-oriented numeric kernel.

The MC64 source boundary was clarified again after accepting that existing
implementations can be reused when they are license-compatible with LGPL KLS.
The rule is license-based rather than authorship-based. The pinned SPRAL
submodule remains the retained vendored source because its BSD-3-Clause license
permits source and binary redistribution with preserved notices, while HSL's
current no-cost licence is personal-use only and does not allow redistribution
in source or binary form. A small BSD `mc64` translation of part of SPRAL is
compatible as a reference, but it does not improve integration over the already
pinned SPRAL C/Fortran interface. The current CMake guard therefore remains
correct: bundled SPRAL checks for redistribution-compatible BSD text, and a
system SPRAL/MC64-style library requires the builder to explicitly set
`KLS_SYSTEM_SPRAL_LGPL_COMPATIBLE=ON` after verifying the selected library.
This keeps MC64-quality preprocessing on the roadmap without admitting HSL
MC64 or restricted solver-tree copies into the LGPL distribution.

The latest discarded prototypes further support the larger-engine diagnosis.
A serial single-block mapped refactor variant reused cached `L`/`U` pointer
indices inside `kls_single_block_mapped_refactor`. It looked attractive as a
row/segment hot-loop cleanup, but same-session A/B checks regressed key
low-work rows (`OPF_3754` about `1.05x`, `xingo_afonso_itaipu` about `1.11x`,
and `ww_vref_6405` about `1.08x` versus the committed baseline), so the patch
was removed. An EGraph wait-loop CPU-relax prototype was also rejected. It
helped some large ASIC samples in one form but regressed `G2_circuit` and
`onetone2`, and narrowing the gate still left mixed results. These are
bookkeeping or spin-wait effects, not the missing CKTSO-scale mechanism.

A current four-thread medium-paper decomposition refresh now shows KLS only
about `1.03x` slower than the saved CKTSO artifact in geomean on the 90 common
completed rows, but the largest losses remain multi-x: `ASIC_320k`,
`ASIC_320ks`, `onetone2`, `ASIC_100ks`, and `G2_circuit` are all above
`2.6x` SPICE-cycle ratio and dominated by repeated refactor work. The refreshed
`bench/suitesparse_cktso_gap_manifest.txt` now tracks the top 40 current
losses plus the timeout rows `mac_econ_fwd500` and `ss1`, so future changes
are checked against the matrices that still expose the large gap instead of
only the older pre-EGraph focus rows.

The BTF EGraph numeric path was then given the same branch-light unscaled
column specialization that single-block EGraph already used. The retained
dispatch is structural: it only applies to unscaled BTF refactors whose largest
block is at least 30k columns, leaving the small compact BTF path on the older
generic kernel after an ungated probe regressed `coupled`. In a same-session
large-BTF A/B with three passes and ten refactors per sample, the retained
candidate improved the nine-row focused geomean by about `1.2%` versus the
committed baseline. Repeated refactor improved on the intended large-block
rows: `ASIC_320k` `0.1033s` to `0.0990s`, `ASIC_320ks` `0.0852s` to
`0.0813s`, `trans5` `0.00665s` to `0.00624s`, and `LeGresley_87936`
`0.00458s` to `0.00420s`. The important unchanged-path guard set
(`G2_circuit`, `onetone1`, `onetone2`, `rajat28`) stayed neutral-to-slightly
positive in a two-pass comparison. This is a real hot-loop cleanup for the
existing EGraph BTF refactor, but it is still an incremental KLU-storage
optimization rather than CKTSO's missing row-oriented pivoting-tail engine.

That BTF-specialized kernel was tightened once more by removing validation
branches already guaranteed by the EGraph dispatcher, refactor map builder, and
LU-pointer-cache builder. The hot path still checks Offx bounds while writing
off-block values, but trusts the already validated block metadata and
topological U pattern, matching the single-block EGraph kernel's trust model.
Against the previous retained specialization, the same nine-row BTF focus
improved by about `1.8%` geomean with no failed rows. The refactor medians
improved on the intended large-block cases: `ASIC_100k` `0.04441s` to
`0.04317s`, `ASIC_100ks` `0.04496s` to `0.04354s`, `ASIC_320k` `0.09904s`
to `0.09701s`, `ASIC_320ks` `0.08130s` to `0.07976s`, and `ASIC_680ks`
`0.04886s` to `0.04740s`. A broader attempt to hoist the
`shared->check_pivots` branch out of every L-update loop was rejected: it
helped some ASIC samples but regressed low-work/scaled guards such as
`LeGresley_87936`, `onetone2`, and `rajat28` under repeated same-session
probes.

Two follow-up BTF bookkeeping probes were also rejected. First, a cached
block-local row map avoided `global_row - block_start` inside the BTF EGraph
kernel, but the extra map memory did not pay back: the nine-row BTF focus
regressed by about `0.5%` geomean versus the retained validation-trim commit,
with clear losses on `LeGresley_87936` and `ckt11752_tr_0`. Second, removing
the remaining per-column LU-pointer-cache validation from the BTF specialized
kernel also regressed the same focus set by about `0.7%` geomean and slowed the
large ASIC rows. The retained BTF specialization therefore keeps the pointer
cache guard and computes block-local row indices directly. The broader signal
is unchanged: small KLU-storage bookkeeping trims are now yielding mixed or
single-digit effects, while the CKTSO gap still requires a row/segment numeric
engine with pivoting-tail restart semantics.

A corresponding single-block EGraph pointer-cache validation trim was also
rejected. Removing the per-column `refactor_lu_pointer_count`/pointer guard
helped `G2_circuit` slightly in a same-session two-pass focus run, but
regressed the Onetone rows by about `10-11%` in SPICE-cycle time and moved the
four-row guard geomean about `4.9%` slower. This confirms the remaining
single-block gap is not a simple pointer-cache guard branch; the retained
single-block kernel should keep its current validation shape until KLS owns a
different row/segment numeric representation.

The MC64 import policy was converted into an executable repository audit after
accepting that existing MC64-style code can be used when it is LGPL-compatible.
`scripts/audit_license_boundary.py` now checks that KLS remains LGPL, that the
only bundled MC64-adjacent implementation used by the build is SPRAL's
BSD-licensed scaling subset, and that first-party/vendor source outside the
vetted SPRAL path does not contain restricted HSL/MC64 markers. This does not
close a performance gap, but it prevents future matching/scaling work from
accidentally crossing the licensing boundary while keeping compatible existing
code available for KLS.

The refactor schedule diagnostics now also record dependency root columns,
leaf columns, and maximum successor fanout. These counters are computed while
KLS already walks the numeric U pattern to build exact EGraph levels. They are
intended to quantify the opportunity for a successor-driven ready scheduler or
row/segment engine on the CKTSO-gap rows, rather than to tune on matrix names
or infer the gap from total edge count alone.

A bounded successor-ready queue was then prototyped for modest EGraph pipeline
tails by retaining reverse U-pattern edges and running eligible pipeline
columns only after their pipeline predecessors completed. The design was
rejected. Same-session A/B against commit `ba8758c` with three passes and ten
refactors per pass showed broad regressions: `onetone1` cycle median
`6.06s -> 18.97s`, `onetone2` `1.60s -> 3.52s`, `ASIC_320k`
`13.62s -> 29.74s`, and `ASIC_320ks` `11.43s -> 25.02s`. Refactor medians
regressed similarly (`ASIC_320k` `0.097s -> 0.261s`, `ASIC_320ks`
`0.079s -> 0.216s`). This means the missing CKTSO mechanism is not a simple
reverse-edge ready queue layered on KLU column storage; KLS needs a different
row/segment numeric representation or separator/private-pipeline engine where
successor scheduling does not add another high-overhead synchronization layer.

The opposite scheduling probe was also rejected: disabling partial EGraph
pipeline tails and forcing the remaining exact-EGraph tail through level
barriers while preserving all-pipeline shapes. Against commit `eac177b`, the
intended ASIC rows regressed sharply in same-session tests with three passes
and ten refactors: `ASIC_320k` cycle median `13.58s -> 30.24s` and refactor
median `0.0976s -> 0.265s`; `ASIC_320ks` cycle median `11.43s -> 25.75s`
and refactor median `0.0797s -> 0.224s`. This confirms that KLS's current
pipeline tail is necessary on these rows, even though it remains much slower
than CKTSO. The remaining gap is therefore not fixed by choosing between
coarse level barriers and KLS's current fetch-and-wait tail; it points back to
the numeric representation and update granularity inside the heavy tail
columns.

A simpler scheduler-overhead probe then tried leasing eight consecutive EGraph
pipeline positions per atomic claim. This was also rejected. The pipeline column
array is not a strict dependency topological order inside every small range; a
worker can lease a range, start an earlier position, and then wait for a
dependency that sits later in the same leased range and therefore cannot be
published by another worker. The focused KLU-first run
(`build/kls_egraph_chunk8_klu_first_gap5_t4_r3_timeout120.jsonl`) exposed the
problem immediately: `ASIC_320k` grew to `30.05s`, `ASIC_320ks` to `24.65s`,
`ASIC_100ks` to `12.84s`, `G2_circuit` to `69.87s`, and `onetone2` timed out
at the 120 s harness limit. The code was removed. Any future batching of the
EGraph tail has to lease dependency-ready sets or provably independent chunks,
not just contiguous positions from the current fetch-and-wait worklist.

All-pipeline single-block scheduling was also probed by replacing natural
column fetch order with exact EGraph level order for huge single-block
refactors. This was meant to expose more independent G2/mc2depi work without
adding a reverse-edge ready queue. It was rejected immediately on the intended
guard: with three passes and five refactors, `G2_circuit` regressed from a
`25.39s` cycle median and `0.220s` refactor median to `43.57s` and `0.404s`.
The natural all-pipeline order is therefore retained for huge single-block
cases. The result further narrows the scheduler diagnosis: changing the order
in which KLS feeds KLU-column kernels is not enough; the missing improvement
needs different per-column update granularity or storage.

The EGraph diagnostics were therefore extended again to report maximum
per-column work and maximum pipeline-tail column work. These fields quantify
when the outer schedule has exposed enough independent columns but one or a
few KLU-column updates still dominate the tail. They are intended to guide the
row/segment numeric-engine work rather than another column-ordering probe.

The immediate hot-loop follow-up was to centralize the repeated sparse
scatter update `x[row] -= L(row,j) * U(j,k)` and unroll it four ways in a
single inline helper used by the serial mapped, BTF, and EGraph refactor
paths. This does not change ordering, matching, pivot policy, or scheduler
semantics; it only reduces overhead in the per-column update loop exposed by
the new maximum-column-work diagnostics. Same-session A/B against commit
`d213eb8` retained the helper: on six heavy CKTSO-gap rows with three passes
and ten refactors per pass, SPICE-cycle geomean improved from `7.77s` to
`7.48s` with no row slower; refactor medians improved on `G2_circuit`,
`onetone1`, `onetone2`, `ASIC_100ks`, `ASIC_320k`, and `ASIC_320ks`. A
one-pass top-20 CKTSO-gap guard also completed all rows and improved geomean
from `3.64s` to `3.55s`; a follow-up three-pass loss check showed the
apparent `rajat28`, `onetone2`, and `ASIC_320ks` losses were noise or reversed
with longer repeats, while tiny `gemat12` remained cycle-noisy despite a
faster measured refactor. This is a small retained KLU-storage kernel cleanup,
not the missing CKTSO row/segment engine.

A wider scatter-unroll follow-up was rejected. Replacing the retained four-way
helper with an eight-way loop, and then with an eight-way-only-for-long-segments
variant, produced mixed same-session results: the six-row heavy focus was only
about `0.6%` positive, the top-20 CKTSO-gap guard was only about `0.3%`
positive, and repeated checks still showed row-level losses such as `rajat28`
and `ASIC_320ks` in refactor time. More importantly, this direction is
CPU-code-generation-specific rather than a solver algorithm improvement. KLS
therefore keeps the simpler four-way scatter helper and leaves further
progress to row/segment numeric storage, pivoting-tail restart, and general
cost-model work.

A branch-light scaled BTF EGraph column specialization was also rejected. The
prototype mirrored the retained unscaled BTF specialization but divided mapped
entries directly by `Rs[oldrow]`, bypassing the generic scaled value helper for
scaled multi-block EGraph refactors with large blocks. It was structurally
sound, built cleanly, and passed the smoke/license tests, but same-session
three-pass checks showed the wrong tradeoff: the scaled-focused eight-row
geomean was only about `0.9%` positive, while the hard scaled Rajat rows
regressed (`rajat20` cycle `3.55s -> 3.59s`, `rajat28` `3.54s -> 3.57s`).
The small `ckt11752_dc_1` win was not enough to justify a duplicate scaled BTF
kernel that worsens the remaining CKTSO-gap rows. KLS therefore keeps the
generic scaled EGraph value path until the larger row/segment numeric engine
can improve scaled and unscaled hard rows together.

A fallback row-segment compact-scratch probe was also rejected. Dense-eligible
segments already keep their compact dense/trailing panel path, but the broader
prototype packed every executable exact-U-pattern multirow segment into
worker-local row-major scratch so later rows could read internal and trailing
U values without KLU double-pointer indirection. It built cleanly and stayed
valid, but the intended four-thread `G2_circuit` guard moved the gated row
path back to about `0.797s` repeated refactor, versus about `0.739s` before
the probe and about `0.233s` for the default column EGraph path in the same
session. The result reinforces the current diagnosis: shallow packing around
KLU storage is not the CKTSO/SubtreeLU mechanism; KLS needs persistent
row/segment numeric storage and scheduler semantics that make segment updates
the native operation rather than a copied side path.

The existing block-restart fallback was then tightened without changing its
numerical semantics. `kls_pivot_restart_rejected_block` previously allocated a
temporary `n`-entry symbolic inverse row map and made multiple full passes over
it before calling the KLU block kernel. It now reuses `numeric->Pinv` as
restart scratch, fills it once from the symbolic row permutation, and rebuilds
the accepted numeric inverse after the repaired block permutation is published.
Focused stressed restart checks on `coupled`, `onetone2`, and `hvdc1` kept
valid residuals and the same restart diagnostics. This reduces current
full-block fallback overhead, but it is still a cleanup around KLU's block
repair path rather than the missing CKTSO pivoting-tail factorization.

The tail-restart diagnostic was then corrected to retain an executable
pivoting-tail plan, not only the old single ETree successor path. For
prefix-current rejects, KLS now seeds the `fast_rejected_pivoting_tail_*` plan
from the checked-refactor unfinished-node bitmap when it is available, and
falls back to the rejected block suffix otherwise; all-current or
unknown-refresh rejects seed the plan from the rejected pivot. The strict
tail-restart saved-work estimate uses this pivoting-tail scope. A focused
`coupled` diagonal-stress probe reported a prefix-current reject at pivot `24`,
a suffix/pivoting-tail scope of `11293` columns and about `1.193e7` work,
versus the older single ETree path of `8194` columns and about `1.191e7` work;
the full repaired-block work remained about `3.650e8`. This keeps the
opportunity estimate tied to the CKTSO restart set and gives the future
row/segment numeric kernel a concrete ordered worklist to consume.

That pivoting-tail plan is now checked as a real scheduler contract instead of
only a count/work estimate. KLS records the first and last global rows in the
retained plan, whether the plan includes the rejected pivot, and whether the
stored order is topologically safe with respect to the ordered-block ETree
parent links. The stressed `coupled` case is still classified as a
root-of-block reject with no reusable serial prefix, but its retained tail plan
is now verified as the ordered worklist a CKTSO-style row/segment tail kernel
would need to consume.

The strict tail-restart readiness gate was then strengthened to compare the
fallback repair against the old block pivot order and to replay KLU's live
prefix bookkeeping through the rejected pivot. The validator rebuilds local
final-to-live row maps, replays KLU `P`/`Pinv` pivot logging, unfinalizes stored
`L` row indices back to local row coordinates, and reconstructs prefix
`Lpend` pruning boundaries from the already-pruned KLU columns. A new smoke
case covers a strict-ready nonzero-prefix repaired block. The same focused
stressed `coupled` probe still reported one strict-ready reject, with about
`3.531e8` saved work over full-block repair. This removes another entry-state
ambiguity before broadening the serial pivoting tail path and before building
the full CKTSO-style pipelined row/segment tail kernel.

The executable serial pivoting-tail restart was then broadened to use the
validated prefix-current proof directly instead of requiring the checked
row-refactor tail-candidate diagnostic. That means default mapped checked
refactors can now run the local tail restart when the rejected pivot is not the
block root, the old prefix is replayable, and the retained pivoting-tail scope
is topological. A new smoke case keeps the default mapped-refactor distinction
explicit: with the checked row-refactor environment disabled, the block repair
uses one serial tail restart, preserves the prefix, and produces a bounded
solve. The benchmark summarizer now reports executed tail restarts separately
from remaining strict-ready opportunities and no longer classifies a missing
row-tail diagnostic as the first hard blocker.

The same serial suffix restart was then opened to unscaled all-current
postcheck rejects. Those failures have already refreshed the block with the
old pivot order, but the prefix before a non-root rejected pivot is still
replayable by the same live-prefix validator. KLS now tries the suffix restart
before full block repair in that state too. The `fast_repaired_tail_restart_*`
work counters now use the actual serial suffix work KLS executes, while the
`fast_rejected_pivoting_tail_*` fields remain the CKTSO-style ETree-tail
diagnostic target for a future row/segment scheduler.

The serial tail retry then stopped copying the full off-block row/value arrays
into scratch storage. The trial factor only needs mutable `Offp` offsets while
constructing columns; after a repair is accepted, KLS rebuilds `Offi`/`Offx`
from the final `Pinv` anyway. This removes a full `nzoff` copy from every
serial tail attempt without changing pivot choices or the published numeric
state.

Accepted serial tail retries then stopped doing a full global off-diagonal
rebuild when the already-refreshed prefix can be preserved. KLS now rebuilds the
off-diagonal suffix from the rejected column through the end of the matrix from
the rebuilt `Pinv`, preserving only the proven-good prefix `Offp` count and
falling back to the full rebuild if the suffix shape proof fails. The BTF
off-block tail smoke changes suffix off-block values between the base and
repaired matrices, so stale values in that suffix are observable in the solve.

The KLS-owned pivoted block restart was then opened to root-of-block rejects.
The serial suffix path still requires a reusable prefix and therefore does not
count root rejects as `fast_tail_restarts`, but `local_reject == 0` can now run
the same sparse pivoting kernel over the whole rejected block before considering
the KLU block-kernel fallback. `fast_kls_block_restarts` records these KLS-owned
block repairs separately from total block restarts, and the prior-pivot smoke
case now requires a root reject to use that path. This is a direct step toward a
KLS-owned factorization kernel, but not the CKTSO ETree-descendant pipelined
tail executor.

KLS now also gives the independent row-up-looking factorization kernel a chance
when a checked reject cannot be recovered by the exact in-place tail repair.
This moves that recovery branch toward the papers' "switch from invalid checked
refactorization to pivoting factorization" rule: exact block/tail repair remains
first, then KLS rebuilds the rejected BTF block with its row-first
dynamic-column-pivot executor, then it tries the quality-checked KLS-first
whole-numeric rebuild before the KLU block-kernel fallback. The block-local
executor now covers both unscaled and KLU row-scaled repair states by consuming
input-row `Rs` during row construction and permuting it back to pivot order
after an accepted scaled block. The scaled preserved-column refresh now uses
the same convention for masked ETree-tail attempts before this broader rebuild
is considered. The whole-numeric checked-reject recovery now covers unscaled
and KLU row-scaled states; it bypasses the normal automatic KLS-first cost gate
but still honors an explicit
`KLS_ENABLE_KLS_FIRST_FACTOR=0` disable. Successful row-first block repairs are
reported as `fast_kls_block_restarts`; successful whole-numeric recoveries are
reported as `fast_kls_rebuild_restarts`. This still stops short of CKTSO
Algorithm 5's ETree-descendant pipelined tail scheduler in place, but removes
another KLU-kernel step from the checked-reject fallback ladder.

The row-refactor benchmark controls were then made explicit after the
checked-row/refactor ambiguity above was found. `kls_bench` now accepts
`--row-refactor env|off|refactor|checked|all`, `run_bench_suite.py` forwards
the same option, and stats report both retained row-pattern shape and actual
last numeric row-kernel execution. A focused four-thread `G2_circuit` check
confirmed the distinction: `--row-refactor checked` ran one checked row
fast-factor pass but reported `row_refactor_last_run=0` after the subsequent
unchecked refactor, while `--row-refactor refactor` reported
`row_refactor_last_run=1` and `row_refactor_last_checked=0` for the final
refactor. In that same noisy sample the explicit row refactor was slower than
the default EGraph refactor, so this is retained as measurement cleanup and
benchmark reproducibility, not as a default policy change.

The EGraph refactor worker then stopped re-running the unscaled single-block
and large-BTF kernel dispatch inside every column. KLS now records the proven
column-kernel kind in the shared worker state before launching the solver-owned
pool, and the hot cluster/pipeline loops call the selected kernel directly. A
same-session three-pass A/B guard with five refactors per sample improved
median repeated refactor on `G2_circuit` by about `1.4%`, `ASIC_320k` by about
`1.2%`, and `onetone2` by about `6.0%`. This is retained as branch cleanup for
the existing exact-EGraph refactor consumer, not as evidence that dispatch
cleanup can close the CKTSO-scale row/segment gap.

The experimental parallel row-refactor pattern then started retaining a group
execution kind: single-row, generic multi-row, or dense multi-row. The row
scheduler consumes that persistent metadata directly and reuses each group's
precomputed trailing length, so it no longer probes the dense-group path for
groups that were already classified as generic. Same-session A/B checks of the
explicit `--row-refactor refactor` path showed `G2_circuit` essentially neutral
to improved at `0.981x` new/base median repeated-refactor time, `OPF_10000`
neutral at `1.000x`, and `xingo_afonso_itaipu` improved to `0.838x`. This is
retained as row-engine groundwork for the future CKTSO-style row/segment
numeric kernel; the row path remains explicit or environment-gated and is still
not a default solver policy.

The next row-scheduler metadata step retained each row group's external
dependency rows. During the dynamic pipeline tail, generic-only row patterns
can now wait once on those external rows before entering the group kernel,
instead of testing predecessor completion inside every row update. The broader
version that also pre-waited dense row-segment groups was rejected because it
regressed `G2_circuit`; dense segments keep their old in-kernel wait behavior.
With the generic-only gate, the same-session Release A/B on the explicit
`--row-refactor refactor` path improved `G2_circuit` to `0.969x` new/base
median repeated-refactor time, kept `OPF_10000` essentially neutral at
`1.003x`, and improved `xingo_afonso_itaipu` to `0.929x`. This is retained as
CKTSO-style scheduler metadata for the experimental row path, not as a claim
that the default KLU-storage solver has closed the CKTSO row-engine gap.

Two follow-up row-group scheduler metadata shortcuts were rejected. First,
retaining a precomputed per-group work estimate avoided repeated scans in
pipeline-scope reporting and level-slice construction, but same-session Release
A/B on the explicit row-refactor path regressed `OPF_10000` repeated refactor
to `1.157x` new/base and `xingo_afonso_itaipu` to `1.076x`, despite a small
`G2_circuit` win at `0.975x`. Second, retaining the generic-only prewait
eligibility byte avoided recomputing the prewait predicate in the pipeline
loop, but regressed `OPF_10000` to `1.110x` and `xingo_afonso_itaipu` to
`1.067x`, while `G2_circuit` was neutral at `1.001x`. Both probes are too
shallow: the remaining CKTSO gap still points to native row/segment numeric
storage and pivoting-tail restart semantics, not more cached scalar scheduler
decisions around the current KLU-backed row experiment.

The retained row-group reverse graph is now consumed by the experimental row
pipeline tail as a bounded successor-ready queue. Queue setup marks the groups
remaining after the cluster split, counts only tail-local predecessor groups,
starts from zero-predecessor tail roots, and enqueues successor groups when
their tail predecessor count reaches zero. A strengthened smoke fixture uses
three dense row groups with a lower coupling, proving a nonzero group edge and
queued execution while preserving a bounded residual. This is closer to the
SubtreeLU private/pipeline scheduler contract than the previous level-list
cursor, but it still runs the current KLU-backed row numeric updates and is not
the missing CKTSO row/segment kernel by itself.

The queued row tail then stopped allocating and publishing the per-row
completion bitmap for unchecked refactors. The queue already carries the
tail-local group predecessor counts, so the bitmap is only needed by checked
row fast-factor runs that may need prefix-current reject validation. Smoke
coverage now proves the split: the checked queued pass uses the bitmap, while
the later unchecked queued pass reports no last-run bitmap use. This removes a
KLU-wrapper bookkeeping artifact from the experimental row scheduler without
changing default KLS policy.

The row-pattern build then added per-row input-cleanup metadata. For each row,
KLS checks whether every input column is already cleared naturally by the row
numeric pass through an `L` dependency, the pivot, or a `U` entry. Rows with
full coverage skip the old final input-column cleanup loop. The smoke fixture
now proves zero cleanup rows on its dependent row-group pattern, and benchmark
JSON exposes cleanup rows and entries so broad runs can confirm whether the
optimization is structural on larger circuit matrices. This removes another
current-row-kernel bookkeeping pass without tuning on matrix names.

The dense row-segment path then narrowed deferred value scatter to checked
pivot probes. Previously, the presence of any dense row segment forced the
entire row refactor to store into row mirrors and copy all row-major `L`/`U`
values back into KLU storage at the end. The native dense group kernel now also
supports unchecked operation by scattering each completed dense row directly
after its in-group update is finalized. Smoke coverage proves the checked pass
still reports deferred scatter while the subsequent unchecked queued refactor
does not. This removes a broad post-pass copy from the experimental row
refactor without changing default KLS policy.

The row ready queue was then made work-aware using the existing group work
estimate. The retained successor graph is sorted once when the row pattern is
built, and each run sorts the initial tail-ready groups before publishing them
to worker threads. This follows SubtreeLU's queue-balancing motivation of
scheduling heavier refactor work earlier while preserving dependency readiness;
it is still running the experimental KLS row layer, not a full separator-tree
private/pipeline scheduler.

The row-pattern builder now retains the row-to-group map and uses it when a
checked row fast-factor pass rejects a pivot. KLS records the conservative
group-tail restart scope reachable from the rejected row's group: number of
groups, covered rows, and retained group-work estimate. This is a planning and
diagnostic bridge to CKTSO's "rows that need to be recomputed with pivoting"
step; it does not yet execute a parallel pivoting row-tail kernel.

The row ready queue workspace then moved from per-run heap allocation to
solver-owned storage tied to the row-pattern lifetime. The ready group array,
ready-slot atomics, tail-local predecessor counters, and tail-membership bitmap
are now resized only when the retained row task graph grows. This follows the
SPICE repeated-refactor requirement that scheduler metadata be retained rather
than rebuilt from scratch on every Newton step. The same chunk also made
checked queued rejects deterministic by publishing completed generic rows
inside multirow groups and by refreshing any missing rows before the rejected
pivot before the existing prefix proof runs.

Full-graph queued row runs then started consuming cached root groups through a
private-root cursor before falling back to the shared ready queue for
successor-released groups. The completing worker also keeps one ready successor
as a local continuation and spills the rest to the shared queue; benchmark JSON
now reports `row_refactor_last_local_ready_groups` and
`row_refactor_local_ready_group_count`, so broad runs can distinguish "ready
queue active" from actual worker-local continuation use. Focused probes showed
the path is heavily exercised: `G2_circuit` reported 64,962 local continuations
in the last row-refactor run and `rajat24` reported 212,824.

Cluster-prefix row tails now use the same private first-wave treatment when
the initially ready tail groups can be copied into per-thread queues, leaving
later successor releases on the shared/local pipeline queue. This is still a
row-group DAG scheduler, not SubtreeLU's full separator-tree partitioner, but
it removes avoidable shared-queue traffic from the private/pipeline boundary.
Benchmark JSON reports `row_refactor_last_private_ready_groups` and
`row_refactor_private_ready_group_count` for this path.

A follow-up attempt to replace the private-root cursor with static strided
root ownership was tested and rejected. Although this looked closer to a
SubtreeLU private-queue shape, it removed dynamic balancing from the first wave
of root tasks. In same-session focused probes it left `rajat24` essentially
flat (`0.27955s` versus the prior `0.27924s` average refactor) but regressed
`G2_circuit` from the prior `0.423s` range to `0.740s`. KLS therefore keeps the
dynamic private-root cursor until a true separator-tree or FLOP-balanced
private-queue partitioner exists.

The tail-restart summarizer then started reporting serial-suffix overcompute:
for executed local tail restarts it now totals the extra suffix columns and work
above the retained CKTSO-style pivoting-tail plan. This does not change solver
behavior, but it makes broad benchmark output rank the cases where replacing
the conservative suffix fallback with a real pipelined pivoting-tail executor
should remove the most work.

Those overcompute counters are now also published by the solver and benchmark
JSON as `fast_repaired_tail_restart_overcompute_columns` and
`fast_repaired_tail_restart_overcompute_work`, with the standalone summarizer
kept backward-compatible for older JSONL runs. This keeps the CKTSO-tail
executor target visible in every benchmark row instead of requiring a separate
postprocessing calculation.

The retained CKTSO-style pivoting-tail plan then gained an explicit executor
shape classification. Benchmark JSON now reports whether the ordered
ETree-descendant restart set is contiguous, whether it exactly matches the
contiguous suffix that KLS's current serial restart can execute, how many
columns are gaps inside the retained plan, and how much extra suffix work the
current serial executor would do above that plan. These fields do not execute
the non-contiguous tail; they make broad paper-suite runs identify the cases
where KLS needs the real CKTSO pipelined row-tail executor instead of another
safe KLU-style suffix cleanup.

KLS then started consuming a conservative executable subset of that retained
tail plan. When the ETree-descendant plan begins at the rejected pivot and its
last planned column is before the full block suffix end, the KLS-owned pivoted
block repair first tries to refactor only the contiguous envelope covering that
plan and preserve the later columns. Internal envelope gaps that are not in the
ETree-tail plan are locked to their old pivot rows, so a gap column cannot
introduce a new pivot change merely because KLS is using a serial envelope
rather than CKTSO's true non-contiguous pipeline. If pivoting selects a row
whose old pivot position lies outside the envelope, if a locked gap changes
pivot, if preserved `L` rows cannot be remapped to the new final pivot order,
or if the speculative attempt fails, KLS restores the original block pointer
metadata and falls back to the existing serial suffix restart. This is still
not the full CKTSO non-contiguous pipelined tail executor, but it is now an
executable tail-envelope approximation rather than only a diagnostic.

The pivoting-tail plan then stopped treating every prefix-current checked
reject as a full suffix when a narrower restart seed is available. KLS first
uses the retained row-refactor tail for checked row-major rejects; otherwise,
when the checked worker bitmap can identify unfinished nodes, it records those
unfinished local columns as `fast_rejected_pivoting_tail_seed_columns` and
closes only that seed set through the ordered-block ETree. The unfinished
bitmap seed is retained separately from the mutable planning queue, so it
remains available if the row-tail seed is rejected by the topological closure
check. This is still diagnostic/planning infrastructure, but it matches CKTSO's
restart-point determination more closely and exposes the true non-suffix
worklist that a pipelined pivoting-tail executor should consume.

The checked refactor-pool path then stopped discarding all unfinished-work
information when it can prove the prefix is current only at BTF-block
granularity. Those rejects now record the rejected block suffix as a
conservative unfinished seed before the ETree tail plan is built, and
benchmark output reports it separately as
`fast_rejected_pivoting_tail_block_seed_columns`. This does not make the pool
path a full CKTSO Algorithm 5 executor: row-tail seeds still take precedence
when they are valid, and the executable repair remains the existing serial
tail envelope. It does, however, preserve the paper's "unfinished guessed
EGraph nodes seed the ETree tail" information for the pool path instead of
collapsing directly to an untagged suffix fallback.

The same block-suffix seed now covers serial prefix-current checked rejects.
Serial row-major checked refactors, serial mapped checked refactors, and the
serial scaled checked continuation all have a precise interrupted set: the
rejected pivot through the end of its diagonal block. KLS now records that set
with `fast_rejected_pivoting_tail_block_seed_columns` before building the
ordered-ETree closure. A validated retained row-tail still wins when it gives a
narrower repair with preserved-gap refresh; otherwise these paths feed the same
conservative CKTSO Algorithm 5 planning fallback as the pool path instead of
clearing the seed and rediscovering an anonymous suffix later.

The retained CKTSO-style tail plan now also materializes the ETree forest
shape that a true pipelined pivoting-tail executor would consume. Once the
seeded worklist has been closed through ordered-block ETree parents and
validated as topological, KLS retains the tail-local parent, child-count, and
leaf-to-root level arrays and reports retained ETree edges, roots, leaves,
maximum fanout, level count, and maximum ready-level width as
`fast_rejected_pivoting_tail_etree_edges`,
`fast_rejected_pivoting_tail_etree_roots`,
`fast_rejected_pivoting_tail_etree_leaves`, and
`fast_rejected_pivoting_tail_etree_max_fanout`,
`fast_rejected_pivoting_tail_etree_levels`, and
`fast_rejected_pivoting_tail_etree_max_width`. Smoke tests require the retained
tail to satisfy the forest invariant `edges + roots == columns` and to expose a
nonempty ready-level structure. This is still scheduler state rather than the
full CKTSO tail executor, but it fills the paper-level representation gap
between a flat suffix/envelope retry and Algorithm 5's ETree-descendant ready
worklist.

The executable exact-mask repair paths now consume that forest representation
as a validation gate. Before the serial narrow-tail retry or the threaded
row-first active-mask repair accepts a retained tail mask, KLS verifies the mask
against the retained scheduler arrays and checks that its roots, leaves, edges,
maximum fanout, levels, and maximum width match the retained forest counters. If
the proof fails, KLS falls back to the existing wider repair path rather than
treating a generic range mask as a CKTSO-style ETree tail. Accepted repairs
report `fast_repaired_tail_restart_etree_mask` in addition to the older
exact-mask counter.

The serial suffix tail retry then stopped copying and mutating a private
`Offp` array. Tail-column construction now supports a discard-only off-block
mode, which is valid because accepted local repairs rebuild `Offi`/`Offx` from
the final `Pinv` before publishing numeric state. This removes one more
whole-matrix scratch allocation from the executable tail-restart path while
preserving the current conservative suffix semantics.

Unchecked row refactors then stopped publishing every completed row back into
KLU's column storage immediately. The KLS-owned row-major `L`/`U` mirrors now
remain authoritative across repeated unchecked row refactors, and a dirty flag
forces a single publish before a KLU solve, transpose solve, later
`kls_factor`, or non-row refactor fallback. This is a direct SPICE-cycle
optimization for the row/segment engine: repeated Newton refactors no longer
pay KLU scatter traffic on every step when the next step can consume the
row-major mirrors directly.
A focused four-thread `G2_circuit` row-refactor probe with three repeated
refactors reported `0.337s` average refactor and a `0.047s` solve that included
the delayed publish, with valid residuals. This is a row-engine throughput step,
not yet the missing exact pivoting-tail executor.

A guarded row-major solve then consumed those dirty row mirrors directly for
unscaled, normal-orientation, single-block, non-transpose solves with no
external KLS row/column scaling or permutation. It does not replace KLU's
general triangular solve; it only skips the lazy publish when the KLS row
mirrors are already authoritative. Focused checks kept residuals valid:
`G2_circuit` with `--row-refactor all` reported three row solves and reduced
the non-transpose solve average to about `0.020s` before the benchmark's
transpose solve forced the expected publish, while an unscaled `nxp1`
row-refactor probe used the row solve once with about `0.029s` solve time.

The pivoting-tail repair path then moved one step closer to CKTSO Algorithm 5.
When the retained ETree-descendant tail is non-contiguous, the preferred KLS
block repair now tries to copy already-finished gap columns from the current
numeric LU stream and numerically refactor only the marked tail columns. A gap
column is copied only if its pivot row is still unchanged, its L rows are still
unpivoted, and its U predecessors do not include a recomputed tail column. If
that copy is impossible because the gap depends on recomputed tail state, KLS
now promotes the gap into the active tail mask and recomputes it with normal
pivoting instead of forcing the old pivot inside a non-tail gap or immediately
falling back to the wider serial restart. Structural/mapping copy failures still
fall back. Benchmark JSON reports copied gap columns/work as
`fast_repaired_tail_restart_skipped_columns` and
`fast_repaired_tail_restart_skipped_work`, and
`fast_repaired_tail_restart_exact_mask` identifies successful repairs whose
executed columns exactly match the ETree-derived pivoting-tail mask; promoted
gap repairs deliberately clear that flag and show up as tail overcompute. This
is still a serial conservative subset of CKTSO's pipelined pivoting-tail
executor, not the full non-contiguous parallel row-tail algorithm, but successful
masked cases are no longer blended with serial gap refactors and dependency
blocked internal gaps no longer force an all-or-nothing suffix fallback.

The row-up first-factor dynamic column-pivot path then stopped relabeling
swapped columns by scanning every previously emitted U entry. U entries now keep
per-column linked lists, so a dynamic column exchange only touches entries in
the two swapped columns. This directly targets the CKTSO Algorithm 1 pivoting
step in the KLS-owned first-factor scaffold. Focused checks stayed valid: the
small dynamic-pivot smoke still exercises the row-up path, `nxp1` KLS-first
initial factor time dropped from about `6.37s` to about `3.72s` with 494
dynamic column pivots, and `rajat24` moved from about `47.2s` to about
`44.2s` with 2355 dynamic column pivots. The modest `rajat24` change confirms
that relabeling was not the dominant missing paper mechanism there; KLS still
needs a parallel row-up/EGraph first-factor executor rather than only cheaper
serial pivot bookkeeping.
The refactor gap remains; this is a storage-ownership bridge toward the
row/segment engine, not the missing CKTSO pivoting-tail executor.

The experimental row-refactor numeric path was then opened to KLU row-scaled
single-block factors. The row path now recomputes KLU row scales before the
serial or threaded row pass, divides input entries by the unpermuted row scale
using the same fixed input-position mapping as the EGraph refactor path, and
permutes `Rs` back to pivot order only after an accepted pass. Checked rejects
leave `Rs` in input-row order for the existing scaled repair/tail machinery, and
the row-tail candidate diagnostic uses the same scaled input loader. A smoke
case covers both one-thread and four-thread scaled row refactors, validates the
constructed solution, and now confirms that the scaled dirty row mirrors remain
authoritative through the normal solve by dividing the RHS through pivot-order
`Rs`, matching KLU's `P*(R\b)` solve setup. This is a general row/segment-engine
coverage step, not the full CKTSO ETree-descendant pivoting-tail executor.

The dirty row-major solve was then extended to transpose solves for the same
single-block eligibility class. It loads `Q' * b`, scatters along row-major `U`
for the `U'` solve, scatters backward along row-major `L` for the `L'` solve, and
writes through `Pnum`, dividing by pivot-order `Rs` when KLU row scaling is
active. The scaled row-refactor smoke now solves both forward and transpose
systems from dirty row mirrors and checks that no lazy publish occurs. This keeps
the benchmark-style transpose validation from forcing KLU column storage after a
row refactor, while external KLS row/column scaling and permutation still use the
publish-and-KLU fallback.

The row-refactor pattern and numeric pass were then opened from single-block
factors to BTF diagonal blocks. The retained refactor map now supplies only
diagonal-block input entries to the row-major update, translates KLU's local
block `L`/`U` row indices into global row order, and refreshes BTF off-block
`Offx` values separately from the original input positions. A smoke case covers
a reducible BTF matrix whose off-block entry changes across refactor, verifies
that the row refactor ran, and initially confirmed that BTF solves published the
dirty row mirrors before using KLU's general triangular solve. This is still a
prerequisite for CKTSO-style row/segment tails rather than the full
pivoting-tail executor.

The dirty row-major solve was then extended across BTF factors. Forward solves
now process BTF blocks in KLU order from last to first, solve each diagonal
block from KLS-owned row-major `L`/`U` mirrors, and subtract refreshed off-block
`Offx` columns from earlier block rows. Transpose solves process blocks from
first to last, apply `Offx'`, and then run the row-major `U'`/`L'` triangular
passes inside the block. The BTF row-refactor smoke now checks both forward and
transpose solves and verifies that dirty row mirrors remain authoritative
without publishing back to KLU storage. This is a storage-ownership step toward
the CKTSO row/segment engine, not the full pivoting-tail executor.

The dirty row-major solve then switched from one-RHS-at-a-time sparse traversal
to KLU-style chunks of up to four right-hand sides. The same BTF and
single-block row-mirror paths now load, update, and store a small RHS batch
while traversing each `L`, `U`, and `Offx` pattern once per chunk. The BTF smoke
now solves five forward and five transpose right-hand sides for both unscaled and
KLU row-scaled factors. This is a structure-adaptive solve-path cleanup enabled
by KLS-owned row storage; it does not replace the missing pivoting-tail
factorization engine.

The auto symbolic selector then started comparing SCOTCH as a guarded
nested-dissection candidate instead of keeping it explicit-only. KLS only pays
for this comparison on large single-block patterns with high estimated symbolic
work, and only keeps SCOTCH when its symbolic fill/work score is at least a
material win over the current AMD/COLAMD/METIS candidate. This follows the
CKTSO/SubtreeLU combined-ordering motivation without adding matrix-name
tuning, and it still leaves the row/segment pivoting-tail executor as the main
missing CKTSO-scale mechanism.

The checked row-tail diagnostic was then tightened for BTF rejects. When the
retained KLS row-successor graph can build a tail from the rejected global row
and every row in that tail stays inside the rejected BTF block, KLS now reports
that row-owned tail directly instead of first falling back to the older KLU
numeric `L` scan. If the retained graph is missing or crosses the block
boundary, the KLU scan remains the conservative fallback. This does not execute
CKTSO's pipelined pivoting tail, but it makes the available row/segment
worklist explicit for multi-block checked-row failures.

The KLS-owned pivoted block kernel was then wired into an experimental first
factorization scaffold behind `KLS_ENABLE_KLS_FIRST_FACTOR=1`. The scaffold
allocates the KLU-compatible numeric object itself, handles singleton BTF blocks
through `Udiag`/`Pnum`, runs the KLS-owned pivoted kernel for multi-column BTF
blocks, rebuilds `Pinv` and `Offp/Offi/Offx`, seeds KLS-owned row-major value
mirrors from the accepted numeric object, and reports
`initial_factor_path:"kls_first"` when it succeeds. The smoke suite covers a
matrix with a 2-column BTF block so this path cannot pass by singleton handling
alone, and it now requires subsequent forward and transpose solves to consume the
seeded row mirrors. The same KLS-owned pivoted tail now consumes the retained
factor-order input map when that map is valid for the current BTF block, so it
can iterate the already split in-block input slice instead of remapping every
original CSC row for each tail column. `kls_stats`, `kls_bench`, and the gap
decomposition script expose this through `kls_tail_last_mapped_columns` and
`kls_tail_mapped_column_count`; the KLS-first smoke case requires the mapped
tail to be exercised.

The scaffold then gained a more direct CKTSO Algorithm 1 bridge instead of only
wrapping KLU-style columns. Under the same `KLS_ENABLE_KLS_FIRST_FACTOR=1` gate,
eligible no-scale or KLU row-scaled matrices first try a KLS-owned sparse
row-major up-looking first factor with pivot checks over every BTF diagonal
block: each factor row scatters the permuted in-block input row, divides by
KLU's input-row `Rs` when row scaling is active, applies already computed
row-major U updates, checks the diagonal against the remaining U-row maximum,
records row-major L/U entries, and only then packs the accepted local block
factors into KLU-compatible numeric storage for the existing solve/refactor API.
After acceptance, scaled factors permute `Rs` through `Pnum` to preserve KLU
solve semantics. The bridge also implements the direct Algorithm 1 pivot
exchange for this eligibility class: when the diagonal fails the threshold
against the largest active U-tail entry, KLS swaps the current block-local
factor column with that entry, updates previously computed row-major U column
labels, publishes the accepted `Q` permutation, and continues rather than
falling straight back to KLU-compatible block tails. `kls_stats`, `kls_bench`,
and the gap
decomposition script report this through
`kls_first_last_row_uplooking_columns` and
`kls_first_row_uplooking_column_count`, and dynamic exchanges through
`kls_first_last_dynamic_column_pivots` and
`kls_first_dynamic_column_pivot_count`; smoke tests require both the no-exchange
row-major path, a weak-diagonal dynamic column-pivot case, a reducible two-block
BTF case with an off-diagonal coupling, and explicit scaled cases. This is the
first direct row-oriented first-factor kernel in KLS. It is still limited: it
falls back to the older KLS/KLU-compatible first-factor paths on unrecoverable
pivot rejection or unsupported scaling/static-pivoting state. The full CKTSO
production target still needs ETree cluster/pipeline scheduling for first
factorization and the ETree-descendant pivoting-tail executor.
When those mirrors were seeded by `kls_first`, unchecked `kls_refactor` now
automatically attempts the existing KLS-owned row-major refactor path even when
`KLS_ENABLE_ROW_REFACTOR=0`, so the scaffold drives the next SPICE-style numeric
update through KLS row storage instead of immediately returning to KLU's refactor
kernel. Checked fast-factor `kls_factor` calls can use the same row-major
ownership when the mirrors remain current and
`KLS_ENABLE_CHECKED_ROW_REFACTOR=1` requests the checked row executor. The smoke
suite verifies automatic unchecked updates with both row-refactor env gates off,
then explicitly enables the checked gate for checked row update and
forward/transpose solve coverage.
The guarded row-major solve is now allowed when analysis selected internal
transpose orientation as well; solve dispatch already passes the required
internal `kernel_transpose` flag, so auto-oriented benchmark runs can consume
KLS-owned mirrors instead of publishing back to KLU solely because the stored
pattern is transposed. `kls_bench` and `run_bench_suite.py` expose
`--kls-first-factor env|off|on` so this experimental KLS-owned first-factor
path can be compared reproducibly across manifest chunks.
After a KLS-first factor has established automatic row-major ownership,
successful checked fast-factor pivot repairs now reseed the row-major mirrors
from the repaired numeric object when the repair path had to rebuild local LU
storage and drop stale mirror metadata. The guarded solve/refactor lifecycle
therefore stays on KLS-owned row storage after a local repivot, matching the
CKTSO paper's row-major fast-factor/recompute flow more closely even though the
full pipelined ETree pivoting-tail executor is still not implemented.
The automatic KLS-first row-refactor handoff then gained a static work gate:
explicit row-refactor controls still force the row engine for experiments, but
the automatic path now compares the retained row/group work estimate with the
exact EGraph refactor work estimate and skips row refactor when the row plan is
already more expensive. This follows the NICSLU/SubtreeLU recommendation to
select parallel kernels from structure and FLOP evidence rather than matrix
names, and prevents KLS-first scaffolding from replacing a cheaper existing
EGraph refactor with a slower row-major mirror update on hard rows.
The gate is now visible in `kls_stats` and benchmark JSON through
`row_refactor_total_group_work`, `row_refactor_auto_enabled`,
`row_refactor_auto_values_ready`, `row_refactor_auto_work_allowed`, and
`row_refactor_auto_should_run`, so future KLS-first comparisons can explain
whether the row-major path was skipped because the paper-style work model
rejected it. Follow-up diagnostics split that skipped state into
`row_refactor_auto_lower_bound_work`,
`row_refactor_auto_lower_bound_rejected`,
`row_refactor_auto_pattern_build_failed`, and
`row_refactor_auto_value_copy_failed`, so the CKTSO-gap runs can distinguish a
cheap lower-bound rejection from an inability to build or populate the
row-major refactor mirrors.
A 2026-06-28 focused top-ten CKTSO-gap rerun with these diagnostics
(`build/kls_auto_reject_diag_gap10_t4_r1_ref3_timeout120.jsonl`) produced a
`3.998s` geomean and showed the row-refactor model recommended and attempted
the handoff on 9 of 10 matrices, but every attempted handoff was rejected by the
lower-bound work scan. The lower-bound estimates were only about `1.00x` to
`1.01x` above the exact EGraph dependency work, and both
`row_refactor_auto_pattern_build_failed` and
`row_refactor_auto_value_copy_failed` stayed zero. This makes the current gap
more concrete: the default slow cases are not reaching BLAS-backed or native
row panels at all; KLS still lacks a row-major refactor algorithm whose
structural work model clearly beats the EGraph scatter path on these patterns.
A forced row-refactor probe on the same top-ten slice
(`build/kls_forced_row_refactor_gap10_t4_r1_ref3_timeout120.jsonl`) confirmed
that simply bypassing this gate is not a plausible CKTSO-gap fix. It completed
without numerical failures but regressed the geomean to `12.086s`; every row was
slower, with `2.57x` to `3.83x` cycle slowdowns against the guarded default.
The forced path did execute the KLS row engine and activated native/compact row
panel counters on the large ASIC/G2 cases, so the loss is not an admission bug.
It is evidence that the current row-major engine is still a scalar
KLU-compatible scaffold rather than the production row/supernode numeric object
described by CKTSO/SubtreeLU. A separate default-EGraph ready-queue probe
(`build/kls_egraph_ready_queue_gap10_t4_r1_ref3_timeout120.jsonl`) also
regressed to `7.272s` geomean, showing that replacing the existing
cluster/pipeline EGraph tail with dynamic ready-column release is not the
missing large lever either.
A same-session 20-matrix CKTSO-gap probe compared this gated KLS-first mode
with `--kls-first-factor off` at 4 threads and 3 refactors. The geomean ratio
was about 1.01x, with wins on `ASIC_100k`, `rajat15`, `onetone1`, and
`transient`, but regressions on `G2_circuit`, `HTC_336_4438`, `ASIC_100ks`,
and `Raj1`. This confirms that the current scaffold should remain
experiment-gated; the papers point to the full row-major fast factorization and
pivoting-tail machinery as the missing general algorithm, not to blindly
enabling KLS-first on every pattern.
The row-major solve experiment was then decoupled from the KLS-first factor
experiment. `KLS_ENABLE_ROW_SOLVE_FROM_NUMERIC=1` now seeds KLS row-major
`L`/`U` solve mirrors from successful ordinary numeric factors/refactors when
the factor has no external KLS row/column permutation or scaling. `kls_bench`
and `run_bench_suite.py` expose this as `--row-solve env|off|on`, and a smoke
test forces ordinary `klu_first` factorization while verifying both forward and
transpose solves use the KLS row-solve path. This keeps CKTSO's row-oriented
solve idea measurable without implying that KLS is delegating to another solver
as a backend.
The row-solve path then gained a single-RHS scalar loop and cached solve
metadata validation. This matches the row-by-row triangular solve form in
CKTSO's Algorithm 2 for the common SPICE one-right-hand-side case, while
preserving the existing four-RHS batched path for wider solves. On the local
five-matrix smoke suite with `--row-solve on`, the geomean moved from about
`0.0161s` before the scalar path to about `0.0123s` after scalar solve plus
validation caching; the `--row-solve off` reference remained about `0.00946s`,
so row solve remains experiment-gated rather than default.
The solve-only seed then stopped building the full row-refactor group/segment
scheduler. It now builds just the row-major `L`/`U` pattern and value mirror
needed by triangular solve, while KLS-first and auto row-refactor paths still
request the full row-refactor pattern. On the same five-matrix smoke suite, the
`--row-solve on` geomean moved again to about `0.0111s`, and solve-only rows
report zero row-refactor groups/segments while still reporting ready row solve
mirrors.
The lean solve setup now also records CKTSO-style triangular partition
diagnostics without allocating the row-refactor scheduler: lower/upper dense
tail start, rows, and entries, plus the fixed eight trapezoid slices used by
the CKTSO paper. The dense tail criterion follows the paper's setup rule: at
least 70% of triangular entries and at least 300,000 entries in the suffix.
KLS now materializes those slice boundaries internally and reports the maximum
per-slice entry count for lower and upper triangular solves, so the parallel
rectangular-slice executor can gate on measured slice balance.
It also precomputes the CKTSO row segmentation step for dense-tail rows: lower
rows split at the slice start, and upper rows split at the slice end, yielding
rectangular versus within-slice triangular entry counts for both factors.
For one-RHS, single-block normal solves, KLS now runs the rectangular part of
each trapezoid slice in the persistent worker pool, with worker barriers around
the sequential triangular piece. It only enables that executor when the
rectangular entries reach the paper's 300,000-entry dense-tail work scale; the
`bcircuit` structure falls back to scalar row solve, while `G2_circuit` crosses
the gate and records parallel slice runs. Rectangular rows inside each slice are
partitioned by accumulated rectangular nonzeros, following CKTSO's thread
workload assignment rule, and diagnostics report the max per-thread rectangular
entries for lower and upper factors.
The next CKTSO Section V gap was the sparse triangular block before the dense
tail: a local `G2_circuit` probe showed that KLS could parallelize about 3.27M
rectangular entries but still left roughly 9.9M lower/upper prefix plus
triangular-piece entries serial. KLS now levelizes those sparse lower/upper
prefixes, runs wide levels in CKTSO-style cluster mode, and lets thread 0 solve
the remaining narrow levels sequentially. The cluster cutoff uses the same
`#threads * 2` width rule already used in the CKTSO-inspired row-refactor
cluster/pipeline split, and the JSON/text diagnostics report sparse level
counts, cluster levels, max widths, and sparse-level run counts.
Validation then showed that the straightforward pthread-barrier executor still
does not close the solve-time gap by itself: `bcircuit` and `G2_circuit` are
too small or synchronization-heavy, and `G3_circuit` ran the sparse-level path
but moved from about `0.199s` scalar row solve to about `0.213s` parallel row
solve. KLS therefore keeps the implementation but only activates it when the
parallelizable work is a substantial share of total solve work and the average
work per synchronization reaches the paper's 300,000-entry dense-tail scale.
These fields are visible in `kls_stats`, `kls_bench` JSON/text, and the gap
decomposition script, giving the parallel triangular-solve path a
structure-based gate instead of a matrix-name heuristic.
That cost model now also controls seeding row-solve structure from ordinary
numeric storage. Earlier `--row-solve on` probes copied `L`/`U` values into
row-major mirrors after each successful ordinary factor/refactor even when the
parallel row-solve executor never ran, adding an `O(nnz(L+U))` tax to the
slow repeated-refactor cases. The adaptive seed first builds the cheap
CKTSO-style partition diagnostics and only prepares row-solve values when the
predicted parallel row-solve work is large enough. On the local top-12 hard-case probe,
forced row-solve seeding had a `3.1772s` geomean, adaptive seeding had a
`3.0613s` geomean with no failures, and explicit `--row-solve off` was
`3.0426s`. This is a useful cleanup, but it is not the missing CKTSO numeric
engine.
The ordinary-numeric row-solve seed no longer performs that value copy even
when the gate accepts. It now keeps the row solve's structural metadata and
reads factor values through the same retained `double *` slots already used to
scatter KLS row-refactor mirrors back into KLU storage. A temporary
`4000 x 4000` lower-triangular probe with `--row-solve on` crossed both
CKTSO-style dense-tail gate terms and reported `row_refactor_last_row_solve=1`
with no row-refactor groups, proving the solve can run directly from ordinary
numeric storage. This removes an `O(nnz(L+U))` setup copy from the experiment,
but it still leaves the larger paper gap in the factor/refactor numeric engine.

The CKTSO Section V solve executor then gained the transpose-side equivalent
instead of using row-scatter loops for transpose solves. KLS now builds retained
transposed row views for `U^T` and `L^T`; each entry points back to the
authoritative KLU-compatible factor slot, so dirty row mirrors and ordinary
numeric value pointers keep the same ownership rule. Single-block, single-RHS
transpose solves use the same sparse-prefix levelization, dense-tail trapezoid
slices, rectangular/triangular split points, worker-pool dispatch, and
work/share gates as normal solves, but execute as dependency-gather triangular
solves to avoid concurrent scatter races. A `4000 x 4000` natural upper
triangular smoke fixture now solves `A^T x=b` through this path and requires a
recorded parallel row-solve run. This closes the obvious normal-vs-transpose
triangular-solve mismatch in the CKTSO Section V implementation; it does not
complete the separate compact-panel transpose group executor.

The same rule now applies to KLS-first row-refactor mirrors. KLS-first still
builds row-refactor metadata so `row_refactor_total_group_work` and
`row_refactor_auto_work_allowed` stay visible, but it copies values into
row-major mirrors only when the retained row-work model allows the row engine
to run. On the local top-six hard-case slice, this removed an unused setup copy
on the rows where `row_refactor_auto_work_allowed=0`; the same-session default
probe moved from `7.7455s` to `7.5856s` geomean. A top-12 guard with
`--row-solve on` completed all rows at `3.0847s`, within noise of the previous
adaptive `3.0613s` run and still faster than the older forced-seed `3.1772s`
baseline.
The row-major solve eligibility was then widened for static-matching factors.
The pre/post solve wrapper already maps normal solves through `R * P * b` and
returns `C * y`, and maps transpose solves through `C * b` and `P' * R * y`.
That means dirty KLS-owned row mirrors and adaptively seeded mirrors from
ordinary numeric storage can be used without publishing to KLU just because
the analysis selected the MC64-adjacent external row permutation or matching
equilibration. A smoke case now refactors a scaled pre-static permutation
matrix and requires both normal and transpose solves to consume the dirty
row-major mirrors.
The row-refactor scheduler was also corrected to preserve CKTSO's
cluster/pipeline split. It now runs the wide cluster prefix chosen by the
paper's `2 * threads` width rule before building the successor-ready queue for
the remaining narrow tail; the previous full-graph ready queue sent very wide
ASIC/G2 levels through hundreds of thousands of tiny atomic tasks. A forced
top-six row-refactor probe improved only from `11.4507s` to `11.3434s`
geomean, and remained much slower than the EGraph path. This confirms that the
current row engine still lacks the paper's production row/segment numeric
kernel and that automatic row-refactor activation should remain cost-gated.
The ready-queue row-refactor tail now also gives each worker a private initial
ready-root range balanced by cached row work before dependent successors enter
the shared queue. This is a small SubtreeLU-style private/pipeline scheduling
step that reduces root-queue contention without changing the row numeric
formulas or claiming the full separator-tree scheduler.

That scaffold was then extended to KLU row-scaled first factors. It computes
`Rs` in input-row order before constructing singleton and multi-column BTF
blocks, rebuilds off-block values before scale permutation, and finally permutes
`Rs` through `Pnum` for solve semantics. The smoke suite now runs the same
2-column BTF case with no scaling and with KLU max-row scaling.

The KLS-first pivoted-block path no longer allocates a fake LU payload before
root block factorization. The shared pivoted-block kernel now accepts an empty
block when the rejected prefix is zero and estimates its own initial LU memory;
non-root tail restarts still require reusable prefix LU. This removes a
wrapper-style placeholder from first factorization while preserving the existing
block-repair semantics.

Re-reading the local CKTSO, NICSLU, and SubtreeLU references leaves one clear
large missing part for the slow cases: KLS still does not own a complete
row/segment-oriented numeric factorization and refactorization engine. CKTSO's
documented advantage is not only METIS-style ordering or MC64-style matching; it
uses row-major sparse LU, EGraph/ETree cluster-pipeline scheduling, fast
factorization with pivot checks, and ETree-descendant tail factorization when a
pivot check fails. SubtreeLU pushes the same direction through separator-tree
private/pipeline queues and row/supernode updates. KLS has implemented the
ordering, LGPL-compatible matching/scaling boundary, diagnostics, and some
row-solve/refactor scaffolding, but the core repeated-iteration speedup in the
papers comes from doing numeric update in that owned row/segment schedule rather
than repeatedly adapting KLU-compatible numeric storage.

## Recommended General Work

1. Build a KLS-owned row/segment-oriented numeric engine instead of adding more
   wrapper-level gates around KLU storage. It should preserve enough row-major
   `L`/`U` access to run CKTSO-style no-pivot fast factorization with pivot
   checks, identify the ETree-descendant restart tail after a failed check, and
   later support SubtreeLU-style separator-tree private/pipeline queues. The
   current KLS-owned block-local pivot restart, suffix restart, and large
   single-block no-pivot cluster/pipeline refactor are useful precursors, but
   the target is a pivoting tail restart inside large blocks.
2. Continue turning matching/scaling into a production MC64-equivalent stage,
   but keep it inside the LGPL-compatible boundary: use the BSD-licensed SPRAL
   scaling submodule, system SPRAL, or independent KLS code, not HSL MC64 or
   restricted MC64 copies from other solver trees. The retained SPRAL path now
   helps large weak-diagonal dominant-BTF cases and avoids replacing no-BTF
   ordering wins, but `pre2` still times out, so matching quality alone is not
   the remaining CKTSO-scale gap.
3. Finish the structure-adaptive triangular solve only after the LU storage
   owned by KLS exposes row-oriented or segmented access cheaply enough that
   solve setup and value mirrors do not dominate repeated refactors.
4. Use static symbolic and numeric-cost models to decide whether a parallel
   kernel should run, so KLS avoids matrix-name-specific tuning. The NICSLU
   R1/R2 counters are now reported and can seed KLS-owned row/segment refactor
   preparation for scheduled numeric states. Benchmark artifacts distinguish
   the model recommendation, attempted preparation, and accepted preparation
   with `row_refactor_auto_model_recommended`,
   `row_refactor_auto_model_attempted`, and
   `row_refactor_auto_model_accepted`; the remaining work is to turn the same
   row/segment substrate into the full pivoting factor and tail-restart
   executor described by CKTSO/SubtreeLU.

The automatic KLS-first row-refactor handoff now applies that last rule before
building the full row-refactor pattern: it scans retained `L`/`U` structure and
input entries to form a conservative lower bound on row-update work, and skips
row metadata setup when even that bound exceeds the exact EGraph refactor work.
On the current 20-row CKTSO-gap focus guard with four threads and three
refactors, this pre-gate skipped row metadata on 11 rows, left nine ASIC/G2/DC
style rows to the exact post-build cost gate, and improved geomean cycle time
by about 2.6% versus the previous copy-gate artifact. This is useful static
cost discipline, but it also confirms the paper reading above: the worst
ASIC/G2/mc2 rows still need the larger CKTSO/SubtreeLU row/segment numeric
engine rather than another wrapper-level policy tweak.

The next forced-row measurements made that distinction sharper. On the top six
current CKTSO-gap rows, explicit `--row-refactor refactor` produced an
`11.31s` geomean cycle time versus `7.55s` for the auto-gated default, with
ASIC/G2 rows regressing by roughly `1.36x-2.20x`. Two same-session experiments
were rejected. First, specializing row input and off-block refresh for unscaled
values regressed the forced-row geomean to `11.75s`; it helped `ASIC_320k` by
about `2%` but lost `7%-11%` on `ASIC_320ks`, `ASIC_100ks`, and `G2_circuit`.
Second, running unchecked row refactors as a full retained group-DAG ready queue
instead of preserving the CKTSO-style cluster prefix regressed geomean to
`11.59s`. The full ready queue did improve `onetone2`, but the large rows
created tens of thousands of local ready continuations, so scheduler churn
outweighed barrier savings. Both tests were reverted. The retained lesson is
that the hard gap is inside the row/segment numeric representation and update
kernel, not a missing unscaled input branch or a blanket full-ready scheduler.

The saved large-paper reconnaissance pair tells the same story at larger scale.
With failures scored at `120s`, KLS was slightly ahead of the saved CKTSO large
artifact in geomean (`0.986x` candidate/reference), because it wins
`TSOPF_FS_b39_c30` and `rajat29` and both solvers struggle on `Hamrle3`.
However, the median ratio was still `1.18x` against KLS and the material losses
were refactor-heavy rows: `nxp1` (`1.67x`), `G3_circuit` (`1.56x`),
`rajat30` (`1.42x`), and `ASIC_680k` (`1.37x`). This keeps the next useful
implementation target aligned with the papers: persistent compact row/segment
numeric storage, batched trailing updates, and pivot-aware tail restart, rather
than another scheduling-only or input-copy-only tweak.
The phase split on those same artifacts rules out triangular solve and
rejected-pivot repair as the main explanation for these rows. Each of the four
losses reported `fast_block_restarts=0` and `fast_tail_restarts=0`. KLS solve
time was already faster than CKTSO on all four, while numeric factor/refactor
was slower: `nxp1` refactor `0.3065s` versus `0.1552s`, `G3_circuit`
`11.03s` versus `7.07s`, `rajat30` `0.2243s` versus `0.1375s`, and
`ASIC_680k` `0.0499s` versus `0.0309s`. Initial factorization was also much
slower on the same cases (`3.13s` versus `0.63s`, `44.16s` versus `5.28s`,
`2.21s` versus `0.58s`, and `1.84s` versus `0.16s`). This is why the local
CKTSO and SubtreeLU papers now point KLS toward the row-major up-looking
factor/refactor kernel with supernode/segment updates first; the full
ETree-descendant pivoting-tail restart remains necessary for stressed
repivoting cases, but it is not what explains these no-reject large losses.

As a small row-engine step after re-reading the CKTSO refactorization section,
the serial row-refactor path now executes the retained row-group processor with
a single local worker instead of keeping a separate scalar row loop. This makes
the dense/generic row-segment kernels available to one-thread refactors and to
single-thread checked fast-factor attempts, preserving the same prefix-scatter,
singular, and rejected-pivot bookkeeping. The smoke suite now asserts that the
dense checked-row prefix repair is a serial checked row-refactor over a dense
segment with deferred scatter. This is still only an incremental move toward
the paper target: KLS remains tied to KLU-compatible value storage for these
rows, so the larger missing item is still persistent compact row/segment
numeric storage and batched trailing updates.

The next CKTSO/SubtreeLU-aligned row-engine step made the compact scratch-panel
path active for unchecked dense row groups when the estimated internal dense and
shared trailing update work crosses a structural gate. The dense-group
dispatcher now tries the compact row-major panel before the native direct
row-mirror kernel for those unchecked groups, falls back to native if the
scratch panel cannot be allocated, and reports
`row_refactor_last_compact_dense_panel` plus
`row_refactor_compact_dense_panel_count`. A generated 48-by-48 dense
unchecked smoke case asserts that the compact-panel path is selected and keeps a
small residual.

A follow-up fixed the checked compact-panel semantics instead of keeping checked
dense groups gated to the native path. The failed trial showed that computing
the whole dense panel before testing any checked pivot can touch suffix state
that CKTSO-style prefix/tail repair expects to remain row-ordered. KLS now uses
a row-ordered compact checked path: after gathering the panel it updates one
row, publishes that row into the row-major mirrors, checks the pivot, and stops
immediately on reject. The dense checked-prefix smoke case now asserts both the
compact-panel marker and the tail repair residual. This is still not the full
pipelined CKTSO tail executor, but it closes a concrete semantic gap for
compact row-major segment updates with pivot-aware prefix repair.

The compact-panel gate was then tightened from a tiny absolute work floor into
a structural arithmetic-intensity check. A dense group now uses worker-local
compact scratch only when its estimated internal dense plus shared-trailing
update work is large enough overall and large enough per copied panel entry.
This avoids treating shallow panel packing as a SubtreeLU/CKTSO mechanism when
the copied scratch has too little update work to amortize it. Focused forced-row
checks stayed valid: `G2_circuit` still exercised 662 dense segments but compact
panel uses dropped from 1280 to about 500 across two refactors, while `mc2depi`
dropped from 3623 to about 1340 compact uses and improved in the same-session
sample. The retained lesson is still that KLS needs persistent compact
row/segment numeric storage, but the scratch bridge now follows a broader
work-per-byte rule instead of a low absolute threshold.

The compact-panel gate is now visible in benchmark artifacts. `kls_stats` and
`kls_bench` report compact-panel eligible dense-group count, eligible rows,
estimated update work, and copied panel entries in addition to the last-run
execution counter. This separates three paper-relevant cases in future suite
runs: no executable dense row segments, dense segments that are too shallow for
compact scratch, and dense segments whose arithmetic intensity is high enough to
exercise the SubtreeLU-style compact-panel bridge. The synthetic dense checked
and unchecked smoke cases now assert those eligibility counters, while small
`add20` still reports zero eligible compact panels.

The gap decomposition CSV now carries the same compact-panel eligibility and
execution counters. This keeps the medium/large CKTSO comparison workflow
aligned with the row-segment diagnostics, so future slow-row reviews can see
whether a refactor-heavy loss has no dense segment work, dense work rejected by
the compact arithmetic-intensity gate, or compact-panel execution that is still
too slow because KLS lacks persistent row/segment numeric storage.

Re-reading the local CKTSO Section IV algorithms also exposed one remaining
serial bridge in the existing block-repair path: after a prefix-current BTF
block repair, KLS recomputed all later diagonal blocks serially even though
those BTF blocks are independent. The solver-owned refactor pool can now start
at a later BTF block, marks earlier blocks as already current for prefix
classification, and the unscaled prefix-current fast-factor restart path tries
that pool before falling back to the old serial suffix. Benchmark stats report
`fast_repaired_parallel_tail_blocks`, and a reducible smoke fixture asserts that
an early repaired block continues over later BTF singleton blocks through the
pool. This follows CKTSO's "continue tail work in parallel" direction only at
BTF-block granularity; it does not implement Algorithm 5's ETree-descendant
row-tail factorization with pivoting, nor SubtreeLU's separator-tree
private/pipeline row queues.

The compact dense row-segment bridge then moved one step from transient scratch
toward persistent row/segment storage. When a dense row group passes the same
compact-panel arithmetic-intensity gate, the row-refactor pattern now allocates
a solver-owned compact panel slice for that group and reports the retained
groups/entries as
`row_refactor_compact_dense_panel_persistent_groups` and
`row_refactor_compact_dense_panel_persistent_entries`. The compact kernel uses
that retained slice before falling back to worker-local scratch, reports actual
retained-slice execution as
`row_refactor_last_compact_dense_panel_persistent` and
`row_refactor_compact_dense_panel_persistent_run_count`, and the dense
checked/unchecked smoke cases assert both execution and retained panel
consumption.
This still repacks current input values on each numeric pass and scatters back
to row mirrors, so it is not the full CKTSO/SubtreeLU row-major numeric storage
model; it does make the compact row-segment value lifetime solver-owned rather
than worker-scratch-owned.

The row solve path now consumes those retained compact panels directly when a
dense group has a complete valid prefix. Row-aware solve value accessors read
in-group `L`, in-group `U`, and trailing-panel entries from the solver-owned
row-major compact panel before falling back to the older sparse row-value
mirrors, and the public stats report last/cumulative compact-panel solve
values. This is still a storage-boundary step rather than the full paper row
numeric engine, but it removes another forced round trip through scattered
KLU-style value arrays from the rows that already have compact SubtreeLU-shaped
storage.

The one-RHS row solve then moved from per-entry compact-panel reads to
group-level compact-panel execution for complete dense groups. The normal
forward solve validates the dense group's lower suffix layout and walks the
row-major lower triangle directly; the backward solve validates the dense and
shared trailing U layout and walks the retained row-major upper/trailing panel
before dividing by the published pivots. Both paths fall back to the previous
sparse row loops if the panel or row layout is not exact, and stats now report
last/cumulative compact group-solve rows and panel entries. This is still a
serial one-RHS kernel, but it is a direct step from "compact panel as storage"
to "compact panel as triangular-solve executor."

The same compact group-solve executor now covers the normal multi-RHS row
solve chunks. KLS validates the retained dense group layout once per group and
updates the chunk's row-major work slab across all RHS columns before falling
back to sparse row loops for non-exact layouts. The compact-panel smoke fixture
now solves three RHS after the one-RHS solve and requires fresh last-run
compact group-solve rows, so the coverage is tied to actual multi-RHS executor
use rather than cumulative one-RHS statistics.

The serial row-refactor transpose solve now consumes the same retained compact
groups instead of walking every in-panel value through row-aware scalar
accessors. The `U^T` phase validates the dense upper panel plus its common
trailing strip, solves the transposed lower triangular panel in row-major form,
and scatters only the external trailing rows. The `L^T` phase validates the
lower-panel suffix, solves the transposed unit-upper panel backward, and
scatters only external lower dependencies. One-RHS and four-RHS chunks share
the same layout guards and fall back to the old scalar loops on non-exact
groups. This closes the serial solve-side panel-executor asymmetry; the
remaining transpose solve gap is the generic parallel transposed-row view,
which is still source-position based rather than a compact-panel group
executor.

The refactor scalar fallback now follows the same storage rule for single-row
producer dependencies. If a completed dense producer row belongs to a retained
compact panel, scalar consumers validate that panel's dense/trailing U-row
layout and apply the update from the row-major panel before falling back to the
scattered row mirror. New stats report last/cumulative scalar rows and entries
served from compact panels, and a focused smoke fixture forces the one-row
producer case that the multi-row compact-supernode update intentionally skips.
This does not add the missing coarse BLAS row-panel engine, but it makes
retained compact panels a broader numeric source for the existing row executor.

KLS now retains METIS `NodeNDP` separator-tree queue metadata instead of
discarding it after ordering. For accepted METIS symbolic analyses, the
`NodeNDP` size tree is converted to a postorder private/pipeline component
sequence: leaf domains are private work components, internal separators are
pipeline components, and accepted ordering positions keep a component map for a
future numeric queue consumer. The stats and benchmark output report analyzed
rows, component counts, private/pipeline row totals, and max component sizes.
A focused `nxp1` analyze-only probe with four threads, METIS ordering, and BTF
disabled reported 414604 analyzed rows, seven components, four private leaf
components, three pipeline separator components, and a 414147/457
private/pipeline row split. This closes the direct state-retention gap against
SubtreeLU's separator-tree setup, but the gap that matters for CKTSO-scale slow
cases remains the numeric consumer: pivot-constrained separator work queues,
FLOP-balanced refactor queues, and row/supernode update kernels are still not
implemented.

The retained separator map is now consumed by one numeric path: the experimental
row-refactor ready queue. When the retained `NodeNDP` map covers the full
factor order and every initial ready row group stays inside one separator
component, KLS assigns those initial ready groups to per-thread private queues
by separator component before dependent successors enter the shared/local ready
queue. Benchmark stats report
`row_refactor_last_separator_private_queue`,
`row_refactor_separator_private_queue_run_count`,
`row_refactor_last_separator_private_components`, and
`row_refactor_separator_private_component_count`. A forced `G2_circuit`
METIS/KLS-first/row-refactor probe with four threads selected the
separator-private path, using four retained components for six initial private
groups, with a valid residual. This is a real SubtreeLU-aligned consumer of the
retained separator tree, but it is still narrower than the paper: it does not
constrain pivot search inside separator-tree subdomains, does not partition all
work by separator FLOP balance, and does not implement the pivoting
first-factor or checked-tail kernel that the slow CKTSO-gap rows still require.

KLS now consumes the same retained separator map more directly for no-pivot row
refactorization through a SubtreeLU Algorithm 6-style queue builder. When the
separator tree covers the row-refactor block, KLS computes retained row-group
work per separator component and subtree, repeatedly moves the dominant
subtree root component to a pipeline queue while returning its child subtrees
to the candidate set, then greedily assigns the remaining subtrees to private
thread queues by work. Row-refactor groups that cross separator-component
boundaries cannot be private-subtree work, so KLS classifies them as pipeline
groups instead of rejecting the whole separator queue. If any private subtree
group depends on a pipeline group, KLS promotes the dependent group into the
pipeline closure and keeps the separator schedule rather than falling back to
the generic ready queue. Benchmark stats report
`row_refactor_last_separator_flop_queue`,
`row_refactor_separator_flop_queue_run_count`,
`row_refactor_last_separator_flop_components`,
`row_refactor_last_separator_flop_private_groups`, and
`row_refactor_last_separator_flop_pipeline_groups`, plus closure promotions via
`row_refactor_last_separator_flop_closure_groups` and
`row_refactor_separator_flop_closure_group_count`. A forced `G2_circuit`
METIS/no-BTF/no-scale/no-pivot row-refactor probe with four threads selected
this queue on all three refactors, using seven separator components, 69,402
private groups and five pipeline groups per refactor, with residual
`3.32e-16` and a repeated-refactor average around `0.20s`. This fills the
paper's FLOP-balanced separator-queue mechanism for the eligible no-pivot row
engine; the remaining SubtreeLU gap is applying the same separator partition
inside checked-tail/pivoting kernels, BTF forests, and deeper supernodal BLAS
numeric storage.

KLS-first row-up factorization now builds and conservatively consumes the
corresponding first-factor queue shape from the retained separator map. For
each covered symbolic block, KLS now first applies the retained separator tree
as an Algorithm 6-style split/collapse: it computes component and subtree work
from permuted block-local row input counts, repeatedly moves the dominant
subtree root into the factor-order pipeline queue, returns child subtrees to
the private candidate set, and greedily assigns remaining subtrees to private
threads by work. That replaces the older direct leaf/private split when the
partition is dependency-safe. The numeric consumer validates the private phase
before threading it: if a private row would read another private thread's
mutable column domain, KLS discards the partitioned queue and falls back to the
legacy retained-component private/pipeline queue. When validation passes, each
thread factors its assigned private rows into local row-up `L`/`U` entries. KLS
now first remaps the block-local row/column order and published `Pnum` to that
private-then-pipeline queue before numeric assembly; this fixes the earlier
unsafe version that tried to execute queue order over the old KLU numeric row
indices. KLS then merges private entries and consumes the pipeline rows with an
Algorithm 3-style atomic task counter. Each pipeline worker now keeps local
scratch, applies all ready private-row predecessors against an immutable
private-U snapshot before waiting for its ordered publish slot, and can consume
already-published earlier pipeline-prefix rows while waiting for intervening
pipeline rows that it does not depend on. Published private and pipeline rows
are also grouped into conservative row-supernodes when adjacent row-major `U`
patterns have the full triangular extension and identical trailing pattern.
Partial pipeline-row updates consume such consecutive ready predecessors, or
the already-finished prefix of a retained predecessor supernode, through one
guarded supernode-run executor. The executor keeps the triangular in-run
discovery behavior, so an update from row `k` can create row `k+1` as the next
dependency before the run proceeds. It also validates the same
dense-prefix/common-trailing shape described for supernodes, solves the
in-supernode prefix in the current sparse row workspace, accumulates the shared
trailing contribution in worker-local scratch, and scatters each trailing column
once. This matches the Algorithm 4 update control flow and scalar compact
numeric shape while still using the current KLS row-entry storage rather than
production BLAS panel storage. The publish step can now accept a
scoped dynamic column exchange for a safely pre-updated pipeline row: while
holding the ordered pipeline lock, it swaps the committed U entries and the
phase-local prefix snapshot, publishes the pivot row, increments a
column-order epoch, and wakes workers. Any speculative suffix row that began
under the older epoch clears its local sparse workspace and recomputes the same
row under the new column order instead of forcing a whole suffix relaunch.
This is still not production BLAS-backed supernode storage or the CKTSO
checked-tail pivoting executor. Benchmark stats expose the last
planned queue through
`kls_first_last_separator_queue`,
`kls_first_last_separator_queue_private_components`,
`kls_first_last_separator_queue_pipeline_components`,
`kls_first_last_separator_queue_private_rows`,
`kls_first_last_separator_queue_pipeline_rows`,
`kls_first_last_separator_queue_nonempty_threads`, and
`kls_first_last_separator_queue_max_thread_rows`, with the structural work
range reported by `kls_first_last_separator_queue_min_thread_work` and
`kls_first_last_separator_queue_max_thread_work`. Algorithm 6-style partition
use is reported through `kls_first_last_separator_queue_partitioned`,
`kls_first_separator_queue_partitioned_count`, and
`kls_first_last_separator_queue_split_components`; scheduled consumption is
reported by `kls_first_last_separator_queue_executed`,
`kls_first_separator_queue_executed_run_count`,
`kls_first_last_separator_queue_executed_private_rows`, and
`kls_first_last_separator_queue_executed_pipeline_rows`; threaded private
execution is reported by `kls_first_last_separator_queue_parallel_private`,
`kls_first_separator_queue_parallel_private_run_count`,
`kls_first_last_separator_queue_parallel_private_rows`, and
`kls_first_last_separator_queue_parallel_private_threads`; threaded pipeline
queue consumption is reported by
`kls_first_last_separator_queue_parallel_pipeline`,
`kls_first_separator_queue_parallel_pipeline_run_count`,
`kls_first_last_separator_queue_parallel_pipeline_rows`, and
`kls_first_last_separator_queue_parallel_pipeline_threads`; private-predecessor
pipeline pre-updates are reported by
`kls_first_last_separator_queue_pipeline_partial`,
`kls_first_separator_queue_pipeline_partial_run_count`,
`kls_first_last_separator_queue_pipeline_partial_rows`, and
`kls_first_last_separator_queue_pipeline_partial_threads`; waiting-prefix
pipeline pre-updates are reported by
`kls_first_last_separator_queue_pipeline_wait_partial`,
`kls_first_separator_queue_pipeline_wait_partial_run_count`,
`kls_first_last_separator_queue_pipeline_wait_partial_rows`, and
`kls_first_last_separator_queue_pipeline_wait_partial_deps`; scalar
row-supernode partial updates are reported by
`kls_first_last_separator_queue_pipeline_supernode_update`,
`kls_first_separator_queue_pipeline_supernode_update_run_count`,
`kls_first_last_separator_queue_pipeline_supernode_update_groups`, and
`kls_first_last_separator_queue_pipeline_supernode_update_rows`; pipeline
pivot-tail restarts are reported by
`kls_first_last_separator_queue_pipeline_pivot_tail`,
`kls_first_separator_queue_pipeline_pivot_tail_run_count`, and
`kls_first_last_separator_queue_pipeline_pivot_tail_rows`, with restart and
actually serialized pivot-row counts reported by
`kls_first_last_separator_queue_pipeline_pivot_restarts`,
`kls_first_separator_queue_pipeline_pivot_restart_count`, and
`kls_first_last_separator_queue_pipeline_pivot_serial_rows`. Prefix
private-predecessor pre-updates, waiting-prefix pre-updates, and scalar
row-supernode updates completed before a pipeline pivot restart remain counted
in the same counters after the epoch retry, while
`kls_first_last_separator_queue_pipeline_pivot_serial_rows` remains reserved
for the older external serialized fallback. After a successful dynamic pivot
inside the separator pipeline, KLS now rebuilds the phase-local panel cache
from the post-exchange column order over the whole committed prefix instead of
disabling the cache for the rest of the phase; suffix rows that restart under
the new epoch can still consume validated dense/common-tail producer panels.
The rebuild is reported through
`kls_first_last_separator_queue_pipeline_prefix_panel_rebuild`,
`kls_first_separator_queue_pipeline_prefix_panel_rebuild_count`, and
`kls_first_last_separator_queue_pipeline_prefix_panel_rebuild_rows`. Smoke
tests cover the queue shape on a 30,000-row METIS-ordered tridiagonal
KLS-first factor and the epoch recovery path on a 30,600-row METIS-ordered grid
with separator-band weak diagonals. The queue-shape test requires more than
one private worker thread, verifies that all planned pipeline rows are consumed
by the guarded pipeline executor, verifies that the private-predecessor partial
pre-update path touches all pipeline rows, and now requires cached
panel-backed producer updates in the separator pipeline. The epoch test
requires a separator-pipeline pivot restart, a prefix-panel rebuild whose
covered rows exceed the private-prefix rows, cached separator-pipeline panel
updates after that restart, zero serialized pivot rows, and a clean solve
residual.
This closes Algorithm 3's queue shape and adds race-free scalar Algorithm 4
partial-update and scoped-pivot retry steps for dependency-safe retained
separator queues.

For KLS-first blocks without a retained separator order, the row-up first factor
now still uses the restartable Algorithm 5-style row pipeline when multiple
threads are available. The executor builds a natural factor-order row list,
lets workers perform partial dependency updates and ordered publication, and
uses the same epoch retry path for dynamic column exchanges. These runs are
reported separately from separator queues through
`kls_first_last_row_pipeline`, `kls_first_row_pipeline_run_count`,
`kls_first_last_row_pipeline_rows`, `kls_first_last_row_pipeline_threads`,
`kls_first_last_row_pipeline_partial`,
`kls_first_row_pipeline_partial_run_count`,
`kls_first_last_row_pipeline_partial_rows`, and
`kls_first_last_row_pipeline_partial_threads`. Dynamic pivot epochs in this
generic row pipeline are now reported separately through
`kls_first_last_row_pipeline_pivot_tail`,
`kls_first_row_pipeline_pivot_tail_run_count`,
`kls_first_last_row_pipeline_pivot_tail_rows`,
`kls_first_last_row_pipeline_pivot_restarts`,
`kls_first_row_pipeline_pivot_restart_count`,
`kls_first_last_row_pipeline_pivot_serial_rows`,
`kls_first_last_row_pipeline_prefix_panel_rebuild`,
`kls_first_row_pipeline_prefix_panel_rebuild_count`, and
`kls_first_last_row_pipeline_prefix_panel_rebuild_rows`, while separator
pipeline epochs stay on the separator-prefixed counters. Smoke coverage
verifies both a clean natural-order execution with no pivot epoch and a
two-thread weak-pivot natural-order execution that rebuilds the committed-prefix
panel cache without setting separator-pipeline counters. This does not add
CKTSO's ETree-descendant pivot-tail executor, but it removes the prior drop to
a purely serial KLS-first row factor whenever the separator queue was
unavailable and makes the generic Algorithm 5-style epoch behavior measurable.

The KLS-first row up-looking dynamic column pivot selector also now consumes
the retained separator map when it is available for the full factor order. On a
weak pivot, KLS now applies the SubtreeLU Algorithm 4 `N'` rule directly:
`N'` is the last retained factor row of the current collapsed component, the
diagonal is compared only with the strongest candidate in `i+1..N'`, and
outside-domain candidates no longer force a separator-scoped rejection through
the old global row maximum. If no safe scoped candidate exists for a required
pivot, KLS rejects the row-up attempt and lets the existing KLS-first fallback
path use the pivoted block kernel. This matches SubtreeLU's requirement that
private-mode pivots do not cross retained component domains. Stats report exact
separator-domain pivots through
`kls_first_last_separator_dynamic_column_pivots` and
`kls_first_separator_dynamic_column_pivot_count`, component-extent pivots
through `kls_first_last_separator_extent_dynamic_column_pivots` and
`kls_first_separator_extent_dynamic_column_pivot_count`, would-have-crossed
fallback candidates through `kls_first_last_separator_dynamic_column_fallbacks`
and `kls_first_separator_dynamic_column_fallback_count`, and strict rejections
through `kls_first_last_separator_dynamic_column_rejects` and
`kls_first_separator_dynamic_column_reject_count`. The existing smoke suite
still exercises KLS-first dynamic column pivoting. This should be treated as a
guarded Algorithm 4 pivot-domain step rather than evidence that the full
SubtreeLU private/pipeline pivoting executor is implemented.

The mapped EGraph fast-refactor path now mirrors the CKTSO checked-reject
prefix recovery that already existed in the row-refactor path. When a checked
EGraph worker rejects a pivot and the pipeline completion bitmap exists, KLS
serially recomputes any unfinished columns before the rejected pivot with the
same mapped EGraph column kernel, marks those columns done, and only then
classifies the reject as prefix-current or unknown. This fills a concrete
CKTSO Section IV gap: a rejected fast factorization no longer loses the valid
prefix just because some earlier EGraph tasks were not scheduled before the
stop flag. Benchmark JSON now reports
`fast_rejected_prefix_refresh_columns` and
`fast_rejected_prefix_refresh_count` so slow-case reruns can distinguish true
unknown-prefix failures from recovered prefix-current ETree-tail candidates.
The implementation is still a serial prefix recovery feeding the existing
serial pivoting-tail/block-repair path; it is not yet CKTSO's full parallel
ETree-scheduled pivoting tail executor.

The compact dense row-refactor kernel then removed a redundant panel sweep in
the unchecked refactor path. Previously the compact dense group first factored
the dense intra-segment panel and then made a second pass over the same
row/dependency pairs to update the trailing panel. The unchecked path now
updates the trailing panel while applying each dependency, matching the
checked compact kernel's dataflow and moving the implementation closer to the
supernode-style dense update direction in CKTSO/SubtreeLU. On the forced
row-refactor focus probe, `G2_circuit` repeated refactor improved from about
`0.507s` to `0.404-0.446s`, `onetone2` improved from about `0.0348s` to
`0.0329s`, and a repeated `ASIC_100ks` probe reported about `0.116s` versus
the previous one-pass `0.148s`, with residuals unchanged. Default automatic
selection remains cost-gated because the row engine is still not broadly
faster than the mapped EGraph refactor.

KLS now maps another explicit SubtreeLU Algorithm 5 detail into the
row-refactor scheduler. When large tail row groups with downstream successors
and width at least `2 * threads` dominate the tail by both group count and row
count, the scheduler prepares a row-dependency ready queue. Dense-group rows
are marked complete immediately after each row is numerically stored and, for
checked runs, after the pivot check passes. Consumers can therefore wait on the
exact finished row prefix instead of the whole supernode. This matches the
paper's large-unfinished-supernode split for no-pivot refactorization and the
checked row fast-factor/refactor path; `KLS_ENABLE_PARTIAL_SUPERNODE_PIPELINE=0`
remains a hard disable. Direct same-tree A/B probes on `coupled` and
`G2_circuit` showed earlier broad selectors were slower than the existing
ready queue on the current KLS row kernel, so the strict structural selector is
kept. It is still not the full CKTSO pivoting-tail factorization.

KLS now consumes retained compact dense row panels as external supernode
update sources during row refactorization. Once a compact dense group finishes
successfully in the current numeric pass, a per-group valid marker lets later
rows recognize a contiguous suffix of two or more dependencies from that same
group and update from the contiguous panel directly instead of expanding each
predecessor row through the sparse row mirror. This implements the
SubtreeLU/CKTSO supernode
update idea more directly than the earlier compact-panel work, which only used
the panel while factoring the producer group itself. The path preserves the
same multiplier checks and row-prefix publication rules as the scalar row
kernel, and falls back to scalar updates when a producer panel is absent,
stale, unfinished, or belongs to the current group. Focused forced-row probes
showed the path active with valid residuals: `coupled` used 1,770 compact
supernode updates over 62,162 predecessor rows, `G2_circuit` used 119,130
updates over 6,220,352 rows, and `ASIC_100ks` used 23,209 updates over
915,915 rows. This is still not full BLAS-backed supernodal factorization or
the CKTSO pivoting-tail executor, but it moves the current row engine from
mere compact storage toward actually consuming supernodes in later updates.

The partial supernode pipeline now closes a direct Algorithm 5-style gap in
that consumer path. Dense and generic producer groups publish the valid prefix
of their retained compact panel as each row is completed, and the ready-queue
scheduler releases successor groups when the exact row-dependency predecessor
count reaches zero. A consumer whose dependency suffix ends at the published
prefix consumes the compact panel immediately; consumers are no longer
speculatively enqueued at the `width - threads` split point. New diagnostics report
`row_refactor_last_compact_supernode_partial_update` and the corresponding
partial-update count, row, and entry totals. The smoke fixture forces a
240-row producer with a 2-thread split and requires a 238-row partial compact
supernode update with a residual-clean solve. This is still not CKTSO's full
pivoting-tail executor or a production BLAS supernodal numeric object, but it
fills the specific missing "start consumers from a finished producer prefix"
semantic instead of treating every compact supernode as all-or-nothing.

The speculative split-point release was removed after direct forced-row
inspection showed it was the large-case pathology, not a scheduler deadlock.
With the older code, `onetone2` timed out under a 35 s forced-row probe while
workers were sampled inside dense group processing and prefactor waits. After
keeping only exact row-dependency release, the same command completed in
`0.203s`, matching the partial-pipeline-disabled control (`0.199s`) with clean
residuals and zero row-prefactor counters. Larger spot checks also completed
cleanly: `G2_circuit` in `1.028s`, `ASIC_100ks` in `0.435s`, and `ASIC_320k`
in `0.508s`. This keeps the paper's row-prefix publication machinery but
rejects the current KLS implementation of early consumer enqueue as too
expensive on dense row groups.

The scheduler now preserves SubtreeLU Algorithm 6's separator private/pipeline
partition when applying row-dependency release. A new row-dependency
ready-queue preparation path can reuse the separator private-group mask, so
private separator groups still run through the retained FLOP-balanced queue
while pipeline groups may be released by exact row-predecessor completion when
the strict large-unfinished-supernode test passes. Checked refactorization uses
the same strict selector when `KLS_ENABLE_CHECKED_ROW_REFACTOR=1`: producer rows
are not released until their pivot checks complete, and the smoke fixture
validates that explicit checked path. A
broader "any large supernode" selector was tested and rejected: on the six-row forced
row-refactor focus probe it became much slower than the separator-only run and
was interrupted after exceeding the normal short-run envelope. With the strict
selector, `build/kls_head_gap_focus6_forced_row_sep_alg5_strict_t4_r3_timeout120.jsonl`
completed all six rows without timeout, but the focus matrices stayed on the
separator queue (`row_refactor_last_partial_supernode_pipeline=0`) and reported
a `15.29s` SPICE-cycle geomean. This closes the scheduling semantic gap without
pretending it is a performance win on the current focus subset.

KLS then removed the remaining dominance selector from the ordinary SubtreeLU
Algorithm 5 row-tail path. The paper's condition is local to each unfinished
dependent supernode: if the supernode has at least `2P` rows, a consumer may
use the completed `k:(k'-P)` prefix instead of waiting for the whole producer.
KLS now selects the row-dependency queue whenever such a dependent
dense/generic producer exists in the pipeline tail; the large producer no
longer needs to dominate the number of tail groups or tail rows. Separator
Algorithm 6 private/pipeline queues still keep their stricter wrapper guard
until private-queue row-dependency release is made fully safe. Private queue groups
are protected from speculative enqueueing and receive normal row-edge
decrements. A new non-dominant-tail smoke fixture places one large dependent
producer behind many small dependent groups and requires the partial-supernode
pipeline by default, covering the paper condition that the large producer need
not dominate the tail.

The compact supernode consumer now applies the paper's matrix-vector update
shape for producer trailing panels. For a later row that
depends on a completed dense producer suffix, KLS still computes and checks the
`L` multipliers in row order, but it no longer scatters each producer row's
trailing contribution directly into the sparse work vector. Instead, it stores
the multipliers in worker scratch, accumulates the producer trailing panel into
a contiguous temporary vector, and scatters that vector once to the sparse
columns. This is an in-KLS BLAS-shaped `gemv` dataflow, not an external BLAS
dependency and not CPU-specific tuning. The smoke fixture now separates a dense
producer, lower consumer rows, and later trailing columns so it asserts this
path. Focused forced-row probes showed the accumulator active with valid
residuals: `coupled` accumulated about 1,770 updates over 62,000 rows and
6.6 million trailing entries with residual `8.97e-16`; `G2_circuit`
accumulated about 119,000 updates over 6.2 million rows and 2.0 billion
trailing entries with residual `3.32e-16`; and `ASIC_100ks` accumulated about
23,000 updates over 916,000 rows and 355 million trailing entries with residual
`1.63e-15`. This fills a
direct SubtreeLU supernode-update detail. The exact high-volume diagnostic
counts can vary slightly in threaded runs because these counters are telemetry,
not synchronization state. Full BLAS-backed supernodal
factorization, separator FLOP-balanced queues, and CKTSO's pipelined pivoting
tail executor remain open.

The default compact-supernode trailing accumulator was then tightened to avoid
an extra producer-suffix pass. The scalar suffix solve still uses sparse `x`
because that was faster on ASIC/G2-style probes, but KLS now accumulates the
contiguous trailing vector while each `L` multiplier is already live and only
delays the final scatter. Focused forced-row probes stayed residual-clean and
showed the default `gemv` path active with `trsv=0`: `coupled` refactor about
`0.0050s`, `ASIC_100ks` about `0.0866s`, and `G2_circuit` about `0.238s` in
same-session one-pass samples. This is a production-path cleanup of the
SubtreeLU-shaped update, not another env-only experiment.

The scalar compact-supernode `trsv` path is now the default algorithmic branch
for every ready producer run with at least two rows, instead of only an explicit
experiment or a work-threshold decision. With
`KLS_ENABLE_COMPACT_SUPERNODE_TRSV` unset, KLS copies the producer suffix from
sparse `x` into contiguous worker scratch, solves it with the producer dense
upper panel, and feeds it to the existing trailing-panel accumulator whenever a
ready compact producer run exists. Setting the variable to `0` disables the
automatic choice, and `1` still forces it for probes. Stats report
`row_refactor_last_compact_supernode_trsv`,
`row_refactor_compact_supernode_trsv_count`,
`row_refactor_compact_supernode_trsv_rows`, and
`row_refactor_compact_supernode_trsv_entries`, and the compact-panel smoke
fixture now verifies the unset automatic mode. The partial-prefix smoke now
also forces the same compact `trsv` executor while using SubtreeLU Algorithm 5's
large-supernode split: KLS solves only the completed producer prefix in
contiguous scratch, scatters dense updates that target unfinished producer rows
back to the sparse work row, and later consumes the tail through the normal
row-done waits. This closes the earlier gap where partial-prefix consumers
could start early but had to use scalar per-row producer updates. On the current
forced-row top slice, this
reduced the completed-row geomean from 20.5 s to 18.1 s; `G2_circuit`'s
explicit row-refactor average moved to about 0.236 s with a clean residual.
The default CKTSO-gap focus still keeps row refactor gated off on the large
ASIC/G2 rows, so this is a direct SubtreeLU supernode-kernel step rather than a
claim that the row engine is ready to replace the exact EGraph path.

KLS now also has an opt-in CBLAS supernode experiment. Configure with
`-DKLS_ENABLE_CBLAS_SUPERNODE=ON` and set
`KLS_ENABLE_CBLAS_SUPERNODE=1` at runtime to use standard CBLAS calls for the
row-major compact supernode updates. The consumer-side update now follows the
SubtreeLU text directly for a completed producer supernode: CBLAS `dtrsv`
solves the producer suffix multipliers and CBLAS `dgemv` applies the retained
trailing panel to the current sparse work row, with the same row-order
multiplier checks before publication. This consumer-side BLAS path is gated by
a structural work-per-copied-entry rule and explicit minimum row/vector/panel
dimensions, because ungated and weakly gated probes issued thousands of tiny
CBLAS calls and regressed badly despite producing valid residuals. The first
large ASIC check also exposed an unfair build artifact: `build-cblas` had been
configured without `CMAKE_BUILD_TYPE=Release`, so even
`KLS_ENABLE_CBLAS_SUPERNODE=0` measured about `100s`. After reconfiguring the
CBLAS tree as Release, the final guarded forced-row probe on
`ASIC_320k`/`ASIC_320ks` measured `38.79s`/`42.87s`, matching the CBLAS-disabled
Release control (`38.84s`/`42.15s`) and the non-CBLAS Release forced-row build
(`39.51s`/`42.51s`). The default auto row-refactor path with CBLAS enabled
remained at `15.81s`/`12.83s` on the same two matrices. Larger gated probes did
exercise the external CBLAS consumer and stayed residual-clean, but they still
did not beat the scalar row kernel in same-session samples: `G2_circuit` was
about `2.50s` versus `0.326s`, and `ASIC_100ks` was about `0.181s` versus
`0.156s`.

The June 27, 2026 CBLAS rerun confirms that "too many small BLAS calls" is not
the dominant current CKTSO gap. The default `build/` tree still has
`KLS_ENABLE_CBLAS_SUPERNODE=OFF`, so the top default CKTSO-gap run cannot enter
external CBLAS at all. In the explicit CBLAS build, the identical-binary
runtime comparison on forced row refactor for `ASIC_320k`/`ASIC_320ks` measured
`42.71s`/`45.39s` with `KLS_ENABLE_CBLAS_SUPERNODE=0` and
`42.90s`/`45.42s` with `KLS_ENABLE_CBLAS_SUPERNODE=1`, a `1.0026x`
candidate/reference ratio with no >2% wins or losses. The external BLAS entry
points are already guarded by minimum row/vector/panel dimensions and
multi-million-operation work thresholds; the cheap two-row threshold belongs to
the native KLS compact triangular-run kernel, not CBLAS. This keeps BLAS as an
opt-in large-case experiment while the remaining large gap stays with the
checked row/scheduler algorithm rather than with BLAS call granularity.

A June 28, 2026 same-binary rerun keeps that conclusion current after the
compact-`trsv` auto-off change. `build-cblas` was rebuilt successfully, then
the top-five CKTSO-gap matrices were run with the runtime CBLAS experiment off
and on:
`build-cblas/kls_cblas_threshold_verify_off_gap5_t4_r1_ref3_timeout120.jsonl`
reported a `7.78748s` geometric mean, while
`build-cblas/kls_cblas_threshold_verify_on_gap5_t4_r1_ref3_timeout120.jsonl`
reported `7.84702s`. Both runs had `build_has_cblas=true`, but all recorded
external EGraph CBLAS and row-supernode counters were zero on every matrix, so
the default path did not enter BLAS at all. The forced row-refactor top-two
check likewise did not reveal a missing small-case guard:
`KLS_ENABLE_CBLAS_SUPERNODE=0` measured `30.8739s` geomean in
`build-cblas/kls_cblas_threshold_verify_forced_off_gap2_t4_r1_ref3_timeout120.jsonl`,
and `KLS_ENABLE_CBLAS_SUPERNODE=1` measured `30.6218s` in
`build-cblas/kls_cblas_threshold_verify_forced_on_gap2_t4_r1_ref3_timeout120.jsonl`.
The retained row-supernode shape counters stayed unchanged. The current code
already guards CBLAS with 512-row/vector or 512-width panel minima plus
multi-million-operation work thresholds; adding another small-case BLAS guard
would not address the observed CKTSO gap.

KLS then added the next, more paper-faithful batch shape: if an unchecked dense
consumer group, or a contiguous row subrange inside it, has the same ordered
list of completed dense producer suffixes as its external dependency pattern,
KLS gathers those consumer rows into the full consumer panel, solves each
producer's consumer-row multipliers as one right-side triangular solve, and
applies each producer trailing panel as one batched panel update before the
consumer panel is factored. The planner now lets an earlier producer update
later producer multiplier columns before those later suffixes are solved, so a
consumer group can batch across multiple completed producer supernodes instead
of requiring a single external producer, and no longer needs every row in the
dense consumer group to share the same external pattern. The original CBLAS
implementation remains an optional backend, but the same planner now has an
LGPL KLS-owned scalar executor for the triangular solves and trailing panel
updates. The smoke fixtures construct two dense producers feeding one dense
consumer and a separate one-producer case where only a 60-row suffix of an
80-row dense consumer group is batchable; both fixtures force the CBLAS runtime
flag off and require batch counters in normal and no-METIS builds. KLS also
reports batch-pattern, batch-candidate, and work-rejected counters so real
matrices can distinguish "no common row subrange" from "candidate too small for
the structural batch gate." A local generated `onetone2_mwmatch` probe after
adding those counters reported 62 compact panels and 17,995 scalar
compact-supernode updates over 417,868 dependency rows, but zero batch patterns
and zero batch candidates. This proves the direct SubtreeLU update shape is
executable in KLS for a broader producer/consumer pattern, while also showing
that the next real-matrix gap is heterogeneous-row producer batching or deeper
row/segment planning, not merely lowering a backend work threshold.

KLS then widened that consumer-side batch shape to the next direct paper case:
contiguous unchecked rows may share one completed dense producer supernode even
when each row starts at a different suffix of that producer. The dense
consumer-panel executor packs those ragged suffixes into a KLS-owned scalar
workspace, solves each row's producer suffix against the retained dense panel,
and applies the producer trailing panel into the current dense consumer panel
before the internal mini-solve. A second executor handles independent
`GROUP_BATCH` rows by keeping each row's pivot and U entries in a temporary row
panel, applying the same ragged producer update across the row batch, and then
publishing the completed row-major values. New smoke fixtures force both
shapes with `KLS_ENABLE_CBLAS_SUPERNODE=0`, so this is not a BLAS backend
dependency. A forced-row top-five CKTSO-gap probe still reported zero batch
candidates on `ASIC_320k`, `ASIC_320ks`, `onetone2`, `ASIC_100ks`, and
`G2_circuit` while retaining hundreds of thousands to millions of scalar
compact-supernode dependency rows. That narrows the remaining paper gap again:
the slow cases are not waiting for identical-pattern or single-producer ragged
suffix batching; they need broader multi-producer, non-contiguous row-panel
numeric execution over the row/segment storage described by CKTSO/SubtreeLU.

KLS also has an opt-in scalar executor for the next contiguous independent-row
paper shape: when `KLS_ENABLE_MULTI_PRODUCER_SUPERNODE=1` is set, a
`GROUP_BATCH` subrange whose full `L` row is an ordered sequence of completed
dense producer suffixes can use a temporary row panel. It solves each producer
suffix against the retained compact producer panel, applies that producer's
trailing panel into later producer multipliers, the pivot, or the row's `U`
workspace, and then publishes the independent rows as a batch. A smoke fixture
sets `KLS_ENABLE_MULTI_PRODUCER_SUPERNODE=1` and validates two dense producer
groups feeding the same independent row batch without using CBLAS. The top-five
CKTSO-gap probe still reported zero batch candidates with this path both
disabled and enabled on `ASIC_320k`, `ASIC_320ks`, `onetone2`, `ASIC_100ks`,
and `G2_circuit`; the scalar compact-supernode update counters remained large.
This rules out contiguous independent multi-producer rows as the large missing
slow-case feature and leaves the direct paper gap at a more structural level:
KLS needs a broader CKTSO/SubtreeLU row/segment panel executor that can batch
heterogeneous, non-contiguous, and internal rows rather than only completed-
producer suffixes exposed by the current row-major storage.

KLS now records that row/segment structure explicitly during row-refactor
symbolic setup. After dense group and compact-panel detection, each row's `L`
pattern is scanned once for maximal contiguous references to completed compact
dense producer groups. The compact-supernode update kernel consumes this
persistent dense-producer run plan before falling back to the old runtime scan.
`kls_stats` and `kls_bench` expose the total planned producer runs, rows,
dependency rows, max runs per row, full-suffix producer runs, multi-run rows,
and fragmented rows. The two-producer smoke fixture requires these counters,
which means KLS now has a concrete symbolic substrate for the CKTSO/SubtreeLU
"for each contributing supernode, update the current row/range" algorithm
instead of rediscovering one suffix at a time inside the numeric kernel. The
remaining gap should be attacked by turning the fragmented and multi-run rows
seen on paper matrices into a broader row/segment panel executor. A forced-row
top-five CKTSO-gap probe after this change reported 4,057-59,571 planned dense
producer runs per matrix over 1,178-9,789 rows, with max 8-27 runs per row and
zero current batch candidates; fragmented rows equaled producer-run rows on all
five matrices. That is a direct paper-algorithm gap: the supernode producer
ranges are present, but KLS still lacks the row-panel executor that can consume
them together with the scalar/non-producer pieces of the same row.

KLS retains that first fragmented dense-consumer executor as an opt-in path under
its structural and work gates. When `KLS_ENABLE_MULTI_PRODUCER_SUPERNODE=1`
is set for an A/B run, unchecked compact dense consumer groups can batch
consecutive rows whose external prefixes mix scalar dependencies with multiple
completed dense producer suffixes. The executor uses the retained compact
dense input panel when possible, processes scalar gaps only up to the next
dense producer run, then solves and applies each planned producer suffix as a
coarse row-panel stage across the batch before moving to the next scalar gap.
Producer trailing updates can feed later external multipliers or the current
dense/trailing panel, and the completed consumer panel is still left for the
normal internal dense-group factor step. Smoke fixtures set
`KLS_ENABLE_MULTI_PRODUCER_SUPERNODE=1`, construct one scalar external row, two
dense producers, and one dense consumer block, and vary producer suffix starts
by row so the older exact-pattern batch path cannot explain the result.

This closes only the structural-dispatch gap, not the performance gap. On the
same five forced-row CKTSO-gap cases, the path fired on all matrices
(`352-2,136` batches and `246,898-1,754,124` batched dependency rows), proving
that the hard rows are dense-consumer fragmented producer rows. The staged
row-panel and compact-direct-input version improved the five-case probe with
this executor enabled from about `33.2s` to about `32.1s`, and improved the
focused `ASIC_320k` refactor probe from about `0.347s` to about `0.319s`. But
the same five-case scalar baseline without this executor was about `13.5s`,
and the focused baseline was about `0.223s`. Keeping this executor opt-in
preserves the paper-algorithm dispatch scaffold without treating the scalar
staging as SubtreeLU's production supernodal kernel. The remaining paper-level
step is still native row/segment panel storage and updates.

The fragmented dense-consumer executor now also prebuilds a bounded symbolic
target map for the active producer run's trailing updates. Each mapped target
is classified once as a later external multiplier, current dense-panel entry,
pivot, current trailing-panel entry, or no-op before the numeric panel update is
applied. This removes the numeric apply loop's destination rediscovery while
keeping memory bounded by `batch_rows * max_run_trailing_len`, not by all runs
in the fragmented row. A focused `ASIC_320k` probe was mixed (`~0.326s`
refactor versus `~0.319s` before the target map), but the same five forced-row
CKTSO-gap probe improved slightly from about `32.1s` to about `31.9s`
geomean. This is useful as a storage/scatter scaffold, but it confirms that the
main remaining paper gap is still native row/segment panel storage and a
production blocked update kernel, not symbolic target lookup alone.

The scalar gaps around fragmented producer runs now also use a worker-owned
stamped row-target map for the current row's external, dense, pivot, and
trailing destinations. This removes the repeated binary searches that dominated
the forced-row profile while keeping the original sorted-search path as an
allocation fallback. A June 28, 2026 forced-row `ASIC_100ks` probe improved the
no-multi-producer path from the earlier `~0.1067s` release artifact to
`0.0936s`; the multi-producer path was still slower at `0.2157s`. On the first
five CKTSO-gap forced-row matrices, the current opt-in multi-producer executor
measured `28.90s` cycle geomean, while leaving it off measured `16.13s`. That is
a structural staging penalty, not a BLAS or CPU-specific threshold, so the
fragmented multi-producer executor remains opt-in until it has native
row/segment panel storage.

The same executor now keeps its fragmented batch metadata in worker-owned
index/byte scratch, and keeps the per-batch pivot vector in the existing double
workspace, instead of allocating row offsets, run offsets, target positions,
target kinds, and pivots for each small batch. This is a production-storage
cleanup rather than a numeric algorithm change. The focused `ASIC_320k` probe
stayed about `0.326s`, but the five forced-row CKTSO-gap geomean improved
slightly again from about `31.9s` to about `31.7s`. A separate temporary
rectangular run-panel plus portable blocked multiply experiment was tested and
rejected in the same session because it regressed the five-case geomean to
about `32.7s`; copying tiny batches into a denser temporary panel is not enough.
The remaining path needs broader row-panel batches and native row/segment
storage that avoids the copy rather than just a local blocked multiply.

KLS now retains the dense-producer target maps symbolically with the producer
run metadata instead of rebuilding them inside each fragmented dense-consumer
batch. Each producer trailing column is classified once as no-op, later
external multiplier, current dense-panel entry, current pivot, or current
trailing-panel entry; the numeric executor stores the accepted global run IDs
for the batch and reads those persistent classifications directly. A synthetic
smoke fixture covers the retained-map path with two dense producers, a dense
consumer block, and producer trailing columns that update the consumer trailing
panel. This also fixed a signed-size guard that had been discarding otherwise
valid retained maps on 64-bit platforms. On `ASIC_320k`, the retained map count
is now nonzero (`3,205,236` entries: `1,742,094` no-op, `316,396` external,
`305,279` dense, `4,477` pivot, and `836,990` trailing), and the focused
four-thread refactor probe improved to about `0.307s`. The same five forced-row
CKTSO-gap probe improved from about `31.7s` to about `30.3s` geomean, still far
from the no-opt-in baseline near `13.5s`. This closes the symbolic target-map
gap from the paper audit but reinforces that the large remaining gap is the
native row/segment panel representation and production blocked update executor,
not another small matcher or CPU-specific kernel tweak.

The same fragmented dense-consumer executor now participates in checked row fast
factorization. The prior guard skipped this batch path whenever pivot checks were
active, even though the CKTSO/SubtreeLU row-supernode update idea still applies
to the already-completed external producer runs before the consumer group's
internal pivot-checked panel step. KLS now performs the same row-wise multiplier
rejection checks before publishing external `L` multipliers from scalar gaps and
producer suffix solves, then lets the existing checked compact-panel factor code
handle internal multipliers and pivots. The retained target-map smoke fixture now
runs the two-producer fragmented dense-consumer shape through checked
factorization with CBLAS disabled and requires the fragmented batch counters,
retained producer target counts, and residual-clean solve. This closes a direct
checked-mode paper coverage gap while leaving the larger native row/segment
numeric engine and CKTSO pipelined pivoting tail open.

The exact-pattern and ragged dense-consumer batch executors now have the same
checked-mode treatment. These executors cover the cleaner SubtreeLU-style cases
where a dense consumer row range depends on completed dense producer suffixes
with either identical external patterns or one ragged producer suffix. They were
previously disabled whenever pivot checks were active, so checked fast
factorization fell back to scalar external updates even though the later
consumer-panel factor step already performs checked internal multiplier and pivot
tests. KLS now runs the exact/ragged producer-suffix batch under checked row fast
factorization and applies row-wise multiplier rejection before each external
`L` value is published. The existing exact and ragged smoke fixtures each run a
fresh checked factorization with CBLAS disabled and require the compact
supernode batch counters plus a residual-clean solve. This fills another direct
paper-coverage hole without claiming to solve the remaining production
row/segment panel engine or CKTSO pivoting-tail executor.

The independent `GROUP_BATCH` row-supernode executors now have the same checked
coverage. The ragged single-producer suffix executor, the contiguous
multi-producer dense-suffix executor, and the fragmented scalar-gap plus
multi-producer executor no longer skip checked row fast factorization. Each
path rejects unsafe `L` multipliers before publication, tests the completed
row pivot against the row's pivot/U maximum before publishing `Udiag`, and only
then marks the independent row done. The smoke suite now runs checked
factorizations for ragged, contiguous multi-producer, and fragmented independent
row batches with CBLAS disabled and requires the compact-supernode batch and
dense-producer counters. This closes the checked-mode dispatch gap for the
paper row-supernode batch shapes while leaving the larger native row/segment
numeric engine and CKTSO pivoting-tail executor open.

An exact-match requirement was deliberately kept for fragmented dense-producer
batches. A common-prefix widening experiment was tried and rejected because it
increased the work on `ASIC_320k` (`~0.387s` focused refactor versus the prior
`~0.329s`) by admitting extra scalar suffix work into the batch. That result is
consistent with the paper lesson: broadened dispatch is not enough when the
underlying row-panel storage and update kernel are still the limiting pieces.

Three further local row-panel shortcuts were tested after retained target maps
and rejected. Streaming producer trailing updates through a one-row scratch
vector avoided the `batch_rows * trailing_len` temporary, but changed the access
pattern without reducing arithmetic: the five forced-row CKTSO-gap geomean
regressed from about `30.27s` to `30.55s`. Computing only retained-map targets
whose class was not no-op also regressed (`32.59s` geomean), and the aggregate
GEMV counters stayed essentially unchanged on the executed batches, showing
that the global no-op target count was not the active-batch bottleneck. A
portable batched TRSV over equal producer suffixes passed smoke but regressed
the focused `ASIC_320k` refactor from about `0.307s` to `0.314s`, because the
fragmented batches are too small and ragged for that loop interchange to
amortize its overhead. A temporary retained dense-input scatter-map experiment
for compact dense rows also regressed (`32.82s` five-case geomean), so repeated
input destination lookup alone was not the large missing paper mechanism. KLS
later kept the cleaned retained-target compact loader for storage-path
completeness and coverage, not as a claimed performance lever. These rejected
variants narrow the next useful implementation target: KLS needs a native
row/segment numeric representation that stores and updates producer/consumer
panels in the execution order directly, rather than more symbolic shortcuts
around the current KLU-shaped row mirrors.

The unchecked producer-panel refactor experiment also uses a blocked panel
algorithm: scalar code factors each diagonal block, `dtrsm` solves the
below-panel multiplier block, and `dgemm` updates both the dense right panel
and the shared trailing panel.
The scalar compact kernel remains the default because same-session `onetone2`
forced-row probes still favored it: the blocked CBLAS panel path was about
`0.072s` versus `0.029s` with the scalar fallback, with both runs
residual-clean. This improved on the earlier per-row CBLAS attempt (`0.219s`)
but confirms the paper gap more precisely: KLS needs batched row-group/supernode
consumer updates and a deeper blocked row-major numeric layout, not BLAS calls
wrapped around each current row in the existing compact group shape. The CBLAS
probe is now restricted to large dense/vector shapes; small checked cases fall
through to the scalar compact kernel instead of paying BLAS call overhead.

KLS then filled a narrower but direct CKTSO tail-restart semantic gap for
prefix-current masked tails. A prefix-current root reject normally means every
later unstarted column belongs to the unfinished seed set, so the safe tail
remains a full suffix. When the reject-only ETree closure is shorter, however,
KLS can now
prove that preserved columns outside that closure have no U dependency on tail
columns, refresh those preserved columns with the existing mapped no-pivot
column kernel, and then execute the shorter pivoted tail envelope. The refresh
proof now works for single-block, BTF, unscaled, and KLU row-scaled diagonal
blocks by using the mapped global-column dispatcher, requiring retained
block/off-diagonal metadata before refreshing omitted block-local columns, and
feeding scaled attempts with input-row `Rs`. The smoke fixtures cover both a
weak root in a 2-by-2 dependent part plus an independent trailing singleton,
and a non-root weak pivot whose independent changed suffix column must be
refreshed before KLS preserves it outside the ETree-derived tail. Both assert
that the retained pivoting-tail plan is shorter than the block suffix and that
KLS counts the repair as `fast_tail_restarts=1`.
This is still a serial envelope rather than CKTSO's parallel Algorithm 5
executor, but it directly applies the paper's distinction between unfinished
EGraph nodes and ETree-descendant pivoting-tail work instead of treating all
safe prefix-current rejects as whole-block repairs.

The first-factor path was then moved closer to the paper by adding an explicit
KLS-owned row-up-looking scaffold. `KLS_ENABLE_KLS_FIRST_FACTOR=1` forces large
eligible symbolic blocks through that scaffold before falling back to KLU, while
`0` remains a hard off switch. The default cold first factor has since been
moved back to the accepted KLU/static path for broad large cases: focused
CKTSO-gap checks showed that automatic KLS-first cold starts spent substantially
longer in the initial factor without reducing repeated EGraph refactor time
enough to pay for it. The remaining automatic cold-start exception is
structural rather than matrix-named: low-density moderate dominant-BTF matrices
with thousands of fringe blocks and an AMD/auto ordering can still use the
row-up scaffold, because disabling it sent `onetone2` into the 120 s
first-cycle timeout while the row-up path completed. Denser many-BTF relatives
such as `onetone1`, and larger dominant-block cases such as `rajat28`, stay on
the KLU/static path because the row-up exception regressed them in focused
checks. More importantly for static-pivoting and recovery states in that same
low-density many-BTF class, an accepted pre-static row-matching candidate can
still replay through the KLS row-up factor in automatic mode unless
`KLS_ENABLE_KLS_FIRST_FACTOR=0` hard-disables it; this keeps those static-match
states from returning to a KLU numeric that can stall the first SPICE cycle. If
the replay fails, KLS restores the accepted trial numeric; if it succeeds,
`initial_factor_path` reports `kls_first` and row-major mirrors are prepared as
for an ordinary KLS-first factor. This does not remove the KLU trial used to
score static-pivot candidates, and it does not implement CKTSO's parallel
row-up/EGraph first-factor executor, but it keeps the row-up engine available
for low-density moderate many-BTF cold starts, static-match recovery, and explicit
production-candidate experiments rather than as a default cold replacement for
broad CKTSO-gap runs.

The row-up first factor now also keeps its row-oriented product alive for the
next phase. While packing the accepted factors into the KLU-compatible numeric
object for fallback solves and pivot repairs, KLS records the same row entries,
factor-order input positions, and numeric value pointers into the shared
row-refactor CSR/group metadata finisher. A successful direct handoff is
reported by `kls_first_last_row_refactor_seeded_rows` and
`kls_first_row_refactor_seeded_row_count`; smoke coverage requires it for the
ordinary row-up, dynamic-column-pivot, BTF, large-auto, and pre-static replay
first-factor cases. This still is not CKTSO's fully row-major primary numeric
object, because BTF off-diagonal refresh and fallback coherence still require
the retained column-oriented refactor map, but it removes the previous
pack-then-reconstruct step for KLS-first row-refactor mirrors.

The pivot-repair path now attempts the retained non-contiguous CKTSO-style
ETree tail mask even when the pivoting tail reaches the end of the BTF block.
Columns outside the unfinished ETree closure are preserved only when the
existing structural dependency proof can lock their old pivots; otherwise KLS
falls back to the full suffix restart. This fills a direct semantic gap in the
current tail executor without pretending to implement CKTSO's full parallel
pipelined tail factorization.

The compact dense row-segment refactor path then removed another KLU-storage
adaptation step. For dense groups with retained compact panel slices, KLS now
loads current input values for in-panel dense columns and shared trailing
columns directly into that compact panel, leaving the sparse work vector only
for external dependency columns and their update deltas. Checked and unchecked
compact dense smoke cases require this direct-input path, and benchmark JSON
reports the last-run and cumulative direct compact-input row counts. This is a
small but direct move toward the SubtreeLU/CKTSO row/segment numeric-storage
model: dense panel values live in the row-segment panel for the numeric pass
instead of being staged through a KLU-shaped sparse accumulator first.

The same direct-input idea now covers native dense row segments that do not use
the retained compact-panel path. For dense groups below the compact arithmetic
intensity gate, KLS writes current in-group lower, diagonal, upper, and shared
trailing input values directly into the row-segment `L`/`U` mirrors and keeps
the sparse work vector for external dependencies and their update deltas. A
46-by-46 dense smoke fixture is deliberately large enough to form a dense
segment and deliberately below the compact-panel gate; it now requires native
dense direct-input rows and zero compact-panel executions. This removes another
KLU-shaped staging step from the default row-segment kernel while preserving
the conservative compact-panel gate.

The retained SubtreeLU separator-tree queue is now consumed by the checked
fast-factor row path as well as the unchecked no-pivot refactor path. Previously
the Algorithm 6-style FLOP-balanced private/pipeline queue was explicitly
disabled when pivot checks were active, leaving checked fast factorization on
the older generic ready queue even when a METIS `NodeNDP` separator map covered
the block. The checked row scheduler now keeps the same separator private
groups and pipeline groups while relying on the existing earliest-reject stop,
done bitmap, and prefix-refresh validation before any CKTSO tail repair is
accepted. A generated 30,600-row sparse-grid smoke fixture uses METIS,
no BTF, no scaling, and checked row fast factorization, and requires
`row_refactor_last_separator_flop_queue=1` with retained separator components,
private groups, pipeline groups, and a clean solve residual. This directly
fills the paper gap of applying SubtreeLU's separator private/pipeline
partition inside the pivot-aware fast path, although it is still not CKTSO's
full ETree-descendant pivoting-tail executor.

A current CKTSO-gap forced-row probe showed the closure rule enabling the
separator FLOP queue on `ASIC_320k`, `ASIC_320ks`, `ASIC_100ks`, and `rajat28`,
where earlier builds rejected the schedule because pipeline groups released
private successors. The effect was modestly positive for the ASIC cases and
neutral to slightly negative for `rajat28`, so this is best understood as
closing a paper-algorithm coverage gap rather than as the remaining large
performance lever by itself.

The separator map itself is now BTF-aware instead of depending on a unique
block-size match after KLU analysis. KLU invokes the METIS user-order callback
once per non-tiny BTF diagonal block, so KLS records each retained `NodeNDP`
component tree with that callback ordinal. After symbolic analysis returns the
accepted `R` block boundaries, KLS stitches the captured local separator trees
into one global postorder forest and synthesizes private components for blocks
that did not run `NodeNDP`. Artificial pipeline nodes connect the BTF roots so
the existing SubtreeLU Algorithm 6-style queue builder can consume one global
component map. A duplicate-size BTF smoke fixture builds two independent
30,600-row sparse-grid blocks; it requires the separator range to cover all
61,200 rows, reports at least 15 global separator components, and runs checked
row fast factorization through the separator FLOP queue with a clean solve.
This closes the earlier duplicate-block ambiguity and lets retained separator
metadata survive BTF forests, while still leaving production supernodal BLAS and
CKTSO's pivoting-tail executor as larger remaining paper gaps.

The KLS-first automatic selector now refuses one class where the implementation
was ahead of the paper coverage: very large scaled matrices whose accepted
symbolic state is one BTF block. The CKTSO paper treats this class with a
parallel row-up/ETree task factorization. KLS therefore keeps the automatic
skip for this class unless the accepted analysis retained a global separator
private/pipeline row queue that the KLS-first row-up executor can consume.
On the CKTSO-gap focus subset, `rajat24` was
the clear failure mode: before this guard it selected `kls_first`, spent about
40.1 s in initial factorization, and modeled at 66.1 s versus CKTSO's 4.68 s.
After the guard, default automatic mode reports
`kls_first_auto_skipped_scaled_single_block_count=1`, keeps the KLU/static first
factor, and models at 10.8 s. The 12-row focus geomean moved from 3.78 s
to 3.34 s, reducing the CKTSO ratio from 2.96x to 2.62x. This is not the final
paper algorithm; it prevents unpartitioned scaled single-block systems from
masking the actual missing piece, which is CKTSO's parallel row-up first factor
and the matching row-oriented refactor/solve engine.

The partial-supernode pipeline split is now a structural automatic row-engine
choice instead of an explicit experiment. With
`KLS_ENABLE_PARTIAL_SUPERNODE_PIPELINE` unset, row refactorization uses the
SubtreeLU/CKTSO split when the pipeline tail is dominated by large producer
groups; setting the variable to `0` remains a hard disable. The same
row-prefix publication is also allowed in checked row fast-factor/refactor
runs: a producer row is marked done only after its row values and pivot check
complete, and consumers still wait on each dependency row before applying the
partial compact-supernode update. This fills a direct Algorithm 5 coverage gap
for prefix-safe producer/consumer overlap, while the full CKTSO pivoting-tail
executor and production row-major first factor remain open.

A broader attempt to enable the row-dependency partial-supernode queue whenever
the structural partial-supernode test passed, even when the separator-private
side dominated the separator-pipeline side, was rejected. It matches the paper's
per-unfinished-supernode wording more directly, but it broke the current KLS
private/pipeline ready-queue invariant: smoke hung because the row-dependency
queue may start with no ready pipeline group while private work still owns the
unlocking dependencies. The current dominance guard therefore remains a
correctness guard, not a performance threshold. Closing this gap needs a queue
that lets private-group completion release row-dependency pipeline groups
without starving the ready queue, not merely relaxing the selector.

The row-level release path now matches the group-level private-queue invariant:
when a completed row drops a private successor group's dependency count to zero,
KLS leaves that successor for its private owner instead of also inserting it into
the global ready queue. This fixes a real mixed-scheduler hazard exposed while
auditing the failed relaxation. Re-running the broader selector relaxation after
this fix still timed out `test_checked_separator_flop_ready_queue`; the
debugger-owned interrupt showed workers in compact/dense group numeric
processing, not in a private-ready wait. Thus duplicate private enqueueing was a
bug, but it was not the only reason the dominance guard is still needed.

KLS then split the mixed queue's dependency accounting by owner: private groups
keep group-level predecessor counts and are released by private group
completion, while pipeline groups keep row-level predecessor counts and can be
released by completed private rows. This is closer to SubtreeLU Algorithm 6 plus
Algorithm 5 than the earlier all-row-dependency accounting, because private
subtrees remain private and only pipeline consumers participate in row-prefix
release. It passes smoke under the existing dominance guard. However, using that
hybrid queue to relax the dominance guard for unchecked row refactor still timed
out `ASIC_320k` at 120 s before producing a JSON record, so the broad automatic
selector remains rejected. The retained value is the corrected private/pipeline
accounting for the already-accepted mixed queue, not a default expansion to the
ASIC/G2 slow cases.

A follow-up unified row-dependency fallback was also rejected. The experiment
kept the mixed private/pipeline queue guarded, discarded the separator-private
phase when that guard failed, and scheduled every row group through the
row-level dependency queue so private completion could no longer starve pipeline
groups. That fixed the specific mixed-queue progress hazard, but it destroyed
the separator-private locality that Algorithm 6 exists to preserve. When enabled
for checked runs, `test_checked_separator_flop_ready_queue` no longer completed
within the 90 s smoke limit; a debugger-owned interrupt showed all workers in
row-group numeric processing rather than in the old private wait. Restricting
the fallback to unchecked row refactor made smoke pass, but the focused
`ASIC_320k` forced-row run timed out at 120 s before writing a JSON record,
versus the retained separator-private queue's roughly 40-45 s forced-row
controls. This rules out an all-groups row-dependency queue as the missing
CKTSO/SubtreeLU scheduler step. The remaining paper gap is narrower: keep the
FLOP-balanced private queues and add a release protocol where private completed
rows can unlock pipeline consumers without moving all private work into the
global ready queue.

The column EGraph refactor schedule now retains exact consecutive
supernode-candidate ranges instead of only counting them. Benchmark JSON and
gap decomposition output report the number of retained candidates and, when
enabled, the number of natural-order EGraph pipeline tasks and columns that
were coarsened. A guarded execution path exists behind
`KLS_ENABLE_EGRAPH_SUPERNODE_TASKS=1`: it only leases candidates whose width is
within the CKTSO-style `2 * threads` task-width bound, and still dispatches and
marks every column through the existing dependency and pivot checks. This is
off by default because the direct experiment showed that scheduler coarsening
without the matching supernodal numeric panel/update kernel is not the missing
CKTSO lever. On `G2_circuit` with four threads, default EGraph refactor stayed
near 0.201 s while the env-enabled coarsened path leased 6,991 tasks covering
23,328 columns and slowed to about 0.209 s. The 11 completed rows of the
CKTSO-gap focus subset were 1.0067x slower than the previous KLS run when the
task path was enabled; the `transient` timeout observed during that run was
also reproduced on the committed `d08ab35` baseline and is a pre-existing
EGraph pipeline flake. This narrows the remaining paper gap: the useful next
step is not more task scheduling, but production supernodal numeric storage
and BLAS-style panel/trailing updates for these retained ranges.

A June 28, 2026 default-path synchronization shortcut was also rejected. The
experiment retained each column's EGraph level and skipped atomic `done` waits
for dependencies known to have completed in the clustered barrier phase. It
preserved residual correctness, but the five-row CKTSO-gap focus artifact
`build/kls_skip_cluster_wait_focus5_t4_r1_ref3_timeout120.jsonl` regressed
geomean to 7.77s versus 7.59s for
`build/kls_input_auto_focus5_t4_r1_ref3_timeout120.jsonl`; four of five rows
slowed, including `G2_circuit` where there were no clustered levels to skip.
The local code was reverted. This reinforces the same paper gap as the
supernode-task experiment: scheduler-only EGraph changes are not enough without
the matching production row/segment or supernodal numeric object.

The intermittent 120 s `transient` timeout was traced to the clustered
EGraph/row-refactor barrier protocol, not to the supernode-range metadata. A
worker could observe `stop` after one level barrier and leave the clustered
phase while another worker had already entered the next level's barrier,
leaving the remaining workers asleep in `pthread_barrier_wait` and the main
thread waiting on the pool completion condition. Clustered column and row
refactor workers now drain all remaining level barriers collectively after a
stop; once stopped, they skip numeric work but still rendezvous with peers
until the level loop is complete. This preserves the existing fast
fetch-and-wait EGraph pipeline while removing the mismatched-barrier timeout
class. A 20-run `transient` stress loop with four threads completed without
timeout in default mode, with ready-queue columns reported as zero and refactor
time staying near the previous 0.026-0.030 s range.

A true column ready-queue scheduler is also available behind
`KLS_ENABLE_EGRAPH_READY_QUEUE=1`. It builds in-tail predecessor counts and
successor lists from the retained EGraph schedule and only dispatches columns
whose pipeline predecessors have completed. On `transient` it exercised all
1,037 pipeline columns and completed reliably, but repeated refactor time rose
to about 0.066-0.072 s, so it remains an off-by-default paper probe rather
than the production path. The ready-queue successor graph is now also built
only when that probe is enabled; a gated default 12-row CKTSO-gap focus run
completed with zero failures and a 3.216 s geomean, effectively unchanged from
the saved 3.206 s KLS reference while still about 2.52x slower than the saved
CKTSO reference on the same common rows. This confirms that the large CKTSO gap
is not closed by stricter task readiness alone; the next algorithmic gap is
still the numeric supernodal/panel update engine.

KLS now has a first column-EGraph numeric supernode update probe behind
`KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1`. It reuses the retained consecutive
supernode candidates, detects a contiguous dependency run, validates the dense
internal L block, solves that dependency vector through a worker-local dense
triangular panel, and only commits to the grouped path when the producer L
columns also share one trailing row list so the update can be accumulated in a
single worker-local vector before one scatter. Benchmark JSON reports
`refactor_last_supernode_update_runs`, rows, entries, and cumulative
`refactor_supernode_update_*` totals. This fills the direct SubtreeLU/CKTSO
TRSV-plus-matrix-vector update semantics in the mapped EGraph refactor kernels,
but it is still off by default because it rebuilds the dense/trailing panel for
each consumer instead of publishing a persistent producer panel.

The same-session six-row CKTSO-gap focus check makes that missing storage piece
clear. Default KLS on the first six focus rows
(`build/kls_default_same_focus6_t4_r3_timeout120.jsonl`) had an 8.24 s
SPICE-cycle geomean. The common-trailing supernode update probe
(`build/kls_egraph_supernode_updates_common_focus6_t4_r3_timeout120.jsonl`)
completed the same rows but rose to 20.36 s. It did exercise real common
trailing panels: `G2_circuit` reported 424,395 grouped updates, 5,541,687 rows,
and 1.51e9 update entries; `ASIC_100ks` reported 160,646 updates and 3.06e8
entries. The large gap therefore is not a missing task queue or lack of
supernode detection anymore; it is the paper's production compact supernode
storage/publish step, so KLS can build a producer panel once and let many
consumers reuse it instead of reconstructing the panel at every dependency run.

KLS now has that persistent producer-panel step for the column EGraph probe.
When `KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1` is set, the mapped EGraph refactor
builds a compact panel cache for retained supernode candidates whose L columns
have a dense internal block and one common trailing row list. Producer columns
publish their finalized L values into the cache before the existing
`pipeline_done` release, so downstream consumers reuse the panel after the
normal dependency wait. On the same six-row CKTSO-gap focus subset, the cached
path (`build/kls_egraph_supernode_panel_cache_focus6_t4_r3_timeout120.jsonl`)
improved the opt-in probe from 20.36 s to 10.43 s geomean, but remained slower
than the same-session default KLS geomean of 8.29 s
(`build/kls_default_panel_cache_base_focus6_t4_r3_timeout120.jsonl`). The
counters show the same high-volume updates are now covered without per-consumer
panel reconstruction: `G2_circuit` still applied 424,395 grouped updates and
1.51e9 entries in the last refactor. This fills the direct storage/publish
piece from the papers but does not close the slow-case gap by itself. The
remaining direct paper gap is compact or batched numerical kernels over these
panels, not more supernode detection.

The column EGraph cache now also consumes published supernode prefixes and
suffixes, matching the partial-publication idea already present in the
row-refactor compact panel path. The cached consumer uses `col_id` rather than
start-only lookup, solves the available dependency subrun inside the retained
dense panel, scatters any in-panel rows beyond the published prefix into the
current column work vector, and accumulates the shared trailing panel once. This
removes another direct paper gap: current columns inside a retained supernode
can reuse already published producer columns instead of rebuilding a partial
panel. On the same six-row focus subset, the corrected prefix/suffix path
(`build/kls_egraph_supernode_panel_prefix_correct_focus6_t4_r3_timeout120.jsonl`)
improved the persistent-panel probe from 10.43 s to 9.44 s geomean. The default
KLS column EGraph path remains faster at 8.29 s geomean, so this is still kept
behind `KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1`. The next gap remains a
coarser compact/batched numeric update kernel; per-consumer scalar loops over
the retained panel are not enough.

The generic mapped EGraph kernel now consumes the same cached compact panels
before falling back to its scalar dependency loop. This closes a coverage gap
left by the earlier single-block and large unscaled-BTF specializations: scaled
BTF and smaller BTF refactors can use retained producer panels when
`KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1` is set. A targeted scaled-BTF
`ckt11752_dc_1` probe reported 313 last-run compact supernode dependency
updates over 6,555 rows and 415,380 entries with the opt-in flag, while the
same command without the flag reported zero such updates. This is a direct
SubtreeLU/CKTSO compact-panel coverage step, but it still does not implement
the coarser batched numeric kernel needed to make the probe a default win.

The cached EGraph-panel consumer now has the same optional CBLAS shape as the
row-supernode experiment for sufficiently large retained panels. In CBLAS
builds with both `KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1` and
`KLS_ENABLE_CBLAS_SUPERNODE=1`, an eligible cached dependency run solves the
unit-diagonal internal panel with `dtrsv`, applies any dense in-panel suffix
with `dgemv`, and applies the shared trailing rows with one more `dgemv`
before the existing scatter. Benchmark artifacts report that BLAS-taken subset
through `refactor_last_supernode_cblas_update_*` and cumulative
`refactor_supernode_cblas_update_*` counters. This is still deliberately behind
the opt-in supernode-update probe; it makes the compact panel mathematically
closer to SubtreeLU's BLAS update shape, but does not replace the remaining
need for broader batched producer/consumer kernels or CKTSO's pivoting tail
executor.

The cached EGraph-panel consumer now also has a portable blocked panel update
when CBLAS is unavailable or `KLS_ENABLE_CBLAS_SUPERNODE=0`. Instead of walking
each producer row all the way through the dense suffix and common trailing rows,
the fallback solves a block of the unit-diagonal internal panel, applies the
remaining in-panel suffix as a dense block update, and accumulates the shared
trailing rows once per block before the existing scatter. Benchmark JSON and
stats expose this through `refactor_last_supernode_blocked_update_runs`, rows,
entries, and cumulative `refactor_supernode_blocked_update_*` counters. The new
smoke fixture forces the public mapped-EGraph path with retained producer
panels, disables CBLAS, and requires the blocked counters to match the total
cached supernode updates with a valid residual. This fills the portable
SubtreeLU/CKTSO compact-panel arithmetic gap for the column EGraph probe; the
remaining larger gap is still broader batched producer/consumer scheduling and
CKTSO's pivoting-tail executor.

A small-fragment guard for the cached EGraph-panel consumer was tested and
rejected. The probe required at least 32 contiguous producer rows and reused the
row-supernode `32768` work / `8x copied entries` threshold before entering the
cached blocked/CBLAS panel path, leaving smaller fragments on the scalar
dependency loop. This did reduce the opt-in update count substantially: on the
five current CKTSO-gap focus rows
(`build/kls_egraph_supernode_guarded_gap5_t4_r3_timeout120.jsonl`), the average
grouped update grew from roughly 7-14 rows in the earlier prefix/suffix probe to
about 24-73 rows, and the blocked subset averaged about 112-244 rows. It still
lost badly. The guarded opt-in geomean was 12.05 s versus 7.96 s for current
default KLS on the same five rows, and it was also slower than the previous
unguarded prefix/suffix opt-in artifact's 9.89 s common-row geomean. `ASIC_320k`
rose from 13.62 s default to 18.04 s guarded opt-in, and `G2_circuit` rose from
23.23 s to 34.35 s. This rules out "too many small BLAS/panel calls" as the main
remaining cause for these slow rows: once small calls are removed, the retained
panel cache build plus the remaining large blocked updates are still more
expensive than the scalar EGraph refactor. The missing paper-level piece remains
a production coarse/batched supernodal numeric executor, not a simple
small-call threshold.

A current top-ten CKTSO-gap retest after the exact row-dependency release fix
keeps the same conclusion. With `KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1`,
`build/kls_egraph_supernode_updates_gap10_after_exact_release_t4_r1_ref3_timeout120.jsonl`
completed all ten rows but measured a `4.7239s` geomean, versus `4.3966s` for
the default KLS artifact
`build/kls_current_gap10_after_exact_release_t4_r1_ref3_timeout120.jsonl`.
The probe applied real cached-panel work, for example `G2_circuit` reported
493 last-refactor grouped updates over 65,484 rows and 17.48M entries, and
`ASIC_320ks` reported 561 updates over 49,454 rows and 10.57M entries. It was
still slower than default on eight of ten rows, while only `gemat12` and
`rajat28` won. Two local blocked-kernel rewrites were also rejected rather than
committed: a row-contiguous traversal
(`build/kls_egraph_supernode_updates_rowmajor_gap10_t4_r1_ref3_timeout120.jsonl`)
and a tiled target accumulator
(`build/kls_egraph_supernode_updates_tiled_gap10_t4_r1_ref3_timeout120.jsonl`).
The tiled version improved some high-work refactor timings such as
`G2_circuit` and `ASIC_320ks` relative to the old opt-in blocked kernel, but
its top-ten geomean was still `4.7505s`, slightly worse than the old opt-in
blocked-kernel geomean and `1.0805x` slower than default KLS. The source was
therefore left on the existing portable blocked kernel.

The column EGraph ready-queue probe was also rerun after the exact-release
change. With `KLS_ENABLE_EGRAPH_READY_QUEUE=1`,
`build/kls_egraph_ready_queue_gap10_after_exact_release_t4_r1_ref3_timeout120.jsonl`
measured a `7.5726s` top-ten geomean, `1.7224x` slower than default. The queue
exercised the full recorded pipeline tails on most rows, for example 2,410
ready columns on `ASIC_320k` and 87,190 on `rajat28`, and slowed those rows by
roughly 2x and 1.7x respectively. The apparent `G2_circuit` win is not enough
to change policy because that run reported zero ready-queue columns and zero
ready-queue runs. This keeps the ready queue as an off-by-default diagnostic
and reinforces that the missing CKTSO-scale improvement is not stricter column
readiness, but a production row/supernode numeric executor with enough coarse
work to amortize its scheduling and panel-cache costs.

An interrupted `gdb` run on `G2_circuit` with a long refactor repeat confirmed
where the current default path is spending time. All four worker threads were
inside `kls_egraph_refactor_single_unscaled_column`, and the sampled instruction
addresses mapped into the inlined `kls_scatter_subtract` loop. The main thread
was only waiting for the worker pool to finish. That rules out metadata rebuild,
barrier wait, and ready-queue absence as the immediate hot-path explanation for
the current top-gap rows: the active kernel is still the sparse column scatter.

The shared scatter kernel now skips the row walk when the dependency scale is
exactly zero. This is a narrow arithmetic no-op, but it directly targets the
sampled hot loop and benefits SPICE-style repeated values that create exact
zero intermediate U entries. On the exact-release top-ten CKTSO-gap focus,
`build/kls_zero_scatter_skip_gap10_t4_r1_ref3_timeout120.jsonl` improved the
KLS geomean from `4.3966s` to `4.1762s`, a `1.0528x` speedup, with seven wins
over 2% and no losses over 2%. The largest wins were `rajat28` at `0.882x`,
`ASIC_320k` at `0.932x`, `ASIC_100ks` at `0.934x`, and `G2_circuit` at
`0.944x`. A broader first-30 focus run
(`build/kls_zero_scatter_skip_gap30_t4_r1_ref3_timeout120.jsonl`) completed
without failures and was slightly faster than the older scatter-unroll common
20-row artifact (`0.9857x` geomean), though it had mixed per-row movement. The
change is therefore retained as a real default hot-path cleanup, while the
remaining common top-ten gap is still `2.4380x` versus the saved CKTSO
four-thread reference.

Thread-count sensitivity reinforces the same point. On the exact-release
top-ten focus, one-thread KLS
(`build/kls_threads1_gap10_after_exact_release_r1_ref3_timeout120.jsonl`) had a
`9.5199s` geomean and two-thread KLS
(`build/kls_threads2_gap10_after_exact_release_r1_ref3_timeout120.jsonl`) had a
`6.5228s` geomean, both slower than the four-thread default artifact's
`4.3966s`. Eight threads
(`build/kls_threads8_gap10_after_exact_release_r1_ref3_timeout120.jsonl`)
improved all ten rows and reached `3.0889s`; sixteen threads
(`build/kls_threads16_gap10_after_exact_release_r1_ref3_timeout120.jsonl`)
improved the geomean further to `2.8736s` but regressed `rajat28`, `onetone1`,
and `gemat12` relative to eight threads. Even sixteen-thread KLS remained
`1.6776x` slower than the saved four-thread CKTSO reference on the common ten
rows. This means KLS is not losing because the four-thread default is
over-parallelized; extra workers help the large scatter-bound cases but do not
replace CKTSO's row/supernodal numeric executor.

A local attempt to reduce pipeline atomic traffic by leasing four columns per
worker was rejected. The patch preserved dependency waits and in-chunk order,
but it let a worker hold later columns while waiting on an early dependency in
the same leased chunk. The result
(`build/kls_pipeline_chunk4_gap10_t4_r1_ref3_timeout120.jsonl`) regressed the
top-ten geomean to `7.1159s`, `1.6185x` slower than the default one-column
lease. `G2_circuit` refactor time rose from about `0.213s` to `0.549s`, and the
large ASIC rows roughly doubled. The source was restored to one-column pipeline
leasing. Any future scheduler batching must be dependency-aware at the ready
task level or coupled to a real supernode task, not a blind contiguous lease.

The same conclusion now has structural evidence from the `L` scatter pattern.
KLS has an opt-in diagnostic,
`KLS_ENABLE_REFACTOR_L_PATTERN_STATS=1`, which records adjacent row-index runs
and contiguous suffixes while building the refactor LU pointer cache. On the
KLU-first CKTSO-gap path
(`build/kls_lpattern_gap5_auto_klufirst_t4_r1_timeout120.jsonl`), the four
completed slow rows had very little direct-contiguous surface: adjacent-run
coverage was 5.53% for `ASIC_320k`, 6.61% for `ASIC_320ks`, 2.58% for
`ASIC_100ks`, and 1.84% for `G2_circuit`; contiguous-suffix coverage was only
1.93%, 2.56%, 0.54%, and 0.09%, with max run length at most 5. `onetone2`
timed out in that forced KLU-first diagnostic mode, but the representative
auto/KLS-first diagnostic
(`build/kls_lpattern_gap5_auto_t4_r1_timeout120.jsonl`) showed a different
pattern family: 63-85% adjacent-run coverage but only 3.58-14.21%
contiguous-suffix coverage, and it changed the initial factor path to
`kls_first`. This rules out a BLAS-threshold guard or suffix-only direct scatter
as the broad CKTSO closer for the current KLU-first losses. A future coarse
row/supernodal engine would need to exploit general adjacent run blocks and
producer/consumer batching, not just avoid small BLAS calls.

A direct CKTSO Algorithm 5-style scalar prefactor was also tested in the
column EGraph refactor path and rejected. The probe scanned later dependencies
while an earlier predecessor was unfinished, consumed only already-published
later dependencies, and required the skipped earlier producers not to update
the candidate dependency slot. This matched the paper's "use finished
predecessors while waiting" semantics, but the safety scan was a poor fit for
the column-packed hot path. With KLS-first auto, the run also made the initial
factor choose the slower `kls_first` path
(`build/kls_egraph_prefactor_gap5_t4_r3_timeout120.jsonl`). A clean
KLU-first check (`build/kls_egraph_prefactor_klu_first_gap5_t4_r3_timeout120.jsonl`)
was still not acceptable: `ASIC_320k` was roughly neutral, `ASIC_320ks`,
`ASIC_100ks`, and `G2_circuit` regressed, and `onetone2` timed out at the
120 s harness limit. Even an off-by-default guarded version disturbed the hot
path enough to reproduce the timeout
(`build/kls_egraph_prefactor_default_off_gap5_t4_r3_timeout120.jsonl`), so the
code was removed rather than retained as another runtime flag. This rules out
that small column-EGraph Algorithm 5 gap as the next likely CKTSO closer; the
same idea remains in the row-major row-refactor executor where the dependency
and update metadata are already row-oriented.

The KLS-owned row-up first factor now executes independent BTF diagonal blocks
in parallel when KLS-first factorization is requested with more than one
thread. Each worker uses private row-up scratch, performs dynamic column pivots
inside its assigned block, and then commits the block's KLU-compatible numeric
storage and row-refactor seed data under a short mutex. This closes the direct
paper gap where the KLS-first bridge was purely serial across BTF blocks even
though CKTSO/SubtreeLU treat first factorization as a row/task parallel
problem. The smoke suite now forces a two-block BTF row-up case with two
threads and requires `kls_first_last_parallel_btf_blocks=2`. This is still not
the full CKTSO ETree pipeline inside one large diagonal block; it is the safe
BTF-level parallel part of that first-factor algorithm.

The automatic scaled single-block guard now recognizes the separator row-up
executor added after the original guard. Very large scaled single-block states
are still skipped in automatic mode when they lack retained global separator
private/pipeline work, but a METIS `NodeNDP` analysis that preserves both
private and pipeline separator rows may now enter KLS-first automatically. This
fills a direct coverage gap against the papers' row-up task-factorization
precondition without enabling the old single-block bridge broadly. The smoke
suite keeps both sides covered: a natural-order 150,000-row scaled diagonal
remains on `klu_first` with `kls_first_auto_skipped_scaled_single_block=1`,
while a 150,000-row scaled METIS/no-BTF tridiagonal uses `kls_first` and
requires the separator private and pipeline phases to execute.

The row-refactor generic segment path now direct-loads raw input rows into
retained sparse row-major segment storage when the existing row pattern can
represent the row exactly: external dependencies remain in the work vector,
internal sparse `L` entries are placed in their row mirror slots, the pivot is
placed in `Udiag`, and in-segment/trailing `U` entries are placed in the row
mirror before numeric updates. Checked runs compute the CKTSO-style row pivot
test from the direct-loaded U entries plus residual updates, so this is a
storage/algorithm bridge rather than a pivot-policy shortcut. That path now
also keeps a retained per-input-entry target map for supported dense and sparse
segment rows. The symbolic row-pattern builder classifies each raw input entry
once as an external work-vector value, row-major `L` slot, pivot, or
in-segment/shared trailing `U` slot; repeated refactors then direct-load
through those retained destinations instead of searching the segment structures
inside the numeric loop. The dense direct-input smoke case now requires target
loading over the dense segment, and the generic sparse smoke case uses a full
upper shared U segment with only sparse subdiagonal L dependencies and requires
both `row_refactor_last_sparse_segment_direct_input_rows` and
`row_refactor_last_segment_target_input_rows` to cover the whole segment. This
narrows the SubtreeLU/CKTSO row-major storage gap for generic segments, but it
still leaves the larger paper work item open: a production compact/batched
row-major numeric engine and CKTSO's pipelined pivoting-tail executor.

The retained segment plan now also owns residual cleanup destinations for
supported target rows that still need post-row input cleanup. During symbolic
row-pattern construction, KLS records the external work-vector positions from
the retained per-input-entry target map in a compact row pointer/list. Numeric
row cleanup consumes that retained list first and only falls back to the old raw
input-column scan when no retained cleanup list exists. The diagnostics now
separate planned cleanup coverage
(`row_refactor_segment_input_cleanup_rows/entries`) from cleanup actually
performed through retained targets
(`row_refactor_last_segment_target_cleanup_rows/entries` and cumulative
companions). This removes another raw KLU/CSC-style rediscovery step from the
row-segment loop, but it remains scaffolding for the paper algorithm rather
than the full production row-major engine.

The independent row-batch executors now also direct-load refactor input rows
into their native batch panels. Ragged, multi-producer, and fragmented
multi-producer batches validate that each raw input entry maps exactly to the
batch's dependency vector, pivot, or row-major U workspace, then populate those
buffers without first staging through the sparse work vector. Rows that do not
fit the current batch shape still fall back to the old residual-safe path. The
smoke probes for all three independent batch shapes now require
`row_refactor_last_batch_direct_input_rows` to cover the executed batch rows.
This fills a direct storage-path gap versus the papers' row-major numeric
updates; it is not yet a replacement for CKTSO's full pipelined pivot-tail
factorization.

The retained SubtreeLU Algorithm 6 queue splitter then gained the missing
private-leaf guard. A focused checked METIS grid exposed that the earlier
FLOP-balance loop could keep promoting separator candidates after it reached
indivisible retained components, leaving no private subtree candidates and
falling back to the generic row-DAG queue. The splitter now only promotes a
candidate component into the pipeline when it has at least one positive-work
child subtree to return to the private candidate set. The checked separator
grid and duplicate-BTF separator forest smoke fixtures both select the
separator FLOP queue again, preserving the paper's private/pipeline structure
instead of turning unbalanced leaves into separator work.

KLS tested the next direct CKTSO Algorithm 5 gap inside the checked
row-refactor pipeline: while a row blocked on an unfinished predecessor, the
worker scanned later dependencies and consumed any already-finished predecessor
whose value was provably final. The scalar and dense-producer-run variants
preserved pivot-reject order and reported through `row_refactor_prefactor_*`
and `row_refactor_prefactor_supernode_*` counters, but the broader forced-row
probes showed that this wait-time prefactor path could dominate runtime on
large dense groups. The default row-refactor partial queue now avoids that
path by using exact row-dependency release; the prefactor counters are retained
as diagnostics but the smoke fixture expects them to remain zero for the exact
partial queue. This leaves CKTSO's full ETree-descendant pivoting-tail restart
open rather than preserving a slow partial emulation.

The KLS-first separator-pipeline factorization now has the analogous grouped
executor for SubtreeLU Algorithm 4's ready-supernode branch. Earlier code
identified published row-major `U` supernode runs and counted them, but consumed
the rows by repeatedly calling the scalar dependency updater. The new
`kls_row_first_partial_apply_supernode_run` helper owns the run: it applies the
first ready dependency, lets that triangular prefix update discover the next
in-run row if needed, pops each newly ready row in order, and reports a
supernode update only when more than one row is actually consumed. It then moved
one step closer to the paper's `trsv`/`gemv` shape: when the published rows
prove a dense upper-triangular prefix plus one common trailing row list, KLS
performs the in-supernode solve over the current row and accumulates the common
trailing update in reusable worker scratch before a single scatter over the
tail. If that compact shape is not validated, the helper falls back to the old
scalar row-by-row run. This fills a numeric-kernel gap in the KLS-owned
first-factor pipeline without pretending that the full paper supernodal storage
layer is finished; the values still flow through KLS row-entry storage rather
than a persistent BLAS panel.

That KLS-first ready-supernode selector now also releases the longest
already-finished prefix of a retained row-supernode instead of requiring the
whole supernode to be below the current ready frontier. The compact helper
treats later rows of the same producer supernode as part of the common tail,
so they remain pending dependencies in the current row while the ready prefix
still contributes immediately. This maps the paper's partial producer-release
idea into the KLS-first pipeline without adding a new scheduling threshold or a
matrix-specific rule.

The same KLS-first compact/prefix supernode executor is now used by ordinary
row-up dependency loops: natural/serial rows, separator private rows, serial
pivot-tail rows, and the BTF-level parallel first-factor worker. A new
conservative run-bound check extends a dependency only through adjacent
row-major `U` rows that are already owned by the same private worker, already
published in the shared ordered suffix, or already finished inside the current
BTF worker block, and that validate as the same supernode pattern. The row then
consumes the run through the shared triangular-prefix/trailing-update helper
instead of re-entering the scalar dependency loop for each producer row. This
closes the earlier mismatch where SubtreeLU-style supernode updates existed in
the separator pipeline but not in the private-mode or BTF-parallel row factor
paths. Stats now expose actual use through
`kls_first_last_row_supernode_update`,
`kls_first_last_row_supernode_update_groups`, and
`kls_first_last_row_supernode_update_rows`, separate from the older
separator-pipeline-only counter.

KLS now fills more of that persistent-panel storage gap. KLS-first row-up
producers publish a completed row-supernode as a dense upper-triangular panel
plus common trailing column list when the next row proves that the supernode has
ended. This covers parallel BTF workers plus the ordinary private and serial
row-up loops. Later rows try that cache before falling back to row-entry
validation, and benchmark output reports actual use through
`kls_first_last_row_supernode_panel_update`,
`kls_first_row_supernode_panel_update_run_count`,
`kls_first_last_row_supernode_panel_update_groups`, and
`kls_first_last_row_supernode_panel_update_rows`. At the start of each ordered
KLS-first separator pipeline phase, the
stable private-prefix `U` snapshot is scanned for validated row-supernodes and
published as phase-local dense upper-triangular panels plus a common trailing
column list. Pipeline workers try that cache before revalidating row entries,
so repeated consumers no longer rebuild the same private-prefix panel. A
dynamic column exchange rebuilds that phase cache from the post-exchange
private-prefix `U` snapshot; row-up producer panel caches are reset because
their open producer tail layout is tied to the old column order. Benchmark
output reports actual panel-backed use through
`kls_first_last_separator_queue_pipeline_supernode_panel_update`,
`kls_first_separator_queue_pipeline_supernode_panel_update_run_count`,
`kls_first_last_separator_queue_pipeline_supernode_panel_update_groups`, and
`kls_first_last_separator_queue_pipeline_supernode_panel_update_rows`. In
CBLAS builds, the existing `KLS_ENABLE_CBLAS_SUPERNODE=1` runtime gate can
consume eligible cached KLS-first panels with CBLAS `dtrsv` over the
non-unit upper-triangular panel and CBLAS `dgemv` for dense suffix/common-tail
updates; otherwise the cached panel uses the scalar in-panel solver. The cached
consumer no longer requires the whole requested dependency run to fit inside the
same panel: it consumes the published prefix and leaves the remaining suffix in
the dependency heap, matching the papers' private/pipeline rule that completed
producer prefixes can be used before the producer tail is available. The
ordinary row-up path now also publishes a compact-validated producer prefix
lazily when a consumer first uses it, so a later consumer can reuse the same
dense/common-tail panel even before the complete producer supernode is known.
The same lazy prefix publication now applies inside separator private/pipeline
row-up paths after the compact prefix has been validated under the scoped pivot
order; dynamic pivots reset or rebuild the affected caches from the committed
post-exchange prefix so stale column layouts are not reused. This is still not
the full paper storage layer: KLS does not proactively maintain mutable
open-supernode panels, and CKTSO's full ETree-descendant pivoting-tail
scheduler remains open. Separator pipeline
pivot-tail rows that
are serialized after a restart now use the same row-up producer panel cache for
completed-supernode publication, so the restarted suffix no longer loses those
completed panels just because a pipeline phase fell back to a serial pivot row.

The BTF-level parallel first-factor worker no longer carries its own older
row-up numeric loop. It now calls the shared `kls_row_first_factor_one_row`
executor used by the serial, private, separator, and pivot-tail row paths, then
performs only the worker-specific row pointer continuation, panel publication,
separator reject accounting, and block commit. This removes a remaining
algorithmic fork: dynamic column exchanges, scoped separator pivot decisions,
owned-supernode run detection, and cached panel consumption now have one
implementation across the KLS-first row-major executors. The paper gap that
remains is the larger one already noted above: KLS still needs the full CKTSO
ETree-descendant pivoting-tail scheduler and a production compact/batched
row-major numeric storage layer.

The retained row-refactor input-target map now covers the scalar/generic row
path instead of only dense and sparse segment kernels. Pattern construction
marks single-row groups whose raw input entries can be represented by the row's
external dependency work-vector slots, pivot, and retained row-major `U` tail.
The scalar row executor consumes that map by loading pivots and `U` entries
directly into row storage, leaving only true predecessor dependencies in the
work vector, and then combines those direct values with dependency updates
before the row-wise pivot check and publish. A new smoke fixture uses a
tridiagonal chain where all rows have retained targets but only a small subset
uses the sparse-segment direct-input kernel, proving the generic row path uses
the retained destinations too. This is another storage-path step toward the
CKTSO/SubtreeLU row-major numeric engine; it still does not replace the larger
compact/batched engine or CKTSO's ETree-descendant pivoting-tail scheduler.

The compact dense row-refactor panel loader now consumes the same retained
per-input-entry target map. When the symbolic row target data is available, raw
input values are placed directly into the compact row's external work-vector
slots, dense `L` panel, pivot, dense in-panel `U`, or trailing `U` panel instead
of rescanning trailing columns during numeric refactorization. The compact-panel
smoke coverage now requires both compact direct-input rows and retained target
consumption, so this closes another storage-path mismatch against the
CKTSO/SubtreeLU row-major algorithm description. The larger open paper gap is
unchanged: KLS still needs a production compact/batched row-major numeric engine
and CKTSO's pipelined ETree-descendant pivoting-tail executor.

The checked-reject KLS-owned block repair now attempts the restartable
Algorithm 5-style row pipeline first when multiple threads are requested, before
the older serial pivot-tail repair and KLU fallback. This path still works at
BTF-block repair scope rather than CKTSO's exact ETree-descendant tail mask, but
it removes a direct mismatch in the fast-reject recovery flow: a failed checked
factor can continue through the row-up pipeline and, if needed, serialize only
dynamic-pivot epochs inside that pipeline. Benchmark and
decomposition output expose the behavior through
`fast_kls_block_restart_last_row_pipeline`,
`fast_kls_block_restart_row_pipeline_count`,
`fast_kls_block_restart_last_row_pipeline_rows`,
`fast_kls_block_restart_last_row_pipeline_threads`,
`fast_kls_block_restart_last_row_pipeline_prefix_rows`,
`fast_kls_block_restart_last_row_pipeline_suffix_rows`,
`fast_kls_block_restart_last_row_pipeline_gap_rows`,
`fast_kls_block_restart_last_row_pipeline_pivot_tail_rows`,
`fast_kls_block_restart_last_row_pipeline_pivot_restarts`,
`fast_kls_block_restart_last_row_pipeline_supernode_update_groups`,
`fast_kls_block_restart_last_row_pipeline_supernode_update_rows`,
`fast_kls_block_restart_last_row_pipeline_supernode_panel_update_groups`, and
`fast_kls_block_restart_last_row_pipeline_supernode_panel_update_rows`; the
prior-pivot smoke fixture now requires the two-thread repair to use that
pipeline.

KLS now narrows that block-repair pipeline when the retained pivoting-tail plan
is an exact topological suffix ending at the end of the rejected BTF block. In
that case, the repair seeds the preserved prefix from the previous LU into
row-major entries, marks those prefix rows finished, and starts the restartable
row pipeline at the rejected suffix rather than at row zero. The activation is
deliberately conservative: if the prefix has non-identity row order or the
retained plan is a non-contiguous/shorter ETree mask, KLS leaves the repair on
the existing serial exact-mask path where preserved later rows are explicitly
validated. The BTF suffix smoke case now requires a nonzero preserved-prefix
count, showing that the multi-thread repair is no longer only a full-block
pipeline in exact-suffix cases.

The same row-pipeline repair now handles the next conservative CKTSO tail case:
a contiguous, topological pivoting-tail envelope that stops before the end of
the rejected BTF block. KLS builds the retained tail mask, refreshes preserved
non-tail block columns just like the serial exact-mask path, seeds both the
preserved prefix and suffix into row-major storage, and runs the restartable
pipeline only on the active envelope. During that envelope, dynamic pivot
selection is bounded to active tail columns; if a row would need a preserved
suffix column to satisfy the pivot check, the pipeline attempt fails and the
existing serial/KLU fallback ladder remains responsible. The root-independent
tail smoke fixture now runs with two threads and requires two active pipeline
rows plus one preserved suffix row for both unscaled and scaled repairs.

KLS now fills the next direct CKTSO tail gap by admitting non-contiguous
topological pivoting-tail masks into the same restartable row pipeline. The
pipeline now consumes the retained ETree-descendant tail worklist as its
compact topological row order, keeps preserved gap and suffix rows seeded from
the prior LU, and applies readiness through an active-rank map instead of raw
row-number thresholds. Suffix and contiguous retained tails also enter through
that worklist path instead of the earlier broad interval shortcut, so successful
KLS-owned block repairs now execute CKTSO Algorithm 5's retained restart-node
sequence inside the row-first pipeline before using the existing fallback ladder
for unsafe pivot exchanges.
Public and benchmark statistics now report
`fast_kls_block_restart_last_row_pipeline_gap_rows`,
`fast_kls_block_restart_last_row_pipeline_etree_tail`,
`fast_kls_block_restart_row_pipeline_etree_tail_count`,
`fast_kls_block_restart_last_row_pipeline_etree_tail_rows`,
`fast_kls_block_restart_last_row_pipeline_etree_tail_gap_rows`, and
`fast_kls_block_restart_last_row_pipeline_etree_tail_exact_mask`; the later
gapped-mask prefactor executor adds
`fast_kls_block_restart_last_row_pipeline_etree_ready`,
`fast_kls_block_restart_row_pipeline_etree_ready_count`,
`fast_kls_block_restart_last_row_pipeline_etree_ready_rows`,
`fast_kls_block_restart_last_row_pipeline_etree_ready_threads`,
`fast_kls_block_restart_last_row_pipeline_etree_prefactor`,
`fast_kls_block_restart_row_pipeline_etree_prefactor_count`,
`fast_kls_block_restart_last_row_pipeline_etree_prefactor_rows`,
`fast_kls_block_restart_last_row_pipeline_etree_prefactor_threads`,
`fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_rows`, and
`fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_deps`. Descendant
pivot epochs inside that ETree-ready branch reuse
`fast_kls_block_restart_last_row_pipeline_pivot_tail_rows` and
`fast_kls_block_restart_last_row_pipeline_pivot_restarts`; the non-contiguous
gap smoke fixture requires a preserved gap row and no serial tail restart when
two threads are enabled.

The retained-tail row-pipeline repair now also validates SubtreeLU's component
pivot domain before accepting the active ETree mask. If separator analysis
covers the repaired BTF block and the retained ETree tail is a suffix-exact
restart, every retained active row must map to a valid collapsed-separator
component and dynamic pivot choice remains bounded by that component extent
intersected with the active tail mask; otherwise KLS falls back to the existing
separator-queue, serial, or KLU repair ladder. A direct non-suffix experiment on
the METIS tridiagonal fast-repair smoke case was rejected: the active tail
covered only 5,731 of 30,000 rows, passed the local row pipeline, but left later
weak pivots outside the retained ETree tail and fell through to KLU. KLS now
pre-validates active-mask repairs with non-recording numeric and rowwise-U pivot
checks before committing them, and keeps non-suffix separator-covered rejects on
the full separator queue path. Benchmark output records this direct
paper-algorithm overlap with
`fast_kls_block_restart_last_row_pipeline_separator_tail_scope`,
`fast_kls_block_restart_row_pipeline_separator_tail_scope_count`, and
`fast_kls_block_restart_last_row_pipeline_separator_tail_scope_rows`.

For separator-covered rejects where the repaired block can use the whole
retained separator tree, KLS now tries the Algorithm 6-style split more
directly instead of immediately reducing the repair to the older retained
ETree-tail active mask. The block-local row and column order is remapped into
private rows followed by collapsed separator-pipeline rows, the private
ownership proof is checked against the remapped row dependencies, validated
private rows execute in parallel worker-local entries, and the pipeline roots
finish through the restartable pivot-capable row pipeline. If the private proof
or private phase is not safe, the same remapped block remains eligible for the
ordinary full row-pipeline repair rather than accepting unsafe private work.
The smoke suite now forces this with a METIS-ordered tridiagonal fast-factor
reject and requires the partitioned separator queue counters. Benchmark output
records the accepted full-block separator repair through
`fast_kls_block_restart_last_row_pipeline_separator_queue`,
`fast_kls_block_restart_row_pipeline_separator_queue_count`,
`fast_kls_block_restart_last_row_pipeline_separator_private_rows`,
`fast_kls_block_restart_last_row_pipeline_separator_pipeline_rows`,
`fast_kls_block_restart_last_row_pipeline_separator_private_threads`,
`fast_kls_block_restart_last_row_pipeline_separator_partitioned`, and
`fast_kls_block_restart_last_row_pipeline_separator_split_components`.

The serial fallback now consumes the same retained topological tail envelope
instead of forcing the failed pivot to be the restart boundary. When the
retained ETree-descendant plan is non-suffix, KLS builds the exact active mask
over that plan's first/last columns, reconstructs the current prefix only before
the retained tail begin, and records/refreshes accepted off-diagonal data from
that actual begin. This matters for unfinished parallel fast-factor states whose
retained descendants can begin before the pivot that finally failed. The
one-thread non-contiguous smoke variant requires the serial exact-mask restart,
while the two-thread variant still requires the row-pipeline path. This closes a
serial wrapper mismatch; the remaining paper gap is still CKTSO's full
multi-task ETree-descendant pivoting-tail scheduler and production row/segment
numeric storage.

The masked row-pipeline repair now makes the same ready-supernode branch
mask-aware. Packed active tail rows no longer force singleton producer rows:
the phase initializes supernode metadata from already-preserved rows, each
ordered active-row publish rebuilds the completed local supernode map across
preserved gaps and active rows, and readiness checks use the active-rank map to
stop a run at unfinished active rows while still admitting ready preserved gap
rows. The phase-local dense/common-tail panel cache now uses the same completed
row map in masked tails: preserved-row panels are built at phase start, active
row publication rebuilds the cache across preserved gaps and completed active
rows, and dynamic column exchanges refresh against the current row-major `U`
before speculative suffix rows retry. This closes the direct Algorithm 4/5
panel-use gap for completed dependencies in non-contiguous tail masks, while
still stopping at unfinished active rows through the active-rank guard. The fast
block-repair pipeline reports this through
`fast_kls_block_restart_last_row_pipeline_supernode_update_groups`,
`fast_kls_block_restart_last_row_pipeline_supernode_update_rows`,
`fast_kls_block_restart_last_row_pipeline_supernode_panel_update_groups`, and
`fast_kls_block_restart_last_row_pipeline_supernode_panel_update_rows`; the
non-contiguous gap smoke fixture now requires one grouped update and one
panel-backed update over two producer rows.

The same KLS-owned row-first tail repair now keeps the accepted off-diagonal
refresh inside the proven tail scope. After accepting a row-first repaired
block, KLS compares the repaired prefix `Pnum` and block-local `Q` order against
the saved prefix. When that prefix is unchanged, it refreshes `Offp`/`Offi`/
`Offx` only from the rejected pivot onward; if the proof or suffix refresh
fails, it falls back to the existing full offdiag rebuild. This applies the
same CKTSO-tail storage boundary that the serial tail restart already used to
the parallel row-first repair path. The non-contiguous mask and BTF suffix smoke
fixtures now require `fast_repaired_last_offdiag_suffix_refresh=1`, one suffix
refresh, and zero full offdiag refreshes.

The retained row-refactor input-target map now covers independent batch groups
as well as dense, sparse, compact-panel, and scalar/generic rows. Batch-group
symbolic setup classifies each raw input entry once as a batch dependency-vector
slot, pivot, or row-major `U` workspace slot. Ragged, contiguous multi-producer,
and fragmented multi-producer batch executors consume those destinations before
falling back to the older per-refactor pattern search. The existing batch smoke
fixtures now require `row_refactor_last_segment_target_input_rows` to cover the
executed batch rows in both checked and unchecked modes. This removes another
KLU-shaped rediscovery step from the KLS row/segment numeric path, while the
larger paper gaps remain the production compact/batched row-major numeric engine
and CKTSO's full ETree-descendant pivoting-tail scheduler.

The exact-pattern compact dense batch executor now consumes the same retained
compact-panel input targets before applying producer batches. Exact dense
consumer row batches no longer stage raw input through the sparse work vector
only to fill the compact panel; raw external dependencies still enter `x`,
while in-panel `L`, pivot, dense-`U`, and trailing-`U` values are direct-loaded
into the retained panel and combined with later residual updates. The exact and
subrange compact-supernode smoke probes require compact direct-input and
retained target counters to cover executed batch rows. This closes another
storage-path mismatch; the larger open items remain production batched
row-major kernels and the full CKTSO tail scheduler.

The same retained compact-panel input contract now covers ragged and fragmented
dense-consumer batches. The single-producer ragged path and the fragmented
multi-producer dense-group path both direct-load raw `L`, pivot, dense-`U`, and
trailing-`U` entries into the retained row-major panel, then merge any residual
work-vector value for that slot before applying producer suffix updates. Their
smoke fixtures now require compact direct-input rows and retained target-input
rows to cover the executed compact-supernode batch rows, including the checked
fragmented target-map case. This removes another KLU-shaped staging rule from
the SubtreeLU-style compact batch executors without claiming to close the
remaining production-kernel and CKTSO-tail work.

Compact dense row-refactor groups now have a portable blocked panel factor path
even when CBLAS is disabled. After the retained panel is loaded, KLS factors
diagonal panel blocks row-by-row, solves the below-block multipliers against the
just-factored upper block, and applies one dense/trailing right-looking update
before publishing the row-major values. The checked/pivoting variant preserves
the existing multiplier and pivot acceptance tests, publishes each accepted
diagonal row as soon as that row is final, and leaves fast-reject prefix
semantics intact. This fills a direct SubtreeLU/CKTSO algorithm gap: compact
supernodal panel arithmetic no longer exists only as a CBLAS probe or as per-row
scalar updates. Benchmark JSON and stats expose this through
`row_refactor_last_compact_dense_panel_blocked`,
`row_refactor_compact_dense_panel_blocked_run_count`,
`row_refactor_compact_dense_panel_blocked_rows`, and
`row_refactor_compact_dense_panel_blocked_entries`; both unchecked and checked
dense compact smoke coverage force CBLAS off and now require the blocked path.
The larger open gaps remain the production checked/pivoting tail executor and
broader row-major numeric engine.

The blocked dense-panel path is now promoted behind the native row-panel
selector rather than being treated only as an unnamed compact-panel probe.
`KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR` is tri-state: unset keeps the structural
work gate, `1` forces retained panel use for eligible dense row groups, and `0`
uses the direct row-major dense-group fallback. New counters
`row_refactor_native_row_panel_*` distinguish selected native panels, blocked
panel executions, unsupported-group fallbacks, and checked-pivot rejects. This
does not add CPU-specific tuning; it exposes the paper-shaped retained-panel
executor as the default candidate for CKTSO-gap refactor cases while preserving
the scalar fallback.

KLS then removed a non-paper threshold from the scalar compact-supernode update
selector. SubtreeLU's row update branch treats a ready supernode as a triangular
solve plus trailing update; the previous default only used KLS's contiguous
worker-scratch `trsv` when a work-per-copied-entry gate said it would amortize.
That broader default was rechecked after the CKTSO-gap rerun because the forced
row path was executing tens of thousands of compact triangular solves with
average producer runs below 80 rows. On the top-ten CKTSO-gap forced-row
comparison, the unset default measured `12.2795s` SPICE-cycle geomean while
`KLS_ENABLE_COMPACT_SUPERNODE_TRSV=0` improved to `11.2955s`
(`build/kls_native_trsv_on_forced_row_gap10_t4_r1_ref3_timeout120.jsonl`,
`build/kls_native_trsv_off_forced_row_gap10_t4_r1_ref3_timeout120.jsonl`).
A follow-up large-run gate at the CBLAS vector scale still did not beat the
full off-switch: it measured `11.6484s` geomean and the few remaining 512+
row triangular solves still hurt rows such as `onetone1`
(`build/kls_native_trsv_largegate_forced_row_gap10_t4_r1_ref3_timeout120.jsonl`).
The final decision used a higher-repeat top-five guard on the dominant slow
rows: the auto-off default measured `16.5528s`, while force-on measured
`16.8222s`
(`build/kls_native_trsv_autoff_forced_row_gap5_t4_r1_ref20_timeout120.jsonl`,
`build/kls_native_trsv_forceon_forced_row_gap5_t4_r1_ref20_timeout120.jsonl`).
The unset default therefore keeps the triangular part on the scalar dependency
walk while retaining the contiguous trailing accumulation. This rejects the
small-BLAS hypothesis as a CKTSO-gap closer for the current row scaffold rather
than the paper algorithm itself. `KLS_ENABLE_COMPACT_SUPERNODE_TRSV=1` remains
the explicit force-on coverage/probe mode, and the compact-panel smoke test uses
that force-on path when it requires the
`row_refactor_last_compact_supernode_trsv` counters to fire.

KLS-first row-up and pivot-tail pipeline phases now also have a portable
cached-supernode panel executor. When a published U-row run has been retained as
a dense panel, the non-CBLAS path solves all ready predecessor multipliers as one
upper-triangular run and applies the common dense/trailing update in aggregate
before falling back to the older per-row panel walk. This moves the
SubtreeLU/CKTSO row-first path closer to the papers' supernode update branch
without claiming the full pipelined ETree-descendant tail factorization is
complete.

The row-first supernode consumer now publishes a compact-validated producer run
before its first consumer falls back to the scalar/compact walk. That makes the
same first consuming row use the cached dense/common-tail panel solve/update
path when the ready run already has the SubtreeLU-style panel shape; the new
single-consumer KLS-first smoke fixture requires exactly one panel-backed update
over the two ready producer rows. This fills the first-consumer half of the
paper ready-supernode branch. The remaining row-major gap is still larger:
production coarse supernode storage/kernels and CKTSO's full pivoting
ETree-tail executor are not complete.

The row-first failed-pivot repair then stopped requiring the retained
pivoting-tail plan to start exactly at the rejected pivot before the
restartable pipeline could use it. KLS now derives a local repair envelope from
the retained topological tail's first and last columns, validates that the
envelope contains the rejected pivot, and either runs the contiguous envelope or
uses the active tail mask for non-contiguous gaps. The reusable-prefix proof,
off-diagonal suffix refresh, and repaired-tail work accounting all start at the
actual retained-tail begin rather than blindly at the rejected pivot. The
non-contiguous gap smoke fixture now requires an exact repaired tail mask, zero
overcompute columns, and suffix-only off-diagonal refresh. This is a direct
CKTSO Algorithm 5 gap closure for the current KLS-owned row-first repair path;
it still is not the full parallel ETree-descendant pivoting-tail scheduler.

KLS now fills another CKTSO Algorithm 5 semantic gap in how that retained tail
is seeded. The CKTSO paper says that, after a pivot-check interruption, the
tail starts from the unfinished guessed-EGraph nodes and their ETree closure.
KLS previously refreshed missing prefix rows first and then built the seed from
the remaining suffix or a narrower row-tail heuristic, which could erase the
actual interrupted unfinished set. The checked row and mapped EGraph reject
paths now snapshot the whole rejected block's unfinished done-bitmap before
prefix refresh, guarantee that the rejected pivot is in the seed, allow seed
columns before the rejected pivot, and prefer that seed before row-tail/suffix
fallbacks. The parallel checked-row smoke fixture now requires the pivoting-tail
seed to come from that unfinished set, with row-tail seeding bypassed. This
still leaves the larger CKTSO executor gap open: the retained worklist is used
by KLS's guarded row-first repair path, not by a full production tail scheduler.

KLS now consumes that retained worklist more directly whenever the retained
topological ETree tail has more than one active row. The row-first repair first
factors the boundary pivot row with the existing pivot-capable row kernel, then
removes that completed boundary from the active mask and runs the remaining
active descendants through the retained row pipeline with an active-rank map
over the ETree tail. If the retained mask covers the whole local repair
envelope, preserved-column refresh is now accepted as a successful no-op, so a
full unfinished block does not fall back to the generic full-block row pipeline
before the retained ETree executor can run. The worker applies already finished
dependencies before the row is publishable, waits for earlier active tail rows,
and then applies the skipped dependencies before pivot check and publication. If
a descendant row still needs a dynamic column exchange, the ETree-prefactor
phase treats that as a restartable pivot epoch: it serializes the pivot row with
the same
row-up-looking kernel, resets the row-up producer panel cache for the new column
order, and resumes the remaining active descendants under the same retained
ETree rank map instead of falling back out of the ETree-ready executor. The
smoke fixtures now require the ETree-ready path for both a contiguous
checked-row unfinished tail and a non-contiguous gap tail; the gapped fixture
also requires `fast_kls_block_restart_last_row_pipeline_etree_prefactor=1`, and
a separate ETree-ready descendant pivot fixture requires both the ETree-ready
counters and nonzero
`fast_kls_block_restart_last_row_pipeline_pivot_tail_rows`/restart counters.
This fills a concrete CKTSO Algorithm 5 prefactor/postfactorization executor
gap for active unfinished descendants, while singleton tails, suffix-shaped
repairs without retained descendant work, and BTF suffix repairs intentionally
remain on the ordered pivot-capable row pipeline until the ETree scheduler
handles those cases without changing pivot semantics.

KLS now also fills the NICSLU Algorithm 4 performance-model gap. After a
numeric factorization, it walks the actual U-dependency graph in factor order,
assigns each column to the thread with the smallest current `END(p)`, charges
`2*nnz(L(:,i))` work for each dependency update, `nnz(L(:,k))` work for
normalization, and a unit sync cost per dependency, then reports the predicted
finish time and flop-only speedup. The model is exposed through
`parallel_task_flow_*` stats and benchmark JSON, and can seed the same
row/segment metadata preparation path as the earlier R1/R2 NICSLU suitability
counters. This is still a model and policy input, not NICSLU's full
ETree/EScheduler-guided factorization executor.

KLS-first pivoting row-up factorization now consumes the retained separator
tree through the same Algorithm 6-style split/collapse idea used by the
row-refactor separator queue. The first-factor planner computes block-local
row-input work per retained component, splits dominant subtrees into pipeline
roots plus private child-subtree candidates, assigns those candidate subtrees
to private thread queues by work, and validates private ownership against the
original row dependencies before remapping the block. Invalid partitions fall
back to the older retained-component queue instead of forcing the generic row
pipeline. Benchmark JSON and smoke coverage expose the path through
`kls_first_last_separator_queue_partitioned`,
`kls_first_separator_queue_partitioned_count`, and
`kls_first_last_separator_queue_split_components`. This closes the direct
KLS-first gap where the pivoting first factor had Algorithm 3 execution but not
the paper's Algorithm 6 separator split. Separator-covered fast-factor block
repair now uses the same split for full-block rejected repairs, with dedicated
`fast_kls_block_restart_*_separator_queue` counters. General refactor pivoting,
the complete CKTSO checked-tail scheduler around guessed ETree interruption,
and production coarse supernode storage remain open.

The default EGraph refactor path now keeps the opt-in supernode numeric update
probe out of the scalar dependency loop unless
`KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1` is set. Earlier CBLAS probes on the
current CKTSO-gap focus rows recorded zero CBLAS update counters, because the
existing BLAS gates already require large row/update shapes. The measurable
default overhead was instead the disabled supernode hook itself: every U
dependency called the probe only to return immediately. Hoisting the
`supernode_numeric_updates` guard around the dependency probe and publish hook
keeps the experimental supernode/CBLAS behavior available when requested, while
the normal scalar EGraph kernel stays branch-light. A top-ten CKTSO-gap
repeat-3 check improved the KLS cycle geomean from `4.3115s` to `4.2172s`
against `build/kls_noauto_klsfirst_gap10_t4_r3_timeout120.jsonl`, with seven
wins over 2%, two ties, and one small absolute loss on `gemat12`. This is a
cleanup around a disabled paper-level experiment, not the missing CKTSO-scale
storage/executor change: the same top rows remain about `2.54x` slower than the
saved CKTSO medium artifact.

A follow-up CBLAS guard check confirmed that small BLAS calls are not the
default CKTSO-gap cause. KLS's CBLAS paths are still build-time optional,
runtime opt-in through `KLS_ENABLE_CBLAS_SUPERNODE=1`, and gated by structural
row/panel sizes plus minimum estimated work before calling `dtrsv`, `dgemv`,
`dtrsm`, or `dgemm`. A same-binary current-source top-ten guard with the CBLAS
runtime gate on versus off moved from `4.8698s` to `4.7526s` geomean on one
repeat-1/refactor-3 pass, but the small `gemat12` row still regressed and an
older repeat-3 artifact comparing the optional CBLAS build against the normal
default remained worse (`4.5385s` versus `4.2172s`). That makes CBLAS useful as
an opt-in probe, not a default policy. The production default should keep using
KLS-owned scalar/blocked panel kernels until a broader row-major storage engine
can feed BLAS-size work without extra staging overhead.

A fresh same-binary spot check on June 27, 2026 confirms the guard is active in
current source. On `onetone2` with the production row refactor disabled, the
CBLAS-enabled build reported zero `refactor_*_cblas_update_*` counters whether
`KLS_ENABLE_CBLAS_SUPERNODE` was `0` or `1` (`0.01426s` versus `0.01356s`
refactor time). On the forced row-refactor path with the known partial
supernode pipeline disabled, `gemat12` had no compact supernode updates and
`onetone2` used the large compact update path without a meaningful CBLAS
penalty (`0.15971s` off versus `0.15829s` on, 4,191 compact-supernode updates
covering 172M update entries). The remaining slow forced-row issue is therefore
not small BLAS calls; it is still the row ready-queue/partial-supernode pipeline
pathology exposed by the default forced-row timeout.

The same conclusion held on the current top-ten CKTSO-gap focus after the exact
row-dependency release fix. A separate `build-cblas` binary configured with
`-DKLS_ENABLE_CBLAS_SUPERNODE=ON` was run with `OPENBLAS_NUM_THREADS=1` to keep
BLAS internal threading out of the comparison. With the runtime gate off,
`build/kls_cblas_build_gate_off_gap10_t4_r1_ref3_timeout120.jsonl` measured a
`4.2390s` geomean; with `KLS_ENABLE_CBLAS_SUPERNODE=1`,
`build/kls_cblas_build_gate_on_gap10_t4_r1_ref3_timeout120.jsonl` measured
`4.2850s`, a `1.0108x` slowdown. More importantly, every default EGraph row in
that run reported zero `refactor_last_supernode_cblas_update_*` and
`refactor_supernode_cblas_update_*` counters even when the runtime gate was on.
So the proposed "BLAS only for large cases" guard is already present for the
production path: the CKTSO-gap losses are occurring before KLS reaches any
default BLAS call. The open paper gap remains durable supernode/panel storage
and an executor that creates BLAS-sized reusable work without per-consumer
staging, not an unguarded small-BLAS threshold.

The EGraph producer-panel cache now also records how many retained panels are
structurally consumable by the actual U dependency streams. The build step
scans the same contiguous dependency-run conditions used by the cached
consumer, exposes `refactor_supernode_panel_count` and
`refactor_supernode_panel_used_count`, and disables panels that no downstream
consumer can use. The June 28, 2026 top-ten opt-in run
(`build/kls_egraph_supernode_pruned_gap10_t4_r1_ref3_timeout120.jsonl`) did
not find dead cache storage: every built panel was usable (`6913/6913` across
the ten matrices), and the geomean was `4.8018s`. The default control
(`build/kls_default_after_panel_prune_gap10_t4_r1_ref3_timeout120.jsonl`)
reported zero panel-cache counters and a `4.0077s` geomean, matching the prior
guarded default within run noise. This rules out unused retained-panel
selection as the large CKTSO gap. The remaining paper-level missing part is
the production row/supernode numeric layout and scheduler that amortizes panel
publication and consumes many rows per stored panel, not another BLAS-size
threshold or a dead-panel prune.

The same conclusion held after splitting the lightweight scalar EGraph
supernode-run consumer away from the retained-panel cache as a local probe.
With scalar supernode updates enabled by default but without building the
cached panel object, the current top-ten CKTSO-gap focus regressed to
`6.8753s` geomean (`build/kls_scalar_supernode_auto_gap10_t4_r1_ref3_timeout120.jsonl`)
against `4.2018s` with `KLS_ENABLE_EGRAPH_SCALAR_SUPERNODE_UPDATES=0`
(`build/kls_scalar_supernode_off_gap10_t4_r1_ref3_timeout120.jsonl`). Adding a
minimum four-column run precheck still measured `6.6764s` versus a disabled
control at `4.1949s`
(`build/kls_scalar_supernode_guarded_auto_gap10_t4_r1_ref3_timeout120.jsonl`,
`build/kls_scalar_supernode_guarded_off_gap10_t4_r1_ref3_timeout120.jsonl`).
The guarded probe recorded zero accepted grouped updates on `ASIC_320k`,
`ASIC_320ks`, `onetone2`, `ASIC_100ks`, `G2_circuit`, and `rajat28`, yet still
roughly doubled several refactor times just by testing candidate starts that
failed the common-trailing validation. Accepted work on the remaining rows was
too small to help: `onetone1` had one run over 67 rows and 2,278 entries,
`transient` had 13 runs over 52 rows and 130 entries, and `rajat24` had nine
runs over 36 rows and 90 entries. This rejected the idea that the existing
column-storage scalar supernode detector only needed to be decoupled from the
panel cache.

The L-column segment statistics tell the same story for a direct segment
scatter shortcut. With `KLS_ENABLE_REFACTOR_L_PATTERN_STATS=1`, `onetone2`
reported 49,349 adjacent-run entries out of 537,026 L entries with max run
length 10; `G2_circuit` reported 118,454 out of 6,446,730 with max run length
5; and `ASIC_100ks` reported 40,190 out of 1,556,952 with max run length 4.
Those runs are too short and sparse to plausibly close a 2-3x CKTSO gap via a
specialized scatter loop. The paper-aligned missing piece remains a production
row/segment-oriented numeric object and executor, not more discovery on top of
KLU-compatible column storage.

Callgrind then isolated the remaining hot default path on `onetone2`: the
release-with-debug sample charged about 35.7% of instructions to
`kls_egraph_refactor_btf_unscaled_column`, with the indirect
`kls_scatter_subtract` update lines dominating inside that kernel. The scalar
scatter loop now uses an eight-entry unroll before the existing four-entry tail
loop. This is not a CPU-specific backend substitution; it reduces loop overhead
in KLS's own sparse scatter kernel. Two top-ten CKTSO-gap repeat-3 guards moved
from the retained `4.2172s` geomean to `4.0588s` and `4.1732s`, and a top-20
one-pass guard completed all rows. The improvement is useful but still small
relative to CKTSO: the focused common-row gap remains about `2.51x`, so the
larger paper gap is still coarse row/supernode storage and scheduling, not this
scatter loop alone.

Several adjacent probes were rejected in the same inspection pass. Retaining
per-column cluster level and skipping waits for cluster dependencies regressed
the top-ten guard by about 1%. Reusing reciprocal pivots inside the EGraph
kernel regressed by about 0.4%. Hoisting the supernode-update branch outside the
dependency loop was unstable, with one apparent win and one repeat loss. Forcing
the KLS-first path timed out on `rajat24` and was much worse than the auto
policy. These failures are kept out of the default path because they do not
close the algorithmic gap described in the CKTSO/SubtreeLU papers.

A June 28, 2026 same-shape rerun confirmed that rejection on the current
source. The stable default top-ten CKTSO-gap run
(`build/kls_current_default_gap10_t4_r1_ref3_timeout120.jsonl`) completed all
ten matrices with a `4.1846s` SPICE-cycle geomean. Forcing the current
KLS-first scaffold
(`build/kls_current_klsfirst_gap10_t4_r1_ref3_timeout120.jsonl`) completed
only nine matrices, timed out `rajat24` at 120s, and was `1.184x` slower than
the default on the nine common successes. The loss was cold-factor dominated:
`ASIC_320k` initial factor time moved from `1.242s` to `3.987s`,
`ASIC_320ks` from `0.994s` to `2.452s`, `G2_circuit` from `1.055s` to
`2.649s`, and `transient` from `1.005s` to `3.011s`, while repeated EGraph
refactor times barely moved. The same default run remains `2.443x` slower than
the retained CKTSO top-ten reference
(`build/cktso_paper_medium93_t4_timeout120.jsonl`). This narrows the useful
next step: not KLS-first policy gating or small BLAS thresholds, but the larger
paper gap in production row/supernode numeric storage and the repeated
refactor executor fed by that storage.

The checked row fast-factor executor is now explicit opt-in for the same
reason. A current top-ten CKTSO-gap A/B with
`KLS_ENABLE_CHECKED_ROW_REFACTOR=0`
(`build/kls_checked_row_off_gap10_t4_r1_ref3_timeout120.jsonl`) reduced
repeated checked `kls_factor` time sharply on matrices that had been attempting
and then abandoning the checked row path: `ASIC_320k` moved from `0.5398s` to
`0.1090s`, `ASIC_320ks` from `0.5760s` to `0.0880s`, `ASIC_100ks` from
`0.4459s` to `0.0581s`, and `G2_circuit` from `1.1145s` to `0.2444s`;
the common-row repeated-factor geomean moved to `0.4517x` of the
zero-scatter-skip baseline.
However the project SPICE-cycle score is dominated by initial factorization and
99 no-pivot refactors, so the same source-default opt-in run
(`build/kls_checked_row_optin_gap10_t4_r1_ref3_timeout120.jsonl`) was neutral
against the zero-scatter-skip baseline: `4.1888s` versus `4.1762s` geomean.
That keeps checked row execution available for targeted paper tests and future
row-major work, while preventing the default production `kls_factor` path from
paying row-pattern setup for an executor that is not yet a broad CKTSO-gap win.

The follow-up BLAS-threshold question exposed a benchmark fairness issue rather
than a new solver algorithm: `cktso_compare` drives CKTSO with 32-bit CSC input
arrays, while `kls_bench` was always calling the public KLS API through
`KLS_INDEX_INT64`. `kls_bench` now defaults to `--input-index auto`, converts
the parsed MatrixMarket structure to 32-bit CSC arrays when the order and nnz
fit, reports `requested_input_index` and `input_index_bytes` in JSON/text, and
falls back to 64-bit for larger inputs. `run_bench_suite.py` forwards
`--input-index auto|32|64` so current paper-gap reruns can compare the fair
32-bit-input path against the old forced-64 harness. This does not close the
paper algorithm gap by itself; it removes an input-width artifact before
judging the remaining EGraph/row-supernode executor losses.

The rerun confirmed that size-gated BLAS and input width are not the large
missing piece. On the five hardest CKTSO-gap rows, auto 32-bit input measured
`7.5867s` geomean versus `7.7692s` for forced 64-bit input
(`build/kls_input_auto_focus5_t4_r1_ref3_timeout120.jsonl`,
`build/kls_input_64_focus5_t4_r1_ref3_timeout120.jsonl`), while the saved CKTSO
artifact is `3.0464s`. On the top-ten gap set, auto measured `4.0602s` versus
`4.0805s` forced-64, still `2.37x` the saved CKTSO geomean
(`build/kls_input_auto_gap10_t4_r1_ref3_timeout120.jsonl`,
`build/kls_input_64_gap10_t4_r1_ref3_timeout120.jsonl`). Every top-ten auto
row reported `input_index_bytes=4`, `build_has_cblas=false`, and zero CBLAS
supernode-update counters. The remaining ratios are repeated-refactor dominated
(`2.08x` to `3.39x` CKTSO refactor ratios on the slowest refactor rows), so the
large gap is still the paper-level numeric executor/storage issue rather than
small BLAS calls or CSC input-width overhead.

A same-source CBLAS-enabled check kept the BLAS hypothesis bounded. With
`-DKLS_ENABLE_CBLAS_SUPERNODE=ON` and `KLS_ENABLE_CBLAS_SUPERNODE=1`,
`ASIC_680k` using the comparable METIS/KLU-first path measured `0.8070s`
initial factor, `0.0659s` repeated factor, and `0.0551s` refactor versus
`0.8188s`, `0.0648s`, and `0.0576s` for the normal build
(`build-cblas/kls_asic680k_metis_klufirst_cblas_t4_r1_ref1.json`,
`build/kls_asic680k_metis_klufirst_baseline_t4_r1_ref1.json`). Both rows
reported zero CBLAS supernode update, GEMV, and TRSV counters. The forced
KLS-first AMD `ASIC_680k` stress path also reported zero CBLAS counters with
the runtime gate both on and off, and the matching CBLAS-enabled `pre2` forced
KLS-first run still timed out at 120s with an empty JSON file
(`build-cblas/kls_pre2_klsfirst_cblas_factor_timeout120.json`). This makes the
proposed "BLAS only for large cases" rule already true for the tested paths;
the unresolved `pre2` and CKTSO-gap losses are not caused by unguarded small
BLAS calls.

The follow-up `pre2` stack inspection isolated the next concrete overhead
inside the row-first storage bridge. In the default BTF-parallel KLS-first
factor path, three workers completed the small BTF blocks and exited; the
surviving worker was serially factoring the 629k-row dominant block, with stack
samples in row-first supernode panel-cache append/growth. A no-BTF sample
reached the restartable intra-block row pipeline, but most workers were waiting
on the pipeline condition while the active worker rebuilt or grew the same
panel-cache storage. KLS now lets `realloc` preserve the row-first panel-cache
metadata, dense values, tail values, and tail-column buffers instead of
allocating fresh storage, zeroing double buffers, copying the old contents in
KLS, and freeing the old allocation at each growth step. Same-session smoke and
CTest passed; `ASIC_680k` forced KLS-first moved from `12.1468s` to `12.1071s`
initial factor and from `0.1160s` to `0.0998s` repeated factor, while the
normal METIS/KLU-first row moved from `0.8188s`/`0.0648s`/`0.0576s`
initial/factor/refactor to `0.7139s`/`0.0577s`/`0.0410s` in a one-pass check
(`build/kls_asic680k_realloc_panel_t4_r1_ref0.json`,
`build/kls_asic680k_metis_klufirst_realloc_panel_t4_r1_ref1.json`). The
bounded `pre2` forced KLS-first factor run still timed out at 120s with an
empty JSON file
(`build/kls_pre2_realloc_panel_klsfirst_factor_timeout120.json`), so this
narrows avoidable storage-growth overhead but does not replace the missing
dominant-block row/ETree executor.

A same-source CBLAS guard rerun confirmed that adding a "large cases only"
BLAS guard would not address the current slow case, because that guard is
already present and the tested paths do not enter CBLAS. The CBLAS-capable
`ASIC_680k` METIS/KLU-first control reported zero CBLAS update counters with
`KLS_ENABLE_CBLAS_SUPERNODE=0` and with `=1`
(`build-cblas/kls_asic680k_cblas_guard_off_current_t4_r1_ref1.json`,
`build-cblas/kls_asic680k_cblas_guard_on_current_t4_r1_ref1.json`). The
matching forced `pre2` run with `KLS_ENABLE_CBLAS_SUPERNODE=1` still timed out
at 120s and left an empty JSON file
(`build-cblas/kls_pre2_cblas_guard_current_t4_factor_timeout120.json`). The
useful code cleanup from this pass is exact rather than heuristic: the
row-first pivot-restart path no longer calls `kls_row_first_supernodes_reset`
immediately before `kls_row_first_pipeline_rebuild_prefix_panel_cache`, because
the rebuild helper resets the same prefix before rebuilding the cache.

A later CBLAS fallback audit did find one exact policy hole, but not the
CKTSO-gap cause. Checked compact dense panels used to skip the native blocked
kernel whenever `KLS_ENABLE_CBLAS_SUPERNODE=1`, even for panels too small for
any guarded CBLAS row update. KLS now only suppresses the checked native blocked
kernel when at least one row in the panel can pass the same structural and
minimum-work CBLAS gate. The normal and CBLAS smoke suites pass with a small
checked-panel fixture forcing the runtime CBLAS gate on and requiring native
blocked-panel stats. A current same-binary top-ten CBLAS rerun still reports
zero CBLAS and compact blocked-panel counters in both gate states
(`build-cblas/kls_cblas_on_checked_guard_gap10_t4_r1_ref3_timeout120.jsonl`,
`build-cblas/kls_cblas_off_checked_guard_gap10_t4_r1_ref3_timeout120.jsonl`),
with geomeans `4.0457s` on and `4.0592s` off. The remaining top-ten signal is
therefore still the column/EGraph repeated-refactor path and missing row-major
executor coverage, not small BLAS calls.

KLS now exposes row-first panel-cache staging counters so the next
SubtreeLU/CKTSO storage work can be separated from actual cached-panel
consumption. Public stats and `kls_bench` JSON report
`kls_first_row_panel_cache_build_count`,
`kls_first_row_panel_cache_build_panels`,
`kls_first_row_panel_cache_build_entries`,
`kls_first_row_panel_cache_append_count`,
`kls_first_row_panel_cache_append_panels`, and
`kls_first_row_panel_cache_append_entries`. Synthetic smoke coverage now
requires nonzero append counters when KLS-first panel updates are expected and
nonzero build counters when a separator pivot epoch rebuilds the prefix cache.
Focused checks showed that forced KLS-first `transient` exercised the cache
heavily (`3,814` appended panels and `1,075,649` stored entries with `97,678`
cached panel update groups in
`build/kls_transient_panel_cache_stats_klsfirst_t4_r1_ref0.json`), while the
forced KLS-first `ASIC_680k` stress control reported zero panel-cache activity
(`build/kls_asic680k_panel_cache_stats_klsfirst_t4_r1_ref0.json`). The forced
KLS-first `pre2` factor probe still timed out at 120s and left an empty JSON
file (`build/kls_pre2_panel_cache_stats_klsfirst_factor_timeout120.json`), so
the new counters do not close that gap directly; they make the next large-case
diagnosis explicit about whether time is going into full prefix rebuilds,
incremental published-panel staging, or later cached-panel consumption.

The June 28, 2026 top-ten CKTSO-gap rerun again rejected "too many small BLAS
calls" as the active explanation. The normal build is compiled without CBLAS,
reported `build_has_cblas=false`, and recorded zero CBLAS update runs, rows,
and entries while measuring a `4.1533s` geomean
(`build/kls_default_group_shape_current_gap10_t4_r1_ref3_timeout120.jsonl`).
The same source built with CBLAS reported `4.1118s` with
`KLS_ENABLE_CBLAS_SUPERNODE=0` and `4.0363s` with
`KLS_ENABLE_CBLAS_SUPERNODE=1`; both runs still recorded zero CBLAS update
runs, rows, entries, and blocked-panel updates
(`build-cblas/kls_cblas_guard_off_group_shape_current_gap10_t4_r1_ref3_timeout120.jsonl`,
`build-cblas/kls_cblas_guard_on_group_shape_current_gap10_t4_r1_ref3_timeout120.jsonl`).
The external BLAS paths are therefore already build-time optional, runtime
opt-in, and size/work gated for the current slow cases.

A same-commit CBLAS guard rerun after the BTF EGraph no-wait split kept that
conclusion unchanged. The accepted non-CBLAS top-ten run measured a `3.9789s`
geomean with `build_has_cblas=false` and zero CBLAS/blocked-panel counters
(`build/kls_btf_nowait_split_gap10_t4_r1_ref3_timeout120.jsonl`). Rebuilding
the CBLAS tree at commit `af1be09` and running the same suite with
`OPENBLAS_NUM_THREADS=1` measured `4.0990s` with
`KLS_ENABLE_CBLAS_SUPERNODE=0` and `4.1731s` with
`KLS_ENABLE_CBLAS_SUPERNODE=1`
(`build-cblas/kls_cblas_guard_off_current_gap10_t4_r1_ref3_timeout120.jsonl`,
`build-cblas/kls_cblas_guard_on_current_gap10_t4_r1_ref3_timeout120.jsonl`).
Both CBLAS-capable runs still recorded zero CBLAS update runs, rows, entries,
and zero native blocked-panel update runs on all ten matrices. This rejects the
"small BLAS calls" hypothesis for the active default CKTSO-gap loss; another
BLAS-size guard would not change the executed code path.

A June 28, 2026 same-binary top-five rerun of the current source reached the
same conclusion for the user's proposed large-case BLAS guard. With
`OPENBLAS_NUM_THREADS=1`, the CBLAS-capable binary measured `1.4959s` geomean
with `KLS_ENABLE_CBLAS_SUPERNODE=0` and `1.4259s` with
`KLS_ENABLE_CBLAS_SUPERNODE=1`
(`build-cblas/kls_cblas_guard_current_off_gap5_t4_r1_ref3_timeout120.jsonl`,
`build-cblas/kls_cblas_guard_current_on_gap5_t4_r1_ref3_timeout120.jsonl`).
Both runs reported `build_has_cblas=true`, but every matrix recorded zero
external CBLAS update runs, rows, and entries, and zero row-refactor compact
GEMV counters. The small variations between the two timings are therefore not
caused by dispatching small CBLAS calls. The active slow path is still missing
broader paper-style row/supernode executor coverage, not another BLAS size
threshold.

KLS now also reports row-refactor group shape in public stats, `kls_bench`
JSON/text output, and `scripts/decompose_solver_gap.py`. The forced row-engine
control with the restored `KLS_ROW_REFACTOR_BATCH_MAX_ROWS=16` measured a
`11.8996s` geomean and showed that singleton groups dominate: `1,107,031` of
`1,187,615` groups were single-row groups, with only `29,059` contiguous batch
groups covering `323,633` rows
(`build/kls_forced_row_group_shape_gap10_t4_r1_ref3_timeout120.jsonl`). The
rejected `64`- and `256`-row batch probes stayed in the same poor band
(`12.1113s` and `11.8470s` post-counter geomeans), so widening the current
contiguous-batch threshold is not a paper-level fix. The next credible gap is
still a durable row/supernode numeric executor that creates reusable
BLAS-sized work and coarser same-level/ETree tasks, not another small-BLAS
guard.

The direct follow-up to the BLAS-threshold concern again showed that small BLAS
calls are not active on the current forced-row slow set. A CBLAS-capable
top-five forced-row run with `OPENBLAS_NUM_THREADS=1` and
`KLS_ENABLE_CBLAS_SUPERNODE=1` measured `21.9360s` geomean, versus `21.3937s`
for the first five rows of the non-CBLAS group-shape control
(`build-cblas/kls_forced_row_cblas_on_focus5_t4_r1_ref3_timeout120.jsonl`,
`build/kls_forced_row_group_shape_gap10_t4_r1_ref3_timeout120.jsonl`). The
same CBLAS run reported zero refactor CBLAS update runs on all five rows. Its
blocked native compact-panel counters were still active, so adding another
"large cases only" BLAS guard would not change the exercised path.

The supported fix from the same pass was instead to keep the paper-style
EGraph refactor enabled for medium-heavy dominant-BTF cases with many fringe
blocks. The earlier `5000` block cap excluded `HTC_336_4438` even though its
largest block covered about 87% of the matrix and the refactor work was
dominated by that block. Raising the cap to `50000` changed that row from the
generic BTF pool path to `egraph`: the same-session old-cap control measured
`22.2447s` with `0.1728s` average refactors, while the candidate measured
`10.6333s` with `0.0559s` average refactors
(`build/kls_medium_manyblock_control_HTC_t4_r1_ref5_timeout120.jsonl`,
`build/kls_medium_manyblock_egraph_HTC_t4_r1_ref5_timeout120.jsonl`). On the
current top-20 CKTSO-gap manifest the geomean moved from the accepted
`3.7633s` artifact to `3.7127s`
(`build/kls_btf_nowait_split_gap20_t4_r1_ref2_timeout120.jsonl`,
`build/kls_medium_manyblock_egraph_gap20_t4_r1_ref2_timeout120.jsonl`). The
remaining `HTC_336_4438` gap then shifts from refactor scheduling to initial
factorization/row-up coverage, which is consistent with the paper-level
executor gap rather than with small BLAS calls.

Two scheduler probes also narrowed the paper gap. Disabling the separator-FLOP
queue measured `11.6238s` geomean against the `11.8996s` control, with mixed
per-matrix movement rather than a clear replacement
(`build/kls_forced_row_no_sepflop_probe_gap10_t4_r1_ref3_timeout120.jsonl`).
Forcing the row-dependency release wrapper whenever any separator pipeline
work exists regressed to `12.1698s`
(`build/kls_forced_row_sep_rowdep_probe_gap10_t4_r1_ref3_timeout120.jsonl`).
Those results reject another broad queue switch as the direct CKTSO-gap fix.

KLS now reports separator-FLOP private-queue balance directly. The diagnostic
top-ten forced-row run measured `11.3918s` geomean, which is in the same noisy
band as the prior forced-row controls, but its queue shape is more important
than the timing
(`build/kls_forced_row_sep_balance_gap10_t4_r1_ref3_timeout120.jsonl`). On
`G2_circuit`, the separator-FLOP queue used all four private threads with
17,238 to 17,540 groups and `2.950e8` to `3.176e8` estimated work per nonempty
private queue, so the remaining loss there points back to numeric layout and
compact update throughput. On `ASIC_320k`, `ASIC_320ks`, `ASIC_100ks`, and
`rajat28`, the same Algorithm 6-style queue used only one nonempty private
thread despite 242,653, 218,928, 54,214, and 33,912 private groups
respectively. That is a concrete paper-aligned gap: KLS retained the
private/pipeline structure, but its separator-private component generation or
queue partitioning is not exposing independent private subtrees on several
slow SPICE matrices.

Three follow-up probes rejected narrower explanations for the single-private
queue collapse. Cutting row-refactor groups at retained separator-component
boundaries reduced closure groups on the ASIC and `rajat28` rows, but still
left the same one nonempty private thread and regressed the top-ten forced-row
geomean to `11.6704s`
(`build/kls_forced_row_sep_component_groups_gap10_t4_r1_ref3_timeout120.jsonl`).
Forcing the separator-FLOP split loop to keep splitting until it had at least
one candidate per worker did not change the reported queue shapes at all and
measured `11.5917s`
(`build/kls_forced_row_sep_min_candidates_gap10_t4_r1_ref3_timeout120.jsonl`).
Finally, dropping the static private thread slices for collapsed separator
queues and letting workers dynamically pull private groups through the exact
predecessor counters regressed to `11.9151s`
(`build/kls_forced_row_sep_dynamic_private_gap10_t4_r1_ref3_timeout120.jsonl`).
Together these results point below the current row-group scheduler policy: the
retained separator tree or numeric work model is not exposing multiple useful
private subtrees for those matrices, and a generic dynamic/private fallback is
not enough to recover CKTSO-level performance.

The production EGraph dependency loop now also skips producer-column metadata
loads when the dependency value is exactly zero. The earlier zero-scatter guard
avoided the arithmetic row walk, but the hot loops still fetched the producer
`L` pointer, values pointer, and length before entering the guarded scatter
helper. The current cleanup keeps the stored `U` value and workspace clearing
unchanged while bypassing those metadata loads in the serial mapped, parallel
mapped, and threaded EGraph refactor kernels. A top-ten CKTSO-gap repeat-3 run
measured `4.0034s` and `4.0203s` geomeans across two passes, versus the saved
current default `4.1533s`
(`build/kls_zero_dep_metadata_skip_gap10_t4_r1_ref3_timeout120.jsonl`,
`build/kls_zero_dep_metadata_skip_gap10_t4_r1_ref3_timeout120_run2.jsonl`,
`build/kls_default_group_shape_current_gap10_t4_r1_ref3_timeout120.jsonl`).
The common-row KLS/CKTSO geomean ratio moved from about `2.43x` to `2.34x`.
A broader top-20 repeat-2 slice completed without failures at `3.8061s`
geomean (`build/kls_zero_dep_metadata_skip_gap20_t4_r1_ref2_timeout120.jsonl`).
This is a useful production-kernel cleanup for the currently dominant sparse
scatter path, but it is still a small effect relative to CKTSO; it does not
replace the missing paper-level row/supernode numeric executor.

A more aggressive dynamic-dependency variant was rejected. The probe moved the
pipeline wait until after reading the current work-vector dependency value and
skipped the wait when an unchecked no-pivot dependency was exactly zero, on the
theory that a zero multiplier does not need the producer `L` column. It stayed
correct on smoke tests but regressed the top-ten CKTSO-gap repeat-3 geomean to
`4.1126s` and `4.1145s` across two passes
(`build/kls_zero_dep_wait_skip_gap10_t4_r1_ref3_timeout120.jsonl`,
`build/kls_zero_dep_wait_skip_gap10_t4_r1_ref3_timeout120_run2.jsonl`).
That result keeps the structural pipeline wait order intact; the retained
zero-dependency improvement is only the producer metadata/scatter bypass after
the dependency is ready.

A zero-`L` store division probe was also rejected. It changed the EGraph
`L`-column store helper to write exact zero directly for nonzero pivots instead
of dividing `0.0 / pivot`, while preserving the old zero-pivot behavior. The
extra branch regressed the top-ten CKTSO-gap repeat-3 geomean to `4.1316s`
(`build/kls_zero_l_store_div_skip_gap10_t4_r1_ref3_timeout120.jsonl`), so the
branch was removed. The profitable exact-zero handling remains at the producer
dependency level, where it avoids whole `L`-column metadata/scatter work.

A current-source policy probe also rejected disabling static pivoting to help
the small initial-factor-dominated `gemat12` row. `--no-static-pivoting`
improved `gemat12` from roughly `0.059s` to `0.049s`, but it broke the general
SPICE set: the top-ten geomean regressed to `9.3922s`, with `onetone2`
at `31.3s`, `onetone1` at `95.2s`, `rajat28` at `20.5s`, and `rajat24` at
`32.4s`
(`build/kls_no_static_pivot_gap10_t4_r1_ref3_timeout120.jsonl`). Static
pivoting therefore remains a necessary broad policy despite the small-case
overhead; closing the `gemat12` gap needs a cheaper small-factor path, not a
global static-pivot disable.

The current source was also rechecked against the two larger paper-aligned
opt-ins after the zero-dependency metadata improvement. EGraph cached
supernode updates with `KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1` still regressed
the top-ten CKTSO-gap set to `4.6533s` geomean versus about `4.00s-4.02s` for
the default (`build/kls_egraph_supernode_current_gap10_t4_r1_ref3_timeout120.jsonl`).
The probe built and consumed panels on every EGraph row, but all such rows
regressed; `gemat12` was the only win and it reported zero supernode panels, so
that movement is not evidence for enabling the panel path. Forced
row-solve-from-numeric also stayed off-policy: `--row-solve on` measured
`4.2294s` geomean and `G2_circuit` rose to `31.5s`
(`build/kls_row_solve_on_current_gap10_t4_r1_ref3_timeout120.jsonl`). These
reruns keep both features as correctness/coverage scaffolding for the future
row/supernode engine rather than current default performance levers.

The EGraph supernode-panel scaffold now keeps its structural panel cache across
repeated fixed-pattern refactors. The old builder freed and rebuilt the panel
map, trailing-row descriptors, dense storage, and trailing storage every time
`KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1` was used, even though schedule and
LU-pointer invalidation already free the cache when the structural state
changes. The new reuse check validates the existing schedule/pointer boundary
and refreshes the enabled-panel count before returning; published producer rows
still clear their own dense/trailing rows before refilling numeric values. On
the same top-ten CKTSO-gap opt-in benchmark, the geomean improved from
`4.6533s` to `4.5505s`
(`build/kls_egraph_supernode_current_gap10_t4_r1_ref3_timeout120.jsonl`,
`build/kls_egraph_supernode_cache_reuse_gap10_t4_r1_ref3_timeout120.jsonl`).
Panel and update counts stayed identical (`6913` panels and `41603` update
runs), which confirms this is a staging-overhead reduction rather than a new
numeric algorithm. A default control stayed residual-clean with zero supernode
panel/update counters at `4.1230s`
(`build/kls_cache_reuse_default_gap10_t4_r1_ref3_timeout120.jsonl`), so the
policy remains unchanged: cached panels are closer to the paper storage model
but still slower than the scalar EGraph default until the row/supernode executor
creates larger reusable work.

The default BTF EGraph refactor now separates the CKTSO-style cluster loop from
the waited pipeline loop inside the hot unscaled column kernel. Cluster columns
are launched only after their predecessor levels have completed, so the inner
dependency walk no longer tests `wait_for_dependencies` on every `U` entry; the
pipeline path still performs the same per-dependency wait before consuming a
producer column. This does not change ordering, pivot checks, or the opt-in
supernode update path. The top-ten CKTSO-gap repeat-3 guards moved from the
saved zero-dependency cleanup runs at `4.0034s` and `4.0203s` geomean to
`3.9789s` and `3.9548s`
(`build/kls_btf_nowait_split_gap10_t4_r1_ref3_timeout120.jsonl`,
`build/kls_btf_nowait_split_gap10_t4_r1_ref3_timeout120_run2.jsonl`). A
broader top-20 repeat-2 guard completed all rows at `3.7633s`, versus the
saved `3.8061s`
(`build/kls_btf_nowait_split_gap20_t4_r1_ref2_timeout120.jsonl`,
`build/kls_zero_dep_metadata_skip_gap20_t4_r1_ref2_timeout120.jsonl`). This is
another small production-kernel cleanup in the dominant scalar EGraph path; the
paper-level gap remains the durable row/supernode numeric executor that can
replace these scalar dependency walks with coarser reusable work.

The analogous single-block EGraph split was tested and rejected. Splitting both
single-block kernels into waited/no-wait loops produced correct tests but did
not generalize: two top-ten repeat-3 runs measured `3.9864s` and `3.9661s`
geomean (`build/kls_single_nowait_split_gap10_t4_r1_ref3_timeout120.jsonl`,
`build/kls_single_nowait_split_gap10_t4_r1_ref3_timeout120_run2.jsonl`) versus
the accepted BTF-only controls at `3.9789s` and `3.9548s`, and the broader
top-20 repeat-2 guard regressed from `3.7633s` to `3.8135s`
(`build/kls_single_nowait_split_gap20_t4_r1_ref2_timeout120.jsonl`). A
scaled-only isolation improved `Raj1` in a focused run but still regressed the
top-ten guard to `4.0511s`
(`build/kls_scaled_single_nowait_split_Raj1_t4_r1_ref5_timeout120.jsonl`,
`build/kls_scaled_single_nowait_split_gap10_t4_r1_ref3_timeout120.jsonl`).
This keeps the accepted loop split limited to the BTF kernel where the measured
gain was stable; the single-block cases need the larger row/supernode executor
rather than more branch-shape variants.

A higher-repeat current-source rerun tightened the CKTSO-gap ranking after the
medium-heavy BTF EGraph gate. The earlier top-20 guard used only two refactors
per matrix and over-weighted timing noise on short cases; rerunning the same
manifest with `--refactor-repeat 20` completed all rows at `3.3573s` geomean
(`build/kls_current_gap20_t4_r1_ref20_timeout120.jsonl`), versus `3.7127s` for
the repeat-2 artifact
(`build/kls_medium_manyblock_egraph_gap20_t4_r1_ref2_timeout120.jsonl`). This
does not prove parity with CKTSO: the stable decomposition still shows
`ASIC_320k` and `ASIC_320ks` at `2.97x` CKTSO cycle time, with no-pivot EGraph
refactors at about `3.30x` and `3.39x` CKTSO refactor time. The broader top-20
gap is now cleaner: `Raj1` was the only repeat-20 regression versus the previous
artifact, while `mc2depi` and `G2_circuit` were effectively unchanged. The
stable evidence keeps the remaining large lever focused on the
CKTSO/SubtreeLU row/supernode numeric executor rather than small BLAS-call
guards, KLS-first broad auto policy, or noisy short-repeat samples.

The same repeat-20 baseline also keeps cached EGraph supernode updates
off-policy. The best cached-panel opt-in still regressed every common top-ten
EGraph row against the stable default: `ASIC_320k` rose from `13.43s` to
`14.45s`, `ASIC_320ks` from `11.30s` to `12.50s`, `G2_circuit` from `20.01s`
to `23.98s`, and `rajat28` from `3.41s` to `4.88s`
(`build/kls_egraph_supernode_cache_reuse_gap10_t4_r1_ref3_timeout120.jsonl`).
The opt-in path executed thousands of panel updates on those rows, so this is
not an inactive feature or BLAS-threshold artifact. It is the wrong granularity
for the current column EGraph executor until the row-major/supernodal storage
path can amortize panel construction and consumption across coarser tasks.

A small-work guard was added to the opt-in cached EGraph panel update path to
test the remaining "small BLAS call" hypothesis without enabling any system
BLAS dependency. CBLAS calls were already build-time optional, runtime gated,
and size/work gated; the normal benchmark build still reports
`build_has_cblas=false`. The new native cached-panel gate skips retained-panel
updates before workspace allocation unless the dependency run has at least 16
contiguous rows, at least 512 estimated triangular/update operations, and at
least 8 estimated operations per staged panel/trailing entry. This removed the
obvious tiny-update pathology: on `ASIC_320k`, the opt-in path dropped from
`502` blocked cached updates over `1006` rows to zero blocked cached updates
(`2` scalar supernode updates remained). It improved the cached-panel opt-in
top-ten geomean only from `4.5505s` to `4.5126s`, still `1.1717x` slower than
the stable default top-ten common rows at `3.8514s`
(`build/kls_egraph_supernode_cache_reuse_gap10_t4_r1_ref3_timeout120.jsonl`,
`build/kls_egraph_supernode_workgate_gap10_t4_r1_ref3_timeout120.jsonl`,
`build/kls_current_gap20_t4_r1_ref20_timeout120.jsonl`). The guard is a useful
scaffold cleanup and stays behind `KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1`, but
it confirms that small BLAS/panel-call overhead is not the large missing piece;
large retained-panel updates on `ASIC_320ks` and `G2_circuit` still regress
because the current executor has not yet moved to the paper-level
row/supernode numeric algorithm.

Two direct scheduler-only probes were also rejected on June 28, 2026. First,
the forced row-refactor path was given shared/private-queue and level-sliced
private-work variants to test whether the observed private-thread collapse was
mostly queue policy. On `ASIC_100ks`, the control forced-row refactor was
`0.284688s`; the shared/private-queue gate slowed it to `0.294783s`. A
separate level-sliced private pass moved a `0.251362s` control to `0.284191s`.
Both variants were residual-clean but slower. Second, a level-order EGraph
supernode-task lease grouped contiguous `refactor_level_cols` under the
existing opt-in task gate. On default EGraph `ASIC_100ks`, it formed 20 tasks
over 140 columns across the run and stayed residual-clean, but repeated
refactor time moved from `0.0432357s` to `0.0480641s`. These failures reinforce
the paper diagnosis: queue reshaping without the production row/supernode
numeric storage and executor does not close the CKTSO gap.

A direct EGraph panel-publication storage probe was also rejected. The probe
filled the opt-in retained supernode panel row inside
`kls_egraph_store_l_column_from_workspace`, eliminating the separate
post-store scan of the just-written L column. That looked closer to the
paper's durable numeric-panel storage, but on the top-five CKTSO-gap EGraph
rows it made the opt-in path materially worse: the same-binary default measured
`8.4375s` geomean, while `KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1` measured
`9.7611s`
(`build/kls_inline_panel_default_gap5_t4_r1_ref3_timeout120.jsonl`,
`build/kls_inline_panel_optin_gap5_t4_r1_ref3_timeout120.jsonl`). The older
work-gated opt-in artifact was `8.4655s` on the same five rows, so the inline
publication was not just still off-policy; it also regressed the opt-in
scaffold. The likely cause is that panel fill pushed extra conditionals into
the hot L-store loop for every published panel column. The probe was removed:
the useful paper-level direction is still native row/supernode storage that
avoids KLU-column staging altogether, not fusing another consumer into the
current column-store loop.

The BTF EGraph hot loop was also given a branch-light default probe that split
the dependency walk into separate `supernode_numeric_updates` on/off loops. The
intent was to remove a disabled opt-in branch from every U-entry in the normal
scalar EGraph path without changing ordering or the opt-in supernode update
semantics. It was correct but not a stable CKTSO-gap improvement: repeat-20
top-five improved only because `onetone2` moved down, while `ASIC_320k`,
`ASIC_100ks`, and `G2_circuit` were slightly worse; the broader repeat-20
top-ten measured `3.8603s` geomean versus the saved `3.8514s` top-ten baseline
(`build/kls_btf_supernode_branch_split_gap5_t4_r1_ref20_timeout120.jsonl`,
`build/kls_btf_supernode_branch_split_gap10_t4_r1_ref20_timeout120.jsonl`,
`build/kls_current_gap20_t4_r1_ref20_timeout120.jsonl`). The split was removed.
This keeps the diagnosis unchanged: the large gap is not a disabled-branch
artifact in the scalar EGraph loop; KLS still needs the paper-level
row/supernode numeric executor.

Disabling row-refactor auto preparation outright was also checked after the
same decomposition showed many default EGraph rows paying a rejected lower-bound
model attempt. `KLS_ENABLE_ROW_REFACTOR=0` correctly suppressed the model and
lower-bound counters, but it did not give a stable policy win: the repeat-3
top-ten moved only from `4.0924s` to `4.0702s` geomean with mixed row-level
refactor noise, and the repeat-20 top-five moved from the saved `7.3959s` to
`7.3566s` only because `onetone2` improved while `ASIC_320k`, `ASIC_320ks`,
`ASIC_100ks`, and `G2_circuit` were flat to worse
(`build/kls_row_auto_disabled_gap10_t4_r1_ref3_timeout120.jsonl`,
`build/kls_row_auto_disabled_gap5_t4_r1_ref20_timeout120.jsonl`,
`build/kls_current_gap20_t4_r1_ref20_timeout120.jsonl`). The default keeps the
existing guarded auto-preparation semantics; the missing lever remains making
the prepared row/supernode path actually faster, not hiding its preparation.

A direct hot-metadata alias probe was also rejected. The probe cached
`numeric->Llen` and `numeric->Ulen` aliases inside the refactor LU pointer cache
and rewired the specialized EGraph refactor kernels to consume those aliases in
the repeated-refactor dependency loops. This was a narrow attempt to reduce
KLU-owned metadata touches without changing numeric semantics, but the rerun was
unambiguously worse: the full CKTSO-gap manifest completed 39 of 42 matrices,
with `HTC_336_9129`, `mac_econ_fwd500`, and `ss1` timing out at 120 seconds
(`build/kls_refactor_len_alias_gap10_t4_r1_ref20_timeout120.jsonl`,
`build/kls_refactor_len_alias_gap10_t4_r1_ref20_timeout120.failures`). On the
20 common rows against the stable repeat-20 baseline, the probe had no wins and
regressed geomean cycle time by `1.412x`; `ASIC_320k`, `G2_circuit`, and
`onetone1` moved from `13.43s`, `20.01s`, and `4.97s` to `19.65s`, `27.96s`,
and `9.21s`. The code was removed. Together with the CBLAS-on/off artifacts
that recorded zero CBLAS update calls on these hard rows, this rejects the
"small BLAS call" and hot-length-metadata hypotheses as primary causes. The
remaining gap still points at the paper-level row/supernode numeric executor
and storage model.

A dependency-aware small-batch grouping probe was also rejected. The probe kept
the old eight-row minimum for root/no-dependency batches but allowed adjacent
same-level dependent rows to form two-row-or-larger batch groups, trying to give
the existing row/supernode executor more consumer rows without touching BLAS
thresholds. The smoke suite passed, but the forced row-refactor top-ten CKTSO-gap
run regressed from `11.8996s` to `14.5117s` geomean against the saved
group-shape baseline
(`build/kls_forced_row_group_shape_gap10_t4_r1_ref3_timeout120.jsonl`,
`build/kls_dep_small_batch_gap10_t4_r1_ref3_timeout120.jsonl`). The probe did
increase batch exposure (`row_refactor_group_batch_count` from `29059` to
`72324`, and `row_refactor_group_batch_rows` from `323633` to `435106` across
the ten rows), but it still reported zero
`row_refactor_compact_supernode_batch_candidate_count` and zero compact
supernode batches. It also recorded zero CBLAS update runs with
`build_has_cblas=false`. The code was removed. This rejects "just make smaller
same-level batches" along with the small-BLAS-call hypothesis; KLS needs
producer-compatible row-panel/supernode batch construction, not more generic
batch groups.

The independent-row multi-producer scaffold then gained retained dense-producer
target maps, but still remains opt-in. Previously those maps only described
dense current groups, so the independent fragmented/exact row batches still did
per-update searches from dense-producer trailing columns into each consumer
row's later L, pivot, or row-major U workspace. KLS now builds the independent
target slices only when `KLS_ENABLE_MULTI_PRODUCER_SUPERNODE=1` is set and the
independent multi-producer executors carry retained run IDs to consume those
slices. A first implementation moved the opt-in top-ten geomean from
`25.4948s` to `22.9011s`; after moving target-slice validation out of the
per-update loop, the same forced row-refactor top-ten opt-in rerun measured
`20.7554s` geomean, with 10 wins over the prior multi-producer artifact
(`build/kls_multi_producer_current_gap10_t4_r1_ref3_timeout120.jsonl`,
`build/kls_multi_targetmap_validated_gap10_t4_r1_ref3_timeout120.jsonl`).
Batch coverage was effectively unchanged (`6864` compact supernode batches and
`2717328408` batch entries), while the independent target maps eliminated the
retained `NONE` target slots in the opt-in run. With the gate off, the top-ten
target-map counters match the saved dense-only baseline exactly, so this is not
a default policy change. The opt-in path is still much slower than the default
forced row-refactor baseline (`20.7554s` versus `11.8996s` geomean, `1.744x`
candidate/reference), so the remaining paper gap is still the lower-overhead
row-panel storage/update kernel, not just retained target lookup or small BLAS
thresholding.

A follow-up worker-scratch probe was rejected. It moved the independent
multi-producer batch metadata (`run_groups`, per-row run offsets/suffixes,
run lengths, run ids, L offsets, and U offsets) into the reusable worker index
workspace and combined pivots, multipliers, and U workspaces into one reusable
worker supernode workspace. This removed many per-batch heap allocations and
looked closer to a persistent row-panel executor, but it regressed the same
opt-in forced row-refactor top-ten run from `20.7554s` to `22.2698s` geomean
with the same `6864` compact supernode batches and `2717328408` batch entries
(`build/kls_multi_targetmap_validated_gap10_t4_r1_ref3_timeout120.jsonl`,
`build/kls_multi_worker_scratch_gap10_t4_r1_ref3_timeout120.jsonl`). The code
was removed. The result narrows the gap further: ordinary heap allocation in
these independent batches is not the main missing piece; the larger issue is
still the arithmetic/data-layout granularity of the row-panel update itself.

A row-major mapped-update probe was rejected for the same reason. In the exact
independent multi-producer branch, retained target maps remove the need to
search each dense-producer trailing column, so the probe inverted the mapped
trailing update loop to process all target columns for one consumer row before
moving to the next row. This looked more cache-local for the consumer
multiplier/pivot/U workspace, but the same opt-in top-ten run regressed from
`20.7554s` to `23.4976s` geomean with identical compact-supernode batch
coverage (`6864` batches and `2717328408` entries) and zero CBLAS calls
(`build/kls_multi_rowmajor_target_gap10_t4_r1_ref3_timeout120.jsonl`). The
code was removed. Together with the worker-scratch probe, this points away from
simple staging/locality rearrangements inside the current per-consumer row
kernel and toward a more substantial packed row-panel update object.

A narrower packed exact-update probe was also rejected. It kept the original
column-major scatter order for retained target maps, but materialized each
exact independent multi-producer run's `batch_rows x trailing_len` update
panel in worker supernode workspace before scattering, mirroring the dense
fragmented batch kernel's packed update object. The top-ten opt-in run
regressed from `20.7554s` to `23.3306s` geomean with the same `6864`
compact-supernode batches, the same `2717328408` batch entries, and zero CBLAS
calls
(`build/kls_multi_packed_exact_update_gap10_t4_r1_ref3_timeout120.jsonl`).
The code was removed. This makes the remaining direction more specific:
packing updates inside the existing row-batch executor is not enough; the data
layout likely needs to avoid the current per-consumer row workspace/scatter
contract rather than adding another temporary panel to it.

KLS now removes another direct mismatch with CKTSO/SubtreeLU's
prefactor/postfactor row algorithm. The row pipeline used to stop trying
retained compact-supernode producer updates after a row had prefactored any
finished dependency while waiting for an earlier dependency; the remaining
post-wait dependencies then fell back to scalar row updates even when a valid
producer panel existed. The compact-supernode update helper is now
applied-mask aware: it accepts the row's already-consumed dependency bitmap,
rejects only candidate runs that would overlap prefactored entries, and remains
available for later postfactor dependencies. This matches Algorithm 5's intent
that both prefactorization and postfactorization can consume supernodes while
preserving the existing scalar fallback. The normal build and smoke suite pass,
but the current hard forced-row top-five repeat-20 guard did not exercise this
case: all rows in
`build/kls_postprefactor_mask_default_gap5_t4_r1_ref20_timeout120.jsonl`
reported zero `row_refactor_prefactor_*` counters. This is therefore a
semantic gap closure, not evidence of a CKTSO-gap speedup on the present focus
set.

A related policy probe was rejected in the same pass. Making
`KLS_ENABLE_PARTIAL_SUPERNODE_PIPELINE` opt-in by default looked promising on a
noisy top-eight repeat-3 forced-row run, but the higher-repeat top-five guard
rejected it: default-off measured `19.1953s` geomean in
`build/kls_partial_pipeline_default_off_gap5_t4_r1_ref20_timeout120.jsonl`,
while forcing the paper partial-prefix release back on measured `16.9455s` in
`build/kls_partial_pipeline_forced_on_gap5_t4_r1_ref20_timeout120.jsonl`. The
default remains unchanged; partial-prefix release stays on unless explicitly
disabled with `KLS_ENABLE_PARTIAL_SUPERNODE_PIPELINE=0`.

A current same-source BLAS-threshold rerun keeps the "small BLAS call" concern
closed for the default CKTSO-gap path. The CBLAS-capable Release tree was
rebuilt at the current source, then the top-five gap manifest was run with
`OPENBLAS_NUM_THREADS=1` and the runtime gate off/on. The first pass measured
`7.71968s` geomean with `KLS_ENABLE_CBLAS_SUPERNODE=0` in
`build-cblas/kls_cblas_guard_off_current_blashyp_gap5_t4_r1_ref3_timeout120.jsonl`
and `8.89658s` with `KLS_ENABLE_CBLAS_SUPERNODE=1` in
`build-cblas/kls_cblas_guard_on_current_blashyp_gap5_t4_r1_ref3_timeout120.jsonl`.
A second noisy pass measured `10.4417s` on and `9.86407s` off in
`build-cblas/kls_cblas_guard_on_current_blashyp2_gap5_t4_r1_ref3_timeout120.jsonl`
and
`build-cblas/kls_cblas_guard_off_current_blashyp2_gap5_t4_r1_ref3_timeout120.jsonl`.
All four runs used `last_refactor_path=egraph`, reported zero
`row_refactor_run_count`, zero native blocked-panel update runs, and zero
`refactor_supernode_cblas_update_*` runs/rows/entries on every matrix. The
existing BLAS policy is already build-time optional, runtime opt-in, and
guarded by 512-row/vector or 512-width panel minima plus multi-million-work
thresholds. Another "BLAS only for large cases" guard would not change the
executed code path; the remaining gap is still the production row/supernode
numeric storage and executor that would create reusable BLAS-sized work.

A schedule-side row-refactor lower-bound cache was also rejected. The idea was
to reuse the U predecessor counts already computed by the EGraph schedule
builder so the automatic row-refactor model would not rescan the same
dependency graph merely to reject the row path by lower-bound work. This is a
reasonable NICSLU-style policy cleanup, but it is not the CKTSO-scale missing
numeric executor, and the measured result did not justify keeping it. The
top-ten focus run regressed to `5.69422s` geomean in
`build/kls_cached_row_lower_bound_gap10_t4_r1_ref3_timeout120.jsonl` versus
the saved default `4.18460s`
(`build/kls_current_default_gap10_t4_r1_ref3_timeout120.jsonl`); a top-five
confirmation also regressed to `11.2459s` in
`build/kls_cached_row_lower_bound_gap5_rerun_t4_r1_ref3_timeout120.jsonl`.
The lower-bound values and rejection decisions matched the saved default
artifacts, so the change only moved bookkeeping around and did not alter the
default EGraph path selection. The source was reverted. This keeps the
diagnosis focused on production row/supernode storage and arithmetic, not on
auto-model bookkeeping.

The natural all-pipeline single-block EGraph path was also tested with a true
successor ready queue, because `G2_circuit` is the one current top-gap row
where the NICSLU task-flow model reports a strong parallel recommendation
(`parallel_task_flow_speedup=3.94669`) while the previous
`KLS_ENABLE_EGRAPH_READY_QUEUE=1` probe reported zero ready-queue columns.
The probe allowed the existing EGraph ready-queue machinery to use natural
global column IDs instead of requiring `refactor_level_cols`. It was residual
clean and exercised the intended graph: with the gate on, `G2_circuit`
reported `refactor_last_ready_queue_columns=150102` and
`refactor_ready_queue_run_count=21` in
`build/kls_natural_ready_queue_g2_on_t4_r1_ref20_timeout120.jsonl`. It was much
slower than the same-session natural-order control
`build/kls_natural_ready_queue_g2_off_t4_r1_ref20_timeout120.jsonl`: SPICE
cycle `52.6120s` versus `29.6737s`, average refactor `0.488389s` versus
`0.261114s`, and average fast-factor refactor `0.426144s` versus `0.235491s`.
The code was removed. This closes the specific "ready queue was never tried on
G2" hole, and again points away from column-scheduler policy toward a different
numeric representation/update granularity for the single-block CKTSO gap.

A more literal NICSLU Algorithm 4-style static task-flow assignment was then
prototyped behind `KLS_ENABLE_EGRAPH_TASK_FLOW=1` for the same huge
single-block EGraph case. The prototype used the existing task-flow model to
assign each column to the thread with the earliest predicted finish time, kept
each thread's assigned columns in natural dependency order, and reused the
existing pipeline-done waits for cross-thread dependencies. This avoided the
global successor ready queue, but it still ran the current EGraph per-column
numeric kernel. On `G2_circuit`, the same-source control measured `26.7858s`
in `build/kls_task_flow_g2_off_t4_r1_ref20_timeout120.jsonl`, while the
static task-flow executor measured `38.2774s` in
`build/kls_task_flow_g2_on_t4_r1_ref20_timeout120.jsonl`. The average refactor
time regressed from `0.232429s` to `0.341268s`, with the model still reporting
`parallel_task_flow_speedup=3.94669` and residuals unchanged. Both runs used
`build_has_cblas=false` and recorded zero CBLAS and supernode-update counters,
so this result does not support the "too many small BLAS calls" hypothesis.
The source was removed. The remaining gap is the paper-side executor/data
structure combination that makes scheduled tasks coarse and cache/BLAS-friendly,
not another thin scheduling wrapper around the current scalar EGraph column
kernel.

A single-row row-refactor allocator cleanup was tried and rejected after the
same BLAS-threshold discussion. The grouped row executors already reuse
`worker->byte_workspace` for their dependency-applied bitmaps, so the
single-row executor was changed to use the same per-worker workspace instead of
per-row `calloc/free`. This is source-structure clean, but the same-source
forced-row native-off benchmark regressed on every top-five gap matrix:
`23.5832s` geomean for the reverted control in
`build/kls_control_forced_row_native_off_gap5_t4_r1_ref3_timeout120.jsonl`
versus `29.3486s` for the workspace patch in
`build/kls_applied_workspace_forced_row_native_off_gap5_t4_r1_ref3_timeout120.jsonl`.
Both runs used `--row-refactor refactor`, `KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=0`,
four threads, and three refactor repeats, and both were residual-clean. The
source was reverted. This closes that allocator-churn hypothesis; it does not
change the earlier conclusion that the default CKTSO-gap path is not entering
CBLAS at all and needs a coarser paper-style numeric executor rather than
another small-case BLAS guard.

A fresh paper-medium rerun on June 28, 2026 keeps that conclusion broad, not
top-five-specific. The current default build completed 89 of the 93 public
medium matrices in
`build/kls_current_paper_medium_t4_r1_ref3_timeout120.jsonl`, with timeouts on
`mac_econ_fwd500`, `ss1`, and `HTC_336_9129`, plus the existing singular
`bips07_1998` setup failure. Against the saved CKTSO artifact, common
successful rows measured `0.258758s` KLS geomean versus `0.253981s` CKTSO
geomean, a `1.0188x` KLS/CKTSO ratio, with 38 wins and 47 losses over 2%.
The top losses remain refactor-heavy: `ASIC_320ks` and `ASIC_320k` are about
`3.0x` slower, followed by `gemat12`, `rajat03`, `ASIC_100ks`, `onetone2`,
`rajat25`, `rajat28`, `rajat20`, and `onetone1`. The broad artifact reports
`build_has_cblas=false`, zero CBLAS update runs, zero EGraph supernode update
runs, zero blocked-panel update runs, and zero row-refactor runs on every
successful row, while still discovering about 977k refactor supernode-candidate
rows with 41.6M trailing entries. This rules out unguarded small BLAS calls as
the broad CKTSO-gap cause; the missing paper-level piece is still the
production row/supernode numeric representation and executor that can turn
those candidates into coarse reusable update work.

The EGraph supernode probe was split to isolate that durable-storage question.
`KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=cached` now builds and publishes the
persistent retained panels but skips the older per-consumer temporary panel
reconstruction fallback used by `=1`. Panel pruning also applies the same
cached-update work gate as the runtime consumer, so opt-in runs do not retain
structurally reachable panels that cannot pass the cached update threshold.
This confirmed the fallback was extra cost but did not make the cached consumer
production-ready: on the refreshed top-ten focus set, the default control
measured `2.23151s` geomean in
`build/kls_cached_workprune_control_gap10_t4_r1_ref3_timeout120.jsonl`,
cached-only measured `2.38338s` in
`build/kls_cached_workprune_gap10_t4_r1_ref3_timeout120.jsonl`, and full
supernode updates measured `2.52051s` in
`build/kls_supernode_full_workprune_gap10_t4_r1_ref3_timeout120.jsonl`. A
higher-repeat `onetone2`/`rajat25` check also rejected the apparent noisy
`onetone2` win: cached-only measured `1.60269s`/`3.50812s` versus the default
`1.42445s`/`3.35944s`. The retained `cached` mode is therefore an experimental
diagnostic for the SubtreeLU-style persistent-panel path, not a default speed
policy. The remaining missing algorithm is coarser producer/consumer batching
or row-major supernode execution that amortizes panel publication across larger
work units.

A direct default EGraph L-index cache hoist was tested and rejected. The patch
passed the already-built `refactor_l_indices32` pointer into the hot scatter
helper so each dependency update would not reload it through the solver object.
That is a harmless storage-access cleanup, but it is not the paper-level
executor change and did not improve the hard refactor rows. With the patch, the
two ASIC top-gap rows measured `12.5128s` geomean in
`build/kls_lindex_cached_scatter_gap2_t4_r1_ref20_timeout120.jsonl`; the
same-session reverted control measured `12.5034s` in
`build/kls_lindex_cached_scatter_control_gap2_t4_r1_ref20_timeout120.jsonl`.
The source stayed reverted. This keeps attention on coarser row/supernode
execution rather than per-dependency pointer plumbing.

The BLAS small-case guard hypothesis was rechecked with a CBLAS-capable build
on the same top-ten CKTSO-gap focus set. `build-cblas` was configured with
`-DKLS_ENABLE_CBLAS_SUPERNODE=ON`, then run once with
`KLS_ENABLE_CBLAS_SUPERNODE=0` and once with `=1`. The disabled control measured
`2.50754s` geomean in
`build-cblas/kls_cblas_off_gap10_t4_r1_ref3_timeout120.jsonl`; the enabled run
measured `2.57907s` in
`build-cblas/kls_cblas_on_gap10_t4_r1_ref3_timeout120.jsonl`. Both artifacts
reported `build_has_cblas=true` on all ten rows, but the enabled run still
reported zero `refactor_supernode_cblas_update_run_count`, zero compact
supernode GEMV/TRSV counts, and zero blocked dense-panel CBLAS runs on every
row. The current BLAS path is already compile-time optional, runtime gated, and
large-work gated, so adding another small-case guard cannot close the observed
CKTSO gap. The slow rows are not paying small BLAS-call overhead; they are not
using BLAS at all.

The cached EGraph supernode experiment was then instrumented and tightened at
the producer-probe boundary. Before the precheck, the cached top-five run in
`build/kls_cached_probe_diag_gap5_t4_r1_ref3_timeout120.jsonl` measured
`1.63104s` geomean and showed millions of per-dependency cached probes for
only thousands of retained-panel hits: for example, `ASIC_320ks` had
6,347,752 probes, 18,548 panel hits, 2,596 contiguous runs, and 1,196 accepted
updates. The source now skips the heavier cached-panel probe unless the
producer column is actually mapped to a retained panel, and benchmark JSON
reports `refactor_supernode_cached_probe_*` counters for attempts, panel hits,
contiguous runs, work-gate acceptance, and applied updates. The same top-five
diagnostic improved to `1.48966s` geomean in
`build/kls_cached_probe_precheck_gap5_t4_r1_ref3_timeout120.jsonl`, and the
top-ten cached run measured `2.30940s` in
`build/kls_cached_probe_precheck_gap10_t4_r1_ref3_timeout120.jsonl` versus the
older cached-only `2.38338s`. This is useful executor cleanup for the
SubtreeLU-style path, but it still trails the default top-ten control
(`2.23151s`), so cached EGraph supernode updates remain opt-in. The counters
make the next gap clearer: panel discovery is not enough; KLS needs a coarser
producer/consumer executor that amortizes the remaining panel publication and
application work across larger row batches.

A same-session scale and row-refactor control sweep also did not reveal a
promotable policy fix. On the top-ten CKTSO-gap focus set, fresh auto scale
measured `2.20572s` geomean in
`build/kls_scale_auto_gap10_t4_r1_ref3_timeout120.jsonl`; forced scale `-1`,
`0`, `1`, and `2` measured `2.16409s`, `2.17264s`, `2.35985s`, and
`2.51107s` in the corresponding `build/kls_scale*_gap10_t4_r1_ref3_timeout120`
artifacts. The best forced scale was only about 2% faster than auto and still
left the ASIC and Rajat hard rows far behind CKTSO, so this is not the missing
paper-level algorithm. Forcing the current row-refactor executor was worse:
`--row-refactor refactor` measured `6.91960s` geomean and
`--row-refactor all` measured `5.87436s` in the top-ten focus artifacts, while
`--row-refactor checked` mostly declined the row path and measured `2.29421s`.
Those runs confirm the current auto model is right to keep these rows on the
EGraph path until the row/supernode executor is made coarser.

An adaptive cached-panel disabling patch was also tried and rejected. The idea
was to mark retained panels that actually produced cached updates in one clean
refactor and disable the rest for later refactors, matching the repeated-SPICE
case where the sparsity pattern is fixed. On the top-ten cached focus run,
however, the patch regressed to `2.38995s` geomean in
`build/kls_cached_panel_adapt_gap10_t4_r1_ref3_timeout120.jsonl` from the
precheck-only cached `2.30940s`. A longer top-five `refactor-repeat=20` check
also did not justify it: default measured `1.38623s` in
`build/kls_default_panel_adapt_gap5_t4_r1_ref20_timeout120.jsonl`, while the
adaptive cached path measured `1.42182s` in
`build/kls_cached_panel_adapt_gap5_t4_r1_ref20_timeout120.jsonl`. The source
was reverted. This rejects another per-panel bookkeeping fix; the remaining
path still needs larger row/supernode work units rather than pruning more
single-consumer panel metadata.

A direct SubtreeLU cluster-mode supernode-task probe was also rejected. The
paper points out that CKTSO only uses supernodes in pipeline mode, so KLS tried
an opt-in extension of `KLS_ENABLE_EGRAPH_SUPERNODE_TASKS=1` that could pull a
consecutive EGraph supernode prefix forward from later cluster levels when all
outside U-predecessors were already marked done. This was dependency-safe and
active, but it moved work to earlier workers without changing the scalar
numeric executor: the top-ten CKTSO-gap control measured `2.33474s` geomean in
`build/kls_cluster_supernode_tasks_off_gap10_t4_r1_ref3_timeout120.jsonl`,
while the gated cluster-supernode run measured `2.42881s` in
`build/kls_cluster_supernode_tasks_on_gap10_t4_r1_ref3_timeout120.jsonl`. The
enabled run formed `90,924` supernode tasks covering `229,939` columns, yet it
regressed all substantive rows (`ASIC_320ks` `1.015x`, `ASIC_320k` `1.017x`,
`ASIC_100ks` `1.008x`, `onetone2` `1.063x`, `rajat28` `1.092x`, `rajat20`
`1.104x`, `onetone1` `1.048x`). The source was reverted. This confirms that
merely grouping cluster scheduling tasks is not the missing SubtreeLU piece;
KLS still needs the row-major supernode numeric update/storage that makes the
grouped task perform less scalar work.

The native row-panel selector is now opt-in rather than the default row
refactor subpath. This is a reversal of the earlier "structural work gate by
default" policy because current same-source probes show that KLS's retained
compact panel still adds panel-copy/update overhead before the full paper
numeric engine exists. On the forced-row top-five CKTSO-gap focus, the old
default measured `4.60205s` geomean in
`build/kls_forced_row_native_default_current_gap5_t4_r1_ref3_timeout120.jsonl`,
while `KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=0` measured `4.38362s` in
`build/kls_forced_row_native_off_current_gap5_t4_r1_ref3_timeout120.jsonl`.
The checked-row analogue also favored native-off slightly (`1.48273s` versus
`1.49574s`) in
`build/kls_checked_row_native_off_current_gap5_t4_r1_ref3_timeout120.jsonl` and
`build/kls_checked_row_native_default_current_gap5_t4_r1_ref3_timeout120.jsonl`.
Unset `KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR` now keeps the scalar row-major
dense-group path; `auto` restores the structural work-gated selector for
experiments, and `1` still forces retained compact panels for coverage. This
does not reject row-major panels as the final paper direction; it rejects
promoting the current scalar blocked-panel scaffold before it reduces enough
numeric work.

A current CBLAS-enabled top-ten CKTSO-gap rerun also rejects promoting the
existing BLAS-backed supernode bridge as the missing paper algorithm. The build
configured successfully with `-DKLS_ENABLE_CBLAS_SUPERNODE=ON` in
`build-cblas` and passed `ctest --test-dir build-cblas --output-on-failure`.
With `OPENBLAS_NUM_THREADS=1`, the CBLAS-capable default measured `2.23148s`
geomean in
`build-cblas/kls_cblas_build_default_ob1_gap10_t4_r1_ref3_timeout120.jsonl`,
essentially matching the non-CBLAS top-ten control, but all
`refactor_supernode_cblas_update_*` counters were zero. Enabling the current
EGraph supernode paths did not help: `KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1`
with `KLS_ENABLE_CBLAS_SUPERNODE=1` measured `2.41350s` in
`build-cblas/kls_cblas_egraph_supernode_gap10_t4_r1_ref3_timeout120.jsonl`,
and the narrower cached-only mode measured `2.34182s` in
`build-cblas/kls_cblas_egraph_cached_gap10_t4_r1_ref3_timeout120.jsonl`.
Those runs did build and use panel metadata (`3386` panels, `4654` cached
blocked updates), but still took zero CBLAS updates because the eligible
producer/consumer shapes did not pass the large-work BLAS gate. This confirms
that simply compiling CBLAS or enabling the current cached-panel consumer does
not close the CKTSO refactor gap; KLS still lacks the coarser row-major
supernode numeric executor described by CKTSO/SubtreeLU, where panel packing
and triangular/update work replace enough scalar dependency traversal to pay
for the staging cost.

A direct producer-side EGraph supernode panel-factor probe was also tried and
rejected. The experiment added an opt-in
`KLS_ENABLE_EGRAPH_SUPERNODE_PANEL_FACTOR=1` path that reused retained L-panel
metadata, handled suffix panels, and tested natural, ready-queue, non-ready,
and cluster scheduling insertion points under the no-pivot/unscaled cases.
The code built and passed the smoke tests while present, but it never reached a
numeric panel-factor application on the two slow ASIC rows. The final relaxed
diagnostic run in
`build/kls_panel_factor_diag_probe_gap2_t4_r1_ref3_timeout120.jsonl` recorded
`2,555,032` producer attempts, `2,776` containing-panel hits, and `3,051`
admitted suffixes covering `162,132` rows, yet still produced `0` applied
updates and `0` update rows. The ready-queue forcing path was separately
rejected: enabling only `KLS_ENABLE_EGRAPH_READY_QUEUE=1` measured `27.6083s`
geomean on the two-row ASIC probe in
`build/kls_ready_queue_probe_gap2_t4_r1_ref3_timeout120.jsonl`, versus
`13.0489s` for the no-env probe in
`build/kls_default_probe_gap2_t4_r1_ref3_timeout120.jsonl`.

The conclusion is stronger than "BLAS is too small" or "scheduling is too
fine-grained." The existing KLU-compatible column-oriented U/L storage and
consumer-side cached panels do not expose the row-major supernode factor/update
object assumed by CKTSO/SubtreeLU. The source was reverted. The next paper-gap
implementation should build a true row-major supernode symbolic and numeric
storage layer, then run panel triangular solve and update work from that
storage, instead of wrapping the current scalar column executor or forcing it
through a different queue.

A consumer-side structural-direct EGraph supernode probe was also rejected.
The experiment kept retained panel maps only as shape metadata, skipped
persistent value-panel publication, and read completed producer `L` columns
directly from the current KLU-compatible numeric object while accumulating one
common trailing workspace per accepted run. It built, passed smoke coverage,
and did real grouped work on `ASIC_320ks`, `ASIC_100ks`, and `rajat03`, but
same-binary top-five CKTSO-gap timing regressed: the default measured
`1.39030s` geomean in
`build/kls_struct_direct_default_gap5_t4_r1_ref3_timeout120.jsonl`, while
`KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=direct` measured `1.62085s` in
`build/kls_struct_direct_gap5_t4_r1_ref3_timeout120.jsonl`, with all five rows
slower. `ASIC_320k` formed `401` retained panels but accepted zero structural
runs and still slowed, showing that shape-only retained panel metadata is not
enough. The source was reverted. This reinforces the same paper-level target:
KLS needs a true row-major/supernodal numeric object and batched executor, not
another structural wrapper around KLU column storage.

A successor-amortized native row-panel auto selector was also rejected. The
paper motivation was reasonable: SubtreeLU's compact row-major supernode
storage is supposed to amortize panel packing through later supernode updates,
so `KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=auto` was prototyped to require an
accepted dense row group to have a downstream row-group successor in addition
to the existing arithmetic-intensity gate. Same-session forced-row controls on
the top-ten CKTSO-gap set showed that the retained native panel path is still
mixed rather than promotable: scalar row-major with
`KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=0` measured `6.68536s` geomean in
`build/kls_current_forced_row_native_off_gap10_t4_r1_ref3_timeout120.jsonl`,
while the existing auto selector measured `6.75167s` in
`build/kls_current_forced_row_native_auto_gap10_t4_r1_ref3_timeout120.jsonl`.
The auto path did real paper-shaped work (`onetone1` alone recorded `57`
native panels, `47` blocked panels, and `19,608` compact supernode updates),
but it still regressed that row by `1.083x` and `ASIC_100ks` by `1.014x`,
offsetting wins on `ASIC_320k`, `onetone2`, `rajat25`, and `rajat20`.

The successor-only prototype itself failed smoke coverage before benchmarking:
the dense checked-prefix repair fixture lost its retained compact panel
(`compact=0`) and the subrange batched compact-supernode fixture lost its
batched update counters. That proves row-group successors are not the right
proxy for panel usefulness in KLS's current row engine. Retained panels also
feed same-group/subrange batch executors and checked prefix repair, not only
later group-DAG successors. The source was reverted. A future selector would
need explicit panel-consumer metadata from the compact batch builders and the
checked-prefix repair path; a simple downstream-successor guard is a structural
under-approximation, not a valid paper-gap closure.

The explicit-consumer selector was then prototyped and rejected as well. The
current build and saved gap artifacts already rule out the user's small-BLAS
hypothesis for this path: `build/CMakeCache.txt` has
`KLS_ENABLE_CBLAS_SUPERNODE=OFF`, the recent top-gap JSON rows report
`build_has_cblas=false`, and all CBLAS update counters are zero. The prototype
therefore left the existing BLAS large-case gates unchanged and instead tried
to make native compact row panels conditional on symbolic consumers. It first
counted dense-producer runs, then had to add scalar compact-panel consumers
and consumer-side batched dense-group prefixes to preserve the smoke fixtures;
with those three consumer classes included, `cmake --build build -j$(nproc)`
and `ctest --test-dir build --output-on-failure` passed.

The focused top-ten CKTSO-gap benchmark still rejected the approach. On the
same source, `KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=0` measured `7.83104s`
geomean in
`build/kls_consumer_guard_native_off_gap10_t4_r1_ref3_timeout120.jsonl`, while
the consumer-aware `auto` path measured `7.62489s` in
`build/kls_consumer_guard_native_auto_gap10_t4_r1_ref3_timeout120.jsonl`. That
looked like a local auto/off win, and it did correctly remove the solve-only
`rajat03` native panels (`3` to `0`). But compared with the saved current auto
artifact, most hard rows slowed by about `1.14x` to `1.23x`, with refactor
averages about `1.15x` to `1.26x` higher despite identical native-panel counts
on the ASIC, onetone, and rajat hard rows. The selector's complete metadata
scan cost more than the avoided panels, so the source was reverted. This keeps
the diagnosis pointed at a coarser row-major supernode numeric executor, not at
another panel-retention policy or small-BLAS guard.

A follow-up scalar BTF EGraph cached-index split was also rejected. The measured
hot path on the top CKTSO-gap rows is
`kls_egraph_refactor_btf_unscaled_column`, and those rows already report
`build_has_cblas=false`, zero CBLAS counters, and active 32-bit L/U/map caches.
The prototype split the BTF unscaled dependency loop into a direct 32-bit cached
case and left the existing supernode update path untouched. It built cleanly and
passed `ctest --test-dir build --output-on-failure`, but the same top-ten
CKTSO-gap benchmark moved only from `2.17494s` to `2.14354s` geomean in
`build/kls_i32_btf_fastpath_gap10_t4_r1_ref3_timeout120.jsonl`. Refactor-time
geomean improved by only about `1.1%`, while `rajat25`, `rajat28`, `rajat03`,
and `onetone1` regressed; `rajat25` was the clearest failure at `1.218x`
overall and `1.372x` refactor time.

That is not a robust paper-level fix and is likely CPU/compiler-sensitive. The
source was reverted. The useful conclusion is that small BLAS calls and another
scalar-loop specialization are not the missing CKTSO/SubtreeLU mechanism; KLS
still needs a true row-major supernode numeric object whose triangular solves
and trailing updates replace enough scalar dependency traversal to pay for the
extra storage.

The large paper-recon comparison was refreshed with the local CKTSO library
instead of relying only on the older saved artifact. `build-cktso` was
configured with
`-DKLS_BUILD_CKTSO_COMPARE=ON -DCKTSO_ROOT=/home/ubuntu/sources/kls-workspace/cktso`.
KLS was rerun on `bench/suitesparse_paper_large_recon_manifest.txt` with
`--threads 4 --repeat 1 --refactor-repeat 1 --timeout 120`, producing
`build/kls_fresh_large_recon_t4_r1_ref1_timeout120.jsonl` plus the skip-four
continuation
`build/kls_fresh_large_recon_t4_r1_ref1_timeout120_skip4.jsonl`.
CKTSO was rerun with the matching repeat settings into
`build/cktso_fresh_large_recon_t4_r1_ref1_timeout120.jsonl`.

The completed common large rows are now close on geomean but still mixed:
KLS/CKTSO SPICE-cycle ratios were about `0.25x` on `TSOPF_FS_b39_c30`,
`1.34x` on `nxp1`, `1.07x` on `G3_circuit`, `1.29x` on `ASIC_680k`,
`1.04x` on `rajat29`, and `1.50x` on `rajat30`, for a completed-row geomean
near `0.95x`. This is not evidence of general parity because KLS still timed
out on `pre2` and `Hamrle3`, while the same-session CKTSO run completed
`pre2` and timed out only on `Hamrle3`.

The fresh `pre2` isolation is consistent with the older diagnosis. CKTSO
completed direct `pre2` in the same environment with about `3.89s` analysis,
`7.03s` initial factor, `6.08s` refactor, and a valid residual. KLS
`--analyze-only` completed in about `9.7s` under auto/AMD and reported a
629628-row dominant block, about `61.1M` estimated entries in each factor, and
about `2.08e11` estimated flops; `--refactor-repeat 0` still timed out at
130s. AMD, METIS, SCOTCH, COLAMD, and natural analyze-only probes completed,
but prior and current numeric probes still time out. The actionable gap for
large cases therefore remains cold first-factor numeric machinery, especially
the CKTSO/SubtreeLU row-up-looking dominant-block executor and ETree-tail
pipeline, not another refactor-only scalar split, small-BLAS guard, or simple
ordering selector.

The small-BLAS guard hypothesis was rechecked against the current source while
filling one direct scheduling gap from that diagnosis. The existing CBLAS paths
remain build-time optional, runtime opt-in, and structurally gated at 512-scale
vectors/panels plus multi-million-work thresholds; the normal build used for the
gap runs has `KLS_ENABLE_CBLAS_SUPERNODE=OFF`, and the forced dominant-block
probes below still report `build_has_cblas=false` and zero CBLAS counters.

KLS now lets the BTF-parallel KLS-first row-up factorization route a large
dominant BTF block through the existing intra-block row pipeline while leaving
the fringe blocks on the BTF worker queue. The selector is structural: more
than one BTF block, a largest block of at least 30,000 rows, and at least 75%
matrix coverage by that largest block. The same patch also propagates the
BTF-worker pipeline counters into the public KLS-first stats, because the
previous parallel-BTF success path could run the row pipeline without exposing
that fact.

This is a paper-aligned scaffold improvement, but not the missing CKTSO-scale
fix. `cmake --build build -j$(nproc)` and
`ctest --test-dir build --output-on-failure` passed. Forced KLS-first probes
showed that the dominant-block row pipeline activates:
`ASIC_320k` reported `kls_first_row_pipeline_run_count=2` and
`kls_first_parallel_btf_block_count=798` in
`build/kls_asic320k_dominant_btf_pipeline_stats_t4_factor_timeout45.json`;
`rajat29` reported `kls_first_row_pipeline_run_count=1` and
`kls_first_parallel_btf_block_count=14307` in
`build/kls_rajat29_dominant_btf_pipeline_stats_t4_factor_timeout70.json`.
Both still fell back to KLU on the follow-up factor in the benchmark harness,
and their forced KLS-first cold factors remained slower than the accepted
default policy (`ASIC_320k` about `3.30s`; `rajat29` about `12.61s`). The
current-source `pre2` forced KLS-first factor-only probe still timed out at
130s and left an empty artifact:
`build/kls_pre2_dominant_btf_pipeline_stats_t4_factor_timeout130.json`.

The conclusion is narrower and clearer: guarding BLAS for large cases is
already true for the tested paths, and simply exposing the existing row
pipeline inside the dominant BTF block is not enough. The remaining large gap
is the paper-level numeric representation/executor itself: CKTSO's production
row/supernode storage, checked pivoting-tail scheduler, and coarse updates must
replace more scalar row/panel-cache growth work before `pre2`-class cases can
approach CKTSO.

The next focused rerun kept the small-BLAS conclusion unchanged and removed
one direct serial rebuild from the large whole-block row pipeline. After a
dynamic column pivot in a full-block pipeline of at least 32,768 rows, KLS now
keeps the scalar U entries and prefix supernode metadata correct but drops the
speculative dense prefix panel cache instead of rebuilding it under the
pipeline lock. This is correctness-preserving because the dense panels are only
an acceleration layer over the already swapped U entries; small pipelines keep
the old rebuild path and the existing dynamic-pivot smoke coverage.

`cmake --build build -j$(nproc)` and
`ctest --test-dir build --output-on-failure` passed. The retained change is not
a CKTSO-gap closer. It gave a repeatable modest improvement on the forced
`ASIC_320k` probe (`3.30s` to about `2.97s` initial factor in
`build/kls_asic320k_final_drop_prefix_cache_t4_factor_timeout45.json`) with
the same residual, but `rajat29` was effectively neutral/noisy
(`12.61s` baseline, one `6.83s` transient run, and a final `12.78s` run in
`build/kls_rajat29_final_drop_prefix_cache_t4_factor_timeout70.json`). The
same `pre2` factor-only probe still timed out at 130s and left an empty
artifact, `build/kls_pre2_drop_prefix_cache_t4_factor_timeout130.json`.

GDB sampling after the cache-drop change confirms that the old
`kls_row_first_supernode_panel_cache_build` sample is no longer the only
bottleneck. The updated `pre2` run instead showed workers serialized under the
pipeline lock in `kls_row_first_partial_apply_one_dep`, with another sample in
`kls_row_first_supernodes_reset`. A stronger experiment that disabled all
prefix supernode acceleration after the first large-pipeline pivot was rejected:
it avoided the reset but pushed `rajat29` back near the old slow path. A
prefix-only reset experiment was also rejected for the same reason. The useful
conclusion is that `pre2` needs the paper-level coarse row/supernode executor,
not just less panel-cache rebuilding and not a small-case BLAS guard.

A follow-up unlocked-snapshot experiment was also rejected. The prototype
copied large ready scalar U rows, then ready compact supernode U runs, while
holding the row-pipeline lock and applied the copied updates outside the lock,
discarding and retrying the row if a later dynamic pivot changed the pipeline
epoch. This was correctness-preserving in smoke tests but not a useful
executor replacement. `cmake --build build -j$(nproc)` and
`ctest --test-dir build --output-on-failure` passed, but the focused probes
were mixed: scalar snapshots moved `ASIC_320k` from about `2.97s` to
`3.00s` and `rajat29` from about `12.78s` to `12.12s`; adding compact
supernode-run snapshots moved `ASIC_320k` to about `3.03s` and `rajat29` to
about `12.61s`. Both `pre2` snapshot variants still timed out at 130s and
left empty artifacts,
`build/kls_pre2_snapshot_scalar_t4_factor_timeout130.json` and
`build/kls_pre2_snapshot_supernode_t4_factor_timeout130.json`.

GDB sampling of the snapshot-supernode variant showed the bottleneck simply
moved to `kls_row_first_partial_apply_supernode_run` and worker lock waits;
the extra copy did not create CKTSO-style coarse tasks, it mostly duplicated
memory traffic. That rejects the lightweight snapshot route and strengthens
the next implementation target: build a real production row/supernode numeric
object with scheduled coarse triangular solves/trailing updates, instead of
copying KLU-compatible row fragments around the existing scalar executor.

A selective large-pipeline panel invalidation increment was then retained. The
small-BLAS hypothesis still did not need a source change: the normal build has
CBLAS disabled, and the optional CBLAS paths already require 512-row/vector or
512-width panel minima plus multi-million-work gates before calling BLAS. The
source change instead makes the large dynamic-pivot row pipeline preserve its
existing row-supernode panel cache and mark only panels touching the two
swapped local columns inactive, rather than dropping the whole prefix cache.
Inactive panels are skipped by the cached executor, and the compact/scalar
fallback still rebuilds a valid panel when it next proves the shape.

`cmake --build build -j$(nproc)`,
`ctest --test-dir build --output-on-failure`,
`cmake --build build-cblas -j$(nproc)`, and
`ctest --test-dir build-cblas --output-on-failure` passed. In the matching
forced KLS-first/no-fast-factor first-factor probes, `ASIC_320k` improved from
about `2.97s` to about `2.20s` and `rajat29` from about `12.78s` to about
`4.63s`, with matching residuals, in
`build/kls_asic320k_selective_panel_invalidate_nofast_t4_factor_timeout45.json`
and
`build/kls_rajat29_selective_panel_invalidate_nofast_t4_factor_timeout70.json`.
The production default path was neutral within repeated-run noise on
`ASIC_320k` (`13.61s` saved current policy versus a rerun at `13.67s`) and
slightly noisy on `rajat29` (`5.52s` saved large-recon policy versus `5.82s`).
The hard `pre2` forced KLS-first factor probe still timed out at 130s and left
an empty
`build/kls_pre2_selective_panel_invalidate_t4_factor_timeout130.json`, so this
is a useful storage-level cleanup but still not the missing paper-scale
row/supernode executor.

A `pre2` METIS forced KLS-first rerun exposed one large non-numeric overhead in
that path. METIS analysis builds a retained separator forest for `pre2`
(`659033` analyzed rows, `58569` components, and a `629628`-row dominant BTF
block in `build/kls_pre2_metis_analyze_current_t4.json`), but the first
bounded numeric probe timed out before useful JSON. An interrupting GDB run
showed the process still in `kls_record_first_separator_queue_plans`, repeatedly
building an `n`-row symbolic map for every one of the `29282` BTF blocks just
to record queue diagnostics. KLS now keeps that diagnostic all-block behavior
for ordinary BTF counts, but when a huge BTF forest has a dominant large block
it records only the dominant or otherwise large executable blocks. The executor
still rebuilds the actual selected block plan before running it.

`cmake --build build -j$(nproc)` and
`ctest --test-dir build --output-on-failure` passed. The same interrupting GDB
probe moved past the diagnostic loop into the real private row factorization:
active worker samples were in `kls_row_first_supernode_panel_cache_append` and
`kls_row_first_partial_apply_supernode_run_cached`, with the parent waiting in
`kls_row_first_run_parallel_private_phase`
(`build/pre2_metis_gdb_after_queue_diag_gate.txt`). The fix therefore removes a
large setup artifact from the `pre2` METIS experiment, but it does not close the
CKTSO gap: `build/kls_pre2_metis_queue_diag_gate_t4_factor_timeout90.json` and
`build/kls_pre2_metis_queue_diag_gate_t4_factor_timeout130.json` both remain
empty after timing out. Focused completed-row checks stayed residual-clean:
`ASIC_320k` still exercised the separator queue with a valid residual in
`build/kls_asic320k_queue_diag_gate_nofast_t4_factor.json`, and `rajat29`
remained on the non-separator KLS-first path with a valid residual in
`build/kls_rajat29_queue_diag_gate_nofast_t4_factor.json`. The remaining
`pre2` gap is now clearer: after setup reaches the paper-aligned separator
private phase, KLS is still dominated by scalar/panel-cache row-supernode
numeric work rather than by the diagnostic queue planner.

The next paper-aligned separator-tree probe tested whether KLS was retaining too
shallow a METIS `NodeNDP` tree for the KLS-first separator executor. SubtreeLU
Algorithm 2 sets the minimum depth to `log2(P)` but continues partitioning until
subdomains are below `max(200, N/1000)`. KLS previously passed exactly the
numeric thread count as `NodeNDP` `npes`, so four-thread `pre2` retained only
four leaves and a `217028`-row largest private component. Analyze-only sweeps
showed that increasing the retained leaves did not materially change METIS
analysis time on `pre2` but strongly improved queue shape: 32 leaves reduced the
largest private component to `44839` rows, 256 leaves to `7826` rows, and the
paper-threshold final 999-leaf tree to `3718` rows.

An unconditional deeper METIS tree was rejected for the production/default path:
completed large METIS/KLU rows kept valid residuals but showed mixed timing
regressions (`G3_circuit` and `rajat30` slowed modestly), because the
permutation itself changes. The first forced-KLS-first version was also too
aggressive in the paper threshold's 200-row saturated regime: the 30.6k-row
weak-pivot separator smoke fixture exposed a residual around `7e-7`, matching
the SubtreeLU warning that fine-grained separator partitioning can restrict
pivoting choices in factorization. The retained change is therefore scoped to
forced KLS-first METIS analyses with at least 200k rows, where `N/1000` rather
than the hard 200-row floor controls the target leaf size. Default METIS
analysis remains at the old thread-count tree unless
`KLS_ENABLE_KLS_FIRST_FACTOR=1` is set.

Current-source checks passed `cmake --build build -j$(nproc)` and
`ctest --test-dir build --output-on-failure`. Paired `pre2` analyze-only probes
confirm the gate: `--kls-first-factor off` retained the old 4-leaf tree with
`217028` largest private rows in
`build/kls_pre2_metis_analyze_default_ndp_t4_final.json`, while
`--kls-first-factor on` retained 999 leaves with `3718` largest private rows in
`build/kls_pre2_metis_analyze_klsfirst_deep_ndp_t4_final.json`. A completed
forced-KLS-first `ASIC_320k` probe remained residual-clean with the deeper tree
in `build/kls_asic320k_klsfirst_deep_ndp_final_t4_factor.json`. The change still
does not close the hard `pre2` numeric gap:
`build/kls_pre2_metis_klsfirst_deep_ndp_final_t4_factor_timeout120.json` is
empty after the 120s timeout. The remaining bottleneck is still the actual
row/supernode numeric executor, especially worker-local panel-cache copying and
scalar/panel-cache updates, not the retained separator-tree depth alone.

A follow-up on the BLAS hypothesis confirmed again that small BLAS calls are not
the active slow-case mechanism. The normal build remains CBLAS-free unless the
optional CMake flag and `KLS_ENABLE_CBLAS_SUPERNODE=1` runtime gate are both
used, and the CBLAS path already has 512-scale plus work/copy gates. The current
`pre2` samples were in native row-first storage, not in CBLAS. The accepted
change instead addresses the sampled native storage bridge directly: a missed
row-supernode panel cache append is now bounded to at most 4M retained entries
per worker before it falls through to the direct compact/scalar update. Already
cached panels and completed-run panel publishing are unchanged. This avoids the
observed first-use path growing a transient worker panel cache from about 32 MB
to 64 MB before updating the current row.

The final same-source checks passed `cmake --build build -j$(nproc)` and
`ctest --test-dir build --output-on-failure`. `ASIC_320k` forced KLS-first
remained residual-clean in
`build/kls_asic320k_missed_panel_cap_final_t4_factor.json` with
`2.20894768s` initial factor, `2.61575448e-15` relative residual, and
`2520001` row panel-cache append entries. The hard `pre2` forced KLS-first
factor run still timed out at 120s with an empty
`build/kls_pre2_metis_missed_panel_cap_t4_factor_timeout120.json`, so the cap is
a storage-copy cleanup rather than the missing CKTSO-scale executor. The GDB
sample moved from `realloc` inside `kls_row_first_supernode_panel_cache_append`
to `kls_row_first_u_rows_fit_supernode`/private row work in
`build/pre2_metis_missed_panel_cap_gdb_interrupt.txt`, which keeps the next
large gap focused on repeated structural supernode validation and scalar
row-update work.

An owner-local private supernode-end cache was also prototyped and rejected. It
cached adjacent-row fit decisions for private separator workers, but a completed
`ASIC_320k` check slowed from the missed-panel-cap control's `2.21615628s`
initial factor to about `2.438s` while `pre2` still timed out at 120s
(`build/kls_asic320k_private_supernode_cache_t4_factor.json`,
`build/kls_asic320k_private_supernode_cache_norebuild_t4_factor.json`, and
`build/kls_pre2_metis_private_supernode_cache_t4_factor_timeout120.json`). The
prototype also exposed that exact adjacent-row validation still has to scan long
tails at least once per candidate pair, so this is not the direct paper-level
replacement for the missing coarse supernode numeric executor.

The narrower retained follow-up memoizes only the adjacent-row structural
predicate inside each private row worker. The cache is tri-state per adjacent
pair, is filled lazily by the existing `kls_row_first_u_rows_fit_supernode`
predicate, and is dropped after a dynamic pivot instead of being rebuilt. The
predicate also now rejects unequal tail lengths before entering the elementwise
tail comparison. This keeps behavior tied to the existing structural test while
avoiding repeated long scans for the same private-worker row pair.

`cmake --build build -j$(nproc)`, `ctest --test-dir build --output-on-failure`,
and `git diff --check` passed. A same-session `ASIC_320k` forced KLS-first
control built from commit `2773488` measured `2.49741041s` initial factor in
`build/kls_asic320k_cap_baseline_same_session_t4_factor.json`; the current
memoized build measured `2.29454936s` and `2.32290521s` in
`build/kls_asic320k_fit_memo_t4_factor.json` and
`build/kls_asic320k_fit_memo_after_baseline_t4_factor.json`, all with
`2.61575448e-15` relative residual. The hard `pre2` forced KLS-first factor run
still timed out at 120s with an empty
`build/kls_pre2_metis_fit_memo_t4_factor_timeout120.json`, but the interrupting
GDB sample moved from `kls_row_first_u_rows_fit_supernode` to
`kls_row_first_partial_apply_supernode_run_cached` / `kls_row_first_heap_pop`
inside the private row worker
(`build/pre2_metis_fit_memo_gdb_interrupt.txt`). The remaining visible gap is
therefore the cached row-supernode update itself and its heap/state movement,
not repeated adjacent-tail validation.

The next retained row-supernode executor cleanup removes another transient heap
operation from validated same-run updates. When a cached or compact
row-supernode run introduces a later dependency that is still inside the same
run, KLS no longer pushes that dependency into the row dependency heap just to
pop it before applying the next row of the run. A small consume helper still
pops dependencies that were already present in the heap and still rejects a
smaller unexpected root, so the change is limited to dependencies created by
the current validated run.

The same build and smoke checks passed again. On `ASIC_320k` forced KLS-first,
the post-memoization control in
`build/kls_asic320k_fit_memo_after_baseline_t4_factor.json` measured
`2.32290521s` initial factor, while the same-run heap-elision build measured
`2.27882813s` with the same `2.61575448e-15` relative residual in
`build/kls_asic320k_same_run_heap_elide_t4_factor.json`. The hard `pre2`
forced KLS-first run still timed out at 120s
(`build/kls_pre2_metis_heap_elide_t4_factor_timeout120.json`), but the GDB
sample no longer showed `kls_row_first_heap_pop`; it moved into
`kls_row_first_partial_apply_supernode_run` itself
(`build/pre2_metis_heap_elide_gdb_interrupt.txt`). The remaining visible gap is
now inside the cached supernode arithmetic/state update body rather than
avoidable heap maintenance for same-run internal dependencies.

The follow-up BLAS-size question was rerun before changing the executor. The
CBLAS-capable top-five CKTSO-gap probe measured `1.40471058s` geomean with
`KLS_ENABLE_CBLAS_SUPERNODE=0` and `1.52710910s` with
`KLS_ENABLE_CBLAS_SUPERNODE=1`
(`build-cblas/kls_cblas_guard_off_current_top5_t4_r1_ref3_timeout120.jsonl`
and `build-cblas/kls_cblas_guard_on_current_top5_t4_r1_ref3_timeout120.jsonl`).
Both runs reported `build_has_cblas=true` but zero CBLAS, compact GEMV/TRSV,
blocked-panel, and KLS-first panel-cache counters on every row. The existing
512-scale and work-per-copy CBLAS guards are therefore not the active
slow-case mechanism; the normal build is not calling BLAS in these paths.

The retained change instead batches row-supernode executor storage and avoids
tiny retained-panel construction. KLS now reserves worker entry arrays by the
whole append count, copies worker-local entries in bulk when no column linked
list is needed, and stores the solved consecutive L multipliers from cached
CBLAS/native panel runs with one consecutive append. Missed row-panel cache
materialization now uses the same cached-supernode structural/work gate as the
EGraph cache before retaining a panel, so width-2 style producer runs stay on
the native compact/scalar path instead of paying dense-panel copy and lookup
costs. Dynamic-pivot panel-cache reset also keeps the allocated cache and
clears only retained panel rows, avoiding an O(block rows) reinitialization of
`panel_id_by_row` on large blocks.

The verification set passed `cmake --build build -j4`,
`ctest --test-dir build --output-on-failure`,
`cmake --build build-cblas -j4`,
`ctest --test-dir build-cblas --output-on-failure`, and `git diff --check`.
`ASIC_320k` forced KLS-first stayed residual-clean: the previous accepted
heap-elision artifact measured `2.27882813s` initial factor with
`2.61575448e-15` relative residual
(`build/kls_asic320k_same_run_heap_elide_t4_factor.json`), while the current
work-gated/batched run measured `0.960101378s` initial factor and
`1.50913325e-15` relative residual
(`build/kls_asic320k_panel_clear_batch_l_t4_factor.json`). The top-five normal
CKTSO-gap focus run completed with no failures at `1.44303173s` geomean in
`build/kls_panel_clear_batch_l_gap5_t4_r1_ref3_timeout120.jsonl`.

The hard `pre2` forced KLS-first METIS factor probe still timed out at 120s
with an empty
`build/kls_pre2_metis_panel_clear_batch_l_t4_factor_timeout120.json`. The GDB
sequence is useful: before the panel gate, workers were still in
`kls_row_first_supernode_panel_cache_append` and
`kls_row_first_supernode_panel_cache_init_rows`
(`build/pre2_metis_batch_l_append_gdb_interrupt.txt` and
`build/pre2_metis_panel_workgate_batch_l_gdb_interrupt.txt`). After the reset
cleanup, the visible stack moved to `kls_row_first_partial_apply_one_dep` and
`kls_row_first_partial_apply_supernode_run_compact`
(`build/pre2_metis_panel_clear_batch_l_gdb_interrupt.txt`). That narrows the
remaining CKTSO-scale gap to the native scalar/compact row-dependency executor,
not BLAS calls, panel-cache growth, or full-row cache resets.

The same small-BLAS-guard hypothesis was rerun on the current source before
adding any new BLAS policy. The CBLAS-capable top-five CKTSO-gap focus measured
`1.52323303s` geomean with `KLS_ENABLE_CBLAS_SUPERNODE=0` and `1.54428317s`
with `KLS_ENABLE_CBLAS_SUPERNODE=1`
(`build-cblas/kls_cblas_guard_off_shape_top5_t4_r1_ref3_timeout120.jsonl` and
`build-cblas/kls_cblas_guard_on_shape_top5_t4_r1_ref3_timeout120.jsonl`).
Both artifacts reported `build_has_cblas=true` but zero CBLAS, compact
GEMV/TRSV, blocked-panel, and KLS-first panel-cache counters on every row. KLS
already has the requested "BLAS only for large cases" behavior on the exercised
paths: CBLAS calls require the runtime gate plus 512-scale shape tests and
work/copy thresholds before entering BLAS. The current slow case is not small
BLAS dispatch.

The retained executor change instead reuses validated row-supernode shape
information. Private and ready row-supernode runs are formed only after adjacent
rows have passed the exact supernode predicate, so the compact consumer no
longer rescans internal and tail column identities for those trusted runs.
Bounds, pivot, length, and fallback checks remain in place. `cmake --build
build -j4`, `ctest --test-dir build --output-on-failure`, `cmake --build
build-cblas -j4`, `ctest --test-dir build-cblas --output-on-failure`, and
`git diff --check` passed. `ASIC_320k` forced KLS-first stayed residual-clean in
the current shape-trust run
(`build/kls_asic320k_trusted_shape_compact_final_t4_factor.json`) with
`0.999742177s` initial factor and `1.50913325e-15` relative residual; an
earlier same-change repeat measured `0.963871569s` in
`build/kls_asic320k_trusted_shape_compact_run2_t4_factor.json`. `pre2` still
timed out at 120s
(`build/kls_pre2_metis_trusted_shape_compact_t4_factor_timeout120.json`), but
the GDB sample moved from compact shape validation to compact arithmetic with
`shape_known=1` (`build/pre2_metis_trusted_shape_compact_gdb_interrupt.txt`).
That leaves the next paper-level issue in numeric row-supernode arithmetic and
static private work distribution, not in BLAS call size.

Two larger-looking follow-ups were rejected before the retained cleanup. First,
a dependency-weighted separator-private work model added an ordered lower-edge
update proxy on top of the existing `row_input_count^2` component weight. This
matched the papers' work-balance language more closely, but it mostly added
planning cost on the current top-five rows, where the separator queue was not
active, and did not close `pre2`: the forced/default top-five probes measured
about `1.740s`/`1.742s` geomean
(`build/kls_depwork_gap5_t4_r1_ref3_timeout120.jsonl` and
`build/kls_depwork_default_gap5_t4_r1_ref3_timeout120.jsonl`), while
`build/kls_pre2_metis_depwork_t4_factor_timeout120.json` remained empty after
120s. Second, delaying compact-run L publication to append the whole
consecutive multiplier segment at once broke a separator pivot-epoch smoke
case with a huge residual, so compact runs must keep publishing each multiplier
before the rest of the row machinery can observe the updated state.

The retained storage cleanup keeps that publication order. Scalar, cached, and
compact row-supernode run loops already reserve the whole L segment before the
loop; they now write each dependency with the reserved-entry helper instead of
re-entering the growable append helper. This removes redundant capacity growth
checks without changing when each L entry becomes visible. `cmake --build build
-j4`, `ctest --test-dir build --output-on-failure`, `cmake --build build-cblas
-j4`, `ctest --test-dir build-cblas --output-on-failure`, and `git diff
--check` passed. `ASIC_320k` forced KLS-first/no-fast remained residual-clean
in `build/kls_asic320k_reserved_l_run_nofast_t4_factor.json`
(`2.35403384s` initial factor, `2.04287279e-15` relative residual). The hard
`pre2` probe still timed out at 120s
(`build/kls_pre2_metis_reserved_l_run_t4_factor_timeout120.json`), so this is
not the missing CKTSO-scale executor. The focused top-five default rerun
completed with no failures and `1.39715503s` geomean in
`build/kls_reserved_l_run_gap5_rerun_t4_r1_ref3_timeout120.jsonl`, after a
noisier first pass measured `1.59194508s` in
`build/kls_reserved_l_run_gap5_t4_r1_ref3_timeout120.jsonl`.

The next `pre2` inspection confirmed that a small-BLAS guard is not the active
problem. The CBLAS-capable top-five probes already reported zero CBLAS,
compact GEMV/TRSV, blocked-panel, and row-panel counters with the runtime gate
both off and on, so KLS is already guarding optional BLAS behind large-shape
and work/copy thresholds on these paths. A GDB interrupt of the current
`pre2` METIS KLS-first factor run instead showed one active worker in
`kls_row_first_supernodes_reset` while other pipeline workers waited on the
pipeline condition (`build/pre2_metis_reserved_l_run_gdb_interrupt.txt`).

KLS now avoids that full reset for active-rank pipeline rows that did not
perform a dynamic column pivot. Such rows publish only their completed row's
supernode state with `kls_row_first_supernodes_publish_row`; pivoted rows still
take the full reset and private panel-cache rebuild because column exchanges
can invalidate trusted row-supernode ranges. `cmake --build build -j4`,
`ctest --test-dir build --output-on-failure`, `cmake --build build-cblas -j4`,
`ctest --test-dir build-cblas --output-on-failure`, and `git diff --check`
passed. A fresh `ASIC_320k` forced KLS-first/no-fast run remained
residual-clean in
`build/kls_asic320k_active_rank_nonpivot_publish_nofast_t4_factor.json`
(`2.12963904s` initial factor, `2.04287279e-15` relative residual).

This retained change does not close `pre2`: the matching 120s factor-only
probe still left an empty
`build/kls_pre2_metis_active_rank_incremental_supernode_t4_factor_timeout120.json`,
and the post-change interrupt still sampled `kls_row_first_supernodes_reset`
through the dynamic-pivot branch
(`build/pre2_metis_active_rank_incremental_supernode_gdb_interrupt.txt`). A
stronger prototype that skipped the pivot reset by invalidating only affected
panels and forcing later compact shape revalidation was rejected: it moved the
`pre2` sample to `kls_row_first_partial_apply_supernode_run`, but a bounded
`rajat29` no-fast probe failed setup as singular
(`build/kls_rajat29_active_rank_shape_invalidate_nofast_t4_factor_timeout90.json`).
Restoring the full pivot reset made the same `rajat29` probe complete with
`1.06772992e-13` relative residual in
`build/kls_rajat29_active_rank_nonpivot_publish_nofast_t4_factor_timeout90.json`.
The remaining paper-level gap is therefore still the checked pivoting
row-supernode executor and coarse numeric storage, not a BLAS threshold.

A fresh CBLAS build check on `ASIC_320k` reconfirmed the same conclusion after
the active-rank reset change. With `OPENBLAS_NUM_THREADS=1` and
`KLS_ENABLE_CBLAS_SUPERNODE=0`,
`build-cblas/kls_asic320k_cblas_off_current_t4_factor.json` measured
`2.31691659s` initial factor and `1.50913325e-15` relative residual. The
same binary with `KLS_ENABLE_CBLAS_SUPERNODE=1` in
`build-cblas/kls_asic320k_cblas_on_current_t4_factor.json` measured
`2.30408845s` initial factor with the same residual. Both artifacts reported
`build_has_cblas=true`, but all external CBLAS, compact GEMV/TRSV, blocked
panel, and KLS-first row-panel counters stayed at zero. This directly rejects
the current "guard BLAS for only large cases" hypothesis: the guard already
requires 512-scale vector/panel shapes plus multi-million-work thresholds, and
the tested gap path is not entering BLAS at all.

The follow-up `rajat29` rerun found a fallback hygiene issue rather than a
small-BLAS issue. Earlier failed KLS-first probes restored only
`status`/`numerical_rank`/`singular_col`/`noffdiag` before entering the KLU
fallback. KLS now restores the full KLU common state after failed row-up-looking
and pivoted-block first-factor attempts, so failed experimental paths cannot
leave changed scale, tolerance, memory-growth, or allocation counters behind.
`cmake --build build -j4`, `ctest --test-dir build --output-on-failure`, and
`git diff --check` passed. Forced KLS-first/no-fast `rajat29` probes completed
at 2, 3, and 4 threads in
`build/kls_rajat29_common_restore_retry_nofast_t2_factor_timeout90.json`,
`build/kls_rajat29_common_restore_retry_nofast_t3_factor_timeout90.json`, and
`build/kls_rajat29_common_restore_retry_nofast_t4_factor_timeout90.json`
(`4.76019409s`, `5.530147s`, and `5.78516555s` initial factor, each with
`9.86042515e-12` relative residual). A matching `ASIC_320k` no-fast probe
completed in
`build/kls_asic320k_common_restore_retry_nofast_t4_factor_timeout60.json`
(`2.09658007s`, `2.02952977e-15` residual). These runs still report
`last_factor_path="klu_fallback"`, zero KLS-first row-up-looking work, and
`build_has_cblas=false`, so this is a correctness cleanup for fallback recovery,
not the paper-level CKTSO gap closer. A broader symbolic `Q`/`Lnz` restore
prototype was rejected because it made the threaded `rajat29` repro fail setup
as singular again.

The same fallback-cleanup source was also checked in the CBLAS-enabled tree.
`cmake --build build-cblas -j4` and
`ctest --test-dir build-cblas --output-on-failure` passed. With
`OPENBLAS_NUM_THREADS=1`, `ASIC_320k` forced KLS-first/no-fast measured
`2.31416104s` with `KLS_ENABLE_CBLAS_SUPERNODE=0` and `2.29213644s` with
`KLS_ENABLE_CBLAS_SUPERNODE=1` in
`build-cblas/kls_asic320k_common_restore_cblas_off_t4_factor.json` and
`build-cblas/kls_asic320k_common_restore_cblas_on_t4_factor.json`. Both runs
reported `build_has_cblas=true` and zero external CBLAS update counters, so the
current-source BLAS guard conclusion is unchanged.

The row-first supernode diagnostics were strengthened because the volatile
`kls_first_last_row_supernode_*` fields can be cleared by later fallback or
refactor bookkeeping even when the KLS-first attempt did use the row-supernode
executor. KLS now reports cumulative
`kls_first_row_supernode_update_groups`/`rows` and the panel-backed subset in
`kls_first_row_supernode_panel_update_groups`/`rows`. On a focused
`ASIC_320k` forced KLS-first factor-only run,
`build/kls_asic320k_rowstats_final_t4_factor.json` measured `2.25258371s`
initial factor with `1.50913325e-15` relative residual and reported `93365`
row-supernode groups over `978512` rows. The panel-backed subset was `42656`
groups over `779428` rows, about `18.27` rows per panel group. This rejects the
idea that the completed `ASIC_320k` row-first path is dominated only by
width-2/3 producer runs; the next CKTSO/SubtreeLU gap remains the arithmetic
and synchronization cost of the existing row-supernode executor.

A direct portable cached-panel row-major accumulation prototype was also
rejected. It accumulated dense-suffix and tail updates into worker scratch with
producer-row-major loops before scattering, matching the shape of a native
`gemv`, but on the same `ASIC_320k` controls it regressed initial factor time
from `2.30488623s` to `2.38630384s` for forced KLS-first and from
`2.29399475s` to `2.38004146s` for forced KLS-first/no-fast, with unchanged
residuals and nearly identical group counts. The source was reverted, leaving
only the cumulative diagnostics.

The small-BLAS hypothesis was rechecked once more on the current CBLAS-capable
top-five CKTSO-gap focus. The default runtime switch comparison measured
`1.49142089s` geomean with `KLS_ENABLE_CBLAS_SUPERNODE=0` and `1.39728323s`
with `KLS_ENABLE_CBLAS_SUPERNODE=1` in
`build-cblas/kls_cblas_guard_off_current_top5_t4_r1_ref3_p1_timeout120.jsonl`
and
`build-cblas/kls_cblas_guard_on_current_top5_t4_r1_ref3_p1_timeout120.jsonl`,
but both runs reported zero EGraph panel, blocked, cached-probe, and CBLAS
update counters. Forcing the only relevant bridge with
`KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=cached` and
`KLS_ENABLE_CBLAS_SUPERNODE=1` measured `1.48961453s` in
`build-cblas/kls_cblas_cached_current_top5_t4_r1_ref3_p1_timeout120.jsonl`;
that run built panels and applied native blocked updates, but still reported
zero CBLAS updates. This leaves no useful "only use BLAS for large cases" patch
to retain: the large-case CBLAS gates already exist, and the active CKTSO-gap
path either does not enter the panel bridge or uses native blocked panel work.
The remaining paper-level gap is still the row/supernode numeric executor and
its coarse producer/consumer scheduling, not another BLAS threshold.

A post-rebuild same-binary rerun kept that conclusion. After rebuilding both
`build/` and `build-cblas/`, the top-five CKTSO-gap focus with
`OPENBLAS_NUM_THREADS=1` measured `1.66028753s` geomean with
`KLS_ENABLE_CBLAS_SUPERNODE=0` in
`build-cblas/kls_cblas_guard_fresh_off_gap5_t4_r1_ref3_timeout120.jsonl` and
`1.40658849s` with `KLS_ENABLE_CBLAS_SUPERNODE=1` in
`build-cblas/kls_cblas_guard_fresh_on_gap5_t4_r1_ref3_timeout120.jsonl`.
Every row in both files reported `build_has_cblas=true` but zero numeric
CBLAS, EGraph-panel, blocked-panel, and row-panel-cache update counters. The
short-run timing difference is therefore not evidence of BLAS work being used;
it is run noise or unrelated scheduling variance. KLS should not add another
small-BLAS guard here because the existing large-shape guard is already
stricter than the proposed fix and the active path does not reach it.

The next scalar row-supernode executor probes were also rejected. A hot-stack
sample of current forced-METIS `pre2` again found one pipeline worker in
`kls_row_first_partial_apply_one_dep` while peer workers waited on the pipeline
condition, so KLS tried splitting the scalar U-row update into dependency and
trailing-column loops. That preserved residuals but regressed the same-session
top-five CKTSO-gap focus from `1.46905424s` to `1.52124999s` geomean in
`build/kls_scalar_dep_control_gap5_t4_r1_ref3_timeout120.jsonl` and
`build/kls_scalar_dep_split_gap5_t4_r1_ref3_timeout120.jsonl`; the source was
reverted.

A broader synchronization prototype then moved a pipeline row's final
dependency application and pivot check outside the ordered mutex once the row
reached its commit slot, disabling shared panel-cache mutation for that
outside-lock final pass. This matched the observed wait shape more directly,
but the evidence was too weak for a default change: three-pass top-five medians
were effectively neutral (`1.41735016s` control versus `1.41398595s` candidate
in `build/kls_pipeline_unlock_control_gap5_t4_r1_ref3_p3_timeout120.jsonl` and
`build/kls_pipeline_unlock_final_gap5_t4_r1_ref3_p3_timeout120.jsonl`), while
forced `ASIC_320k` KLS-first factor-only regressed from `3.92563111s` to
`4.01844172s` SPICE-cycle median in
`build/kls_asic320k_pipeline_unlock_control_t4_factor_p3.jsonl` and
`build/kls_asic320k_pipeline_unlock_final_t4_factor_p3.jsonl`. An interrupt of
the unresolved `pre2` run moved the sampled active worker to
`kls_row_first_supernodes_reset`, but that stack movement did not come with a
completed timing win, so the source was reverted.

A safer active-rank pivot-reset cleanup was also rejected. Instead of skipping
the full pivot reset, it reset supernode metadata only for rows already marked
done and made panel-cache rebuild scans skip not-done rows before reading
supernode metadata. This remained residual-clean on smoke tests but regressed
the three-pass top-five focus to `1.49098394s` in
`build/kls_done_row_reset_gap5_t4_r1_ref3_p3_timeout120.jsonl`, and forced
`ASIC_320k` KLS-first factor-only to `4.13749417s` in
`build/kls_asic320k_done_row_reset_t4_factor_p3.jsonl`. The source was
reverted. The retained conclusion is narrower: the current gap is not the cost
of initializing not-done supernode rows, and lock-scope changes need a true
row/supernode numeric representation that avoids shared panel-cache mutation
rather than temporarily bypassing it.

A retained-input batch-row probe was also rejected. The row-refactor symbolic
builder already creates target maps for `KLS_ROW_REFACTOR_GROUP_BATCH`, but the
process-row loader intentionally declined those rows because L-input targets in
that path must remain in the scalar workspace for existing dependency and
compact-supernode update code. A corrected probe allowed batch rows through the
target map while mapping L targets back into `x[dep]`; this raised last-run
targeted rows substantially (`ASIC_320ks` about 217k to 321k and `ASIC_100ks`
about 68k to 99k), but it did not make the executor faster. Same-session
forced-row top-five checks with `KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=0`
regressed `ASIC_320ks` from `0.2901s` to `0.2934s`, `ASIC_320k` from
`0.3109s` to `0.3613s`, `gemat12` from `0.000907s` to `0.000988s`,
`rajat03` from `0.00227s` to `0.00290s`, and `ASIC_100ks` from `0.1439s` to
`0.1529s` average refactor time. A longer two-row repeat-20 check likewise
regressed `ASIC_320ks` from `0.2474s` to `0.2667s` and `ASIC_100ks` from
`0.1210s` to `0.1410s`. The source was reverted. The useful conclusion is that
retained input-target coverage is not the missing paper lever by itself; batch
rows still need a coarser producer/consumer numeric kernel instead of more
target-map admission into the scalar row loop.

A direct EGraph supernode wait-tail probe was rejected after profiling the
current default CKTSO-gap path. `perf` was unavailable in the container
(`perf_event_paranoid=4`), so a long `ASIC_320k` GDB interrupt sampled the
default four-thread refactor with two workers inside
`kls_scatter_subtract_i32` from `kls_egraph_refactor_btf_unscaled_column`, two
workers waiting in `kls_egraph_refactor_wait_done` on the same predecessor
column, and the parent blocked in `kls_egraph_mapped_refactor`. The paper
motivation was reasonable: exact adjacent supernode chains give a transitive
done relation, so a consumer could wait for the last column in a consecutive
dependency run and then consume the scalar L scatters without per-column waits.
The prototype did exactly that without enabling the slower cached-panel
numeric path, and it was residual-clean on smoke and a focused `ASIC_100ks`
run, but it regressed `ASIC_100ks` refactor average from the saved current
`0.0382s` to `0.0441s` while recording about `101824` wait runs over `966694`
columns in the last run. The source was reverted. The useful conclusion is
that tail waiting removes overlap: consumers can start scattering the first
ready predecessor while later columns in the same supernode chain are still
finishing. This also explains why the earlier producer-side supernode-task
artifact (`KLS_ENABLE_EGRAPH_SUPERNODE_TASKS=1`) was not a general fix. Closing
the CKTSO gap still requires a real coarse row/supernode numeric task, not only
coarser waits over the existing scalar scatter loop.

A retained diagnostic now measures that missing row/supernode numeric target
directly instead of inferring it from rejected panel probes. With
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_STATS=1`, schedule construction builds
a temporary map from every column inside a retained EGraph supernode to that
supernode's end, then scans the actual U dependency streams for contiguous
consumer runs that a SubtreeLU-style row-major panel could consume. This is
only a counter path; it does not change the default executor. On `ASIC_100ks`,
the current top-gap default refactor had `1,556,952` dependency edges and the
diagnostic found `160,299` supernode-shaped consumer runs covering `1,148,099`
dependency rows, with max width `282`, `306,222,715` covered L-entry visits,
and `19,687,996` internal dense entries
(`build/kls_supernode_consumer_asic100ks_t4_r1_ref3.json`). On `ASIC_320k`, it
found `72,446` runs covering `939,167` rows, max width `537`,
`303,544,843` L-entry visits, and `40,628,085` internal entries
(`build/kls_supernode_consumer_asic320k_t4_r1_ref1.json`). The top-five focus
artifact `build/kls_supernode_consumer_gap5_t4_r1_ref3_timeout120.jsonl` kept
residuals clean and showed the same pattern on `ASIC_320ks`, `rajat03`, and
`ASIC_100ks`; `gemat12` remains a first-factor loss with no EGraph refactor
edges.

The comparison against the existing cached-panel executor is the important
part. On the same `ASIC_100ks` source,
`KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=cached` accepted only `267` cached updates
over `39,436` rows and `7,400,641` blocked update entries
(`build/kls_supernode_cached_asic100ks_current_t4_r1_ref3.json`), while the
new consumer-run diagnostic saw more than a million supernode-shaped
dependency rows. The missing paper algorithm is therefore not absent supernode
structure and not just a wait or BLAS threshold. It is the durable
row-major/supernode numeric object that can consume these runs without
requiring the current common-trailing cached-panel shape or
rebuilding/scattering through KLU-compatible column storage for each consumer.

A direct ragged EGraph supernode consumer was then tried and rejected. The
prototype added an opt-in `KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=ragged` mode,
skipped the common-trailing panel cache, validated each retained contiguous
run, solved the internal triangular part in worker scratch, and scattered each
producer column's ragged tail directly from the KLU-compatible L columns. It
was residual-clean and applied substantial grouped work, which confirms the
consumer-run diagnostic: with a 16-row gate, `ASIC_100ks` applied `14,193`
runs over `558,741` rows and `213,794,108` update entries, but refactor time
was `0.0934s` versus the same-session default `0.0413s` in the top-five screen
(`build/kls_ragged_supernode_gap5_t4_r1_ref3_timeout120.jsonl` and
`build/kls_ragged_control_gap5_t4_r1_ref3_timeout120.jsonl`). Tightening the
gate to 64-row runs still regressed a serial same-matrix check:
`ASIC_100ks` default refactor was `0.03796s`, while ragged-wide mode applied
`2,074` runs over `210,696` rows and `77,517,080` entries but took `0.06573s`
(`build/kls_ragged64_control_asic100ks_t4_r1_ref3.json` and
`build/kls_ragged64_asic100ks_t4_r1_ref3.json`). The source was reverted.
This is a useful negative result: simply replacing scalar scatter with a
per-consumer ragged validation plus branchy direct L-column scatter does not
close the paper gap. The durable row-major/supernode numeric object needs to
precompute the producer-panel row-major representation and feed coarse
triangular/update kernels without rediscovering or branching over the KLU
column layout for every consumer.

A follow-up durable ragged-panel cache was implemented and rejected on
2026-06-28. This opt-in prototype used
`KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=ragged-panel`, built one cached
row-major internal block plus per-producer-row tail slices for each retained
EGraph supernode candidate, published each completed L column into that cache,
and consumed only large contiguous U dependency runs. This removed the previous
per-consumer L-column validation/scatter discovery but still kept the tail
updates as scalar ragged row scatters. It built and passed smoke tests, but the
top-five CKTSO-gap screen regressed from `1.6463s` to `1.7533s` geomean
(`build/kls_ragged_panel_control_gap5_t4_r1_ref3_timeout120.jsonl` versus
`build/kls_ragged_panel_gap5_t4_r1_ref3_timeout120.jsonl`). The mode did reach
substantial work: `ASIC_320ks` applied `2,216` runs over `339,862` rows and
`137,375,654` entries while refactor time moved from `0.08496s` to `0.09879s`;
`ASIC_320k` applied `2,294` runs over `265,289` rows and `130,486,393` entries
while refactor moved from `0.12431s` to `0.13253s`; and `ASIC_100ks` applied
`1,391` runs over `154,944` rows and `60,814,817` entries while refactor moved
from `0.05702s` to `0.06315s`. A repeat-20 `ASIC_100ks` amortization check
confirmed that this was not just first-use cache construction:
`build/kls_ragged_panel_control_asic100ks_t4_r1_ref20.json` measured
`0.04439s` refactor, while
`build/kls_ragged_panel_asic100ks_t4_r1_ref20.json` measured `0.05449s` after
`29,211` grouped runs, `3,253,824` rows, and `1,277,111,157` update entries.
The build used here had `build_has_cblas=false`, and earlier CBLAS-enabled
guards reported zero CBLAS counters on these same default slow paths, so this is
not evidence for a small-BLAS threshold problem. The source was reverted. The
remaining paper gap is more specific: KLS needs a production row/segment or
supernodal numeric engine that batches/shared-tail updates into coarse kernels,
not merely a prepacked scalar ragged-tail cache.

A retained cleanup now removes one piece of avoidable work from the existing
row-first cached-panel executor. `kls_row_first_partial_apply_supernode_run_cached`
previously reserved and zeroed fallback trailing scratch before trying the
cached CBLAS/native row-major panel path, even though successful cached-panel
runs allocate their own right-hand-side scratch and return before the fallback
buffer is used. The reserve/memset is now deferred until after the cached
paths decline the run, preserving all existing validation and fallback
semantics while avoiding dead scratch clearing on accepted cached panels.
`cmake --build build -j4`, `ctest --test-dir build --output-on-failure`,
`cmake --build build-cblas -j4`, `ctest --test-dir build-cblas
--output-on-failure`, and `git diff --check` passed. The same-session top-five
default CKTSO-gap focus stayed residual-clean and moved from `1.5319s` to
`1.5045s` geomean
(`build/kls_defer_cached_scratch_control_gap5_t4_r1_ref3_timeout120.jsonl`,
`build/kls_defer_cached_scratch_gap5_t4_r1_ref3_timeout120.jsonl`), although
that default path reported no row-panel-cache activity and should be treated
mostly as a regression guard. The focused `ASIC_320k` forced KLS-first factor
probe stayed residual-clean at `2.1189s` initial factor with `42,665`
row-supernode panel groups covering `779,443` rows
(`build/kls_asic320k_defer_cached_scratch_t4_factor.json`). The hard forced
METIS `pre2` factor probe still timed out at 120s with an empty JSON file
(`build/kls_pre2_metis_defer_cached_scratch_t4_factor_timeout120.json`), so
this is a retained executor hygiene improvement, not the missing CKTSO-scale
row/supernode numeric engine.

A follow-up GDB interrupt of the current forced-METIS `pre2` run still sampled
one active pipeline worker in `kls_row_first_supernodes_reset`, with peer
pipeline workers waiting on the condition variable and the parent joining
(`build/pre2_current_defer_scratch_gdb_interrupt.txt`). That made private
panel-cache rebuild storage reuse look attractive: keep the dynamic-pivot
rebuild semantics, but avoid freeing/reallocating `panel_id_by_row`, dense
panel, tail, and panel metadata arrays every time
`kls_row_first_supernode_panel_cache_build` is called. The prototype was
rejected for correctness. Reusing arrays, then adding back dense zeroing to
match the old `calloc`, and then adding full `panel_id_by_row` clearing still
made the sensitive forced KLS-first/no-fast `rajat29` probe fail setup as
singular with empty JSON files
(`build/kls_rajat29_reuse_panel_build_nofast_t4_factor_timeout90.json`,
`build/kls_rajat29_reuse_panel_build_zero_nofast_t4_factor_timeout90.json`,
`build/kls_rajat29_reuse_panel_build_fullclear_nofast_t4_factor_timeout90.json`).
The source was reverted. The result reinforces that the dynamic-pivot rebuild
path has hidden state-ordering requirements; fixing the `pre2` reset sample
needs a designed pivot-aware row/supernode representation, not a storage-reuse
shortcut around the current cache rebuild.

A more paper-aligned active-rank pivot prototype was also rejected. Instead of
resetting supernode metadata inside the pipeline worker after a dynamic pivot,
the worker aborted the active-rank phase before appending/exchanging the pivot
row and asked the existing restartable suffix logic to refactor that row
serially. This matches CKTSO's preference for restarting the dependent tail
after a pivot more directly than doing a full reset while peer workers wait,
and `ASIC_320k` remained residual-clean
(`build/kls_asic320k_active_pivot_abort_t4_factor.json`, `2.3335s` initial
factor, `1.509e-15` relative residual). The same change made the sensitive
forced KLS-first/no-fast `rajat29` probe fail setup as singular with an empty
JSON file
(`build/kls_rajat29_active_pivot_abort_nofast_t4_factor_timeout90.json`), so
the source was reverted. The existing active-rank pipeline restart path is not
yet semantically equivalent to the in-worker pivot path; a correct CKTSO-style
pivot-tail restart needs explicit state transfer for the active-rank
dependency/panel state, not just an earlier `KLS_ROW_FIRST_PIPELINE_FAIL_PIVOT`.

The retained follow-up is diagnostic rather than another small BLAS policy
change. KLS now reports active-rank pivot reset and panel-rebuild counters via
`kls_first_active_rank_pivot_reset_count`,
`kls_first_active_rank_pivot_reset_rows`,
`kls_first_active_rank_pivot_panel_rebuild_count`, and
`kls_first_active_rank_pivot_panel_rebuild_rows`, so future `pre2`-style
samples can distinguish full reset/rebuild cost from ordinary row-panel cache
activity. The current top-five CKTSO-gap rerun with this instrumentation stayed
residual-clean, reported zero active-rank pivot resets on all five matrices,
and measured `1.5481s` geomean in
`build/kls_active_rank_pivot_stats_gap5_t4_r1_ref3_timeout120.jsonl`. A focused
`ASIC_320k` KLS-first factor probe also reported zero active-rank reset and
panel-rebuild counters (`build/kls_asic320k_active_rank_pivot_stats_t4_factor.json`,
`2.4361s` initial factor, `1.509e-15` relative residual). The sensitive
`rajat29` no-fast pivot probe failed setup as singular in this diagnostic build,
so it was not used as a pass criterion.

The same rerun closes the small-BLAS guard hypothesis for the current focus
set. The CBLAS-capable binary already requires the runtime
`KLS_ENABLE_CBLAS_SUPERNODE=1` gate plus 512-scale vector/panel or minimum-work
tests before calling BLAS. With `OPENBLAS_NUM_THREADS=1`, the CBLAS-capable
top-five control measured `1.4519s` geomean with
`KLS_ENABLE_CBLAS_SUPERNODE=0` and `1.5900s` with
`KLS_ENABLE_CBLAS_SUPERNODE=1`
(`build-cblas/kls_cblas_guard_active_rank_off_gap5_t4_r1_ref3_timeout120.jsonl`
and `build-cblas/kls_cblas_guard_active_rank_on_gap5_t4_r1_ref3_timeout120.jsonl`).
Both runs had `build_has_cblas=true`, but every matrix reported zero CBLAS
update runs, rows, and entries. Adding a stricter "large cases only" BLAS guard
would therefore not affect the observed CKTSO-gap path; the active issue remains
the native row-first/pivot-aware numeric executor shape.

The next profiled hot-loop cleanup was also rejected as too small and mixed to
keep. `perf record` was unavailable on this host because `perf_event_paranoid`
is set to `4`, so Callgrind was used on `ASIC_100ks` instead. That mixed cold
factor/refactor profile put `kls_egraph_refactor_btf_unscaled_column` at
`37.58%` of sampled instructions, behind the cold KLU kernel at `41.50%`. A
prototype hoisted the cached 32-bit `L` row-index table lookup out of each
EGraph scatter helper call. It built and passed the normal smoke suite, but a
same-commit top-five comparison was noise-level and matrix-mixed: the clean
`c5b3f8a` control measured `1.44549s` geomean
(`build/kls_baseline_c5b3f8a_gap5_t4_r1_ref3_timeout120.jsonl`) and the
prototype measured `1.44269s`
(`build/kls_cached_i32_gap5_t4_r1_ref3_timeout120.jsonl`), while
`ASIC_320ks` and `rajat03` regressed. Focused repeats were also mixed
(`ASIC_100ks` baseline refactor `0.04017s` versus prototype `0.05558s` and
`0.03743s`; `ASIC_320k` baseline `0.10861s` versus prototype `0.10458s` and
`0.10562s`). The source was reverted. This keeps the next work item at the
algorithm level: a production CKTSO/SubtreeLU-style row/supernode numeric
executor, not a scalar EGraph scatter cleanup.

The auto-ordering policy now avoids one redundant analysis pass on the large
full/nearly-full diagonal METIS-start class. Fresh controls showed that
`--ordering auto` already selected the same METIS/BTF symbolic as forced
METIS on `ASIC_320ks`, but spent about twice the analysis time because the
auto METIS-start branch also tried a guarded no-BTF symbolic and then kept the
BTF result. Forced no-BTF on the same matrix had similar analysis cost but
worse repeated-refactor work (`0.10059s` versus `0.08390s` refactor average in
`build/asic320ks_metis_nobtf_ref3.json` and
`build/asic320ks_metis_ref3.json`). KLS now skips that no-BTF retry when the
large diagonal METIS-start structural predicate fired, preserving the retained
BTF symbolic. The focused analyze-only check moved `ASIC_320ks` auto analysis
from `2.01145s` to `1.01080s`
(`build/asic320ks_auto_analyze.json`,
`build/asic320ks_auto_skip_nobtf_analyze.json`) with the same METIS/BTF
selection; matching post-change controls measured `1.00920s` on `ASIC_320k`
and `0.35153s` on `ASIC_100ks`
(`build/asic320k_auto_skip_nobtf_analyze.json`,
`build/asic100ks_auto_skip_nobtf_analyze.json`). This improves production
auto-analysis overhead for the current top-gap ASIC shape, but it does not
change the main refactor diagnosis: the repeated numeric path still needs the
paper-level row/supernode executor.

The next post-change selector checks rejected three tempting but narrower
shortcuts. First, forcing the existing row-refactor engine on the top-five
CKTSO-gap focus set measured `4.26764s` geomean versus `1.58809s` for the
current auto control
(`build/kls_forced_row_after_skip_nobtf_gap5_t4_r1_ref3_timeout120.jsonl`,
`build/kls_auto_skip_nobtf_gap5_t4_r1_ref3_timeout120.jsonl`). The forced row
path did run in parallel and built row work, but every focus row regressed,
confirming that the row lower-bound/model gate is protecting the scalar KLU
storage path rather than hiding a ready win. Second, enabling the column EGraph
successor ready queue measured `2.52798s` geomean on the same focus set
(`build/kls_egraph_ready_queue_gap5_t4_r1_ref3_timeout120.jsonl`). It recorded
ready-queue columns on the EGraph rows, but `ASIC_320ks`, `ASIC_320k`, and
`ASIC_100ks` regressed sharply because the queue still feeds the same
per-column scalar scatter executor. Third, a small static-pivot policy check on
`gemat11`/`gemat12` again showed why disabling static pivoting is not a broad
fix. `--no-static-pivoting` improved `gemat12` (`0.05018s` versus `0.05661s`)
but regressed `gemat11` (`0.05351s` versus `0.03990s`), worsening the pair
geomean from `0.04753s` to `0.05182s`
(`build/kls_gemat_static_current_t4_r1_ref20_timeout120.jsonl`,
`build/kls_gemat_nostatic_t4_r1_ref20_timeout120.jsonl`). These checks leave
the paper-aligned target unchanged: KLS needs a production row-major/supernode
numeric representation and coarse update kernels, not just different selector
policy around the current scalar executors.

A narrower cached-panel cleanup was retained after the next current-source
rerun. The previous per-panel adaptive pruning was rejected because it added
bookkeeping and regressed the cached opt-in path, but the current top-ten probe
showed a simpler no-work case: `ASIC_320k`, `rajat24`, `transient`, and
`rajat28` built retained EGraph panels and paid cached-probe attempts while
applying zero cached updates. KLS now records a pass-level cached-only guard:
when `KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=cached` completes a valid numeric pass
with cached-probe attempts but zero work-gate-accepted updates, later passes
with the same retained panel cache skip the cached probe and report
`refactor_supernode_cached_probe_disabled` plus
`refactor_supernode_cached_probe_disable_count`. The focused cached top-ten
artifact after the change
(`build/kls_cached_disable_guard_gap10_t4_r1_ref3_timeout120.jsonl`) shows the
guard firing once on the zero-work rows and reducing their later-pass cached
attempts to zero (`ASIC_320k`, `rajat24`, `transient`, and `rajat28`), while
productive cached rows such as `ASIC_320ks`, `ASIC_100ks`, `G2_circuit`,
`onetone1`, and `onetone2` remain enabled. The same artifact is not evidence
that cached EGraph panels should become the default: it measured `4.41862s`
geomean versus the same-session default top-ten `4.33669s`, and the cached path
still loses the top-ten common set to CKTSO by about `2.47x`. This is retained
only as opt-in row/supernode executor cleanup; the missing paper-level item
remains a coarse row-major/supernode numeric engine rather than another cached
selector rule.

The full EGraph supernode experiment now has the same kind of pass-level
amortization guard for its low-work case. A current-source top-five rerun first
confirmed that cached-only panels should remain opt-in:
`build/kls_current_default_rerun_gap5_t4_r1_ref3_timeout120.jsonl` measured
`1.47590s` geomean, while
`build/kls_current_cached_rerun_gap5_t4_r1_ref3_timeout120.jsonl` measured
`1.52588s`. The full `KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1` bridge measured
`1.47197s` before the guard in
`build/kls_current_egraph_supernode_rerun_gap5_t4_r1_ref3_timeout120.jsonl`,
but `ASIC_320k` did only two tiny supernode updates per pass
(`10` entries) while still probing about `1,600` cached candidates. KLS now
records `refactor_supernode_update_disabled` and
`refactor_supernode_update_disable_count`; after one valid full-supernode pass
with cached-probe attempts but fewer than `512` accepted update entries, later
passes skip the full supernode probe for the retained panel cache. The focused
guard artifact
`build/kls_supernode_lowwork_guard_gap5_t4_r1_ref3_timeout120.jsonl` shows
`ASIC_320k` reducing cumulative cached attempts from `6,388` to `1,597` and
cumulative supernode work from `40` to `10` entries, with the disable counter
set to `1`. Productive rows such as `ASIC_320ks`, `ASIC_100ks`, and `rajat03`
did not disable. A second guarded rerun still remained noisy and mixed
(`1.54805s` geomean in
`build/kls_supernode_lowwork_guard_gap5_rerun2_t4_r1_ref3_timeout120.jsonl`),
so the full supernode bridge remains experimental rather than a new default.
This is retained only to avoid repeated no-amortization work inside the
paper-aligned supernode prototype; the broad CKTSO gap still requires a
production row-major/supernode numeric executor.

KLS now has an opt-in structural cache for the paper-style row-major
U-supernode producer object, separate from the stricter common-trailing L-panel
cache. With `KLS_ENABLE_REFACTOR_U_SUPERNODE_PATTERN=1`, schedule construction
retains each detected adjacent U-pattern supernode's start, local start, width,
and exact right-side successor-column list, and reports
`refactor_u_supernode_pattern_*` counters. The current top-five diagnostic with
consumer stats enabled
(`build/kls_u_supernode_pattern_gap5_t4_r1_ref3_timeout120.jsonl`) measured a
`1.50934s` geomean and confirmed that the new structural object covers the
producer side on the slow ASIC rows: `ASIC_320ks` retained `4,705` patterns,
`13,346` rows, `70,267` right entries, and `158,601` internal entries;
`ASIC_320k` retained `4,819`/`13,557`/`69,170`/`218,207`; and `ASIC_100ks`
retained `8,481`/`24,581`/`153,053`/`108,725`. The same rows still had zero
default retained common-tail panels, while the consumer diagnostic saw
`73551`, `72446`, and `160299` dependency runs respectively. This narrows the
paper gap: KLS now preserves the producer pattern that CKTSO/SubtreeLU would
consume, but it still lacks the production numeric executor that updates
consumer rows from that row-major U object.

KLS now also has an opt-in numeric cache for the same retained U-supernode
producer object. With `KLS_ENABLE_REFACTOR_U_SUPERNODE_VALUES=1`, the EGraph
refactor allocates row-major dense and right-side U value buffers for the
retained pattern and records U values as columns publish. The focused top-five
diagnostic
(`build/kls_u_supernode_values_gap5_t4_r1_ref3_timeout120.jsonl`) stayed
residual-clean and populated substantial numeric storage: across the five rows
it allocated `1,038,387` dense slots and `2,607,110` right-side slots, with
`2,182,336` dense writes and `10,428,440` right-side writes over the repeated
refactors. The slow ASIC rows had no external CBLAS calls in this run, so this
is not a small-BLAS-call artifact. The same opt-in run measured a worse
`1.70716s` geomean, while the gate-off current-source control
(`build/kls_u_supernode_values_default_gap5_t4_r1_ref3_timeout120.jsonl`)
measured `1.41058s` and kept all U-supernode value counters at zero. This
retains the paper-style numeric staging object but confirms the missing
performance piece is still a consumer/executor that uses those cached values
without paying an extra scalar recording pass.

KLS now has an opt-in ragged L-panel executor for the same retained
U-supernode producer ranges. With
`KLS_ENABLE_REFACTOR_U_SUPERNODE_RAGGED_L=1`, schedule construction retains the
U-supernode structure, the refactor path allocates a ragged L-panel object,
producer columns publish fresh L values into that object as each numeric pass
finishes them, and later contiguous dependency runs can consume the cached
internal panel plus per-producer trailing L rows. This is the direct
paper-aligned extension beyond the older common-trailing cached panel: it can
represent nearly all retained U-supernode producers on the ASIC focus rows
instead of requiring identical trailing L row lists. The focused top-five
artifact
(`build/kls_u_supernode_ragged_l_gap5_t4_r1_ref3_timeout120.jsonl`) stayed
residual-clean and built `18,218` valid ragged L panels over `747,694` dense
entries and `2,610,433` trailing entries. The executor did real work:
`637` last-pass runs covered `89,907` producer rows and `18,044,382` update
entries, with `2,548` cumulative runs over the repeated refactors. It still
lost the focused set, measuring a `1.65257s` geomean versus the same-source
gate-off control
(`build/kls_u_supernode_ragged_l_default_gap5_t4_r1_ref3_timeout120.jsonl`) at
`1.33400s`, and `ASIC_100ks` regressed from `5.23s` to `11.34s`. The gate-off
control kept all U-supernode ragged-L counters at zero. This validates the
broader producer/consumer executor mechanically, but shows that simply caching
ragged producer L rows is not enough; the next paper gap is coarser batching or
a policy that only pays the publish/cache cost when producer rows are reused
enough to amortize it.

The ragged L-panel executor now has the same clean-pass amortization guard as
the earlier cached-panel probes. During the first valid numeric pass it marks
only panels that actually feed a dependency run, then disables untouched panels
for later repeated refactors. The focused prune artifact
(`build/kls_u_supernode_ragged_l_prune_gap5_t4_r1_ref3_timeout120.jsonl`)
improved the opt-in geomean from `1.65257s` to `1.44225s` and reduced the
`ASIC_100ks` regression from `11.34s` to `5.70s`. The counters show why:
`18,215` of `18,218` valid ragged panels were pruned after one clean pass,
covering `557,004` dense entries and all `2,610,433` trailing entries that
were not used later. The executor still did the same real work on the remaining
panels (`637` last-pass runs, `89,907` producer rows, and `18,044,382` update
entries), and residuals remained clean. The same-source gate-off control
(`build/kls_u_supernode_ragged_l_prune_default_gap5_t4_r1_ref3_timeout120.jsonl`)
measured `1.36151s` with all ragged-L counters at zero, so the prune is
retained as opt-in executor cleanup rather than promoted to default policy. It
narrows the missing piece further: KLS can now identify and discard unused
ragged panels, but the remaining used panel work still needs a coarser kernel or
an auto policy that predicts when the retained panel will beat scalar EGraph.

A fresh same-binary CBLAS guard check on the current source again rejects the
"small BLAS calls" hypothesis. The CBLAS-capable build was rebuilt, then the
top-five CKTSO-gap focus set was run with `OPENBLAS_NUM_THREADS=1` and the
runtime CBLAS gate off and on. The gate-off artifact
(`build-cblas/kls_cblas_guard_latest2_off_gap5_t4_r1_ref3_timeout120.jsonl`)
measured a `1.40038s` geomean; the gate-on artifact
(`build-cblas/kls_cblas_guard_latest2_on_gap5_t4_r1_ref3_timeout120.jsonl`)
measured `1.47634s`. Every row in both files reported
`build_has_cblas=true` but zero external CBLAS update runs, rows, and entries,
and zero compact-supernode GEMV/TRSV or blocked-panel counters. The active
paths already require the runtime CBLAS gate plus 512-scale vector/panel tests
and multi-million-operation work thresholds before calling BLAS. There is
therefore no useful additional "only use BLAS for large cases" patch for these
slow rows; the remaining gap is still the paper-level row-major/supernode
executor that creates reusable coarse work, not BLAS call granularity.

A direct ragged-panel tail-union accumulator was prototyped and rejected. The
experiment kept the opt-in retained U-supernode ragged L-panel cache, then
added `KLS_ENABLE_REFACTOR_U_SUPERNODE_L_UNION=1` to precompute each panel's
unique trailing row set and map per-row ragged tail entries into that union.
At run time, accepted ragged updates accumulated repeated trailing rows in
worker scratch and scattered each touched union row once. This tested whether
the remaining used ragged panels were losing mainly because they wrote the same
trailing `x` rows repeatedly. Correctness passed (`ctest --test-dir build
--output-on-failure` and opt-in `kls_smoke`), but the focused benchmarks did
not support retaining the code. On the current five-row focus
(`ASIC_320k`, `ASIC_320ks`, `onetone2`, `ASIC_100ks`, `G2_circuit`), ragged
without the union measured `8.94150s` geomean in
`build/kls_u_supernode_l_union_ragged_control_gap5_t4_r1_ref3_timeout120.jsonl`,
while union measured `9.24227s` in
`build/kls_u_supernode_l_union_gap5_t4_r1_ref3_timeout120.jsonl`; the default
control was still better at `8.11402s` in
`build/kls_union_default_control_gap5_current_t4_r1_ref3_timeout120.jsonl`.
On the older ASIC/gemat/rajat focus, union measured `1.46720s` in
`build/kls_u_supernode_l_union_oldgap5_t4_r1_ref3_timeout120.jsonl`, worse
than the saved ragged-prune `1.44225s` and default `1.36151s` controls. The
prototype was reverted. This narrows the missing paper executor again: the
problem is not just duplicate trailing-row writes inside retained ragged
producer panels; KLS needs a coarser multi-consumer/panel task or row-major
numeric object that avoids the current per-current-column execution shape.

KLS now records the missing multi-consumer reuse shape explicitly. The existing
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_STATS=1` diagnostic still scans actual
U-dependency streams for SubtreeLU-style contiguous consumer runs, but it now
also groups those runs by producer supernode start. New stats report how many
producer panels have any consumer run, how many are reused by at least two
consumer runs, the reused run/row totals, max runs and rows per producer
panel, and the L/internal work covered by reused panels. The older
ASIC/gemat/rajat focus in
`build/kls_consumer_reuse_stats_gap5_t4_r1_ref3_timeout120.jsonl` shows that
the hard ASIC rows are not isolated single-consumer opportunities:
`ASIC_320ks` has `73,413` of `73,551` consumer runs and `896,979` of
`897,255` consumer rows inside reused producer panels; `ASIC_320k` has
`72,321`/`72,446` runs and `938,917`/`939,167` rows; `ASIC_100ks` has
`160,281`/`160,299` runs and `1,148,063`/`1,148,099` rows. The current
five-row focus in
`build/kls_consumer_reuse_stats_current_gap5_t4_r1_ref3_timeout120.jsonl`
confirms this is not ASIC-only: `onetone2` has `431,426` of `431,440`
consumer rows in reused panels, and `G2_circuit` has all `5,569,687`
consumer rows in reused panels. This strengthens the next implementation
target: a producer-centered multi-consumer task can potentially amortize the
retained row-major panel publication across hundreds of consumer runs per
panel; further single-current-column tail or scalar scatter variants are the
wrong granularity.

KLS now has an opt-in retained consumer-run plan for that producer-centered
task. With `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN=1`, schedule building
keeps a producer-panel pointer array plus compact `(current column, producer
dependency, run rows)` arrays for the same contiguous consumer runs reported by
the reuse diagnostic. This does not execute a new numeric kernel yet; it closes
the scheduling/data-ownership gap needed before a CKTSO/SubtreeLU-style
producer task can publish one panel and feed many consumers. A JSON smoke on
`rajat03` parsed cleanly and retained `6,598` runs and `41,466` consumer rows
across `257` reused panels. The current top-five gap focus
(`build/kls_consumer_plan_current_gap5_t4_r1_ref3_timeout120.jsonl`) retained
large reusable plans on every row: `ASIC_320k` has `72,446` runs and
`939,167` rows in a `4.31 MB` plan, `ASIC_320ks` has `73,551` runs and
`897,255` rows in `4.34 MB`, `onetone2` has `50,274` runs and `431,440`
rows in `1.50 MB`, `ASIC_100ks` has `160,299` runs and `1,148,099` rows in
`4.64 MB`, and `G2_circuit` has `424,370` runs and `5,569,687` rows in
`11.39 MB`. The same run again reported zero external CBLAS, compact GEMV, and
compact TRSV counters on all five rows, so the small-BLAS-call hypothesis is
not active on these losses; the next useful implementation step is a numeric
executor that consumes this retained producer-panel plan.

The retained plan now also stores each run's producer-panel start column and
producer-panel local offset. The plan executor and producer-panel cache pruning
consume those retained fields instead of rediscovering the producer panel from
the dependency column. This is still scaffolding rather than a speed path: it
turns the plan into a producer-addressable run table needed by a future
producer-centered task, but it does not create the persistent target
accumulators that would let that task update many consumers at once.

The retained consumer plan can now drive an opt-in cached-panel executor with
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_EXEC=1`. The executor builds the
same retained plan, adds a current-column lookup index, uses planned
`(current, producer, run_rows)` entries to select cached producer-panel updates,
and records attempt/hit/applied counters. The producer-panel cache pruning step
also consumes the retained plan directly instead of rediscovering every
contiguous run from U. Correctness remained clean: the `rajat03` smoke applied
`44` planned runs and `1,310` rows on the first refactor, then disabled the
executor after the low-coverage clean pass; the final relative residual was
about `3.9e-15`.

The focused top-five result rejects this cached-panel executor as the missing
CKTSO mechanism. The unguarded plan executor was correct but slower
(`8.54s` geomean in
`build/kls_consumer_plan_exec2_current_gap5_t4_r1_ref3_timeout120.jsonl`).
After adding the structural coverage guard, the same focus set measured
`7.67s` in
`build/kls_consumer_plan_exec_guard_current_gap5_t4_r1_ref3_timeout120.jsonl`,
slightly slower than the same-binary no-env control at `7.62s` in
`build/kls_consumer_plan_exec_guard_default_current_gap5_t4_r1_ref3_timeout120.jsonl`.
The counters explain why the path is insufficient: `ASIC_320ks` applied only
`48,841` of `897,255` retained plan rows, `onetone2` applied `26,108` of
`431,440`, `ASIC_100ks` applied `39,436` of `1,148,099`, and `G2_circuit`
applied `65,265` of `5,569,687`; `ASIC_320k` found hits but no planned run
passed the cached-panel executor's shape/work gates. This is a direct
paper-gap result, not a small tuning result: KLS can now retain and execute
planned producer-panel runs, but the existing completed-panel cache covers only
a tiny fraction of the reusable work. Closing the CKTSO gap requires a broader
row-major or persistent-consumer accumulator that can consume the retained plan
beyond the current equal-trailing cached-panel shape.

A fresh same-source CBLAS-capable rerun confirms that "guard BLAS for only large
cases" is already handled and is not the active loss mode. The current code's
CBLAS call sites require both the runtime `KLS_ENABLE_CBLAS_SUPERNODE` gate and
large shape/work thresholds before calling BLAS. With `OPENBLAS_NUM_THREADS=1`,
the top-five focus measured `1.47818s` geomean with the runtime CBLAS gate off
in
`build-cblas/kls_cblas_guard_rerun_current3_off_gap5_t4_r1_ref3_timeout120.jsonl`
and `1.43948s` with it on in
`build-cblas/kls_cblas_guard_rerun_current3_on_gap5_t4_r1_ref3_timeout120.jsonl`.
Every matrix in both runs reported `build_has_cblas=true`, but all external
CBLAS update, compact-supernode GEMV, and compact-supernode TRSV counters were
zero. The small timing movement is therefore ordinary run noise or indirect
layout noise, not BLAS work. The slow rows stayed on the EGraph refactor path,
and the row-refactor model again rejected automatic handoff through the
lower-bound work gate. Further BLAS-size threshold changes would not close the
CKTSO gap on these cases; the missing mechanism is still the paper-level
producer/consumer row-major numeric executor.

The June 29, 2026 follow-up rerun reached the same conclusion on the current
source after the object-workspace cleanup. A freshly rebuilt CBLAS-capable tree
was run on the first five CKTSO-gap matrices with `OPENBLAS_NUM_THREADS=1`,
four KLS threads, one initial factor, and three refactors. With the runtime
CBLAS gate disabled,
`build-cblas/kls_cblas_guard_current_off_gap5_t4_r1_ref3_timeout120.jsonl`
measured a `2.67328s` SPICE-cycle geomean; with
`KLS_ENABLE_CBLAS_SUPERNODE=1`,
`build-cblas/kls_cblas_guard_current_on_gap5_t4_r1_ref3_timeout120.jsonl`
measured `4.14533s`. Both artifacts reported `build_has_cblas=true`, but every
row still had zero external CBLAS update runs, rows, and entries. The default
non-CBLAS build in
`build/kls_default_guard_current_gap5_t4_r1_ref3_timeout120.jsonl` reported
`build_has_cblas=false` and zero CBLAS/compact-panel counters. This means the
requested "use BLAS only for large cases" policy is already true in the source:
CBLAS requires the build option, the runtime gate, 512-scale shape checks, and
multi-million-work gates, and the current hard-gap rows do not enter those call
sites at all.

After the low-work fast-factor guard, the same current-source check was repeated
with the CBLAS-capable binary on `LeGresley_87936`, `trans4`, `ACTIVSg2000`,
`onetone2`, and `rajat28`. Two order-reversed passes with
`OPENBLAS_NUM_THREADS=1` and `KLS_ENABLE_CBLAS_SUPERNODE=0/1` kept every row on
`last_refactor_path=egraph` and reported zero CBLAS update runs, rows, and
entries in both gate settings. The runtime geomean moved from `2.78635s` with
the gate off to `2.39433s` with the gate on, but no BLAS call was taken, so that
movement is ordinary run/cache noise. There is no useful source change in
another small-BLAS threshold guard; the active gap remains outside the optional
CBLAS kernels.

The retained-plan ragged-L pruning/prefilter direction was also tested and
rejected. The combined ragged-L plus retained consumer-plan run measured
`1.51766s` geomean on the top-five gap focus in
`build/kls_ragged_l_plan_probe_gap5_t4_r1_ref3_timeout120.jsonl`, already
slower than the same-binary default control at `1.39148s` in
`build/kls_default_control_probe_gap5_t4_r1_ref3_timeout120.jsonl`. Moving the
retained plan into a post-build prune step measured `1.54823s` in
`build/kls_ragged_l_plan_preprune_gap5_t4_r1_ref3_timeout120.jsonl`, and
skipping the second prune scan only improved that to `1.51137s` in
`build/kls_ragged_l_plan_preprune2_gap5_t4_r1_ref3_timeout120.jsonl`.
Filtering panels before metadata construction was much worse: the focused run
measured `3.69737s` in
`build/kls_ragged_l_plan_prefilter_gap5_t4_r1_ref3_timeout120.jsonl`, and the
counters showed that ragged update rows dropped to zero on the ASIC rows after
the filter. This rejects static plan pruning as the clear missing CKTSO idea:
it can remove some publication scans, but it does not create the reusable
row-major/persistent accumulator that the papers rely on to feed many consumers
from one producer panel.

A focused rerun after adding retained producer-panel starts/offsets confirmed
that this metadata alone is not the missing mechanism. With
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_EXEC=1`, the top-five focus
measured `3.95505s` geomean in
`build/kls_consumer_plan_panel_fields_gap5_t4_r1_ref3_timeout120.jsonl`,
versus the same-binary no-env control at `3.18461s` in
`build/kls_consumer_plan_panel_fields_default_gap5_t4_r1_ref3_timeout120.jsonl`.
Coverage remained the same narrow shape as before: `ASIC_320ks` retained
`897,255` plan rows but applied `48,841`, `ASIC_320k` retained `939,167` and
applied none, `rajat03` retained `41,466` and applied `1,310`, and
`ASIC_100ks` retained `1,148,099` and applied `39,436`; the executor disabled
itself after the low-coverage clean pass on those rows. The retained
producer-panel fields are useful for the next producer-task implementation, but
they do not change the current conclusion: KLS needs a coarser row-major or
persistent-consumer accumulator, not more cached-panel lookup metadata.

The independent row multi-producer executor then had its per-accepted-batch
allocation path removed. Both contiguous and fragmented independent-row
multi-producer kernels now use the existing worker index workspace for run
groups, row offsets, per-run offsets/suffixes/lengths/ids, and U offsets; the
worker double workspace for multipliers, pivots, and U scratch; and the worker
object/byte workspace for run panel descriptors. This matches the storage
direction implied by CKTSO/SubtreeLU better than repeated `malloc`/`calloc`
setup around the row-panel executor, while preserving the same scalar arithmetic
and pivot checks.

The focused forced-row top-five rerun is coverage-limited, not a gap-closing
result. The same-session control
`build/kls_worker_scratch_forced_control_gap5_t4_r1_ref3_timeout120.jsonl`
measured `8.13943s` geomean, while
`build/kls_worker_scratch_multiproducer_gap5_t4_r1_ref3_timeout120.jsonl` with
`KLS_ENABLE_MULTI_PRODUCER_SUPERNODE=1` measured `5.69134s` geomean. However,
all five matrices still reported zero
`row_refactor_compact_supernode_batch_count` and zero candidate counts; the
observed wins therefore cannot be credited to accepted independent
multi-producer batches. The synthetic smoke fixtures remain the correctness
coverage for this storage cleanup. The paper-level conclusion is unchanged:
the current hard rows need the broader native row/segment producer-consumer
numeric engine, not another small gate around the independent-batch scaffold.

The same worker-object storage was then separated from the byte workspace and
applied to the dense-consumer fragmented producer executor's run-panel
descriptors. This matters because the dense path can need both live panel
descriptors and byte scratch for fallback target maps; using the byte workspace
for both would be unsafe. The worker now owns a dedicated object buffer that is
freed with the other row-refactor worker scratch buffers, and the dense
fragmented path no longer `calloc`/`free`s descriptor arrays per accepted batch.

This is still a storage cleanup, not the missing paper kernel. With current
selector policy, unset `KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR` keeps the scalar
dense-group path and therefore accepts zero dense-fragmented batches on the old
top-five dense-fragmented focus. Forcing the current compact selector with
`KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=auto` and
`KLS_ENABLE_MULTI_PRODUCER_SUPERNODE=1` exercises the same batch shape, but the
overall scaffold remains slow (`62.4818s` geomean in
`build/kls_dense_fragmented_object_workspace_auto_top5_t4_r1_ref2_timeout120.jsonl`).
An isolated same-mode `ASIC_320k` comparison against the previous commit shows
the descriptor storage change itself is not the regression source:
`515f202` measured `72.3339s` in
`build/kls_prev515_dense_fragmented_auto_asic320k_t4_r1_ref2_timeout120.jsonl`,
while the object-workspace build measured `69.4805s` in
`build/kls_object_workspace_dense_fragmented_auto_asic320k_t4_r1_ref2_timeout120.jsonl`;
both runs had identical dense-fragmented coverage (`70` compact panels,
`352` accepted batches, `246,898` batched dependency rows). The remaining gap
is therefore still the native row/segment numeric representation and production
blocked update, not descriptor allocation.

The next broad-suite timeout fix targets the policy path around
`HTC_336_9129`. The refreshed medium comparison had scored this row as a
120-second timeout while CKTSO solved it in `8.08898s`; focused probes showed
that explicit AMD/no-BTF/unscaled factorization completed cleanly, but auto
policy could spend the whole cap in static-matching or METIS-promotion retry
territory before settling. KLS now recognizes a narrow low-work direct-AMD
class before those expensive numeric probes: auto ordering, AMD-selected,
single no-BTF block, at least 200k rows, symbolic work below `2e8`, average
stored degree at most `3.40`, and a 40-60% structural diagonal fraction. For
that class, KLS skips pre-static row matching, the SPRAL Hungarian numeric
trial, and numeric METIS promotion, letting the already-low-work AMD/no-BTF
factor run directly.

This is a failure-removal policy, not the missing row/supernode kernel. The
targeted default-auto suite row now completes:
`build/kls_lowwork_direct_amd_htc9129_suite_t4_r1_ref3_timeout120.jsonl`
reports `HTC_336_9129` at `16.2184s` with AMD, no BTF, scale `-1`, no static
matching, and a clean EGraph refactor path; the refreshed medium run
`build/kls_lowwork_direct_amd_paper_medium_t4_r1_ref3_timeout120.jsonl`
measured the same row at `7.64591s`, slightly faster than the saved CKTSO
artifact for that matrix. The denser `HTC_336_4438` row stayed on its existing
AMD/BTF/static-SPRAL shape and did not match the new sparse low-work guard.
The smoke suite still passes. The full medium run improved the failure count
from four rows to three (`mac_econ_fwd500`, `ss1`, and singular
`bips07_1998` remain), but KLS is still not generally on par with CKTSO:
with failures scored at 1000 seconds, the refreshed artifact measured
`0.63138s` geomean versus CKTSO's saved `0.33559s`. The next gap is therefore
still the remaining timeout rows and the broad slow repeated-refactor cases,
not this now-removed HTC auto-policy failure.

The next failure-removal change addresses `ss1`, a large, very-low-degree,
full-diagonal matrix. Earlier probes showed that default auto selected METIS
with BTF, creating `95,977` BTF blocks and timing out under the full
`repeat=1, refactor_repeat=3, timeout=120` benchmark loop. The same METIS
ordering without BTF completed under the wall-time cap and made repeated
refactors much cheaper (`~11-12s` instead of `~57-60s` in focused probes).
KLS now starts the existing structural class
`is_large_very_low_degree_full_diagonal_pattern` with `use_btf=0` directly,
instead of running a score-based no-BTF retry. This is intentionally structural:
the paper-medium set has three matching matrices (`G2_circuit`, `mc2depi`,
and `ss1`), and the two collateral rows already behaved as single METIS
blocks even when BTF was requested.

Focused verification on the current source:
`build/kls_lowdegree_nobtf_G2_circuit_current_t4_r1_ref3_timeout120.jsonl`
completed `G2_circuit` with METIS/no-BTF at `44.3114s`,
`build/kls_lowdegree_nobtf_mc2depi_current_t4_r1_ref3_timeout120.jsonl`
completed `mc2depi` with METIS/no-BTF at `191.031s`, and
`build/kls_lowdegree_nobtf_ss1_current_t4_r1_ref3_timeout120.jsonl`
completed `ss1` with METIS/no-BTF at `1125.86s`. The `ss1` row is therefore
no longer a timeout, but it is still far slower than CKTSO's saved
`348.787s` SPICE-cycle estimate. A projection artifact that replaces only
the three affected rows in the previous medium run,
`build/kls_lowdegree_nobtf_medium_projection_t4_r1_ref3_timeout120.jsonl`,
scores one remaining missing candidate row (`mac_econ_fwd500`) and measures
KLS at `0.58340s` geomean versus CKTSO at `0.30764s` on the common scored
set. This keeps the main conclusion unchanged: the BLAS-size guard hypothesis
is not active on these rows (all focused rows reported zero CBLAS calls), and
the remaining large gap is repeated-refactor numeric work rather than small
external BLAS calls.

The remaining medium timeout row, `mac_econ_fwd500`, was then traced on the
current source. A temporary stage trace showed that default auto reached
`maybe_select_pre_static_row_match` after METIS/BTF analysis and did not leave
that large SPRAL/static pre-factor candidate before a 100-second cap. Disabling
static pivoting confirmed that the base METIS/BTF path still completed its
first factor, but then spent about 43 seconds in the auto pivot-tolerance
trial before entering the repeated factor. The retained fix is therefore a
structural policy correction for the existing
`is_large_sparse_diagonal_low_degree_pattern` class: skip the large pre-static
match, skip the post-factor SPRAL numeric trial, skip the redundant auto
pivot-tolerance trial, and choose `1e-5` as the initial pivot tolerance for
that class. Explicit probes rejected a fast-factor exception for the same
dominant-BTF shape because it timed out under a 140-second cap, and rejected
lower tolerance choices as an accuracy tradeoff (`1e-6` was faster but had a
larger relative residual; `3e-6` was worse still).

The focused standard medium-row rerun now completes:
`build/kls_sparse_diag_trials_default_mac_r1_ref3_timeout120.jsonl` reports
`mac_econ_fwd500` with METIS/BTF, scale `1`, selected pivot tolerance `1e-5`,
no selected static/SPRAL match, analysis `11.8659s`, initial factor `33.4531s`,
repeat factor `28.5324s`, refactor `8.8153s`, and SPICE-cycle estimate
`930.249s`; the wall-clock wrapper completed in `101.37s` under the 120-second
cap. The residual was valid but not CKTSO-class (`relative_residual_l2`
`2.07e-7`), so this is not a license to lower pivot tolerance further without
more accuracy evidence. A projection replacing this row in the previous
medium projection,
`build/kls_sparse_diag_trials_medium_projection_t4_r1_ref3_timeout120.jsonl`,
has zero scored missing KLS rows versus the saved CKTSO medium artifact, but
still measures KLS at `0.58294s` geomean versus CKTSO at `0.30764s`
(`1.895x` candidate/reference). The next paper-level gap is therefore no
longer timeout coverage on the medium set; it is the broad repeated-factor and
refactor numeric gap on rows such as `ACTIVSg2000`, `LeGresley_87936`,
`trans4`, `rajat24`, `HTC_336_4438`, and the ASIC family.

The next retained improvement targets a low-work no-pivot repeated-factor
policy failure. In the projected medium comparison, `ACTIVSg2000`,
`ACTIVSg10K`, and `ACTIVSg70K` all used the checked `kls_fast_refactor` factor
path even though their numeric work was only about `2.0e6`, `4.4e6`, and
`9.2e6` flops and there were no off-diagonal pivots. Focused default versus
`--no-fast-factor` probes showed that the checked fast-factor scaffold was
more expensive than rebuilding the KLU numeric object and keeping the EGraph
refactor schedule for these shapes. The retained gate therefore skips
fast-factor repair for two narrow structural cases: single-block/no-BTF,
no-static, no-pivot matrices with `1e6-1e7` factor flops and at least 10k
rows, and small dominant-BTF/no-pivot matrices with 128-512 BTF blocks,
`1e6-3e6` factor flops, 3k-10k rows, and 90-97% largest-block coverage. Guard
probes kept `powersim`, `OPF_10000`, `G2_circuit`, `rajat03`, `hcircuit`, and
`memplus` outside the new skip or on their previous fast path.

Focused artifacts on the retained source report `ACTIVSg2000` at `0.658s`,
`ACTIVSg10K` at `1.581s`, and `ACTIVSg70K` at `1.146s`
(`build/kls_lowwork_skip_fast_activsg2000_t4_r1_ref3.json`,
`build/kls_lowwork_skip_fast_activsg10k_t4_r1_ref3.json`, and
`build/kls_lowwork_skip_fast_activsg70k_t4_r1_ref3.json`). Replacing those
three rows in the current medium projection gives
`build/kls_lowwork_skip_fast_activs_medium_projection_t4_r1_ref3_timeout120.jsonl`,
with KLS geomean `0.56554s` versus CKTSO `0.30764s` (`1.838x`). This is a
real movement from the previous `1.895x` projection, but it also confirms the
dominant remaining loss is still the broader CKTSO numeric kernel gap, not
just low-work fast-factor dispatch.

The BLAS-small-case hypothesis was checked against the current build before
changing policy: this build reports `build_has_cblas=false`, and the previous
focused rows reported zero CBLAS update calls. The retained change is
therefore diagnostic rather than another executor toggle. Public stats and
`kls_bench --json` now classify retained cached-supernode probe misses into
shape/materialization, sparse-stream contiguity, work gate, and workspace
pressure. This directly tests the SubtreeLU/CKTSO gap around persistent
supernode numeric panels and producer/consumer reuse without making a
case-specific timing change.

Focused runs with
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_EXEC=1`, 4 threads, one factor
repeat, and three refactor repeats were saved as
`build/kls_retained_plan_diag_asic100ks_t4_r1_ref3.json` and
`build/kls_retained_plan_diag_onetone2_t4_r1_ref3.json`. On `ASIC_100ks`,
the cached probe attempted 5468 runs, hit a panel every time, found 2369
contiguous streams, and applied only 267 runs / 39436 rows. The misses were
2750 shape/materialization rejects, 349 stream rejects covering 61106 rows,
2102 work-gate rejects covering 4666 rows, and zero workspace rejects. On
`onetone2`, the probe attempted 9309 runs and applied 235 runs / 26108 rows;
the misses were 1241 shape rejects, 6837 stream rejects covering 939932 rows,
996 work-gate rejects covering 3052 rows, and zero workspace rejects. The
large remaining gap is therefore not small BLAS dispatch or workspace
allocation. It is the paper-level missing piece: make retained numeric panels
usable for non-contiguous consumer streams and for the currently
shape-rejected producer/consumer cases, then revisit the work gate once that
coverage is materially higher.

A direct partial-prefix experiment was rejected. Letting the cached executor
consume only the contiguous prefix of a non-contiguous probe did reduce stream
rejects and increased applied cached rows on `onetone2`, but the focused
retained-plan runs got slower, so this is not the missing paper mechanism.
The retained source instead skips cached probes for panel ids whose panel
width has already been pruned to zero. This does not change any accepted
numeric update, but it removes invalid-panel probe overhead exposed by the
new counters. With
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_EXEC=1`,
`build/kls_valid_panel_probe_asic100ks_t4_r1_ref3.json` kept the same 267
applied cached runs / 39436 rows while reducing cached probe attempts from
5468 to 734 and shape rejects from 2750 to 15. The matching `onetone2` run
kept 235 applied runs / 26108 rows while reducing attempts from 9309 to 7280
and shape rejects from 1241 to 22. The remaining misses are still the real
paper gap: non-contiguous consumer streams (`onetone2`) and low-coverage
retained plan application, not invalid-panel probing.

A stricter gathered cached-panel solve was also tested and rejected. The
prototype solved the whole retained panel slice in worker scratch and would
commit only if every internal U row missing from the fixed KLU-compatible U
pattern computed exactly zero. This preserves correctness, because committing
nonzero values for missing U rows would require extra U storage and solve
support. Focused runs with both
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_EXEC=1` and
`KLS_ENABLE_REFACTOR_GATHERED_SUPERNODE=1` on `ASIC_100ks` and `onetone2`
kept the same applied cached rows and the same stream-reject counts as the
valid-panel-probe baseline. The non-contiguous streams therefore are not
recoverable inside the current KLU-compatible U pattern by a guarded gathered
panel consumer. Closing this paper gap requires changing the native numeric
object/U storage, or routing these cases through a row-major supernodal
refactor that owns the additional internal U entries.

A current-source CBLAS-capable rerun confirms that adding another "large cases
only" BLAS guard would be a no-op for the focused retained-panel slow cases.
`build-cblas` was rebuilt and run with `OPENBLAS_NUM_THREADS=1`,
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_EXEC=1`, and
`KLS_ENABLE_CBLAS_SUPERNODE=0/1` on `ASIC_100ks` and `onetone2`
(`build-cblas/kls_cblas_guard_plan_{off,on}_{asic100ks,onetone2}_t4_r1_ref3.json`).
Both matrices reported `build_has_cblas=true`, but every run had zero
`refactor_supernode_cblas_update_*` calls. The retained cached-panel coverage
was unchanged by the runtime CBLAS gate: `ASIC_100ks` kept 734 cached probes,
267 blocked updates / 39436 rows, 347 stream rejects / 61101 rows, and zero
workspace rejects; `onetone2` kept 7280 cached probes, 235 blocked updates /
26108 rows, 6837 stream rejects / 939932 rows, and zero workspace rejects.
The existing source already requires the runtime CBLAS gate plus 512-scale
shape checks and multi-million-operation work thresholds, and these focused
paths fall through to the in-KLS blocked cached-panel kernel. The active
paper gap remains retained row/segment numeric storage for non-contiguous
consumer streams, not small BLAS dispatch.

The retained consumer-plan executor now stays plan-only when
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_EXEC=1` is the gate that enables
supernode numeric updates. Previously the executor used the retained plan for
known producer/consumer runs, then still tried opportunistic full-panel probes
for non-plan dependencies; on the focused rows those extra probes did not
produce accepted updates and dominated the stream-reject bucket. The new guard
returns to scalar handling before recording a cached-panel probe when no
retained consumer-plan run exists. Explicit
`KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1` still keeps the broader opportunistic
cached-panel experiment.

Focused plan-only reruns stayed residual-clean and preserved all accepted
cached rows while removing non-plan stream rejects. On `onetone2`,
`build/kls_plan_only_probe_onetone2_t4_r1_ref3.json` kept 235 accepted cached
updates / 26108 rows, but cached attempts dropped from 7280 to 255 and stream
rejects dropped from 6837 / 939932 rows to zero. On `ASIC_100ks`,
`build/kls_plan_only_probe_asic100ks_t4_r1_ref3.json` kept 267 accepted
updates / 39436 rows, while cached attempts dropped from 734 to 281 and stream
rejects dropped from 347 / 61101 rows to zero. The remaining misses are 20
work-gate rejects on `onetone2` and 14 on `ASIC_100ks`; the large row-count
gap is still the missing native row/segment storage, but the plan executor no
longer pays for non-plan KLU-stream probes that cannot implement that storage.

A fresh native row-panel policy probe on the current source shows why the
existing selector still needs a better structural signal before promotion.
Forced row-refactor runs on the first five CKTSO-gap manifest rows compared
`KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=0`, `auto`, and `1` with
`--row-refactor refactor`, one factor repeat, and two refactor repeats
(`build/kls_native_policy_{off,auto,on}_gap5_t4_r1_ref2_timeout120.jsonl`).
The same-source geomeans were `10.3462s`, `10.0424s`, and `9.4441s`,
respectively, and all rows stayed residual-clean. `auto` improved the large
ASIC rows (`ASIC_320ks` `0.9065x`, `ASIC_320k` `0.9157x`, `ASIC_100ks`
`0.9720x` versus native-off) but regressed `rajat03` by `1.629x`; forced `1`
was best geomean in this pass but still regressed `rajat03` and `ASIC_100ks`.
This is promising coverage movement, not enough evidence for a default change.

A narrower auto rule that rejected unchecked pure-internal dense panels
(`trailing_len == 0`) was prototyped and rejected by smoke coverage. It would
avoid the `rajat03` auto regression shape, but the checked-row prefix repair
fixture and the partial-prefix fixture both require zero-trailing compact
panels to publish the expected retained panel state. That means trailing
presence alone is not the right reuse proxy. The next native selector needs to
distinguish no-reuse pure panel factorization from zero-trailing panels that
feed checked-prefix or partial-prefix publication.

The follow-up implementation changed native row-panel `auto` from a fixed
compact-work threshold into a reuse-retiring policy. In `auto`, KLS now tries
native row-panel execution for width > 1, but after a successful unchecked
refactor it disables native row panels for later refactors of the same
row-refactor pattern when that run used native panels and produced neither
compact-supernode updates nor a partial-supernode pipeline. Checked refactors
and partial-prefix/pipeline runs are intentionally exempt, preserving the smoke
coverage that rejected the trailing-only rule. The benchmark output now exposes
`row_refactor_native_row_panel_auto_disabled` and
`row_refactor_native_row_panel_auto_disable_count`.

This is also the current answer to the small-BLAS hypothesis. The source
already gates external CBLAS behind the build option, the runtime
`KLS_ENABLE_CBLAS_SUPERNODE=1` switch, 512-scale row/panel tests, and
multi-million-operation work thresholds. The focused runs below used the
default build (`build_has_cblas=false`), so a stricter "BLAS only for large
cases" guard cannot explain or fix these timings.

Validation after the selector change:

- `cmake --build build -j4` completed.
- `ctest --test-dir build --output-on-failure` passed both tests.
- Native-off top-five forced row-refactor control:
  `build/kls_native_adaptive_off_gap5_t4_r1_ref3_timeout120.jsonl`,
  geomean `10.0890s`.
- Broad native-auto:
  `build/kls_native_adaptive_broad_auto_gap5_t4_r1_ref3_timeout120.jsonl`,
  geomean `8.1896s`, five wins over native-off, geomean ratio `0.8117`.
- Forced native rerun:
  `build/kls_native_adaptive_broad_on_gap5_t4_r1_ref3_timeout120.jsonl`,
  geomean `8.4649s`; broad auto is within normal run noise and slightly
  faster in this pass. An earlier forced-native artifact at `4.3748s` did not
  reproduce and should be treated as an outlier.

The top-five broad-auto rows did not trip auto-disable: ASIC rows produced
compact-supernode updates, while `gemat12` and `rajat03` used the partial
pipeline. A controlled `rajat03` run with
`KLS_ENABLE_PARTIAL_SUPERNODE_PIPELINE=0`
(`build/kls_native_adaptive_broad_auto_no_partial_rajat03_t4_r1_ref3_timeout120.jsonl`)
confirmed the retirement path: after three refactors it reported
`row_refactor_native_row_panel_auto_disabled=1`,
`row_refactor_native_row_panel_auto_disable_count=1`, six native panels / 200
rows from the first run, and zero compact-supernode or partial-pipeline reuse.

A same-source rerun then checked whether the remaining top-ten gap could be
explained by small BLAS calls. It cannot: the default benchmark binary again
reported `build_has_cblas=false`, and every focused row reported zero external
CBLAS update counters. The optional CBLAS paths are already large-case gated
behind the build option, the runtime `KLS_ENABLE_CBLAS_SUPERNODE=1` switch,
512-scale shape checks, and multi-million-operation work thresholds, so an
additional "BLAS only for large cases" guard would be a no-op for these rows.

The retained code change is instead a narrow structural auto-row-refactor
escape hatch for small dominant-BTF cases. KLS now lets the row-refactor
pattern pass the lower-bound and work-ratio gates when the matrix is
moderately sized, has no off-diagonal pivots, and has one dominant BTF block;
that shape is where the EGraph refactor was losing to the forced row path.
For that internally selected shape, native row-panel `auto` is enabled without
changing the global `KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR` default. Diagnostic
save/restore now also carries the internal native-panel auto flag with the
row-refactor auto state.

Validation for the small dominant-BTF gate:

- `cmake --build build -j4` completed.
- `ctest --test-dir build --output-on-failure` passed both tests.
- Single-pass top-ten candidate:
  `build/kls_small_dombtf_rowauto_gap10_t4_r1_ref3_timeout120.jsonl`,
  geomean `6.3066s`; `rajat03` switched to `last_refactor_path=row_refactor`,
  used 18 native panels / 600 rows, and reported zero CBLAS calls.
- Same-binary row-refactor-disabled control:
  `build/kls_rowauto_disabled_samebin_gap10_t4_r1_ref3_timeout120.jsonl`,
  geomean `4.5419s`; this run had large unchanged-row timing swings, so it was
  not used alone for the policy decision.
- Pass-3 median candidate:
  `build/kls_small_dombtf_rowauto_gap10_t4_r1_ref3_pass3_timeout120.jsonl`,
  geomean `6.0088s`; `rajat03` measured `0.7027s`.
- Pass-3 median row-refactor-disabled control:
  `build/kls_rowauto_disabled_samebin_gap10_t4_r1_ref3_pass3_timeout120.jsonl`,
  geomean `5.9559s`; `rajat03` measured `3.4663s`.

This is a local improvement, not a CKTSO-gap closer. Against
`build/cktso_paper_medium93_t4_timeout120.jsonl`, the pass-3 candidate remains
`6.409x` slower on the common top-ten set, with `rajat03` still `13.98x`
slower than CKTSO despite the row-path win. The missing paper-scale piece is
therefore still the durable row/segment numeric storage and executor that
creates reusable coarse work broadly, not BLAS dispatch size.

A follow-up fast-factor policy check found a broader repeated-refactor issue.
The benchmark's SPICE-cycle score ignores repeated `kls_factor` time but runs
one `kls_factor` before the timed refactor loop, so the numeric state left by
KLS fast-factor repair can still affect the 99-refactor projection. On the
current top-ten focus, a same-source `--no-fast-factor` pass-3 control
(`build/kls_no_fast_factor_current_gap10_t4_r1_ref3_pass3_timeout120.jsonl`)
measured `4.7280s` geomean versus `6.0088s` for the previous default
(`build/kls_small_dombtf_rowauto_gap10_t4_r1_ref3_pass3_timeout120.jsonl`).
The improvement was not clean enough to turn fast factor off globally:
`onetone2` regressed in that pass, and direct paired samples showed large run
variance on rows whose final factor/refactor paths were otherwise identical.

The retained change is therefore a structural guard, not a global default flip.
When an existing EGraph refactor schedule has at least `1e8` dependency-work
units and at least 95% of that work is in the pipeline region, KLS now skips
the in-place fast-factor repair attempt and falls back to rebuilding the KLU
numeric object. This leaves cold first factors unchanged, keeps the small
dominant-BTF `rajat03` row-refactor path eligible for KLS fast factor, and
keeps lower-work cases such as `onetone2` outside the new guard. The guard is
recorded internally as `KLS_FAST_FACTOR_FAIL_PIPELINE_REFACTOR_GUARD` when
stats are read immediately after `kls_factor`.

Validation for the high-work pipeline guard:

- `cmake --build build -j4` completed.
- `ctest --test-dir build --output-on-failure` passed both tests.
- Single-pass top-ten guard:
  `build/kls_fast_pipeline_guard_gap10_t4_r1_ref3_timeout120.jsonl`,
  geomean `4.3822s`, nine wins over the previous single-pass default
  (`6.3066s`), and one noisy `gemat12` loss.
- Pass-3 top-ten guard:
  `build/kls_fast_pipeline_guard_gap10_t4_r1_ref3_pass3_timeout120.jsonl`,
  geomean `4.6536s`, eight wins over the previous pass-3 default
  (`6.0088s`), with `onetone2` and `rajat25` as noisy losses.
- First-20 focus screen:
  `build/kls_fast_pipeline_guard_gap20_t4_r1_ref3_timeout120.jsonl`,
  geomean `3.5419s`; no rows failed under the 120-second cap.

Against CKTSO, the pass-3 top-ten guard remains `4.964x` slower, and the
single-pass first-20 focus remains `4.027x` slower. This confirms the guard is
worth retaining as repeated-refactor policy cleanup, but it still does not
replace the missing CKTSO/SubtreeLU-style row/segment numeric executor.

A June 29 current-source rerun after the fast-factor guard gives the next
baseline for that executor work. The fresh top-ten default artifact
(`build/kls_current_baseline_gap10_t4_r1_ref3_timeout120.jsonl`) measured
`3.5428s` geomean with no failures, improving over the saved guard artifact
but still `3.779x` slower than `build/cktso_paper_medium93_t4_timeout120.jsonl`
on the same ten rows. Nine of the ten rows were still dominated by repeated
refactor time. The slow EGraph rows again reported
`row_refactor_auto_model_recommended=1`,
`row_refactor_auto_model_attempted=1`, and
`row_refactor_auto_lower_bound_rejected=1`, with lower-bound work close to the
exact EGraph work. Earlier forced-row artifacts and the current focused probes
therefore still apply: bypassing the gate is not the missing fix because the
current row-major executor does more scalar work than the EGraph path.

Current-source opt-in checks also reject promoting the existing supernode
prototypes. With `KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=cached`, the top-five
artifact `build/kls_current_cached_egraph_gap5_t4_r1_ref3_timeout120.jsonl`
measured `2.3725s` geomean, while the same-source default top-five control
`build/kls_current_baseline_gap5_t4_r1_ref3_timeout120.jsonl` measured
`1.5925s`. The consumer-plan executor
(`build/kls_current_consumer_plan_gap5_t4_r1_ref3_timeout120.jsonl`) measured
`2.4539s`. On the top-ten set, full EGraph supernode updates
(`build/kls_current_full_egraph_supernode_gap10_t4_r1_ref3_timeout120.jsonl`)
measured `4.1474s`, and the U-supernode ragged-L path
(`build/kls_current_u_ragged_gap10_t4_r1_ref3_timeout120.jsonl`) measured
`3.8816s`, both behind the default. These results keep the existing
paper-aligned supernode paths opt-in: they validate retained producer/consumer
metadata, but they do not yet amortize enough scalar work to be production
defaults.

A thread-count sweep reinforces that the retained EGraph path is parallel-work
limited but also has a high-thread synchronization boundary. On the same
top-ten slice, one, four, eight, sixteen, and thirty-two KLS threads measured
`6.0093s`, `3.5428s`, `2.9934s`, `2.3307s`, and `18.0589s` geomean in the
`build/kls_current_baseline_gap10_t{1,4,8,16,32}_r1_ref3_timeout120.jsonl`
artifacts. Sixteen threads still loses to CKTSO's four-thread artifact by
`3.193x`, and the thirty-two-thread collapse is hardware/synchronization
sensitive, so KLS did not add a CPU-dependent hard cap. The useful conclusion
is narrower: the current scalar EGraph executor can use more parallelism, but
it still needs coarser row/supernode numeric work per requested thread to close
the CKTSO gap fairly.

The only retained source cleanup from this pass removes inactive work from the
default EGraph hot loops. U-supernode numeric-value recording is an opt-in
paper prototype behind `KLS_ENABLE_REFACTOR_U_SUPERNODE_VALUES=1`, yet the
default scalar EGraph dependency loops still called the recorder for every
dependency and pivot, where it immediately returned because the value cache was
disabled. KLS now makes that recorder explicitly inline and guards those call
sites on `shared->u_supernode_values`; the opt-in value-cache smoke still
passes with `KLS_ENABLE_REFACTOR_U_SUPERNODE_VALUES=1`. The guarded default
artifact
`build/kls_no_u_value_record_default_gap10_t4_r1_ref3_pass3_timeout120.jsonl`
measured `3.4767s` pass-3 geomean with zero U-supernode value writes and no
failures, versus `4.6536s` for the previous saved pass-3 guard artifact. The
single-pass A/B was mixed (`3.6672s` after the guard versus `3.5428s` for the
fresh pre-change control), so the cleanup is not a standalone speed claim. It
is still `3.708x` slower than CKTSO on the same top-ten rows. This is retained
as low-risk cleanup on the active EGraph path, not as a claim that hot-loop
cleanup replaces the missing producer-centered row/supernode executor.

A fresh same-binary CBLAS gate check on June 29, 2026 reconfirms that "use
BLAS only for large cases" is already the current policy and is not active on
the CKTSO-gap focus rows. The normal `build` tree has
`KLS_ENABLE_CBLAS_SUPERNODE=OFF`; the separate `build-cblas` tree was rebuilt
with `KLS_ENABLE_CBLAS_SUPERNODE=ON` and run with `OPENBLAS_NUM_THREADS=1`.
With the runtime gate off,
`build-cblas/kls_cblas_gate_off_gap10_t4_r1_ref3_current_timeout120.jsonl`
measured `3.4664s` geomean. With `KLS_ENABLE_CBLAS_SUPERNODE=1`,
`build-cblas/kls_cblas_gate_on_gap10_t4_r1_ref3_current_timeout120.jsonl`
measured `3.6931s` geomean. Both artifacts reported
`build_has_cblas=true`, but every row had zero
`refactor_last_supernode_cblas_update_*` and zero
`refactor_supernode_cblas_update_*` counters. Eight rows still used the EGraph
refactor path, one used the mapped path, and one used row refactor. Therefore
an additional small-case BLAS guard would be a no-op for this slow set; the
gap remains in creating and scheduling reusable coarse row/supernode numeric
work rather than in external BLAS call granularity.

The next paper-aligned EGraph supernode cleanup moves the cached-panel pruning
decision before numeric panel allocation. Previously
`kls_build_refactor_supernode_panel_cache` allocated dense and trailing value
storage for every structural candidate and then pruned to the panels that the
U dependency stream or retained consumer plan could actually consume. On rows
such as `ASIC_320k`, this meant hundreds of candidate panels could be staged
even when zero cached updates would run. The builder now records lightweight
candidate metadata first, marks only candidates whose dependency run passes
the same contiguity and work tests used by the executor, and allocates dense
and trailing numeric storage only for those marked candidates. If an opt-in
supernode mode requests the cache but the marked set is empty, KLS frees the
empty cache and disables that opt-in mode for later refactors of the same
numeric object.

Validation kept the change scoped to the experimental supernode paths:
`cmake --build build -j4`, `ctest --test-dir build --output-on-failure`,
`KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1 ./build/kls_smoke`,
`KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=cached ./build/kls_smoke`, and
`git diff --check` all passed. The single-pass top-ten artifacts show the
intended staging effect. Default mode stayed free of supernode-cache work in
`build/kls_used_panel_cache_default_gap10_t4_r1_ref3_timeout120.jsonl`
(`3.6069s` geomean). Cached mode in
`build/kls_used_panel_cache_cached_gap10_t4_r1_ref3_timeout120.jsonl`
reduced no-use rows to zero retained panels and one disable, but still measured
`3.6520s`, so it remains opt-in. Full supernode mode improved the focused
top-ten artifact to `2.6895s` in
`build/kls_used_panel_cache_full_gap10_t4_r1_ref3_timeout120.jsonl`, with
usable rows retaining exactly one panel instead of hundreds of staged
candidates; however, the broader top-20 promotion check rejected making this
the default. `build/kls_used_panel_cache_full_gap20_t4_r1_ref3_timeout120.jsonl`
measured `2.2795s`, while the same-source default
`build/kls_used_panel_cache_default_gap20_t4_r1_ref3_timeout120.jsonl`
measured `2.0286s`. The retained value is therefore narrower: the
SubtreeLU/CKTSO-style cached-panel scaffold now avoids paying for panels that
cannot be consumed, but the scalar update executor still lacks enough
arithmetic intensity to be a general default.

The same source was then rerun with the retained supernode consumer-plan
executor after the used-panel cache change. The plan-only path avoids the
ordinary cached-probe stream rejects by looking up retained producer/consumer
runs directly, and it now benefits from the same pre-allocation panel pruning.
It improved relative to the older pre-pruning consumer-plan artifact but still
did not clear the broader default gate: the top-ten plan-only artifact
`build/kls_used_panel_cache_consumer_plan_gap10_t4_r1_ref3_timeout120.jsonl`
measured `2.4451s`, and the top-20 artifact
`build/kls_used_panel_cache_consumer_plan_gap20_t4_r1_ref3_timeout120.jsonl`
measured `2.2768s`. Because the same-source default top-20 artifact above
measured `2.0286s`, the consumer-plan executor also remains opt-in. The
remaining gap is the executor arithmetic and writeback cost after a plan hit,
not just rejected-probe overhead or unused panel staging.

After the focused row-refactor loss shifted to `coupled` and `rajat03`, KLS
also reran the small-BLAS hypothesis on those rows specifically. The current
source already requires a CBLAS-capable build, `KLS_ENABLE_CBLAS_SUPERNODE=1`,
512-scale row/vector shape checks, and multi-million-operation work thresholds
before dispatching to external BLAS. A rebuilt `build-cblas` tree was run with
`OPENBLAS_NUM_THREADS=1` on just those two matrices. With the runtime gate off,
`build/kls_cblas_gate_off_coupled_rajat03_t4_r1_ref3_timeout120.jsonl`
measured `0.4246s` geomean; with the gate on,
`build/kls_cblas_gate_on_coupled_rajat03_t4_r1_ref3_timeout120.jsonl`
measured `0.4388s`. Both artifacts reported `build_has_cblas=true`, but the
external CBLAS counters remained zero. `coupled` used 81 native row panels over
2421 rows, while `rajat03` used 18 panels over 600 rows, so the native panels
average only about 30 to 33 rows, far below the CBLAS threshold. The compact
batched-supernode counters also stayed at zero on both rows. An additional
"BLAS only for large cases" guard would therefore not change these losses; the
active row-refactor gap is still in the scalar native/compact panel executor
and value writeback, not in small external BLAS dispatch.

The native row-panel off control itself was then corrected so
`KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=0` cannot be promoted back to structural
auto mode by `row_refactor_auto_native_row_panel`. Unset still keeps the
current structural default, `auto` requests the auto gate, and `1` forces
eligible retained panels. With the corrected explicit off state,
`build/kls_native_forced_off_coupled_rajat03_t4_r1_ref3_pass3_timeout120.jsonl`
measured `0.4508s` geomean on `coupled`/`rajat03`; native-panel counters were
zero as intended. The same-source default rerun in
`build/kls_default_after_native_env_fix_coupled_rajat03_t4_r1_ref3_pass3_timeout120.jsonl`
measured `0.4158s` geomean and still reported zero external CBLAS calls. This
rejects a scalar/off guard for these small panels; the retained native panel
path is not the main reason for the CKTSO gap on the focused pair.

The same focused rows were then checked against the SubtreeLU Algorithm 6
private/pipeline queue idea. The default artifacts showed
`separator_analyzed_rows=0` for both `coupled` and `rajat03`; KLS's retained
`METIS_NodeNDP` separator map is currently limited to 30,000+ row threaded
METIS analyses, while these row-refactor losses have only 11,341 and 7,602
rows. A temporary source probe lowered `KLS_METIS_NDP_MIN_ROWS` to 4,096 and
rebuilt KLS. This did activate retained separator metadata and the
FLOP-balanced separator queue for METIS-ordered `rajat03` and forced-METIS
`coupled`, but it was a loss: auto ordering on the two-row focus measured
`0.5380s` geomean in
`build/kls_ndp4096_auto_coupled_rajat03_t4_r1_ref3_pass3_timeout120.jsonl`,
and forced METIS measured `0.6240s` in
`build/kls_ndp4096_metis_coupled_rajat03_t4_r1_ref3_pass3_timeout120.jsonl`.
For those runs, the separator queue turned almost the whole row-refactor group
graph into pipeline work (`row_refactor_group_pipeline_groups` equaled all
groups), increasing both analysis and repeated refactor time. The threshold
change was reverted; medium-matrix separator retention is not a general gap
closer in the current row executor.

The competing partial-supernode pipeline policy was also probed because the
default row-refactor rows showed partial-supernode activity but no retained
separator queue. With `KLS_ENABLE_PARTIAL_SUPERNODE_PIPELINE=0`, the focused
two-row median artifact
`build/kls_sepqueue_no_partial_coupled_rajat03_t4_r1_ref3_pass3_timeout120.jsonl`
improved `coupled` but slightly worsened `rajat03`, for a `0.4221s` geomean
versus `0.4482s` in the same-source default artifact
`build/kls_default_coupled_rajat03_t4_r1_ref3_pass3_timeout120.jsonl`.
The broader top-20 guard rejected making that policy a default:
`build/kls_ndp4096_no_partial_gap20_t4_r1_ref3_timeout120.jsonl` measured
`2.0638s`, worse than the saved same-source default top-20 `2.0286s`. The
retained conclusion is that KLS should not broadly disable the partial
supernode pipeline; the remaining CKTSO gap needs a better row/supernode
numeric executor, not just a different queue-selection switch.

The current source also rechecked whether forcing the KLS-owned first-factor
bridge could improve the large EGraph refactor rows enough to justify replacing
the accepted KLU first factor. It did not. With
`KLS_ENABLE_KLS_FIRST_FACTOR=1`,
`build/kls_forced_first_large3_t4_r1_ref3_timeout120.jsonl` measured
`7.4156s` on `ASIC_100ks`, `14.0002s` on `ASIC_320k`, and `22.7419s` on
`G2_circuit`, versus the saved default values of `5.4732s`, `13.0825s`, and
`21.1155s`. The forced path changed `last_factor_path` to `kls_first`, but
initial factor time rose sharply (`ASIC_100ks` from `0.6480s` to `2.0697s`,
`ASIC_320k` from `1.2720s` to `2.1234s`, and `G2_circuit` from `1.0330s` to
`2.0720s`) while repeated EGraph refactor time stayed flat or slightly worse.
The KLS-first automatic selector therefore remains conservative for these
large CKTSO-gap rows. The missing piece is still the production
row/supernode numeric executor that lowers repeated refactor cost, not merely
seeding the existing EGraph path from the current KLS-first bridge.

A fresh large-paper reconnaissance on June 29, 2026 keeps the same conclusion
with current source. `build/kls_current_large_recon_t4_r1_ref1_timeout120.jsonl`
completed six of the eight selected large rows at `39.0334s` geomean and timed
out on `pre2` and `Hamrle3` under the 120s process cap. Against the fresh CKTSO
large-recon artifact, charging failures as 1000s, KLS measured `1.060x` slower:
it won `TSOPF_FS_b39_c30` (`82.99s` versus `334.19s`), tied CKTSO's `Hamrle3`
timeout, but lost `pre2` because CKTSO completed it (`509.50s` projected
cycle) while KLS timed out. The completed large losses are now concentrated in
`rajat30` (`37.33s` versus `22.15s`), `nxp1` (`36.17s` versus `25.34s`),
`ASIC_680k` (`7.10s` versus `5.89s`), `G3_circuit` (`809.24s` versus
`749.95s`), and the near-tie `rajat29`.

Focused current-source probes on the completed-but-slow `nxp1`/`rajat30` pair
rejected the obvious paper-path bypasses. Forcing row refactor in
`build/kls_forced_row_large_nxp1_rajat30_t4_r1_ref1_timeout120.jsonl`
measured `151.35s` and `131.69s`, far behind the default `36.17s` and
`37.33s`; this validates the existing row lower-bound rejection on these large
single-block cases. Full EGraph supernode updates in
`build/kls_full_supernode_large_nxp1_rajat30_t4_r1_ref1_timeout120.jsonl`
measured `48.18s` and `40.56s`. The retained consumer-plan executor in
`build/kls_consumer_plan_large_nxp1_rajat30_t4_r1_ref1_timeout120.jsonl`
was mixed, losing `nxp1` (`49.34s`) while slightly improving `rajat30`
(`36.32s`), for a net geomean loss. These results keep the current scalar
supernode prototypes opt-in even on large rows. The next implementation target
is still the production row/supernode numeric representation that can store
and update the extra non-contiguous row-major state directly, plus a separate
first-factor path for `pre2`; the current toggles do not close either gap.

A follow-up `pre2` rerun added explicit dominant-BTF first-factor coverage
diagnostics instead of inferring path selection from analysis-only output. The
new stats record whether the KLS-first row-up-looking dominant-BTF pipeline was
used, which BTF block it targeted, how many rows that block contained, and
whether retained separator metadata covered it; `KLS_TRACE_KLS_FIRST_FACTOR=1`
prints the same decision before a long factor run can time out. A fast
`transient` forced-first probe confirmed the counters are wired, reporting the
178,823-row dominant block with no separator coverage under AMD. On `pre2`, the
default AMD forced-first run emitted
`block=13843 rows=629628 n=659033 nblocks=29282 separator=0` and timed out at
the 120s cap. The forced METIS run emitted the same block and row count with
`separator=1`, but it also timed out at 120s. The local CKTSO comparison binary
completed the matching `pre2` run in about 21s wall time, with
`analysis_seconds=3.657411`, `initial_factor_seconds=6.308607`,
`factor_seconds_avg=5.083979`, and `refactor_seconds_avg=4.893165`. This
rejects both the timeout-limit and ordering-only explanations for `pre2`: KLS
is entering the correct dominant block, and METIS can provide separator
coverage, but the current KLS numeric executor still does not implement the
coarse CKTSO/SubtreeLU row/supernode work inside that block.

The next `pre2` probe sampled the METIS path under `gdb` and found a more
specific serialization point inside that numeric executor. The interrupted run
had one pipeline worker in `kls_row_first_supernodes_reset()` while the other
pipeline workers waited on the same mutex. That reset was rebuilding the
completed-prefix row-supernode map after a pivot; in the separator-queue case,
the completed prefix is almost the whole 629,628-row dominant block. KLS now
keeps the small-prefix behavior but, once the completed prefix reaches the
existing 32k-row prefix-cache rebuild cutoff, invalidates the speculative
row-supernode accelerator and continues with scalar dependency updates instead
of rescanning the huge prefix under the pipeline lock. A fast forced-first
`transient` run still completed with a valid residual (`relative_residual_l2`
about `3.44e-13`) and comparable timing. On `pre2`, the METIS forced-first
120s run now emitted the dominant-block trace a second time, showing the first
factor completed and the measured factor started; before this change the same
bounded run timed out before that point. A 240s `pre2` METIS factor-only run
still timed out during the measured factor, however, and the new stack sample
moved to `kls_row_first_partial_apply_one_dep()` with other pipeline workers
waiting. The default AMD forced-first path still timed out before a second
trace. This narrows the next CKTSO/SubtreeLU gap further: after avoiding the
large prefix supernode reset, KLS still serializes too much dependency-row
numeric update work behind the pipeline mutex.

The follow-up retained-panel path keeps completed dense/common-tail row-up
panels alive after speculative row-supernode metadata is disabled, invalidates
stale panels across dynamic column pivots, and lets later pipeline rows consume
those panels before falling back to scalar row-entry streaming. This is a
correct SubtreeLU-shaped bridge but not enough for the hard first-factor case.
On forced-METIS `pre2`, a capped trace still reached the pivot tail at row
274,430 with about 15.3B scalar U entries after applying about 260k
panel-backed update groups over roughly 2.0M rows. The new scalar-stream split
shows why this is not just a missing ready-supernode dependency run: about
5.0B of those scalar entries update still-pending dependency rows, while about
10.3B update the current row's output/trailing pattern. On
`Freescale/transient`, the repeated pass is much smaller but shows the same
split direction after retained panels: about 1.47M internal entries and 2.57M
output/trailing entries. The next paper-level implementation should therefore
target producer-to-current-row output streaming through a broader row-major
numeric object, producer/output accumulator, or comparable coarse supernode
executor; a dependency-only grouping change or another BLAS threshold would
leave most of the measured `pre2` scalar stream intact.

A same-binary CBLAS guard rerun confirmed that this conclusion still applies
to the current top CKTSO-gap medium rows. `build-cblas` was rebuilt and run
with `OPENBLAS_NUM_THREADS=1` on the first five
`bench/suitesparse_cktso_gap_manifest.txt` entries. With
`KLS_ENABLE_CBLAS_SUPERNODE=0`,
`build-cblas/kls_cblas_guard_current_off_gap5_t4_r1_ref3_timeout120.jsonl`
measured a `1.76679s` SPICE-cycle geomean; with
`KLS_ENABLE_CBLAS_SUPERNODE=1`,
`build-cblas/kls_cblas_guard_current_on_gap5_t4_r1_ref3_timeout120.jsonl`
measured `1.61226s`. Both artifacts reported `build_has_cblas=true`, but all
five rows had zero `refactor_last_supernode_cblas_update_*` and zero
`refactor_supernode_cblas_update_*` counters. The timing variation is therefore
not caused by executing small BLAS kernels. The current source already
implements the proposed "BLAS only for large cases" policy for these paths;
tightening the threshold again would be a no-op on the focused losses.

A fresh current-source rerun rejected the earlier small-dominant-BTF
row-refactor escape hatch. That escape let row-refactor run when its cheap
lower-bound or retained group-work estimate was up to `1.35x` the exact EGraph
dependency work, on the assumption that row-major native/partial-supernode
updates would amortize the extra scalar work. The current top-gap artifacts show
the opposite for the retained KLU-compatible row scaffold: `rajat03` and
`coupled` both switched to `last_refactor_path=row_refactor` under that escape
and lost badly even though `coupled` reported compact-supernode activity. A
same-session top-twenty control measured `2.12037s` geomean
(`build/kls_current_control_gap20_t4_r1_ref3_timeout120.jsonl`); forcing
row-refactor off measured `2.05920s`
(`build/kls_row_refactor_off_gap20_t4_r1_ref3_timeout120.jsonl`) while moving
`rajat03` and `coupled` back to EGraph.

KLS now applies the row-refactor lower-bound gate uniformly: if the cheap
lower-bound work exceeds the exact EGraph dependency work, row metadata setup is
rejected even for small dominant-BTF matrices; if the full retained row-pattern
work is built, it must be no larger than the EGraph dependency work before
automatic row-refactor execution is allowed. This follows the paper-level
principle more directly: the current row-major scaffold should only replace the
exact EGraph refactor when it demonstrably reduces work or creates broad enough
coarse kernels, not merely because the shape is small.

Validation for the strict row lower-bound gate:

- `cmake --build build -j4` completed.
- `ctest --test-dir build --output-on-failure` passed both tests.
- Top-ten CKTSO-gap focus:
  `build/kls_row_lb_strict_gap10_t4_r1_ref3_timeout120.jsonl`, geomean
  `2.31589s`, versus the same top-ten slice from the fresh top-twenty control
  at `2.36727s`.
- Top-twenty CKTSO-gap focus:
  `build/kls_row_lb_strict_gap20_t4_r1_ref3_timeout120.jsonl`, geomean
  `1.92561s`, `0.9081x` of the fresh KLS control geomean. It had 11 wins, 3
  ties, and 6 losses over 2% versus
  `build/kls_current_control_gap20_t4_r1_ref3_timeout120.jsonl`.
- The intended rows now reject at the lower-bound gate and stay on EGraph:
  `rajat03` reports `row_refactor_auto_lower_bound_work=2110490` versus
  `refactor_dependency_work=2018944`, `row_refactor_auto_lower_bound_rejected=1`,
  and `last_refactor_path=egraph`; `coupled` reports
  `12186383` versus `11933363`, lower-bound rejected, and EGraph.

Against `build/cktso_paper_medium93_t4_timeout120.jsonl`, the strict top-twenty
run is still `2.1896x` slower than CKTSO on the common rows. This change removes
a self-inflicted selector loss, but it does not close the main CKTSO gap: the
remaining losses are still the EGraph scalar update kernels and the missing
production row/supernode numeric executor.

An EGraph tail-wait experiment tested whether CKTSO's cluster/pipeline split was
leaving avoidable dependency waits in KLS. The idea was valid: cluster levels
complete under barriers before the pipeline tail starts, so tail columns do not
need to spin on predecessors that belong to the completed cluster prefix. Three
variants were tried: an ungated pipeline-column flag, a boundary-edge gated flag,
and a simpler split-prefix flag enabled only when `cluster_levels > 0`.

The experiment was rejected as a source change. The ungated run
`build/kls_tail_wait_flag_gap20_t4_r1_ref3_timeout120.jsonl` measured
`1.87993s` geomean versus the strict row lower-bound baseline at `1.92561s`, but
the wins included all-pipeline matrices where the flag cannot skip any waits,
so that result was not a defensible algorithmic signal. The boundary-gated run
`build/kls_tail_wait_flag_diffgate_gap20_t4_r1_ref3_timeout120.jsonl` measured
`1.90098s` geomean, still had 7 losses over 2%, and remained `2.1615x` slower
than CKTSO. The simpler split-prefix run
`build/kls_tail_wait_splitflag_gap20_t4_r1_ref3_timeout120.jsonl` measured
`1.92149s`, effectively the same as strict KLS. The source was therefore left at
the strict row lower-bound gate rather than retaining a noisy scheduler tweak.

This also keeps the BLAS hypothesis bounded. The measured retained build reports
`build_has_cblas=false`, and the earlier same-binary CBLAS checks reported zero
CBLAS update counters even when CBLAS was enabled. The proposed "use BLAS only
for large cases" policy is already implemented for the optional CBLAS paths and
is not active in the current CKTSO-gap losses. The main gap is still coarse
numeric work aggregation inside the dominant blocks, not BLAS thresholding or
tail-wait bookkeeping.

The direct follow-up to the BLAS-threshold question repeated the same-binary
check on the current rebuilt CBLAS tree. With `OPENBLAS_NUM_THREADS=1` and the
runtime CBLAS gate off,
`build-cblas/kls_cblas_rerun_off_gap5_t4_r1_ref3_timeout120.jsonl` measured a
`1.43193s` top-five geomean. With `KLS_ENABLE_CBLAS_SUPERNODE=1`,
`build-cblas/kls_cblas_rerun_on_gap5_t4_r1_ref3_timeout120.jsonl` measured
`1.36371s`. The row-level movement was mixed (`ASIC_100ks` regressed while
`rajat03` improved), and every row in both artifacts reported
`build_has_cblas=true` with zero `refactor_last_supernode_cblas_update_*` and
zero cumulative `refactor_supernode_cblas_update_*` counters. Therefore the
runtime gate did not execute BLAS at all; adding another small-case BLAS guard
would not change this focused path.

A retained-plan lookup probe was also rejected. The experiment sorted each
consumer-plan current-column run list by dependency and replaced the executor's
linear run search with a binary lookup. It preserved smoke correctness
(`cmake --build build -j 4`, `ctest --test-dir build --output-on-failure`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_EXEC=1 ./build/kls_smoke` passed),
but it did not improve the focused benchmark:
`build/kls_consumer_plan_binlookup_gap5_t4_r1_ref3_timeout120.jsonl` measured
`1.42891s` versus the same-source default
`build/kls_binlookup_default_gap5_t4_r1_ref3_timeout120.jsonl` at `1.42766s`.
The plan executor still applied only the already-known limited cached-panel
subset before auto-disable, so lookup overhead is not the large paper gap. The
source change was reverted; the next useful target remains a broader
producer/output accumulator or row-major numeric object that can consume the
retained plan beyond the current common-tail cached-panel shape.

A follow-up U-supernode ragged-L run-length probe tested a more direct retained
plan use without the failed static pre-prune. The experiment built the retained
consumer plan whenever `KLS_ENABLE_REFACTOR_U_SUPERNODE_RAGGED_L=1` was set
and let the ragged-L executor cap each candidate producer run to the plan's
actual contiguous `(current, dependency)` run length before validating the U
stream. This preserved correctness (`cmake --build build -j 4`,
`ctest --test-dir build --output-on-failure`, and
`KLS_ENABLE_REFACTOR_U_SUPERNODE_RAGGED_L=1 ./build/kls_smoke` passed), but the
focused top-ten CKTSO-gap benchmark regressed:
`build/kls_ragged_planrun_gap10_t4_r1_ref3_timeout120.jsonl` measured
`3.52997s` geomean versus the same-source retained default
`build/kls_current_retained_gap10_t4_r1_ref3_timeout120.jsonl` at `2.28216s`,
with no wins over 2%. The counters showed no new broad coverage: the path
still updated the same limited ragged-L shapes (`ASIC_320ks` cumulative
`146,769` rows, `ASIC_100ks` `118,584`, `onetone2` `78,600`) and still applied
zero ragged-L rows on `ASIC_320k` and `rajat28`, while paying the retained-plan
construction cost. The source change was reverted. The useful conclusion is
that the retained plan's run length metadata alone does not create the missing
paper mechanism; KLS still needs the persistent producer/output accumulator or
row-major numeric object that lets one producer panel feed many consumers.

The small-BLAS guard hypothesis was then checked against a new opt-in default
EGraph stream diagnostic rather than another threshold change. KLS now reports
`refactor_stream_dependency_entries`, `refactor_stream_pivot_entries`, and
`refactor_stream_output_entries` when
`KLS_ENABLE_REFACTOR_STREAM_STATS=1` is set. The counters are populated while
building the exact EGraph dependency schedule and are otherwise inert; they
classify each producer L-column entry used by a U dependency as still-pending
dependency work, the current pivot row, or output/trailing work.

Validation for this diagnostic:

- `cmake --build build -j 4` completed.
- `ctest --test-dir build --output-on-failure` passed both tests.
- `KLS_ENABLE_REFACTOR_STREAM_STATS=1` on the current top-ten CKTSO-gap focus
  produced `build/kls_stream_stats_gap10_t4_r1_ref3_timeout120.jsonl` with a
  `2.29825s` SPICE-cycle geomean. That is `1.007x` of the retained default
  `build/kls_current_retained_gap10_t4_r1_ref3_timeout120.jsonl`, so the
  diagnostic did not materially perturb the path.
- The same diagnostic run remains `2.451x` slower than
  `build/cktso_paper_medium93_t4_timeout120.jsonl` on the common rows.

The stream split is the useful finding. Across the nine EGraph rows in that
top-ten set, KLS counted about `932.1M` dependency-side L entries, `8.8M` pivot
entries, and `921.4M` output/trailing entries. Output/trailing work is therefore
about `49.5%` of the measured scalar stream, while pivot entries are only about
`0.47%`. The hard ASIC rows are almost exactly balanced:
`ASIC_320ks` is `139.2M` dependency versus `139.3M` output/trailing entries,
and `ASIC_320k` is `170.1M` versus `169.8M`. This confirms that the current
loss is not caused by small external BLAS calls: the CBLAS artifacts still show
zero BLAS update counters, and the scalar EGraph stream itself is split across
dependency and output/trailing work. A dependency-only grouping tweak or another
BLAS threshold guard would leave roughly half of the stream untouched. The more
direct paper-aligned target remains a persistent producer/output accumulator,
row-major numeric object, or broad supernode executor that can reuse one
producer panel across many consumers and scatter/output rows without rebuilding
the same stream each time.

The follow-up ragged-L miss instrumentation explains why the existing
U-supernode ragged executor is not that missing mechanism. KLS now reports
last-pass ragged-L probe attempts, panel misses, short rejects, stream rejects,
and work rejects, plus an auto-disable flag/counter for the ragged-L executor.
The counters are only active when
`KLS_ENABLE_REFACTOR_U_SUPERNODE_RAGGED_L=1` requests that experimental path.

On the focused top-ten CKTSO-gap run before the auto-disable guard,
`build/kls_ragged_miss_stats_gap10_t4_r1_ref3_timeout120.jsonl` measured a
`3.02317s` geomean versus the retained default's `2.28216s`. The miss counters
showed the real cause: `9,238,545` of `9,268,826` last-pass ragged probes were
panel misses (`99.7%`). Short rejects, stream rejects, and work rejects were
tiny by comparison (`107`, `27,449`, and `724`). The executor did build many
candidate panels, but after the first pass it pruned almost all unused panels
and then kept probing the invalidated panel map on later SPICE refactors.

KLS therefore now auto-disables the ragged-L executor for a numeric object when
a clean pass applies less than one quarter as many ragged update rows as it
probes. This keeps the experimental path from repeatedly paying panel-miss
overhead while preserving the default EGraph path. Validation:

- `cmake --build build -j 4` completed.
- `ctest --test-dir build --output-on-failure` passed both tests.
- Final-source default top-ten run
  `build/kls_default_after_ragged_guard_gap10_t4_r1_ref3_timeout120.jsonl`
  measured `2.14452s` geomean and remains `2.287x` slower than CKTSO on the
  same ten rows, so the default path was not regressed but the CKTSO gap remains.
- Final-source ragged-L run
  `build/kls_ragged_autodisable_final2_gap10_t4_r1_ref3_timeout120.jsonl`
  measured `2.60687s`, `0.862x` of the raw ragged miss run but still
  `1.216x` slower than the default final-source EGraph run. It disabled
  ragged-L on 9 of 10 rows after one sparse-use pass and preserved
  `357,981` cumulative ragged update rows over `2,001` runs.

This is a useful guardrail, not the CKTSO closer. It directly rejects a
paper-shaped but too-narrow executor once measured coverage is sparse. The
large consumer plan still exposes the important opportunity (`6.3M` run rows
and about `1.75B` L entries on the same focus set), but the current ragged-L
executor realizes only a small fraction of it. The next implementation should
consume those retained runs with a broader row-major producer/output object
instead of probing one narrow panel map from every scalar dependency.

A same-source CBLAS guard rerun on June 29, 2026 confirms that the current
focused loss is still not caused by small BLAS calls. `build-cblas` was rebuilt
as a Release tree with `-DKLS_ENABLE_CBLAS_SUPERNODE=ON` and run on the saved
top-five gap manifest with `OPENBLAS_NUM_THREADS=1`. With
`KLS_ENABLE_CBLAS_SUPERNODE=0`,
`build-cblas/kls_cblas_small_guard_latest_off_gap5_t4_r1_ref3_timeout120.jsonl`
measured `1.38157s` geomean. With `KLS_ENABLE_CBLAS_SUPERNODE=1`,
`build-cblas/kls_cblas_small_guard_latest_on_gap5_t4_r1_ref3_timeout120.jsonl`
measured `1.32706s` geomean, but every row still reported zero
`refactor_last_supernode_cblas_update_*` and zero cumulative
`refactor_supernode_cblas_update_*` counters. The timing difference is
therefore run noise or unrelated branch effects, not evidence that BLAS
granularity is active. The existing CBLAS call sites already require the
runtime gate plus 512-scale row/vector or panel checks and multi-million-work
thresholds; another "large only" guard would not affect these rows.

The next direct supernode-gap prototype builds the U-supernode ragged-L panel
pattern from the retained consumer plan instead of from the narrower
`supernode_pipeline_end` candidate map. Set
`KLS_ENABLE_REFACTOR_U_SUPERNODE_PLAN_PATTERN=1` together with
`KLS_ENABLE_REFACTOR_U_SUPERNODE_RAGGED_L=1` to select this path. The schedule
now forces the retained consumer plan when the plan-pattern probe is requested,
then populates the existing ragged-L panel storage from each planned producer
panel's maximum retained run width. The executor also caps a ragged update to
the exact retained `(current, dependency)` plan run when the plan-derived
pattern is active. Existing default and old ragged modes still use the original
`supernode_pipeline_end` source.

Validation for the implementation:

- `cmake --build build -j4` completed.
- `ctest --test-dir build --output-on-failure` passed both tests.
- On the saved top-five CKTSO-gap focus, the current default measured
  `1.53003s` geomean in
  `build/kls_planpattern_default_gap5_t4_r1_ref3_timeout120.jsonl`, the old
  ragged path measured `1.61840s` in
  `build/kls_planpattern_oldragged_gap5_t4_r1_ref3_timeout120.jsonl`, and the
  plan-derived ragged path measured `1.63377s` in
  `build/kls_planpattern_ragged_gap5_t4_r1_ref3_timeout120.jsonl`.
- On the top-ten CKTSO-gap focus,
  `build/kls_planpattern_ragged_gap10_t4_r1_ref3_timeout120.jsonl` measured
  `2.68365s` geomean, `1.251x` slower than
  `build/kls_default_after_ragged_guard_gap10_t4_r1_ref3_timeout120.jsonl`.

This rejects a tempting direct reading of the retained-plan opportunity. The
top-ten plan-pattern run exposed `750,247` retained consumer runs over
`6,285,128` planned rows and about `1.754B` planned L entries, but it still
applied only `2,001` ragged updates over `357,981` rows before disabling on
nine of ten matrices. The retained plan therefore is not missing merely because
the panel map starts from `supernode_pipeline_end`; the larger missing
algorithm is the CKTSO/SubtreeLU-style batching step that groups many planned
small consumer runs against persistent producer/output storage. Single-run
ragged triangular updates remain too sparse and too probe-heavy even when their
panel ranges come from the broad consumer plan.

The retained consumer-plan diagnostics now distinguish single-run sparsity from
panel-level batching opportunity. `kls_stats` and `kls_bench` report planned
runs below the ragged single-run row floor, plus reused-panel batches whose
aggregate row count clears that floor:
`refactor_supernode_consumer_plan_small_run_*`,
`refactor_supernode_consumer_plan_batch_*`, and
`refactor_supernode_consumer_plan_batch_small_run_*`.

On the top-ten CKTSO-gap focus,
`build/kls_consumer_plan_batch_stats_gap10_t4_r1_ref3_timeout120.jsonl`
measured `2.26037s` geomean with plan construction enabled. It counted
`750,247` retained plan runs over `6,285,128` rows. Of those, `599,264` runs
and `1,826,985` rows are below the current ragged single-run row floor, which
explains why single-run ragged execution has poor coverage. But reused panels
contain `725,372` runs and `6,235,295` rows, including `574,389` small runs and
`1,777,152` small-run rows. This narrows the implementation target: the paper
gap is a panel-level grouped consumer executor that batches many small planned
runs against one producer panel/output accumulator, not a different single-run
threshold or a broader panel start map.

A follow-up retained-plan diagnostic checks whether that batching opportunity
is only cross-current producer-panel reuse, or whether enough work is visible
inside each current column to justify a smaller column-grouped executor first.
KLS now also reports current-column retained-plan distribution through
`refactor_supernode_consumer_plan_column_count`,
`refactor_supernode_consumer_plan_max_column_*`, and
`refactor_supernode_consumer_plan_column_batch_*`. The counters are computed
from the same retained plan and do not change numeric execution.

Validation:

- `cmake --build build -j4` completed.
- `ctest --test-dir build --output-on-failure` passed both tests.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN=1 ./build/kls_smoke` passed.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN=1` on the top-ten CKTSO-gap
  focus produced
  `build/kls_consumer_plan_column_batch_stats_gap10_t4_r1_ref3_timeout120.jsonl`
  with a `2.26381s` geomean and no failed matrices.

The new split is actionable. Across the same top-ten focus, the retained plan
still has `750,247` runs over `6,285,128` rows, with `599,264` small runs over
`1,826,985` rows. Producer-panel batches cover `725,372` runs and
`6,235,295` rows. Current-column batches cover almost the same numeric surface:
`35,873` current columns contain `670,653` retained runs and `6,106,008` rows,
including `519,716` small runs and `1,648,550` small-run rows. The largest
single current column batch is also substantial (`3,450` runs and `9,421`
planned rows on the Rajat rows). This means the next paper-aligned prototype
does not have to start with a fully cross-current producer scheduler. A
current-column grouped consumer executor can first consume many retained runs
for one active `x` workspace, then later promote the same retained producer
storage to cross-current reuse if the column-grouped path proves insufficient.
The existing executor does not do this yet: it still applies or rejects one
planned producer run at a time, so the broad column-batch surface is currently
only measured, not exploited.

An immediate column-batch gate prototype was tried and reverted. The probe made
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_COLUMN_EXEC=1` retain cached
producer panels for current columns whose total retained-plan rows cleared the
existing cached-supernode row floor, then allowed those planned small runs to
enter the cached panel path. It preserved correctness
(`cmake --build build -j4`, `ctest --test-dir build --output-on-failure`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_COLUMN_EXEC=1 ./build/kls_smoke`
passed), but it did not execute any cached updates on the focus rows. The
top-ten default measured `2.35504s` in
`build/kls_columnbatch_default_gap10_t4_r1_ref3_timeout120.jsonl`, the existing
plan executor measured `2.37290s` in
`build/kls_columnbatch_planexec_gap10_t4_r1_ref3_timeout120.jsonl`, and the
column gate measured `2.48611s` in
`build/kls_columnbatch_columnexec_gap10_t4_r1_ref3_timeout120.jsonl`. Both
plan-enabled runs still reported zero cached-probe attempts, zero applied
cached rows, and plan execution disabled on nine of ten rows.

The failed probe clarifies the missing implementation boundary. A column-batch
policy alone cannot help while the cached producer-panel object is still built
from the narrower `supernode_pipeline_end` candidate map. The next real
implementation must first materialize producer panels from the retained plan's
panel starts and column-batch extents, then run a grouped current-column
executor over that storage. Merely relaxing the individual-run work gate around
the old panel cache adds retained-plan construction overhead without creating
the missing producer storage.

The next diagnostic makes that storage boundary explicit. KLS now reports how
many retained-plan producer-panel starts could be materialized as the existing
strict cached-panel object, and how many of those pass the current
common-trailing predicate:
`refactor_supernode_consumer_plan_cached_panel_*` and
`refactor_supernode_consumer_plan_strict_cached_panel_*`. These counters are
computed when the cached panel builder runs and survive the old panel cache's
auto-disable path, so they describe the missing materialization opportunity
even when no cached updates execute.

Validation:

- `cmake --build build -j4` completed.
- `ctest --test-dir build --output-on-failure` passed both tests.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_EXEC=1 ./build/kls_smoke`
  passed.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_EXEC=1` on the top-ten
  CKTSO-gap focus produced
  `build/kls_plan_strict_cached_panel2_gap10_t4_r1_ref3_timeout120.jsonl`
  with a `2.29734s` geomean and no failed matrices.

The strict cached-panel shape is far too narrow to close the gap. The retained
plan still contains `750,247` runs over `6,285,128` rows, and current-column
batches cover `6,106,008` rows. The plan-derived producer starts expose
`39,658` possible panel objects over only `113,830` panel rows, and the
existing common-trailing cached-panel predicate accepts just `2,798` of them
over `7,168` rows with `14,203` trailing value entries. The benchmark still
reports zero actual cached panels, zero cached-probe attempts, and plan
execution disabled on nine of ten matrices. Therefore the next implementation
should not try to stretch the strict cached-panel object. It needs a ragged
retained-plan producer panel or current-column grouped accumulator that can
consume non-common trailing patterns directly.

The next retained-plan diagnostic measures that grouped-accumulator direction
directly. With
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_OUTPUT_STATS=1`, KLS now builds the
retained consumer plan, scans each current column's planned producer runs, and
counts output entries at or after the current pivot against the unique local
rows those entries would touch. The counters are reported as
`refactor_supernode_consumer_plan_deferred_*` for all current columns and
`refactor_supernode_consumer_plan_batch_deferred_*` for the current-column
batches that pass the existing ragged-supernode row gate. This is deliberately a
diagnostic scan, not production execution.

Validation:

- `cmake --build build -j4` completed.
- `ctest --test-dir build --output-on-failure` passed both tests.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_OUTPUT_STATS=1 ./build/kls_smoke`
  passed.
- A no-diagnostic same-binary top-ten CKTSO-gap baseline produced
  `build/kls_default_after_outputdiag_gap10_t4_r1_ref3_timeout120.jsonl` with a
  `4.63823s` geomean and no failed matrices.
- The output-stats diagnostic produced
  `build/kls_plan_output_defer_stats_gap10_t4_r1_ref3_timeout120.jsonl` with a
  `6.99864s` geomean and no failed matrices. The slower time is instrumentation
  cost from scanning retained-plan output structure.

This rules out the small-BLAS explanation for the current default gap. The
baseline run was built without CBLAS support and reported zero external CBLAS
updates and zero compact-supernode GEMV activity. The diagnostic run likewise
reported zero CBLAS updates. Across the top-ten focus, the retained plan still
covered `6,285,128` run rows and `6,106,008` current-column batch rows, but the
deferred output side exposed `875,343,104` raw entries collapsing to only
`6,257,172` unique output rows. The batchable subset was almost the same:
`873,846,433` raw entries collapsed to `5,748,636` unique rows. Per-matrix
collapse ratios ranged from about `21.9x` on `rajat03` to about `197.5x` on
`onetone1`, with the large ASIC and Rajat rows mostly above `100x`.

The direct paper-aligned next step is therefore a current-column grouped output
accumulator or ragged retained-plan producer-panel executor. It should consume
all retained runs for one current column, accumulate into a sparse/dense touched
row set once, and write each output row once before the pivot/store step. A
BLAS-size guard remains correct policy for optional external CBLAS, but it is
not the missing large mechanism on these slow cases.

That grouped-output direction now has retained shape-target substrate behind
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_SHAPE_TARGETS=1`. Instead of only
counting current-column deferred rows, KLS stores a compact target-row map for
each reusable producer shape and reports raw per-current publish entries versus
retained target rows. Correctness passed `cmake --build build -j2`,
`ctest --test-dir build --output-on-failure`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_SHAPE_TARGETS=1 ./build/kls_smoke`.
Focused probes completed on `ASIC_100ks` and `ASIC_320k`: `ASIC_100ks` reported
`8,434` target groups, `152,496` target runs, `133,615,959` raw entries,
`159,004` target rows, and max `611` rows/group; `ASIC_320k` reported `4,521`
target groups, `68,258` target runs, `109,078,358` raw entries, `75,161`
target rows, and max `508` rows/group. This is still not the executor, but it
is the retained sparse publish surface the CKTSO/SubtreeLU-style grouped
executor needs.

The follow-up execution prototype rejected the scalar grouped-accumulator half
of that direction. An opt-in current-column retained-plan accumulator was
implemented first as a full sparse output accumulator and then narrowed to only
accumulate the pivot/output value. The smoke suite stayed residual-clean, and
the top-three probes measured `2.83539s` for the full sparse accumulator in
`build/kls_plan_accum_gap3_t4_r1_ref3_timeout120.jsonl`, `2.80389s` after a
column-batch gate in
`build/kls_plan_accum_batchgate_gap3_t4_r1_ref3_timeout120.jsonl`, and
`2.38787s` for the pivot-only variant in
`build/kls_plan_accum_pivot_gap3_t4_r1_ref3_timeout120.jsonl`. That apparent
win was against a stale baseline.

The decisive same-binary top-ten comparison was worse on every row. The default
run in `build/kls_plan_accum_default_gap10_t4_r1_ref3_timeout120.jsonl` measured
`2.28458s` geomean, while the pivot-only accumulator in
`build/kls_plan_accum_pivot_gap10_t4_r1_ref3_timeout120.jsonl` measured
`3.42396s`, a `1.499x` regression. Per-matrix ratios ranged from `1.017x` on
`gemat12` to `2.029x` on `onetone1`; the large ASIC rows regressed by `1.268x`
to `1.833x`. The accumulator artifact did execute large retained-plan work on
the hard rows, for example `ASIC_320ks` applied `617,135` plan rows and
`185,394,303` entries, but still reported zero CBLAS and zero compact-GEMV
counters. The default comparison likewise reported zero CBLAS counters.

The direct conclusion is that another "use BLAS only for large cases" guard is
not an actionable patch here: the optional CBLAS paths already require the
build flag, runtime gate, 512-scale shape checks, and minimum-work checks, and
the current slow default path does not enter CBLAS. It is also not enough to
wrap the existing scalar EGraph replay in a write-combining accumulator. The
paper gap now points more specifically to a real retained producer-panel
numeric object with row-major or panel-major values, grouped dependency scans,
and fewer scalar L-entry probes before output accumulation, rather than a
post-hoc accumulator around the current per-entry replay.

The next rerun inspected whether that retained producer-panel executor can be a
simple exact-shape batch first, rather than a fully ragged batch from day one.
KLS now reports
`refactor_supernode_consumer_plan_shape_batch_*` counters from the retained
consumer plan. For each producer panel it groups planned runs by exact
`(panel_offset, run_rows)` and counts only groups with at least two consumers.
These counters are diagnostic; the executor is still the existing scalar EGraph
path.

Validation:

- `cmake --build build -j4` completed.
- `ctest --test-dir build --output-on-failure` passed both tests.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN=1 ./build/kls_smoke` passed.
- The current top-five auto-order rerun with the diagnostic produced
  `build/kls_plan_shape_batch_gap5_t4_r1_ref3_timeout120.jsonl` with a
  `7.74922s` geomean and no failures. The same-source auto-order control was
  `7.58543s`; CKTSO on the same five was `2.93787s`.
- The broader top-ten manifest produced
  `build/kls_plan_shape_batch_gap10_t4_r1_ref3_timeout120.jsonl` with a
  `4.03514s` geomean and no failures.

The BLAS-small-case hypothesis remains rejected by the new artifacts. The
top-ten run was built without CBLAS and reported zero external CBLAS updates
and zero compact-supernode GEMV activity on every row. The top-five run likewise
reported zero CBLAS and compact GEMV counters. Therefore an additional "BLAS
only for large cases" guard would not change these measured paths.

The exact-shape batching signal is strong. Across the top-ten manifest, the
retained plan covered `11,689,643` run rows, the existing broad batch gate
covered `11,619,009` rows, and exact `(panel_offset, run_rows)` groups covered
`10,069,373` rows in `45,742` groups and `1,016,724` runs. Among matrices with
nonzero retained-plan work, exact-shape coverage ranged from `58.4%`
(`transient`) to `91.2%` (`ASIC_100ks`) of planned rows, with an `80.9%`
geometric-mean coverage. The exact-shape groups still include many short runs:
`2,370,852` exact-shape rows were below the current ragged-supernode row gate,
so the production executor should batch large exact shapes first and keep a
cheap scalar/ragged fallback for small shapes.

This gives a clearer paper-aligned next implementation target than another
threshold tweak: build a retained producer-panel numeric object for exact-shape
consumer groups, apply all consumers of one producer panel through one grouped
scan/update, and leave unmatched or short groups on the current scalar path.
That directly attacks the missing coarse panel reuse while preserving the
already-correct BLAS size guards.

KLS now retains those exact-shape groups as internal plan metadata rather than
only counting them. The retained consumer plan stores, for every exact
`(panel_start, panel_offset, run_rows)` group with at least two consumers, a
group pointer, the group run ids, and the producer shape fields. The grouping
uses the same definition as the diagnostic counters above and is built only
when the retained consumer plan is requested. This is still not the grouped
numeric executor; it is the executor-ready index that removes the need to
rediscover shape batches from the run stream.

Validation:

- `cmake --build build -j4` completed.
- `ctest --test-dir build --output-on-failure` passed both tests.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN=1 ./build/kls_smoke` passed.
- `KLS_ENABLE_REFACTOR_U_SUPERNODE_PLAN_PATTERN=1
  KLS_ENABLE_REFACTOR_U_SUPERNODE_RAGGED_L=1 ./build/kls_smoke` passed.
- The top-five retained-plan run
  `build/kls_shape_groups_metadata_gap5_t4_r1_ref3_timeout120.jsonl` measured
  `7.58465s` geomean with no failures, effectively neutral against the
  same-source auto-order control at `7.58543s`.
- The broader top-ten retained-plan run
  `build/kls_shape_groups_metadata_gap10_t4_r1_ref3_timeout120.jsonl` measured
  `4.13033s` geomean with no failures. The previous shape-counter diagnostic
  was `4.03514s`; this run adds group metadata allocation and some timing noise,
  not a new numeric path.

The retained groups cover the same top-ten opportunity as the counter-only
probe: `45,742` exact-shape groups, `1,016,724` grouped runs, and
`10,069,373` grouped run rows out of `11,689,643` planned rows. The stored group
metadata increased retained-plan bytes by `9,597,608` across the top-ten
manifest. This keeps the memory cost modest relative to the matrix sizes and
makes the next implementation target concrete: a producer-panel grouped
executor can iterate `shape_group_ptr[group]..shape_group_ptr[group+1]`, gather
the current columns from the stored run ids, and apply one producer shape to
many consumers with a scalar/ragged fallback for unmatched groups.

The retained-plan U-supernode ragged-L experiment now prefilters dependency
probes through the retained consumer plan before touching the cached panel
state. The implementation records whether each current column's retained runs
are dependency-sorted, gives each EGraph worker a cursor for that sorted run
stream, and uses the retained `(panel_start, panel_offset, run_rows)` metadata
instead of rediscovering the planned run twice. This is still the scalar
ragged-L executor; it does not yet use the retained exact-shape groups as a
grouped numeric kernel.

Validation:

- `cmake --build build -j4` completed.
- `ctest --test-dir build --output-on-failure` passed both tests.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN=1 ./build/kls_smoke` passed.
- `KLS_ENABLE_REFACTOR_U_SUPERNODE_PLAN_PATTERN=1
  KLS_ENABLE_REFACTOR_U_SUPERNODE_RAGGED_L=1 ./build/kls_smoke` passed.
- On `onetone2` with one refactor repeat, the prefiltered path reported
  `50,274` U-supernode L probe attempts, `1,885` panel misses, `48,146`
  work rejects, `26,200` applied rows, and `4.37e-16` relative residual. The
  earlier unfiltered retained-pattern ragged-L probe attempted `528,565`
  dependencies on the same matrix, so this removes roughly `90.5%` of raw
  scalar probes before panel/cache checks.
- The final top-five CKTSO-gap focus run
  `build/kls_planpattern_ragged_prefilter_gap5_t4_r1_ref3_timeout120.jsonl`
  measured `8.50835s` geomean with no failures. The same-code default rerun
  `build/kls_current_auto_rerun_gap5_t4_r1_ref3_timeout120.jsonl` measured
  `7.56138s` geomean. Both runs were built without CBLAS and reported zero
  CBLAS updates.

This confirms that retained-plan filtering is useful but insufficient. The
experimental ragged-L path improved over the previous plan-pattern ragged-L
top-five run (`9.01978s` geomean), but it remains slower than the default KLU
column replay. On repeated refactors, KLS prunes U-supernode L panels that were
not actually applied; later iterations therefore turn most retained-run probes
into cheap panel misses instead of small-work rejects. That behavior reinforces
the same conclusion as the output-collapse and exact-shape diagnostics: the
missing paper-aligned mechanism is not a smaller BLAS threshold. Optional CBLAS
is already behind build/runtime gates plus 512-scale and minimum-work guards,
and it is absent from these measurements. The next plausible gap closer is the
grouped producer-panel executor over the retained exact-shape groups, so one
producer scan feeds many consumers without one scalar L-entry probe per planned
dependency.

The cached retained-plan executor now shares the same sorted current-column run
cursor used by the ragged-L retained-pattern path. The worker cursor resets
when a current column changes, when the stored cursor leaves the current
column's run slice, or when dependency order moves backward, and it can return
the matched run id for a future grouped executor. This removes the old linear
scan through one current column's retained runs for every planned cached-panel
probe, but it does not change the cached-panel numeric shape.

Validation:

- `cmake --build build -j4` completed.
- `ctest --test-dir build --output-on-failure` passed both tests.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_EXEC=1 ./build/kls_smoke`
  passed.
- `KLS_ENABLE_REFACTOR_U_SUPERNODE_PLAN_PATTERN=1
  KLS_ENABLE_REFACTOR_U_SUPERNODE_RAGGED_L=1 ./build/kls_smoke` passed.
- The top-five plan-exec focus improved only slightly, from `8.13022s` in
  `build/kls_plan_exec_precursor_gap5_t4_r1_ref3_timeout120.jsonl` to
  `8.10120s` in
  `build/kls_plan_exec_cursor_gap5_t4_r1_ref3_timeout120.jsonl`. The same-code
  default run
  `build/kls_default_after_plan_cursor_gap5_t4_r1_ref3_timeout120.jsonl`
  measured `7.67658s`, so the retained cached-panel executor remains slower
  than default.
- One-refactor samples show the narrow-coverage cause directly: `ASIC_320ks`
  made `3,784` retained-plan attempts, hit `317` plan runs, applied `299`
  cached-panel runs over `48,841` rows, and then disabled the executor;
  `onetone2` made `7,280` attempts, hit `255`, applied `235` over `26,108`
  rows, and also disabled the executor. Residuals remained clean.

This confirms that lookup overhead is not the large missing paper mechanism.
The executor still materializes only the strict cached-panel/common-trailing
shape, so most retained rows never become executable numeric work. The next
gap-closing implementation remains the broader materialized producer-panel or
current-column grouped accumulator that can consume retained plan runs whose
trailing patterns are ragged rather than common.

A fresh CBLAS size-guard rerun on the current source confirms that the
"small BLAS calls" hypothesis is already covered by the implementation and is
not the active CKTSO-gap cause. The optional CBLAS call sites require the CBLAS
build, `KLS_ENABLE_CBLAS_SUPERNODE=1`, 512-scale vector/panel shape checks
(`2048` rows on the row-first cached path), and minimum work/copy thresholds
before calling `dtrsv`, `dgemv`, `dtrsm`, or `dgemm`.

Validation:

- `cmake --build build-cblas -j4` completed.
- `ctest --test-dir build-cblas --output-on-failure` passed both tests.
- With `OPENBLAS_NUM_THREADS=1` and `KLS_ENABLE_CBLAS_SUPERNODE=0`, the
  top-five CKTSO-gap focus
  `build-cblas/kls_cblas_sizeguard_current_off_gap5_t4_r1_ref3_timeout120.jsonl`
  measured `1.38984s` geomean.
- With the same binary and `KLS_ENABLE_CBLAS_SUPERNODE=1`,
  `build-cblas/kls_cblas_sizeguard_current_on_gap5_t4_r1_ref3_timeout120.jsonl`
  measured `1.38296s` geomean, a `0.995x` on/off ratio.
- Every row in both artifacts reported zero
  `refactor_supernode_cblas_update_*`,
  `row_refactor_compact_supernode_gemv_*`, and
  `row_refactor_compact_supernode_trsv_*` counters.

Therefore adding another "use BLAS only for large cases" switch would be a
no-op for the focused slow rows. The current gap is still in the missing
paper-style coarse producer/consumer numeric executor, not in dispatching small
external BLAS kernels.

KLS now also reports whether those exact retained producer-panel groups exist
inside a single current column. The new
`refactor_supernode_consumer_plan_column_shape_batch_*` counters group each
current column's retained runs by exact
`(panel_start, panel_offset, run_rows)` and count only groups with at least two
runs. This deliberately includes `panel_start`, because the question is whether
one active EGraph `x` workspace can reuse one producer-panel numeric object
without a cross-current scheduler.

Validation:

- `cmake --build build -j4` completed.
- `ctest --test-dir build --output-on-failure` passed both tests.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN=1 ./build/kls_smoke` passed.
- `KLS_ENABLE_REFACTOR_U_SUPERNODE_PLAN_PATTERN=1
  KLS_ENABLE_REFACTOR_U_SUPERNODE_RAGGED_L=1 ./build/kls_smoke` passed.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN=1` on the top-ten CKTSO-gap
  focus produced
  `build/kls_column_shape_batch_gap10_t4_r1_ref3_timeout120.jsonl` with a
  `2.13626s` geomean and no failed matrices.

The result closes the current-column exact-shape branch. Across the top-ten
focus, the retained plan covered `6,285,128` rows, broad current-column batches
covered `6,106,008` rows, and producer-side exact-shape groups still covered
`5,394,956` rows in `39,077` groups. The new current-column exact producer
shape counters were zero on every row: zero groups, zero runs, zero rows, and
zero small-run rows. That is consistent with the plan structure: within one
current column, an exact `(panel_start, panel_offset)` identifies one dependency
run, so repeated exact producer-panel shape is fundamentally cross-current.

The next executor should therefore stop trying to stay inside the current
worker's single `x` workspace for the exact-shape path. To exploit the retained
groups already stored in `shape_group_ptr`, KLS needs a producer-panel grouped
task that can gather several current columns, apply one materialized producer
panel to multiple current workspaces, and then publish/merge the resulting
updates before those current columns complete. A current-column executor may
still help for ragged output accumulation, but it cannot consume the retained
exact producer-panel groups measured above.

The row-refactor work model now treats dense row groups as a single
supernode-panel update rather than counting their internal dependencies twice.
For `KLS_ROW_REFACTOR_GROUP_DENSE`, `kls_row_refactor_compute_group_work`
keeps input, U storage, and external dependency work, skips scalar dependency
charges for prior rows inside the same dense group, and then adds the dense
panel update cost once. This matches the SubtreeLU/CKTSO paper intent more
closely: the dense supernode update replaces the internal row-by-row walk.

Validation on the top-five CKTSO-gap focus:

- `cmake --build build -j4` completed.
- `ctest --test-dir build --output-on-failure` passed both tests.
- The default/env run
  `build/kls_dense_group_work_defaultenv_gap5_t4_r1_ref3_timeout120.jsonl`
  measured `1.38939s` geomean and kept row refactor off on all five matrices.
  The auto model still recommended the ASIC and `rajat03` cases, but the
  pre-pattern lower-bound gate rejected them.
- The forced row run
  `build/kls_dense_group_work_forcedrow_gap5_t4_r1_ref3_timeout120.jsonl`
  measured `5.09228s` geomean. It remained much slower than the column path,
  but its reported row work dropped as expected: `ASIC_320ks` from about
  `1.188x` to `1.011x` of dependency work, `ASIC_320k` from about `1.083x` to
  `1.009x`, and `ASIC_100ks` from about `1.085x` to `1.008x`.

This is a useful correction but not a gap closer. The forced row path now looks
nearly equal in modeled arithmetic on the ASIC cases while still running
roughly `3.5x` to `6.8x` slower in refactor time. That points back to staging,
workspace traffic, ready-queue overhead, and missing coarse producer/consumer
execution, not to small BLAS calls or scalar overcounting in the model.

KLS now also measures whether the retained exact producer-panel groups are
simple enough for a batch-start executor. The new
`refactor_supernode_consumer_plan_first_dep_shape_batch_*` counters keep only
exact-shape retained runs whose producer panel is the first actual U dependency
of the current column. That subset could be executed from freshly scattered
current-column workspaces without first replaying earlier dependencies.

Validation:

- `cmake --build build -j4` completed.
- `ctest --test-dir build --output-on-failure` passed both tests.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN=1 ./build/kls_smoke` passed.
- The top-ten retained-plan focus
  `build/kls_first_dep_shape_batch_gap10_t4_r1_ref3_timeout120.jsonl`
  completed with no failed matrices and measured `2.20643s` geomean.

The result rules out the smaller batch-start shortcut. The same run reported
`39,077` exact-shape retained groups covering `713,396` runs and `5,394,956`
rows, but every `first_dep_shape_batch_*` counter was zero on every matrix.
The next executor therefore must include the harder paper-style stage:
materialize several current-column workspaces, advance each one through its
earlier dependency stream to the retained `(panel_start, panel_offset)`, apply
the shared producer panel across that gathered batch, and then finish or
publish those current columns without double-processing them in the normal
EGraph scheduler.

A current CBLAS-enabled build confirms that small BLAS calls are already
guarded out of the focused slow paths. The CBLAS call sites require explicit
`KLS_ENABLE_CBLAS_SUPERNODE=1` plus minimum row/batch/panel sizes and estimated
work before issuing `dtrsv`, `dgemv`, `dtrsm`, or `dgemm`; otherwise they fall
back to the KLS scalar/compact kernels. With the refreshed `build-cblas`
binary, both
`build/kls_cblas_guard_current_off_gap10_t4_r1_ref3_timeout120.jsonl` and
`build/kls_cblas_guard_current_on_gap10_t4_r1_ref3_timeout120.jsonl`
completed the top-ten CKTSO-gap focus with no failures. The runtime-gate-on
run reported zero `refactor_supernode_cblas_update_*` counters on every row;
its cycle geomean was `2.52445s` versus `2.64174s` with the runtime gate off.
The follow-up `Freescale/transient` checks also reported zero CBLAS update
counters, with `4.31849s` cycle time on and `4.45335s` off. The remaining gap
therefore is not unguarded small BLAS dispatch. It remains the larger
paper-level executor problem: KLS must avoid repeatedly streaming long
producer rows and instead apply retained producer panels to gathered current
workspaces.

KLS now measures the missing advance-to-offset stage for those retained
producer-panel groups. The new
`refactor_supernode_consumer_plan_shape_batch_advance_*` counters scan each
exact-shape retained group, locate the producer dependency in the current
column's numeric U stream, and estimate the earlier dependency work needed to
advance a gathered current-column workspace to that retained
`(panel_start,panel_offset)` before applying the shared producer panel.

Validation:

- `cmake --build build -j4` completed with only the pre-existing long JSON
  format-string warning in `bench/kls_bench.c`.
- `ctest --test-dir build --output-on-failure` passed both tests.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN=1 ./build/kls_smoke` passed.
- The top-ten CKTSO-gap retained-plan focus
  `build/kls_shape_batch_advance_gap10_t4_r1_ref3_timeout120.jsonl`
  completed with no failed matrices and measured `2.23896s` geomean.

The counters now cover the same exact-shape retained groups as the existing
shape-batch metadata: `39,077` groups, `713,396` runs, and `5,394,956` retained
run rows across the top-ten focus. Advancing those gathered current workspaces
to the retained offsets requires `312,069,123` earlier U dependencies and about
`2.55859e10` estimated scalar update work, averaging about `437` earlier
dependencies and `35,865` scalar-work units per retained run. The largest
single retained run is preceded by `32,090` dependencies and about `678k`
estimated scalar-work units.

This rules out a lightweight panel-only executor as a CKTSO-gap closer. The
paper-shaped grouped executor must claim a batch of current columns, scatter
their inputs into several workspaces, advance each workspace through its earlier
dependency stream, apply the shared retained producer panel, and then finish or
publish those current columns so the normal EGraph loop does not process them a
second time. If the advance work is added on top of the existing per-column
factor loop, it will be too expensive; it has to replace that part of the
normal numeric path for the claimed columns.

KLS now has the scheduler substrate needed for that replacement path. The
opt-in `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_CLAIMS=1` path allocates a
per-column claim-generation map next to the existing pipeline-done map. Normal
EGraph workers check the map before factoring each column; if a future grouped
retained-plan producer claims a column, the worker waits for that column's
pipeline-done generation and skips the duplicate scalar factorization. Ready
queue mode also publishes successors after a claimed column finishes, so the
dependency graph can continue without refactoring the column twice.

This chunk deliberately does not claim columns yet. It only makes the normal
scheduler safe for a grouped producer to own columns in a later patch and
reports claim/skip/wait counters through `kls_stats` and `kls_bench`.

Validation:

- `cmake --build build -j4` completed with only the pre-existing long JSON
  format-string warning in `bench/kls_bench.c`.
- `ctest --test-dir build --output-on-failure` passed both tests.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN=1 ./build/kls_smoke` passed.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN=1
  KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_CLAIMS=1 ./build/kls_smoke`
  passed.
- The top-ten retained-plan focus with claims enabled completed with no failed
  matrices in `build/kls_claim_substrate_gap10_t4_r1_ref3_timeout120.jsonl`
  and measured `2.27105s` geomean.

The benchmark reported zero claimed columns, zero claim skips, and zero claim
waits, as expected for a dormant substrate. The same run used the normal
non-CBLAS build (`build_has_cblas=false`) and reported zero CBLAS updates, so
the latest check also confirms that this slow-path work is not being explained
by small external BLAS calls. The remaining paper gap is still the grouped
numeric producer that actually claims retained-plan current columns, advances
their workspaces through earlier dependencies, applies the shared retained
producer panel, publishes the finished columns, and avoids the old scalar path
for those claimed columns.

KLS now also retains the executable U-stream location for every retained
consumer-plan run. Each run stores the producer dependency position and block
start in the current column's numeric stream when the retained plan is built,
so the future grouped producer does not need to rediscover those offsets while
holding a batch of gathered current-column workspaces. The new
`refactor_supernode_consumer_plan_positioned_*` counters report how much of the
retained plan has a reusable location.

Validation:

- `cmake --build build -j4` completed with only the pre-existing long JSON
  format-string warning in `bench/kls_bench.c`.
- `ctest --test-dir build --output-on-failure` passed both tests.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN=1 ./build/kls_smoke` passed.
- The top-ten retained-plan focus
  `build/kls_plan_positions_gap10_t4_r1_ref3_timeout120.jsonl` completed with
  no failed matrices and measured `2.28316s` geomean.

The retained-position coverage is complete on the focus set: `750,247` planned
runs and all `6,285,128` planned rows also reported positioned locations. The
same run used the normal non-CBLAS build (`build_has_cblas=false`) and reported
zero CBLAS update runs/rows. Therefore the user's proposed large-case BLAS
guard is already satisfied for the active code path and would be a no-op on
this benchmark. The remaining paper gap is still the grouped numeric producer:
use these retained positions to advance claimed current-column workspaces once,
apply the shared retained producer panel across the batch, and skip the old
per-column scalar replay for those claimed columns.

The retained-position arrays are now consumed by the opt-in plan executor. When
the builder verifies that the stored U positions are monotone within each
current column, the retained-plan cursor keys directly on the current U-stream
position and still checks that the dependency column matches before accepting a
run. If a future matrix lacks complete or monotone retained positions, the
executor falls back to the previous dependency-key lookup.

Validation:

- `cmake --build build -j4` completed.
- `ctest --test-dir build --output-on-failure` passed both tests.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_EXEC=1 ./build/kls_smoke`
  passed.
- `KLS_ENABLE_REFACTOR_U_SUPERNODE_PLAN_PATTERN=1
  KLS_ENABLE_REFACTOR_U_SUPERNODE_RAGGED_L=1 ./build/kls_smoke` passed.
- The top-ten retained-plan executor run
  `build/kls_planexec_poscursor_gap10_t4_r1_ref3_timeout120.jsonl` completed
  with no failed matrices and measured `2.29693s` geomean.
- The same-source default control
  `build/kls_poscursor_default_gap10_t4_r1_ref3_timeout120.jsonl` completed
  with no failed matrices and measured `2.12339s` geomean.
- Same-source cached-panel and ragged-L toggles remained slower:
  `build/kls_poscursor_cached_gap10_t4_r1_ref3_timeout120.jsonl` measured
  `2.45566s`, and
  `build/kls_poscursor_ragged_gap10_t4_r1_ref3_timeout120.jsonl` measured
  `3.07614s`.

This improves the opt-in plan-executor mechanics but does not close the CKTSO
gap. The plan executor applied only `179,892` rows out of `6,285,128` planned
rows and disabled itself on 9 of the 10 matrices. The same-source default KLS
run is still `2.2648x` slower than the saved CKTSO top-ten subset
(`2.12339s` versus `0.937565s` geomean), and all runs above reported zero CBLAS
updates. The next required paper-level step is still broader grouped numeric
execution over retained producer panels, not more retained-run lookup tuning.

KLS now has an opt-in retained exact-shape producer L-cache behind
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_GROUP_CACHE=1`. This is the first
direct storage step for the paper-aligned grouped producer-panel executor:
instead of using the strict common-trailing cached-panel predicate, it validates
each retained exact-shape group as a dense internal L subpanel plus ragged
per-row trailing rows. The cache records structural coverage, dense slots,
trailing row slots, covered retained runs/run rows, and a conservative storage
cap result. It does not yet publish values or execute the grouped update, so
timing with the gate enabled measures staging overhead rather than the final
algorithm.

Validation:

- `cmake --build build -j4` completed with only the existing long JSON
  format-string warning in `bench/kls_bench.c`.
- `ctest --test-dir build --output-on-failure` passed both tests.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_GROUP_CACHE=1 ./build/kls_smoke`
  passed.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN=1
  KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_GROUP_CACHE=1 ./build/kls_smoke`
  passed.
- The top-ten opt-in cache run
  `build/kls_group_l_cache_gap10_t4_r1_ref3_timeout120.jsonl` completed with
  no failed matrices and measured `2.32137s` geomean.
- The same-source default control
  `build/kls_group_l_cache_default_gap10_t4_r1_ref3_timeout120.jsonl`
  completed with no failed matrices and measured `2.17693s` geomean.

The new cache validates `38,782` of `39,077` retained exact-shape groups,
covering `704,933` of `713,396` grouped runs and `5,314,536` of `5,394,956`
grouped run rows. The structural object needs `112,662` unique group rows,
`1,523,436` dense slots, `6,341,846` ragged trailing slots, and
`115,230,013` bytes on the focus set; no matrix hit the storage cap. The
default control kept the new cache disabled (`group_l_built=0`) and again used
the normal non-CBLAS build with zero CBLAS calls. This strengthens the current
gap diagnosis: a high-coverage retained producer-panel object is feasible, but
KLS still needs the value-publish and grouped executor that consumes these
panels instead of replaying the scalar current-column path.

The retained exact-shape producer L-cache now has a numeric publication layer.
The cache builder also creates a producer-column map from each completed L
column to the retained exact-shape groups and local rows that contain it. Under
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_GROUP_CACHE=1`, the EGraph
refactor publishes the freshly computed L values into the cached dense internal
rows and ragged trailing rows at the same column-completion points used by the
existing U-supernode and cached-panel publishers. The cache remains opt-in and
is still not consumed by a grouped executor.

Validation:

- `cmake --build build -j4` completed. The benchmark executable now reports
  two long JSON format-string warnings in `bench/kls_bench.c`; they are
  non-fatal and come from the expanded diagnostic JSON fields.
- `ctest --test-dir build --output-on-failure` passed both tests.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_GROUP_CACHE=1 ./build/kls_smoke`
  passed.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN=1
  KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_GROUP_CACHE=1 ./build/kls_smoke`
  passed.
- The top-ten opt-in publish run
  `build/kls_group_l_publish_gap10_t4_r1_ref3_timeout120.jsonl` completed with
  no failed matrices and measured `2.26451s` geomean.
- The same-source default control
  `build/kls_group_l_publish_default_gap10_t4_r1_ref3_timeout120.jsonl`
  completed with no failed matrices and measured `2.25914s` geomean.

The publish run wrote `1,716,223` dense cached values and `18,534,182` ragged
trailing cached values across the top-ten focus, while retaining the same
structural coverage as the prior cache run: `38,782` cached exact-shape groups,
`704,933` covered grouped runs, and `5,314,536` covered grouped run rows. It
recorded `45` cumulative group invalidations across repeated refactors, all
left disabled for fallback instead of being used by future grouped execution.
The default control kept the cache disabled (`group_l_built=0`) with zero group
L writes and zero CBLAS calls. This turns the high-coverage retained producer
object from structural-only storage into numeric storage; the remaining
paper-level gap is now the executor that claims compatible current columns,
uses these published panels, and skips the old scalar replay for those claimed
columns.

The retained group-L cache now also has an opt-in single-current-column
executor gate behind
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_GROUP_EXEC=1`. The executor adds a
run-to-group map, uses the retained exact-shape dense/ragged L values when a
current U dependency stream exactly matches a retained run, and reports
separate structural executable-run and applied-update counters. It deliberately
reuses the existing native ragged-supernode work gate and adds a structural
minimum of `1024` executable retained runs before allocating/publishing numeric
group-L values for exec-only runs.

Validation:

- `cmake --build build -j4` completed, with the same non-fatal long JSON
  format-string warnings in `bench/kls_bench.c`.
- `ctest --test-dir build --output-on-failure` passed both tests.
- The current-build default top-ten control
  `build/kls_group_exec_default_current_gap10_t4_r1_ref3_timeout120.jsonl`
  completed with no failed matrices and measured `2.18820s` geomean.
- The guarded group-exec top-ten run
  `build/kls_group_exec_minguard_gap10_t4_r1_ref3_timeout120.jsonl` completed
  with no failed matrices and measured `2.22099s` geomean.

The important result is negative and algorithmic: the top-ten focus has large
exact-shape cache coverage, but the single-run executor finds no worthwhile
per-current updates. Only `onetone1` reports any structurally executable
single-run groups (`63` runs, `34,650` rows), below the amortization guard, and
all matrices report zero applied group-L update rows. This explains why the
ungated probe regressed badly: it performed retained-plan lookups without
skipping scalar replay. The paper-aligned missing piece is therefore not more
single-run lookup tuning; it is a true grouped/batched executor that claims and
updates many same-shape consumers together, amortizing the retained producer
panel over a group as in the paper algorithms.

KLS now distinguishes standalone retained group-L executor eligibility from
aggregate batch-candidate eligibility. The prior gate asked whether each
single current-column run was large enough to justify a retained exact-L
update. That is the wrong diagnostic for the paper-shaped grouped executor:
the paper opportunity is many small same-shape consumers sharing one retained
producer panel. The group-L cache therefore now reports
`refactor_supernode_consumer_plan_group_l_batch_candidate_*` counters for
exact-L-valid groups with at least two matching retained runs. These counters
are explicitly named as candidates; they do not claim that a batched executor
has been implemented yet.

Validation:

- `cmake --build build -j4` completed, with the same non-fatal long JSON
  format-string warnings in `bench/kls_bench.c`.
- `ctest --test-dir build --output-on-failure` passed both tests.
- The opt-in group-L cache/executor metadata run
  `build/kls_group_l_batch_candidates_gap10_t4_r1_ref3_timeout120.jsonl`
  completed the top-ten CKTSO-gap focus with no failed matrices and measured
  `2.24087s` geomean.
- The same-source default control
  `build/kls_group_l_batch_candidates_default_gap10_t4_r1_ref3_timeout120.jsonl`
  completed with no failed matrices and measured `2.25333s` geomean; because
  the retained group-L cache was not requested, it reported zero batch
  candidate runs.

The retained group-L batch-candidate diagnostic now also measures whether a
candidate has enough grouped update work to pay for scalar prefix advancement
to the retained producer panel. With
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_GROUP_CACHE=1`, the current
top-ten CKTSO-gap run
`build/kls_group_l_batch_payoff_gap10_t4_r1_ref3_timeout120.jsonl` completed
with no failed matrices and measured `2.58472s` geomean. The same-source
default control
`build/kls_group_l_batch_payoff_default_gap10_t4_r1_ref3_timeout120.jsonl`
completed with no failed matrices and measured `2.19009s` geomean, with zero
group-L cache counters because the cache was not requested.

Across the cache-enabled top ten, all `38,782` batch candidates and `704,933`
candidate runs had computable advance costs, covering `1,454,726,482`
candidate update entries. The scalar prefix work to reach those panels was
about `25,156,176,333` estimated work units. Only `1,830` candidates,
`5,997` runs, and `3,141,804` update entries passed the simple payoff filter
where grouped update entries are at least the scalar advance work; their
advance work was only about `2,118,881`. On the hard ASIC rows, the payoff
surface was especially small: `ASIC_320ks` kept `288,654` of `211,576,799`
candidate entries, and `ASIC_320k` kept `22,566` of `239,333,804` candidate
entries.

This directly addresses the small-BLAS-guard hypothesis. Optional external
CBLAS is already gated by large size/work thresholds and these default focused
runs do not enter the external CBLAS counters. The retained-panel loss is
dominated by the cost of advancing current columns to the shared producer
panel, not by BLAS dispatch granularity after reaching the panel. A plausible
paper-aligned next executor must therefore batch or otherwise reduce the
advance stage itself before applying the retained L panels; merely replaying
the existing scalar prefix and adding another BLAS threshold cannot close the
ASIC gap.

The opt-in metadata run reports only `63` standalone executable group-L runs
and `34,650` standalone executable rows, but `38,782` aggregate batch-candidate
groups, `704,933` candidate runs, `5,314,536` candidate run rows, and
`1,454,726,482` candidate update entries, with a maximum of `766` runs in one
exact-shape group. This closes the diagnostic gap left by the single-run gate:
the retained exact-L storage has enough aggregate same-shape work to matter,
but KLS still needs the executor that owns multiple current columns and applies
one producer panel over the group instead of replaying those runs one column at
a time.

The next opt-in probe deliberately executed those aggregate candidates through
the existing single-current-column group-L executor behind
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_GROUP_BATCH_EXEC=1`. The cache
keeps retained exact-L values when aggregate batch candidates clear the
structural threshold, maps each retained run to its exact-L group, and relaxes
the standalone per-run work gate only for that explicit batch-exec experiment.
The hot run lookup now uses the same worker cursor as the other retained-plan
executors, so the probe is testing the algorithmic granularity rather than a
known linear lookup artifact.

Validation:

- `git diff --check` passed.
- `cmake --build build -j4` completed.
- `ctest --test-dir build --output-on-failure` passed both tests.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_GROUP_BATCH_EXEC=1 ./build/kls_smoke`
  passed.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_GROUP_EXEC=1 ./build/kls_smoke`
  passed.
- The top-ten CKTSO-gap focus
  `build/kls_group_l_batch_exec_cursor_gap10_t4_r1_ref3_timeout120.jsonl`
  completed with no failed matrices and measured `3.63452s` geomean. This is
  slightly better than the same aggregate executor before the cursor lookup
  cleanup (`3.76413s`), but still much slower than the same-source metadata-only
  cache run (`2.24087s`) and default control (`2.25333s`).

The aggregate executor consumed the full candidate surface: `38,782` candidate
groups, `704,933` candidate runs, and `5,314,536` candidate run rows became
`2,111,342` applied group-L update runs, `15,672,408` applied update rows, and
`4,273,553,183` applied update entries across the three repeated refactors.
That is useful negative evidence. KLS can now retain and replay the exact-L
groups, but replaying them one current column at a time is slower than the
normal path. The paper-level gap is therefore the true grouped producer-panel
executor: claim compatible current columns, advance each workspace once to the
retained producer offset, apply the shared exact-L panel over the group, and
skip the old scalar/current-column replay for those claimed updates. The same
run was built without CBLAS support (`build_has_cblas=false`) and reported zero
`refactor_supernode_cblas_update_*` calls, so an additional "use BLAS only for
large cases" guard cannot affect this result.

KLS now has default-off ownership infrastructure for that future grouped
executor behind
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_SHAPE_CLAIMS=1`: an inverse
`refactor_level_pos` map, distinct pipeline ownership tags for ordinary leased
columns versus shape-claimed columns, and a guarded shape-group owner probe
that only activates when the retained group-L run map is present. This is an
ownership substrate, not the final paper algorithm.

Validation:

- `git diff --check` passed.
- `cmake --build build -j4` completed.
- `ctest --test-dir build --output-on-failure` passed both tests.
- `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_SHAPE_CLAIMS=1 ./build/kls_smoke`
  passed.
- The default top-three control after the patch
  `build/kls_default_after_shape_claims_gap3_t4_r1_ref3_timeout120.jsonl`
  completed with no failed matrices and measured `2.05866s` geomean:
  `ASIC_320ks=10.8233s`, `ASIC_320k=12.6456s`, `gemat12=0.063746s`.
- The opt-in shape-claim probe
  `build/kls_shape_claims_gap3_t4_r1_ref3_timeout120.jsonl` timed out on both
  hard ASIC cases at the 120s per-matrix limit and only completed `gemat12` at
  `0.0647132s`.

This is stronger negative evidence against scalar future-column replay. The
ownership guard is safe enough for smoke tests and leaves default performance
unchanged, but using it to claim future columns and dispatch their complete
scalar EGraph refactors serializes the slow ASIC cases instead of closing the
gap. The missing CKTSO-paper-sized step remains the same and is now sharper:
once columns are safely owned, KLS must execute a true grouped producer-panel
kernel over multiple current workspaces. A BLAS size guard is still orthogonal:
these default and failed probe runs do not exercise external CBLAS update
counters, so small-case BLAS overhead is not the main cause of the ASIC loss.

A first BTF-unscaled gathered-workspace executor prototype was implemented and
rejected before commit. The prototype added
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_GROUP_L_CLAIM_EXEC=1`, claimed
future columns from one retained exact group, scattered each claimed current
column into a private block-sized workspace, advanced each workspace to the
retained producer offset through scalar dependencies, applied the shared
retained exact-L panel, then finished the scalar suffix and marked the claimed
columns done. It passed `cmake --build build -j4`,
`ctest --test-dir build --output-on-failure`, `git diff --check`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_GROUP_L_CLAIM_EXEC=1
./build/kls_smoke`, but the focused ASIC measurements were strongly negative.

Artifacts:

- Unguarded active grouped-claim run
  `build/kls_group_l_claim_exec_active_gap3_t4_r1_ref3_timeout120.jsonl`:
  `ASIC_320ks` and `ASIC_320k` both timed out at 120s; only `gemat12`
  completed (`0.0878407s`).
- With a 128-dependency prefix-advance guard,
  `build/kls_group_l_claim_exec_guard128_gap3_t4_r1_ref3_timeout120.jsonl`:
  the same two ASIC cases timed out at 120s; `gemat12` completed
  (`0.0729079s`).
- With a hard cap of one claimed column per refactor,
  `build/kls_group_l_claim_exec_cap1_gap3_t4_r1_ref3_timeout120.jsonl`
  completed but regressed badly: `ASIC_320ks=17.5253s`,
  `ASIC_320k=21.3527s`, `gemat12=0.0660925s`, geomean `2.91356s`.
  The two ASIC rows each claimed only one column and applied only one retained
  group-L update (`9` rows / `467` entries on `ASIC_320ks`, `7` rows / `329`
  entries on `ASIC_320k`), yet refactor time rose from the same-source default
  `0.081893s` to `0.144222s` and from `0.097155s` to `0.180003s`.

This rejects the naive gathered-workspace implementation, not the paper
algorithm itself. The expensive part is advancing and finishing claimed current
columns as isolated scalar columns with full block-sized workspace traffic.
Even a single claimed column can be an expensive current column, and the small
retained exact-L update it reaches does not amortize that cost. A viable
CKTSO-style executor must make the advance stage itself a batched/coarse task:
select current columns by estimated advance cost, advance multiple workspaces
without full block memset per column, apply substantial retained panels, and
merge/publish results without making other workers wait on one stolen heavy
column. The failed prototype was reverted from `src/kls.c`.

A safer follow-up keeps the shape-claim ownership probe from firing unless the
retained group-L executor is active and the retained group has measured
candidate payoff. The cache now stores a separate validity bit and payoff bit:
plain valid groups can still be retained for diagnostics or future executor
work, but `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_SHAPE_CLAIMS=1` no
longer claims columns for groups whose retained-L update work would not repay
the estimated prefix advance.

Validation:

- `build/kls_shape_claims_payoff_gap3_t4_r1_ref3_timeout120.jsonl` completed
  without failures but regressed the top-three probe: `ASIC_320ks=14.9006s`,
  `ASIC_320k=13.4147s`, `gemat12=0.0601335s`, geomean `2.29069s`.
  The two ASIC rows claimed future columns (`78` and `2` in the last refactor)
  while applying zero retained group-L updates, so the probe was still only
  rescheduling scalar current-column work.
- After requiring the retained update executor,
  `build/kls_shape_claims_payoff_execguard_gap3_t4_r1_ref3_timeout120.jsonl`
  completed with no failures and no shape claims: `ASIC_320ks=11.3446s`,
  `ASIC_320k=13.5294s`, `gemat12=0.0598688s`, geomean `2.09454s`.
  The same rows still reported payoff-positive retained groups
  (`288,654` and `22,566` retained update entries), but zero group-L update
  executions because the executor itself was not enabled.
- A CBLAS-capable top-ten control already has "use BLAS only for large cases"
  behavior: the source requires the runtime `KLS_ENABLE_CBLAS_SUPERNODE=1`
  gate plus 512/2048-scale row or panel dimensions and multi-million estimated
  work. `build/kls_cblas_guard_current_off_gap10_t4_r1_ref3_timeout120.jsonl`
  and `build/kls_cblas_guard_current_on_gap10_t4_r1_ref3_timeout120.jsonl`
  both report zero `refactor_supernode_cblas_update_*` counters on all
  top-ten CKTSO-gap rows, including `ASIC_320ks` and `ASIC_320k`.

This narrows the actionable gap again: another small-case BLAS guard would be
dead code for the current slow rows. The missing paper-aligned piece is still a
real grouped/coarse numeric executor that performs the prefix advance and
producer-panel application as a batch, not a standalone shape claim or a
different threshold around an inactive CBLAS call site.

The active-rank row pipeline now releases the pipeline mutex during the
commit-cursor dependency drain for large active-rank phases as well. To avoid
racing waiting workers that can still mutate the shared completed-panel cache,
the unlocked cursor uses a shallow workspace with `supernode_panel_cache=NULL`;
it can still consume read-only row-supernode metadata and compact/scalar
published-U rows, and any scratch reallocation is synchronized back to the
worker workspace after the drain. This is a direct concurrency fix for the
sampled `kls_row_first_partial_apply_one_dep` lock hold, not another BLAS
threshold.

Validation:

- `cmake --build build -j4` completed.
- `ctest --test-dir build --output-on-failure` passed both tests.
- `git diff --check` passed.
- `build/kls_asic320k_unlocked_active_rank_t4_factor.json` completed with
  `initial_factor_seconds=2.16228979` and `relative_residual_l2=2.61575448e-15`.
- `build/kls_rajat29_unlocked_active_rank_nofast_t4_factor.json` completed
  under forced KLS-first/no-fast with `initial_factor_seconds=5.13465759` and
  `relative_residual_l2=9.86042515e-12`, so it did not reproduce the
  correctness failures seen in earlier active-rank pivot/restart prototypes.
- `build/kls_pre2_unlocked_active_rank_trace70.stderr` reached the
  `pivot-tail` trace at row 274,430 within an 80s process cap. The previous
  same-source control trace `build/kls_pre2_current_trace60.stderr` only
  reached row 196,608 before timing out. The new trace still reported
  `15.379B` scalar U entries, with `10.364B` output/trailing entries.
- `build/kls_pre2_unlocked_active_rank_t4_factor_timeout130.json` remained
  empty after the 130s timeout, although stderr printed the dominant-BTF trace
  twice. The change improves first-pass concurrency but does not close the
  hard measured-factor timeout.
- The focused top-five CKTSO-gap reruns completed without failures:
  `build/kls_unlocked_active_rank_gap5_t4_r1_ref3_timeout120.jsonl` measured
  `1.41674462s` geomean, and
  `build/kls_unlocked_active_rank_gap5_rerun_t4_r1_ref3_timeout120.jsonl`
  measured `1.36983167s` geomean.

This retained change removes a real serialized section from the large
active-rank pipeline, but the trace confirms the remaining CKTSO gap is still
the paper-sized grouped producer/output numeric executor rather than a mutex
alone.

A follow-up wait-snapshot prototype tried to let workers waiting behind the
commit cursor drain ready active-rank dependencies against a snapshot of the
published U prefix. It preserved correctness on the sensitive forced
KLS-first probes, but was rejected before commit: the top-five CKTSO-gap
geomean regressed twice (`1.4503s` and `1.5076s` in
`build/kls_wait_snapshot_guard_gap5*_t4_r1_ref3_timeout120.jsonl`), and
`pre2` still reached the same pivot-tail point as the retained active-rank
commit-drain unlock (`274,430/629,628`) before the timeout. The retained
source therefore stays with the committed cursor-only unlock.

The CBLAS large-case guard was tightened independently. The optional
KLS-first cached-panel CBLAS path now requires at least 512 dense/tail update
columns in addition to the existing 2048 producer-row, 50M-work, and
work-per-copied-entry thresholds. This directly implements the conservative
"BLAS only for large useful updates" policy without changing the default
non-CBLAS build. Same-binary CBLAS-capable forced `ASIC_320k` checks completed
with clean residuals (`2.61575448e-15`): `KLS_ENABLE_CBLAS_SUPERNODE=0`
measured `2.5094s`, and `=1` measured `2.4994s`. The final non-CBLAS top-five
CKTSO-gap check also completed without failures at `1.4545s` geomean in
`build/kls_cblas_update_guard_final_gap5_t4_r1_ref3_timeout120.jsonl`. The
guard is retained as a safe policy bound rather than as evidence that BLAS
thresholding closes the paper gap.

The next retained cleanup splits the unscaled BTF EGraph column kernel into a
plain scalar loop when all optional supernode, group-L, ragged-L, and
U-supernode value hooks are disabled. The default CKTSO-gap path was paying
those disabled hook branches on every U-stream dependency even in non-CBLAS
runs. The hook-capable loop is unchanged for experimental modes, while the
plain path still uses the same exact EGraph schedule and KLU-compatible
column storage. Validation passed `cmake --build build -j4`,
`ctest --test-dir build --output-on-failure`, and `git diff --check`.
Focused top-five reruns measured `1.4181s` and `1.4180s` geomean in
`build/kls_plain_btf_branch_gap5*_t4_r1_ref3_timeout120.jsonl`; the top-ten
continuation comparison measured `2.1172s` geomean versus `2.1949s` for
`build/kls_current_gap10_continuation_t4_r1_ref3_timeout120.jsonl`. This is
kept as low-risk hot-path cleanup only. The remaining slow rows still run the
scalar EGraph executor, so the paper-sized gap remains the missing
row/segment-oriented grouped numeric engine and checked pivoting-tail
scheduler, not BLAS granularity or disabled-hook branch overhead.

A follow-up payoff-gating probe tested whether the existing single-current
group-L executor could be salvaged by consuming only groups whose retained
update entries exceed the modeled scalar prefix advance. It could not. The
first version executed only payoff-positive groups but still published exact-L
values for every valid retained group; the top-ten CKTSO-gap run regressed to
`4.2323s` geomean in
`build/kls_group_l_payoff_exec_gap10_t4_r1_ref3_timeout120.jsonl`. A second
version added an explicit standalone-executable flag and skipped producer
publish work for groups that were neither standalone-executable nor
payoff-positive under batch execution; it improved the failed probe to
`3.9022s` geomean in
`build/kls_group_l_payoff_publishskip_gap10_t4_r1_ref3_timeout120.jsonl`, but
that was still slower than the prior all-batch opt-in run (`3.6345s`) and far
slower than the current default (`2.1172s`). The source experiment was
reverted. The result is useful negative evidence: filtering a one-current-
column replay path cannot close the gap, even when the retained panel work is
payoff-positive. The needed paper-scale executor must own multiple compatible
current workspaces, advance them as a batch to the retained producer offset,
and apply the shared producer panel once across that group.

The retained consumer-plan diagnostics now also measure exact non-empty prefix
sharing inside each retained producer-shape group. The new
`refactor_supernode_consumer_plan_prefix_advance_batch_*` counters hash and
collision-check the actual earlier U-dependency stream before the retained
panel dependency, so they answer whether a direct "same prefix, then shared
producer panel" executor has enough structural surface. The top-ten CKTSO-gap
diagnostic with `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN=1` completed in
`build/kls_prefix_advance_diag_gap10_t4_r1_ref3_timeout120.jsonl`; the matching
group-cache diagnostic completed in
`build/kls_prefix_advance_groupcache_gap10_t4_r1_ref3_timeout120.jsonl`.
Across the ten rows, exact shape groups still represented `713,396` runs and
`25.586B` modeled prefix-advance work, but identical non-empty prefix batches
covered only `36,122` runs, `1.213M` dependencies, and `136.4M` work
(`0.53%` of advance work). The ASIC rows, which dominate the CKTSO gap, had
almost no coverage (`0.04%`, `0.01%`, and `0.05%` work share for
`ASIC_320ks`, `ASIC_320k`, and `ASIC_100ks`). `onetone2` and `onetone1` had
larger prefix sharing (`7.7%` and `2.3%` work share), but not enough to explain
the broad slow set. This rules out an exact identical-prefix batch executor as
the main gap closer. The next paper-aligned implementation should therefore
batch less-identical row/segment streams with shared producer panels and
coarse output accumulation, not another exact-shape or exact-prefix replay.

A direct Algorithm 5 payoff audit then separated whole producer panels from
payoff-positive subsets inside those panels. The first top-ten CKTSO-gap
diagnostic,
`build/kls_alg5_payoff_diag_gap10_t4_r1_ref3_timeout120.jsonl`, measured
`2.1783s` geomean and found only `3` whole panels where replaying all selected
producer work paid off. In contrast, selected subsets inside panels covered
`594` groups and `13,980` runs, with about `251.4M` modeled update work versus
`132.9M` modeled prefix-advance work. This is the clearest paper-aligned
surface so far: the executor should not replay whole producer panels, but it
should batch the payoff-positive current subsets that share a retained
producer panel.

The scalar scheduling shortcut was rejected. A claim-based prototype that tried
to steal those Algorithm 5 payoff columns without a new batched numeric kernel
regressed to `2.2933s` geomean in
`build/kls_alg5_payoff_claims_gap10_t4_r1_ref3_timeout120.jsonl`, and the
claim counters stayed at zero. That path was removed. The retained source now
keeps only an opt-in selected-run substrate behind
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_PLAN=1`; ordinary
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN=1` diagnostics still report the
same Algorithm 5 payoff counters without allocating dead execution-plan arrays.
Current reruns on the same source reported the same `594` groups and `13,980`
runs in
`build/kls_alg5_payoff_diag_rerun_gap10_t4_r1_ref3_timeout120.jsonl` and
`build/kls_alg5_payoff_plan_gated_gap10_t4_r1_ref3_timeout120.jsonl`; both were
non-CBLAS builds with zero CBLAS update counters. The gated retained-plan run
was slower (`3.9534s` versus `3.5523s` diagnostic-only), confirming that
retaining the selector is useful only as a substrate for the future batched
executor, not as a performance change by itself. This also keeps the latest
BLAS-size-guard hypothesis bounded: the focused slow paths are not dispatching
external BLAS at all, so another small-case BLAS threshold cannot close this
gap.

The next implementation tried the direct ragged producer-panel variant of that
Algorithm 5 subset. `KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_EXEC=1`
now builds the retained consumer plan, keeps the payoff-positive Algorithm 5
subset, materializes only those selected prefix extents as the existing
plan-derived ragged U-supernode L pattern, and clamps the ragged-L executor to
selected runs. Correctness passed `ctest --test-dir build --output-on-failure`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_EXEC=1 ./build/kls_smoke`,
and the same smoke with `KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN=1`.

The focused top-ten result rejects this as the missing CKTSO mechanism. The
new opt-in numeric path measured `4.9043s` geomean in
`build/kls_alg5_payoff_ragged_exec_gap10_t4_r1_ref3_timeout120.jsonl`, while
the same-source default control measured `2.1604s` in
`build/kls_alg5_payoff_ragged_exec_default_gap10_t4_r1_ref3_timeout120.jsonl`.
The run was still a normal non-CBLAS build with zero CBLAS update counters.
It built `594` Algorithm 5 payoff patterns over `11,650` producer rows and did
execute the ragged path (`567` updates, `228,450` rows, and `89.2M` update
entries), but seven rows disabled the ragged executor after low useful
coverage. The ASIC rows built payoff patterns but reported no cumulative
ragged updates before disable; `onetone1` executed `546` updates over `223,131`
rows and still remained slower. This means the payoff-positive ragged producer
panel alone is not enough. It removes some producer L scans, but it still
advances and executes one current column at a time. The remaining paper gap is
more specific: KLS needs a multi-current Algorithm 5 executor or persistent
consumer accumulator so one producer-prefix panel and one prefix-advance phase
feed many current workspaces, rather than another per-current replay path.

The follow-up multi-current diagnostic makes that conclusion sharper and moves
the next work away from BLAS guards. `build/kls_alg5_multicurrent_diag_gap10_t4_r1_ref3_timeout120.jsonl`
reported the same `594` payoff-positive Algorithm 5 producer-panel subsets and
`13,980` selected runs, but those runs also correspond to `13,980` distinct
current columns. Every payoff subset was multi-current, with max per-panel
current counts matching max selected-run counts: `196` on `ASIC_320ks`, `233`
on `ASIC_320k`, `153` on `onetone2`, and `322` on `onetone1`. That means the
paper-aligned opportunity is not repeated work inside one current column. The
missing executor must gather many current-column workspaces for one retained
producer prefix, advance them once to the selected offset, and apply/publish the
shared prefix update as a batch.

KLS now retains the corresponding batch-substrate map explicitly. The
Algorithm 5 payoff selector keeps a compact group pointer over payoff-positive
producer panels, the producer-panel start for each group, and a run-to-group
map for selected consumer-plan runs. The old ragged payoff executor now accepts
only grouped runs, so single payoff-positive runs are no longer treated as
executable batch work. `build/kls_alg5_groupmap_diag_gap10_t4_r1_ref3_timeout120.jsonl`
reproduced the same `594` groups, `13,980` selected runs/current columns,
`670,921` prefix rows, and modeled work `251.4M/132.9M` update/advance totals.
This still is not the final numeric batch executor, but the scheduler and
kernel no longer have to infer sibling current columns by rescanning the whole
consumer plan.

The retained group map now also carries the execution metadata needed by that
future batch kernel. Each selected run stores its exact Algorithm 5 prefix-row
count, each payoff group stores total and max prefix rows, and stats report how
many grouped runs already have retained dependency positions. The focused
artifact `build/kls_alg5_execmeta_plan_gap10_t4_r1_ref3_timeout120.jsonl`
verified that the grouped metadata exactly matches the selected payoff surface:
`594` groups, `13,980` selected runs/current columns, `670,921` retained
group-prefix rows, and `13,980` positioned grouped runs. The largest retained
per-run prefix has `546` rows. The remaining missing piece is therefore the
numeric multi-current workspace/publish kernel itself, not plan discovery,
prefix-length recovery, or dependency-position lookup.

The next executable slice moved the opt-in ragged Algorithm 5 path onto that
retained workspace layout. Under
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_EXEC=1`, a selected run now
validates its group/current slot, allocates the group's compact current
workspace extent, and executes through the mapped current-slot offset rather
than a standalone run-sized scratch array. Correctness passed
`cmake --build build -j4`, `ctest --test-dir build --output-on-failure`,
`./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_PLAN=1 ./build/kls_smoke`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_EXEC=1 ./build/kls_smoke`.

The focused top-ten CKTSO-gap check
`build/kls_alg5_mapped_workspace_exec_gap10_t4_r1_ref3_timeout120.jsonl`
completed without failures at `2.9601s` geomean, while the same-source default
control `build/kls_mapped_workspace_default_gap10_t4_r1_ref3_timeout120.jsonl`
measured `2.2298s`. The mapped path did execute numeric ragged updates on
`rajat25` (`7` updates, `1,773` rows, `363,412` entries) and `onetone1`
(`182` updates, `74,377` rows, `29.37M` entries), but the ASIC rows still only
populated the retained multi-current workspace map. This confirms the remaining
paper gap is not BLAS dispatch and not the workspace-addressing substrate; the
missing CKTSO-scale mechanism is the actual multi-current executor that gathers
and advances many current workspaces for one retained producer group before
publishing them.

KLS now retains Algorithm 5 prefix-advance counts inside the payoff descriptor
instead of only retaining current-slot workspace geometry. Each selected run
stores the number of producer dependencies that must be advanced before the
selected prefix update, and public stats report the total advance dependencies,
maximum per-run advance dependencies, and zero-advance selected runs across all
payoff groups. The focused top-five CKTSO-gap artifact
`build/kls_alg5_advance_descriptor_gap5_t4_r1_ref3_timeout120.jsonl` completed
with no failed passes at `1.391166153602165` geomean. It reports
`190,994` retained advance dependencies for `ASIC_320ks` with zero
zero-advance selected runs, and `194,219` retained advance dependencies for
`ASIC_320k` with one zero-advance selected run. `ASIC_100ks` still has
`32,849` retained advance dependencies and only five zero-advance selected
runs. This makes the next paper gap concrete: a first-dependency shortcut is
too small for the slow cases, and BLAS dispatch is not the first explanation to
pursue. The missing CKTSO/SubtreeLU mechanism is a grouped Algorithm 5 executor
that claims current workspaces, batches prefix advancement across those
workspaces, and only then applies and publishes the shared producer-prefix
update.

The next implementation added the first shared-work scheduler slice for that
descriptor. `KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_QUEUE=1` now builds
the retained Algorithm 5 payoff groups as a self-contained request, allocates a
bounded solver-owned queue, publishes current-column candidates when a
producer-prefix trigger column is marked done, and lets at most half the worker
threads claim queued columns on pop before dispatching them through the normal
EGraph dependency-checked column executor. This deliberately does not require
the older `KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_EXEC=1` ragged-L
numeric path; the queue uses the paper payoff descriptor as scheduling metadata.
Correctness passed `cmake --build build -j4`,
`ctest --test-dir build --output-on-failure`, `./build/kls_smoke`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_QUEUE=1 ./build/kls_smoke`.

The bounded claim-on-pop queue is safe but not yet the missing CKTSO-scale win.
Same-build 4-thread controls measured `0.0408858587s` refactor on
`ASIC_100ks`, `0.0802955404s` on `ASIC_320ks`, and `0.000656136405s` on
`rajat03`. The queue probes measured `0.077946977s`, `0.13615748s`, and
`0.000625623402s` respectively, with `13`, `29`, and `4` queued claimed columns
and matching residuals. A more aggressive claim-on-publish variant was also
tried because claim-on-pop loses most candidates to the ordinary pipeline, but
it hit the `40s`/`60s` timeouts on `ASIC_100ks`/`ASIC_320ks`; that variant was
rejected. The result sharpens the remaining gap: the paper mechanism is not
just earlier column ownership. KLS still needs the true multi-current numeric
batch that advances many current workspaces for one retained producer-prefix
group and publishes the shared update without converting those current columns
back into ordinary scalar EGraph column tasks.

A follow-on queue pass tested the safer half of claim-on-publish: preclaim a
candidate only when all U predecessors are already complete in the current
EGraph generation. Keeping not-ready candidates in the same queue still exposed
the old scheduler hazard: `ASIC_320ks` hit the `60s` timeout because pipeline
workers can wait on preclaimed columns before any worker drains the queue. The
kept implementation therefore makes the queue ready-only and lets a worker that
encounters a claimed column drain payoff-queue work while it waits. Correctness
passed `cmake --build build -j4`, `ctest --test-dir build --output-on-failure`,
and `KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_QUEUE=1 ./build/kls_smoke`.
Focused probes then completed:
`build/asic100ks_alg5_queue_readyclaim_probe.json` measured `0.0362660967s`
refactor with `0` queued claims, `build/asic320ks_alg5_queue_readyclaim_probe.json`
measured `0.0792330637s` with `2` queued claims, and
`build/rajat03_alg5_queue_readyclaim_probe.json` measured `0.00091396682s` with
`0` queued claims. This confirms the ready-only queue is a safe diagnostic but
not the CKTSO-closing mechanism; the useful paper gap remains the grouped
multi-current numeric executor, not more scalar column stealing.

KLS now records the retained Algorithm 5 target/accumulator surface needed by
that grouped executor. Each selected payoff run stores a target-entry count
covering dense suffix updates inside the retained producer pattern plus
best-effort L-trailing entries outside that pattern. Public stats and
`kls_bench` report the total target entries, maximum per group, and maximum per
run. The implementation deliberately leaves the existing scalar ragged-L
executor unchanged; these counters are descriptor substrate for the paper
algorithm, not a BLAS or scalar-scheduler tuning change. Correctness passed
`cmake --build build -j4`, `ctest --test-dir build --output-on-failure`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_PLAN=1 ./build/kls_smoke`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_EXEC=1 ./build/kls_smoke`. The
top-five CKTSO-gap descriptor probe
`build/kls_alg5_target_surface_gap5_t4_r1_ref3_timeout120.jsonl` completed with
no failures and `1.7929971541164775s` geomean. It reported large retained target
surfaces on the ASIC EGraph rows: `ASIC_320ks` has `129` groups, `3,596`
current workspaces, `157,067` prefix workspace rows, `46,839,997` target
entries, and `190,994` advance dependencies; `ASIC_320k` has `121` groups,
`3,755` current workspaces, `146,027` prefix workspace rows, `48,674,072`
target entries, and `194,219` advance dependencies; `ASIC_100ks` has `102`
groups, `1,112` current workspaces, `20,859` prefix workspace rows,
`4,668,124` target entries, and `32,849` advance dependencies. This makes the
next paper gap more concrete: the missing executor must retain and publish tens
of millions of accumulator updates for grouped current workspaces, not just
claim earlier scalar columns or adjust BLAS thresholds.

The retained target surface is now addressable rather than only counted. The
Algorithm 5 payoff descriptor stores a target offset for each selected run,
each group's retained producer-pattern width, and each group's target-entry
span, and the opt-in mapped path validates that a selected run's target slice
fits inside its group span before using the descriptor. This still leaves the
numeric executor unchanged, but it removes another rediscovery step from the
future grouped multi-current accumulator. Correctness passed
`cmake --build build -j4`, `ctest --test-dir build --output-on-failure`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_PLAN=1 ./build/kls_smoke`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_EXEC=1 ./build/kls_smoke`.
The top-five descriptor probe
`build/kls_alg5_addressable_target_gap5_t4_r1_ref3_timeout120.jsonl` completed
with no failures and `1.7738882990342122s` geomean. It reports retained
pattern-width sums and target entries on the ASIC EGraph rows:
`ASIC_320ks` has pattern-width sum/max `2,425/264` and `46,839,997` target
entries, `ASIC_320k` has `2,604/535` and `48,674,072`, and `ASIC_100ks` has
`1,535/157` and `4,668,124`. `gemat12` and `rajat03` used the mapped path and
therefore correctly report zero Algorithm 5 payoff groups. The next CKTSO-paper
gap remains the actual grouped numeric kernel: allocate current workspaces and
target accumulators from these retained offsets, batch prefix advancement, and
publish the accumulated dense-suffix/L-trailing updates.

KLS now makes the other half of that descriptor addressable as well: selected
Algorithm 5 payoff runs retain a prefix-advance offset, and each payoff group
retains the total advance-dependency span that owns those offsets. The scalar
opt-in path validates both target-span and advance-span ownership before it
accepts a selected run, so the future grouped executor can consume the same
descriptor without rediscovering either slice. This is still descriptor
substrate, not a BLAS or scalar-scheduler tuning change. Correctness passed
`cmake --build build -j2`, `ctest --test-dir build --output-on-failure`,
`./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_PLAN=1 ./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_EXEC=1 ./build/kls_smoke`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_QUEUE=1 ./build/kls_smoke`.
The top-five descriptor probe
`build/kls_alg5_advance_spans_gap5_t4_r1_ref3_timeout120.jsonl` completed with
no failures and `1.8396695264324108s` geomean. The hard ASIC EGraph rows still
show nonzero retained advance spans alongside the existing target spans:
`ASIC_320ks` has `129` groups, `190,994` advance deps, and `46,839,997` target
entries; `ASIC_320k` has `121`, `194,219`, and `48,674,072`; `ASIC_100ks` has
`102`, `32,849`, and `4,668,124`. `gemat12` and `rajat03` remain mapped-path
rows with zero Algorithm 5 payoff groups.

KLS now separates the retained Algorithm 5 target-entry work volume from the
exact accumulator slots a grouped executor would need to publish. Selected
payoff runs retain a target-slot offset, each group retains its target-slot
span, and `group_target_cols` stores those slots as compact local workspace row
indices. Dense suffix slots are retained directly; L-trailing rows are stamped
per selected run so duplicated rows become one accumulator/publish position for
that run. The existing scalar opt-in path validates target-entry, target-slot,
and advance-span ownership before accepting a selected run, but the numeric
kernel remains unchanged. This is paper-aligned descriptor work, not BLAS
threshold tuning.

Correctness passed `cmake --build build -j2`,
`ctest --test-dir build --output-on-failure`, `./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_PLAN=1 ./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_EXEC=1 ./build/kls_smoke`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_QUEUE=1 ./build/kls_smoke`.
The top-five CKTSO-gap probe
`build/kls_alg5_target_slots_gap5_t4_r1_ref3_timeout120.jsonl` completed with
no failures and `1.8057280532263165s` geomean. The hard ASIC EGraph rows show
why slot-indexed accumulation is the right next step: `ASIC_320ks` has
`46,839,997` target entries but only `784,673` target slots
(`67,651` max group, `492` max run), `ASIC_320k` has `48,674,072` entries but
`971,718` slots (`96,854` max group, `598` max run), and `ASIC_100ks` has
`4,668,124` entries but `188,230` slots (`22,332` max group, `615` max run).
The next CKTSO-paper gap remains the true grouped multi-current executor:
allocate current workspaces and slot-indexed accumulators from these retained
offsets, batch prefix advancement across the group, then publish the accumulated
dense-suffix and L-trailing updates without replaying each current column as a
scalar EGraph task.

The retained group descriptor now also owns the inverse map needed by that
executor: `group_current_run_ptr` and `group_current_runs` list the selected
Algorithm 5 payoff runs attached to each retained current slot. The scalar
opt-in path validates that its selected run is present in the slot's retained
run list before using the descriptor. This removes the last per-column
consumer-plan lookup a group-triggered executor would otherwise need when it
starts from a payoff group rather than from an individual current column. It
also handles the conservative case where a retained current owns more than one
selected run in the same group. Correctness passed `cmake --build build -j2`,
`ctest --test-dir build --output-on-failure`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_PLAN=1 ./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_EXEC=1 ./build/kls_smoke`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_QUEUE=1 ./build/kls_smoke`.
The focused top-five CKTSO-gap probe
`build/kls_alg5_current_runs_gap5_t4_r1_ref3_timeout120.jsonl` also completed
with no failures and `1.8909105332569978s` geomean. The hard ASIC EGraph rows
still retain the same large grouped surfaces: `ASIC_320ks` has `129` groups,
`3,596` currents, `784,673` target slots, and `190,994` advance dependencies;
`ASIC_320k` has `121`, `3,755`, `971,718`, and `194,219`; `ASIC_100ks` has
`102`, `1,112`, `188,230`, and `32,849`. This is still descriptor substrate,
not a BLAS or CPU-specific tuning change; the remaining paper gap is the
numeric grouped executor that uses the retained current-run map, current
workspace offsets, advance spans, and target slots together.

KLS now fills a narrower CKTSO Algorithm 5 execution semantic in the EGraph
numeric kernels under the modeled-work auto gate. A blocked pipeline column can
consume later already-finished scalar U predecessors before the current
predecessor is done, but only when a structural safety scan proves that no
earlier unapplied predecessor can still write the candidate workspace entry.
The kernel tracks consumed positions in an applied bitmap, marks ordinary
scalar and batched dependencies as they are consumed, and disables further
batched dependency runs for that column after the first true out-of-order
prefactor. This implements Algorithm 5's skip-unfinished prefactor idea without
turning the retained payoff descriptor into a misleading "backend" or a
CPU-specific micro-tuning knob.

Correctness passed `cmake --build build -j2`,
`ctest --test-dir build --output-on-failure`, `./build/kls_smoke`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_EXEC=1 ./build/kls_smoke`.
The top-five CKTSO-gap probe with the prefactor path
`build/kls_pref_update_only_gap5_t4_r1_ref3_timeout120.jsonl` completed with
no failures and `1.3938181393489573s` geomean versus the current no-new-flag
control `build/kls_current_default_gap5_t4_r1_ref3_timeout120.jsonl` at
`1.4204760684642577s`. The last-refactor prefactor counters confirm that the
path fired on the hard ASIC rows: `ASIC_320ks` used `61` columns / `125` deps,
`ASIC_320k` used `66` / `81`, and `ASIC_100ks` used `5` / `425`. Combining it
with the existing scalar ragged payoff executor remained worse:
`build/kls_alg5_pref_update_gap5_t4_r1_ref3_timeout120.jsonl` measured
`1.8658851843625461s`. The lesson is that the direct Algorithm 5 prefactor
semantic helps modestly, but the large CKTSO gap still points at the missing
grouped multi-current numeric executor rather than BLAS thresholds or the
current scalar payoff replay.

The Algorithm 5 prefactor flag is now opportunistic rather than scheduler
shaping. The earlier implementation allocated `pipeline_done` solely because
`KLS_ENABLE_EGRAPH_ALGORITHM5_PREF_UPDATE=1` was set, which made the experiment
measure both the skip-unfinished prefactor semantic and a changed EGraph
scheduler. KLS now only uses pipeline completion state that the selected
refactor path already needed for other reasons. Inside a BTF EGraph column, the
applied bitmap is allocated lazily only after a dependency actually blocks; when
allocated, it marks the already-consumed prefix as applied before scanning later
published predecessors. That prefix initialization is part of the algorithmic
guard, because without it the out-of-order scan can treat ordinary consumed
predecessors as still pending.

Correctness passed `cmake --build build -j2`,
`ctest --test-dir build --output-on-failure`, and `./build/kls_smoke`. The final
top-five CKTSO-gap probe
`build/kls_pref_update_lazy_prefix_gap5_t4_r1_ref3_timeout120.jsonl` completed
with no failures and `1.3868464599815646s` geomean. The 20-matrix CKTSO-gap
slice `build/kls_pref_update_lazy_prefix_gap20_t4_r1_ref3_timeout120.jsonl`
completed with no failures and `1.895666332219691s` geomean, versus the current
default control `build/kls_current_default_gap20_t4_r1_ref3_timeout120.jsonl`
at `1.9827861293806655s`, the earlier forced-pipeline prefactor probe
`build/kls_pref_update_gap20_t4_r1_ref3_timeout120.jsonl` at
`1.9475523420313228s`, and the no-force/eager-bitmap probe
`build/kls_pref_update_noforce_gap20_t4_r1_ref3_timeout120.jsonl` at
`1.9157238058535386s`. This confirms that the useful part is the paper-level
skip-unfinished prefactor semantic, not changing the scheduler just to enable
the flag.

A prefix-triggered current-replay queue was also tested and rejected rather than
kept. `build/kls_alg5_queue_pref_gap5_t4_r1_ref3_timeout120.jsonl` measured
`2.9581724708498953s` top-five geomean, and combining it with the prefactor
flag in
`build/kls_alg5_queue_pref_plus_pref_update_gap5_t4_r1_ref3_timeout120.jsonl`
measured `2.9219258140638535s`; `ASIC_100ks` alone rose to about `23s`. That
experiment replayed whole scalar current columns early. It did not implement
the paper's grouped multi-current accumulator, so keeping it would be misleading
and would distract from the still-missing large algorithmic piece.

KLS now allocates retained Algorithm 5 payoff runtime state under
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_WORKSPACE=1`. The allocation is
attached to schedule preparation, not to the numeric hot path, and therefore
serves as the persistent current/target substrate for the missing grouped
multi-current executor rather than as a speed claim. The top-five CKTSO-gap
probe `build/kls_alg5_runtime_workspace_gap5_t4_r1_ref3_timeout120.jsonl`
completed with no failures and `1.7766791317372144s` geomean. On the hard ASIC
rows the runtime counters match the retained descriptors exactly:
`ASIC_320ks` has `3,596` currents, `157,067` workspace rows, and `784,673`
target slots; `ASIC_320k` has `3,755`, `146,027`, and `971,718`; `ASIC_100ks`
has `1,112`, `20,859`, and `188,230`. This keeps BLAS out of the first-order
explanation: the remaining gap is still the paper-level grouped numeric
accumulator that uses this state.

The first numeric consumer of that retained state is now implemented as an
opt-in slot accumulator under
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_SLOT_ACCUM=1`. The flag now
triggers the same payoff-plan preparation as the scalar Algorithm 5 executor,
allocates retained current/target arrays, copies each selected current's prefix
values into its retained workspace slice, accumulates dense-suffix and
L-trailing updates into retained target slots, and publishes each target slot
once. Correctness passed `cmake --build build -j2`,
`ctest --test-dir build --output-on-failure`, `./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_SLOT_ACCUM=1 ./build/kls_smoke`,
and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_SLOT_ACCUM=1 KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_WORKSPACE=1 ./build/kls_smoke`.
The top-five CKTSO-gap probe
`build/kls_alg5_slot_accum_runtime_gap5_t4_r1_ref3_timeout120.jsonl` completed
with no failures and `2.1749853505877623s` geomean. It fired on the hard ASIC
EGraph rows: `ASIC_320ks` recorded `3,596` slot-accumulated runs, `157,067`
rows, `46,839,997` target entries, and `784,673` target slots; `ASIC_320k`
recorded `3,719`, `133,371`, `47,336,990`, and `964,970`; `ASIC_100ks`
recorded `1,112`, `20,859`, `4,668,124`, and `188,230`. A same-binary
workspace-only rerun
`build/kls_alg5_runtime_workspace_rerun_gap5_t4_r1_ref3_timeout120.jsonl`
measured `1.804664999716117s`, so the slot path remains experimental rather
than a default speed path. The direct lesson matches the CKTSO paper gap: a
per-current slot accumulator removes repeated target scatter but still replays
each current column separately. Closing the large slow-case gap requires the
actual grouped multi-current numeric executor that advances several retained
current workspaces through the same producer panel together and then publishes
their target slots, not BLAS threshold tuning or CPU-specific changes.

The Algorithm 5 payoff queue now also has an opt-in prefix-triggered prefetch
probe under
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_QUEUE_PREFETCH=1`. The first
implementation preclaimed current columns before all scalar U predecessors were
complete; the focused CKTSO-gap run
`build/kls_queuepref_prefetch_fixed_gap5_t4_r1_ref3_timeout120.jsonl` timed out
on `ASIC_320ks`, `ASIC_320k`, and `ASIC_100ks` at the 120s per-matrix limit.
That failure is useful: without the paper's separate postfactor task graph,
preclaiming an unfinished scalar column can make the ordinary EGraph pipeline
wait on work that is not actually executable yet.

The committed version therefore uses unclaimed prefetch hints. A prefix trigger
may queue up to `8 * threads` not-yet-ready current columns, but a queue consumer
claims and executes such a hint only if the dependencies are complete when the
hint is popped; otherwise the ordinary pipeline remains responsible for the
column. Correctness passed `cmake --build build -j2`,
`ctest --test-dir build --output-on-failure`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_QUEUE=1 ./build/kls_smoke`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_QUEUE_PREFETCH=1 ./build/kls_smoke`.
The same-build top-five CKTSO-gap control
`build/kls_queuepref_default_samebuild_gap5_t4_r1_ref3_timeout120.jsonl`
measured `1.3987571728147734s` geomean; the existing payoff queue
`build/kls_queuepref_queue_samebuild_gap5_t4_r1_ref3_timeout120.jsonl` measured
`1.7752554941984127s`; and the safe prefetch-hint run
`build/kls_queuepref_prefetch_hint_gap5_t4_r1_ref3_timeout120.jsonl` measured
`1.7573032538167448s`. The ASIC rows claimed only `13`, `13`, and `0` columns
respectively under prefetch hints, while the zero-group `rajat03` row still paid
queue/plan overhead and rose from `0.122342586544s` to `0.37090529484s`.
This rejects queue scheduling as the first-order CKTSO-gap fix. The missing
paper-level part remains grouped multi-current numeric execution, not BLAS
thresholding and not earlier scalar-column queue placement.

The cached EGraph supernode consumer now also implements SubtreeLU's
large-dependent-supernode split rule. When `KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES`
is active and dependency waiting is required, a retained cached-panel run with
at least `2 * threads` producer rows may consume the completed prefix ending
`threads` rows before the producer supernode tail, then leave the tail for the
next dependency-run iteration. `KLS_ENABLE_EGRAPH_SUPERNODE_SPLIT=0` disables
the split for same-binary A/B checks; the existing ragged-L Algorithm 5
prefactor split uses the same helper.

Correctness passed `cmake --build build -j2`,
`ctest --test-dir build --output-on-failure`,
`KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1 ./build/kls_smoke`, and
`KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=1 KLS_ENABLE_EGRAPH_SUPERNODE_SPLIT=0 ./build/kls_smoke`.
The focused top-five CKTSO-gap default control
`build/kls_subtree_split_default_gap5_t4_r1_ref3_timeout120.jsonl` measured
`1.419884637649015s` geomean. The cached-supernode path with the split disabled
`build/kls_subtree_split_supernode_off_gap5_t4_r1_ref3_timeout120.jsonl`
measured `1.5582262938112768s`; with the split enabled,
`build/kls_subtree_split_supernode_on_gap5_t4_r1_ref3_timeout120.jsonl`
measured `1.5670838569234982s`. The split reduced cached rows on
`ASIC_100ks`, `ASIC_320ks`, and `rajat03`, but the full cached path still lost
to the default scalar EGraph path. Forced SCOTCH ordering was also rejected on
this focused set: `build/kls_order_scotch_gap5_t4_r1_ref3_timeout120.jsonl`
measured `2.8342883267318753s`, while auto already selected METIS on the hard
ASIC rows. This keeps the next high-value work on the grouped Algorithm 5
multi-current executor, not SCOTCH promotion, BLAS thresholds, or cached-panel
splitting alone.

The retained Algorithm 5 prefix handoff is now real under
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_PREFIX_PREP=1`. A
producer-prefix trigger can mark a current slot as preparing, apply already
published advance dependencies into private scratch, retain the prefix U
coefficients and target deltas, publish the slot as ready, and let the normal
ragged-L column path consume that retained prefix instead of recomputing the
producer update. Current-slot flags are atomic because producer preparation and
scalar column consumption can occur on different EGraph workers. The blocked
smoke fixture now saves/restores this new environment variable so Algorithm 5
experiments do not leak into the cached-supernode stats test.

Correctness passed `cmake --build build -j2`,
`ctest --test-dir build --output-on-failure`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_PREFIX_PREP=1 ./build/kls_smoke`,
and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_PREFIX_PREP=1 KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_EXEC=1 ./build/kls_smoke`.
The top-five CKTSO-gap refactor probe
`build/kls_prefix_prep_final_gap5_t4_r1_ref3_timeout120.jsonl`
completed with no failures. Its `refactor_seconds_avg` geomean was
`0.14802499172433736s`, versus the same-source default
`build/kls_prefix_prep_default_gap5_t4_r1_ref3_timeout120.jsonl` at
`0.05883573675690034s` and the retained-workspace-only check
`build/kls_prefix_prep_workspace_check_gap5_t4_r1_ref3_timeout120.jsonl` at
`0.11620079150786915s`. The path did consume retained prefixes on the ASIC
rows: `ASIC_320k` consumed `1,974` runs / `107,853` rows / `697,212` target
slots, `ASIC_320ks` consumed `2,029` / `127,978` / `593,344`, and
`ASIC_100ks` consumed `354` / `12,206` / `110,744`. This confirms the
paper-level handoff but rejects the per-current retained-prefix replay as a
default speed fix. The remaining gap is the larger grouped executor: prepare
several current workspaces and advance them through the shared producer panel
together, rather than repeating scatter and advance work once per current.

`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_COMPLETE=1` now adds the
next direct paper-level scheduling slice: a retained final-dependency trigger
map for prepared Algorithm 5 current slots. Prefix triggers still prepare
current-state workspaces in grouped batches; additionally, KLS maps each
single-run current slot to its final U dependency and, when that dependency
publishes, claims and dispatches the prepared current column through the
existing direct-prefix completion path. This fixes the earlier structural issue
where grouped completion at the prefix trigger could only claim `3` to `4`
ASIC currents because most prepared currents still had suffix dependencies.

Correctness passed `cmake --build build -j2`,
`git diff --check`, `ctest --test-dir build --output-on-failure`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_COMPLETE=1 ./build/kls_smoke`.
The first full 41-row CKTSO-gap run before final-trigger activation,
`build/kls_alg5_group_complete_gap5_t4_r1_ref3_timeout120.jsonl`, behaved like
the old direct-complete path and had two 120s timeouts. After adding the final
trigger map, the top-five focused run
`build/kls_alg5_group_complete_final_trigger_gap5_t4_r1_ref3_timeout120.jsonl`
completed but remained slower than default: SPICE-cycle geomean `2.4134s`
versus `1.4629s` in `build/kls_current_default_gap5_t4_r1_ref3_timeout120.jsonl`.
The ASIC rows show that the new trigger really fires:
`ASIC_320ks` claimed `289` prepared currents with `103` claim waits,
`ASIC_320k` claimed `337` with `104` waits, and `ASIC_100ks` claimed `75` with
`27` waits. This closes the prefix-trigger timing gap from the papers, but it
also rejects delayed scalar dispatch as the missing performance mechanism. The
next plausible paper-aligned step is not BLAS thresholding; it is to execute
multiple prepared current workspaces inside one grouped numeric task so claimed
columns do not re-enter the ordinary per-column dispatcher and make other
workers wait on scalar follow-on work.

The follow-up wait-free prepared-completion slice keeps that same guarded
Algorithm 5 surface but removes one scalar scheduler artifact. A final-triggered
prepared current has already passed `kls_egraph_refactor_column_dependencies_done`,
so the claim processor now first calls the BTF direct-prefix completion helper
with dependency waits disabled, after doing the same off-block input copy the
normal BTF dispatcher would have done. If the prepared-state shape does not
match, it falls back to the ordinary claimed-column dispatcher.

Correctness passed `cmake --build build -j2`, `git diff --check`,
`ctest --test-dir build --output-on-failure`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_COMPLETE=1 ./build/kls_smoke`.
The same-build top-five control
`build/kls_default_after_waitfree_gap5_t4_r1_ref3_timeout120.jsonl` measured
`1.4558s` SPICE-cycle geomean. The guarded wait-free completion run
`build/kls_alg5_group_complete_waitfree_final_gap5_t4_r1_ref3_timeout120.jsonl`
completed all five rows but still measured `2.3694s`. The ASIC rows did use the
path: `ASIC_320ks` claimed `921` prepared columns and still processed `1,360`
current-state seed runs / `702,618` skipped dependencies / `1,147,422` state
rows; `ASIC_320k` claimed `1,064` with `1,357` / `742,805` / `1,223,719`;
`ASIC_100ks` claimed `215` with `309` / `179,851` / `282,959`. This confirms
the wait checks were not the large gap.

A direct in-place live-state suffix attempt was also tested and rejected before
commit. Mutating the retained current-state buffer as the post-prefix workspace
matched the desired paper direction superficially, but the top-five run failed
the three ASIC rows with singular-matrix errors. The committed code therefore
keeps the conservative private/delta state completion. The missing piece is
still a designed multi-current numeric owner with correct lifetime and rollback
rules, not ad hoc mutation of the retained seed buffers.

The retained current-state lifetime is now represented explicitly: prepared
current slots move from `READY` to `COMPLETING` before the guarded direct final
state completion reads them, restore to `READY` on fallback/error, and clear to
`EMPTY` only after a successful completion. This closes the state ownership gap
that made the rejected in-place mutation unsafe, but it does not by itself
provide the paper's grouped numeric executor. Correctness passed
`cmake --build build -j2`, `git diff --check`,
`ctest --test-dir build --output-on-failure`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_COMPLETE=1 ./build/kls_smoke`.
The guarded top-five run
`build/kls_alg5_group_complete_owned_state_gap5_t4_r1_ref3_timeout120.jsonl`
completed with no failures but measured `2.4617s` SPICE-cycle geomean, while
the same-build default control
`build/kls_default_after_owned_state_gap5_t4_r1_ref3_timeout120.jsonl`
measured `1.4734s`. Retained-state seed runs dropped on the ASIC rows
(`ASIC_320ks` `1360 -> 887`, `ASIC_320k` `1357 -> 1055`, `ASIC_100ks`
`309 -> 85` versus the previous wait-free guarded artifact), confirming that
owning the slot changes lifetime behavior; the speed loss confirms again that
the main gap is not BLAS or scalar wait checks, but the absence of a grouped
multi-current numeric task that advances several owned current workspaces
together.

The next direct live-state retry is now correctness-safe but still rejected as
the large speed fix. `KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_LIVE_STATE=1`
adds a logged in-place retained-state suffix path under the guarded direct
final-state completion: it requires the owned `COMPLETING` slot state, records
the original value of each retained state entry before mutation, restores the
log on fallback/error, and clears the slot only after a successful completion.
Correctness passed `cmake --build build -j2`, `git diff --check`,
`ctest --test-dir build --output-on-failure`, and
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_COMPLETE=1 KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_LIVE_STATE=1 ./build/kls_smoke`.
The top-five CKTSO-gap run
`build/kls_alg5_group_complete_live_owned_gap5_t4_r1_ref3_timeout120.jsonl`
completed with no failures but measured `2.5316s` geomean, versus `2.4617s`
for `build/kls_alg5_group_complete_owned_state_gap5_t4_r1_ref3_timeout120.jsonl`.
The live path did remove retained-row copy accounting on the ASIC rows
(`seed_rows` went to zero while seed-run counts stayed essentially the same),
but refactor time rose by about `3%`. This confirms that the earlier singular
failures were a lifetime/rollback bug and that copied retained rows are not the
primary CKTSO-gap cause. The missing algorithm remains a grouped multi-current
numeric owner that applies the suffix dependencies to several live current
states together.

The next paper-level row-refactor check moved away from BLAS thresholds and back
to SubtreeLU Algorithm 6 scheduling. The forced row-refactor ASIC runs were using
the retained separator flop queue, but the private side reported only one active
private thread: `ASIC_320ks`, `ASIC_320k`, and `ASIC_100ks` each had
`row_refactor_last_separator_flop_private_threads=1` despite a 4-thread run and
hundreds of retained separator components. KLS now assigns the remaining active
separator components, rather than each whole surviving candidate subtree, to the
least-loaded private thread queue. This keeps separator-crossing groups in the
pipeline closure but prevents a single large surviving private subtree from
serializing all private component work.

Correctness passed `cmake --build build -j2`, `git diff --check`, and
`ctest --test-dir build --output-on-failure`. The forced unchecked row-refactor
top-five final-source run
`build/kls_row_sep_component_balance_final_gap5_t4_r1_ref3_timeout120.jsonl`
measured `4.2548s` geomean, a `1.075x` speedup over the prior forced-row artifact
`build/kls_row_refactor_current_gap5_t4_r1_ref3_timeout120.jsonl` at `4.5736s`.
The ASIC rows now report `row_refactor_last_separator_flop_private_threads=4`;
`ASIC_320ks` improved from `32.8883s` to `29.2094s`, `ASIC_320k` from
`37.7034s` to `35.8048s`, and `ASIC_100ks` from `19.7514s` to `17.7677s`.
This is worth retaining because it closes a clear Algorithm 6 implementation
gap, but it is not a CKTSO-gap closer: the same forced row run is still
`2.888x` slower than the default top-five control
`build/kls_default_after_owned_state_gap5_t4_r1_ref3_timeout120.jsonl`.

Two follow-up probes explain why the selector should remain conservative. With
`KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=auto`, the same forced-row top-five
artifact
`build/kls_row_sep_component_balance_native_auto_env_gap5_t4_r1_ref3_timeout120.jsonl`
improved to `4.1008s` and executed compact dense panels on the ASIC rows, but it
was still `2.783x` slower than default and regressed small rows such as
`gemat12`. A finer private group-level balance attempt was rejected and removed:
`build/kls_row_sep_group_balance_gap5_t4_r1_ref3_timeout120.jsonl` balanced the
private work almost perfectly, but worsened the forced-row geomean to `4.4904s`
and slowed all three ASIC rows. The result is a useful negative signal: the
remaining large gap is not raw private queue balance or small-BLAS dispatch, but
the scalar row/segment numeric executor and the missing coarse grouped numeric
owner described by the CKTSO/SubtreeLU papers.

A current-source scheduler interleaving probe also rejects the remaining
private/pipeline handoff as the first-order gap. SubtreeLU says the
refactorization private-to-pipeline barrier can be removed because the sparsity
pattern is reused. KLS already removes the hard level barrier for the retained
separator FLOP queue, but a trial change let each worker non-blockingly consume
one ready separator chain between private groups. Correctness passed
`cmake --build build -j2`, `git diff --check`,
`./build/kls_smoke`, and `ctest --test-dir build --output-on-failure`, but the
forced row-refactor top-five run
`build/kls_row_sep_interleave_gap5_t4_r1_ref3_timeout120.jsonl` measured
`4.2960s` geomean versus `4.2548s` for the retained component-balance baseline.
The ASIC rows slowed by about `1-2%`, and the local-ready counters dropped
rather than exposing more useful overlap. The probe was removed.

A current-source native EGraph supernode A/B likewise keeps BLAS and existing
coarse-update switches out of the lead-cause slot. With
`KLS_ENABLE_EGRAPH_SUPERNODE_UPDATES=cached`, the top-five default-path focus
run `build/kls_egraph_supernode_cached_gap5_t4_r1_ref3_timeout120.jsonl`
measured `1.4449s` geomean versus `1.4624s` for the same current-source default
subset, but the gain came from `rajat03`; the ASIC rows were neutral or slower
(`ASIC_320ks` `1.060x`, `ASIC_320k` `1.001x`, `ASIC_100ks` `1.085x` relative
to default). Full native updates were worse at `1.6025s` in
`build/kls_egraph_supernode_updates_gap5_t4_r1_ref3_timeout120.jsonl`. These
paths are paper-aligned scalar/native supernode experiments, not BLAS calls, but
they still do not provide the missing CKTSO/SubtreeLU grouped multi-current
numeric owner.

A broader current-source sweep confirms that none of the existing scheduler or
coarse-update switches should be promoted as the next default. Forced cached
EGraph supernode updates over the first 20 CKTSO-gap rows completed cleanly in
`build/kls_egraph_supernode_cached_gap20_t4_r1_ref3_timeout120.jsonl`, but the
geomean regressed to `2.1650s` versus `2.0798s` for the current default. Only
`coupled` and `OPF_3754` improved materially; most EGraph rows slowed even when
the cached update counters fired. Forced `KLS_ENABLE_EGRAPH_READY_QUEUE=1` was
much worse on the same 20-row focus:
`build/kls_egraph_ready_queue_gap20_t4_r1_ref3_timeout120.jsonl` measured
`3.3987s`, with large regressions on nearly every row. A forced SCOTCH-ordering
probe was stopped after the first six rows because setup/runtime had already
exceeded the useful probe budget; the completed rows were slower than current
`auto` ordering (`ASIC_320ks` `16.3754s`, `ASIC_320k` `16.8603s`,
`gemat12` `0.0949s`, `rajat03` `0.2090s`, `ASIC_100ks` `6.0265s`,
`onetone2` `7.5626s`). These results reject broad cached-supernode, ready-queue,
and forced-SCOTCH promotion as general CKTSO-gap fixes.

A guarded grouped-terminal Algorithm 5 implementation was also rejected and
removed. The trial added a publisher that consumed already prepared grouped
prefix items directly when the prefix covered all U dependencies, preserving
the existing sorted pipeline release for completed claims. Correctness passed
normal smoke, group-complete smoke, `ctest --test-dir build
--output-on-failure`, and `git diff --check`, but
`build/kls_alg5_group_complete_group_terminal_gap5_t4_r1_ref3_timeout120.jsonl`
regressed the guarded top-five group-complete geomean to `2.5191s` versus
`2.4617s` for
`build/kls_alg5_group_complete_owned_state_gap5_t4_r1_ref3_timeout120.jsonl`.
The terminal surface was tiny on the hard ASIC rows, so the added claim/scan
work outweighed the saved per-column dispatch. The remaining paper gap is still
not terminal publication; it is a real grouped multi-current suffix executor
that advances several live current workspaces together.

To make that remaining gap measurable, KLS now reports grouped Algorithm 5
suffix counters whenever the grouped payoff descriptor is built:
`refactor_supernode_algorithm5_payoff_group_suffix_deps`,
`refactor_supernode_algorithm5_payoff_group_max_suffix_deps`,
`refactor_supernode_algorithm5_payoff_group_suffix_update_entries`, and
`refactor_supernode_algorithm5_payoff_group_max_run_suffix_update_entries`.
These count the U dependencies and `1 + Llen(dep)` update-entry proxy left
after the current advance and prefix work. The focused top-five grouped-prep
run
`build/kls_alg5_group_suffix_stats_gap5_t4_r1_ref3_timeout120.jsonl`
reported large non-terminal suffix surfaces on the hard ASIC rows:
`ASIC_320ks` had `1,212,799` suffix dependencies and `393,218,349` suffix
update entries, `ASIC_320k` had `1,369,088` and `503,781,567`, and
`ASIC_100ks` had `415,992` and `132,383,507`. Across the five-row focus set,
the grouped descriptor exposed `2,997,879` suffix dependencies and
`1,029,383,423` suffix update entries over `8,463` positioned runs. This keeps
BLAS and terminal publication out of the lead-cause slot: the paper-aligned
work still missing from KLS is a grouped suffix/current-state executor that can
amortize this remaining scalar update surface.

A runtime final-trigger batch probe rules out final-trigger completion as the
place to recover that grouped suffix surface. KLS now reports
`refactor_last_supernode_algorithm5_payoff_final_trigger_batches`,
`refactor_last_supernode_algorithm5_payoff_final_trigger_multi_batches`,
`refactor_last_supernode_algorithm5_payoff_final_trigger_claims`, and
`refactor_last_supernode_algorithm5_payoff_final_trigger_suffix_deps`. The
focused grouped-complete run
`build/kls_alg5_group_complete_final_trigger_surface_gap5_t4_r1_ref3_timeout120.jsonl`
measured `2.5581s` geomean, with `412` final-trigger batches, `0`
multi-claim batches, `412` claims, and only `54,739` suffix dependencies reached
through that path. The same rows still expose `2,997,879` planned grouped
suffix dependencies and `1,029,383,423` suffix update-entry proxy work in the
descriptor. So the CKTSO/SubtreeLU paper gap is not a missed final-trigger
batch: it is earlier, where grouped current-state suffix work should be owned
and advanced before individual final-trigger claims become visible.

The grouped suffix surface now also records whether those suffix dependencies
are unique or repeated inside each Algorithm 5 payoff group:
`refactor_supernode_algorithm5_payoff_group_suffix_unique_deps`,
`refactor_supernode_algorithm5_payoff_group_suffix_duplicate_deps`,
`refactor_supernode_algorithm5_payoff_group_suffix_shared_deps`,
`refactor_supernode_algorithm5_payoff_group_suffix_max_dep_fanout`, and
`refactor_supernode_algorithm5_payoff_group_suffix_duplicate_update_entries`.
The focused grouped-prep run
`build/kls_alg5_group_suffix_sharing_gap5_t4_r1_ref3_timeout120.jsonl`
completed all five rows with clean residuals and `1.7751s` geomean. The hard
ASIC rows show the intended next executor target clearly: `ASIC_320ks` had
`1,212,799` suffix deps but only `212,759` group-local unique deps, with
`1,000,040` duplicate deps, `111,618` shared dependency keys, max fanout `173`,
and `348,084,511` duplicate update-entry proxy work. `ASIC_320k` had
`1,369,088` / `207,442` / `1,161,646` / `112,649` / `223` /
`454,748,014`; `ASIC_100ks` had `415,992` / `139,796` / `276,196` /
`68,503` / `79` / `98,483,225`. Across the five-row focus set, `2,437,882` of
`2,997,879` suffix dependencies were duplicate group-local producer columns and
`901,315,750` of the suffix update-entry proxy work was duplicate-after-first
surface. This makes the next paper-level step more concrete: the grouped
multi-current owner should batch suffix advancement by producer dependency
inside a payoff group, not only retain one current state per scalar consumer.

KLS can now retain that producer-keyed suffix surface in the Algorithm 5 payoff
descriptor with
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_SUFFIX_MAP=1`. The internal
map is `group -> suffix producer column -> (run, U position)` and is built
after the existing grouped current-state and trigger maps are complete,
validated against the measured unique-dependency and suffix dependency totals.
It is intentionally descriptor-only: no numeric path uses it yet, because the
previous unsafe suffix-advance probe showed that mutating live current states
outside a clear owner protocol can introduce residual drift. The map gives the
next executor the paper-aligned shape it needs to advance all occurrences of
one producer inside a payoff group through owned current states. The focused
top-five grouped-prep benchmark with the map enabled in
`build/kls_alg5_group_suffix_map_on_gap5_t4_r1_ref3_timeout120.jsonl`
completed with clean residuals, the same `2,997,879` suffix dependencies and
`2,437,882` duplicate producer occurrences, and `1.8473s` geomean versus
`1.8208s` with the map flag off. The descriptor is kept behind its own flag
because it is still preparation for the grouped suffix executor, not a speedup
by itself.

To keep the next refactor work grounded in the repeated numeric pass rather
than the initial factorization, `kls_bench` now has an optional callgrind hook:
set `KLS_BENCH_CALLGRIND_REFACTOR=1` and run under callgrind with collection
disabled at start. The hook starts instrumentation, zeros stats, dumps stats,
and stops instrumentation around only the repeated `kls_refactor` loop. A
focused `ASIC_100ks` profile using
`build/callgrind_asic100ks_refactor.out.1` collected `2,669,453,961`
instruction references in the repeated refactor region. The profile is
decisive: `kls_egraph_refactor_apply_btf_scalar_dep` accounts for
`2,460,170,381` instructions (`92.16%`), called from
`kls_egraph_refactor_btf_unscaled_column`, and the call count matches
`refactor_dependency_edges=1,556,952`. In the same run all default
supernode/panel/Algorithm 5 payoff execution counters were zero. This rejects
another suffix-only or BLAS-threshold tweak as the next first-order fix for
the slow ASIC refactor cases. The missing paper-level executor is still the
coarse grouped numeric owner that reduces those 1.56M scalar BTF dependency
applications, either by batched prefix advancement or by a multi-current
producer-panel kernel.

A direct-default scalar supernode fallback was tested and rejected before
commit. The prototype let the existing direct supernode dependency helper run
without the cached-panel gate and called it from the plain BTF scalar loop. On
the top-five CKTSO-gap focus,
`build/kls_direct_supernode_plain_gap5_t4_r1_ref3_timeout120.jsonl` regressed to
`1.9612s` geomean versus the same-source default artifact
`build/kls_prepared_suffix_default_gap5_t4_r1_ref3_timeout120.jsonl` at
`1.4287s`. It applied only tiny runs on the hard rows: in the last refactor
`ASIC_320ks` applied `488` supernode rows, `ASIC_320k` applied `1,004`, and
`ASIC_100ks` applied `4,025`, out of `1.56M-1.72M` scalar dependencies. Adding
the same 16-row minimum used by cached EGraph supernodes removed all applied
updates but still measured `1.8815s` in
`build/kls_direct_supernode_min16_gap5_t4_r1_ref3_timeout120.jsonl`. The source
change was reverted. This closes the narrow "just enable the direct supernode
helper in the scalar loop" idea: the default gap is not hidden by an inactive
direct helper; the probe cost is visible before any paper-sized batch is formed.

The new opt-in group-L state-plan probe
`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_GROUP_STATE=1` keeps compact
per-run row-state metadata for grouped consumer-plan runs. It first exposed a
too-narrow payoff filter: on the hard `ASIC_320ks` and `ASIC_320k` rows the
old payoff gate retained zero state rows despite `69,341` and `68,260`
advance-valid candidate runs. KLS now marks advance-valid group-L shapes
separately from payoff-valid execution shapes, so the diagnostic can retain
advance-state metadata without enabling those groups in the numeric executor.
The focused top-five probe
`build/kls_group_l_state_plan_gap5_t4_r1_ref3_timeout120_v2.jsonl` then showed
the full advance surface is large: `ASIC_320ks` retained `69,341` state runs,
`28,823,026` rows, max `2,546` rows, and `231MB`; `ASIC_320k` retained
`68,260` / `29,211,275` / `2,658` / `234MB`; `ASIC_100ks` retained
`152,496` / `72,670,355` / `2,588` / `583MB`. The diagnostic itself is
therefore intentionally not a speed path: the same run regressed to `11.7841s`
geomean because it materializes the full surface before numerics. The default
focused control stayed normal at `1.4497s` in
`build/kls_group_l_state_default_gap5_t4_r1_ref3_timeout120.jsonl`. This
confirms the refactor gap is large, structured, and paper-aligned, but the
next implementation should stream or own current states inside a grouped
executor instead of pre-materializing every advance row.

`KLS_ENABLE_REFACTOR_SUPERNODE_CONSUMER_PLAN_GROUP_STATE_FOCUS=1` now makes
that state diagnostic bounded and self-contained. The focus flag triggers the
retained consumer plan, the group-L cache, and the state scan, but only stores
the largest advance-valid groups that fit the current focus cap
(`64` groups and `64MB` of state rows plus run pointers). It still counts the
full candidate surface, so it is a descriptor for a future grouped executor
rather than a speed path.

On `ASIC_100ks` with 4 threads, one factor, and one refactor, the focus probe
`build/kls_group_l_state_focus_asic100ks_t4_r1_ref1.json` reported a clean
residual (`relative_residual_l2=1.9225e-15`), `1,556,952` dependency edges,
`8,434` candidate groups, `152,496` candidate runs, and `72,670,355`
candidate rows. The bounded state retained `16` groups, `7,375` runs, and
`8,228,308` rows at exactly `67,108,864` bytes, with
`refactor_supernode_consumer_plan_group_l_state_storage_limited=1`. The
same-source no-state control
`build/kls_group_l_state_focus_asic100ks_default_t4_r1_ref1.json` kept group-L
state counters at zero and refactored in `0.0390s`, while the focus scan took
`12.1663s`. That cost is expected because focus mode still walks the whole
advance surface; the important result is that KLS can now inspect the largest
grouped refactor states without allocating the previous hundreds of megabytes
of row-state arrays. The next refactor implementation should consume this
shape as an owned/streaming group executor, not by materializing another
diagnostic table.

`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_ADVANCE_PREP=1` now
adds an opt-in grouped pre-prefix advance executor for the Algorithm 5 payoff
path. It seeds each claimed current into its planned sparse current-state rows
once, advances shared pre-prefix producers across the group in producer order,
then reuses the existing prefix application. Correctness smoke tests pass for
the plain flag, direct-prefix final-state, and group-complete combinations.
The focused top-five CKTSO-gap benchmark was a negative result:
`build/kls_alg5_group_advance_prep_gap5_t4_r1_ref3_timeout120.jsonl` measured
`3.3473s` geomean, versus `2.1922s` for the equivalent scalar direct-current
prep control
`build/kls_scalar_direct_current_advance_map_gap5_t4_r1_ref3_timeout120.jsonl`.
The path did run on the ASIC cases (`ASIC_320ks` seeded `4,011` current states
and `3,565,150` sparse rows), but sparse-state row lookup and per-current row
updates dominated the producer reuse. This rejects a sparse binary-search
grouped pre-prefix replay as the missing large CKTSO-paper mechanism. The
next refactor attempt should avoid per-row sparse lookup in the inner numeric
loop, likely by owning a dense or position-coded multi-current state for a
bounded producer-panel kernel.

`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_ADVANCE_PREP_HASH=1`
tests the next narrower version of that idea: the grouped pre-prefix advance
executor builds a compact per-batch row-to-state-position hash table, so the
inner L-row update no longer binary-searches every sparse current-state row.
This correctly activates the same direct-current-state machinery as
`GROUP_ADVANCE_PREP` and passes smoke/ctest. The fixed top-five CKTSO-gap run
`build/kls_alg5_group_advance_prep_hash_gap5_t4_r1_ref3_timeout120_fixed.jsonl`
measured `2.8283s` geomean. That is much better than the binary-search grouped
replay at `3.4233s`, but still worse than scalar direct-current prep at
`2.3865s` and default at `1.4532s` on the same source. Counters confirm the
hash path ran on the ASIC cases (`ASIC_320ks` seeded `4,010` states and
`3,564,401` sparse rows). This narrows the loss: row-position lookup alone was
not the large missing CKTSO mechanism. The remaining overhead is the doubled
current-state seeding and per-current sparse state update volume; the next
attempt should avoid materializing both scalar and grouped current states, or
own a bounded dense producer-panel state instead of hashing sparse rows per
batch.

The retained Algorithm 5 pre-prefix producer map was tested as the next direct
refactor fix and rejected before commit. The prototype kept the same grouped
sparse current-state semantics, but replaced runtime fanout discovery with the
stored `group_advance_map_*` producer-to-run references. Correctness passed
(`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_ADVANCE_PREP=1
./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_ADVANCE_PREP_HASH=1
./build/kls_smoke`, and `ctest --test-dir build --output-on-failure`), but the
focused same-source benchmark regressed:
`build/kls_alg5_group_advance_map_hash_gap5_t4_r1_ref3_timeout120.jsonl`
measured `3.0131s` SPICE-cycle geomean, versus `2.4192s` for scalar
direct-current prep in
`build/kls_mapchange_scalar_direct_current_gap5_t4_r1_ref3_timeout120.jsonl`
and `1.4674s` for default in
`build/kls_mapchange_default_gap5_t4_r1_ref3_timeout120.jsonl`. On
`ASIC_320ks`, refactor time rose to `0.2877s`, worse than the previous hashed
grouped run at `0.2400s`. This shows the retained full-group map is the wrong
runtime granularity: each prefix-trigger batch claims only a small active subset,
so iterating full symbolic group fanout can add work. The CKTSO/SubtreeLU gap
therefore remains the active-batch grouped numeric executor described by the
papers, not merely a faster lookup for the already retained producer map.

An active-batch producer dependency map was also prototyped and rejected before
commit. The source patch added
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_ADVANCE_ACTIVE_MAP=1`
and kept only each claimed current's live next pre-prefix dependency in a
worker-local bucket map, avoiding the old grouped advance loop's repeated scan
over every claimed current. It built cleanly and passed `./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_ADVANCE_ACTIVE_MAP=1
./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_ADVANCE_ACTIVE_MAP=1
KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_ADVANCE_POS=1
./build/kls_smoke`, and `ctest --test-dir build --output-on-failure`.
The focused same-source top-five CKTSO-gap benchmark was negative:
`build/kls_active_map_default_gap5_t4_r1_ref3_timeout120.jsonl` measured
`1.4820s` default geomean,
`build/kls_group_advance_prep_gap5_t4_r1_ref3_timeout120.jsonl` measured
`3.3484s` for the existing grouped-advance prep,
`build/kls_active_map_gap5_t4_r1_ref3_timeout120.jsonl` measured `4.5098s`,
and
`build/kls_active_map_pos_gap5_t4_r1_ref3_timeout120.jsonl` measured
`2.4771s`. Position-coded row lookup helped but was still far slower than
default, while the plain active map was worse than the old grouped loop. This
rejects active dependency bucketing as the clear missing CKTSO/SubtreeLU paper
piece. The dominant loss is still the surrounding grouped sparse-state executor:
it seeds and updates too much per-current state before the producer-panel work
is coarse enough to amortize that cost.

The scalar BTF EGraph dependency hot path now passes the already-built
32-bit L-row alias table directly into `kls_egraph_refactor_apply_btf_scalar_dep`
and into the inline non-wait BTF scatter. This removes a per-dependency lookup
through `solver->refactor_l_indices32` and a helper branch from the path that
callgrind identified as the hard ASIC refactor bottleneck. It is not the missing
CKTSO/SubtreeLU grouped executor, but it is a default-path cleanup on the same
measured scalar dependency stream. Validation passed `cmake --build build -j2`,
`./build/kls_smoke`, and `ctest --test-dir build --output-on-failure`.
The same-source top-five CKTSO-gap baseline before the change measured
`1.4820s` in
`build/kls_active_map_default_gap5_t4_r1_ref3_timeout120.jsonl`; two candidate
runs measured `1.4184s` and `1.4151s` in
`build/kls_btf_i32_direct_gap5_t4_r1_ref3_timeout120.jsonl` and
`build/kls_btf_i32_direct_gap5_t4_r1_ref3_timeout120_r2.jsonl`. A stashed
same-source top-ten baseline measured `2.3084s` in
`build/kls_pre_i32_direct_gap10_t4_r1_ref3_timeout120.jsonl`; the candidate
measured `2.2805s` in
`build/kls_btf_i32_direct_gap10_t4_r1_ref3_timeout120.jsonl`. This is a small
accepted scalar cleanup, not evidence that scalar tuning replaces the still
missing grouped producer/current-state executor.

A follow-up generalization of that cleanup was tested and rejected before
commit. The prototype added a reusable alias-aware scatter helper, marked the
BTF scalar dependency helper always-inline, and threaded local `refactor_l_indices32`
aliases through the single-block, generic, and Algorithm 5 prefix-advance
EGraph scalar loops. It built cleanly and passed `./build/kls_smoke`,
`ctest --test-dir build --output-on-failure`, and `git diff --check`, but the
timing did not justify the added churn. The focused top-five CKTSO-gap run
`build/kls_egraph_alias_direct_gap5_t4_r1_ref3_timeout120.jsonl` measured
`1.4184s`, essentially the same as the retained BTF-only cleanup. On top ten,
one candidate run improved a same-machine stashed baseline from `2.3467s` to
`2.2915s`, but a second candidate pass regressed to `2.3734s` in
`build/kls_egraph_alias_direct_gap10_t4_r1_ref3_timeout120_r2.jsonl`, and the
hard ASIC rows were neutral or slower. The source patch was reverted. This
keeps the accepted direct-i32 change limited to the measured BTF scalar stream
and reinforces that the next useful refactor work is still the grouped
multi-current producer/current-state executor, not broader scalar scatter
plumbing.

The refactor focus was re-checked against the CKTSO and SubtreeLU paper
algorithms after rejecting BLAS-threshold tuning as the primary explanation.
The papers point at row-major up-looking refactorization, Algorithm 5
prefactor/postfactorization, supernode TRSV/GEMV updates, and Algorithm 6
private/pipeline separator scheduling. KLS already has probes or retained
metadata for each of those names, but the forced-row CKTSO-gap evidence shows
they are still wrapped around the old scalar/current-state executor rather than
replacing it. The top-five forced row-refactor run
`build/kls_forced_row_refactor_current_gap5_t4_r1_ref3_timeout120.jsonl`
measured `4.22647s` geomean and put `ASIC_320ks` at `28.977s`
(`refactor_seconds_avg=0.258348`) versus the retained EGraph row in
`build/kls_btf_i32_direct_gap10_t4_r1_ref3_timeout120.jsonl` at `10.862s`
(`refactor_seconds_avg=0.0828143`). The scheduler did use the paper-shaped
separator queue on the ASIC rows, but almost all work became tiny private
groups: `ASIC_320ks` reported `218928` separator-private groups, `125`
pipeline groups, `95` closure promotions, and zero compact-supernode update
rows.

Forcing the more paper-shaped row-panel/supernode knobs on the same worst case
was also negative:
`build/kls_forced_row_refactor_full_panel_asic320ks_t4_r1_ref3_timeout120.jsonl`
measured `37.144s` and `refactor_seconds_avg=0.341396`. It executed
`1,468,143` compact-supernode update rows and `608,619,828` compact-supernode
update entries, including `462,827,352` GEMV entries, `82,867,395` TRSV
entries, `732` compact-supernode batches, and `82,314,839` native-panel blocked
entries. This is more than twice the row-group work estimate
(`285,164,776`) and worse than the scalar forced row-refactor run. The
accepted conclusion is that simply enabling existing compact-panel, native
row-panel, multi-producer, or BLAS-shaped probes will not close the CKTSO gap.
The missing large piece remains a production row-major supernode/current-state
numeric object that avoids duplicated sparse current states and private-group
bookkeeping instead of adding panel work on top of them.

`scripts/decompose_solver_gap.py` now records this distinction directly. It
uses `last_refactor_path` in addition to the first-factor path and reports
`candidate_last_refactor_path`, `row_refactor_panel_overstage_ratio`, and
`row_refactor_private_group_share`. Its `paper_gap_signal` can now classify
forced row probes as `row_refactor_private_scalar_scaffold` or
`row_refactor_panel_overstaged`, while default EGraph rows whose row-refactor
lower bound is already worse are marked `row_refactor_lower_bound_rejected`.
On the forced full-panel `ASIC_320ks` artifact, the signal is
`row_refactor_panel_overstaged` with a panel-overstage ratio of `2.134`.

The refactor focus was rerun on the current tree after the BLAS-threshold
discussion. A same-source top-five control measured `1.46637s` geomean in
`build/kls_current_default_gap5_t4_r1_ref3_timeout120.jsonl`. Forcing the
existing EGraph ready queue with `KLS_ENABLE_EGRAPH_READY_QUEUE=1` measured
`2.49138s` in
`build/kls_ready_queue_probe_gap5_t4_r1_ref3_timeout120.jsonl`, with large
regressions on `ASIC_320ks`, `ASIC_320k`, and `ASIC_100ks`. Explicit
Algorithm 5 prefactor on/off checks were also neutral-to-negative:
`KLS_ENABLE_EGRAPH_ALGORITHM5_PREF_UPDATE=0` measured `1.47113s`
(`build/kls_prefactor_off_probe_gap5_t4_r1_ref3_timeout120.jsonl`) and
`=1` measured `1.48939s`
(`build/kls_prefactor_on_probe_gap5_t4_r1_ref3_timeout120.jsonl`). In the
default hard rows, nearly all modeled refactor work is already in the pipeline
tail (`ASIC_320ks` `273,590,879 / 282,035,457`, `ASIC_320k`
`335,407,299 / 343,575,699`, `ASIC_100ks`
`306,394,308 / 316,640,487`), but all supernode, ragged-L, and Algorithm 5
payoff execution counters remain zero. This rejects ready-queue scheduling and
the scalar prefactor guard as the missing CKTSO/SubtreeLU mechanism on the
current slow rows.

`scripts/decompose_solver_gap.py` now marks that sharper case as
`egraph_scalar_tail_numeric_owner_missing` when refactor work is dominated by
the EGraph pipeline tail, row-refactor auto is lower-bound rejected, and no
grouped supernode/current-state executor counter fires. That label is meant to
keep the next implementation directed at the paper's row-major supernode or
multi-current numeric owner, not at BLAS thresholds, ready queues, or another
thin wrapper around the scalar dependency stream.

The row-refactor path now takes one narrow paper-aligned step toward that
owner: unset `KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR` behaves like `auto` only
when the retained row-refactor pattern already has persistent compact dense
panels and at least `1e7` modeled compact-panel update work. This keeps explicit
`0`, `auto`, and `1` semantics unchanged, skips small pure-panel cases such as
`rajat03`, and lets large ASIC-style row refactors use the retained
row-panel/supernode executor without a manual env override. Same-source
forced row-refactor top-five checks with refactor-repeat 10 measured
`3.20631s` geomean for the structural auto gate versus `3.42501s` with
`KLS_ENABLE_NATIVE_ROW_PANEL_REFACTOR=0`
(`build/kls_structural_native_auto_gap5_t4_r1_ref10_timeout120.jsonl`,
`build/kls_structural_native_auto_off_gap5_t4_r1_ref10_timeout120.jsonl`).
The default selector remained on EGraph/mapped refactor paths at `1.46054s` in
`build/kls_structural_native_auto_default_gap5_t4_r1_ref3_timeout120.jsonl`.
This is worth retaining as a scoped row-panel dispatch fix, but the large CKTSO
gap remains the missing production row-major grouped numeric owner.

The scalar row-refactor row kernel now also uses the retained worker byte
workspace for its per-row dependency-applied bitmap instead of allocating and
freeing that bitmap for every row. This matches the dense/native row-refactor
paths and removes one avoidable source of private-group scaffold overhead.
The focused forced row-refactor top-five run
`build/kls_row_bitmap_workspace_gap5_t4_r1_ref3_timeout120.jsonl` measured
`4.14697s` geomean, compared with the recent forced-row artifacts around
`4.22647s` to `4.5736s`; the same-source default EGraph selector remained
clean at `1.41046s` in
`build/kls_row_bitmap_workspace_default_gap5_t4_r1_ref3_timeout120.jsonl`.
This is a useful refactor-path cleanup, but it is not the paper-level gap
closer: the slow cases still need the production row-major grouped numeric
owner rather than per-row scalar private-group bookkeeping.

The default EGraph refactor loop now also skips the dormant Algorithm 5
direct-prefix seed probe before calling into its helper. The helper already
returned immediately when the opt-in current-state/advance-seed executors were
inactive, but the hard ASIC profile still charged the inactive call chain inside
the repeated scalar dependency stream. The change gates the helper at the four
EGraph refactor loops and preserves the explicit
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_ADVANCE_SEED`
path. Verification passed `./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_ADVANCE_SEED=1 ./build/kls_smoke`,
`ctest --test-dir build --output-on-failure`, and `git diff --check`. Focused
top-five CKTSO-gap timing stayed neutral at the cycle level:
`build/kls_seed_guard_gap5_t4_r1_ref3_timeout120.jsonl` measured `1.41525s`
geomean and repeat
`build/kls_seed_guard_gap5_t4_r1_ref3_timeout120_r2.jsonl` measured
`1.42174s`, versus the retained default artifact at `1.41046s`; refactor-only
geomean was essentially unchanged (`0.00905s` and `0.00913s` versus
`0.00913s`). This is worth keeping as a default hot-loop cleanup, but it does
not change the paper-level conclusion: inactive probes and BLAS thresholds are
not the first-order CKTSO gap; the missing mechanism is still the grouped
row-major/supernode current-state numeric owner.

A bounded dense-current panel inside the Algorithm 5 suffix-group executor was
prototyped and rejected before commit. The experiment added an opt-in
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_SUFFIX_GROUP_DENSE_PANEL`
path that gathered each producer batch's retained sparse current-state values
into a worker-local dense `lcol_len x item_count` panel, applied the shared
producer update there, and scattered the panel back. It passed `./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_SUFFIX_GROUP_ADVANCE=1 ./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_SUFFIX_GROUP_DENSE_PANEL=1 ./build/kls_smoke`,
`ctest --test-dir build --output-on-failure`, and `git diff --check`, but the
focused top-five CKTSO-gap benchmark did not move. Same-source default measured
`7.8135s` SPICE-cycle geomean and `0.06097s` refactor-only geomean in
`build/kls_dense_panel_default_gap5_t4_r1_ref3_timeout120.jsonl`. Existing
suffix-group advance measured `27.0181s` / `0.24840s` in
`build/kls_suffix_group_control_gap5_t4_r1_ref3_timeout120.jsonl`; the dense
panel measured `26.9663s` / `0.24787s` in
`build/kls_suffix_group_dense_panel_gap5_t4_r1_ref3_timeout120.jsonl`. On the
hard ASIC rows, the panel still seeded the same sparse current states
(`ASIC_320k` `1,023,899` rows, `ASIC_320ks` `789,478`, `ASIC_100ks` `75,621`)
and only changed the single-producer update storage. This rejects a bounded
gather/update/scatter panel as the large paper gap. The missing piece has to be
earlier and coarser: avoid materializing the per-current sparse state in this
form, or make the row-major/supernode owner hold the live state across producer
windows instead of copying it around one producer at a time.

A direct retained-current-state seed prototype was also rejected before commit.
The patch changed
`kls_egraph_refactor_seed_algorithm5_payoff_group_item_state` to initialize the
Algorithm 5 retained sparse state directly from the current column's mapped
input entries instead of staging through the worker `x` vector and copying the
retained rows back. Correctness passed `./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_SUFFIX_GROUP_ADVANCE=1 ./build/kls_smoke`,
and `ctest --test-dir build --output-on-failure`, but the hard seeded ASIC
rows did not improve. The binary-search version measured `26.7430s` SPICE-cycle
geomean and `0.24605s` refactor-only geomean in
`build/kls_direct_state_seed_suffix_group_gap5_t4_r1_ref3_timeout120.jsonl`,
versus the previous suffix-group artifact at `27.0181s` / `0.24840s`; however
the seeded rows themselves were neutral or slightly slower
(`ASIC_320k` `0.23594s` vs `0.23583s`, `ASIC_320ks` `0.20829s` vs
`0.20814s`, `ASIC_100ks` `0.09873s` vs `0.09838s`). Replacing the per-entry
binary search with a sorted merge cursor regressed the same top-five suffix
group run to `27.0366s` in
`build/kls_direct_state_seed_merge_suffix_group_gap5_t4_r1_ref3_timeout120.jsonl`.
This rejects the seed-copy staging as the large missing mechanism. The
Algorithm 5 gap is still downstream of, or broader than, the initial seed:
the retained sparse current states and grouped executor ownership model remain
the target.

A branch-light scalar scatter cleanup in the default BTF EGraph refactor was
also tested and rejected before commit. The prototype split the generic
`kls_scatter_subtract*` helpers into checked and already-nonzero variants,
marked `kls_egraph_refactor_apply_btf_scalar_dep` always-inline, and skipped
L-column pointer loads for empty L columns. It passed `cmake --build build -j2`,
`./build/kls_smoke`, `ctest --test-dir build --output-on-failure`, and
`git diff --check`, but a clean same-source A/B against a detached HEAD
baseline showed no robust improvement. On the current top-five CKTSO-gap
manifest, patched runs measured `7.7090s` and `7.7365s` SPICE-cycle geomean in
`build/kls_scalar_scatter_nonzero_gap5_t4_r1_ref3_timeout120.jsonl` and
`build/kls_scalar_scatter_nonzero_gap5_t4_r1_ref3_timeout120_r2.jsonl`, while
the clean baseline worktree measured `7.7531s` and `7.6663s` in
`build/kls_baseline_head_gap5_t4_r1_ref3_timeout120*.jsonl`.
Median refactor-only geomean was effectively identical: `0.06003s` patched
versus `0.06002s` baseline. The patch helped the hard ASIC_320k/ASIC_320ks
median refactor rows by about `1.8%` and `2.2%`, but regressed ASIC_100ks,
G2_circuit, and onetone2 enough to erase the geomean. This keeps the accepted
diagnosis unchanged: the refactor loss is not a missed scalar branch in
`kls_egraph_refactor_apply_btf_scalar_dep`; the paper-sized gap is still the
row-major/supernode current-state numeric owner that reduces the scalar
dependency stream instead of polishing it.

A scheduler-barrier shortcut for the EGraph refactor was tested and rejected
before commit. The paper-level hypothesis was that, after the SubtreeLU-style
cluster-level barriers, a pipeline column should not need to poll
`pipeline_done` for predecessors scheduled before
`refactor_level_ptr[cluster_level_count]`; those predecessors have already
passed the cluster barriers. The prototype routed both
`kls_egraph_refactor_wait_done` and
`kls_egraph_refactor_dependency_done_now` through that schedule-position check.
Correctness passed `cmake --build build -j2`, `./build/kls_smoke`,
`ctest --test-dir build --output-on-failure`, and `git diff --check`, but the
timing was not robust enough to keep. The initial uncached helper measured
`7.5566s`, `7.9116s`, and `7.3799s` SPICE-cycle geomean on the top-five
CKTSO-gap manifest in
`build/kls_cluster_done_skip_gap5_t4_r1_ref3_timeout120*.jsonl`, versus clean
baseline runs at `7.7531s`, `7.6663s`, and `7.6065s`
(`build/kls_baseline_head_gap5_t4_r1_ref3_timeout120*.jsonl` and
`build/kls_baseline_c565_gap5_t4_r1_ref3_timeout120.jsonl`). An isolated
higher-repeat ASIC_100ks check was mildly positive (`5.2809s` versus
`5.3952s` in
`build/kls_cluster_done_skip_asic100ks_t4_r2_ref8_timeout120.jsonl` and
`build/kls_baseline_c565_asic100ks_t4_r2_ref8_timeout120.jsonl`), but the full
run still had unstable ASIC_100ks losses. A tightened version that cached the
pipeline-begin schedule position regressed two consecutive full passes to
`7.8890s` and `8.0413s` in
`build/kls_cluster_done_skip_cached_gap5_t4_r1_ref3_timeout120*.jsonl`.
This rejects redundant cluster-predecessor polling as the clear missing
paper mechanism. The remaining gap is not the cluster/pipeline barrier
bookkeeping; it is still the absent coarse numeric owner that prevents the
pipeline tail from being expressed as hundreds of millions of scalar
dependency updates.

A June 30, 2026 cost-gated grouped pre-prefix advance patch was tested and
rejected before commit. The patch tried to make
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_ADVANCE_PREP=1`
closer to the paper's low-advance Algorithm 5 subset: before claiming a
current for the grouped sparse-state batch, it required the retained producer
prefix work (`target_entries + prefix triangular work`) to exceed the modeled
cost of sparse current-state seeding plus pre-prefix advance dependencies.
It also refused to form a grouped sparse-state batch with only one surviving
current. Correctness passed `cmake --build build -j2`, `git diff --check`,
`./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_ADVANCE_PREP=1 ./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_GROUP_ADVANCE_POS=1 ./build/kls_smoke`,
and `ctest --test-dir build --output-on-failure`.
The performance result was decisively negative. The same-source default
top-five CKTSO-gap run
`build/kls_group_adv_paygate_default_gap5_t4_r1_ref3_p3_timeout120.jsonl`
measured `1.3865s` SPICE-cycle geomean, but the gated grouped-advance run
did not produce the first matrix row after several minutes and was interrupted
with an empty
`build/kls_group_adv_paygate_prep_gap5_t4_r1_ref3_p3_timeout120.jsonl`.
A direct single-pass probe on `ASIC_320ks`
(`repeat=1`, `refactor-repeat=1`) hit a 90 second external timeout with the
same gated grouped-advance flag. This rejects a scalar cost gate around the
existing sparse current-state batch as the missing CKTSO/SubtreeLU mechanism.
The problem is not simply selecting fewer high-advance currents; the retained
executor still materializes and mutates sparse per-current state in a form that
is too expensive before it reaches the coarse producer-panel work.

To test whether a dense current-state window could be that missing owner, KLS
now reports the retained Algorithm 5 current-state row span in addition to the
actual retained row count. The diagnostics are exposed as
`refactor_supernode_algorithm5_payoff_current_state_span_rows` and
`refactor_supernode_algorithm5_payoff_current_state_max_span_rows`.
A one-pass top-three CKTSO-gap plan-only probe
(`build/kls_alg5_state_span_plan_gap3_t4_r1_ref1.jsonl`) measured
`10.5621s`, `12.7739s`, and `0.061058s` on `ASIC_320ks`, `ASIC_320k`,
and `gemat12`, for a `2.0196s` geomean. The two ASIC rows retained
`2,568,461` and `2,925,515` sparse current-state rows, but those rows spanned
`703,336,131` and `781,830,290` dense row slots, with per-current maxima of
about the full block (`320,926`). The compact-state variant
(`build/kls_alg5_state_span_compact_gap3_t4_r1_ref1.jsonl`) reduced actual
state rows to `910,070` and `1,067,018`, but still spanned `572,250,409` and
`698,030,903` row slots while regressing the top-three geomean to `4.6545s`.
This rejects a naive dense-window retained-state executor for the hard ASIC
cases: the row ranges are hundreds of times larger than the useful sparse
state. The paper-aligned target remains a grouped sparse current-state owner
that keeps live producer/current batches together without copying or scanning
each current independently.

The default BTF EGraph scalar tail now has an opt-in producer-run surface
diagnostic behind `KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_STATS=1`. The tracker is
on-the-fly and default-off: it only records consecutive scalar dependencies
whose U positions and block-local producer rows advance together, then reports
last-pass and cumulative run count, row count, L-entry surface, and max run
length through `refactor_*_btf_scalar_run_*` stats and benchmark JSON. It does
not change numeric execution. Validation passed `cmake --build build -j2`,
`./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_STATS=1 ./build/kls_smoke`,
`ctest --test-dir build --output-on-failure`, and `git diff --check`.

The first focused probe on `ASIC_100ks` with 4 threads, one factor, and one
refactor (`build/kls_btf_scalar_run_stats_asic100ks_t4_r1_ref1.json`) stayed
residual-clean (`relative_residual_l2=1.9225e-15`) and reported
`107,917` contiguous scalar producer runs covering `1,245,501` of the
`1,556,952` scalar dependency edges, with `312,716,559` summed L entries and max
run length `462`. The same no-env run
(`build/kls_btf_scalar_run_stats_off_asic100ks_t4_r1_ref1.json`) kept all new
counters at zero. `scripts/decompose_solver_gap.py` now labels this measured
case `egraph_scalar_tail_producer_runs_unowned` instead of the broader
`egraph_scalar_tail_numeric_owner_missing`. This strengthens the refactor
direction: the slow tail has enough contiguous producer-run surface for a
coarse row-major/current-state owner, and the missing piece is not another BLAS
threshold, ready queue, or scalar scatter cleanup.

The first direct BTF scalar producer-run executor was implemented behind
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_EXEC=1` and kept default execution unchanged.
The unsafe wait-mode mutation hazard was fixed before validation: the executor
now waits for every dependency in the run before clearing any live workspace
entry. It also records last-pass and cumulative applied-run counters through
`refactor_*_btf_scalar_run_exec_*`, and the gap script treats nonzero applied
runs as `egraph_btf_scalar_run_executor_active`.

The unguarded direct form was residual-clean but too slow on `ASIC_100ks`: it
covered `46,540` producer runs, `1,100,807` rows, and `310,512,292` L entries,
but raised one-refactor `refactor_seconds_avg` from about `0.0425s` to
`0.1145s`. The reason matches the CKTSO/SubtreeLU paper gap: wrapping KLU's
current scalar dependency stream added an in-run/trailing row-range branch to
every L entry instead of preserving the branch-light trailing scatter.

The committed opt-in executor therefore uses a stricter split: while building
the retained LU pointer cache under the exec flag, KLS records whether all L row
lists are strictly ascending. Only in that case can the executor apply the
in-run triangular prefix locally and send the suffix through the existing
scatter primitive. On the current `ASIC_100ks` focus case that guard correctly
prevents execution: the same-source split benchmark
`build/kls_btf_scalar_run_exec_split_off_asic100ks_t4_r1_ref3.json` measured
`refactor_seconds_avg=0.0394436027`, and
`build/kls_btf_scalar_run_exec_split_on_asic100ks_t4_r1_ref3.json` measured
`0.0393036917` with all `refactor_last_btf_scalar_run_exec_*` counters at zero
and `relative_residual_l2=1.92251861e-15`.

Validation for this chunk passed `git diff --check`, `cmake --build build -j2`,
`./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_EXEC=1 ./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_EXEC=1 KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_STATS=1 ./build/kls_smoke`,
and `ctest --test-dir build --output-on-failure`. The result is not a default
speed path; it narrows the paper-aligned next step to an owned row-ordered or
grouped current-state numeric representation instead of more scalar-loop
tuning.

The next refactor slice added that owned row-ordered representation without
changing the KLU-compatible numeric factor: under
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_EXEC=1`, KLS now builds a sorted L-row mirror
with source positions, refreshes each sorted value vector when the EGraph stores
the corresponding L column, and routes the BTF scalar producer-run executor
through that mirror. This turns the previous no-op split into an actual
triangular-prefix/trailing-scatter executor while keeping the path default-off.
The row-lower invariant is checked when the mirror is built, so the hot executor
does not rescan each in-run prefix just for validation.

Focused measurements show that the representation is correct and much closer to
the paper shape, but still not a default speed path. On `ASIC_100ks`, the
default control
`build/kls_btf_scalar_run_sorted_l_off_asic100ks_t4_r1_ref3.json` measured
`refactor_seconds_avg=0.0444740077`. The sorted-mirror executor with a 4-row
floor covered `46,540` runs, `1,100,807` rows, and `310,512,292` L entries, but
measured `0.050005735`; a 16-row floor covered `12,548` runs, `879,371` rows,
and `297,303,361` entries, but still measured `0.049890663`. On `ASIC_320ks`,
the default control
`build/kls_btf_scalar_run_sorted_l_off_asic320ks_t4_r1_ref1.json` measured
`0.085005967`; the sorted-mirror executor measured `0.10847028` with a 4-row
floor and `0.106261211` with a 16-row floor. All focused sorted-mirror runs were
residual-clean (`relative_residual_l2` about `1.9e-15` to `2.1e-15`).

This rejects sorted row order alone as the clear missing CKTSO/SubtreeLU
mechanism. The next paper-aligned executor has to batch multiple current
workspaces through the retained producer run, so the producer's L suffix is
streamed once across a grouped current-state owner rather than replaying a
separate scalar scatter for each current column.

The follow-on slice added the structural grouped-current descriptor for that BTF
producer-run target. `KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUPS=1` keeps numeric
execution unchanged, scans retained BTF U patterns for contiguous scalar
producer runs, groups those runs by producer start, and retains the multi-current
groups as `(producer start, max rows, current, U position, current rows)` arrays.
The benchmark fields `refactor_btf_scalar_run_group_*` report the total run
surface, the subset with multiple current columns, and the row/L-entry surface
that could be reused if a future live current-state executor streams the producer
L suffix once across the group. The gap decomposition script now labels rows
with nonzero `refactor_btf_scalar_run_group_reused_entries` as
`egraph_scalar_tail_grouped_producer_runs_unowned`, separating this paper-level
missing owner from the earlier single-current producer-run executor.

Focused one-refactor probes confirm that this is the right next executor target,
not just a bookkeeping variant. On `ASIC_100ks`, the cached planner artifact
`build/kls_btf_group_cached_on_asic100ks_t4_r1_ref1.json` stayed
residual-clean (`relative_residual_l2=1.9225e-15`), built `13,116`
multi-current producer groups covering `104,280` current-run memberships, and
reported `299,572,567` reused L-entry reads, with up to `362` currents sharing a
producer start and max run length `462`. The cached planner measured
`refactor_seconds_avg=0.043987771`, in the same noise band as the default
one-refactor control. On `ASIC_320ks`,
`build/kls_btf_group_on_asic320ks_t4_r1_ref1.json` built `33,856`
multi-current groups over `121,701` current-run memberships and reported
`270,529,235` reused L-entry reads, with max current fanout `444` and max run
length `586`. This gives the next implementation a concrete target: a live
multi-current BTF executor should consume these retained groups and update
several current workspaces per streamed producer L suffix.

The next diagnostic slice checked whether the existing EGraph/BTF schedule
already creates a natural rendezvous point for that grouped executor. Under
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_WAIT_STATS=1`, KLS now reuses the
retained grouped producer-run descriptor and records live wait overlap when a
current column reaches a grouped scalar producer run whose producer dependency is
not yet published. This is still diagnostic-only; it does not change numeric
execution.

The result is clear on the current hard SPICE cases. On `ASIC_100ks`,
`build/kls_btf_group_wait_on_asic100ks_t4_r1_ref1.json` built the same `13,116`
multi-current groups and `299,572,567` reusable L-entry reads, but reported
zero grouped waits, zero grouped overlap waits, and max live waiters zero
(`relative_residual_l2=1.92251861e-15`). On `ASIC_320ks`,
`build/kls_btf_group_wait_on_asic320ks_t4_r1_ref1.json` likewise built `33,856`
multi-current groups and `270,529,235` reusable L-entry reads, but also reported
zero grouped waits and max live waiters zero
(`relative_residual_l2=2.08436857e-15`). Enabling scalar producer-run stats at
the same time confirms the scalar-run surface is present: `ASIC_100ks` reported
`107,920` producer runs and `312,769,800` touched L entries, while `ASIC_320ks`
reported `230,896` producer runs and `278,005,604` touched L entries. Both
still had zero grouped live waits.

This rejects a simple wait-rendezvous implementation. The paper-level grouped
current-state owner has to be introduced explicitly in the BTF schedule; it
cannot be obtained by attaching work to the current dependency wait point. The
next refactor implementation should make grouped producer runs schedule-visible
and route multiple current workspaces through a producer-owned update, rather
than relying on incidental concurrent waits.

The follow-on scheduler slice made the grouped producer runs schedule-visible
without changing the numeric kernel. Under
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_CLAIMS=1`, KLS now triggers from a
completed producer column, finds retained BTF scalar-run groups for that
producer, and attempts to claim grouped current columns whose full dependency
lists are already published. The hook also fires from barriered cluster columns,
so it sees producer starts before the pipeline tail begins. Benchmark JSON now
separates producer-trigger surface counters from the subset of current columns
actually claimed and dispatched.

Focused probes show that this simple producer-side claim scheduler is still too
late in the dependency stream. On `ASIC_100ks`,
`build/kls_btf_group_claims2_on_asic100ks_t4_r1_ref1.json` saw the full
schedule-visible grouped surface: `13,116` producer triggers/groups and
`104,280` current memberships, with `299,572,567` reusable L-entry reads. It
claimed zero currents. On `ASIC_320ks`,
`build/kls_btf_group_claims2_on_asic320ks_t4_r1_ref1.json` saw `33,856`
producer triggers/groups and `121,701` current memberships, with `270,529,235`
reusable L-entry reads. It also claimed zero currents. Both runs were
residual-clean (`relative_residual_l2=1.92251861e-15` and
`2.08436857e-15`). Same-binary controls without the flag measured
`0.045162553s` on `ASIC_100ks` and `0.084273274s` on `ASIC_320ks`; the claim
probe measured `0.041062299s` and `0.087357188s`, but with zero actual claimed
currents those timing deltas are noise rather than a numeric improvement.

This narrows the next implementation further. Producer completion is now proven
to expose the grouped surface, but waiting for full current-column readiness
still leaves no work to claim. The missing paper-level executor must own partial
current state before all dependencies are complete, advance those live
workspaces through grouped producer runs, and then hand the completed columns
back to the normal publish path. A scheduler-only current claim does not fill
that gap.

The prefix-readiness probe fills in the next dependency detail. Under
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_PREFIX_STATS=1`, KLS now checks each
producer-triggered grouped current membership for two states: all dependencies
before the contiguous producer run are already published, and the full producer
run is already published. With both prefix stats and the old claim probe enabled,
`ASIC_100ks` saw `104,280` schedule-visible memberships, `44,454` prefix-ready
memberships covering `150,680,054` L entries, but only `246` full-run-ready
memberships covering `2,937` L entries; full-column claims remained zero.
`ASIC_320ks` saw `121,701` memberships, `49,742` prefix-ready memberships
covering `149,632,708` L entries, and only `18` full-run-ready memberships
covering `1,152` L entries; full-column claims again remained zero. Both probes
were residual-clean (`1.92251861e-15` and `2.08436857e-15`).

This rejects producer-start dispatch as the whole refactor answer. The useful
surface is the prefix-ready current state, not an immediately executable full
producer run. The next implementation should create a retained current-state
owner for these prefix-ready memberships and register wakeups on later producer
columns in the contiguous run, so grouped L-stream reuse is applied when the run
actually becomes ready instead of waiting until the whole current column is
done.

The run-end wake probe adds that missing second trigger. Under
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_WAKE_STATS=1`, KLS now builds a
reverse descriptor from the last producer in each grouped contiguous run to the
group memberships that should wake there. At producer-start completion it arms
memberships whose prefix dependencies are already published; at producer-run-end
completion it records the armed memberships whose full run is now ready.
`ASIC_100ks` armed `44,454` memberships covering `150,680,054` L entries and
woke `44,072` memberships covering `150,430,068` L entries. `ASIC_320ks` armed
`49,740` memberships covering `149,632,370` L entries and woke `49,643`
memberships covering `149,422,154` L entries. Both probes were residual-clean
(`1.92251861e-15` and `2.08436857e-15`).

This confirms that the scheduler-side shape needed by the paper-level retained
current-state owner is present: most prefix-ready memberships have a matching
run-end wake. The next numeric implementation should attach a compact retained
workspace to the armed membership and consume those wake-ready batches with one
producer-owned L-stream pass across multiple current states.

The retained-state wake probe now attaches that compact workspace under
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_STATS=1`. The flag implies the
run-end wake probe, plans compact per-membership row maps, materializes the
prefix-ready state at arm time, and advances that state when the contiguous
producer run becomes ready. `ASIC_100ks` planned `58,663,555` compact rows,
materialized `44,334` currents / `13,633,139` rows / `1,409,330` prefix
dependencies, advanced `43,776` currents through `485,305` state rows and
`127,452,176` L entries, and recorded zero state rejects. `ASIC_320ks` planned
`50,818,849` compact rows, materialized `49,765` currents / `15,442,593` rows /
`1,941,488` prefix dependencies, advanced `47,515` currents through `493,567`
state rows and `125,352,798` L entries, and also recorded zero state rejects.
Both runs remained residual-clean (`1.92251861e-15` and `2.08436857e-15`).

This fills the direct paper-algorithm gap that was still missing after the wake
probe: KLS can now retain and update live BTF current state at the right
producer-run boundaries. It is not a speed result yet, because the normal scalar
refactor still executes and the retained state is only diagnostic. The next
refactor step should make this state the numeric owner for those current
columns, avoiding the duplicated scalar prefix/run recomputation and sharing the
producer L stream across the wake-ready batch.

The first guarded numeric consumer for that state is now implemented under
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_EXEC=1`. It stores retained U
coefficients while materializing/advancing compact state, publishes a READY
state only after the run-end advance has completed, arms only the farthest
member per current column, restores that state at dispatch, writes the skipped U
entries, and then returns to the existing scalar suffix, pivot check, and L
storage path. The READY-after-advance ordering matters: the first trial
published READY before advance and produced a `1.53e-8` residual on
`ASIC_100ks`; delaying READY until after advance restored the clean
`1.92251861e-15` residual.

Best-member filtering makes the probe much cheaper but still not competitive.
On `ASIC_100ks`, execution mode materialized `15,124` states / `940,885` rows,
advanced `14,818` states through `10,756,780` L entries, consumed `604`
retained states, skipped `129,480` U dependencies, restored `291,012` state
rows, and refactored in `9.9392s`; the no-flag same-build control refactored in
`0.0475s`. On `ASIC_320ks`, it materialized `16,742` states / `1,168,129`
rows, advanced `16,464` states through `14,290,944` L entries, consumed `998`
states, skipped `193,588` U dependencies, restored `463,173` rows, and
refactored in `8.5472s`; the no-flag control was `0.0808s`. Both execution
probes had zero state-exec rejects and zero state rejects.

The next diagnostic split records why most wake-ready retained states are not
consumed: `refactor_last_btf_scalar_run_group_state_exec_not_ready_currents`
counts dispatch attempts where the selected retained state is not READY yet,
`refactor_last_btf_scalar_run_group_state_exec_dispatch_bypass_currents`
counts target-column dispatches with a selected retained state where the
restore hook was gated off, the matching `_ready_` counter records how often
that bypassed state was already READY,
`refactor_last_btf_scalar_run_group_state_exec_owned_ready_currents` counts
states that become READY while the target current column is already
claimed/leased, and
`refactor_last_btf_scalar_run_group_state_exec_late_ready_currents` counts
states that become READY only after their target current column is already done.
This distinguishes a simple consumer lookup miss from the larger paper-level
scheduling problem of publishing current-state work too late for the existing
scalar EGraph dispatch.

The bypass split makes the dominant large-case miss explicit. On `ASIC_100ks`,
the guarded executor produced `14,818` READY retained states, consumed only
`607`, and saw `14,130` dispatch bypasses where the selected retained state was
already READY; only `1,660` dispatches saw not-ready state, `7` became READY
while already owned, and `74` became READY after done. On `ASIC_320ks`, it
produced `16,441` READY states, consumed `973`, and saw `15,386` READY
bypasses; not-ready, owned-ready, and late-ready were only `1,351`, `18`, and
`64`. Both probes remained residual-clean with zero state rejects. This points
directly at the restore-hook eligibility gate around `plain_scalar_updates`:
the retained BTF state exists, but most target columns bypass the consumer
hook instead of being owned by a grouped current-state executor.

The immediate follow-up removes the `wait_for_dependencies` part of that hook
gate while keeping the plain-scalar and Algorithm 5 seed guards. This makes
already-READY retained state consumable from both wait and no-wait BTF scalar
dispatches. On `ASIC_100ks`, the same guarded executor then consumed `14,697`
states, skipped `302,057` U dependencies, restored `737,313` rows, and had zero
dispatch bypasses, but still refactored in `10.5670s`. On `ASIC_320ks`, it
consumed `16,356` states, skipped `345,690` dependencies, restored `894,007`
rows, and refactored in `8.7802s`. Both probes stayed residual-clean with zero
state rejects. This closes the hook-reachability gap, but it also confirms that
the per-current sparse-state restore is still the wrong CKTSO-gap owner: once
nearly all READY states are consumed, restored-row traffic dominates and the
path remains orders of magnitude slower than the default refactor.

The retained-state advance helper now uses the sorted L-row fact already
available in the LU pointer cache. When the state executor is requested, KLS
records whether L row lists are sorted; if they are, advancing a retained sparse
state walks the sorted retained-state rows and L rows once instead of binary
searching the retained state for every L entry. This is not the missing grouped
producer owner, but it removes a real cost from the sparse-state experiment:
`ASIC_100ks` improved from `10.5670s` to `9.8842s`, and `ASIC_320ks` improved
from `8.7802s` to `8.3343s`, both with clean residuals and zero state rejects.
The remaining gap is still much larger than row-position lookup; the executor
continues to materialize, advance, and restore one sparse current state per
consumer instead of streaming each producer L suffix once across a live group.

This confirms that the retained current-state handoff can be made
correctness-clean, but it also rejects a per-current restored sparse state as
the missing CKTSO-speed mechanism. The next paper-aligned owner needs to keep a
multi-current producer/window batch alive and stream the producer L data once
across that batch, rather than materializing and restoring one sparse state per
eventual current column.

The default BTF EGraph scalar tail now also skips the dormant direct-complete
final-state probe before calling its helper. The helper already returned
immediately unless
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_COMPLETE=1` was
set, but the inactive call still sat in the repeated unscaled BTF column hot
path. Correctness passed `cmake --build build -j2`, `./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_SUPERNODE_ALGORITHM5_PAYOFF_DIRECT_PREFIX_COMPLETE=1 ./build/kls_smoke`,
`ctest --test-dir build --output-on-failure`, and `git diff --check`. Focused
default probes stayed residual-clean: `ASIC_100ks` measured
`refactor_seconds_avg=0.0470649`, `ASIC_320ks` measured
`0.0809238`, and the full 41-row CKTSO-gap manifest completed without failures
at `2.08205s` SPICE-cycle geomean. This is retained as another safe inactive
probe guard, not as a CKTSO-gap explanation; the missing refactor mechanism is
still the grouped live current-state owner described above.

The retained BTF state executor now reports how much scalar suffix remains after
a consumed retained state. The new last/cumulative counters split terminal
currents, remaining suffix dependencies, and remaining suffix L-entry work. This
tests the direct CKTSO Algorithm 5 post-factorization question: if most retained
states are already terminal, should KLS finish those columns directly from the
state map? Focused probes stayed residual-clean, but the counters show that
terminal direct completion is not the large missing piece. On `ASIC_100ks`,
`build/kls_btf_group_state_exec_suffix_asic100ks_t4_r1_ref1.json` consumed
`14,697` retained states, of which `10,085` were terminal, with only `12,566`
remaining suffix dependencies and `390,272` remaining L entries; nevertheless it
still restored `738,957` state rows and refactored in `10.2357s`. On
`ASIC_320ks`,
`build/kls_btf_group_state_exec_suffix_asic320ks_t4_r1_ref1.json` consumed
`16,360` states, `2,728` terminal, with `45,717` remaining dependencies and
`396,029` remaining L entries, while restoring `898,406` rows and refactoring in
`8.7118s`. This confirms the suffix after restore is not the first-order loss;
the paper-level owner still has to avoid materializing/restoring one sparse
state per current and instead keep the multi-current producer/window batch live.

The follow-up implementation now makes those terminal retained states publish
directly from the retained state map instead of copying through the scalar
workspace. The path is still guarded by
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_EXEC=1`; it validates the
retained pivot and L rows, applies the same pivot-reject rule, writes the
retained U coefficients, diagonal, and L multipliers, and uses a merge walk for
sorted retained L rows. The focused probes stayed residual-clean with zero
state-exec rejects. On `ASIC_100ks`,
`build/kls_btf_group_state_exec_terminal_direct_merge_asic100ks_t4_r1_ref1.json`
reduced restored rows from `738,957` to `128,014` and measured `10.1458s`
versus `10.2357s` for the restore path. On `ASIC_320ks`,
`build/kls_btf_group_state_exec_terminal_direct_merge_asic320ks_t4_r1_ref1.json`
reduced restored rows from `898,406` to `434,141` and measured `8.5201s`
versus `8.7118s`. This is worth retaining as a correct Algorithm 5
post-factorization writeback slice, but the timing confirms the same larger
diagnosis: retained-state materialization and producer-run advancement dominate,
so the CKTSO-speed path still needs a live multi-current producer/window owner.

The wake-time terminal publisher now fills the next direct Algorithm 5 gap:
terminal retained states can be claimed and published immediately when the
producer wake fires, even when the future current is still in the
level-synchronous, non-tail-pipeline phase. The first diagnostic version found
that all terminal wake candidates were outside the tail-pipeline position
surface, not blocked by dependency or claim contention. The final guarded
implementation claims those currents through the same `pipeline_claimed` array,
marks them done, and lets the later level worker skip them through the existing
claimed-column wait path. Focused probes stayed residual-clean with zero
state-exec rejects. On `ASIC_100ks`,
`build/kls_btf_group_state_exec_wake_terminal_nonpipeline_asic100ks_t4_r1_ref1.json`
published all `10,104` terminal wake candidates, all nonpipeline, restored
`130,330` rows, and measured `10.0954s`. On `ASIC_320ks`,
`build/kls_btf_group_state_exec_wake_terminal_nonpipeline_asic320ks_t4_r1_ref1.json`
published all `2,684` terminal wake candidates, all nonpipeline, restored
`425,154` rows, and measured `8.5030s`. This closes the obvious terminal
writeback ownership gap but still only moves timing slightly, so the remaining
paper-level loss is earlier: KLS materializes and advances one retained sparse
state per current instead of keeping CKTSO-style grouped current state live
across the producer/window batch.

The retained-state wake path now has an explicit grouped-advance diagnostic,
guarded by
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_ADVANCE_BATCH_STATS=1`, to
check whether several retained states from the same producer group advance
together often enough to justify a direct grouped owner. The diagnostic is
gated separately from
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_EXEC=1`, so normal retained
state timing leaves these counters at zero. Focused normal runs stayed
residual-clean and confirmed the gate: `ASIC_100ks` measured `10.4165s` with
all advance-batch counters at zero, and `ASIC_320ks` measured `8.34268s` with
all advance-batch counters at zero. With the diagnostic enabled, `ASIC_100ks`
reported `1,702` same-group advance batches covering `3,526` currents, only
`46,878` unique producer-stream entries versus `95,062` duplicated per-current
entries, max batch width `3`, and `9.98670s` refactor time. `ASIC_320ks`
reported `981` batches covering `1,967` currents, `215,456` unique entries
versus `455,595` duplicated entries, max batch width `3`, and `8.32107s`.
Both diagnostic probes had clean residuals and zero state-exec rejects. This
shows a real pairwise sharing surface, but it is too small to explain the slow
cases by itself: batched duplicate entries were about `0.8%` of all advanced
entries on `ASIC_100ks` and `3.3%` on `ASIC_320ks`. The larger paper gap is
therefore not simply "batch adjacent wake advances"; it remains the need for a
different live current-state owner that avoids materializing, advancing, and
restoring one sparse state per current across the whole producer/window region.

A workspace lazy-clear variant was also rejected. Removing the full sequential
clear from retained BTF state preparation kept residuals clean but moved cost
into scattered worker-side initialization: `ASIC_100ks` slowed from the
`10.0954s` wake-terminal baseline to `10.2738s`, and `ASIC_320ks` slowed from
`8.5030s` to `8.59457s`. KLS therefore keeps the full clear because it likely
pre-touches the retained-state workspace more effectively than piecemeal lazy
initialization on these cases.

The BTF scalar-run group descriptor now also reports producer-step fanout for
the live grouped-current owner the papers imply. This is different from the
wake-adjacent batch diagnostic above: for each producer row inside a grouped
contiguous run, KLS counts how many current columns would be updated if a
completed-producer hook streamed that producer's L rows once across all live
current states. The counters are descriptor-level and remain gated by
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUPS=1`, so the default no-env run leaves
them at zero. Focused controls stayed residual-clean:
`build/kls_btf_group_producer_step_off_asic100ks_t4_r1_ref1.json` reported no
group descriptor and no producer-step counters. With groups enabled,
`ASIC_100ks` reported `76,244` producer steps, `68,467` multi-current steps,
`1,225,700` active current-step slots, `9,319,210` unique producer L entries,
`308,891,777` duplicate per-current entries, `299,572,567` reusable entries,
max fanout `362`, and a clean residual. `ASIC_320ks` reported `115,602`
producer steps, `98,964` multi-current steps, `1,139,656` active slots,
`4,398,767` unique producer L entries, `274,928,002` duplicate entries,
`270,529,235` reusable entries, max fanout `444`, and a clean residual. The
duplicate/unique ratios are about `33.1x` and `62.5x`, respectively. This
confirms the prior small wake-adjacent result was only measuring a narrow
scheduling coincidence; the real CKTSO-sized surface is a producer-row indexed
live grouped-current executor that streams each completed producer once across
all active currents.

The follow-up refactor chunk now retains that producer-row indexed schedule in
the BTF scalar-run group descriptor. For every indexed completed producer, KLS
stores the dependent group steps and the active grouped-current members for each
step; the descriptor remains tied to the LU pointer cache and is still gated by
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUPS=1`. This is not yet the numeric
executor, but it removes the next structural gap between the diagnostic fanout
and a real completed-producer hook. Focused probes stayed residual-clean. The
default-off `ASIC_100ks` control left the group/index counters at zero and
measured `0.0499743s`. With the descriptor enabled, `ASIC_100ks` built `34,328`
producer index entries for `76,244` producer steps, with at most `24` group
steps on one completed producer, and measured `0.0437695s`. `ASIC_320ks` built
`82,133` producer index entries for `115,602` producer steps, with at most `11`
group steps on one producer, and measured `0.0838740s`. The duplicate/unique
producer-entry ratios remain the key signal: about `33.1x` on `ASIC_100ks` and
`62.5x` on `ASIC_320ks`. The next paper-aligned implementation step is to make
the runtime producer-completion path consume this retained schedule while
keeping the grouped current workspaces live, instead of materializing/restoring
one sparse current state per eventual consumer.

The retained schedule is now consumed by an opt-in producer-step retained-state
advance path,
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_STEP_ADVANCE=1`. Instead of
waiting until the last producer in a contiguous run and then replaying the whole
run for one retained sparse current state, KLS tracks each materialized
membership's next producer row and advances it from the completed-producer hook
through the producer-step index. The path uses the existing wake-state CAS for
per-member ownership, can finish a state at run-end without replaying already
advanced rows, and remains default-off. Correctness stayed clean:
`./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_STEP_ADVANCE=1 ./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_STEP_ADVANCE=1 KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_EXEC=1 ./build/kls_smoke`,
and `ctest --test-dir build --output-on-failure` all passed.

The focused ASIC probes show this is a useful paper-aligned substrate but not
the CKTSO-gap closer. With step advance but no retained-state execution,
`ASIC_100ks` stayed residual-clean and consumed the whole producer-step index:
`34,328` triggers, `76,244` steps, `396,292` advanced current events,
`479,962` advanced rows, `125,069,484` L entries, `44,115` ready retained
states, and zero rejects. It still measured `11.6130s` because it materialized
`44,115` sparse current states with `13,361,860` retained rows. `ASIC_320ks`
was similar: `82,133` triggers, `115,602` steps, `440,257` advanced current
events, `541,432` rows, `144,652,075` entries, `49,635` ready states, zero
rejects, and `10.2718s`, after materializing `15,281,862` retained rows. With
`STATE_EXEC=1` on `ASIC_100ks`, step advance reduced the retained work surface
to `15,284` materialized states, `900,370` retained rows, `75,626` advanced
rows, and `11,390,181` step-advance entries, while consuming `15,193` states
and restoring only `133,941` rows. The run still measured `10.5671s`, about the
same order as the prior wake-terminal retained-state executor. This rejects
"incrementally advance one sparse retained current state per member" as the
large missing CKTSO mechanism. The next implementation should use the same
producer-step schedule to own a grouped live workspace and stream each producer
L column once across multiple current states, rather than maintaining thousands
of separate sparse retained states.

The producer-step path now has that first producer-column batched slice: for a
completed producer, KLS claims all eligible retained BTF scalar-run memberships
across the producer's retained step range, captures each current state's U
coefficient, and streams the producer L column once across the claimed retained
states. The scalar per-member helper remains as a fallback and catch-up path.
Correctness stayed clean on `./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_STEP_ADVANCE=1 ./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_STEP_ADVANCE=1 KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_EXEC=1 ./build/kls_smoke`,
and `ctest --test-dir build --output-on-failure`. The default `ASIC_100ks`
control still left the step and batch counters at zero and measured
`0.0474679s`.

The focused probes show the batching is active and useful substrate, but still
not the large CKTSO gap closer. On `ASIC_100ks`, the producer batch covered
`60,384` retained steps, `314,957` current updates, and `65,762,438` update
entries with zero rejects; refactor time moved from the previous step-only
`11.6130s` to `11.2642s`, with a clean `1.92e-15` relative residual. On
`ASIC_320ks`, it covered `61,619` retained steps, `332,333` current updates,
and `72,724,209` entries with zero rejects; refactor time moved from the
previous step-only `10.2718s` to `9.78565s`, with a clean `2.08e-15` relative
residual. This directly fills the "stream one producer across multiple current
states" surface, but the retained-state materialization remained enormous
(`13,493,577` rows on `ASIC_100ks`, `15,218,280` on `ASIC_320ks`). The next
paper-aligned gap is therefore not the producer trigger itself; it is replacing
the sparse retained-state objects with a grouped live workspace/supernodal
window owner so the producer batch does not pay per-current retained row maps.

Rechecking the selected-large timeout cases keeps `pre2` as the largest
actionable KLS-vs-CKTSO gap. The saved large recon sidecars show KLS timing out
on `pre2` and `Hamrle3` under 120s, while CKTSO times out only on `Hamrle3`;
the same-session CKTSO `pre2` rerun completes the full compare in `21.50s`
wall time with `initial_factor_seconds=6.210766` and
`refactor_seconds_avg=4.922633`. KLS now makes the dominant-BTF parallel
KLS-first path fall back to the same legacy count-based separator queue used by
the serial path when a METIS partitioned queue fails private-ownership
validation. On forced-METIS `pre2`, the 60s trace improved from the old
131,072-row checkpoint with about `4.06e9` scalar published-U entry touches to
about `1.75e8` touches at the same checkpoint and reached the pivot-tail
restart near row `593,557`. However,
`build/kls_pre2_metis_fallback_t4_r1_ref0_timeout130.json` remained empty
after a non-traced 130s factor-only timeout, and the traced 130s run was still
inside the 629,628-row pipeline after the 589,824-row checkpoint. This confirms
the timeout policy is not the issue; the remaining clear missing piece is the
paper-level coarse row/supernode first-factor executor for the dominant block,
especially the scalar dependency drain and pivot-tail continuation.

A July 1, 2026 timeout-pair recheck keeps that diagnosis. CKTSO completed
`pre2` in `build/cktso_recheck_timeout_pair_t4_r1_ref1_timeout120.jsonl` with
`analysis_seconds=4.584171`, `initial_factor_seconds=7.237705`,
`factor_seconds_avg=4.191527`, `refactor_seconds_avg=5.262766`, and
`solve_seconds_avg=0.138894`, while `Hamrle3` timed out. The matching KLS run
`build/kls_recheck_timeout_pair_t4_r1_ref1_timeout120.jsonl` emitted no rows
because both `pre2` and `Hamrle3` timed out. A forced-METIS KLS-first `pre2`
trace reached the 589,824-row checkpoint with `804,575,231` scalar published-U
entry touches, then spent the cap on rows `598596`-`598872`; the 15 committed
rows over the 5M trace threshold were all dynamic-pivot rows and summed
`84,654,869` scalar U-entry touches. An attempted serial pivot-storm drain was
not retained: it triggered after only 5,024 committed rows and still timed out
at 120s. The clear missing mechanism is therefore not a timeout limit, small
BLAS threshold, or simple pivot-storm serial fallback; it is still the
CKTSO/SubtreeLU-style coarse row/supernode numeric executor that avoids
replaying long published-U streams for adjacent pivoting rows.

KLS now keeps row-first cached supernode panels across dynamic column pivots
when both exchanged columns are only trailing columns of the panel, while still
invalidating panels whose dense block would be changed. This fills a direct
paper-aligned retention gap: a tail-only column exchange does not invalidate
the supernode's triangular block, so the cached panel can continue to serve
later row updates after its common tail column labels are exchanged. The target
`pre2` result is positive but partial. A 90s forced-METIS KLS-first long-row
trace with this change,
`build/kls_pre2_panel_tail_exchange_trace_t4_r1_ref0_timeout90.stderr`,
reached long rows through about `599710` instead of about `598872`, and the
long-row sample's panel-backed update rows rose from `2,940` to `92,353`.
The low-noise 130s trace
`build/kls_pre2_panel_tail_exchange_trace_t4_r1_ref0_timeout130.stderr` still
timed out at the 589,824-row checkpoint, so the remaining gap is still the
coarser CKTSO/SubtreeLU numeric executor rather than this one cache-retention
detail. Validation included `./build/kls_smoke`,
`ctest --test-dir build --output-on-failure`, a no-failure top-five CKTSO-gap
medium run at `1.6325s` geomean, and a completed forced KLS-first
`transient` pivot-plus-panel run with `2.75e-13` relative residual.

A same-binary timeout-pair refresh after the tail-only panel exchange confirms
that `pre2` remains the largest clean KLS-vs-CKTSO timeout gap. CKTSO completed
`pre2` in `build/cktso_timeout_pair_refresh_t4_r1_ref1_timeout120.jsonl` with
`analysis_seconds=4.347466`, `initial_factor_seconds=8.525571`,
`factor_seconds_avg=4.311943`, `refactor_seconds_avg=3.348520`, and
`solve_seconds_avg=0.134543`, while `Hamrle3` still timed out. The matching KLS
pair `build/kls_timeout_pair_refresh_t4_r1_ref1_timeout120.jsonl` emitted no
rows because both `pre2` and `Hamrle3` timed out. A fresh forced-METIS
KLS-first `pre2` trace,
`build/kls_pre2_timeout_pair_refresh_trace_t4_r1_ref0_timeout90.stderr`,
entered the 629,628-row dominant BTF block with separator coverage and reached
the 589,824-row checkpoint before the 90s cap. At that checkpoint it had
`8,252,845` scalar dependencies and `804,367,990` scalar published-U entry
touches; only `114,634,524` were internal dependency-row touches, while
`689,733,466` were current-row output/trailing touches. The same-run GDB sample
in `build/pre2_timeout_pair_refresh_gdb_run.txt` captured the active pipeline
worker in `kls_row_first_partial_apply_one_dep()` called from
`kls_row_first_partial_apply_ready()`, while the peer pipeline workers waited
on the row-pipeline condition/mutex and the parent joined the pipeline phase.
This keeps the next implementation target at the paper-level coarse
row/supernode numeric executor for producer-to-current output streaming, not
timeout policy, `Hamrle3`, ordering alone, or BLAS dispatch.

An exact current-row output-suffix grouping experiment was rejected as too
narrow for the `pre2` pivot tail. The prototype split scalar dependency rows
into immediate internal updates and a deferred output suffix, then grouped
adjacent ready dependencies only when both rows had the same large ordered
`col >= current_row` suffix. It passed `./build/kls_smoke` and
`ctest --test-dir build --output-on-failure` after adding a pivot-safe guard
that makes unsorted dynamic-pivot rows fall back to scalar. The target trace
`build/kls_pre2_output_run_guarded_trace_t4_r1_ref0_timeout90.stderr` showed
the path activating, but only for `980` groups, `3,841` rows, and `2,161,382`
U entries by the 589,824-row checkpoint. The same checkpoint moved the wrong
way versus `build/kls_pre2_timeout_pair_refresh_trace_t4_r1_ref0_timeout90.stderr`:
scalar published-U touches increased from `804,367,990` to `804,557,892`, and
output touches increased from `689,733,466` to `689,906,083`. This rules out
a single-current exact-suffix accumulator as the clear missing paper mechanism.
The first-factor gap needs a broader producer-indexed or multi-current numeric
owner that reuses producer output streams across several current rows despite
dynamic pivot ordering, rather than another exact-tail cache inside one current
row.

A fresh timeout-pair recheck after rejecting the exact-suffix experiment again
isolates `pre2` as the largest clean KLS-vs-CKTSO timeout gap. CKTSO completed
`pre2` in `build/cktso_timeout_recheck_pair_t4_r1_ref1_timeout120.jsonl` with
`analysis_seconds=4.272073`, `initial_factor_seconds=8.18128`,
`factor_seconds_avg=5.73322`, `refactor_seconds_avg=4.233533`,
`solve_seconds_avg=0.146472`, and `spice_cycle_seconds=446.22032`; CKTSO timed
out only on `Hamrle3`. The matching KLS pair
`build/kls_timeout_recheck_pair_t4_r1_ref1_timeout120.jsonl` emitted no rows,
with both `pre2` and `Hamrle3` recorded in the failure sidecar as 120s
timeouts. A short-start probe reached the forced-METIS KLS-first dominant BTF
row-pipeline start for `pre2` at about `16.5s`, so KLS setup/ordering is slower
than CKTSO's analysis but is not the full 120s failure. The 90s factor-only
trace `build/kls_pre2_timeout_recheck_trace_t4_r1_ref0_timeout90.stderr`
entered the same 629,628-row dominant BTF block and reached the 589,824-row
checkpoint with `8,252,734` scalar dependencies, `804,607,861` scalar U-entry
touches, `114,660,977` internal touches, and `689,946,884` output/trailing
touches. The gap is therefore still dominated by first-factor row-pipeline
published-U/output streaming in the dominant BTF block; the actionable missing
mechanism remains a CKTSO/SubtreeLU-style coarse producer-to-multiple-current
numeric owner for pivoting rows, with setup/ordering as a secondary gap.

KLS now has the first safe slice of that producer-owned numeric path in the
first-factor row pipeline: after a row commits and publishes its U row, the
commit cursor can claim waiting partial rows whose next scalar dependency is
that producer, then stream the producer U row once across the claimed current
rows. The claim is guarded by the pipeline order epoch, requires the producer
to be the heap root for each target row, and skips targets where the existing
supernode or cached-panel executor can own a wider run. The path is on by
default only for coarse producers with at least `1024` saved stream entries,
and can be disabled with `KLS_DISABLE_ROW_PIPELINE_PRODUCER_BATCH=1`. The
unthresholded version was too eager: on the forced KLS-first spot
`G2_circuit`/`transient` run it measured `10.7099s` geomean versus `10.1031s`
with batching disabled. A `512` saved-stream threshold was still slightly
negative (`12.0067s` versus `11.8314s`). The retained `1024` threshold is
neutral on that forced spot (`11.2926s` versus `11.3075s`) with clean
residuals. On `pre2`, the retained threshold activates only for `340` coarse
producer batches by the 589,824-row checkpoint in
`build/kls_pre2_producer_batch_threshold_trace_t4_r1_ref0_timeout90.stderr`:
it replaces `585,522` target U updates with `195,174` producer-stream reads,
and scalar published-U touches move from `804,607,861` in
`build/kls_pre2_timeout_recheck_trace_t4_r1_ref0_timeout90.stderr` to
`804,044,491`. This is a correct paper-aligned substrate, but it also proves
the current pipeline's live current-row window is too small to capture the
large missing reuse surface; the next gap-closing step is a broader live-row
window or grouped current-state owner, not lowering the threshold again.

An experimental bounded live-row lookahead substrate now exists, but it is not
ready for default use. It is gated behind both
`KLS_ENABLE_EXPERIMENTAL_ROW_PIPELINE_LOOKAHEAD=1` and
`KLS_ROW_PIPELINE_LOOKAHEAD=<slots>` so ordinary runs, and even accidental
`KLS_ROW_PIPELINE_LOOKAHEAD` settings, stay on the safe worker-only pipeline.
The substrate reserves extra future row states, initializes them to the current
pipeline epoch, lets committed producers update them through the producer-batch
path, and lets workers adopt the oldest unclaimed state later. This directly
tests the broader live-current-window idea that the `pre2` trace called for.
The result is mixed and therefore intentionally not default: with the
experimental path effectively enabled during development, `G2_circuit`
improved from the prior forced KLS-first `13.7s` range to `9.20s`-`9.79s`,
but `transient` timed out at 120s with both 2 and 8 lookahead slots. The normal
guarded build still passes `./build/kls_smoke`, `ctest --test-dir build
--output-on-failure`, and the forced KLS-first `G2_circuit`/`transient` spot
completed at `8.3534s` geomean in
`build/kls_lookahead_guard_default_forced_spot_t4_r1_ref0_timeout120.jsonl`.
This confirms the live-row window is a real paper-aligned lever, but the
current implementation needs eligibility control or cheaper catch-up before it
can be used generally.

The lookahead substrate now has that first eligibility control: by default it
only reserves rows with at most `4` raw input entries
(`KLS_ROW_PIPELINE_LOOKAHEAD_MAX_INPUTS=0` disables this cap), and it prepares
lookahead rows by loading only the raw input pattern under the commit mutex.
Older dependency catch-up is left to the adopting worker, so the lookahead
window no longer performs large scalar drains while publishing a producer.
This made the experimental smoke check complete under a 30s cap, but it did
not make the policy generally safe: with
`KLS_ENABLE_EXPERIMENTAL_ROW_PIPELINE_LOOKAHEAD=1` and
`KLS_ROW_PIPELINE_LOOKAHEAD=8`, the forced KLS-first spot still completed
`G2_circuit` at `10.9122s` but timed out `transient` at 120s in
`build/kls_lookahead8_inputonly_forced_spot_t4_r1_ref0_timeout120.jsonl`.
The safe default path remained unaffected: the same forced spot without the
experimental flag completed both matrices at `8.9437s` geomean in
`build/kls_lookahead_inputonly_default_forced_spot_t4_r1_ref0_timeout120.jsonl`,
and `KLS_ROW_PIPELINE_LOOKAHEAD=8` alone stayed inert. The remaining problem is
therefore not just expensive lookahead preparation; it is that the current
lookahead adoption policy can reserve or prioritize rows that are bad for some
pipelines. The next implementation should add shape/owner eligibility before
adoption, or make lookahead rows advisory rather than consuming the global row
cursor.

The lookahead path has now been made advisory rather than cursor-consuming:
workers first claim the normal global row position, then adopt a prepared
lookahead state only when it exactly matches that position. This fixes the
pathological scheduling failure but does not make lookahead a default win. The
explicit experimental run
`build/kls_lookahead8_advisory_forced_spot_t4_r1_ref0_timeout120.jsonl`
completed both forced KLS-first spot matrices with clean residuals
(`G2_circuit=10.9588s`, `transient=12.2114s`), while the guarded default
artifact
`build/kls_lookahead_inputonly_default_forced_spot_t4_r1_ref0_timeout120.jsonl`
remained faster (`G2_circuit=10.6131s`, `transient=7.5370s`). The result keeps
lookahead as an experimental substrate: exact-position advisory adoption is
safer, but it still needs a stronger owner/shape filter before it can help the
large timeout cases.

A current timeout-case recheck keeps the largest clean KLS-vs-CKTSO gap on
`pre2`, and narrows it to numeric first factor rather than the benchmark
refactor loop. The direct KLS factor-only probes
`build/kls_pre2_factoronly_current_t4_r1_ref0_timeout120.json`,
`build/kls_pre2_factoronly_metis_current_t4_r1_ref0_timeout120.json`, and
`build/kls_pre2_factoronly_scotch_current_t4_r1_ref0_timeout120.json` all
timed out at 120s without JSON, while the direct `Hamrle3` factor-only probe
also timed out. CKTSO's refreshed workspace artifact
`build/cktso_timeout_pair_refresh_t4_r1_ref1_timeout120.jsonl` completed
`pre2` with `analysis_seconds=4.347466`,
`initial_factor_seconds=8.525571`, `factor_seconds_avg=4.311943`,
`refactor_seconds_avg=3.348520`, and `solve_seconds_avg=0.134543`; CKTSO only
timed out on `Hamrle3`. KLS analyze-only on `pre2` is slower but not the main
120s failure: AMD took `21.398362s`, METIS took `18.6061678s`, and SCOTCH took
`20.9142238s`, and all three retained the same `629628`-row dominant block.
The current 75s traced AMD factor-only run
`build/kls_pre2_factoronly_trace_current_t4_r1_ref0_timeout75.stderr` reached
the dominant row-pipeline start but no 65,536-row progress checkpoint; the
thread `/proc` sample in `build/kls_pre2_proc_sample_current.txt` showed the
process had entered the threaded factor phase after setup. This makes ordering
a secondary issue for this case: METIS and SCOTCH did not rescue the timeout.
The actionable largest gap remains the first-factor row-pipeline executor for
the dominant BTF block, especially the published-U/output streaming and
producer-to-multiple-current reuse that the CKTSO/SubtreeLU-style algorithms
avoid doing one current row at a time.

The experimental lookahead path now has a producer-stop catch-up slice. Before
the completed-producer batch scan tests an unclaimed lookahead state, KLS
advances that partial row through ready dependencies strictly before the
producer and stops with the producer still at the heap root. This directly
matches the paper-level producer-owned update idea: the prefix is made ready,
but the just-published U row is still streamed once across the batch rather
than scalar-applied independently per current row. The effect is visible on
`pre2`: at the 589,824-row checkpoint, the no-lookahead retained producer
batch trace `build/kls_pre2_producer_batch_threshold_trace_t4_r1_ref0_timeout90.stderr`
had `340` producer batches, `1020` targets, `195174` streamed producer U
entries, and `804044491` scalar U-entry touches. With
`KLS_ENABLE_EXPERIMENTAL_ROW_PIPELINE_LOOKAHEAD=1` and
`KLS_ROW_PIPELINE_LOOKAHEAD=8`,
`build/kls_pre2_lookahead8_catchup_trace_t4_r1_ref0_timeout90.stderr`
increased that to `1967` producer batches, `11617` targets, and `704485`
streamed producer U entries, while scalar U-entry touches fell to
`802908719`. That is real movement in the intended direction but still far
too small for the CKTSO gap. Pushing to `32` lookahead slots grew the target
surface further (`3652` batches and `32104` targets by the 524,288-row
checkpoint in
`build/kls_pre2_lookahead32_catchup_trace_t4_r1_ref0_timeout90.stderr`) but
lost row progress by the same 90s cap, so the next step should not simply
increase slots. The forced KLS-first `G2_circuit`/`transient` spot stayed
residual-clean with the 8-slot catch-up path but was not a robust runtime win:
`build/kls_lookahead8_catchup_forced_spot_t4_r1_ref0_timeout120.jsonl`
measured `11.2465s` geomean, while the exact-code rerun
`build/kls_lookahead8_catchup_rerun_forced_spot_t4_r1_ref0_timeout120.jsonl`
measured `12.0360s`, both slower than the guarded default `8.9437s` geomean.
This keeps catch-up lookahead experimental and points to a persistent
grouped-current owner, not a larger pool of independent future-row states, as
the remaining paper-aligned mechanism.

A follow-up producer-directed raw-input scan was tested and rejected before it
was retained in source. The idea was to fill experimental lookahead slots with
future rows whose raw input already referenced the just-completed producer,
then run the existing producer-batch update on those states. This is still not
the needed grouped-current owner. In
`build/kls_pre2_lookahead8_producerdirect_trace_t4_r1_ref0_timeout90.stderr`,
the run reached only the 393,216-row checkpoint by the same 90s cap, versus
589,824 rows for
`build/kls_pre2_lookahead8_catchup_trace_t4_r1_ref0_timeout90.stderr`.
At that lower checkpoint it had `1132` producer batches and `6095` targets,
so it did find some producer-specific states, but the scan/fill policy spent
too much time on poorer future rows and reduced overall row progress. This
rejects raw-input future-row scanning as the next CKTSO-gap closer; the next
implementation should keep a producer-indexed grouped-current state live
without repeatedly searching future rows or materializing many independent
full workspaces.

A fresh timeout recheck on the current tree shifts the immediate `pre2`
culprit earlier inside that first-factor pipeline. A direct CKTSO run on the
same workspace completed `pre2` with `analysis_seconds=5.304347`,
`initial_factor_seconds=9.175596`, `refactor_seconds_avg=6.084399`, and
`solve_seconds_avg=0.140293`. KLS analyze-only on the same command shape
completed with `analysis_seconds=8.33114405` in
`build/kls_pre2_timeout_recheck_analyze_current_t4.json`, so ordering/analyze
is not enough to explain the 120s failure. The current KLS traced factor run
`build/kls_pre2_timeout_recheck_trace_current_t4_r1_ref0_timeout120.stderr`
entered the `629628`-row dominant row-pipeline block but emitted no 65,536-row
progress checkpoint before the timeout. A debugger-owned interrupt of the same
AMD/auto command
`build/pre2_timeout_recheck_auto_gdb_interrupt.txt` sampled one active worker
inside `kls_row_first_supernode_panel_cache_build()` /
`kls_row_first_supernode_panel_cache_free()` while the other row-pipeline
workers waited on the pipeline mutex/condition. This makes the largest current
gap more specific than generic scalar published-U streaming: early first-factor
prefix/pivot panel-cache rebuild work is serialized under the row-pipeline
commit lock before useful row progress. The next fix should remove or sharply
bound those early panel-cache rebuilds on large dominant BTF blocks, or make
the rebuild/exchange state incremental and outside the global pipeline mutex,
before returning to broader grouped-current producer ownership.

The direct large-prefix-cache bypass prototype was tested and rejected. The
first version skipped full prefix panel-cache rebuilds for all blocks above
`32768` rows after a dynamic pivot and only cleared active cached panels. It
remained residual-clean on `ASIC_320k`, but it pushed that probe from the
saved `~2s` initial-factor range to `11.2391421s` in
`build/kls_asic320k_large_prefix_clear_t4_factor.json`, so the guard was too
broad. Narrowing the bypass to blocks at least `16 * 32768` rows left
`ASIC_320k` below the cutoff but still did not improve the hard case:
`build/kls_pre2_large_prefix_clear_trace_t4_r1_ref0_timeout120.stderr` reached
only the row-pipeline start. The matching interrupt
`build/pre2_large_prefix_clear_gdb_interrupt.txt` showed the sampled worker had
moved from panel-cache build/free into `kls_row_first_supernodes_reset_prefix`,
still under the row-pipeline lock. A stronger huge-block variant then disabled
prefix supernodes immediately after the pivot and skipped wait-time dependency
drains under the lock. It passed `ctest --test-dir build --output-on-failure`
and kept the sensitive `rajat29` no-fast probe residual-clean
(`build/kls_rajat29_large_prefix_disable_nofast_t4_factor_timeout90.json`),
but `build/kls_pre2_large_prefix_disable_nowait_trace_t4_r1_ref0_timeout75.stderr`
still did not reach the first 65,536-row checkpoint. The short low-threshold
long-row trace
`build/kls_pre2_large_prefix_disable_nowait_longrow_t4_timeout30.stderr`
proved rows were committing below that coarse interval (`completed=21520` by
the cap), but many committed rows were still applying roughly `3.3M` scalar U
entries apiece. The last interrupt
`build/pre2_large_prefix_disable_nowait_gdb_interrupt.txt` sampled the active
worker in `kls_row_first_partial_apply_ready_until()` /
`kls_row_first_partial_apply_one_dep()` with peer workers waiting. The source
was restored. This sequence rules out cache rebuild/reset alone as the
CKTSO-scale gap closer; after that serialization is removed, the same `pre2`
region is dominated by repeated scalar dependency/output application. The next
retained implementation needs the paper-level grouped-current/producer-output
executor, not another prefix-cache invalidation shortcut.

A rebuilt-binary timeout recheck again isolates `pre2` as the largest clean
KLS-vs-CKTSO timeout gap. With `repeat=1`, `refactor-repeat=1`, 4 threads, and
the same 120s wall cap over the known timeout pair, CKTSO completed `pre2` in
`build/cktso_timeout_pair_current_t4_r1_ref1_timeout120.jsonl` with
`analysis_seconds=4.652263`, `initial_factor_seconds=6.003679`,
`refactor_seconds_avg=3.714327`, `solve_seconds_avg=0.133762`, and
`spice_cycle_seconds=391.750515`; CKTSO timed out only on `Hamrle3`. KLS
emitted no successful row in
`build/kls_timeout_pair_current_t4_r1_ref1_timeout120.jsonl`; the sidecar shows
both `pre2` and `Hamrle3` timed out. Analyze-only is not the largest part:
fresh `pre2` analyze-only completed in `9.72s` wall / `8.05047587s` reported
analysis for auto ordering
(`build/kls_pre2_timeout_pair_current_analyze_auto_t4.json`) and `5.46s` wall /
`3.80298773s` reported analysis for explicit AMD
(`build/kls_pre2_timeout_pair_current_analyze_amd_t4.json`). A forced KLS
first-factor trace reached the dominant `629628`-row BTF row-pipeline block
but not the first 65,536-row progress checkpoint before a 75s cap
(`build/kls_pre2_timeout_pair_current_forced_trace_t4_r1_ref0_timeout75.stderr`).
The only committed progress shown there was around rows `14054` and `14070`,
with individual rows applying about `1.05M-1.12M` scalar U entries and
roughly `0.90M-0.96M` of those as output entries. A gdb-owned rerun
(`build/pre2_timeout_pair_current_gdb_run_under.txt`) sampled one pipeline
worker in `kls_row_first_supernode_panel_cache_build()` while peer row-pipeline
workers waited on the condition variable and the outer KLS first-factor worker
joined that phase. The largest clean gap is therefore KLS initial numeric
factorization of the dominant `pre2` BTF block: first serialized panel-cache
rebuild/reset work in the row-pipeline commit path, then repeated scalar
producer-output replay across adjacent current rows. The timeout limit and
ordering/analyze phase are secondary for this case.

The July 1, 2026 timeout recheck also rejects a narrower same-output-tail
batching attempt. The prototype grouped consecutive ready dependency rows only
when their published-U tail columns at or beyond the current row were exactly
identical and long enough to amortize one output scatter. It built cleanly,
passed `ctest --test-dir build --output-on-failure`, and stayed residual-clean
on a forced KLS-first `transient` factor probe, but it did not touch the active
`pre2` loss. The 75s trace
`build/kls_pre2_output_tail_group_trace_t4_r1_ref0_timeout75.stderr` reached
only row `22501` of the `629628`-row dominant BTF block, versus `22419` in the
saved baseline `build/kls_pre2_gaprefresh_trace_t4_r1_ref0_timeout75.stderr`.
Aggregated long-row counters recorded just `45` grouped scalar runs, `102`
grouped rows, and `192701` grouped U entries, while the same long rows still
performed `158704950` scalar U-entry touches and `138095349` scalar
trailing/output touches. The hot long-row lines at the timeout tail still had
`scalar_runs=0`. The source was restored, and the evidence narrows the missing
paper-level owner again: it cannot require identical current-row output tails.
It needs to be a CKTSO/SubtreeLU-style producer/panel-to-many-current executor
that owns nonidentical sparse current states and streams published output work
once across that live group.

The producer-batch miss counters now make that conclusion more concrete. The
opt-in `KLS_TRACE_ROW_PIPELINE=1` long-row trace reports producer probe counts,
ready roots, underfilled probes, and rejection reasons. In the 75s default
`pre2` trace
`build/kls_pre2_producer_miss_default_trace_t4_r1_ref0_timeout75.stderr`, KLS
again reached only row `22540` of the `629628`-row dominant BTF block. The
logged long rows applied `413850748` scalar U entries, including `362332719`
scalar output entries. Producer batching fired, but its scale was tiny:
`112` batches, `287` targets, `210797` streamed producer U entries, and
`539981` target U entries. The rejection split rules out out-of-order
dependency handling as the first fix: `producer_reject_not_root=7`,
`producer_reject_not_ready=0`, `producer_reject_supernode=0`, and
`producer_reject_cached_panel=0`. The bigger misses were absent or unavailable
states: `producer_reject_bad_state=390`, `producer_reject_dep_absent=120`,
and `192` underfilled probes from `912` worker probes. The eight-slot
lookahead check
`build/kls_pre2_producer_miss_lookahead8_trace_t4_r1_ref0_timeout75.stderr`
did not change the conclusion. It reached only row `22501`, and although it
added `841` lookahead probes, it still applied only `36` producer batches and
`284012` target U entries against `138562400` scalar output entries in the
logged long rows; `producer_reject_bad_state=903` dominated the lookahead
misses. The next paper-aligned implementation should therefore create a
persistent producer-indexed live-current owner, with current states retained
and advanced outside the commit-lock producer scan, rather than increasing
lookahead slots or trying to apply producer dependencies out of order.

A July 1, 2026 producer-indexed lookahead prototype tested the most direct
small version of that idea and was rejected before commit. It built a
block-local reverse map from each completed producer row to future row
positions whose raw input contained that producer, then filled lookahead
states from that producer list before the normal producer-batch scan. The path
stayed residual-clean on forced KLS-first `Freescale/transient`
(`build/kls_transient_producer_index_t4_r1_ref0_timeout120.json`), but it did
not move the hard timeout case. The 75s `pre2` trace
`build/kls_pre2_producer_index_trace_t4_r1_ref0_timeout75.stderr` reached only
row `22535` of the `629628`-row dominant BTF block, compared with `22540` for
the default trace, and logged `454599862` scalar U touches with `398381325`
scalar output touches. Producer batching remained too small: `114` batches,
`277` targets, and `520283` target U entries. A same-session CKTSO rerun,
`build/cktso_pre2_recheck_producer_index_t4_r1_ref1_timeout120.json`,
completed the same matrix in `21.15s` wall time, with
`initial_factor_seconds=6.157833` and `relative_residual_l2=1.01703404e-16`.
This narrows the retained implementation target: KLS needs the paper-level
live grouped-current numeric owner that keeps several current states in one
producer-owned batch and streams output work once across that batch. A
producer-to-future-row index feeding the existing per-current sparse states is
not enough.

A retained-state producer-step row-union diagnostic now measures the live
workspace opportunity directly instead of inferring it from producer-entry
counts. The new JSON fields
`refactor_last_btf_scalar_run_group_state_step_batch_state_rows`,
`refactor_last_btf_scalar_run_group_state_step_batch_unique_state_rows`, and
their cumulative counterparts are populated only when
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_ADVANCE_BATCH_STATS=1` is used
with retained producer-step advance. They merge the sorted sparse retained row
sets for the members that actually advanced in each producer-step batch and
count duplicate state rows versus exact unique global rows. The code builds
cleanly, and `ctest --test-dir build --output-on-failure` still passes both
tests.

The diagnostic confirms that the missing mechanism is a grouped live workspace,
not another small threshold around the current per-member sparse state. On
`ASIC_100ks`, `build/kls_btf_group_state_step_union_asic100ks_t4_r1_ref1.json`
was residual-clean (`relative_residual_l2=2.65029351e-15`) and reported
`13016768` duplicate producer-step batch state rows but only `3241691` unique
global state rows, a `4.02x` collapse opportunity. On `ASIC_320ks`,
`build/kls_btf_group_state_step_union_asic320ks_t4_r1_ref1.json` was also
residual-clean (`2.08436857e-15`) and reported `61569635` duplicate state rows
versus `6419228` unique rows, a `9.59x` opportunity. The same runs remained
slow because the opt-in retained-state experiment still materializes tens of
millions of per-member sparse states (`69160671` on `ASIC_100ks` and
`56280350` on `ASIC_320ks`) before the diagnostic can stream producer batches.

A fresh same-tree `pre2` factor-only recheck keeps the largest clean timeout
gap in the same place. The command writing
`build/kls_pre2_factoronly_rowunion_tree_t4_r1_ref0_timeout120.json` timed out
at the 120s cap with no JSON output. The current saved CKTSO pair artifact,
`build/cktso_timeout_pair_current_t4_r1_ref1_timeout120.jsonl`, completed
`pre2` with `analysis_seconds=4.652263`, `initial_factor_seconds=6.003679`,
`refactor_seconds_avg=3.714327`, and `solve_seconds_avg=0.133762`, while
`build/kls_timeout_pair_current_t4_r1_ref1_timeout120.failures` records KLS
timeouts on both `pre2` and `Hamrle3`. Analyze-only KLS is not large enough to
explain the miss (`3.80298773s` on forced AMD ordering). The next
paper-aligned implementation target should therefore replace per-current
retained sparse states with a grouped live-current workspace/window owner that
streams one producer L column across many current states and writes the
nonidentical output/trailing rows from that owner.

The first retained-value storage cut is now implemented as an opt-in substrate,
not a performance win. `KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_COMPACT_VALUES=1`
lets the retained BTF scalar-run executor allocate value storage only for
memberships belonging to current columns that have a selected retained executor
state; the structural sparse row descriptors remain complete. The first
best-member-only variant was rejected because producer-step advance also uses
member positions from the step index, causing many compact value lookups to
reject. The retained implementation now broadens compact storage to all
memberships for selected current columns, which preserves the producer-step
path.

Validation for the compact value substrate:

- `cmake --build build -j4` completed.
- `ctest --test-dir build --output-on-failure` passed both tests.
- `KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_EXEC=1
  KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_STEP_ADVANCE=1
  KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_COMPACT_VALUES=1
  build/kls_smoke` passed.
- On `ASIC_100ks`, the same-binary full-value baseline
  `build/kls_btf_state_exec_full_values_asic100ks_t4_r1_ref1.json` was
  residual-clean with `refactor_seconds_avg=18.8164198`; compact values in
  `build/kls_btf_state_exec_compact_values2_asic100ks_t4_r1_ref1.json` were
  also residual-clean with zero state rejects but measured `18.8269919`
  (`1.0006x` of full).
- On `ASIC_320ks`,
  `build/kls_btf_state_exec_full_values_asic320ks_t4_r1_ref1.json` measured
  `15.876602`, while
  `build/kls_btf_state_exec_compact_values_asic320ks_t4_r1_ref1.json` measured
  `15.8905815`, again residual-clean with zero state rejects (`1.0009x`).

This result is useful mainly because it proves the retained executor can be
decoupled from absolute per-member value offsets without losing correctness.
It does not close the CKTSO gap by itself: the selected-current compact store
still preserves per-membership sparse states and therefore leaves the core
producer-to-many-current owner missing. The next implementation should use this
indirection to replace those member-local value slices with an actual grouped
workspace/window, rather than expecting compact allocation alone to move time.

A current timeout-pair recheck at commit `e37f923` again isolates `pre2` as the
largest clean KLS-vs-CKTSO gap. With the same 4-thread, one-factor,
one-refactor, 120s-per-matrix cap over the large recon timeout pair, CKTSO
completed `pre2` in
`build/cktso_timeout_pair_e37f923_t4_r1_ref1_timeout120.jsonl` with
`analysis_seconds=3.586958`, `initial_factor_seconds=6.326116`,
`refactor_seconds_avg=4.875776`, `solve_seconds_avg=0.121033`, and
`spice_cycle_seconds=504.718198`; CKTSO timed out on `Hamrle3`. The matching
KLS run, `build/kls_timeout_pair_e37f923_t4_r1_ref1_timeout120.jsonl`, emitted
no rows, and its sidecar records 120s timeouts on both `pre2` and `Hamrle3`.
That makes `pre2` the actionable comparator because CKTSO solves it inside the
same cap.

The current `pre2` KLS analysis cost is not the timeout-sized gap. Analyze-only
KLS completed with auto/AMD-selected ordering in `10.5177002s` and explicit AMD
in `7.37985796s`
(`build/kls_pre2_analyze_auto_e37f923_t4.json` and
`build/kls_pre2_analyze_amd_e37f923_t4.json`). Both analyze rows show the same
dominant BTF block: `nblocks=29282`, `max_block=629628`, and estimated flops
`2.07660366e+11`. A 75s forced first-factor trace,
`build/kls_pre2_timeout_gap_trace_e37f923_t4_r1_ref0_timeout75.stderr`, reached
the dominant row-first pipeline (`block=13843`, `rows=629628`) but produced no
completed JSON before the cap.

A debugger-owned live sample of the same factor-only `pre2` window,
`build/pre2_timeout_gap_gdb_pty_e37f923.txt`, refines the hot point: the main
thread is waiting in `kls_factor` through
`kls_try_first_factor_row_uplooking_blocks_impl`, one row-first worker is
waiting in `kls_row_first_run_parallel_pipeline_phase` /
`kls_row_first_run_restartable_pipeline_suffix`, one active worker is in
`kls_row_first_supernode_panel_cache_build`, and the other row-pipeline workers
are blocked on the pipeline condition variable. This single sample is not a
full profile, but it is consistent with the previous scalar-output traces: the
gap is still cold numeric factorization inside the dominant BTF row/supernode
executor, with serialized panel-cache rebuild/suffix-pipeline work blocking
parallel progress. The paper-aligned next target remains a coarser CKTSO /
SubtreeLU-style row/supernode numeric owner that preserves and streams producer
panel work across the suffix instead of repeatedly rebuilding or replaying it
under a one-worker bottleneck.

The pivot path now removes that rebuild/reset bottleneck from the ordered
row-pipeline critical section. Dynamic-pivot epochs no longer rebuild prefix
row-supernode metadata or the prefix panel cache while peer workers wait.
Instead, KLS disables stale prefix supernode metadata and preserves only cached
panels that survive the pivot by exchanging tail columns or deactivating panels
whose dense block was affected. This applies to generic row-pipeline pivots,
separator-queue pipeline pivots, and active-rank suffix pivots; future panels
can still be appended as rows publish, so correctness falls back to scalar or
newly appended panels rather than relying on stale prefix metadata.

Validation for the no-rebuild pivot path:

- `cmake --build build -j4` completed.
- `ctest --test-dir build --output-on-failure` passed both tests after the
  smoke diagnostics were updated to require zero prefix rebuild counters on the
  pivot-epoch paths.
- `KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_EXEC=1
  KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_STEP_ADVANCE=1
  KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_STATE_COMPACT_VALUES=1
  build/kls_smoke` passed.
- A forced KLS-first `transient` check,
  `build/kls_transient_no_prefix_rebuild_t4_r1_ref0.json`, was residual-clean
  (`relative_residual_l2=4.12024204e-13`) with
  `initial_factor_seconds=2.03343866`, `factor_seconds_avg=1.15802924`, and
  zero row-pipeline/active-rank prefix panel rebuild counters.

On `pre2`, this is a real progress fix but not the CKTSO-gap closer. The old
75s trace at `e37f923`,
`build/kls_pre2_timeout_gap_trace_e37f923_t4_r1_ref0_timeout75.stderr`, only
printed the dominant-pipeline start before the timeout. The intermediate
no-active-rank-rebuild trace still did the same. After removing both prefix
rebuild and prefix reset from the ordered lock,
`build/kls_pre2_no_prefix_rebuild_trace_t4_r1_ref0_timeout75.stderr` reached
long-row commits around `completed=45171/629628`. A debugger-owned sample,
`build/pre2_no_prefix_rebuild_gdb_pty.txt`, now lands in
`kls_row_first_partial_apply_one_dep` rather than
`kls_row_first_supernode_panel_cache_build` or
`kls_row_first_supernodes_reset_prefix`; peer workers are still waiting on the
ordered pipeline lock. The non-traced 120s factor-only run
`build/kls_pre2_no_prefix_rebuild_factor_t4_r1_ref0_timeout120.json` still
timed out with no JSON output.

The remaining large gap is therefore sharper: after removing the full prefix
rebuild/reset lock bottleneck, the hard `pre2` rows still spend their time in
single-row scalar published-U replay. The 75s trace shows repeated long rows
with roughly `5M` to `10M` scalar U entries each, dominated by output/trailing
entries, while producer batching usually covers only a few thousand target
entries per row. The next paper-aligned step should move that output/trailing
update stream into a grouped producer/panel-to-many-current owner instead of
trying to recover more prefix cache state after pivots.

Auto ordering now starts or promotes to METIS for large weak-diagonal
moderate-degree matrices and for high-work dominant-BTF symbolic candidates
whose METIS trial builds a retained separator queue. This fills a direct
CKTSO-paper ordering gap without using a matrix-name exception: the direct
start requires at least 500k rows, 6-12 average structural entries per column,
bounded row/column degree, nonempty rows, and 20%-60% diagonal coverage; the
post-AMD/COLAMD promotion requires a 95%+ dominant BTF block, thousands of BTF
blocks, and very high estimated factor work. On `pre2`,
`build/kls_pre2_auto_direct_metis_skipretry_analyze_t4.json` selected
`ordering=metis`, retained `separator_pipeline_rows=4126`, and reduced auto
analysis to `12.8011618s` by skipping the unhelpful no-BTF retry for this
direct-start class. The matching 45s trace
`build/kls_pre2_auto_direct_metis_trace45.stderr` reached
`599203/629628` rows, compared with `36176/629628` for the AMD-auto trace
`build/kls_pre2_timeout_recheck_current_trace45.stderr`. This is substantial
progress, but the untraced 120s factor-only run
`build/kls_pre2_auto_direct_metis_t4_r1_ref0_timeout120.json` still timed out
with no JSON. The remaining `pre2` gap is still numeric pivot-tail execution:
late rows continue to replay multi-million-entry scalar U streams with only
small producer batches.

A same-commit timeout-pair recheck at `6633f97` keeps that conclusion and makes
`pre2` the largest clean KLS-vs-CKTSO timeout gap. With the same two-matrix
large-recon slice, four threads, `repeat=1`, `refactor-repeat=1`, and a 120s
per-matrix cap, CKTSO completed `pre2` in
`build/cktso_timeout_pair_6633f97_t4_r1_ref1_timeout120.jsonl`
(`analysis_seconds=3.618970`, `initial_factor_seconds=6.303858`,
`factor_seconds_avg=5.097129`, `refactor_seconds_avg=4.922565`,
`solve_seconds_avg=0.121867`) and timed out only on `Hamrle3`. KLS emitted no
rows in `build/kls_timeout_pair_6633f97_t4_r1_ref1_timeout120.jsonl`; the
sidecar records 120s timeouts on both `pre2` and `Hamrle3`. Current KLS
analyze-only on `pre2` still selects `ordering=metis` and finishes in
`11.6469298s` with `separator_pipeline_rows=4126`
(`build/kls_pre2_timeout_pair_6633f97_analyze_auto_t4.json`), so setup and
ordering are not the timeout-sized part.

The matching forced KLS-owned decomposition trace
`build/kls_pre2_timeout_pair_6633f97_forced_trace45.stderr` reached the
`589824/629628` checkpoint with `869321669` scalar U-entry touches, split into
`125359072` internal dependency touches and `743962597` trailing/output touches.
After that checkpoint, the 45s cap logged `2359` long rows through row `599195`;
those rows alone replayed `6410013403` scalar U entries, `4263724314` of them
trailing/output entries. Dynamic pivots accounted for `998` of those rows
(`42.3%`) and had zero producer-batch targets. Non-pivot rows did trigger the
retained producer path, but only for `3817` targets and `8001601` target U
entries, with `16813` panel-backed rows total. That is orders of magnitude
smaller than the scalar trailing/output stream. The largest gap is therefore
still the paper-level numeric owner for the pivoting dominant-BTF tail: a
producer/panel-to-many-current row/supernode executor that avoids replaying
published U rows through the scalar current-row path.

`scripts/summarize_row_pipeline_trace.py` now turns those row-pipeline stderr
files into repeatable CSV or JSON summaries, including the last progress event,
long-row count, pivoted long rows, scalar U entries, trailing/output entries,
producer target entries, and panel update rows. Applying it to same-commit
probes rejects the obvious policy-only alternatives. Eight-slot lookahead
(`build/kls_pre2_6633f97_lookahead8_trace45.stderr`) increased checkpoint
producer target U entries to `35781784` but advanced only to long row `598218`
in the 45s cap; 32-slot lookahead
(`build/kls_pre2_6633f97_lookahead32_trace45.stderr`) raised checkpoint target
entries to `94687787` but only reached long row `594654`. Explicit scale `2`
(`build/kls_pre2_6633f97_scale2_forced_trace45.stderr`) also reached fewer
tail rows than the auto-scale forced trace. Lower pivot tolerance reduces some
traced scalar work but not the timeout: `pivot_tol=1e-4`
(`build/kls_pre2_6633f97_tol1e4_forced_trace45.stderr`) still timed out in
`build/kls_pre2_6633f97_tol1e4_forced_t4_r1_ref0_timeout120.json`, and
`pivot_tol=1e-5` similarly timed out in
`build/kls_pre2_6633f97_tol1e5_forced_t4_r1_ref0_timeout120.json`. This keeps
the next code target on a real grouped-current numeric owner, not larger
lookahead, stronger scaling, or lower pivot tolerance around the existing
per-current scalar state.

A fresh current-tree timeout-pair recheck at `b6fb241` keeps the same largest
gap. With the same two-matrix large-recon slice, four threads, `repeat=1`,
`refactor-repeat=1`, and a 120s per-matrix cap, KLS emitted no rows in
`build/kls_timeout_pair_b6fb241_t4_r1_ref1_timeout120.jsonl`; the sidecar
records 120s timeouts on both `pre2` and `Hamrle3`. CKTSO completed `pre2` in
`build/cktso_timeout_pair_b6fb241_t4_r1_ref1_timeout120.jsonl`
(`analysis_seconds=4.515627`, `initial_factor_seconds=8.042882`,
`factor_seconds_avg=6.024751`, `refactor_seconds_avg=5.303302`,
`solve_seconds_avg=0.141433`) and timed out only on `Hamrle3`. The KLS failure
diagnostic for `pre2` completed analyze-only in `14.0594538s`, selected
`ordering=metis`, and again reported a 629,628-row dominant BTF block with
`separator_pipeline_rows=4126`, so ordering/analyze remains a secondary cost
rather than the timeout-sized gap.

The matching current-tree forced KLS-owned trace,
`build/kls_pre2_b6fb241_forced_trace45.stderr`, reached the same
`589824/629628` checkpoint. The trace summarized `869372808` scalar U-entry
touches, `744003478` of them trailing/output touches. Producer batching
covered only `912507` target U entries at that checkpoint, so scalar output was
`815.34x` larger than the producer-target U stream and producer targets covered
only `0.105%` of scalar U work. The 45s cap also logged `2132` long rows
through row `598963`; those long rows replayed `5129278578` scalar U entries,
`3574480212` of them trailing/output entries, with only `7403498`
producer-target U entries. Pivot rows were `874` of those long rows (`40.99%`),
and long-row scalar output was still `482.81x` larger than producer-target U
work. This recheck points to the same missing paper algorithm: a grouped
producer/panel-to-many-current numeric owner for the pivoting dominant-BTF
tail, not a wider ordering package, BLAS threshold, or larger lookahead around
the existing scalar current-row replay.

A July 1, 2026 current-`HEAD` rerun at `879d478` reproduces the timeout split
and keeps `pre2` as the clean KLS-only large timeout. With the same two-matrix
large-recon slice, four threads, `repeat=1`, `refactor-repeat=1`, and a 120s
per-matrix cap, KLS emitted no rows in
`build/kls_timeout_pair_head_t4_r1_ref1_timeout120.jsonl`; the sidecar records
120s timeouts on both `pre2` and `Hamrle3`. CKTSO completed `pre2` in
`build/cktso_timeout_pair_head_t4_r1_ref1_timeout120.jsonl`
(`analysis_seconds=4.225813`, `initial_factor_seconds=7.914385`,
`factor_seconds_avg=5.161273`, `refactor_seconds_avg=5.269733`,
`solve_seconds_avg=0.148217`) and timed out only on `Hamrle3`. The KLS
failure diagnostic for `pre2` completed analyze-only in `12.7715154s`, selected
`ordering=metis`, and retained the same 629,628-row dominant BTF block with
`separator_pipeline_rows=4126`, so analysis and ordering remain secondary to
numeric execution.

The matching 90s forced KLS-owned trace,
`build/kls_pre2_head_forced_trace90.stderr`, reached the `589824/629628`
checkpoint (`93.68%` complete). It recorded `869303361` scalar U-entry touches
and `743949878` trailing/output touches, while the retained producer-batch path
covered only `936240` target U entries. Overall scalar output was therefore
`794.61x` larger than producer-target U work, with producer targets covering
only `0.108%` of scalar U touches. The long-row tail logged `2083` long rows
through row `598912`; those rows replayed `4843807560` scalar U entries,
`3415853832` of them trailing/output entries, but only `7187407`
producer-target U entries. Pivot rows were `850` of those long rows
(`40.81%`), and long-row scalar output was still `475.26x` larger than
producer-target U work. This fresh rerun rejects the same policy-only fixes and
keeps the next implementation target on the paper-level grouped
producer/panel-to-many-current numeric owner for the pivoting dominant-BTF tail.

A post-pivot producer-lookahead diagnostic now tests one bounded slice of that
owner without enabling the broader experimental lookahead path. With
`KLS_ENABLE_ROW_PIPELINE_PIVOT_LOOKAHEAD=1`, a dynamic-pivot commit clears stale
lookahead, advances the order epoch, fills fresh lookahead states under the
post-pivot column order, and immediately streams the just-published pivot row
through the existing producer-batch kernel. The default diagnostic window uses
one slot per worker; `KLS_ROW_PIPELINE_PIVOT_LOOKAHEAD=<slots>` can widen it up
to the existing row-pipeline lookahead cap. Correctness passed
`./build/kls_smoke`,
`KLS_ENABLE_ROW_PIPELINE_PIVOT_LOOKAHEAD=1 ./build/kls_smoke`,
`KLS_ENABLE_ROW_PIPELINE_PIVOT_LOOKAHEAD=1 KLS_ROW_PIPELINE_PIVOT_LOOKAHEAD=8 ./build/kls_smoke`,
and `ctest --test-dir build --output-on-failure`.

The `pre2` trace response confirms that post-pivot producer ownership is the
right direction but not enough. At the same `589824/629628` checkpoint, the
four-slot diagnostic
`build/kls_pre2_pivot_lookahead_trace90.stderr` raised producer target U
entries to `2112190` and reduced the overall scalar-output/producer-target
ratio to `351.93x`; the eight-slot diagnostic
`build/kls_pre2_pivot_lookahead8_trace90.stderr` raised target U entries to
`3601125` and reduced the ratio to `206.30x`; the sixteen-slot diagnostic
`build/kls_pre2_pivot_lookahead16_trace90.stderr` raised target U entries to
`6252151` and reduced the ratio to `118.72x`. Eight and sixteen slots also
advanced the traced long-tail maximum to about row `608k`, versus `598912` for
the default trace. However, non-traced 125s probes with four and eight slots
still timed out
(`build/kls_pre2_pivot_lookahead_t4_r1_ref0_timeout125.json` and
`build/kls_pre2_pivot_lookahead8_t4_r1_ref0_timeout125.json`). This keeps the
post-pivot path opt-in: it proves that pivot rows were a real missing producer
surface, but closing the CKTSO gap still requires a persistent grouped-current
numeric owner that avoids materializing and replaying independent sparse
current states.

The follow-up pivot-row split in `scripts/summarize_row_pipeline_trace.py`
confirms why the default-capped diagnostic did not close the timeout. With
`KLS_ROW_PIPELINE_PIVOT_LOOKAHEAD=8`, only `14` of `1826` logged pivot long rows
actually applied a producer batch, and pivot-row producer target U entries were
just `77502` against `3.870788451e9` pivot scalar output entries. Disabling the
lookahead raw-input cap for the same pivot-only diagnostic
(`KLS_ROW_PIPELINE_LOOKAHEAD_MAX_INPUTS=0`,
`build/kls_pre2_pivot_lookahead8_maxinputs0_trace90.stderr`) made `916` of
`921` logged pivot long rows batch the pivot producer, raised pivot producer
target U entries to `15566778`, and reduced the pivot scalar-output/producer
target ratio from `49944x` to `109.82x`. The non-traced 125s run still timed
out in `build/kls_pre2_pivot_lookahead8_maxinputs0_t4_r1_ref0_timeout125.json`,
so broader pivot eligibility is not sufficient by itself.

A sparse-state preservation prototype was rejected before commit. It tried to
carry unclaimed lookahead states across a pivot by applying the column exchange
to each saved sparse workspace and rebuilding its local dependency heap, rather
than clearing and refilling from raw input. Correctness passed the normal smoke
test and the opt-in preserve smoke test, but the focused `pre2` traces were
negative. With the normal input cap,
`build/kls_pre2_pivot_lookahead8_preserve_trace90.stderr` increased scalar
touches to `982777460` and only batched `281` of `1095` logged pivot long rows,
versus `868249990` scalar touches for the capped refill path. With the cap
disabled,
`build/kls_pre2_pivot_lookahead8_maxinputs0_preserve_trace90.stderr` increased
scalar touches to `1139176876` and worsened the pivot scalar-output/producer
target ratio to `147.59x`, versus `859854990` scalar touches and `109.82x` for
the no-cap refill path. This rejects preserving independent sparse lookahead
states as the persistent owner. The next implementation target remains a
coarser grouped-current numeric owner that advances a batch as one object,
rather than carrying many separately caught-up sparse row workspaces across
pivots.

A grouped lookahead catch-up prototype was also tested and rejected before
commit. It tried to batch shared ready roots across unclaimed lookahead states
before the normal current-producer scan, directly testing whether the missing
owner could be approximated by synchronizing existing sparse lookahead rows
around earlier completed producers. The trace moved in the right direction but
not nearly enough: against the no-cap pivot refill control
`build/kls_pre2_pivot_lookahead8_maxinputs0_trace90.stderr`, the grouped
catch-up trace
`build/kls_pre2_group_catchup_pivot8_maxinputs0_trace90.stderr` reduced overall
scalar U touches only from `859854990` to `859296862` and improved the overall
scalar-output/producer-target ratio only from `28.44x` to `26.58x`. Pivot
long-row coverage rose from `916/921` to `1108/1113`, and the pivot
scalar-output/producer-target ratio improved from `109.82x` to `83.68x`, but
both traces still stopped at the same `589824/629628` progress checkpoint. The
non-traced 125s factor-only probe
`build/kls_pre2_group_catchup_pivot8_maxinputs0_t4_r1_ref0_timeout125.json`
also timed out with no JSON row. This rejects grouped catch-up of separately
owned sparse lookahead states as the large CKTSO-gap closer. The required
algorithm is still a production grouped-current or row/supernode numeric owner
that owns the current batch and its target map, rather than repeatedly catching
up independent sparse current states.

A CKTSO static-pivoting shortcut was tested and rejected before commit. CKTSO's
preprocessing uses maximum-weight matching to reduce dynamic pivots, so the
trial allowed the existing large SPRAL pre-static candidate to accept its
matched symbolic pattern without first requiring a successful trial numeric
factorization. This specifically tested whether `pre2` is timing out because
KLS never reaches the post-factor SPRAL acceptance gate. Correctness smoke
passed for the normal binary and for the opt-in
`KLS_ENABLE_PRE_STATIC_SPRAL_SYMBOLIC_ONLY=1` path, but the focused trace was
negative. The 45s `pre2` trace
`build/kls_pre2_prestat_spral_symbolic_trace45.stderr` reached only
`393216/629628` rows, while the retained METIS-start trace reaches the late
tail in the same window. At 90s,
`build/kls_pre2_prestat_spral_symbolic_trace90.stderr` reached the same
`589824/629628` checkpoint as `build/kls_pre2_head_forced_trace90.stderr`, but
long-row scalar U work worsened from `4.843807560e9` to `8.057231222e9` and
the long-row scalar-output/producer-target ratio worsened from `475.26x` to
`578.88x`. This rejects accepting a large SPRAL static match without numeric
proof as the clear `pre2` gap closer. Static pivoting remains useful through
the existing guarded numeric-trial paths, but the timeout-sized gap is still
the row-pipeline numeric owner for the pivoting tail.

A current timeout-pair recheck at `54d839d` again isolates `pre2` as the
largest clean KLS-only timeout. With four threads, `repeat=1`,
`refactor-repeat=1`, and a 120s per-matrix cap over the `pre2`/`Hamrle3`
slice, CKTSO completed `pre2` in
`build/cktso_timeout_pair_54d839d_t4_r1_ref1_timeout120.jsonl`
(`analysis_seconds=4.252091`, `initial_factor_seconds=7.752507`,
`factor_seconds_avg=5.061632`, `refactor_seconds_avg=5.312715`,
`solve_seconds_avg=0.160632`) and timed out on `Hamrle3`. The matching KLS run
`build/kls_timeout_pair_54d839d_t4_r1_ref1_timeout120.jsonl` emitted no rows;
its sidecar records 120s timeouts on both matrices. The `pre2` KLS
analyze-only diagnostic completed in `12.5914825s`, selected `ordering=metis`,
and reported the same `629628`-row dominant BTF block, so ordering/setup is not
the largest part of this gap. `Hamrle3` remains a shared hard case because
CKTSO also times out there; `pre2` is the actionable comparator because CKTSO
finishes the same case inside the cap.

The matching current `pre2` factor-only trace
`build/kls_pre2_timeout_pair_54d839d_trace_t4_r1_ref0_timeout75.stderr`
reached the familiar `589824/629628` progress checkpoint before the 75s cap.
At that checkpoint KLS had replayed `869303159` scalar U entries, including
`743949331` trailing/output entries, while the retained producer-batch path
covered only `957366` target U entries; scalar output was therefore `777.08x`
larger than producer-target U work. The long-row tail logged `2687` rows
through row `599516`, replaying `8421328975` scalar U entries and
`5208144838` trailing/output entries, with only `8876221` producer-target U
entries. Dynamic-pivot long rows were `1141` of those rows and still had zero
producer rows in the baseline trace. This pins the largest clean CKTSO gap on
numeric first-factor row-pipeline replay in the pivoting dominant-BTF tail, not
refactorization, solve time, ordering, or the timeout limit.

A wider pivot-lookahead no-cap probe was also rejected as a substitute for the
live grouped-current owner. With `KLS_ROW_PIPELINE_PIVOT_LOOKAHEAD=32` and
`KLS_ROW_PIPELINE_LOOKAHEAD_MAX_INPUTS=0`,
`build/kls_pre2_pivot_lookahead32_maxinputs0_trace90.stderr` greatly increased
producer coverage at the `589824/629628` checkpoint: producer target U entries
rose to `114915909` and the scalar-output/producer-target ratio fell to
`6.28x`, versus `28.44x` for the eight-slot no-cap control and `794.61x` for
the default trace. Long-row scalar replay also dropped to `1776215140` entries
with a `21.60x` long-row output/target ratio. However, the trace still stopped
at the same checkpoint and the matching non-traced 125s factor-only probe
`build/kls_pre2_pivot_lookahead32_maxinputs0_t4_r1_ref0_timeout125.json`
timed out with no JSON row. This rejects keeping a much larger pool of
independent sparse lookahead states as the CKTSO-gap closer; KLS needs a
production grouped-current numeric owner that owns the batch and target map,
not more advisory sparse rows around the existing replay path.

A July 1, 2026 timeout-pair refresh on the same `pre2`/`Hamrle3` large slice
keeps `pre2` as the largest clean KLS-only timeout gap. With four threads,
`repeat=1`, `refactor-repeat=1`, and a 120s per-matrix cap, CKTSO completed
`pre2` in
`build/cktso_timeout_pair_readyroot_recheck_t4_r1_ref1_timeout120.jsonl`
(`analysis_seconds=4.544027`, `initial_factor_seconds=7.761729`,
`factor_seconds_avg=5.175383`, `refactor_seconds_avg=3.904002`,
`solve_seconds_avg=0.142911`, `spice_cycle_seconds=413.093054`) and timed out
only on `Hamrle3`. The matching KLS run
`build/kls_timeout_pair_readyroot_recheck_t4_r1_ref1_timeout120.jsonl` emitted
no rows; its sidecar records 120s timeouts on both matrices. The `pre2`
analyze-only diagnostic in that sidecar completed in `13.1605754s`, selected
`ordering=metis`, and again reported the `629628`-row dominant BTF block. This
keeps the largest actionable gap inside numeric first factorization of the
pivoting dominant block; `Hamrle3` remains shared-hard because CKTSO also times
out there.

A ready-root batching prototype was rejected before commit. It kept the
commit-drain lock and, whenever a current row popped a ready scalar dependency,
tried to batch that producer across active workers and unclaimed lookahead
states with the same root. This directly tested whether the missing
CKTSO/SubtreeLU grouped-current owner could be approximated by synchronizing
existing sparse states at ready roots. The same-checkpoint `pre2` result was
negative: `build/kls_pre2_baseline_trace45.stderr` and
`build/kls_pre2_readyroot_trace45.stderr` both reached `458752/629628` rows
under a 45s cap, but the prototype increased scalar U work from `682232871` to
`682418775`, reduced producer target U entries from `854772` to `591684`, and
worsened the scalar-output/producer-target ratio from `682.25x` to `985.93x`.
It also inflated underfilled probes from `264766` to `2427994`. This rejects
ready-root catch-up over independently owned sparse states as the missing
CKTSO-scale mechanism. The next paper-aligned implementation still needs a
real live grouped-current numeric owner that advances one batch and target map
coarsely instead of repeatedly probing and catching up separate row states.

An unlocked wait-partial prototype was also rejected before commit. It targeted
one plausible CKTSO/SubtreeLU scheduling gap: KLS advances waiting current rows
under the shared row-pipeline mutex, while the papers describe current rows
using finished predecessors as they become available. The prototype let
active-rank wait drains run outside the mutex with a snapshot of published U
storage, hid the mutable sparse state from producer-batch probing, and made
commits wait before any U reallocation or pivot column exchange that could
invalidate unlocked readers. Correctness checks passed
(`./build/kls_smoke`,
`KLS_ENABLE_ROW_PIPELINE_UNLOCK_WAIT_PARTIAL=1 ./build/kls_smoke`,
`ctest --test-dir build --output-on-failure`, plus residual-clean forced
`add20` and `bcircuit` probes), but the `pre2` timing evidence was not strong
enough to keep the code.

The 45s trace was only mildly positive:
`build/kls_pre2_unlock_wait_trace45.stderr` reached the same
`458752/629628` checkpoint as `build/kls_pre2_baseline_trace45.stderr`,
reduced scalar U touches from `682232871` to `682160113`, and raised producer
target U entries from `854772` to `931806` (`682.25x` to `625.77x` scalar
output per producer-target U entry). The longer 75s trace rejected it as a
CKTSO-gap closer:
`build/kls_pre2_unlock_wait_trace75.stderr` reached the same
`589824/629628` checkpoint as the saved baseline
`build/kls_pre2_timeout_pair_54d839d_trace_t4_r1_ref0_timeout75.stderr`, but
total scalar U touches were essentially unchanged (`869321969` versus
`869303159`), producer target U entries were slightly lower (`936438` versus
`957366`), and the long-row log advanced only to row `595435` versus `599516`
for the baseline. The no-trace factor-only probe
`build/kls_pre2_unlock_wait_factor_t4_r1_ref0_timeout125.json` still timed out
with no JSON row. This rejects mutex release around existing independent sparse
current states as the clear missing mechanism; the larger paper-level gap
remains a live grouped-current numeric owner that changes the state ownership
and producer-stream reuse, not just where the current scalar drain holds the
mutex.

The first-factor row-pipeline trace now has the matching producer-batch
state-union diagnostic. For each successful producer batch, KLS counts the
sum of target sparse-state pattern rows and the exact unique column rows across
those targets, reporting them as `producer_state_rows` and
`producer_unique_state_rows` in both progress and long-row trace lines. The
summary helper also reports the state-rows-per-unique ratio. This is trace-only
and does not affect default execution.

The first `pre2` trace with the new fields,
`build/kls_pre2_state_union_trace45.stderr`, confirms that the active producer
batches have real grouped-workspace reuse but the live target window is still
too narrow. At the `589824/629628` checkpoint, successful batches covered
`1160592` target sparse-state rows but only `415135` unique rows, a `2.80x`
collapse opportunity. In the logged long rows, the matching split was
`2121701` state rows versus `827050` unique rows (`2.57x`). However, producer
target U entries were still only `921402` against `743978570` scalar output
entries, about `807x` smaller than the scalar output stream. This sharpens the
next implementation target: grouped state storage is useful once targets are
present, but the CKTSO-scale gap still requires a larger live grouped-current
owner that brings far more current rows into producer-owned batches.

A fresh timeout-case recheck at `88aaf4b` confirms the largest gap component
without relying on older artifacts. The same two-row large slice, four threads,
`repeat=1`, `refactor-repeat=1`, and 120s per-matrix cap produced
`build/cktso_timeout_pair_88aaf4b_t4_r1_ref1_timeout120.jsonl` for CKTSO and
`build/kls_timeout_pair_88aaf4b_t4_r1_ref1_timeout120.jsonl` for KLS. CKTSO
completed `pre2` (`analysis_seconds=4.162433`,
`initial_factor_seconds=5.589610`, `factor_seconds_avg=3.541263`,
`refactor_seconds_avg=4.520806`, `solve_seconds_avg=0.129361`) and timed out
only on `Hamrle3`. KLS emitted no rows; its failure sidecar records 120s
timeouts on both matrices. The KLS `pre2` analyze-only diagnostic completed in
`13.0557551s`, selected METIS, and reported the same `629628`-row dominant BTF
block, so analyze/order setup remains bounded far below the timeout.

The matching capped `pre2` factor-only trace,
`build/kls_pre2_timeout_pair_88aaf4b_trace75.stderr`, reached the same
`589824/629628` checkpoint before the 75s cap. At that point KLS had replayed
`869301906` scalar U entries, including `743948197` trailing/output entries,
while successful producer batches covered only `944607` target U entries. The
scalar-output/producer-target ratio is therefore `787.57x`. The new state-union
diagnostic also reports `1190540` producer target state rows but only `425873`
unique rows (`2.80x` reuse), which means grouped state storage is not the
missing scale by itself; KLS is failing to make enough live currents targets of
producer-owned numeric updates. This keeps the paper-level target on a
persistent grouped-current row/supernode numeric owner for first-factor
pivoting tails, rather than another timeout-policy, ordering-package, BLAS, or
small lock-refactor change.

An active-current catch-up batching probe was added behind
`KLS_ENABLE_ROW_PIPELINE_ACTIVE_CATCHUP_BATCH=1` and rejected as a default
gap-closer. The idea directly tested whether KLS was missing CKTSO/SubtreeLU
producer reuse because active current rows were not being advanced to the
just-completed producer before producer-batch selection. The guarded path
applies ready dependencies below the producer for active rows that already
contain that producer, then reruns the normal producer-batch eligibility check.
The implementation keeps extra trace fields:
`producer_active_catchup_attempts`, `producer_active_catchup_deps`, and
`producer_active_catchup_targets`; the summary helper includes the same fields
for progress and long-row traces.

Correctness checks passed (`./build/kls_smoke`,
`KLS_ENABLE_ROW_PIPELINE_ACTIVE_CATCHUP_BATCH=1 ./build/kls_smoke`,
`ctest --test-dir build --output-on-failure`, plus forced KLS-first `add20`
and `bcircuit` residual probes). The focused `pre2` trace was mildly positive
but far below CKTSO scale. Against
`build/kls_pre2_timeout_pair_88aaf4b_trace75.stderr`,
`build/kls_pre2_active_catchup_trace75.stderr` reached the same
`589824/629628` checkpoint, eliminated `producer_reject_not_root` by creating
`62131` catch-up targets and advancing `79883` scalar dependencies, but raised
producer-target U entries only from `944607` to `1054425`. The scalar
output/producer-target ratio improved only from `787.57x` to `705.53x`. The
matching no-trace probe
`build/kls_pre2_active_catchup_factor_t4_r1_ref0_timeout125.json` still timed
out with no JSON row. This rejects active-state catch-up around the existing
independent sparse row states as the clear missing mechanism; the required
paper-level work remains a persistent grouped-current numeric owner that makes
many more live currents targets of producer-owned row/supernode updates.

The next trace-only slice measures the ready-root candidate surface before the
current producer-batch filters discard it. KLS now reports
`producer_candidate_targets`, `producer_candidate_target_u_entries`,
`producer_underfilled_target_u_entries`, and
`producer_low_saved_stream_target_u_entries` in row-pipeline progress and
long-row trace lines, and the summary helper reports candidate/scalar ratios
and accepted-candidate share. This answers whether the current gap could be
closed by lowering `KLS_ROW_FIRST_PRODUCER_BATCH_MIN_SAVED_STREAM` or admitting
single-target batches.

On `pre2`, the answer is no. The default candidate-surface trace
`build/kls_pre2_candidate_surface_trace75.stderr` reached the same
`589824/629628` checkpoint with `28216007` candidate target U entries before
thresholds, of which only `954507` were accepted, `5489691` were underfilled,
and `21771809` failed the saved-stream gate. That looks like a large discard
ratio, but it is still tiny next to the `743979273` scalar U-output entries:
even accepting every candidate would leave scalar output `26.37x` larger than
candidate producer-target work. With active catch-up enabled,
`build/kls_pre2_candidate_surface_active_catchup_trace75.stderr` eliminated
`not_root` and raised candidate target U entries only to `32567534`; scalar
output was still `22.84x` larger. This rejects threshold tuning and
active-state catch-up as CKTSO-scale fixes. The missing first-factor mechanism
has to create a much larger live target window or persistent grouped-current
owner before producer filtering, not merely accept more of the current
ready-root surface.

The timeout-case refocus now drops `Hamrle3` from the tight tuning loop because
CKTSO also times out there. The interrupted two-row reference run,
`build/cktso_timeout_pair_refocus_t4_r1_ref1_timeout120.jsonl`, completed
`pre2` under the same 120s cap (`analysis_seconds=3.670268`,
`initial_factor_seconds=6.197750`, `refactor_seconds_avg=4.851256`,
`solve_seconds_avg=0.118040`) and timed out only on `Hamrle3`. The matching
KLS `pre2`-only run,
`build/kls_pre2_refocus_t4_r1_ref1_timeout120.jsonl`, emitted no rows; its
failure sidecar records analyze-only success in `11.6621665s`, METIS ordering,
and the same `629628`-row dominant BTF block before the 120s timeout. Thus KLS
spends at least `108.34s` after analysis inside first numeric factorization
without completing the benchmark row, already `17.48x` CKTSO's `pre2` initial
factor time.

KLS also now has a trace-only compact live-window diagnostic for this
first-factor question. `KLS_TRACE_ROW_PIPELINE_COMPACT_WINDOW=1` keeps compact
symbolic sparse states for a bounded current/future row window until their rows
complete, advances those states structurally as producer U rows are published,
and reports `compact_window_*` counters in row-pipeline progress and long-row
trace lines. `KLS_ROW_PIPELINE_COMPACT_WINDOW=<slots>` controls the window
width, and `KLS_ROW_PIPELINE_COMPACT_WINDOW_MAX_ENTRIES=<entries>` caps one
state. The summary helper reports compact-window target/scalar ratios.

The focused `pre2` evidence rejects bounded independent compact states as the
CKTSO-sized missing mechanism. The baseline 45s trace
`build/kls_pre2_refocus_trace45.stderr` reached `589824/629628` rows with
`743927452` scalar U-output entries and only `1006368` producer-target U
entries (`739.22x`). The completed-frontier 64-state compact-window trace,
`build/kls_pre2_compact_window64_completed_trace45.stderr`, raised the modeled
target surface to `9081488` U entries with zero compact overflows, but scalar
output remained `81.92x` larger. This is more surface than the current accepted
producer batches, but still far below CKTSO scale. The paper-aligned next step
therefore remains a coarser grouped-current row/supernode numeric owner, not a
bounded side window of independent sparse symbolic states.

The direct sparse numeric version of that compact window was implemented behind
`KLS_ENABLE_ROW_PIPELINE_COMPACT_EXEC=1` and is also rejected as a default
CKTSO-gap closer. The prototype keeps numeric values and retained L entries in
the compact states, advances unreserved states when producer U rows publish,
and lets a worker claim the prepared sparse state instead of replaying the row
prefix from scratch. Correctness smoke passed, and forced KLS-first `add20` and
`bcircuit` probes stayed residual-clean:
`build/kls_add20_compact_exec_lean_t4_r1_ref0.json` reported
`relative_residual_l2=3.39e-16`, and
`build/kls_bcircuit_compact_exec_lean_t4_r1_ref0.json` reported
`relative_residual_l2=8.17e-17`.

The focused `pre2` timing rejects the shape even after trimming the compact
state storage. Retained compact-state L entries now live in lean `(dep,value)`
arrays, and compact-claim tracing baselines against the worker's current
allocation counters so the earlier hundreds of millions of reported local
reserve-copy entries are no longer misread as real per-row churn. The corrected
16-state compact-exec trace
`build/kls_pre2_compact_exec_lean_w16_tracefix_trace45.stderr` reached the
same `589824/629628` checkpoint as
`build/kls_pre2_compact_exec_control_trace45.stderr` and raised compact target
work to `69009179` U entries, but scalar U output remained `10.83x` larger.
The untraced 16-state factor probe
`build/kls_pre2_compact_exec_lean_w16_factor_t4_r1_ref0_timeout125.json` still
timed out with no JSON row. This means the missing CKTSO/SubtreeLU mechanism is
not "more independent sparse states"; it has to be a coarser grouped-current
owner that streams producers across a batch instead of replaying most
trailing/output work through separately materialized sparse states.

That compact-exec prototype now has the first grouped producer/current update:
when compact numeric execution is enabled, a completed producer row collects
compact states whose root dependency is that producer, pops/appends their L
entries, and streams the published U row once across the compact target list.
The trace records `compact_window_batches` and
`compact_window_stream_u_entries`, and the summary helper reports compact
target/stream reuse. Correctness stayed clean under
`KLS_ENABLE_ROW_PIPELINE_COMPACT_EXEC=1`: default and compact-exec smoke passed,
`ctest --test-dir build --output-on-failure` passed, and forced KLS-first
`add20`/`bcircuit` probes with a 64-state compact window reported
`relative_residual_l2=3.39e-16` and `8.17e-17`.

The grouped compact update improves the rejected compact-exec shape but still
does not close the `pre2` gap. The old 16-state compact-exec trace,
`build/kls_pre2_compact_exec_lean_w16_tracefix_trace45.stderr`, reached
`589824/629628` rows with `872880120` scalar U entries and `69009179` compact
target U entries. The grouped 16-state trace,
`build/kls_pre2_compact_group_w16_trace45.stderr`, reached the same checkpoint
with `651254322` scalar U entries and streamed `13169290` producer U entries
across `69926314` compact target U entries (`5.31x` reuse). A wider 64-state
grouped trace, `build/kls_pre2_compact_group_w64_trace45.stderr`, raised
compact target work to `153055010` U entries with `9.95x` target/stream reuse,
but still reached only the same checkpoint and remained slightly behind the
same-binary no-compact control in scalar output. The no-trace 64-state factor
probe, `build/kls_pre2_compact_group_w64_factor_t4_r1_ref0_timeout125.json`,
timed out with no JSON row. This keeps the grouped compact path opt-in: it is
a better paper-aligned substrate than independent compact states, but the
bounded compact window is still not the missing CKTSO-scale grouped
row/supernode owner.

The compact-exec target lookup is now indexed by each live compact state's root
dependency instead of scanning every compact slot at every producer publication.
On the 64-state `pre2` trace, this preserves the same grouped compact surface
while cutting `compact_window_probes` from `36497216` to `2078045`, matching the
number of actual compact targets. The indexed 45s trace
`build/kls_pre2_compact_index_w64_trace45.stderr` reached the same
`589824/629628` checkpoint with `153075875` compact target U entries and
`15382089` streamed producer U entries (`9.95x` reuse). This removes avoidable
scan overhead in the opt-in compact prototype, but it is not the large missing
paper mechanism: the no-trace indexed run
`build/kls_pre2_compact_index_w64_notrace120.json` still timed out with no JSON
row under the 120s cap.

The compact sparse state now also uses append-only pattern storage with an
open-addressed row-to-position index, instead of keeping each state sorted by
binary-search insertion and shifting values on every new fill. This is still
only an opt-in compact-window substrate, but it attacks the cost exposed by
wider grouped-current probes. The 64-state `pre2` trace stayed effectively
neutral (`build/kls_pre2_compact_hashindex_w64_trace45.stderr` reached the same
`589824/629628` checkpoint and preserved about `153M` compact target U entries).
The wider 256-state run improved from
`build/kls_pre2_compact_index_w256_trace45.stderr`, which reached only
`393216/629628` rows, to
`build/kls_pre2_compact_hashindex2_w256_trace45.stderr`, which reached
`458752/629628` rows with `245453325` compact target U entries and
`12935884` streamed producer U entries (`18.97x` reuse). This confirms that
per-state sparse insertion was a real overhead in wider windows. It still does
not close `pre2`: the 256-state path remains behind the 64-state checkpoint, so
the missing paper-scale mechanism is still a coarser grouped row/supernode
owner, not simply a larger independent sparse window.

A post-pivot compact refill prototype was rejected before commit. It tested
whether compact execution was missing an obvious pivot-row publication step by
clearing stale compact states after a dynamic column pivot, advancing the order
epoch, then filling fresh compact states and streaming the just-published pivot
row through the compact producer batch. The focused 64-state `pre2` trace
`build/kls_pre2_pivot_compact_w64_trace45.stderr` reached the same
`589824/629628` checkpoint as
`build/kls_pre2_compact_hashindex_w64_trace45.stderr`, but compact target U work
stayed flat (`153016248` versus `153032141`), scalar U output rose slightly
(`556025156` versus `555564933`), and normal producer target U entries fell
(`903186` versus `930291`). This rejects "post-pivot compact rehydrate" as the
clear missing mechanism. The useful conclusion is narrower: pivoted rows are a
real missing surface for producer ownership, but bounded compact states still do
not create the paper-level grouped current owner.

An early-abort guard for separator-private dynamic pivots was tested and
rejected before commit. The experiment stopped the private separator phase on
`pre2` after `1024` dynamic column pivots, before the later `pivot-reject`,
then fell back to the guarded full row pipeline without publishing partial
private rows. The same 35s cap gives no evidence of a speedup:
`build/kls_pre2_nocap_current_trace35.stderr` and
`build/kls_pre2_cap1024_current_trace35.stderr` both timed out at the
`589824/629628` checkpoint, while the capped run raised scalar U output from
`534466380` to `728292143`. This rejects "fail separator-private earlier" as
the large missing paper mechanism; the gap remains the first-factor
producer/current grouping problem in the pivoting dominant BTF tail.

A direct ready-supernode producer handoff was also rejected before commit. The
prototype let producer-batch probing consume a ready supernode or cached-panel
run for an active target row when the scalar producer candidate was rejected
only because the target's root dependency was inside that ready run. Default
and opt-in smoke tests passed, and `ctest --test-dir build --output-on-failure`
passed. However, the focused `pre2` trace
`build/kls_pre2_supernode_producer_trace45.stderr` timed out at the same
`589824/629628` progress checkpoint as the current no-compact trace, and the
new path did not trigger at that checkpoint (`producer_supernode_targets=0`,
`producer_supernode_rows=0`). Scalar U output was essentially unchanged versus
`build/kls_pre2_nocap_current_trace35.stderr` (`534573340` versus
`534466380`). This rejects the simple "let producer probing advance one
ready-supernode target" idea as the missing paper mechanism. The useful target
remains a coarser grouped-current owner that makes several live current states
share producer-row or producer-supernode numeric work, not a one-target
catch-up around an already ready run.

A post-pivot exact-supernode-retention prototype was rejected before commit.
The idea followed SubtreeLU's row-by-row supernode detection more directly:
after a dynamic pivot, KLS kept the supernode metadata arrays alive but raised a
validity floor to the pivoted row, so old prefix metadata was ignored and only
newly published post-pivot rows could form exact adjacent U-supernodes. The
prototype built cleanly and passed `./build/kls_smoke` plus
`ctest --test-dir build --output-on-failure`, but it did not activate the
missing hot path on `pre2`. Against
`build/kls_pre2_current_refresh_trace45.stderr`, the prototype trace
`build/kls_pre2_supernode_floor_trace45.stderr` still timed out at the same
`589824/629628` periodic checkpoint; tail inspection reached slightly fewer
rows (`599260` max row and `599257` max completed versus `599283` and `599280`
for the control), scalar U output was effectively unchanged
(`534467391` versus `534467374`), and the long-row tail still reported
`scalar_runs=0`. This rejects "resume exact adjacent supernodes after pivots"
as the clear CKTSO/SubtreeLU gap closer. The remaining paper-level gap is
broader than metadata lifetime: the pivoting tail needs a grouped
row/supernode/current owner that can share producer work across many live
current rows, including rows that do not satisfy KLS's exact adjacent
supernode-shape predicate.

A compact-window pivot-preservation prototype was rejected before commit. This
tested a stronger version of the earlier post-pivot compact refill idea: instead
of clearing prepared compact states after a dynamic column pivot, the prototype
permuted each state accumulator through the column swap, rebuilt its sparse
index and root-dependency heap, conservatively discarded states that had already
consumed either swapped column, and streamed the just-published pivot row through
the preserved compact buckets. Default smoke, opt-in compact-preserve smoke,
`ctest --test-dir build --output-on-failure`, and residual probes stayed clean:
`build/kls_add20_compact_pivot_preserve_t4_r1_ref0.json` reported
`relative_residual_l2=3.38738965e-16`, and
`build/kls_bcircuit_compact_pivot_preserve_t4_r1_ref0.json` reported
`8.16829708e-17`.

The focused `pre2` evidence was not good enough to retain the code. Against the
same-binary compact-exec 64-state control
`build/kls_pre2_compact_w64_current_trace45.stderr`, the preserve trace
`build/kls_pre2_compact_pivot_preserve_w64_trace45.stderr` reached the same
`589824/629628` periodic checkpoint. It raised compact target U work from
`153034087` to `170767920` and reduced compact refills, so the preservation
path did activate. However, scalar U output also rose slightly
(`555604758` to `556795398`), and the no-trace factor probe
`build/kls_pre2_compact_pivot_preserve_w64_factor_t4_r1_ref0_timeout125.json`
still timed out with no JSON row. This rejects "preserve bounded compact
states across dynamic pivots" as the missing CKTSO/SubtreeLU mechanism. The
useful conclusion is that retaining sparse side states across pivots can expose
more producer/target surface, but the surface is still far too bounded and
state-local; the next implementation needs a coarser owner for many current
rows, not better lifetime management for independent compact states.

The BTF scalar-run group descriptor now has an opt-in live-state union
diagnostic,
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_LIVE_STATE_STATS=1`, that measures
the sparse workspace shape needed by the missing grouped-current owner before
the retained-state executor materializes per-current sparse states. For each
multi-current BTF scalar-run group, KLS reuses the existing structural row
collector to count both the duplicated per-member retained rows and the exact
union of local state rows that a single grouped live workspace would need. The
fields are visible through `kls_stats`, benchmark JSON/text, and the gap
decomposition CSV as
`refactor_btf_scalar_run_group_live_state_groups`,
`refactor_btf_scalar_run_group_live_state_currents`,
`refactor_btf_scalar_run_group_live_state_rows`,
`refactor_btf_scalar_run_group_live_state_unique_rows`,
`refactor_btf_scalar_run_group_live_state_reused_rows`,
`refactor_btf_scalar_run_group_live_state_max_currents`,
`refactor_btf_scalar_run_group_live_state_max_rows`, and
`refactor_btf_scalar_run_group_live_state_max_unique_rows`.

Validation stayed clean:
`cmake --build build --target kls_bench kls_smoke -j2`,
`./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_LIVE_STATE_STATS=1 ./build/kls_smoke`,
and `ctest --test-dir build --output-on-failure` all passed. The default
`ASIC_100ks` control
`build/kls_asic100ks_live_state_default_t4_r1_ref1.json` left the descriptor
and live-state counters at zero and measured `initial_factor_seconds=0.752872`
and `refactor_seconds_avg=0.0471574`. With the diagnostic enabled,
`build/kls_asic100ks_live_state_stats_t4_r1_ref1.json` stayed residual-clean
(`relative_residual_l2=1.92251861e-15`) and reported `13,116` live-state
groups covering `104,280` current memberships, `58,663,555` duplicated
per-current state rows, `20,355,030` unique grouped rows, `38,308,525` reusable
rows, max fanout `362`, max duplicated group rows `489,578`, and max unique
group rows `27,386`. `ASIC_320ks` in
`build/kls_asic320ks_live_state_stats_t4_r1_ref1.json` was also residual-clean
(`2.08436857e-15`) and reported `33,856` groups, `121,701` memberships,
`50,818,849` duplicated rows, `18,462,708` unique grouped rows, `32,356,141`
reusable rows, max fanout `444`, max duplicated group rows `501,702`, and max
unique group rows `21,683`.

This fills a diagnostic gap between the paper-level idea and the current
retained-state prototype. The whole-group state collapse is real, about `2.88x`
on `ASIC_100ks` and `2.75x` on `ASIC_320ks`, but it is smaller than the
producer-step row-union collapse observed after materialization. That means the
next implementation should not simply allocate one full union workspace per
static group and expect a CKTSO-sized win. The stronger target is a
producer/window-scoped grouped workspace that keeps only the live row union for
the active producer slice, then publishes or restores the nonidentical current
outputs from that owner.

The current live-workspace follow-up refines that surface to producer-step
granularity. `KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_LIVE_STEP_STATS=1`
counts, for each multi-current producer step, the duplicated retained rows that
the scalar current-state path would carry and the exact unique row union that a
single grouped producer/window workspace would need. The optional
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_LIVE_STEP_PLAN=1` additionally keeps
the per-step unique-row sizing prefix; it is still an allocation substrate only
and leaves numeric execution unchanged. Default runs keep the counters and plan
empty.

Validation passed `cmake --build build --target kls_bench kls_smoke -j2`,
`./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_LIVE_STEP_STATS=1 ./build/kls_smoke`,
`KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_LIVE_STEP_PLAN=1 ./build/kls_smoke`,
and `ctest --test-dir build --output-on-failure`. The default `ASIC_100ks`
control
`build/kls_asic100ks_live_step_default_t4_r1_ref1.json` measured
`initial_factor_seconds=0.673518`, `refactor_seconds_avg=0.0463456`, clean
residual `1.92251861e-15`, and zero live-step counters. With live-step stats
enabled, `build/kls_asic100ks_live_step_stats_t4_r1_ref1.json` reported
`68,467` multi-current live producer steps, `1,217,923` active current
memberships, `1,183,247,933` duplicated state rows, `183,365,884` unique rows,
and `999,882,049` reusable rows, for a `6.45x` row-collapse surface. The max
step reached `362` active currents, `489,578` duplicated rows, and `27,386`
unique rows. `ASIC_320ks` in
`build/kls_asic320ks_live_step_stats_t4_r1_ref1.json` reported `98,964` live
steps, `1,123,018` current memberships, `960,528,608` duplicated rows,
`120,767,097` unique rows, and `839,761,511` reusable rows, for a `7.95x`
collapse surface, with max fanout `444`. The plan variant on `ASIC_100ks`
matched those counters and stayed residual-clean.

This is the strongest paper-aligned refactor signal so far: the useful sharing
surface is much larger at active producer-step scope than at whole static group
scope, while the default KLS path remains unchanged. The next implementation
should therefore build a real producer/window live workspace from this sizing
prefix, including row-position storage and a numeric owner that streams each
producer row once across the active current set. It should not spend more time
on larger whole-group state plans or per-current retained-state replay.

That prefix has now been extended into the first real row descriptor substrate.
With `KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_LIVE_STEP_PLAN=1`, KLS retains
the sorted unique row list for each live producer step when memory allows, and
reports `refactor_btf_scalar_run_group_live_step_stored_rows` plus
`refactor_btf_scalar_run_group_live_step_storage_limited`. Correctness passed
`git diff --check`, `python3 -m py_compile scripts/decompose_solver_gap.py`,
`cmake --build build --target kls_bench kls_smoke -j2`, `./build/kls_smoke`,
both live-step opt-in smoke modes, `ctest --test-dir build --output-on-failure`,
and the decomposition CSV smoke. The default `ASIC_100ks` control
`build/kls_asic100ks_live_step_rows_default_t4_r1_ref1.json` still left the
descriptor inactive (`built=0`, live-step unique/stored rows both zero) and
measured `refactor_seconds_avg=0.0444901` with clean residual.

The full descriptor retained on both hard ASIC probes without hitting the
storage guard. `build/kls_asic100ks_live_step_rows_t4_r1_ref1.json` stored all
`183,365,884` unique live-step rows (`storage_limited=0`) and stayed
residual-clean (`1.92251861e-15`), while
`build/kls_asic320ks_live_step_rows_t4_r1_ref1.json` stored all `120,767,097`
unique rows with `storage_limited=0` and residual `2.08436857e-15`. These
plan-mode runs are intentionally expensive (`22.93s` and `17.62s` average
refactor time) because they still recompute and retain the descriptor rather
than executing through it. The value is not current speed; it is that the
producer/window owner now has the exact row sets it needs on the slow cases.
The next step is to add per-current row positions or a dense workspace mapping
on top of this descriptor and then move one producer-step numeric update into
the grouped owner.

For benchmark selection, matrices where the reference solver times out under
the same limit, such as `Hamrle3` under the current CKTSO cap, should remain in
a separate stress/scalability bucket. They are useful for robustness and
long-run profiling, but they should not drive CKTSO-relative tuning or geomean
gap claims because there is no finite CKTSO timing to close. The primary
gap-closing slice should use cases where CKTSO completes and KLS still loses.

The next row-descriptor chunk connects that retained live-step substrate to
runtime producer batches. `KLS_ENABLE_REFACTOR_BTF_SCALAR_RUN_GROUP_LIVE_STEP_PLAN=1`
now leaves the existing numeric path unchanged but records how many
producer-step retained-state batches exactly cover all memberships of a planned
live step. Benchmark JSON exposes last and cumulative
`refactor_*_live_step_runtime_full_*` counters for full steps, current
memberships, and unique descriptor rows. If these counters are high on the
slow ASIC rows, the next grouped owner can use the retained live-step row list
directly; if they are low, the remaining gap is scheduler coverage before
numeric ownership. This avoids another per-current replay tweak and tests the
paper-level grouped workspace precondition directly.

The first focused `ASIC_100ks` probe with live-step plan plus retained
producer-step advance stayed residual-clean (`1.92251861e-15`) and retained the
full `183,365,884` live-step row descriptor, but runtime exact coverage was
small: `2,723` full steps, `7,301` currents, and `1,625,695` unique descriptor
rows, while the producer batch advanced `59,689` steps and `308,194` current
updates. This makes the next missing large piece clearer: a grouped owner must
either claim/schedule whole live steps more often or handle partial-step
membership directly. The existing per-current retained-state path is not close
enough to the paper algorithm just because the row descriptor exists.
