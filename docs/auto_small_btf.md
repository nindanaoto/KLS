# Bounded small-BTF AUTO selection

This focused change starts at main f3e6c41. It reuses the native factor format
and existing compact-index cache. It does not incorporate the accumulated
forest, register-kernel, or SIMD-batching experiments.

After the ordinary orientation/BTF competition selects an AMD winner, default
componentwise AUTO may make one METIS symbolic proposal. The envelope is
transpose storage, 2..4096 unknowns, at most 65535 entries, multiple BTF blocks,
known full structural rank, finite positive AMD estimates, and estimated work
below 100000. Explicit ordering, scale, backend and tuning-profile choices
retain their existing selection behavior. Competing frame rewrites decline.

The proposal needs at least 20% predicted-fill reduction. A known METIS work
estimate must not increase; only the documented -1 sentinel represents an
unavailable estimate. These are proposal signals, not promises of realized
fill. Internal P/Q/R permutations are owned by the selected symbolic; the
public coordinate frame and canonical input pattern remain unchanged.

A private lease retains the original AMD symbolic/common/separator plus an
exact pattern snapshot. Actual public numeric failure may consume it once
and refactor current A using AMD. Invalid arguments never trigger a retry.
Changing the symbolic, input frame, or analysis revokes the lease. Numeric
caches are discarded before rollback. Proposal time and correctness recovery
are charged to their initiating operation; no tuning solve is replayed.
Rollback also restores public symbolic block, rank, and estimated-work
statistics from AMD. Allocation-failure and singular-retry tests check these
statistics alongside the restored symbolic ownership.

The compact factor/transpose solve uses current values and RHS on every call.
Native factors remain authoritative. Only immutable capability is cached;
current ownership and masked-nearest FP controls remain executor gates. The
ordinary raw entry retains complete admission. Unsupported modes and optional
metadata failures use the original paths. Retained compact metadata and the
saved symbolic/pattern are included in the LLC/32 budget.

## Validation and merge decision

Merged into `main` on 2026-11-04 as `f7b3319f`. The merged production files
match the qualified source, public headers are unchanged, and all 17 new
policy/gate tests pass from the merged tree.

The latest user-approved gate, `auto-btf-v3`, allows a cell to regress against
main in linear and wall time when its paired linear-time medians strictly beat
both CKTSO and SubtreeLU. Otherwise the v2 limits below apply. The inverter
target, overall geometric mean, and correctness remain unchanged. The broad
gate uses the separately approved all-recorded-pairs reassessment below.
The versioned v3 gate and independent auditor agree on both retained SPICE
campaigns, which pass. All ten gate tests pass on host and container. The
user authorizes merging after the required checks pass. All current gates
pass, including the approved broad reassessment described below.

The user-approved revised gate (`auto-btf-v2`) supersedes the historical
`auto-btf` decision policy below. Both serial inverter cells must improve
linear time by at least 2% versus main. A cell where main already wins against
either peer must preserve each such win and stay within 2% of main in linear
and wall time. A cell where main loses to both peers may regress by at most
5% in linear and wall time. The geometric mean of the 28 candidate/main
linear-time cell medians must not exceed 1.0. Existing peer ties are protected
as wins. Classification uses medians of paired main/peer observations, never
a quotient of candidate/peer and candidate/main medians. Correctness and full
grid requirements are unchanged. The original scripts and failed receipts
remain intact; versioned launch/review and independent audit scripts implement
the revised policy. Seven focused tests check both implementations and raw
paired classification. The broad matrix comparison retains its original
correctness coverage and 2% lifecycle threshold; on 2026-11-04 the user
approved using all recorded pairs for the nine investigated cells, with
the existing both-peer-win exception.

The latest scientific commit, 9d8f9350, passes all 51 native, portable, and
ASan/UBSan tests, plus four scoped TSan ownership tests. The revised gate's
seven focused Python tests pass on the host and in the container. Its ordinary
API and retained-execution campaigns each pass 70 observations, including
independent agreement and retained execution in all four inverter cells.
The `latest_validation` receipt identifies this exact frozen source.

Both fresh opposite-order timing campaigns pass correctness on all 462
observations each, with exact independent agreement. Initial serial inverter
ratios versus main are 0.89252 and 0.89364; confirmation ratios are 0.90360 and
0.89953. Their overall linear-time geometric means are 0.98365 and 0.98355.
The initial campaign passes `auto-btf-v2`. Confirmation fails only the wall-time
limit for LLC32 serial `ram2k`: its linear ratio is 1.00263 and wall ratio is
1.03062, with no lost peer win. That v2 failure is retained. The user subsequently approved the v3 peer-win
exception: ram2k confirmation ratios are 0.90259 against CKTSO and 0.44521
against SubtreeLU. Both v3 reassessments pass, and the combined receipt verifies matching inputs,
binaries, and opposite campaign orders. The full
1760-observation broad comparison has completed and its independent audit
agrees exactly. It has no candidate-only correctness failure, with 80 shared
invalid observations retained. Nine qualified matrix/cache groups exceed the
unchanged 1.02 lifecycle ratio limit; the largest is LLC32 `1138_bus` at
1.09224. The original four-pair broad gate fails and is retained as a historical
decision. The approved reassessment below passes. No
complete-corpus geometric mean is claimed from the incomplete correctness coverage.

Follow-up diagnostics preserve every original result. A 216-observation,
12-pair screen found two cells still above 2%; a 324-observation, balanced
three-arm experiment found no consistent benefit from keeping native and
compact refactor routines out of line. Those experimental annotations were
not promoted to production.

The peer benchmark tools now use timing schema 2: summed elapsed public API
calls, matching the KLS driver. Internal library timers remain diagnostic
fields. The previous peer tools summed internal timers and cannot establish
this public-call comparison. Fresh frozen peer tools completed 432 balanced
observations across all nine original over-limit groups; all raw logs,
residuals, timing schemas, and paired ratios were revalidated. Five groups
beat both peers: LLC32 rajat21, adder_dcop_01, case9, powersim, and LLC96
adder_dcop_01. The four others lose to CKTSO.

Pooling all 40 recorded main/candidate pairs per investigated group retains
the original failed campaign and every diagnostic pair. The four groups
without both peer wins then have ratios 1.00642 (LLC32 1138_bus), 1.00849
(LLC96 circuit_3), 1.00989 (LLC96 hcircuit), and 1.00508 (LLC96
legresley_2508). Only case9 remains over 2%, at 1.02121, and beats both
peers. On 2026-11-04 the user approved using all 40 recorded pairs instead
of the frozen four-pair gate. `reassess_auto_broad_pooled.py` validates the
original full-grid audit, exact binary identities, all raw logs and numerical
checks, all launch orders, and all 40 pairs in every original over-limit cell.
It preserves the original failure and passes the approved gate. An independent
reduction agrees exactly. Seven gate tests pass on host and container. The
original 1,760 observations, including 80 shared failures, remain intact.
No complete-corpus aggregate is claimed. The unmodified production sources
match the build qualified in all four profiles; public headers are unchanged.

The previous scientific commit, 74365dfc, passes all 51 native, portable, and
ASan/UBSan tests, plus four scoped TSan ownership tests. Its frozen source and
library receipts are recorded under `current_validation` in the qualification
receipt. Its ordinary API and retained AUTO execution campaigns both pass 70
observations with exact independent agreement. Its first full timing campaign
passes correctness on all 462 observations, but fails the accepted performance
gate: LLC96 eight-thread `dac` has a 1.02568 linear-time ratio versus main and
does not qualify for a peer exception. The serial inverter ratios are 0.87594
and 0.88033. Confirmation and broad comparison were not launched for that
snapshot. The failed campaign is retained without filtering.

`dac` has 6366 unknowns and never enters compact execution. The next revision
moves existing size limits ahead of ownership and FP-environment admission,
including ordinary refactor/solve entries. A 12-pair diagnostic per thread count
measures size-first/previous median ratios of 0.99187 (one thread) and 0.98935
(eight threads), with all 49 observations qualified. These diagnostics support
the rejection-cost change; they do not override the failed gate. Complete fresh
qualification is required for the revised source. The measurements below
describe the earlier 98a3893e snapshot and do not establish current merge
readiness.

The dedicated private test retains exact native-factor and vendor-solve word
and FP-flag oracles. The real one-RHS vendor transpose loops explicitly
round the multiply before subtraction, matching their existing native release
code generation and the compact executor. A register barrier prevents the
sanitizer build from silently contracting those operations; complex and
multiple-RHS walks retain their existing arithmetic.

A captured caller-order inverter fixture exercises positive
selection through the public API, unlike the grid rejection fixtures. Tests
cover changed A/RHS, aliasing, both solve directions, numeric allocation failure,
singular updates, snapshot allocation failure, lease revocation and reanalysis.
The adjacent fixture provenance records source and capture hashes.

The frozen main-based candidate at 98a3893e passes native, portable and
ASan/UBSan suites (49 tests each), plus all four scoped TSan ownership tests.
The unchanged exact oracle also passes against fresh main's unmodified native
vendor archive. Both the ordinary API and retained-execution campaigns pass
70 observations, with exact independent agreement. All four inverter cells
retain AUTO-selected METIS and execute compact factor and solve operations on
successful public calls; mux8 retains no proposal.

Both opposite-order timing campaigns pass: 924 qualified observations with
exact independent agreement on every comparison cell. The serial inverter
linear-time ratios versus fresh main are:

| Campaign | LLC 96 MiB | LLC 32 MiB |
| --- | ---: | ---: |
| Initial | 0.90795 | 0.89920 |
| Confirmation | 0.88952 | 0.88886 |

Serial inverter linear time is 26.3–28.8% lower than SubtreeLU's and
3.7–6.2% higher than CKTSO's in these campaigns. The existing CKTSO gap remains;
this focused change does not complete the broader peer-win objective.

All other SPICE cells remain within the 2% linear/wall non-regression tolerance
in both orders. No peer exception is needed for those campaigns.

The broad experiment stopped during FullChip teardown after 464 recorded
observations. Its strict cleanup rejected a transient process in the completed
session; a later strict identity check removed only the remaining mount helper.
The original campaign and its unrecorded attempted stage are retained. The
independently audited complete prefix contains no candidate-only correctness
loss, but circuit_1 has a 1.02594 median lifecycle ratio, exceeding the 1.02
broad gate. This candidate is therefore **not merge ready**. A production fix
and fresh qualification are required; the failing measurements are not filtered
or accepted through additional repetitions. A subsequent implementation moves
normal-orientation rejection ahead of the added small-BTF admission checks.
The focused native checks pass. A separate 72-observation diagnostic using an
identical benchmark object places both the previous candidate and the new guard
close to main, indicating startup sensitivity rather than proving a systematic
2.6% cost. That screen does not replace the failed acceptance measurements;
the new implementation requires complete fresh qualification.

A separate bounded post-exit wait preserves the strict cleanup's PID, session,
owner and image checks. It never signals an unknown process and never extends
the numerical acceptance deadline. Three real-session tests and five prefix
preservation tests pass, including both container CTest registrations. The new
broad harness persists numeric exit status before teardown. Continuation can
retain an entire immutable complete-matrix prefix only with identical frozen
binaries and inputs, including every failed row. It was not launched for this
candidate because the audited prefix already fails the lifecycle gate.

The broad ledger includes shared 900-second Hamrle3 timeouts. Those pairs supply
no performance ratio; no complete-corpus aggregate is claimed. The
[qualification receipt](evidence/auto_small_btf_qualification.json) records frozen
sources, libraries, receipts, the failed broad gate and remaining work. The
[development record](evidence/auto_small_btf_development.json) preserves the
initial ASan exact-word failure and earlier excluded preparations.

The `auto-btf` campaign policy uses two fresh opposite-order full campaigns.
Both one-thread inverter cells must improve linear time by at least 2% versus
fresh main. Other cells may regress beyond the 2% linear/wall tolerance only
when that same cell remains at least 2% faster in linear time than each of
CKTSO and SubtreeLU. Missing, failed, duplicate, diagnostic or incomplete grids
cannot pass. The independent auditor reconstructs the gate without importing
the primary reducer. Existing non-regressed peer gaps do not block this focused
change; the broader all-56 peer-win objective is not claimed complete.
