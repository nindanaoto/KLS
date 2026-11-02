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
