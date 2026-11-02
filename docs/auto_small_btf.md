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

Native development tests and the initial ten-observation unhinted capture
screen passed. These are development evidence, not complete qualification or
a performance result against main. Full native/portable/ASan/scoped-TSan,
API, retained-execution and broad-corpus qualification remain required.

The `auto-btf` campaign policy uses two fresh opposite-order full campaigns.
Both one-thread inverter cells must improve linear time by at least 2% versus
fresh main. Other cells may regress beyond the 2% linear/wall tolerance only
when that same cell remains at least 2% faster in linear time than each of
CKTSO and SubtreeLU. Missing, failed, duplicate, diagnostic or incomplete grids
cannot pass. The independent auditor reconstructs the gate without importing
the primary reducer. Existing non-regressed peer gaps do not block this focused
change; the broader all-56 peer-win objective is not claimed complete.
