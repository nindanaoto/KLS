# Remove producer rejection-reason diagnostics

Preceding baseline: `688e2ab`; resumed-pass baseline: pushed `372f9d1`.
Retain the SNB alignment fix and previous accepted slimming.

Remove eight diagnostic rejection counters, their aggregation/printing,
the reason-to-counter switch, and the private rejection enum. Only NOT_ROOT
was consumed outside tracing. Preserve it as the scalar admission helper's
optional `not_root_out` boolean: initialize false on entry, set true only
at the original heap-root mismatch, and test it at the unchanged catch-up
decision. The retry passes NULL because its reason is no longer needed.
The supernode helper needs only its existing success/failure return.

All admission predicates, their evaluation order, helper calls, locking,
and numerical update behavior are preserved. The output boolean remains
false for every earlier safety rejection. This removes 152 net C-source
lines without deleting catch-up or producer batching. The summary script
keeps archived rejection-field support and documents its historical status.

## Correctness evidence

Release and ASan/UBSan CTest each pass 6/6. A differential harness extracts
the actual pre-change and post-change scalar/supernode admission bodies,
uses small fixture structs and controlled readiness helpers, and compares
917,504 cases over combinations of null inputs, missing workspace fields,
epoch mismatch, ownership, invalid state, absent dependencies, non-root
dependencies, readiness, supernodes, cached panels and dependency end bounds.
It verifies identical admission results and equivalence of the old NOT_ROOT
reason to the new boolean, including NULL output handling. All pass.
This is predicate/flag coverage with stubbed readiness, not a replacement
for full asynchronous pipeline tests. Harness and log:
`build/prep-trim-repair-dBcoha/check-reject-admission.py` and
`reject-admission-check.log`.

An additional release smoke run enables pipeline tracing, active catch-up,
and supernode producer batching. General trace parsing is retained. Timed
campaigns run only after builds and correctness checks finish. This run
passes, emits 24 general progress records and no retired rejection fields,
and reports up to 160 active-catch-up attempts and 160 successful catch-up
targets. Thus catch-up is actually exercised, not merely enabled; the run
does not establish executed supernode producer batches (batch counts are 0).

## Performance validation

Artifacts: `build/prep-trim-repair-dBcoha/reject-trim-*`. The rotating H100
runner compares the preceding binary (`original`), aligned starting binary
(`start`), and candidate (`fix`), with pinned CPU sets, clean environments,
verified build provenance, source diff, commands/hashes, and every
entrywise-changed refactor checked at a 1e-8 residual limit. No builds/tests
overlap timed runs. This is not a full paper/Xyce or other-machine campaign.

All **435/435 timed launches pass**: 36 onetone1 comparisons (six triplets
on each cache domain), 351 automatic-policy controls across six other
matrices/nine configurations, and 48 forced-KLS-first controls. The recorded
campaign diff matches the final source patch.

Median paired H100 changes versus `688e2ab` / `372f9d1`:

| Matrix / configuration | Incremental | Resumed-pass cumulative |
|---|---:|---:|
| onetone1, CPU 8, one thread | +0.30% | +0.21% |
| onetone1, CPU 0, one thread | -0.01% | +0.08% |
| TSOPF_FS_b9_c6, CPUs 8-15, eight threads | +0.61% | +0.20% |
| TSOPF_FS_b9_c6, CPUs 0-7, eight threads | -0.30% | +0.66% |
| TSOPF_FS_b9_c6, one thread | -0.11% | +0.21% |
| transient, one thread | -0.03% | +0.19% |
| twotone, eight threads | +0.71% | +1.26% |
| ASIC_680ks, eight threads | -0.83% | -1.73% |

Forced-first changes span -1.16% to +0.06% incrementally. No focused case
meets the approximately 3% repeatable lifecycle stop condition. Retain this
chunk. The SNB entry remains 64-byte aligned (0x64dc0, size 0x1931) and its
recovered onetone1 performance is retained. Final Release/sanitizer builds
have no warnings or errors. These observations do not prove zero regression
on every paper matrix or another architecture.

After three resumed-pass cleanup chunks (270 net C-source lines), broaden
the cumulative paper-corpus screen against `372f9d1` before further source
deletions. The existing focused set cannot establish full-corpus neutrality.
Further candidates include producer probe/candidate trace counters, but
their associated operational admission quantities must be preserved.
The broader cleanup goal remains active.
