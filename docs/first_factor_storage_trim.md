# First-factor storage accounting trial

Baseline: `0839ac6`; cumulative reference: pushed `88a177e`.

The trial removes private allocation-growth/copy counters from first-factor
entry storage, their worker snapshots, trace aggregation fields, collector,
and printed fields. Source-wide references show no allocator or numerical
policy consumes those counts. Allocation sizes, copying, ownership, failure
handling, and all other pipeline telemetry are unchanged. Source reduction:
80 net lines, but private structure sizes change and require performance
validation rather than an assumption of identical generated behavior.

Release and ASan/UBSan CTest pass 6/6. Artifacts are under
`build/prep-trim-repair-dBcoha/storage-*`. The three-side comparison rotates
the preceding, starting, and candidate binaries, checks build provenance,
records commands/hashes/source diff, and verifies H100 entrywise refactors
at a 1e-8 residual limit. Timed runs do not overlap builds or tests.

Besides the automatic-policy controls, explicit `KLS_ENABLE_KLS_FIRST_FACTOR=1`
tests cover 1138_bus, circuit_4 and bcircuit; all report `kls_first` initial
factor paths. Their four-triplet screen passes 36/36 launches but flags
1138_bus at +2.97% versus the preceding commit and +3.74% versus the pushed
starting point. This is near the stopping threshold and needs a repeated
check before accepting the chunk.

The repeat uses 24 rotating triplets on each of two eight-core CPU sets
(144/144 valid launches). Median paired lifecycle changes are -1.00% and
+0.30% versus `0839ac6`, and -0.76% and -0.07% versus `88a177e`.
The initial small-sample slowdown did not reproduce, so this chunk is retained.
These forced-first fixtures report zero parallel row-pipeline rows; they
exercise KLS-first but do not establish timed parallel-pipeline coverage.

The automatic-policy campaign passes 282/282 launches across seven matrices
and ten configurations. Incremental lifecycle changes range from -0.53% to
+1.09%; cumulative changes range from -0.78% to +1.37%. TSOPF_FS_b9_c6
changes by -0.53%, +0.05%, and +0.56% for the two eight-thread sets and the
single-thread control, respectively. No notable regression was reproduced.

Exact evidence: `storage-trim-{metadata,summary}.json`, `storage-trim.jsonl`,
`storage-first-recheck-{metadata,summary}.json`, and
`storage-first-recheck.jsonl` under the artifact directory above. This is
focused H100 validation, not a full paper or Xyce benchmark rerun.
