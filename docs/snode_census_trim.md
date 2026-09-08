# Remove hypothetical supernode censuses

Baseline: `55d1f6b`. This chunk removes 174 net source lines from
`kls_numeric_schedule.inc` without changing supernode construction.

Removed: the `KLS_TRACE_RELAXED_SNODE` symmetric-difference census, the
`KLS_TRACE_SNODE` run-width histogram, and its hypothetical subset/padding
simulation. All scan results were local and used only in diagnostic output.
The scans do not select a representation or update numerical data.

Retained: strict run detection, sorting, padded-panel construction, coverage
admission, elapsed-time accounting, preparation-decline diagnostics, and a
basic sorted/covered/accepted preparation summary. Numeric thresholds and
all SNB/row/residual kernels are unchanged.

Release and ASan/UBSan CTest each pass 6/6. Performance artifacts are under
`build/prep-trim-repair-dBcoha/census-trim*`. The comparison alternates the
frozen `trace-trim-kls_bench` (`55d1f6b` equivalent) with the candidate, using
matching build provenance, automatic policies and verified H100 entrywise
refactors. It records commands, source diff, binary/matrix hashes, and phase
timings. Positive paired percentage changes mean slower.

The focused campaign completed 188 launches across seven matrices
and eleven configurations, emphasizing both TSOPF CPU sets and large SNB
serial controls. Builds/tests do not overlap timed comparisons. This is not
a full paper or Xyce rerun.

All **188/188 launches pass**. No notable regression was observed; paired
median lifecycle changes range from -1.76% to +1.04%. TSOPF is +0.51% on
CPUs 8-15 and +0.58% on CPUs 0-7. Serial onetone1 is -0.77%, twotone +0.20%,
transient +0.20%, and ASIC_680ks +0.10%. These are focused measurements,
not a guarantee of identical performance for every input.

The actual supernode code path and numerical safeguards remain intact.
The next independently reviewed low-risk candidate is the diagnostic-only
stored-zero/run/tail-padding block in compact-solve preparation.
