# Resolve an add32 corpus-screen launch confound

The second medium-corpus checkpoint compares `78a4c89` with `372f9d1`:
467 fewer net C-source lines over the resumed pass. Its runner is
`build/resumed-corpus-chhiOV/scan-producer.py`, with artifacts prefixed
`producer-corpus0`. It paused after case 23 (`add32`, CPUs 0-7): 24 cases,
168 valid launches, 162 cases remaining. No source changes were pending.

The initial screen measured +8.35% lifecycle time. Its 12-pair repeat was
+7.13% (12/12 slower) and +17.80% steady refactor time. Median lifecycles
were 5.876 ms baseline and 6.308 ms current. This is a large relative flag
on a short workload and warranted investigation, not automatic acceptance.

## Confirmation and conflicting evidence

- A separate 24-triplet comparison with the preceding `2a35dcd` binary
  measured +0.04% incrementally and +0.79% cumulatively on CPUs 0-7.
  Twelve triplets on CPUs 8-15 measured +0.18% / -0.48%.
- Reusing the exact corpus executable paths and alternating pair order
  reproduced +7.29% over 24 pairs (23/24 slower), followed by +7.27% over
  12 pairs (12/12 slower). Steady refactor changes were +18.53% / +19.67%.
- Under the same paired runner, using the focused runner's byte-identical
  candidate copy instead measured -0.70% over 24 pairs. Both candidate
  copies have SHA-256
  `8dd9eea773608a6dcb526c6c630a5efd369a298a707554db5d0cb6e130ea41b5`.
- Finally, baseline and current were copied sequentially to one shared
  temporary executable path, checking its hash before every launch.
  This equal-path 24-pair comparison measured +0.46%, below the review
  threshold. Frozen input executables were not overwritten.

These checks demonstrate executable-copy/path sensitivity in this setup;
they do not identify its lower-level cause. In particular, the raw +7%
result cannot be attributed confidently to source slimming. Do not call
this a confirmed code regression or conceal the original flag. No numerical
failures occurred. Continue the corpus screen from case 24, retaining the
frozen candidate and treating further flags with the same scrutiny.

## Artifacts

Corpus directory: `build/resumed-corpus-chhiOV/`:
`producer-corpus0-*`, `producer-add32-exact-repeat-*`,
`producer-add32-copy-repeat-*`, `producer-add32-equal-repeat-*`, and their
`scan-add32*.py` runners. The shared temporary executable is
`add32-equal-path-kls_bench`. The 108-launch three-binary control is under
`build/prep-trim-repair-dBcoha/producer-add32-stopcheck*`.
The original screen and these controls total 444 valid launches. Exact
commands, binary hashes, build provenance and empty source diffs are saved.

## Second flag and revised corpus method

The resumed segment `producer-corpus1` completed cases 24-40, then paused
on `hcircuit`, CPUs 8-15. Its initial three-pair screen measured +3.02%;
the 12-pair repeat measured +3.08% lifecycle and +7.19% steady refactor,
with all 12 pairs slower. A separate 24-pair shared-executable-path check
measured only +0.08%, with all launches valid. Its runner and artifacts are
`scan-hcircuit-equal.py` and `producer-hcircuit-equal-repeat-*` in the corpus
directory. This second flag also does not establish a source regression.

Two resolved flags show that separate executable paths are a material
confound for this checkpoint. Restart the full 186-case comparison using
`scan-producer-equal.py`, tag `producer-equal0`, from case 0. This is a
deliberate method correction after terminal review pauses, not a restart
caused by an observation timeout. Rechecking earlier cases is necessary
because launch effects can conceal slowdowns as well as create flags.

Before each sequential launch, copy the chosen frozen input binary to
`corpus-equal-path-kls_bench` and verify its SHA-256. Both sides execute at
that same path with the same arguments. Copies and hashes happen before
the benchmark process starts, outside its lifecycle timer. Input binaries
remain untouched. The manifest, two CPU domains, H100 workload, numerical
verification, limits and repeat rule remain unchanged. The new metadata
records the execution-copy policy. Further source deletions stay paused
until this corrected checkpoint is resolved.
