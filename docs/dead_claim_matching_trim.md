# Dead claim state and explicit matching experiments

Removed the EGraph claim/lease helpers and state: no remaining constructor
allocated claim storage, and both ownership generations were always zero.
Pipeline completion, ordinary dispatch, fused dispatch, and worker affinity
remain unchanged. Removed the explicit matching cycle-path search and its
exclusive helpers, plus the explicit unscaled static-match override. Normal
matching, pair swaps, and production deferred matching equilibration remain.

Replaced the obsolete README consumer-plan/Algorithm 5 payoff instructions
with a retirement notice. Historical reports and public reserved stats remain.

Release and ASan/UBSan CTest pass 5/5 each. This cleanup does not claim to fix
the previously measured TSOPF_FS_b9_c6 one-thread first-refactor regression;
that investigation is a separate change.
