# Remove predicted-factor default-argument wrappers

Preceding baseline: `efc1235`; resumed-pass baseline: `372f9d1`.

Two private wrappers supplied a NULL span argument to the existing suffix
closure and pivoting fill implementations. Each had one call site. Pass NULL
directly at those sites and remove the wrappers and three now-unnecessary
forward declarations. The span-aware callers remain intact. No closure,
pivoting, allocation, recovery, or numerical acceptance logic changes.
This removes 18 net C-source lines.

Release and ASan/UBSan builds succeeded; CTest passes 6/6 in each (3.80 and
15.64 seconds). Build provenance matches the preceding accepted comparison,
and the frozen preceding executable hash matches its recorded candidate hash.

After `objcopy --strip-all --remove-section .note.gnu.build-id`, the complete
before/after executable files are byte-identical. Both SHA-256 hashes are:

`afaf9e992dc53d2c028cbbb8e8664175f2c7fa617f4974c91f407da1258079e5`

Comparison files: `build/predicted-wrapper-rBmshL/before-all` and `after-all`.
Preceding executable: `build/prep-trim-repair-dBcoha/policy-wrapper-accepted-kls_bench`.
The SNB kernel remains at 0x64540, size 0x1931, with 64-byte alignment.

No new timing campaign is needed for this chunk's unchanged runtime image;
the preceding 579-launch results remain the latest focused measurements.
This is build-specific binary equivalence, not a cross-compiler or
cross-machine performance claim. The broader cleanup goal remains active.
