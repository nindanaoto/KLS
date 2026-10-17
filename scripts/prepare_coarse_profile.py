#!/usr/bin/env python3
"""Prepare an isolated diagnostic source archive; never edit solver sources.

Only complete API/engine/job functions get scopes. No column/entry clocks.
The manifest proves removing every insertion recovers the exact Git source.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import tarfile

SCOPES = {
    "kls.c": [
        ("kls_refactor_solve", "caller"),
        ("private_refactor_body", "caller"),
        ("private_factor_body", "caller"),
        ("private_numeric_snapshot_finish", "caller"),
        ("private_caller_rhs", "caller"),
        ("private_solve_body", "caller"),
        ("private_caller_finish", "caller"),
        ("private_componentwise_finish", "caller"),
        ("kls_refactor_solve_impl", "caller"),
        ("kls_lean_btf_map32_refactor", "caller"),
    ],
    "kls_original_operator.inc": [("private_original_operator_build", "caller")],
    "kls_numeric_prepare_repair.inc": [
        ("prepare_numeric_values", "caller"),
        ("kls_sync_authoritative_packed_l_values", "caller"),
        ("kls_sync_authoritative_packed_u_values", "caller"),
    ],
    "kls_egraph_refactor.inc": [
        ("kls_mapped_refactor", "caller"),
        ("kls_egraph_mapped_refactor", "caller"),
        ("kls_egraph_refactor_worker_run", "worker"),
        ("kls_egraph_refactor_pool_run_worker", "worker"),
        ("kls_row_solve_worker_run", "worker"),
        ("kls_try_parallel_row_solve_one_rhs", "caller"),
    ],
    "kls_row_refactor.inc": [
        ("kls_lean_single_block_map32_refactor", "caller"),
        ("kls_lean_row_refactor_numeric", "caller"),
        ("kls_lean_parallel_worker_run", "worker"),
        ("kls_row_refactor_worker_run", "worker"),
        ("kls_row_publish_worker_run", "worker"),
    ],
    "kls_numeric_schedule.inc": [
        ("kls_snb_refactor", "caller"),
        ("kls_snb_refactor_parallel", "caller"),
        ("kls_parallel_refactor", "caller"),
    ],
    "kls_compact_solve.inc": [
        ("kls_i32_solve", "caller"),
        ("kls_native_transpose_solve", "caller"),
    ],
}


def body_start(source, name):
    # Definition signatures start at column zero; reject prototypes explicitly.
    pattern = r"^(?:static\s+)?(?:KLS_ALWAYS_INLINE\s+)?[\w *]+\b" + re.escape(name) + r"\s*\("
    candidates = []
    for match in re.finditer(pattern, source, re.M):
        i, depth = match.end(), 1
        while depth and i < len(source):
            depth += (source[i] == "(") - (source[i] == ")")
            i += 1
        while i < len(source) and source[i].isspace():
            i += 1
        if i < len(source) and source[i] == "{":
            candidates.append(i + 1)
    if len(candidates) != 1:
        raise ValueError(f"Expected one definition of {name}, found {len(candidates)}")
    return candidates[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--revision", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    repo = Path(subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip())
    revision = subprocess.check_output(["git", "rev-parse", args.revision + "^{commit}"], text=True).strip()
    for filename, scopes in SCOPES.items():
        original = subprocess.check_output(
            ["git", "show", f"{revision}:src/{filename}"], text=True)
        for name, _ in scopes:
            body_start(original, name)
    out = args.output.resolve()
    out.mkdir()  # No pooling, overwriting, or reusing a previous preparation.
    source = out / "source"
    source.mkdir()
    archive = subprocess.check_output(["git", "archive", "--format=tar", revision])
    with tarfile.open(fileobj=io.BytesIO(archive)) as stream:
        stream.extractall(source, filter="data")
    records, edits = [], {}
    for filename, scopes in SCOPES.items():
        path = source / "src" / filename
        original = text = path.read_text()
        insertions = []
        for name, role in scopes:
            slot = len(records)
            addition = f"\n  KLS_COARSE_SCOPE({slot}); /* diagnostic: {name} */"
            offset = body_start(text, name)
            text = text[:offset] + addition + text[offset:]
            insertions.append(addition)
            records.append({"slot": slot, "name": name, "role": role, "file": filename})
        if filename == "kls.c":
            for name, anchor in (
                ("ordinary_numeric_job", "    const UF_long nblocks = shared->symbolic->nblocks;"),
                ("ordinary_certificate_job", "    if(shared->solver && shared->solver->componentwise_work) {"),
            ):
                assert text.count(anchor) == 1, (name, text.count(anchor))
                slot = len(records)
                addition = f"\n      KLS_COARSE_SCOPE({slot}); /* diagnostic: {name} */"
                if name == "ordinary_numeric_job":
                    text = text.replace(anchor, addition + "\n" + anchor)
                    insertions.append(addition + "\n")
                else:
                    text = text.replace(anchor, anchor + addition)
                    insertions.append(addition)
                records.append({"slot": slot, "name": name, "role": "worker", "file": filename})
            anchor = '#include "kls/kls.h"'
            addition = '\n#include "kls_coarse_profile.inc"'
            assert text.count(anchor) == 1
            text = text.replace(anchor, anchor + addition)
            insertions.append(addition)
        recovered = text
        for addition in reversed(insertions):
            assert recovered.count(addition) == 1
            recovered = recovered.replace(addition, "", 1)
        assert recovered == original
        path.write_text(text)
        edits[filename] = insertions
    helper = Path(__file__).with_name("coarse_profile.inc")
    names = ",".join(json.dumps(r["name"]) for r in records)
    roles = ",".join(json.dumps(r["role"]) for r in records)
    prefix = (f"#define KLS_COARSE_COUNT {len(records)}\n"
              f"static const char *const kls_coarse_names[] = {{{names}}};\n"
              f"static const char *const kls_coarse_roles[] = {{{roles}}};\n")
    (source / "src/kls_coarse_profile.inc").write_text(prefix + helper.read_text())
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    manifest = {
        "revision": revision, "inverse_exact": True, "scopes": records, "insertions": edits,
        "tool_sha256": {str(p.resolve()): sha(p) for p in (Path(__file__), helper)},
        "source_sha256": {str(p): sha(p) for p in source.rglob("*") if p.is_file()},
        "allocator_wrappers": ["malloc", "calloc", "realloc", "free", "aligned_alloc", "posix_memalign"],
        "limits": "Diagnostic only. No per-column clocks. Worker times overlap and are not caller wall time. Threads counts thread lifetimes entering a scope, not simultaneous occupancy. Allocation requests cover wrapped link boundaries, not live/process memory. No speed claim from instrumented execution.",
    }
    (out / "source-receipt.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Prepared {len(records)} coarse scopes from {revision} at {source}")


if __name__ == "__main__":
    main()
