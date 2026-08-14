#!/usr/bin/env python3
"""Run post-reveal KLS route diagnostics on the frozen validation suite.

These runs explain AUTO regret; they are never part of the confirmatory
solver comparison and must not be used to change the frozen implementation.
"""

from __future__ import annotations

import argparse
import json
import os
import pathlib
import subprocess
import sys
from typing import Any

from run_external_validation import clean_environment, manifest_names


ARMS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("auto", ("--backend", "auto")),
    ("serial_backend", ("--backend", "serial")),
    ("row_off", ("--backend", "auto", "--row-refactor", "off")),
    ("no_fast_factor", ("--backend", "auto", "--no-fast-factor")),
    ("no_btf", ("--backend", "auto", "--no-btf")),
    ("forced_amd", ("--backend", "auto", "--ordering", "amd")),
    ("forced_metis", ("--backend", "auto", "--ordering", "metis")),
)


def locate_matrices(matrix_dir: pathlib.Path, names: list[str]) -> dict[str, pathlib.Path]:
    by_name: dict[str, pathlib.Path] = {}
    duplicates: set[str] = set()
    for path in matrix_dir.rglob("*.mtx"):
        key = path.stem.casefold()
        if key in by_name:
            duplicates.add(key)
        else:
            by_name[key] = path
    if duplicates:
        raise ValueError("duplicate matrix basenames: " + ", ".join(sorted(duplicates)))
    missing = [name for name in names if name.casefold() not in by_name]
    if missing:
        raise ValueError("missing validation matrices: " + ", ".join(missing))
    return {name: by_name[name.casefold()] for name in names}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign", type=pathlib.Path, required=True)
    parser.add_argument("--matrix-dir", type=pathlib.Path, required=True)
    parser.add_argument("--kls-bench", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--passes", type=int, default=7)
    parser.add_argument("--timeout", type=float, default=300.0)
    parser.add_argument(
        "--values",
        choices=("rank-preserving", "entrywise", "localized-entrywise"),
        default="entrywise",
    )
    parser.add_argument("--amplitude", type=float, default=0.001)
    parser.add_argument("--refactor-repeat", type=int, default=20)
    parser.add_argument("--solve-repeat", type=int, default=100)
    parser.add_argument(
        "--matrix",
        action="append",
        default=[],
        help="run only this manifest matrix basename; repeatable",
    )
    parser.add_argument(
        "--arm",
        action="append",
        choices=tuple(arm for arm, _ in ARMS),
        default=[],
        help="run only this diagnostic arm; repeatable",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.threads <= 0 or args.passes <= 0:
        parser.error("--threads and --passes must be positive")
    if args.timeout <= 0 or not 0 < args.amplitude < 1:
        parser.error("timeout and amplitude are invalid")
    try:
        campaign: dict[str, Any] = json.loads(args.campaign.read_text(encoding="utf-8"))
        manifest = pathlib.Path(campaign["manifest"])
        names = manifest_names(manifest)
        if args.matrix:
            requested = {pathlib.Path(name).stem.casefold() for name in args.matrix}
            available = {name.casefold() for name in names}
            missing = sorted(requested - available)
            if missing:
                raise ValueError(
                    "requested matrices are not in the campaign: " + ", ".join(missing)
                )
            names = [name for name in names if name.casefold() in requested]
        matrices = locate_matrices(args.matrix_dir, names)
        expected_binary = campaign["binaries"]["kls"]["sha256"]
        from run_external_validation import sha256
        if sha256(args.kls_bench) != expected_binary:
            raise ValueError("KLS benchmark binary differs from the confirmatory campaign")
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(exc, file=sys.stderr)
        return 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    environment = clean_environment()
    arms = tuple(item for item in ARMS if not args.arm or item[0] in args.arm)
    mode = "w"
    if args.dry_run:
        output = None
    else:
        output = args.output.open(mode, encoding="utf-8")
    try:
        for matrix_index, name in enumerate(names):
            for pass_index in range(args.passes):
                shift = (pass_index + matrix_index) % len(arms)
                order = arms[shift:] + arms[:shift]
                for arm, extra in order:
                    cmd = [
                        str(args.kls_bench), str(matrices[name]), "--json",
                        "--threads", str(args.threads),
                        "--repeat", str(args.solve_repeat),
                        "--factor-repeat", "0",
                        "--refactor-repeat", str(args.refactor_repeat),
                        "--refactor-values", args.values,
                        "--refactor-value-amplitude", str(args.amplitude),
                        "--no-transpose-solve", "--input-index", "64",
                        *extra,
                    ]
                    print(f"{name} pass={pass_index} arm={arm}")
                    if args.dry_run:
                        continue
                    try:
                        proc = subprocess.run(
                            cmd, text=True, capture_output=True, check=False,
                            timeout=args.timeout, env=environment,
                        )
                    except subprocess.TimeoutExpired:
                        row: dict[str, Any] = {"status": "timeout_or_fail"}
                    else:
                        if proc.returncode == 0:
                            try:
                                row = json.loads(proc.stdout)
                            except json.JSONDecodeError:
                                row = {"status": "parse_error"}
                        else:
                            row = {"status": "solver_error", "returncode": proc.returncode}
                    row["matrix"] = str(matrices[name])
                    row["validation_arm"] = arm
                    row["validation_pass"] = pass_index
                    output.write(json.dumps(row, sort_keys=True) + "\n")
                    output.flush()
    finally:
        if output is not None:
            output.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
