#!/usr/bin/env python3
"""Run the frozen, counterbalanced KLS external-validation campaign."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import pathlib
import platform
import re
import shlex
import subprocess
import sys
from typing import Any


PROTOCOL = "kls-external-validation-v1"
WORKLOADS = (
    ("rank_p001", "rank-preserving", 0.001),
    ("entry_p001", "entrywise", 0.001),
    ("entry_p100", "entrywise", 0.1),
    ("localized_p100", "localized-entrywise", 0.1),
)
DEFAULT_THREADS = (1, 4, 8)
SOLVER_SOURCE_PATHS = ("CMakeLists.txt", "include", "src", "third_party")


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def manifest_header(path: pathlib.Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        if not raw.startswith("#"):
            continue
        key, separator, value = raw[1:].strip().partition(":")
        if separator:
            result[key.strip()] = value.strip()
    return result


def manifest_names(path: pathlib.Path) -> list[str]:
    names: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        token = raw.split("#", 1)[0].strip()
        if token:
            names.append(pathlib.Path(token.split()[0]).stem)
    return names


def git_output(root: pathlib.Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(root), *args], text=True
    ).strip()


def assert_frozen_source(root: pathlib.Path, expected_commit: str) -> str:
    commit = git_output(
        root, "log", "-1", "--format=%H", "--", *SOLVER_SOURCE_PATHS
    )
    if commit != expected_commit:
        raise ValueError(
            "solver tree changed after reveal: last source commit "
            f"is {commit}, expected {expected_commit}"
        )
    dirty = subprocess.run(
        [
            "git", "-C", str(root), "diff", "--quiet", "HEAD", "--",
            *SOLVER_SOURCE_PATHS,
        ],
        check=False,
    )
    if dirty.returncode != 0:
        raise ValueError("tracked solver source changed after reveal")
    listing = git_output(root, "ls-tree", "-r", "HEAD", "--", *SOLVER_SOURCE_PATHS)
    return hashlib.sha256((listing + "\n").encode()).hexdigest()


def cpu_model() -> str:
    try:
        for line in pathlib.Path("/proc/cpuinfo").read_text().splitlines():
            if line.lower().startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "unknown"


def clean_environment() -> dict[str, str]:
    blocked_prefixes = ("KLS_", "CKTSO_", "SUBTREELU_", "BENCH_")
    environment = {
        key: value for key, value in os.environ.items()
        if not key.startswith(blocked_prefixes)
    }
    environment.update(
        {
            "OPENBLAS_NUM_THREADS": "1",
            "OMP_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "KLS_BENCH_VERIFY_EACH_REFACTOR": "1",
            "BENCH_VERIFY_EACH_REFACTOR": "1",
        }
    )
    return environment


def parse_threads(value: str) -> tuple[int, ...]:
    threads = tuple(int(item) for item in value.split(",") if item.strip())
    if not threads or any(item <= 0 for item in threads) or len(set(threads)) != len(threads):
        raise argparse.ArgumentTypeError("threads must be unique positive integers")
    return threads


def campaign_description(args: argparse.Namespace, binaries: dict[str, pathlib.Path]) -> dict[str, Any]:
    header = manifest_header(args.manifest)
    if header.get("Protocol") != PROTOCOL:
        raise ValueError(f"manifest does not declare protocol {PROTOCOL}")
    expected_commit = header.get("Frozen KLS commit", "")
    expected_tree = header.get("Frozen solver-tree SHA256", "")
    if not re.fullmatch(r"[0-9a-f]{40}", expected_commit):
        raise ValueError("manifest has no valid frozen KLS commit")
    actual_tree = assert_frozen_source(args.root, expected_commit)
    if actual_tree != expected_tree:
        raise ValueError("frozen solver-tree digest does not match manifest")
    missing = [str(path) for path in binaries.values() if not path.is_file()]
    if missing:
        raise ValueError("missing benchmark binary: " + ", ".join(missing))
    wanted = manifest_names(args.manifest)
    found: dict[str, list[pathlib.Path]] = {name.casefold(): [] for name in wanted}
    for matrix in args.matrix_dir.rglob("*.mtx"):
        key = matrix.stem.casefold()
        if key in found:
            found[key].append(matrix)
    absent = [name for name in wanted if not found[name.casefold()]]
    duplicate = [name for name in wanted if len(found[name.casefold()]) > 1]
    if absent:
        raise ValueError("missing validation matrices: " + ", ".join(absent))
    if duplicate:
        raise ValueError("duplicate validation matrix basenames: " + ", ".join(duplicate))
    return {
        "protocol": PROTOCOL,
        "frozen_kls_commit": expected_commit,
        "frozen_solver_tree_sha256": actual_tree,
        "manifest": str(args.manifest.resolve()),
        "manifest_sha256": sha256(args.manifest),
        "matrix_count": len(wanted),
        "matrix_dir": str(args.matrix_dir.resolve()),
        "threads": list(args.threads),
        "workloads": [
            {"id": ident, "values": values, "amplitude": amplitude}
            for ident, values, amplitude in WORKLOADS
        ],
        "passes": args.passes,
        "solve_repetitions": args.solve_repeat,
        "refactor_repetitions": args.refactor_repeat,
        "timeout_seconds": args.timeout,
        "residual_limit": args.residual_limit,
        "input_index": args.input_index,
        "binaries": {
            name: {"path": str(path.resolve()), "sha256": sha256(path)}
            for name, path in binaries.items()
        },
        "host": {
            "cpu": cpu_model(),
            "platform": platform.platform(),
            "python": platform.python_version(),
        },
        "fixed_environment": {
            "OPENBLAS_NUM_THREADS": "1",
            "OMP_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "KLS_BENCH_VERIFY_EACH_REFACTOR": "1",
            "BENCH_VERIFY_EACH_REFACTOR": "1",
        },
    }


def write_or_verify_metadata(path: pathlib.Path, description: dict[str, Any]) -> None:
    canonical = json.dumps(description, indent=2, sort_keys=True) + "\n"
    if path.exists():
        if path.read_text(encoding="utf-8") != canonical:
            raise ValueError(
                f"{path} describes a different campaign; use a new output directory"
            )
    else:
        path.write_text(canonical, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=pathlib.Path, default=pathlib.Path("."))
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--matrix-dir", type=pathlib.Path, required=True)
    parser.add_argument("--kls-bench", type=pathlib.Path, required=True)
    parser.add_argument("--cktso-compare", type=pathlib.Path, required=True)
    parser.add_argument("--subtreelu-compare", type=pathlib.Path, required=True)
    parser.add_argument("--klu-compare", type=pathlib.Path, required=True)
    parser.add_argument("--output-dir", type=pathlib.Path, required=True)
    parser.add_argument("--threads", type=parse_threads, default=DEFAULT_THREADS)
    parser.add_argument("--passes", type=int, default=4)
    parser.add_argument("--solve-repeat", type=int, default=100)
    parser.add_argument("--refactor-repeat", type=int, default=20)
    parser.add_argument("--timeout", type=float, default=300.0)
    parser.add_argument("--residual-limit", type=float, default=1e-8)
    parser.add_argument("--input-index", choices=("32", "64"), default="64")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--only",
        action="append",
        default=[],
        help="run only configuration ids containing this substring; repeatable",
    )
    args = parser.parse_args()
    args.root = args.root.resolve()
    args.manifest = args.manifest.resolve()
    args.matrix_dir = args.matrix_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    if args.passes <= 0 or args.passes % 4 != 0:
        parser.error("--passes must be a positive multiple of four")
    if args.solve_repeat <= 0 or args.refactor_repeat < 2:
        parser.error("--solve-repeat must be positive and --refactor-repeat at least two")
    if args.timeout <= 0 or args.residual_limit <= 0:
        parser.error("timeout and residual limit must be positive")
    binaries = {
        "kls": args.kls_bench.resolve(),
        "ck": args.cktso_compare.resolve(),
        "st": args.subtreelu_compare.resolve(),
        "klu": args.klu_compare.resolve(),
    }
    try:
        description = campaign_description(args, binaries)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        write_or_verify_metadata(args.output_dir / "campaign.json", description)
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        print(exc, file=sys.stderr)
        return 1

    runner = args.root / "scripts" / "run_paired_suite.sh"
    environment = clean_environment()
    configurations: list[str] = []
    for threads in args.threads:
        for workload_id, values, amplitude in WORKLOADS:
            ident = f"t{threads}_{workload_id}"
            if args.only and not any(token in ident for token in args.only):
                continue
            configurations.append(ident)
            outputs = {
                name: args.output_dir / f"{ident}.{name}.jsonl"
                for name in binaries
            }
            complete = args.output_dir / f"{ident}.complete"
            if complete.exists():
                print(f"skip complete {ident}")
                continue
            cmd = [
                str(runner),
                str(binaries["kls"]),
                str(binaries["ck"]),
                str(args.matrix_dir),
                str(args.manifest),
                str(outputs["kls"]),
                str(outputs["ck"]),
                str(binaries["st"]),
                str(outputs["st"]),
                str(binaries["klu"]),
                str(outputs["klu"]),
            ]
            config_env = dict(environment)
            config_env.update(
                {
                    "THREADS": str(threads),
                    "PASSES": str(args.passes),
                    "COUNTERBALANCE_SIDES": "1",
                    "STOP_INVALID_SIDES": "1",
                    "TIMEOUT": str(args.timeout),
                    "RESIDUAL_LIMIT": str(args.residual_limit),
                    "SOLVE_REPEAT": str(args.solve_repeat),
                    "REFACTOR_REPEAT": str(args.refactor_repeat),
                    "REFACTOR_VALUES": values,
                    "REFACTOR_VALUE_AMPLITUDE": str(amplitude),
                    "FACTOR_REPEAT": "0",
                    "INPUT_INDEX": args.input_index,
                }
            )
            print(f"run {ident}: {shlex.join(cmd)}")
            if args.dry_run:
                continue
            result = subprocess.run(cmd, env=config_env, check=False)
            if result.returncode != 0:
                print(f"configuration failed: {ident}", file=sys.stderr)
                return result.returncode
            complete.write_text("complete\n", encoding="utf-8")

    selection_path = args.output_dir / "selected_configurations.json"
    selection_path.write_text(
        json.dumps(configurations, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
