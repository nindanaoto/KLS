#!/usr/bin/env python3
"""Run a resumable, pass-paired KLS paper campaign.

The output is one append-only JSONL stream.  Every attempted launch carries a
configuration, matrix, pass, solver, and launch-position identity, including
timeouts and invalid solves.  No failed side is silently dropped or retried.
"""

from __future__ import annotations

import argparse
import collections
import datetime as dt
import hashlib
import json
import math
import os
import pathlib
import platform
import shlex
import socket
import subprocess
import sys
import time
from typing import Any


SCHEMA = "kls-paper-observation-v2"
SOLVERS = ("kls", "ck", "st", "klu")
BALANCED_ORDERS = (
    ("kls", "ck", "st", "klu"),
    ("ck", "kls", "klu", "st"),
    ("st", "klu", "kls", "ck"),
    ("klu", "st", "ck", "kls"),
)


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git(root: pathlib.Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(root), *args], text=True, stderr=subprocess.DEVNULL
    ).strip()


def manifest_names(path: pathlib.Path) -> list[str]:
    names: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        token = raw.split("#", 1)[0].strip()
        if token:
            names.append(pathlib.Path(token.split()[0]).stem)
    if not names or len({name.casefold() for name in names}) != len(names):
        raise ValueError(f"{path}: manifest is empty or has duplicate basenames")
    return names


def resolve_matrices(matrix_dir: pathlib.Path, names: list[str]) -> dict[str, pathlib.Path]:
    wanted = {name.casefold(): name for name in names}
    found: dict[str, list[pathlib.Path]] = {key: [] for key in wanted}
    for path in matrix_dir.rglob("*.mtx"):
        key = path.stem.casefold()
        if key in found:
            found[key].append(path.resolve())
    missing = [wanted[key] for key, paths in found.items() if not paths]
    duplicate = [wanted[key] for key, paths in found.items() if len(paths) > 1]
    if missing or duplicate:
        parts = []
        if missing:
            parts.append("missing: " + ", ".join(missing))
        if duplicate:
            parts.append("duplicate: " + ", ".join(duplicate))
        raise ValueError("; ".join(parts))
    return {wanted[key]: paths[0] for key, paths in found.items()}


def clean_environment(extra: dict[str, Any]) -> dict[str, str]:
    blocked = ("KLS_", "CKTSO_", "SUBTREELU_", "BENCH_", "OMP_", "MKL_")
    env = {key: value for key, value in os.environ.items()
           if not key.startswith(blocked)}
    env.update({
        "OPENBLAS_NUM_THREADS": "1",
        "OMP_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "KLS_BENCH_VERIFY_EACH_REFACTOR": "1",
        "BENCH_VERIFY_EACH_REFACTOR": "1",
    })
    env.update({str(key): str(value) for key, value in extra.items()})
    return env


def command_for(solver: str, binary: pathlib.Path, matrix: pathlib.Path,
                config: dict[str, Any]) -> list[str]:
    threads = int(config["threads"])
    systems = int(config.get("lifecycle_systems", 100))
    mode = str(config.get("values", "entrywise"))
    amplitude = str(config.get("amplitude", 0.001))
    if solver == "kls":
        command = [
            str(binary), str(matrix), "--lifecycle-systems", str(systems),
            "--threads", str(threads), "--backend", "auto", "--ordering", "auto",
            "--orientation", "auto", "--scale", "auto", "--input-index",
            str(config.get("input_index", "64")), "--refactor-values", mode,
            "--refactor-value-amplitude", amplitude, "--no-transpose-solve", "--json",
        ]
        if "expected_refactors" in config:
            command += ["--expected-refactors", str(config["expected_refactors"])]
        if "expected_solves" in config:
            command += ["--expected-solves", str(config["expected_solves"])]
        command += [str(item) for item in config.get("kls_args", [])]
        return command
    if solver in ("ck", "st"):
        return [
            str(binary), str(matrix), str(threads), "1", str(max(systems - 1, 0)),
            "0", mode, amplitude, "--lifecycle-systems", str(systems),
        ]
    return [
        str(binary), str(matrix), "--lifecycle-systems", str(systems),
        "--refactor-values", mode, "--refactor-value-amplitude", amplitude, "--json",
    ] + [str(item) for item in config.get("klu_args", [])]


def result_valid(result: dict[str, Any], limit: float, systems: int) -> tuple[bool, str]:
    if result.get("status") not in (None, 0):
        return False, f"solver-status:{result.get('status')}"
    if result.get("verify_each_refactor") is not True:
        return False, "refactors-not-verified"
    if result.get("lifecycle_mode") != "direct" or result.get("lifecycle_systems") != systems:
        return False, "wrong-lifecycle"
    for key in ("relative_residual_l2", "refactor_max_relative_residual"):
        try:
            value = float(result[key])
        except (KeyError, TypeError, ValueError):
            return False, f"missing:{key}"
        if not (math.isfinite(value) and 0.0 <= value <= limit):
            return False, f"residual:{key}={value:.6g}"
    try:
        elapsed = float(result["measured_lifecycle_seconds"])
    except (KeyError, TypeError, ValueError):
        return False, "missing:measured_lifecycle_seconds"
    if not (math.isfinite(elapsed) and elapsed > 0.0):
        return False, "invalid:measured_lifecycle_seconds"
    return True, "valid"


def parse_json_output(stdout: str) -> dict[str, Any]:
    for raw in reversed(stdout.splitlines()):
        try:
            value = json.loads(raw)
        except json.JSONDecodeError:
            continue
        if isinstance(value, dict):
            return value
    raise ValueError("no JSON object in stdout")


def existing_keys(path: pathlib.Path) -> set[tuple[str, str, int, str]]:
    keys: set[tuple[str, str, int, str]] = set()
    if not path.exists():
        return keys
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip():
            continue
        row = json.loads(raw)
        key = (row["config_id"], row["matrix_id"], int(row["pass_id"]), row["solver"])
        if key in keys:
            raise ValueError(f"{path}:{line_number}: duplicate observation {key}")
        keys.add(key)
    return keys


def metadata(root: pathlib.Path, manifest: pathlib.Path, configs: pathlib.Path,
             binaries: dict[str, pathlib.Path], passes: int,
             residual_limit: float) -> dict[str, Any]:
    diff = subprocess.check_output(
        ["git", "-C", str(root), "diff", "--binary", "HEAD"],
    )
    try:
        lscpu = subprocess.check_output(["lscpu", "--json"], text=True)
        lscpu_value: Any = json.loads(lscpu)
    except (OSError, subprocess.CalledProcessError, json.JSONDecodeError):
        lscpu_value = None
    runner = pathlib.Path(__file__).resolve()
    reducer = runner.with_name("reduce_paper_campaign.py")
    scientific_status = git(
        root, "status", "--short", "--untracked-files=all", "--",
        "CMakeLists.txt", "include", "src", "third_party", "bench", "scripts",
    ).splitlines()
    return {
        "schema": "kls-paper-campaign-v2",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "git_commit": git(root, "rev-parse", "HEAD"),
        "git_diff_sha256": hashlib.sha256(diff).hexdigest(),
        "scientific_status": scientific_status,
        "runner_sha256": sha256(runner),
        "reducer_sha256": sha256(reducer),
        "manifest": str(manifest), "manifest_sha256": sha256(manifest),
        "manifest_entries": manifest_names(manifest),
        "configs": str(configs), "configs_sha256": sha256(configs),
        "configuration_spec": json.loads(configs.read_text(encoding="utf-8")),
        "passes": passes, "residual_limit": residual_limit,
        "binaries": {key: {"path": str(path), "sha256": sha256(path)}
                     for key, path in binaries.items()},
        "host": {"hostname": socket.gethostname(), "platform": platform.platform(),
                 "python": platform.python_version(), "lscpu": lscpu_value},
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=pathlib.Path, default=pathlib.Path("."))
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--matrix-dir", type=pathlib.Path, required=True)
    parser.add_argument("--configs", type=pathlib.Path, required=True)
    parser.add_argument("--output-dir", type=pathlib.Path, required=True)
    for solver in SOLVERS:
        parser.add_argument(f"--{solver}-bench", type=pathlib.Path, required=True)
    parser.add_argument("--passes", type=int, default=8)
    parser.add_argument("--timeout", type=float, default=900.0)
    parser.add_argument("--residual-limit", type=float, default=1e-8)
    parser.add_argument("--only", action="append", default=[])
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.passes <= 0 or args.passes % 4:
        parser.error("--passes must be a positive multiple of four")
    if args.timeout <= 0 or args.residual_limit <= 0:
        parser.error("timeout and residual limit must be positive")
    root, manifest, matrix_dir, configs_path, output_dir = (
        args.root.resolve(), args.manifest.resolve(), args.matrix_dir.resolve(),
        args.configs.resolve(), args.output_dir.resolve(),
    )
    binaries = {solver: getattr(args, f"{solver}_bench").resolve() for solver in SOLVERS}
    try:
        configs = json.loads(configs_path.read_text(encoding="utf-8"))
        if not isinstance(configs, list) or not configs:
            raise ValueError("configs must be a nonempty JSON array")
        names = manifest_names(manifest)
        matrices = resolve_matrices(matrix_dir, names)
        ids = [str(config["id"]) for config in configs]
        if len(set(ids)) != len(ids):
            raise ValueError("configuration ids are not unique")
        by_id = {str(config["id"]): config for config in configs}
        for config in configs:
            config_id = str(config["id"])
            if int(config.get("threads", 0)) <= 0:
                raise ValueError(f"{config_id}: threads must be positive")
            if int(config.get("lifecycle_systems", 100)) <= 0:
                raise ValueError(f"{config_id}: lifecycle_systems must be positive")
            baseline = config.get("baseline")
            if baseline is not None:
                if str(baseline) not in by_id:
                    raise ValueError(f"{config_id}: unknown baseline {baseline}")
                if config.get("pair_group") != by_id[str(baseline)].get("pair_group"):
                    raise ValueError(f"{config_id}: baseline must share pair_group")
        output_dir.mkdir(parents=True, exist_ok=True)
        campaign_path = output_dir / "campaign.json"
        description = metadata(root, manifest, configs_path, binaries,
                               args.passes, args.residual_limit)
        if campaign_path.exists():
            previous = json.loads(campaign_path.read_text(encoding="utf-8"))
            for key in ("git_commit", "git_diff_sha256", "scientific_status",
                        "runner_sha256", "reducer_sha256", "manifest_sha256",
                        "configs_sha256", "passes", "residual_limit", "binaries"):
                if previous.get(key) != description.get(key):
                    raise ValueError(f"campaign metadata changed at {key}; use a new output directory")
        else:
            campaign_path.write_text(json.dumps(description, indent=2, sort_keys=True) + "\n")
        observations = output_dir / "observations.jsonl"
        done = existing_keys(observations)
    except (OSError, ValueError, KeyError, json.JSONDecodeError,
            subprocess.CalledProcessError) as exc:
        print(exc, file=sys.stderr)
        return 1

    selected_configs = [
        config for config in configs
        if not args.only or any(token in str(config["id"]) for token in args.only)
    ]
    if not selected_configs:
        print("--only did not select any configuration", file=sys.stderr)
        return 1
    groups: dict[str, list[dict[str, Any]]] = collections.OrderedDict()
    for config in selected_configs:
        config_id = str(config["id"])
        selected = list(config.get("matrices", names))
        selected_solvers = tuple(config.get("solvers", SOLVERS))
        if any(name not in matrices for name in selected):
            raise ValueError(f"{config_id}: unknown matrix in selection")
        if any(solver not in SOLVERS for solver in selected_solvers):
            raise ValueError(f"{config_id}: unknown solver")
        group_id = str(config.get("pair_group", f"config:{config_id}"))
        groups.setdefault(group_id, []).append(config)

    with observations.open("a", encoding="utf-8", buffering=1) as stream:
        for group_id, group_configs in groups.items():
            selections = [list(config.get("matrices", names)) for config in group_configs]
            if any(selected != selections[0] for selected in selections[1:]):
                raise ValueError(f"{group_id}: paired configurations require identical matrix order")
            selected = selections[0]
            if any(name not in matrices for name in selected):
                raise ValueError(f"{group_id}: unknown matrix in selection")
            for matrix_id in selected:
                for pass_id in range(args.passes):
                    offset = pass_id % len(group_configs)
                    config_order = group_configs[offset:] + group_configs[:offset]
                    config_ids = [str(item["id"]) for item in config_order]
                    for config_position, config in enumerate(config_order):
                        config_id = str(config["id"])
                        selected_solvers = tuple(config.get("solvers", SOLVERS))
                        order = [solver for solver in BALANCED_ORDERS[pass_id % 4]
                                 if solver in selected_solvers]
                        for launch_position, solver in enumerate(order):
                            key = (config_id, matrix_id, pass_id, solver)
                            if key in done:
                                continue
                            command = command_for(solver, binaries[solver], matrices[matrix_id], config)
                            cpus = str(config.get("cpus", "")).strip()
                            if cpus:
                                command = ["taskset", "-c", cpus, *command]
                            print(f"[{config_id} {matrix_id} p{pass_id} {solver}] {shlex.join(command)}")
                            if args.dry_run:
                                continue
                            started = dt.datetime.now(dt.timezone.utc).isoformat()
                            start = time.monotonic()
                            timed_out = False
                            try:
                                completed = subprocess.run(
                                    command, env=clean_environment(config.get("env", {})),
                                    text=True, capture_output=True, timeout=args.timeout,
                                    check=False,
                                )
                                returncode, stdout, stderr = (
                                    completed.returncode, completed.stdout, completed.stderr
                                )
                            except subprocess.TimeoutExpired as exc:
                                timed_out, returncode = True, None
                                stdout = exc.stdout if isinstance(exc.stdout, str) else ""
                                stderr = exc.stderr if isinstance(exc.stderr, str) else ""
                            result: dict[str, Any] | None = None
                            valid, reason = False, "timeout" if timed_out else "process-failed"
                            if not timed_out and returncode == 0:
                                try:
                                    result = parse_json_output(stdout)
                                    valid, reason = result_valid(
                                        result, args.residual_limit,
                                        int(config.get("lifecycle_systems", 100)),
                                    )
                                except ValueError as exc:
                                    reason = f"parse-error:{exc}"
                            row = {
                                "schema": SCHEMA, "config_id": config_id,
                                "pair_group": group_id,
                                "config_launch_position": config_position,
                                "config_launch_order": config_ids,
                                "matrix_id": matrix_id, "matrix": str(matrices[matrix_id]),
                                "pass_id": pass_id, "solver": solver,
                                "launch_position": launch_position, "launch_order": order,
                                "started_utc": started,
                                "wall_seconds": time.monotonic() - start,
                                "command": command, "returncode": returncode,
                                "timed_out": timed_out, "valid": valid,
                                "validity_reason": reason, "result": result,
                                "stderr_tail": stderr[-4000:],
                            }
                            stream.write(json.dumps(row, sort_keys=True) + "\n")
                            stream.flush()
                            os.fsync(stream.fileno())
                            done.add(key)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
