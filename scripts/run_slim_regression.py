#!/usr/bin/env python3
"""Paired before/after regression audit using the full paper headline protocol."""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
import pathlib
import platform
import resource
import shutil
import statistics
import subprocess
import time

from run_paper_campaign import (
    clean_environment, command_for, manifest_names, parse_json_output,
    resolve_matrices, result_valid, sha256,
)
from slim_build_provenance import require_matching_builds


def reduce_records(out: pathlib.Path) -> dict:
    records = [json.loads(line) for line in (out / "observations.jsonl").read_text().splitlines()]
    groups: dict[tuple, dict] = {}
    for row in records:
        key = (row["config"], row["matrix"])
        groups.setdefault(key, {})[(row["pass"], row["side"])] = row
    rows = []
    for (config, matrix), group in groups.items():
        ratios = []
        statuses = {side: [] for side in ("before", "after")}
        for p in sorted({key[0] for key in group}):
            pair = [group.get((p, side)) for side in ("before", "after")]
            for side, record in zip(statuses, pair):
                statuses[side].append(record["status"] if record else "pending")
            if all(r and r["status"] == "valid" for r in pair):
                ratios.append(pair[1]["result"]["measured_lifecycle_seconds"] /
                              pair[0]["result"]["measured_lifecycle_seconds"])
        row = {"config": config, "matrix": matrix, "valid_pairs": len(ratios),
               "after_over_before": statistics.median(ratios) if ratios else None,
               "before_statuses": ";".join(statuses["before"]),
               "after_statuses": ";".join(statuses["after"])}
        for side in statuses:
            valid = [r["result"] for (_, s), r in group.items()
                     if s == side and r["status"] == "valid"]
            row[side + "_seconds"] = statistics.median(
                r["measured_lifecycle_seconds"] for r in valid) if valid else None
            row[side + "_paths"] = ";".join(sorted({str(r.get("last_refactor_path")) for r in valid}))
        rows.append(row)
    if rows:
        with (out / "comparison.csv").open("w") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    summary = {}
    for config in sorted({r["config"] for r in rows}):
        selected = [r for r in rows if r["config"] == config]
        valid = [r for r in selected if r["after_over_before"] is not None]
        summary[config] = {
            "matrices_observed": len(selected), "matrices_with_valid_pairs": len(valid),
            "geomean_after_over_before": math.exp(statistics.mean(
                math.log(r["after_over_before"]) for r in valid)) if valid else None,
            "regressions_over_5_percent": sum(r["after_over_before"] > 1.05 for r in valid),
            "improvements_over_5_percent": sum(r["after_over_before"] < 1 / 1.05 for r in valid),
            "worst": sorted(valid, key=lambda r: r["after_over_before"], reverse=True)[:15],
        }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--before", type=pathlib.Path, required=True)
    ap.add_argument("--after", type=pathlib.Path, required=True)
    ap.add_argument("--before-revision", required=True)
    ap.add_argument("--output", type=pathlib.Path, required=True)
    ap.add_argument("--passes", type=int, default=8)
    ap.add_argument("--timeout", type=float, default=900)
    ap.add_argument("--memory-gib", type=int, default=80)
    ap.add_argument("--reduce-only", action="store_true")
    args = ap.parse_args()
    root = pathlib.Path(__file__).resolve().parents[1]
    out = args.output.resolve()
    if args.reduce_only:
        print(json.dumps(reduce_records(out), indent=2)); return
    if args.passes < 1 or args.timeout <= 0 or args.memory_gib <= 0:
        ap.error("passes, timeout, and memory limit must be positive")
    manifests = [root / "bench" / f"suitesparse_paper_{tier}_manifest.txt"
                 for tier in ("medium", "large")]
    names = [name for path in manifests for name in manifest_names(path)]
    assert len(names) == len({n.casefold() for n in names}) == 110
    matrices = resolve_matrices(root / "data/suitesparse", names)
    configs = [c for c in json.loads((root / "bench/paper_campaign_configs.json").read_text())
               if c["id"].startswith("headline_")]
    assert len(configs) == 2
    binaries = {"before": args.before.resolve(), "after": args.after.resolve()}
    try:
        build_provenance = require_matching_builds(binaries)
    except ValueError as error:
        raise SystemExit(str(error)) from error
    env = clean_environment({})
    diff = subprocess.check_output(["git", "diff", "HEAD"], cwd=root)
    metadata = {
        "before_revision": args.before_revision,
        "after_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
        "after_diff_sha256": hashlib.sha256(diff).hexdigest(),
        "binaries": {s: {"path": str(p), "sha256": sha256(p)} for s, p in binaries.items()},
        "build_flags": {s: (p.parent / "CMakeFiles/kls.dir/flags.make").read_text()
                        for s, p in binaries.items()},
        "build_provenance": build_provenance,
        "configs": configs, "manifests": {str(p): sha256(p) for p in manifests},
        "matrices": names, "passes": args.passes, "timeout": args.timeout,
        "memory_gib": args.memory_gib, "kernel": platform.platform(),
        "cpu_affinity": sorted(os.sched_getaffinity(0)), "runner_sha256": sha256(pathlib.Path(__file__)),
        "method": "alternating paired launches; median paired H100 ratios; failures retained",
    }
    out.mkdir(parents=True, exist_ok=True)
    meta_path = out / "metadata.json"
    if meta_path.exists():
        if json.loads(meta_path.read_text()) != metadata:
            raise SystemExit("metadata changed; use a new output directory")
    else:
        meta_path.write_text(json.dumps(metadata, indent=2) + "\n")
        (out / "source.diff").write_bytes(diff)
    for side, binary in binaries.items():
        frozen = out / (side + "-kls_bench")
        if frozen.exists():
            if sha256(frozen) != metadata["binaries"][side]["sha256"]:
                raise SystemExit("frozen binary hash mismatch")
        else:
            shutil.copy2(binary, frozen)
        binaries[side] = frozen
    observations = out / "observations.jsonl"
    done = set()
    if observations.exists():
        done = {(r["config"], r["matrix"], r["pass"], r["side"])
                for r in map(json.loads, observations.read_text().splitlines())}
    limit = args.memory_gib * 1024 ** 3

    def memory_limit() -> None:
        resource.setrlimit(resource.RLIMIT_AS, (limit, limit))

    total = len(names) * len(configs) * args.passes * 2
    with observations.open("a", buffering=1) as stream:
        for config in configs:
            for matrix_index, name in enumerate(names):
                for p in range(args.passes):
                    order = ("before", "after") if (p + matrix_index) % 2 == 0 else ("after", "before")
                    for position, side in enumerate(order):
                        key = (config["id"], name, p, side)
                        if key in done: continue
                        command = ["taskset", "-c", config["cpus"]] + command_for(
                            "kls", binaries[side], matrices[name], config)
                        record = {"config": config["id"], "matrix": name, "pass": p,
                                  "side": side, "position": position, "command": command,
                                  "utc": dt.datetime.now(dt.timezone.utc).isoformat()}
                        start = time.monotonic()
                        try:
                            result = subprocess.run(command, env=env, text=True, capture_output=True,
                                                    timeout=args.timeout, preexec_fn=memory_limit)
                            record.update(exit_code=result.returncode, stderr_tail=result.stderr[-6000:])
                            try:
                                record["result"] = parse_json_output(result.stdout)
                                valid, why = result_valid(record["result"], 1e-8, config["lifecycle_systems"])
                                record["status"] = "valid" if result.returncode == 0 and valid else f"exit:{result.returncode}/{why}"
                            except ValueError:
                                record["status"] = f"exit:{result.returncode}/no-json"
                        except subprocess.TimeoutExpired:
                            record.update(status="timeout", exit_code=None)
                        record["wall_seconds"] = time.monotonic() - start
                        stream.write(json.dumps(record) + "\n")
                        done.add(key)
                        print(f"{len(done)}/{total} {config['id']} {name} p{p} {side}: "
                              f"{record['status']} ({record['wall_seconds']:.2f}s)", flush=True)
                reduce_records(out)
    print(json.dumps(reduce_records(out), indent=2))


if __name__ == "__main__":
    main()
