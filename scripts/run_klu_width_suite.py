#!/usr/bin/env python3
"""Run the vendored KLU32/KLU64 diagnostic as paired width samples."""

from __future__ import annotations

import argparse
import json
import math
import os
import pathlib
import statistics
import subprocess
import sys


def geometric_mean(values: list[float]) -> float:
    positives = [value for value in values if value > 0.0 and math.isfinite(value)]
    if not positives:
        return math.nan
    return math.exp(sum(math.log(value) for value in positives) / len(positives))


def read_manifest(path: pathlib.Path) -> list[str]:
    names: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        name = raw.split("#", 1)[0].strip().lower()
        if name:
            names.append(name)
    return names


def select_matrices(
    matrix_dir: pathlib.Path, manifest: pathlib.Path
) -> list[pathlib.Path]:
    by_stem: dict[str, pathlib.Path] = {}
    duplicates: set[str] = set()
    for matrix in matrix_dir.rglob("*.mtx"):
        stem = matrix.stem.lower()
        if stem in by_stem:
            duplicates.add(stem)
        else:
            by_stem[stem] = matrix
    if duplicates:
        raise ValueError("duplicate matrix stems: " + ", ".join(sorted(duplicates)))
    selected: list[pathlib.Path] = []
    missing: list[str] = []
    for name in read_manifest(manifest):
        matrix = by_stem.get(name)
        if matrix is None:
            missing.append(name)
        else:
            selected.append(matrix)
    if missing:
        raise ValueError("missing manifest matrices: " + ", ".join(missing))
    return selected


def median_field(samples: list[dict[str, object]], width: str, field: str) -> float:
    return statistics.median(float(sample[width][field]) for sample in samples)


def summarize_matrix(matrix: pathlib.Path, samples: list[dict[str, object]]) -> dict[str, object]:
    comparable = [
        sample
        for sample in samples
        if sample.get("same_symbolic")
        and sample["klu32"]["nnz_l"] == sample["klu64"]["nnz_l"]
        and sample["klu32"]["nnz_u"] == sample["klu64"]["nnz_u"]
        and sample["klu32"]["nblocks"] == sample["klu64"]["nblocks"]
    ]
    row: dict[str, object] = {
        "matrix": str(matrix),
        "samples_ok": len(samples),
        "comparable_samples": len(comparable),
        "same_symbolic_all": len(comparable) == len(samples),
        "same_numeric_shape_all": all(
            bool(sample.get("same_numeric_shape")) for sample in samples
        ),
        "samples": samples,
    }
    for width in ("klu32", "klu64"):
        row[width] = {
            "refactor_steady_seconds_median": median_field(
                samples, width, "refactor_steady_seconds_avg"
            ),
            "solve_seconds_median": median_field(samples, width, "solve_seconds_avg"),
            "spice_cycle_seconds_median": median_field(
                samples, width, "spice_cycle_seconds"
            ),
            "max_relative_residual": max(
                float(sample[width]["refactor_max_relative_residual"])
                for sample in samples
            ),
        }
    if comparable:
        ratios = [
            float(sample["klu32"]["refactor_steady_seconds_avg"])
            / float(sample["klu64"]["refactor_steady_seconds_avg"])
            for sample in comparable
        ]
        solve_ratios = [
            float(sample["klu32"]["solve_seconds_avg"])
            / float(sample["klu64"]["solve_seconds_avg"])
            for sample in comparable
        ]
        cycle_ratios = [
            float(sample["klu32"]["spice_cycle_seconds"])
            / float(sample["klu64"]["spice_cycle_seconds"])
            for sample in comparable
        ]
        row["ratio32_over_64"] = {
            "refactor_steady_median": statistics.median(ratios),
            "solve_median": statistics.median(solve_ratios),
            "spice_cycle_median": statistics.median(cycle_ratios),
        }
    return row


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bench", type=pathlib.Path, required=True)
    parser.add_argument("--matrix-dir", type=pathlib.Path, required=True)
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--jsonl", type=pathlib.Path, required=True)
    parser.add_argument("--passes", type=int, default=3)
    parser.add_argument("--repeat", type=int, default=2)
    parser.add_argument("--factor-repeat", type=int, default=0)
    parser.add_argument("--refactor-repeat", type=int, default=7)
    parser.add_argument(
        "--refactor-values",
        choices=["unchanged", "rank-preserving", "entrywise", "localized-entrywise"],
        default="entrywise",
    )
    parser.add_argument("--refactor-value-amplitude", type=float, default=0.001)
    parser.add_argument("--ordering", choices=["amd", "colamd", "natural"], default="amd")
    parser.add_argument("--scale", choices=["-1", "0", "1", "2"], default="2")
    parser.add_argument("--no-btf", action="store_true")
    parser.add_argument("--cpus", default="0")
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument("--skip", type=int, default=0)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()

    if args.passes <= 0 or args.repeat <= 0 or args.refactor_repeat <= 1:
        parser.error("passes/repeat must be positive and refactor-repeat must exceed one")
    try:
        matrices = select_matrices(args.matrix_dir, args.manifest)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    matrices = matrices[args.skip :]
    if args.limit is not None:
        matrices = matrices[: args.limit]
    if not matrices:
        print("matrix selection is empty", file=sys.stderr)
        return 1

    args.jsonl.parent.mkdir(parents=True, exist_ok=True)
    failures_path = args.jsonl.with_suffix(".failures")
    rows: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []
    env = os.environ.copy()
    env["BENCH_VERIFY_EACH_REFACTOR"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["OMP_NUM_THREADS"] = "1"

    with args.jsonl.open("w", encoding="utf-8") as output:
        for matrix_index, matrix in enumerate(matrices):
            samples: list[dict[str, object]] = []
            errors: list[str] = []
            for pass_index in range(args.passes):
                width_order = (
                    "32-first" if (matrix_index + pass_index) % 2 == 0 else "64-first"
                )
                cmd = [
                    "taskset",
                    "-c",
                    args.cpus,
                    str(args.bench),
                    str(matrix),
                    "--repeat",
                    str(args.repeat),
                    "--factor-repeat",
                    str(args.factor_repeat),
                    "--refactor-repeat",
                    str(args.refactor_repeat),
                    "--refactor-values",
                    args.refactor_values,
                    "--refactor-value-amplitude",
                    str(args.refactor_value_amplitude),
                    "--ordering",
                    args.ordering,
                    "--scale",
                    args.scale,
                    "--width-order",
                    width_order,
                    "--json",
                ]
                if args.no_btf:
                    cmd.append("--no-btf")
                try:
                    proc = subprocess.run(
                        cmd,
                        text=True,
                        capture_output=True,
                        check=False,
                        timeout=args.timeout,
                        env=env,
                    )
                except subprocess.TimeoutExpired:
                    errors.append(f"{width_order}: timeout after {args.timeout:g}s")
                    continue
                if proc.returncode != 0:
                    errors.append(
                        f"{width_order}: exit {proc.returncode}: {proc.stderr.strip()}"
                    )
                    continue
                try:
                    sample = json.loads(proc.stdout)
                except json.JSONDecodeError as exc:
                    errors.append(f"{width_order}: invalid JSON: {exc}")
                    continue
                samples.append(sample)
            if not samples:
                failure = {"matrix": str(matrix), "errors": errors}
                failures.append(failure)
                print(f"{matrix.stem}: failed", flush=True)
                continue
            row = summarize_matrix(matrix, samples)
            if errors:
                row["errors"] = errors
            rows.append(row)
            output.write(json.dumps(row, sort_keys=True) + "\n")
            output.flush()
            ratio = row.get("ratio32_over_64", {}).get("refactor_steady_median")
            ratio_text = f"{float(ratio):.4f}" if ratio is not None else "not-comparable"
            print(f"{matrix.stem}: refactor 32/64 {ratio_text}", flush=True)

    with failures_path.open("w", encoding="utf-8") as output:
        for failure in failures:
            output.write(json.dumps(failure, sort_keys=True) + "\n")

    comparable_rows = [row for row in rows if "ratio32_over_64" in row]
    ratios = [
        float(row["ratio32_over_64"]["refactor_steady_median"])
        for row in comparable_rows
    ]
    residual_ok = all(
        float(row[width]["max_relative_residual"]) <= 1.0e-8
        for row in comparable_rows
        for width in ("klu32", "klu64")
    )
    wins = sum(ratio < 0.98 for ratio in ratios)
    losses = sum(ratio > 1.02 for ratio in ratios)
    refactor_geomean = geometric_mean(ratios)
    gate_pass = (
        len(comparable_rows) >= 8
        and refactor_geomean <= 0.90
        and wins > len(comparable_rows) // 2
        and residual_ok
    )
    summary = {
        "matrices_total": len(matrices),
        "matrices_ok": len(rows),
        "matrices_failed": len(failures),
        "matrices_comparable": len(comparable_rows),
        "refactor_ratio32_over_64_geomean": refactor_geomean,
        "wins_over_2pct": wins,
        "losses_over_2pct": losses,
        "residual_gate_pass": residual_ok,
        "backend_implementation_gate_pass": gate_pass,
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
