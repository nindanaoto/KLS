#!/usr/bin/env python3
"""Run kls_bench over a directory of MatrixMarket files and summarize results."""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import subprocess
import sys


def spice_cycle_seconds(row: dict[str, object]) -> float:
    analysis = float(row["analysis_seconds"])
    factor = float(row.get("initial_factor_seconds", row["factor_seconds_avg"]))
    refactor = float(row["refactor_seconds_avg"])
    solve = float(row["solve_seconds_avg"])
    return analysis + factor + solve + 99.0 * (refactor + solve)


def geometric_mean(values: list[float]) -> float:
    positives = [v for v in values if v > 0.0 and math.isfinite(v)]
    if not positives:
        return math.nan
    return math.exp(sum(math.log(v) for v in positives) / len(positives))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kls-bench", type=pathlib.Path, required=True)
    parser.add_argument("--matrix-dir", type=pathlib.Path, required=True)
    parser.add_argument("--jsonl", type=pathlib.Path)
    parser.add_argument("--repeat", type=int, default=5)
    parser.add_argument("--refactor-repeat", type=int, default=5)
    parser.add_argument("--passes", type=int, default=1)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--timeout", type=float)
    parser.add_argument("--ordering", choices=["auto", "amd", "colamd", "natural", "metis"], default="auto")
    parser.add_argument("--orientation", choices=["auto", "normal", "transpose"], default="auto")
    parser.add_argument("--scale", choices=["auto", "-1", "0", "1", "2"], default="auto")
    parser.add_argument("--pivot-tol", type=float, default=None)
    parser.add_argument("--no-btf", action="store_true")
    parser.add_argument("--no-fast-factor", action="store_true")
    parser.add_argument("--no-static-pivoting", action="store_true")
    args = parser.parse_args()

    matrices = sorted(args.matrix_dir.rglob("*.mtx"))
    if not matrices:
        print(f"no .mtx files found under {args.matrix_dir}", file=sys.stderr)
        return 1
    if args.passes <= 0:
        print("--passes must be positive", file=sys.stderr)
        return 1

    rows: list[dict[str, object]] = []
    failures: list[tuple[pathlib.Path, str]] = []
    if args.jsonl is not None:
        args.jsonl.parent.mkdir(parents=True, exist_ok=True)
        out = args.jsonl.open("w", encoding="utf-8")
    else:
        out = None

    try:
        for matrix in matrices:
            samples: list[dict[str, object]] = []
            sample_failures: list[str] = []
            for _ in range(args.passes):
                cmd = [
                    str(args.kls_bench),
                    str(matrix),
                    "--repeat",
                    str(args.repeat),
                    "--refactor-repeat",
                    str(args.refactor_repeat),
                    "--threads",
                    str(args.threads),
                    "--ordering",
                    args.ordering,
                    "--orientation",
                    args.orientation,
                    "--scale",
                    args.scale,
                    "--json",
                ]
                if args.no_btf:
                    cmd.append("--no-btf")
                if args.no_fast_factor:
                    cmd.append("--no-fast-factor")
                if args.no_static_pivoting:
                    cmd.append("--no-static-pivoting")
                if args.pivot_tol is not None:
                    cmd.extend(["--pivot-tol", str(args.pivot_tol)])
                try:
                    proc = subprocess.run(
                        cmd,
                        text=True,
                        capture_output=True,
                        check=False,
                        timeout=args.timeout,
                    )
                except subprocess.TimeoutExpired:
                    sample_failures.append(f"timeout after {args.timeout:g}s")
                    continue
                if proc.returncode != 0:
                    sample_failures.append(proc.stderr.strip())
                    continue
                sample = json.loads(proc.stdout)
                sample["spice_cycle_seconds"] = spice_cycle_seconds(sample)
                samples.append(sample)
            if not samples:
                failures.append((matrix, "; ".join(sample_failures)))
                continue
            samples.sort(key=lambda row: float(row["spice_cycle_seconds"]))
            row = dict(samples[len(samples) // 2])
            cycle_samples = [float(sample["spice_cycle_seconds"]) for sample in samples]
            row["spice_cycle_seconds_samples"] = cycle_samples
            row["spice_cycle_seconds_median"] = row["spice_cycle_seconds"]
            row["passes_ok"] = len(samples)
            row["passes_failed"] = len(sample_failures)
            rows.append(row)
            if out is not None:
                out.write(json.dumps(row, sort_keys=True) + "\n")
                out.flush()
            if args.passes == 1:
                print(f"{matrix}: {row['spice_cycle_seconds']:.6g}s")
            else:
                print(
                    f"{matrix}: {row['spice_cycle_seconds']:.6g}s "
                    f"median of {len(samples)} pass(es)"
                )
    finally:
        if out is not None:
            out.close()

    cycle_values = [float(row["spice_cycle_seconds"]) for row in rows]
    summary = {
        "matrices_total": len(matrices),
        "matrices_ok": len(rows),
        "matrices_failed": len(failures),
        "spice_cycle_geomean_seconds": geometric_mean(cycle_values),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    if failures:
        print("failures:", file=sys.stderr)
        for matrix, err in failures:
            print(f"  {matrix}: {err}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
