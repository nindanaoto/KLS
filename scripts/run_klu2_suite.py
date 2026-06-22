#!/usr/bin/env python3
"""Run klu2_compare over a directory of MatrixMarket files and summarize results."""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import subprocess
import sys


def geometric_mean(values: list[float]) -> float:
    positives = [v for v in values if v > 0.0 and math.isfinite(v)]
    if not positives:
        return math.nan
    return math.exp(sum(math.log(v) for v in positives) / len(positives))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--klu2-compare", type=pathlib.Path, required=True)
    parser.add_argument("--matrix-dir", type=pathlib.Path, required=True)
    parser.add_argument("--jsonl", type=pathlib.Path)
    parser.add_argument("--repeat", type=int, default=5)
    parser.add_argument("--refactor-repeat", type=int, default=5)
    parser.add_argument("--timeout", type=float)
    args = parser.parse_args()

    matrices = sorted(args.matrix_dir.rglob("*.mtx"))
    if not matrices:
        print(f"no .mtx files found under {args.matrix_dir}", file=sys.stderr)
        return 1
    if args.repeat <= 0 or args.refactor_repeat < 0:
        print("repeat arguments are invalid", file=sys.stderr)
        return 1

    rows: list[dict[str, object]] = []
    failures: list[tuple[pathlib.Path, str]] = []
    out = None
    if args.jsonl is not None:
        args.jsonl.parent.mkdir(parents=True, exist_ok=True)
        out = args.jsonl.open("w", encoding="utf-8")

    try:
        for matrix in matrices:
            cmd = [
                str(args.klu2_compare),
                str(matrix),
                str(args.repeat),
                str(args.refactor_repeat),
            ]
            try:
                proc = subprocess.run(
                    cmd,
                    text=True,
                    capture_output=True,
                    check=False,
                    timeout=args.timeout,
                )
            except subprocess.TimeoutExpired:
                failures.append((matrix, f"timeout after {args.timeout:g}s"))
                continue
            if proc.returncode != 0:
                failures.append((matrix, proc.stderr.strip()))
                continue
            row = json.loads(proc.stdout)
            rows.append(row)
            if out is not None:
                out.write(json.dumps(row, sort_keys=True) + "\n")
                out.flush()
            print(f"{matrix}: {float(row['spice_cycle_seconds']):.6g}s")
    finally:
        if out is not None:
            out.close()

    cycle_values = [float(row["spice_cycle_seconds"]) for row in rows]
    print(json.dumps({
        "matrices_total": len(matrices),
        "matrices_ok": len(rows),
        "matrices_failed": len(failures),
        "spice_cycle_geomean_seconds": geometric_mean(cycle_values),
    }, indent=2, sort_keys=True))
    if failures:
        print("failures:", file=sys.stderr)
        for matrix, reason in failures:
            print(f"  {matrix}: {reason}", file=sys.stderr)
    return 0 if rows else 1


if __name__ == "__main__":
    raise SystemExit(main())
