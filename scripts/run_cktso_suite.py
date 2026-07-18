#!/usr/bin/env python3
"""Run cktso_compare over a directory of MatrixMarket files and summarize results."""

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


def read_manifest(path: pathlib.Path) -> list[str]:
    names: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        if line:
            names.append(line.lower())
    return names


def select_matrices(matrix_dir: pathlib.Path, manifest: pathlib.Path | None) -> list[pathlib.Path]:
    matrices = sorted(matrix_dir.rglob("*.mtx"))
    if manifest is None:
        return matrices
    wanted = read_manifest(manifest)
    matrix_by_stem: dict[str, pathlib.Path] = {}
    duplicate_stems: list[str] = []
    for matrix in matrices:
        stem = matrix.stem.lower()
        if stem in matrix_by_stem:
            duplicate_stems.append(stem)
            continue
        matrix_by_stem[stem] = matrix
    if duplicate_stems:
        raise ValueError(
            "duplicate matrix basenames under "
            f"{matrix_dir}: {', '.join(sorted(set(duplicate_stems)))}"
        )
    selected: list[pathlib.Path] = []
    missing: list[str] = []
    for name in wanted:
        matrix = matrix_by_stem.get(name)
        if matrix is None:
            missing.append(name)
        else:
            selected.append(matrix)
    if missing:
        raise ValueError(
            f"manifest entries not found under {matrix_dir}: {', '.join(missing)}"
        )
    return selected


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cktso-compare", type=pathlib.Path, required=True)
    parser.add_argument("--matrix-dir", type=pathlib.Path, required=True)
    parser.add_argument("--manifest", type=pathlib.Path)
    parser.add_argument("--jsonl", type=pathlib.Path)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--repeat", type=int, default=5)
    parser.add_argument(
        "--factor-repeat",
        type=int,
        default=None,
        help="extra Factorize calls after the initial factor (default: follows --repeat)",
    )
    parser.add_argument("--refactor-repeat", type=int, default=5)
    parser.add_argument(
        "--refactor-values",
        choices=("unchanged", "rank-preserving", "entrywise"),
        default="unchanged",
    )
    parser.add_argument(
        "--refactor-value-amplitude",
        type=float,
        default=1.0e-3,
    )
    parser.add_argument("--timeout", type=float)
    parser.add_argument(
        "--skip",
        type=int,
        default=0,
        help="skip this many selected matrices before running",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="run at most this many selected matrices",
    )
    args = parser.parse_args()

    try:
        matrices = select_matrices(args.matrix_dir, args.manifest)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    if not matrices:
        print(f"no .mtx files found under {args.matrix_dir}", file=sys.stderr)
        return 1
    if (
        args.threads <= 0
        or args.repeat <= 0
        or args.refactor_repeat < 0
        or (args.factor_repeat is not None and args.factor_repeat < 0)
        or not math.isfinite(args.refactor_value_amplitude)
        or args.refactor_value_amplitude < 0.0
        or args.refactor_value_amplitude >= 1.0
        or (
            args.refactor_values != "unchanged"
            and args.refactor_value_amplitude <= 0.0
        )
    ):
        print("threads and repeat arguments are invalid", file=sys.stderr)
        return 1
    if args.skip < 0:
        print("--skip must be non-negative", file=sys.stderr)
        return 1
    if args.limit is not None and args.limit < 0:
        print("--limit must be non-negative", file=sys.stderr)
        return 1
    if args.skip:
        matrices = matrices[args.skip :]
    if args.limit is not None:
        matrices = matrices[: args.limit]
    if not matrices:
        print("matrix selection is empty after applying --skip/--limit", file=sys.stderr)
        return 1

    rows: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []
    out = None
    failure_out = None
    if args.jsonl is not None:
        args.jsonl.parent.mkdir(parents=True, exist_ok=True)
        out = args.jsonl.open("w", encoding="utf-8")
        failure_out = args.jsonl.with_suffix(".failures").open("w", encoding="utf-8")

    try:
        for matrix in matrices:
            cmd = [
                str(args.cktso_compare),
                str(matrix),
                str(args.threads),
                str(args.repeat),
                str(args.refactor_repeat),
                str(
                    args.factor_repeat
                    if args.factor_repeat is not None
                    else args.repeat
                ),
                args.refactor_values,
                f"{args.refactor_value_amplitude:.17g}",
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
                failure = {
                    "matrix": str(matrix),
                    "reason": f"timeout after {args.timeout:g}s",
                    "returncode": None,
                }
                failures.append(failure)
                if failure_out is not None:
                    failure_out.write(json.dumps(failure, sort_keys=True) + "\n")
                    failure_out.flush()
                continue
            if proc.returncode != 0:
                failure = {
                    "matrix": str(matrix),
                    "reason": proc.stderr.strip(),
                    "returncode": proc.returncode,
                }
                failures.append(failure)
                if failure_out is not None:
                    failure_out.write(json.dumps(failure, sort_keys=True) + "\n")
                    failure_out.flush()
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
        if failure_out is not None:
            failure_out.close()

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
        for failure in failures:
            print(f"  {failure['matrix']}: {failure['reason']}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
