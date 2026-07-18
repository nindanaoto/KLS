#!/usr/bin/env python3
"""Run the flag-style KLU comparison benchmark over a matrix corpus."""

from __future__ import annotations

import argparse
import json
import math
import pathlib
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
        line = raw.split("#", 1)[0].strip()
        if line:
            names.append(line.lower())
    return names


def select_matrices(
    matrix_dir: pathlib.Path, manifest: pathlib.Path | None
) -> list[pathlib.Path]:
    matrices = sorted(matrix_dir.rglob("*.mtx"))
    if manifest is None:
        return matrices

    matrix_by_stem: dict[str, pathlib.Path] = {}
    duplicate_stems: list[str] = []
    for matrix in matrices:
        stem = matrix.stem.lower()
        if stem in matrix_by_stem:
            duplicate_stems.append(stem)
        else:
            matrix_by_stem[stem] = matrix
    if duplicate_stems:
        raise ValueError(
            "duplicate matrix basenames under "
            f"{matrix_dir}: {', '.join(sorted(set(duplicate_stems)))}"
        )

    selected: list[pathlib.Path] = []
    missing: list[str] = []
    for name in read_manifest(manifest):
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--klu-compare", type=pathlib.Path, required=True)
    parser.add_argument("--matrix-dir", type=pathlib.Path, required=True)
    parser.add_argument("--manifest", type=pathlib.Path)
    parser.add_argument("--jsonl", type=pathlib.Path)
    parser.add_argument("--repeat", type=int, default=1)
    parser.add_argument("--factor-repeat", type=int, default=None)
    parser.add_argument("--refactor-repeat", type=int, default=5)
    parser.add_argument(
        "--refactor-values",
        choices=("unchanged", "rank-preserving", "entrywise"),
        default="unchanged",
    )
    parser.add_argument("--refactor-value-amplitude", type=float, default=1.0e-3)
    parser.add_argument("--passes", type=int, default=1)
    parser.add_argument("--timeout", type=float)
    parser.add_argument("--skip", type=int, default=0)
    parser.add_argument("--limit", type=int)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    factor_repeat = args.repeat if args.factor_repeat is None else args.factor_repeat
    if (
        args.repeat <= 0
        or factor_repeat < 0
        or args.refactor_repeat < 0
        or args.passes <= 0
        or args.skip < 0
        or (args.limit is not None and args.limit < 0)
        or not math.isfinite(args.refactor_value_amplitude)
        or args.refactor_value_amplitude < 0.0
        or args.refactor_value_amplitude >= 1.0
        or (
            args.refactor_values != "unchanged"
            and args.refactor_value_amplitude <= 0.0
        )
    ):
        print("repeat, selection, or refactor-value arguments are invalid", file=sys.stderr)
        return 1

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

    out = None
    failure_out = None
    if args.jsonl is not None:
        args.jsonl.parent.mkdir(parents=True, exist_ok=True)
        out = args.jsonl.open("w", encoding="utf-8")
        failure_out = args.jsonl.with_suffix(".failures").open("w", encoding="utf-8")

    rows: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []
    try:
        for matrix in matrices:
            samples: list[dict[str, object]] = []
            sample_failures: list[str] = []
            for _ in range(args.passes):
                cmd = [
                    str(args.klu_compare),
                    str(matrix),
                    "--repeat",
                    str(args.repeat),
                    "--factor-repeat",
                    str(factor_repeat),
                    "--refactor-repeat",
                    str(args.refactor_repeat),
                    "--refactor-values",
                    args.refactor_values,
                    "--refactor-value-amplitude",
                    f"{args.refactor_value_amplitude:.17g}",
                    "--json",
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
                    sample_failures.append(f"timeout after {args.timeout:g}s")
                    continue
                if proc.returncode != 0:
                    sample_failures.append(proc.stderr.strip())
                    continue
                try:
                    sample = json.loads(proc.stdout)
                    float(sample["spice_cycle_seconds"])
                except (json.JSONDecodeError, KeyError, TypeError, ValueError) as exc:
                    sample_failures.append(f"invalid JSON result: {exc}")
                    continue
                samples.append(sample)

            if not samples:
                failure = {
                    "matrix": str(matrix),
                    "reason": "; ".join(sample_failures),
                    "sample_failures": sample_failures,
                }
                failures.append(failure)
                if failure_out is not None:
                    failure_out.write(json.dumps(failure, sort_keys=True) + "\n")
                    failure_out.flush()
                continue

            samples.sort(key=lambda row: float(row["spice_cycle_seconds"]))
            row = dict(samples[len(samples) // 2])
            row["spice_cycle_seconds_samples"] = [
                float(sample["spice_cycle_seconds"]) for sample in samples
            ]
            row["spice_cycle_seconds_median"] = float(row["spice_cycle_seconds"])
            row["passes_ok"] = len(samples)
            row["passes_failed"] = len(sample_failures)
            rows.append(row)
            if out is not None:
                out.write(json.dumps(row, sort_keys=True) + "\n")
                out.flush()
            print(
                f"{matrix}: {float(row['spice_cycle_seconds']):.6g}s "
                f"median of {len(samples)} pass(es)"
            )
    finally:
        if out is not None:
            out.close()
        if failure_out is not None:
            failure_out.close()

    summary = {
        "matrices_total": len(matrices),
        "matrices_ok": len(rows),
        "matrices_failed": len(failures),
        "spice_cycle_geomean_seconds": geometric_mean(
            [float(row["spice_cycle_seconds"]) for row in rows]
        ),
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
