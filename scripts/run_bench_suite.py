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


def clipped_text(text: str, limit: int = 4000) -> str:
    if len(text) <= limit:
        return text
    return text[:limit] + "...<truncated>"


def append_solver_options(cmd: list[str], args: argparse.Namespace) -> None:
    cmd.extend(
        [
            "--threads",
            str(args.threads),
            "--ordering",
            args.ordering,
            "--orientation",
            args.orientation,
            "--scale",
            args.scale,
            "--input-index",
            args.input_index,
        ]
    )
    if args.no_btf:
        cmd.append("--no-btf")
    if args.no_fast_factor:
        cmd.append("--no-fast-factor")
    if args.no_static_pivoting:
        cmd.append("--no-static-pivoting")
    if args.pivot_tol is not None:
        cmd.extend(["--pivot-tol", str(args.pivot_tol)])
    if args.row_refactor != "env":
        cmd.extend(["--row-refactor", args.row_refactor])
    if args.kls_first_factor != "env":
        cmd.extend(["--kls-first-factor", args.kls_first_factor])
    if args.row_solve != "env":
        cmd.extend(["--row-solve", args.row_solve])
    if args.stress_diagonal_scale is not None:
        cmd.extend(["--stress-diagonal-scale", str(args.stress_diagonal_scale)])
    if args.stress_diagonal_column is not None:
        cmd.extend(["--stress-diagonal-column", str(args.stress_diagonal_column)])


def collect_failure_diagnostic(
    args: argparse.Namespace, matrix: pathlib.Path
) -> dict[str, object] | None:
    if args.failure_diagnostics == "none":
        return None
    cmd = [
        str(args.kls_bench),
        str(matrix),
        "--analyze-only",
        "--json",
    ]
    append_solver_options(cmd, args)
    try:
        proc = subprocess.run(
            cmd,
            text=True,
            capture_output=True,
            check=False,
            timeout=args.failure_diagnostic_timeout,
        )
    except subprocess.TimeoutExpired:
        timeout = args.failure_diagnostic_timeout
        return {
            "kind": "analyze",
            "status": "timeout",
            "timeout_seconds": timeout,
        }
    diagnostic: dict[str, object] = {
        "kind": "analyze",
        "returncode": proc.returncode,
    }
    if proc.returncode != 0:
        diagnostic["status"] = "failed"
        diagnostic["stderr"] = clipped_text(proc.stderr.strip())
        return diagnostic
    try:
        diagnostic["row"] = json.loads(proc.stdout)
        diagnostic["status"] = "ok"
    except json.JSONDecodeError:
        diagnostic["status"] = "invalid-json"
        diagnostic["stdout"] = clipped_text(proc.stdout.strip())
        diagnostic["stderr"] = clipped_text(proc.stderr.strip())
    return diagnostic


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
    parser.add_argument("--kls-bench", type=pathlib.Path, required=True)
    parser.add_argument("--matrix-dir", type=pathlib.Path, required=True)
    parser.add_argument("--manifest", type=pathlib.Path)
    parser.add_argument("--jsonl", type=pathlib.Path)
    parser.add_argument("--repeat", type=int, default=5)
    parser.add_argument("--refactor-repeat", type=int, default=5)
    parser.add_argument("--passes", type=int, default=1)
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
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--timeout", type=float)
    parser.add_argument("--ordering", choices=["auto", "amd", "colamd", "natural", "metis", "scotch"], default="auto")
    parser.add_argument("--orientation", choices=["auto", "normal", "transpose"], default="auto")
    parser.add_argument("--scale", choices=["auto", "-1", "0", "1", "2"], default="auto")
    parser.add_argument("--input-index", choices=["auto", "32", "64"], default="auto")
    parser.add_argument("--pivot-tol", type=float, default=None)
    parser.add_argument(
        "--row-refactor",
        choices=["env", "off", "refactor", "checked", "all"],
        default="env",
    )
    parser.add_argument(
        "--kls-first-factor",
        choices=["env", "off", "on"],
        default="env",
    )
    parser.add_argument(
        "--row-solve",
        choices=["env", "off", "on"],
        default="env",
    )
    parser.add_argument("--stress-diagonal-scale", type=float, default=None)
    parser.add_argument("--stress-diagonal-column", type=int, default=None)
    parser.add_argument("--no-btf", action="store_true")
    parser.add_argument("--no-fast-factor", action="store_true")
    parser.add_argument("--no-static-pivoting", action="store_true")
    parser.add_argument("--require-metis", action="store_true")
    parser.add_argument("--require-scotch", action="store_true")
    parser.add_argument("--require-spral-scaling", action="store_true")
    parser.add_argument(
        "--failure-diagnostics",
        choices=["none", "analyze"],
        default="analyze",
        help="attach an analyze-only diagnostic row when every sample fails",
    )
    parser.add_argument(
        "--failure-diagnostic-timeout",
        type=float,
        default=30.0,
        help="timeout in seconds for failure diagnostics",
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
    if args.passes <= 0:
        print("--passes must be positive", file=sys.stderr)
        return 1
    if args.skip < 0:
        print("--skip must be non-negative", file=sys.stderr)
        return 1
    if args.limit is not None and args.limit < 0:
        print("--limit must be non-negative", file=sys.stderr)
        return 1
    if args.failure_diagnostic_timeout is not None and args.failure_diagnostic_timeout <= 0:
        print("--failure-diagnostic-timeout must be positive", file=sys.stderr)
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
    if args.jsonl is not None:
        args.jsonl.parent.mkdir(parents=True, exist_ok=True)
        out = args.jsonl.open("w", encoding="utf-8")
        failure_out = args.jsonl.with_suffix(".failures").open("w", encoding="utf-8")
    else:
        out = None
        failure_out = None

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
                    "--json",
                ]
                append_solver_options(cmd, args)
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
                missing_features: list[str] = []
                if args.require_metis and not sample.get("build_has_metis"):
                    missing_features.append("METIS")
                if args.require_scotch and not sample.get("build_has_scotch"):
                    missing_features.append("SCOTCH")
                if args.require_spral_scaling and not sample.get("build_has_spral_scaling"):
                    missing_features.append("SPRAL scaling")
                if missing_features:
                    sample_failures.append(
                        "benchmark binary missing required feature(s): "
                        + ", ".join(missing_features)
                    )
                    continue
                sample["spice_cycle_seconds"] = spice_cycle_seconds(sample)
                samples.append(sample)
            if not samples:
                failure = {
                    "matrix": str(matrix),
                    "reason": "; ".join(sample_failures),
                    "sample_failures": sample_failures,
                }
                diagnostic = collect_failure_diagnostic(args, matrix)
                if diagnostic is not None:
                    failure["diagnostic"] = diagnostic
                failures.append(failure)
                if failure_out is not None:
                    failure_out.write(json.dumps(failure, sort_keys=True) + "\n")
                    failure_out.flush()
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
