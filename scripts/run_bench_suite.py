#!/usr/bin/env python3
"""Run kls_bench over a directory of MatrixMarket files and summarize results."""

from __future__ import annotations

import argparse
import json
import math
import os
import pathlib
import subprocess
import sys

from summarize_row_pipeline_trace import summarize as summarize_row_pipeline_trace


def spice_cycle_seconds(row: dict[str, object]) -> float:
    analysis = float(row["analysis_seconds"])
    factor = float(row.get("initial_factor_seconds", row["factor_seconds_avg"]))
    solve = float(row["solve_seconds_avg"])
    first = float(row.get("refactor_first_seconds", row["refactor_seconds_avg"]))
    steady = float(
        row.get("refactor_steady_seconds_avg", row["refactor_seconds_avg"])
    )
    refactor_solve_first = float(
        row.get("refactor_solve_first_seconds", solve)
    )
    refactor_solve_steady = float(
        row.get("refactor_solve_steady_seconds_avg", solve)
    )
    return (
        analysis
        + factor
        + solve
        + first
        + refactor_solve_first
        + 98.0 * (steady + refactor_solve_steady)
    )


def geometric_mean(values: list[float]) -> float:
    positives = [v for v in values if v > 0.0 and math.isfinite(v)]
    if not positives:
        return math.nan
    return math.exp(sum(math.log(v) for v in positives) / len(positives))


def clipped_text(text: str, limit: int = 4000) -> str:
    if len(text) <= limit:
        return text
    return text[:limit] + "...<truncated>"


def completed_text(value: str | bytes | None) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return value


def append_solver_options(cmd: list[str], args: argparse.Namespace) -> None:
    cmd.extend(
        [
            "--threads",
            str(args.threads),
            "--backend",
            args.backend,
            "--ordering",
            args.ordering,
            "--orientation",
            args.orientation,
            "--scale",
            args.scale,
            "--input-index",
            args.input_index,
            "--refactor-values",
            args.refactor_values,
            "--refactor-value-amplitude",
            str(args.refactor_value_amplitude),
        ]
    )
    if args.no_btf:
        cmd.append("--no-btf")
    if args.no_fast_factor:
        cmd.append("--no-fast-factor")
    if args.no_static_pivoting:
        cmd.append("--no-static-pivoting")
    if args.no_transpose_solve:
        cmd.append("--no-transpose-solve")
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


def collect_analyze_failure_diagnostic(
    args: argparse.Namespace, matrix: pathlib.Path
) -> dict[str, object] | None:
    cmd = [
        str(args.kls_bench),
        str(matrix),
        "--analyze-only",
        "--json",
    ]
    if args.factor_repeat is not None:
        cmd.extend(["--factor-repeat", str(args.factor_repeat)])
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
        diagnostic["stderr"] = clipped_text(completed_text(proc.stderr).strip())
        return diagnostic
    try:
        diagnostic["row"] = json.loads(proc.stdout)
        diagnostic["status"] = "ok"
    except json.JSONDecodeError:
        diagnostic["status"] = "invalid-json"
        diagnostic["stdout"] = clipped_text(completed_text(proc.stdout).strip())
        diagnostic["stderr"] = clipped_text(completed_text(proc.stderr).strip())
    return diagnostic


def failure_trace_path(args: argparse.Namespace, matrix: pathlib.Path) -> pathlib.Path:
    trace_dir = args.failure_trace_dir
    if trace_dir is None:
        if args.jsonl is not None:
            trace_dir = args.jsonl.parent / f"{args.jsonl.stem}.failure-traces"
        else:
            trace_dir = pathlib.Path("build") / "failure-traces"
    trace_dir.mkdir(parents=True, exist_ok=True)
    return trace_dir / f"{matrix.stem}.stderr"


def collect_trace_failure_diagnostic(
    args: argparse.Namespace, matrix: pathlib.Path
) -> dict[str, object]:
    trace_path = failure_trace_path(args, matrix)
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
    if args.failure_trace_kls_first_factor != "env":
        cmd.extend(["--kls-first-factor", args.failure_trace_kls_first_factor])

    env = os.environ.copy()
    env["KLS_TRACE_KLS_FIRST_FACTOR"] = "1"
    env["KLS_TRACE_ROW_PIPELINE"] = "1"
    env["KLS_TRACE_ROW_PIPELINE_LONG_ROW_ENTRIES"] = str(
        args.failure_trace_long_row_entries
    )
    env_overrides = {
        "KLS_TRACE_KLS_FIRST_FACTOR": env["KLS_TRACE_KLS_FIRST_FACTOR"],
        "KLS_TRACE_ROW_PIPELINE": env["KLS_TRACE_ROW_PIPELINE"],
        "KLS_TRACE_ROW_PIPELINE_LONG_ROW_ENTRIES":
            env["KLS_TRACE_ROW_PIPELINE_LONG_ROW_ENTRIES"],
    }
    try:
        proc = subprocess.run(
            cmd,
            text=True,
            capture_output=True,
            check=False,
            timeout=args.failure_trace_timeout,
            env=env,
        )
        stdout = completed_text(proc.stdout)
        stderr = completed_text(proc.stderr)
        status = "ok" if proc.returncode == 0 else "failed"
        returncode: int | None = proc.returncode
    except subprocess.TimeoutExpired as exc:
        stdout = completed_text(exc.stdout)
        stderr = completed_text(exc.stderr)
        status = "timeout"
        returncode = None

    trace_path.write_text(stderr, encoding="utf-8")
    diagnostic: dict[str, object] = {
        "kind": "row-pipeline-trace",
        "status": status,
        "timeout_seconds": args.failure_trace_timeout,
        "returncode": returncode,
        "stderr_path": str(trace_path),
        "command": cmd,
        "env": env_overrides,
    }
    if stdout.strip():
        try:
            diagnostic["row"] = json.loads(stdout)
        except json.JSONDecodeError:
            diagnostic["stdout"] = clipped_text(stdout.strip())
    if status == "failed":
        diagnostic["stderr"] = clipped_text(stderr.strip())
    try:
        diagnostic["summary"] = summarize_row_pipeline_trace(trace_path)
    except OSError as exc:
        diagnostic["summary_error"] = str(exc)
    return diagnostic


def collect_failure_diagnostic(
    args: argparse.Namespace, matrix: pathlib.Path
) -> dict[str, object] | None:
    if args.failure_diagnostics == "none":
        return None
    if args.failure_diagnostics == "analyze":
        return collect_analyze_failure_diagnostic(args, matrix)
    if args.failure_diagnostics == "trace":
        return collect_trace_failure_diagnostic(args, matrix)
    return {
        "kind": "combined",
        "analyze": collect_analyze_failure_diagnostic(args, matrix),
        "trace": collect_trace_failure_diagnostic(args, matrix),
    }


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
    parser.add_argument(
        "--analyze-only",
        action="store_true",
        help="run only symbolic analysis and select the median analysis time",
    )
    parser.add_argument(
        "--factor-repeat",
        type=int,
        default=None,
        help="extra full factors after the initial factor (default: bench follows --repeat)",
    )
    parser.add_argument("--refactor-repeat", type=int, default=5)
    parser.add_argument(
        "--refactor-values",
        choices=[
            "unchanged",
            "rank-preserving",
            "entrywise",
            "localized-entrywise",
        ],
        default="unchanged",
    )
    parser.add_argument("--refactor-value-amplitude", type=float, default=0.001)
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
    parser.add_argument(
        "--backend", choices=["auto", "kls", "serial"], default="auto"
    )
    parser.add_argument("--timeout", type=float)
    parser.add_argument(
        "--ordering",
        choices=["auto", "amd", "colamd", "natural", "metis", "scotch", "amf"],
        default="auto",
    )
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
    parser.add_argument("--no-transpose-solve", action="store_true")
    parser.add_argument("--require-metis", action="store_true")
    parser.add_argument("--require-scotch", action="store_true")
    parser.add_argument("--require-spral-scaling", action="store_true")
    parser.add_argument(
        "--failure-diagnostics",
        choices=["none", "analyze", "trace", "all"],
        default="analyze",
        help=(
            "attach diagnostics when every sample fails: analyze-only, "
            "row-pipeline trace, or both"
        ),
    )
    parser.add_argument(
        "--failure-diagnostic-timeout",
        type=float,
        default=30.0,
        help="timeout in seconds for failure diagnostics",
    )
    parser.add_argument(
        "--failure-trace-timeout",
        type=float,
        default=45.0,
        help="timeout in seconds for row-pipeline failure trace diagnostics",
    )
    parser.add_argument(
        "--failure-trace-long-row-entries",
        type=int,
        default=100000,
        help=(
            "KLS_TRACE_ROW_PIPELINE_LONG_ROW_ENTRIES threshold used by "
            "row-pipeline failure trace diagnostics"
        ),
    )
    parser.add_argument(
        "--failure-trace-dir",
        type=pathlib.Path,
        help=(
            "directory for row-pipeline failure trace stderr files; defaults "
            "to a sibling directory named after --jsonl, or build/failure-traces"
        ),
    )
    parser.add_argument(
        "--failure-trace-kls-first-factor",
        choices=["env", "off", "on"],
        default="env",
        help=(
            "trace-only override for --kls-first-factor; use 'on' when a "
            "timeout needs KLS-owned row-pipeline evidence"
        ),
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
    if args.repeat <= 0 or args.refactor_repeat < 0 or (
        args.factor_repeat is not None and args.factor_repeat < 0
    ):
        print("repeat arguments are invalid", file=sys.stderr)
        return 1
    if not 0.0 <= args.refactor_value_amplitude < 1.0 or (
        args.refactor_values != "unchanged"
        and args.refactor_value_amplitude <= 0.0
    ):
        print("--refactor-value-amplitude must be in (0, 1) for changed values", file=sys.stderr)
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
    if args.failure_trace_timeout is not None and args.failure_trace_timeout <= 0:
        print("--failure-trace-timeout must be positive", file=sys.stderr)
        return 1
    if args.failure_trace_long_row_entries < 0:
        print("--failure-trace-long-row-entries must be non-negative", file=sys.stderr)
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
                cmd = [str(args.kls_bench), str(matrix), "--json"]
                if args.analyze_only:
                    cmd.append("--analyze-only")
                else:
                    cmd.extend(
                        [
                            "--repeat",
                            str(args.repeat),
                            "--refactor-repeat",
                            str(args.refactor_repeat),
                        ]
                    )
                    if args.factor_repeat is not None:
                        cmd.extend(["--factor-repeat", str(args.factor_repeat)])
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
                if args.analyze_only:
                    sample["selection_seconds"] = float(sample["analysis_seconds"])
                else:
                    sample["spice_cycle_seconds"] = spice_cycle_seconds(sample)
                    sample["selection_seconds"] = sample["spice_cycle_seconds"]
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
            samples.sort(key=lambda row: float(row["selection_seconds"]))
            row = dict(samples[len(samples) // 2])
            selection_samples = [
                float(sample["selection_seconds"]) for sample in samples
            ]
            row["selection_seconds_samples"] = selection_samples
            row["selection_seconds_median"] = row["selection_seconds"]
            if not args.analyze_only:
                row["spice_cycle_seconds_samples"] = selection_samples
                row["spice_cycle_seconds_median"] = row["spice_cycle_seconds"]
            row["passes_ok"] = len(samples)
            row["passes_failed"] = len(sample_failures)
            rows.append(row)
            if out is not None:
                out.write(json.dumps(row, sort_keys=True) + "\n")
                out.flush()
            if args.passes == 1:
                print(f"{matrix}: {row['selection_seconds']:.6g}s")
            else:
                print(
                    f"{matrix}: {row['selection_seconds']:.6g}s "
                    f"median of {len(samples)} pass(es)"
                )
    finally:
        if out is not None:
            out.close()
        if failure_out is not None:
            failure_out.close()

    selection_values = [float(row["selection_seconds"]) for row in rows]
    summary = {
        "matrices_total": len(matrices),
        "matrices_ok": len(rows),
        "matrices_failed": len(failures),
        (
            "analysis_geomean_seconds"
            if args.analyze_only
            else "spice_cycle_geomean_seconds"
        ): geometric_mean(selection_values),
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
