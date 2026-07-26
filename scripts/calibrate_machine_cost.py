#!/usr/bin/env python3
"""Build and evaluate a persistent, diagnostic KLS machine-cost profile.

The profile contains machine/build-specific rates and uncertainty, never
matrix names or per-matrix decisions.  Collection deliberately uses explicit
ordering/BTF variants so AUTO behavior is not part of the fitted labels.
Nothing in this tool changes solver routing; it is an offline experiment for
deciding whether a calibrated relative-cost model is accurate enough to do so.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
import pathlib
import platform
import re
import shutil
import statistics
import subprocess
import sys
import time
from collections import defaultdict
from typing import Iterable


SAMPLE_SCHEMA = "kls-machine-cost-samples-v1"
PROFILE_SCHEMA = "kls-machine-cost-profile-v1"
REPORT_SCHEMA = "kls-machine-cost-report-v1"
PREDICTION_MODES = ("symbolic-only", "baseline-measured")

VARIANTS: dict[str, dict[str, object]] = {
    "amd_btf": {"ordering": "amd", "btf": True},
    "colamd_btf": {"ordering": "colamd", "btf": True},
    "amf_btf": {"ordering": "amf", "btf": True},
    "amf_no_btf": {"ordering": "amf", "btf": False},
    "metis_btf": {"ordering": "metis", "btf": True},
    "scotch_btf": {"ordering": "scotch", "btf": True},
}
DEFAULT_VARIANTS = tuple(VARIANTS)

NUMERIC_SAMPLE_FIELDS = (
    "n",
    "nnz",
    "analysis_seconds",
    "initial_factor_seconds",
    "factor_seconds_avg",
    "refactor_first_seconds",
    "refactor_steady_seconds_avg",
    "refactor_solve_first_seconds",
    "refactor_solve_steady_seconds_avg",
    "solve_seconds_avg",
    "relative_residual_l2",
    "refactor_max_relative_residual",
    "nblocks",
    "max_block",
    "structural_rank",
    "numerical_rank",
    "offdiag_pivots",
    "nnz_l",
    "nnz_u",
    "estimated_flops",
    "factor_flops",
    "memory_bytes",
    "memory_peak_bytes",
)

TEXT_SAMPLE_FIELDS = (
    "ordering",
    "orientation",
    "requested_scale",
    "scale",
    "requested_btf",
    "btf",
    "initial_factor_path",
    "last_refactor_path",
    "build_has_metis",
    "build_has_scotch",
    "build_has_spral_scaling",
    "build_has_cblas",
)

SYMBOLIC_NUMERIC_FIELDS = tuple(
    f"symbolic_{key}"
    for key in (
        "analysis_seconds",
        "nblocks",
        "max_block",
        "structural_rank",
        "nnz_l",
        "nnz_u",
        "estimated_flops",
    )
)

SYMBOLIC_TEXT_FIELDS = (
    "symbolic_ordering",
    "symbolic_orientation",
    "symbolic_scale",
    "symbolic_btf",
)

STAGES: dict[str, tuple[str, tuple[str, ...]]] = {
    "analysis": ("analysis_seconds", ("nnz_log2_n",)),
    "factor": ("initial_factor_seconds", ("estimated_flops", "factor_entries")),
    "refactor": (
        "refactor_steady_seconds_avg",
        ("estimated_flops", "factor_entries"),
    ),
    "solve": ("refactor_solve_steady_seconds_avg", ("factor_entries",)),
}

FINGERPRINT_ENV = (
    "OMP_NUM_THREADS",
    "OMP_PROC_BIND",
    "OMP_PLACES",
    "OMP_DYNAMIC",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
)


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha256_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def cpu_identity() -> dict[str, str]:
    wanted = {
        "vendor_id",
        "cpu family",
        "model",
        "model name",
        "stepping",
        "microcode",
    }
    identity: dict[str, str] = {}
    cpuinfo = pathlib.Path("/proc/cpuinfo")
    try:
        for raw in cpuinfo.read_text(encoding="utf-8", errors="replace").splitlines():
            if not raw.strip():
                break
            key, separator, value = raw.partition(":")
            normalized = key.strip().lower()
            if separator and normalized in wanted:
                identity[normalized.replace(" ", "_")] = value.strip()
    except OSError:
        pass
    if not identity:
        identity["processor"] = platform.processor() or "unknown"
    return identity


def machine_fingerprint(
    kls_bench: pathlib.Path, threads: int, cpu_list: str | None
) -> dict[str, object]:
    try:
        affinity = sorted(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        affinity = []
    controlled_env = {name: os.environ.get(name) for name in FINGERPRINT_ENV}
    controlled_env["OPENBLAS_NUM_THREADS"] = "1"
    return {
        "kls_bench_sha256": sha256_file(kls_bench),
        "machine": platform.machine(),
        "kernel": platform.release(),
        "cpu": cpu_identity(),
        "logical_cpu_count": os.cpu_count(),
        "parent_affinity": affinity,
        "requested_cpu_list": cpu_list,
        "threads": threads,
        "environment": controlled_env,
    }


def parse_variants(raw: str) -> list[str]:
    variants = [item.strip() for item in raw.split(",") if item.strip()]
    if not variants:
        raise ValueError("at least one variant is required")
    unknown = sorted(set(variants) - VARIANTS.keys())
    if unknown:
        raise ValueError("unknown variant(s): " + ", ".join(unknown))
    if len(set(variants)) != len(variants):
        raise ValueError("variant list contains duplicates")
    return variants


def read_manifest(path: pathlib.Path) -> list[str]:
    names: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        if line:
            names.append(line.lower())
    return names


def select_matrices(
    matrix_dir: pathlib.Path, manifests: list[pathlib.Path]
) -> list[pathlib.Path]:
    by_stem: dict[str, pathlib.Path] = {}
    duplicates: set[str] = set()
    for matrix in sorted(matrix_dir.rglob("*.mtx")):
        stem = matrix.stem.lower()
        if stem in by_stem:
            duplicates.add(stem)
        else:
            by_stem[stem] = matrix
    if duplicates:
        raise ValueError(
            "duplicate matrix basenames under matrix directory: "
            + ", ".join(sorted(duplicates))
        )
    wanted: list[str] = []
    seen: set[str] = set()
    for manifest in manifests:
        for name in read_manifest(manifest):
            if name not in seen:
                wanted.append(name)
                seen.add(name)
    missing = [name for name in wanted if name not in by_stem]
    if missing:
        raise ValueError("manifest matrices not downloaded: " + ", ".join(missing))
    return [by_stem[name] for name in wanted]


def clipped(text: str | bytes | None, limit: int = 2000) -> str:
    if text is None:
        return ""
    if isinstance(text, bytes):
        text = text.decode("utf-8", errors="replace")
    stripped = text.strip()
    return stripped if len(stripped) <= limit else stripped[:limit] + "...<truncated>"


def benchmark_command(
    args: argparse.Namespace,
    matrix: pathlib.Path,
    variant: str,
    analyze_only: bool = False,
) -> list[str]:
    definition = VARIANTS[variant]
    cmd = [
        str(args.kls_bench),
        str(matrix),
        "--json",
        "--threads",
        str(args.threads),
        "--backend",
        "auto",
        "--ordering",
        str(definition["ordering"]),
        "--orientation",
        "normal",
        "--scale",
        args.scale,
        "--input-index",
        "auto",
    ]
    if analyze_only:
        cmd.append("--analyze-only")
    else:
        cmd.extend(
            [
                "--repeat",
                "1",
                "--factor-repeat",
                "0",
                "--refactor-repeat",
                str(args.refactor_repeat),
                "--refactor-values",
                "entrywise",
                "--refactor-value-amplitude",
                str(args.refactor_value_amplitude),
                "--no-transpose-solve",
            ]
        )
    if not bool(definition["btf"]):
        cmd.append("--no-btf")
    if args.cpu_list is not None:
        cmd = ["taskset", "-c", args.cpu_list, *cmd]
    return cmd


def run_benchmark(
    cmd: list[str], timeout: float, child_env: dict[str, str]
) -> tuple[str, dict[str, object] | None, float, str, str, str]:
    start = time.monotonic()
    try:
        proc = subprocess.run(
            cmd,
            text=True,
            capture_output=True,
            check=False,
            timeout=timeout,
            env=child_env,
        )
    except subprocess.TimeoutExpired as exc:
        return (
            "timeout",
            None,
            time.monotonic() - start,
            clipped(exc.stdout),
            clipped(exc.stderr),
            f"timeout after {timeout:g}s",
        )
    if proc.returncode != 0:
        return (
            "failed",
            None,
            time.monotonic() - start,
            clipped(proc.stdout),
            clipped(proc.stderr),
            f"exit status {proc.returncode}",
        )
    try:
        parsed = json.loads(proc.stdout)
    except json.JSONDecodeError as exc:
        return (
            "failed",
            None,
            time.monotonic() - start,
            clipped(proc.stdout),
            clipped(proc.stderr),
            f"invalid benchmark JSON: {exc}",
        )
    return "ok", parsed, time.monotonic() - start, "", clipped(proc.stderr), ""


def compact_sample(
    row: dict[str, object], matrix: pathlib.Path, variant: str, pass_index: int
) -> dict[str, object]:
    sample: dict[str, object] = {
        "record_type": "sample",
        "status": "ok",
        "matrix": str(matrix),
        "matrix_name": matrix.stem,
        "variant": variant,
        "pass": pass_index,
    }
    for key in (*NUMERIC_SAMPLE_FIELDS, *TEXT_SAMPLE_FIELDS):
        if key in row:
            sample[key] = row[key]
    return sample


def attach_symbolic_sample(
    sample: dict[str, object], symbolic: dict[str, object]
) -> None:
    for key in (
        "analysis_seconds",
        "nblocks",
        "max_block",
        "structural_rank",
        "nnz_l",
        "nnz_u",
        "estimated_flops",
        "ordering",
        "orientation",
        "scale",
        "btf",
    ):
        if key in symbolic:
            sample[f"symbolic_{key}"] = symbolic[key]


def collect(args: argparse.Namespace) -> int:
    if args.threads <= 0 or args.passes <= 0 or args.refactor_repeat < 2:
        raise ValueError(
            "threads/passes must be positive and refactor-repeat at least 2"
        )
    if args.timeout <= 0.0:
        raise ValueError("timeout must be positive")
    if not 0.0 < args.refactor_value_amplitude < 1.0:
        raise ValueError("refactor value amplitude must be in (0, 1)")
    if args.cpu_list is not None:
        if shutil.which("taskset") is None:
            raise ValueError("--cpu-list requires taskset")
        if re.fullmatch(r"[0-9,-]+", args.cpu_list) is None:
            raise ValueError(
                "--cpu-list must contain only CPU numbers, commas, and dashes"
            )
    variants = parse_variants(args.variants)
    matrices = select_matrices(args.matrix_dir, args.manifest)
    if args.skip < 0 or args.limit is not None and args.limit < 0:
        raise ValueError("skip and limit must be non-negative")
    matrices = matrices[args.skip :]
    if args.limit is not None:
        matrices = matrices[: args.limit]
    if not matrices:
        raise ValueError("matrix selection is empty")

    args.samples.parent.mkdir(parents=True, exist_ok=True)
    metadata = {
        "record_type": "metadata",
        "schema": SAMPLE_SCHEMA,
        "created_utc": utc_now(),
        "fingerprint": machine_fingerprint(args.kls_bench, args.threads, args.cpu_list),
        "manifests": [
            {"path": str(path), "sha256": sha256_file(path)} for path in args.manifest
        ],
        "settings": {
            "passes": args.passes,
            "timeout_seconds": args.timeout,
            "refactor_repeat": args.refactor_repeat,
            "refactor_value_amplitude": args.refactor_value_amplitude,
            "symbolic_probe": True,
            "scale": args.scale,
            "variants": variants,
        },
    }
    ok = 0
    failed = 0
    total = len(matrices) * len(variants) * args.passes
    completed = 0
    child_env = os.environ.copy()
    child_env["OPENBLAS_NUM_THREADS"] = "1"
    child_env["KLS_BENCH_VERIFY_EACH_REFACTOR"] = "1"
    with args.samples.open("w", encoding="utf-8") as output:
        output.write(json.dumps(metadata, sort_keys=True) + "\n")
        output.flush()
        for matrix in matrices:
            for variant in variants:
                for pass_index in range(args.passes):
                    symbolic_cmd = benchmark_command(
                        args, matrix, variant, analyze_only=True
                    )
                    (
                        symbolic_status,
                        symbolic,
                        symbolic_wall,
                        symbolic_stdout,
                        symbolic_stderr,
                        symbolic_reason,
                    ) = run_benchmark(symbolic_cmd, args.timeout, child_env)
                    if symbolic_status != "ok" or symbolic is None:
                        record = {
                            "record_type": "sample",
                            "status": symbolic_status,
                            "phase": "symbolic",
                            "matrix": str(matrix),
                            "matrix_name": matrix.stem,
                            "variant": variant,
                            "pass": pass_index,
                            "wall_seconds": symbolic_wall,
                            "reason": symbolic_reason,
                            "stdout": symbolic_stdout,
                            "stderr": symbolic_stderr,
                        }
                        failed += 1
                    else:
                        cmd = benchmark_command(args, matrix, variant)
                        status, parsed, wall, stdout, stderr, reason = run_benchmark(
                            cmd, args.timeout, child_env
                        )
                        if status != "ok" or parsed is None:
                            record = {
                                "record_type": "sample",
                                "status": status,
                                "phase": "numeric",
                                "matrix": str(matrix),
                                "matrix_name": matrix.stem,
                                "variant": variant,
                                "pass": pass_index,
                                "wall_seconds": wall,
                                "symbolic_wall_seconds": symbolic_wall,
                                "reason": reason,
                                "stdout": stdout,
                                "stderr": stderr,
                            }
                            attach_symbolic_sample(record, symbolic)
                            failed += 1
                        else:
                            record = compact_sample(parsed, matrix, variant, pass_index)
                            record["wall_seconds"] = wall
                            record["symbolic_wall_seconds"] = symbolic_wall
                            attach_symbolic_sample(record, symbolic)
                            ok += 1
                    output.write(json.dumps(record, sort_keys=True) + "\n")
                    output.flush()
                    completed += 1
                    print(
                        f"[{completed}/{total}] {matrix.stem} {variant}: "
                        f"{record['status']}",
                        flush=True,
                    )
                    if record["status"] != "ok":
                        # Repeating a deterministic setup failure or a full
                        # timeout cannot improve a median and can multiply an
                        # expensive negative by --passes.
                        break
    print(
        json.dumps(
            {
                "attempted": completed,
                "planned": total,
                "ok": ok,
                "failed": failed,
            },
            sort_keys=True,
        )
    )
    return 0 if ok else 2


def load_sample_files(
    paths: Iterable[pathlib.Path],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    metadata: list[dict[str, object]] = []
    records: list[dict[str, object]] = []
    for path in paths:
        seen_metadata = False
        with path.open(encoding="utf-8") as stream:
            for line_number, line in enumerate(stream, 1):
                if not line.strip():
                    continue
                try:
                    row = json.loads(line)
                except json.JSONDecodeError as exc:
                    raise ValueError(
                        f"{path}:{line_number}: invalid JSON: {exc}"
                    ) from exc
                if row.get("record_type") == "metadata":
                    if seen_metadata:
                        raise ValueError(f"{path}:{line_number}: duplicate metadata")
                    if row.get("schema") != SAMPLE_SCHEMA:
                        raise ValueError(
                            f"{path}:{line_number}: unsupported sample schema"
                        )
                    metadata.append(row)
                    seen_metadata = True
                elif row.get("record_type") == "sample":
                    records.append(row)
                else:
                    raise ValueError(f"{path}:{line_number}: unknown record type")
        if not seen_metadata:
            raise ValueError(f"{path}: missing metadata")
    return metadata, records


def compatible_fingerprints(metadata: list[dict[str, object]]) -> dict[str, object]:
    if not metadata:
        raise ValueError("no sample metadata")
    reference = metadata[0].get("fingerprint")
    if not isinstance(reference, dict):
        raise ValueError("sample metadata has no fingerprint")
    for row in metadata[1:]:
        if row.get("fingerprint") != reference:
            raise ValueError(
                "sample files have incompatible machine/build fingerprints"
            )
    return reference


def finite_number(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def collapse_samples(
    records: list[dict[str, object]], residual_threshold: float
) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in records:
        if row.get("status") != "ok":
            continue
        matrix = str(row.get("matrix_name", "")).lower()
        variant = str(row.get("variant", ""))
        residual = finite_number(row.get("relative_residual_l2"))
        refactor_residual = finite_number(row.get("refactor_max_relative_residual"))
        if not matrix or variant not in VARIANTS or residual is None:
            continue
        if residual > residual_threshold:
            continue
        if refactor_residual is not None and refactor_residual > residual_threshold:
            continue
        grouped[(matrix, variant)].append(row)

    collapsed: list[dict[str, object]] = []
    for (matrix, variant), rows in sorted(grouped.items()):
        result: dict[str, object] = {
            "matrix_name": matrix,
            "matrix": rows[0].get("matrix", matrix),
            "variant": variant,
            "passes_ok": len(rows),
        }
        for key in (*NUMERIC_SAMPLE_FIELDS, *SYMBOLIC_NUMERIC_FIELDS):
            values = [finite_number(row.get(key)) for row in rows]
            usable = [value for value in values if value is not None]
            if usable:
                result[key] = statistics.median(usable)
                median = float(result[key])
                if key.endswith("_seconds") or "seconds_" in key:
                    deviations = [abs(value - median) for value in usable]
                    result[f"{key}_relative_mad"] = (
                        statistics.median(deviations) / median if median > 0.0 else 0.0
                    )
        for key in (*TEXT_SAMPLE_FIELDS, *SYMBOLIC_TEXT_FIELDS):
            if key in rows[0]:
                result[key] = rows[0][key]
        collapsed.append(result)
    return collapsed


def feature_value(row: dict[str, object], name: str) -> float:
    n = finite_number(row.get("n")) or 0.0
    nnz = finite_number(row.get("nnz")) or 0.0
    fill = (finite_number(row.get("symbolic_nnz_l", row.get("nnz_l"))) or 0.0) + (
        finite_number(row.get("symbolic_nnz_u", row.get("nnz_u"))) or 0.0
    )
    if name == "nnz_log2_n":
        return nnz * math.log2(max(n, 2.0))
    if name == "estimated_flops":
        return max(
            finite_number(
                row.get("symbolic_estimated_flops", row.get("estimated_flops"))
            )
            or 0.0,
            0.0,
        )
    if name == "factor_entries":
        return max(fill, 0.0)
    raise ValueError(f"unknown feature {name}")


def solve_linear(matrix: list[list[float]], rhs: list[float]) -> list[float] | None:
    n = len(rhs)
    augmented = [list(matrix[row]) + [rhs[row]] for row in range(n)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda row: abs(augmented[row][col]))
        if abs(augmented[pivot][col]) <= 1.0e-18:
            return None
        augmented[col], augmented[pivot] = augmented[pivot], augmented[col]
        scale = augmented[col][col]
        for item in range(col, n + 1):
            augmented[col][item] /= scale
        for row in range(n):
            if row == col:
                continue
            scale = augmented[row][col]
            for item in range(col, n + 1):
                augmented[row][item] -= scale * augmented[col][item]
    return [augmented[row][n] for row in range(n)]


def percentile(values: list[float], fraction: float) -> float:
    if not values:
        return math.nan
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def predict_stage(model: dict[str, object], row: dict[str, object]) -> float:
    coefficients = model.get("coefficients")
    if not isinstance(coefficients, dict):
        raise ValueError("stage model has no coefficients")
    return sum(
        float(coefficient) * feature_value(row, name)
        for name, coefficient in coefficients.items()
    )


def fit_stage(
    rows: list[dict[str, object]], timing_key: str, feature_names: tuple[str, ...]
) -> dict[str, object]:
    examples: list[tuple[list[float], float, float]] = []
    for row in rows:
        timing = finite_number(row.get(timing_key))
        features = [feature_value(row, name) for name in feature_names]
        if (
            timing is None
            or timing <= 0.0
            or not any(value > 0.0 for value in features)
        ):
            continue
        noise = finite_number(row.get(f"{timing_key}_relative_mad")) or 0.0
        examples.append((features, timing, noise))
    if not examples:
        raise ValueError(f"no usable {timing_key} samples")

    active_sets = [tuple(range(len(feature_names)))]
    if len(feature_names) > 1:
        active_sets.extend((index,) for index in range(len(feature_names)))
    candidates: list[tuple[float, int, list[float], list[float]]] = []
    for active in active_sets:
        weights = [1.0] * len(examples)
        coefficients: list[float] | None = None
        for _ in range(6):
            size = len(active)
            normal = [[0.0] * size for _ in range(size)]
            rhs = [0.0] * size
            for weight, (features, timing, _noise) in zip(weights, examples):
                normalized = [features[index] / timing for index in active]
                for left in range(size):
                    rhs[left] += weight * normalized[left]
                    for right in range(size):
                        normal[left][right] += (
                            weight * normalized[left] * normalized[right]
                        )
            diagonal = max((normal[index][index] for index in range(size)), default=1.0)
            for index in range(size):
                normal[index][index] += max(diagonal * 1.0e-12, 1.0e-24)
            coefficients = solve_linear(normal, rhs)
            if coefficients is None or any(value < 0.0 for value in coefficients):
                coefficients = None
                break
            residuals = []
            for features, timing, _noise in examples:
                predicted = sum(
                    coefficient * features[index]
                    for index, coefficient in zip(active, coefficients)
                )
                residuals.append(abs(predicted / timing - 1.0))
            weights = [1.0 if value <= 0.25 else 0.25 / value for value in residuals]
        if coefficients is None:
            continue
        errors = []
        expanded = [0.0] * len(feature_names)
        for index, coefficient in zip(active, coefficients):
            expanded[index] = coefficient
        for features, timing, _noise in examples:
            predicted = sum(
                value * coefficient for value, coefficient in zip(features, expanded)
            )
            errors.append(abs(predicted / timing - 1.0))
        score = statistics.median(errors) + 0.25 * percentile(errors, 0.9)
        candidates.append((score, len(active), expanded, errors))
    if not candidates:
        raise ValueError(f"could not fit a nonnegative {timing_key} model")
    candidates.sort(key=lambda item: (item[0], item[1]))
    _score, _active_count, coefficients, errors = candidates[0]
    noises = [noise for _features, _timing, noise in examples]
    return {
        "timing": timing_key,
        "coefficients": {
            name: coefficient
            for name, coefficient in zip(feature_names, coefficients)
            if coefficient > 0.0
        },
        "samples": len(examples),
        "median_absolute_relative_error": statistics.median(errors),
        "p90_absolute_relative_error": percentile(errors, 0.9),
        "median_timing_relative_mad": statistics.median(noises),
        "p90_timing_relative_mad": percentile(noises, 0.9),
    }


def pair_feature(
    baseline: dict[str, object], candidate: dict[str, object], name: str
) -> float | None:
    if name == "log_estimated_flops_ratio":
        baseline_value = feature_value(baseline, "estimated_flops")
        candidate_value = feature_value(candidate, "estimated_flops")
    elif name == "log_factor_entries_ratio":
        baseline_value = feature_value(baseline, "factor_entries")
        candidate_value = feature_value(candidate, "factor_entries")
    else:
        raise ValueError(f"unknown pair feature {name}")
    if baseline_value <= 0.0 or candidate_value <= 0.0:
        return None
    return math.log(candidate_value / baseline_value)


def fit_log_ratio_stage(
    pairs: list[tuple[dict[str, object], dict[str, object]]],
    timing_key: str,
    feature_names: tuple[str, ...],
) -> dict[str, object]:
    examples: list[tuple[list[float], float]] = []
    for baseline, candidate in pairs:
        baseline_time = finite_number(baseline.get(timing_key))
        candidate_time = finite_number(candidate.get(timing_key))
        features = [pair_feature(baseline, candidate, name) for name in feature_names]
        if (
            baseline_time is None
            or candidate_time is None
            or baseline_time <= 0.0
            or candidate_time <= 0.0
            or any(value is None for value in features)
        ):
            continue
        examples.append(
            (
                [float(value) for value in features],
                math.log(candidate_time / baseline_time),
            )
        )
    if len(examples) < 3:
        raise ValueError(f"too few paired {timing_key} samples")

    active_sets = [tuple(range(len(feature_names))), tuple()]
    if len(feature_names) > 1:
        active_sets.extend((index,) for index in range(len(feature_names)))
    candidates: list[tuple[float, int, float, list[float], list[float]]] = []
    for active in active_sets:
        weights = [1.0] * len(examples)
        solution: list[float] | None = None
        for _ in range(6):
            size = 1 + len(active)
            normal = [[0.0] * size for _ in range(size)]
            rhs = [0.0] * size
            for weight, (features, target) in zip(weights, examples):
                design = [1.0, *(features[index] for index in active)]
                for left in range(size):
                    rhs[left] += weight * design[left] * target
                    for right in range(size):
                        normal[left][right] += weight * design[left] * design[right]
            diagonal = max((normal[index][index] for index in range(size)), default=1.0)
            for index in range(1, size):
                normal[index][index] += max(diagonal * 1.0e-10, 1.0e-18)
            solution = solve_linear(normal, rhs)
            if solution is None or any(value < 0.0 for value in solution[1:]):
                solution = None
                break
            residuals = []
            for features, target in examples:
                predicted = solution[0] + sum(
                    coefficient * features[index]
                    for index, coefficient in zip(active, solution[1:])
                )
                residuals.append(abs(predicted - target))
            weights = [1.0 if value <= 0.25 else 0.25 / value for value in residuals]
        if solution is None:
            continue
        expanded = [0.0] * len(feature_names)
        for index, coefficient in zip(active, solution[1:]):
            expanded[index] = coefficient
        relative_errors = []
        for features, target in examples:
            predicted = solution[0] + sum(
                coefficient * value for coefficient, value in zip(expanded, features)
            )
            relative_errors.append(abs(math.exp(predicted - target) - 1.0))
        score = statistics.median(relative_errors) + 0.25 * percentile(
            relative_errors, 0.9
        )
        candidates.append((score, len(active), solution[0], expanded, relative_errors))
    if not candidates:
        raise ValueError(f"could not fit paired {timing_key} model")
    candidates.sort(key=lambda item: (item[0], item[1]))
    _score, _active_count, intercept, coefficients, errors = candidates[0]
    return {
        "timing": timing_key,
        "intercept": intercept,
        "coefficients": {
            name: coefficient
            for name, coefficient in zip(feature_names, coefficients)
            if coefficient > 0.0
        },
        "samples": len(examples),
        "median_absolute_relative_error": statistics.median(errors),
        "p90_absolute_relative_error": percentile(errors, 0.9),
    }


def fit_pairwise_models(
    collapsed: list[dict[str, object]], baseline: str, minimum_samples: int
) -> dict[str, object]:
    by_matrix: dict[str, dict[str, dict[str, object]]] = defaultdict(dict)
    for row in collapsed:
        by_matrix[str(row["matrix_name"])][str(row["variant"])] = row
    models: dict[str, object] = {}
    for candidate in sorted(VARIANTS):
        if candidate == baseline:
            continue
        pairs = [
            (rows[baseline], rows[candidate])
            for rows in by_matrix.values()
            if baseline in rows and candidate in rows
        ]
        if len(pairs) < minimum_samples:
            continue
        stages: dict[str, object] = {}
        try:
            stages["factor"] = fit_log_ratio_stage(
                pairs,
                "initial_factor_seconds",
                ("log_estimated_flops_ratio", "log_factor_entries_ratio"),
            )
            stages["refactor"] = fit_log_ratio_stage(
                pairs,
                "refactor_steady_seconds_avg",
                ("log_estimated_flops_ratio", "log_factor_entries_ratio"),
            )
            stages["solve"] = fit_log_ratio_stage(
                pairs,
                "refactor_solve_steady_seconds_avg",
                ("log_factor_entries_ratio",),
            )
        except ValueError:
            continue
        if min(int(model["samples"]) for model in stages.values()) < minimum_samples:
            continue
        models[candidate] = {
            "baseline": baseline,
            "candidate": candidate,
            "matrices": len(pairs),
            "stages": stages,
        }
    return models


def fit_profile(
    metadata: list[dict[str, object]],
    records: list[dict[str, object]],
    residual_threshold: float,
    minimum_samples: int,
    baseline: str = "amd_btf",
) -> dict[str, object]:
    fingerprint = compatible_fingerprints(metadata)
    collapsed = collapse_samples(records, residual_threshold)
    by_variant: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in collapsed:
        by_variant[str(row["variant"])].append(row)
    models: dict[str, object] = {}
    for variant, rows in sorted(by_variant.items()):
        if len(rows) < minimum_samples:
            continue
        stages: dict[str, object] = {}
        for stage, (timing_key, feature_names) in STAGES.items():
            stages[stage] = fit_stage(rows, timing_key, feature_names)
        models[variant] = {
            "definition": VARIANTS[variant],
            "matrices": len(rows),
            "stages": stages,
        }
    if len(models) < 2:
        raise ValueError(
            "fewer than two variants have enough valid calibration matrices"
        )
    if baseline not in models:
        raise ValueError(f"baseline variant {baseline} lacks a complete absolute model")
    pairwise_models = fit_pairwise_models(collapsed, baseline, minimum_samples)
    if not pairwise_models:
        raise ValueError("no candidate has enough paired symbolic/numeric samples")
    timing_noise = []
    for model in models.values():
        assert isinstance(model, dict)
        stages = model["stages"]
        assert isinstance(stages, dict)
        for stage_model in stages.values():
            assert isinstance(stage_model, dict)
            timing_noise.append(float(stage_model["p90_timing_relative_mad"]))
    minimum_margin = max(0.03, min(0.20, 4.0 * statistics.median(timing_noise)))
    manifest_hashes = sorted(
        {
            str(manifest["sha256"])
            for row in metadata
            for manifest in row.get("manifests", [])
            if isinstance(manifest, dict) and "sha256" in manifest
        }
    )
    return {
        "schema": PROFILE_SCHEMA,
        "created_utc": utc_now(),
        "fingerprint": fingerprint,
        "calibration": {
            "sample_files": len(metadata),
            "source_manifest_sha256s": manifest_hashes,
            "successful_matrices": len({str(row["matrix_name"]) for row in collapsed}),
            "successful_matrix_variants": len(collapsed),
            "residual_threshold": residual_threshold,
            "minimum_samples_per_variant": minimum_samples,
        },
        "decision": {
            "baseline": baseline,
            "minimum_relative_margin": minimum_margin,
            "horizon_source": (
                "caller-supplied; optional calibrated margins are keyed by horizon"
            ),
            "calibrated_relative_margins": {},
        },
        "models": models,
        "pairwise_models": pairwise_models,
    }


def fit_command(args: argparse.Namespace) -> int:
    if args.minimum_samples < 3:
        raise ValueError("minimum samples must be at least three")
    if args.residual_threshold <= 0.0:
        raise ValueError("residual threshold must be positive")
    metadata, records = load_sample_files(args.samples)
    profile = fit_profile(
        metadata,
        records,
        args.residual_threshold,
        args.minimum_samples,
        args.baseline,
    )
    args.profile.parent.mkdir(parents=True, exist_ok=True)
    args.profile.write_text(
        json.dumps(profile, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    summary = {
        name: {
            "matrices": model["matrices"],
            "stage_p90_relative_error": {
                stage: details["p90_absolute_relative_error"]
                for stage, details in model["stages"].items()
            },
        }
        for name, model in profile["models"].items()
    }
    summary["pairwise"] = {
        name: {
            "matrices": model["matrices"],
            "stage_p90_relative_error": {
                stage: details["p90_absolute_relative_error"]
                for stage, details in model["stages"].items()
            },
        }
        for name, model in profile["pairwise_models"].items()
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"wrote {args.profile}")
    return 0


def load_profile(path: pathlib.Path) -> dict[str, object]:
    profile = json.loads(path.read_text(encoding="utf-8"))
    if profile.get("schema") != PROFILE_SCHEMA:
        raise ValueError(f"{path}: unsupported profile schema")
    return profile


def stage_prediction(
    variant_model: dict[str, object], row: dict[str, object], stage: str
) -> tuple[float, float]:
    stages = variant_model.get("stages")
    if not isinstance(stages, dict) or not isinstance(stages.get(stage), dict):
        raise ValueError(f"profile has no {stage} model")
    model = stages[stage]
    assert isinstance(model, dict)
    predicted = predict_stage(model, row)
    uncertainty = max(
        float(model.get("p90_absolute_relative_error", 0.0)),
        float(model.get("p90_timing_relative_mad", 0.0)),
    )
    return predicted, uncertainty


def predict_pair_stage_ratio(
    stage_model: dict[str, object],
    baseline: dict[str, object],
    candidate: dict[str, object],
) -> tuple[float, float]:
    raw_coefficients = stage_model.get("coefficients")
    if not isinstance(raw_coefficients, dict):
        raise ValueError("pairwise stage model has no coefficients")
    value = float(stage_model.get("intercept", 0.0))
    for name, coefficient in raw_coefficients.items():
        feature = pair_feature(baseline, candidate, str(name))
        if feature is None:
            raise ValueError("candidate symbolic has no positive work estimate")
        value += float(coefficient) * feature
    return math.exp(value), float(stage_model["p90_absolute_relative_error"])


def symbolic_analysis_seconds(row: dict[str, object]) -> float:
    return max(
        finite_number(row.get("symbolic_analysis_seconds", row.get("analysis_seconds")))
        or 0.0,
        0.0,
    )


def predicted_pair_lifecycle(
    pair_model: dict[str, object],
    baseline: dict[str, object],
    candidate: dict[str, object],
    expected_refactors: int,
    expected_solves: int,
    baseline_numeric: dict[str, tuple[float, float]] | None = None,
) -> tuple[float, float, dict[str, float]]:
    raw_stages = pair_model.get("stages")
    if not isinstance(raw_stages, dict):
        raise ValueError("pairwise model has no stages")
    if baseline_numeric is None:
        baseline_numeric = {
            "factor": (
                finite_number(baseline.get("initial_factor_seconds")) or 0.0,
                0.0,
            ),
            "refactor": (
                finite_number(baseline.get("refactor_steady_seconds_avg")) or 0.0,
                0.0,
            ),
            "solve": (
                finite_number(baseline.get("refactor_solve_steady_seconds_avg")) or 0.0,
                0.0,
            ),
        }
    components = {"analysis": symbolic_analysis_seconds(candidate)}
    weighted_error = 0.0
    for stage, multiplier in (
        ("factor", 1),
        ("refactor", expected_refactors),
        ("solve", expected_solves),
    ):
        baseline_time, baseline_uncertainty = baseline_numeric[stage]
        raw_model = raw_stages.get(stage)
        if not isinstance(raw_model, dict):
            raise ValueError(f"pairwise model has no {stage} stage")
        ratio, uncertainty = predict_pair_stage_ratio(raw_model, baseline, candidate)
        contribution = multiplier * baseline_time * ratio
        components[stage] = contribution
        weighted_error += contribution * (baseline_uncertainty + uncertainty)
    total = sum(components.values())
    return total, weighted_error / total if total > 0.0 else math.inf, components


def predicted_baseline_lifecycle(
    variant_model: dict[str, object],
    row: dict[str, object],
    expected_refactors: int,
    expected_solves: int,
) -> tuple[
    float,
    float,
    dict[str, float],
    dict[str, tuple[float, float]],
]:
    numeric: dict[str, tuple[float, float]] = {}
    components = {"analysis": symbolic_analysis_seconds(row)}
    weighted_error = 0.0
    for stage, multiplier in (
        ("factor", 1),
        ("refactor", expected_refactors),
        ("solve", expected_solves),
    ):
        prediction, uncertainty = stage_prediction(variant_model, row, stage)
        numeric[stage] = (prediction, uncertainty)
        contribution = multiplier * prediction
        components[stage] = contribution
        weighted_error += contribution * uncertainty
    total = sum(components.values())
    return (
        total,
        weighted_error / total if total > 0.0 else math.inf,
        components,
        numeric,
    )


def actual_lifecycle(
    row: dict[str, object], expected_refactors: int, expected_solves: int
) -> float:
    analysis = finite_number(row.get("analysis_seconds")) or 0.0
    factor = finite_number(row.get("initial_factor_seconds")) or 0.0
    refactor = finite_number(row.get("refactor_steady_seconds_avg")) or 0.0
    solve = finite_number(row.get("refactor_solve_steady_seconds_avg")) or 0.0
    return analysis + factor + expected_refactors * refactor + expected_solves * solve


def calibrated_margin_key(
    prediction_mode: str, expected_refactors: int, expected_solves: int
) -> str:
    return f"{prediction_mode}:refactors={expected_refactors}:solves={expected_solves}"


def profile_calibrated_margin(
    profile: dict[str, object],
    prediction_mode: str,
    expected_refactors: int,
    expected_solves: int,
    regression_tolerance: float,
) -> float | None:
    decision = profile.get("decision")
    if not isinstance(decision, dict):
        return None
    margins = decision.get("calibrated_relative_margins")
    if not isinstance(margins, dict):
        return None
    entry = margins.get(
        calibrated_margin_key(prediction_mode, expected_refactors, expected_solves)
    )
    if not isinstance(entry, dict):
        return None
    calibrated_tolerance = finite_number(entry.get("regression_tolerance"))
    margin = finite_number(entry.get("relative_margin"))
    if (
        calibrated_tolerance is None
        or margin is None
        or not math.isclose(
            calibrated_tolerance, regression_tolerance, rel_tol=0.0, abs_tol=1.0e-15
        )
    ):
        return None
    return margin


def evaluate_profile(
    profile: dict[str, object],
    metadata: list[dict[str, object]],
    records: list[dict[str, object]],
    expected_refactors: int,
    expected_solves: int,
    baseline: str,
    residual_threshold: float,
    regression_tolerance: float,
    prediction_mode: str = "symbolic-only",
    decision_margin: float | None = None,
) -> dict[str, object]:
    if prediction_mode not in PREDICTION_MODES:
        raise ValueError(f"unknown prediction mode {prediction_mode}")
    if decision_margin is not None and not 0.0 <= decision_margin < 1.0:
        raise ValueError("decision margin must be in [0, 1)")
    sample_fingerprint = compatible_fingerprints(metadata)
    if sample_fingerprint != profile.get("fingerprint"):
        raise ValueError(
            "profile and validation samples have incompatible fingerprints"
        )
    models = profile.get("models")
    if not isinstance(models, dict) or baseline not in models:
        raise ValueError(f"profile has no baseline model {baseline}")
    pairwise_models = profile.get("pairwise_models")
    if not isinstance(pairwise_models, dict):
        raise ValueError("profile has no pairwise models")
    profile_baseline = profile.get("decision", {}).get("baseline")  # type: ignore[union-attr]
    if profile_baseline != baseline:
        raise ValueError(
            f"profile was fitted against {profile_baseline}, not requested {baseline}"
        )
    collapsed = collapse_samples(records, residual_threshold)
    by_matrix: dict[str, dict[str, dict[str, object]]] = defaultdict(dict)
    for row in collapsed:
        variant = str(row["variant"])
        if variant in models:
            by_matrix[str(row["matrix_name"])][variant] = row
    minimum_margin = float(
        profile.get("decision", {}).get("minimum_relative_margin", 0.05)  # type: ignore[union-attr]
    )
    fixed_margin = decision_margin
    margin_source = "command-line override"
    if fixed_margin is None:
        fixed_margin = profile_calibrated_margin(
            profile,
            prediction_mode,
            expected_refactors,
            expected_solves,
            regression_tolerance,
        )
        margin_source = (
            "calibrated profile" if fixed_margin is not None else "model error bound"
        )
    results = []
    false_positives = 0
    decisive = 0
    correct = 0
    for matrix, rows in sorted(by_matrix.items()):
        if baseline not in rows or len(rows) < 2:
            continue
        actuals: dict[str, float] = {
            variant: actual_lifecycle(row, expected_refactors, expected_solves)
            for variant, row in rows.items()
        }
        if prediction_mode == "baseline-measured":
            predictions: dict[str, float] = {baseline: actuals[baseline]}
            uncertainties: dict[str, float] = {baseline: 0.0}
            components: dict[str, dict[str, float]] = {
                baseline: {
                    "analysis": symbolic_analysis_seconds(rows[baseline]),
                    "factor": finite_number(
                        rows[baseline].get("initial_factor_seconds")
                    )
                    or 0.0,
                    "refactor": expected_refactors
                    * (
                        finite_number(rows[baseline].get("refactor_steady_seconds_avg"))
                        or 0.0
                    ),
                    "solve": expected_solves
                    * (
                        finite_number(
                            rows[baseline].get("refactor_solve_steady_seconds_avg")
                        )
                        or 0.0
                    ),
                }
            }
            baseline_numeric = None
        else:
            (
                baseline_prediction,
                baseline_uncertainty,
                baseline_components,
                baseline_numeric,
            ) = predicted_baseline_lifecycle(
                models[baseline],
                rows[baseline],
                expected_refactors,
                expected_solves,
            )
            predictions = {baseline: baseline_prediction}
            uncertainties = {baseline: baseline_uncertainty}
            components = {baseline: baseline_components}
        for variant, pair_model in pairwise_models.items():
            if variant not in rows or not isinstance(pair_model, dict):
                continue
            try:
                prediction, uncertainty, breakdown = predicted_pair_lifecycle(
                    pair_model,
                    rows[baseline],
                    rows[variant],
                    expected_refactors,
                    expected_solves,
                    baseline_numeric,
                )
            except ValueError:
                continue
            if prediction <= 0.0 or not math.isfinite(prediction):
                continue
            predictions[variant] = prediction
            uncertainties[variant] = uncertainty
            components[variant] = breakdown
        if baseline not in predictions or len(predictions) < 2:
            continue
        predicted_best = min(predictions, key=predictions.get)  # type: ignore[arg-type]
        actual_best = min(actuals, key=actuals.get)  # type: ignore[arg-type]
        predicted_advantage = 1.0 - predictions[predicted_best] / predictions[baseline]
        if fixed_margin is None:
            comparison_margin = min(
                0.75,
                max(
                    minimum_margin,
                    uncertainties[predicted_best] + uncertainties[baseline],
                ),
            )
        else:
            comparison_margin = fixed_margin
        selected = baseline
        verdict = "baseline_predicted"
        if predicted_best != baseline:
            if predicted_advantage > comparison_margin:
                selected = predicted_best
                verdict = "candidate_selected"
                decisive += 1
            else:
                verdict = "candidate_uncertain"
        actual_ratio = actuals[selected] / actuals[baseline]
        oracle_ratio = actuals[actual_best] / actuals[baseline]
        if selected == baseline or actual_ratio <= 1.0:
            correct += 1
        if selected != baseline and actual_ratio > 1.0 + regression_tolerance:
            false_positives += 1
        results.append(
            {
                "matrix": matrix,
                "predicted_best": predicted_best,
                "selected": selected,
                "actual_best": actual_best,
                "verdict": verdict,
                "predicted_advantage_over_baseline": predicted_advantage,
                "comparison_margin": comparison_margin,
                "actual_selected_over_baseline": actual_ratio,
                "actual_oracle_over_baseline": oracle_ratio,
                "predictions": predictions,
                "actuals": actuals,
                "prediction_components": components,
            }
        )
    return {
        "schema": REPORT_SCHEMA,
        "created_utc": utc_now(),
        "profile_schema": PROFILE_SCHEMA,
        "horizon": {
            "expected_refactors": expected_refactors,
            "expected_solves": expected_solves,
        },
        "baseline": baseline,
        "prediction_mode": prediction_mode,
        "decision_margin": fixed_margin,
        "decision_margin_source": margin_source,
        "matrices": len(results),
        "decisive_nonbaseline": decisive,
        "safe_decisions": correct,
        "false_positive_regressions": false_positives,
        "regression_tolerance": regression_tolerance,
        "results": results,
    }


def evaluate_command(args: argparse.Namespace) -> int:
    if args.expected_refactors < 0 or args.expected_solves < 0:
        raise ValueError("expected lifecycle counts must be non-negative")
    if args.regression_tolerance < 0.0:
        raise ValueError("regression tolerance must be non-negative")
    profile = load_profile(args.profile)
    metadata, records = load_sample_files(args.samples)
    report = evaluate_profile(
        profile,
        metadata,
        records,
        args.expected_refactors,
        args.expected_solves,
        args.baseline,
        args.residual_threshold,
        args.regression_tolerance,
        args.prediction_mode,
        args.decision_margin,
    )
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    for row in report["results"]:
        print(
            f"{row['matrix']}: predicted={row['predicted_best']} "
            f"selected={row['selected']} actual={row['actual_best']} "
            f"pred_gain={row['predicted_advantage_over_baseline']:.3f} "
            f"margin={row['comparison_margin']:.3f} "
            f"actual/base={row['actual_selected_over_baseline']:.3f}"
        )
    print(
        json.dumps(
            {key: value for key, value in report.items() if key != "results"},
            indent=2,
            sort_keys=True,
        )
    )
    return 0 if report["false_positive_regressions"] == 0 else 2


def calibrate_margin_command(args: argparse.Namespace) -> int:
    if args.expected_refactors < 0 or args.expected_solves < 0:
        raise ValueError("expected lifecycle counts must be non-negative")
    if not 0.0 <= args.regression_tolerance < 1.0:
        raise ValueError("regression tolerance must be in [0, 1)")
    profile = load_profile(args.profile)
    metadata, records = load_sample_files(args.samples)
    report = evaluate_profile(
        profile,
        metadata,
        records,
        args.expected_refactors,
        args.expected_solves,
        args.baseline,
        args.residual_threshold,
        args.regression_tolerance,
        args.prediction_mode,
        0.0,
    )
    unsafe_gains = []
    for row in report["results"]:
        predicted = str(row["predicted_best"])
        if predicted == args.baseline:
            continue
        actuals = row["actuals"]
        actual_ratio = float(actuals[predicted]) / float(actuals[args.baseline])
        if actual_ratio > 1.0 + args.regression_tolerance:
            unsafe_gains.append(float(row["predicted_advantage_over_baseline"]))
    decision = profile.get("decision")
    if not isinstance(decision, dict):
        raise ValueError("profile has no decision metadata")
    minimum_margin = float(decision.get("minimum_relative_margin", 0.05))
    unsafe_floor = max(unsafe_gains, default=0.0)
    margin = max(minimum_margin, unsafe_floor)
    if margin >= 1.0:
        raise ValueError("calibrated decision margin is not finite below one")
    margins = decision.setdefault("calibrated_relative_margins", {})
    if not isinstance(margins, dict):
        raise ValueError("profile calibrated margins are malformed")
    key = calibrated_margin_key(
        args.prediction_mode, args.expected_refactors, args.expected_solves
    )
    margins[key] = {
        "calibrated_utc": utc_now(),
        "calibration_matrices": report["matrices"],
        "prediction_mode": args.prediction_mode,
        "expected_refactors": args.expected_refactors,
        "expected_solves": args.expected_solves,
        "regression_tolerance": args.regression_tolerance,
        "minimum_noise_margin": minimum_margin,
        "largest_unsafe_predicted_gain": unsafe_floor,
        "unsafe_raw_choices": len(unsafe_gains),
        "relative_margin": margin,
    }
    args.output_profile.parent.mkdir(parents=True, exist_ok=True)
    args.output_profile.write_text(
        json.dumps(profile, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "key": key,
                "calibration_matrices": report["matrices"],
                "unsafe_raw_choices": len(unsafe_gains),
                "largest_unsafe_predicted_gain": unsafe_floor,
                "minimum_noise_margin": minimum_margin,
                "relative_margin": margin,
            },
            indent=2,
            sort_keys=True,
        )
    )
    print(f"wrote {args.output_profile}")
    return 0


def synthetic_row(
    matrix: str,
    variant: str,
    n: int,
    nnz: int,
    fill: float,
    flops: float,
    rates: tuple[float, float, float, float, float, float],
) -> dict[str, object]:
    (
        analysis_rate,
        factor_flop,
        factor_fill,
        refactor_flop,
        refactor_fill,
        solve_rate,
    ) = rates
    return {
        "record_type": "sample",
        "status": "ok",
        "matrix": matrix + ".mtx",
        "matrix_name": matrix,
        "variant": variant,
        "pass": 0,
        "n": n,
        "nnz": nnz,
        "analysis_seconds": analysis_rate * nnz * math.log2(n),
        "initial_factor_seconds": factor_flop * flops + factor_fill * fill,
        "refactor_steady_seconds_avg": refactor_flop * flops + refactor_fill * fill,
        "refactor_solve_steady_seconds_avg": solve_rate * fill,
        "relative_residual_l2": 1.0e-14,
        "refactor_max_relative_residual": 1.0e-14,
        "nnz_l": fill / 2.0,
        "nnz_u": fill / 2.0,
        "estimated_flops": flops,
    }


def self_test() -> int:
    fingerprint = {"test": True}
    metadata = [
        {
            "record_type": "metadata",
            "schema": SAMPLE_SCHEMA,
            "fingerprint": fingerprint,
            "manifests": [],
        }
    ]
    records = []
    amd_rates = (2.0e-9, 1.2e-9, 4.0e-8, 4.0e-10, 1.0e-8, 2.0e-8)
    amf_rates = (2.4e-9, 1.0e-9, 3.0e-8, 2.0e-10, 8.0e-9, 1.5e-8)
    for index in range(6):
        n = 2000 * (index + 1)
        nnz = 7 * n
        fill = float((20 + index) * n)
        flops = float((300 + 50 * index) * n)
        records.append(
            synthetic_row(f"m{index}", "amd_btf", n, nnz, fill, flops, amd_rates)
        )
        records.append(
            synthetic_row(
                f"m{index}", "amf_btf", n, nnz, fill * 0.8, flops * 0.7, amf_rates
            )
        )
    profile = fit_profile(metadata, records, 1.0e-8, 3)
    report = evaluate_profile(
        profile, metadata, records, 20, 20, "amd_btf", 1.0e-8, 0.02
    )
    if set(profile["models"]) != {"amd_btf", "amf_btf"}:
        raise AssertionError("self-test profile variants are incomplete")
    if (
        report["matrices"] != 6
        or report["decisive_nonbaseline"] != 6
        or report["false_positive_regressions"] != 0
        or report["prediction_mode"] != "symbolic-only"
    ):
        raise AssertionError("self-test evaluation failed")
    key = calibrated_margin_key("symbolic-only", 20, 20)
    profile["decision"]["calibrated_relative_margins"][key] = {  # type: ignore[index]
        "regression_tolerance": 0.02,
        "relative_margin": 0.5,
    }
    guarded = evaluate_profile(
        profile, metadata, records, 20, 20, "amd_btf", 1.0e-8, 0.02
    )
    if (
        guarded["decision_margin_source"] != "calibrated profile"
        or guarded["decision_margin"] != 0.5
        or guarded["decisive_nonbaseline"] != 0
    ):
        raise AssertionError("self-test calibrated margin was not applied")
    encoded = json.dumps(profile, sort_keys=True)
    if PROFILE_SCHEMA not in encoded or "m0" in encoded:
        raise AssertionError("profile schema leaked calibration matrix identity")
    print("machine cost calibration self-test: passed")
    return 0


def add_collection_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--kls-bench", type=pathlib.Path, required=True)
    parser.add_argument("--matrix-dir", type=pathlib.Path, required=True)
    parser.add_argument("--manifest", type=pathlib.Path, action="append", required=True)
    parser.add_argument("--samples", type=pathlib.Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--cpu-list")
    parser.add_argument("--passes", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument("--refactor-repeat", type=int, default=5)
    parser.add_argument("--refactor-value-amplitude", type=float, default=0.001)
    parser.add_argument(
        "--scale", choices=("auto", "-1", "0", "1", "2"), default="auto"
    )
    parser.add_argument("--variants", default=",".join(DEFAULT_VARIANTS))
    parser.add_argument("--skip", type=int, default=0)
    parser.add_argument("--limit", type=int)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    collect_parser = subparsers.add_parser(
        "collect", help="collect explicit-variant samples"
    )
    add_collection_arguments(collect_parser)

    fit_parser = subparsers.add_parser("fit", help="fit and persist a machine profile")
    fit_parser.add_argument(
        "--samples", type=pathlib.Path, action="append", required=True
    )
    fit_parser.add_argument("--profile", type=pathlib.Path, required=True)
    fit_parser.add_argument("--baseline", choices=tuple(VARIANTS), default="amd_btf")
    fit_parser.add_argument("--residual-threshold", type=float, default=1.0e-8)
    fit_parser.add_argument("--minimum-samples", type=int, default=4)

    evaluate_parser = subparsers.add_parser(
        "evaluate", help="evaluate a profile on separately collected samples"
    )
    evaluate_parser.add_argument("--profile", type=pathlib.Path, required=True)
    evaluate_parser.add_argument(
        "--samples", type=pathlib.Path, action="append", required=True
    )
    evaluate_parser.add_argument("--expected-refactors", type=int, required=True)
    evaluate_parser.add_argument("--expected-solves", type=int, required=True)
    evaluate_parser.add_argument(
        "--baseline", choices=tuple(VARIANTS), default="amd_btf"
    )
    evaluate_parser.add_argument("--residual-threshold", type=float, default=1.0e-8)
    evaluate_parser.add_argument("--regression-tolerance", type=float, default=0.02)
    evaluate_parser.add_argument(
        "--prediction-mode", choices=PREDICTION_MODES, default="symbolic-only"
    )
    evaluate_parser.add_argument("--decision-margin", type=float)
    evaluate_parser.add_argument("--report", type=pathlib.Path)

    margin_parser = subparsers.add_parser(
        "calibrate-margin",
        help="persist an empirically safe margin for one lifecycle horizon",
    )
    margin_parser.add_argument("--profile", type=pathlib.Path, required=True)
    margin_parser.add_argument("--output-profile", type=pathlib.Path, required=True)
    margin_parser.add_argument(
        "--samples", type=pathlib.Path, action="append", required=True
    )
    margin_parser.add_argument("--expected-refactors", type=int, required=True)
    margin_parser.add_argument("--expected-solves", type=int, required=True)
    margin_parser.add_argument("--baseline", choices=tuple(VARIANTS), default="amd_btf")
    margin_parser.add_argument("--residual-threshold", type=float, default=1.0e-8)
    margin_parser.add_argument("--regression-tolerance", type=float, default=0.02)
    margin_parser.add_argument(
        "--prediction-mode", choices=PREDICTION_MODES, default="symbolic-only"
    )

    subparsers.add_parser("self-test", help="run deterministic model tests")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    try:
        if args.command == "collect":
            return collect(args)
        if args.command == "fit":
            return fit_command(args)
        if args.command == "evaluate":
            return evaluate_command(args)
        if args.command == "calibrate-margin":
            return calibrate_margin_command(args)
        return self_test()
    except (OSError, ValueError) as exc:
        print(exc, file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("interrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
