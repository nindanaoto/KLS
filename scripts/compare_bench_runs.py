#!/usr/bin/env python3
"""Compare two JSONL benchmark runs by matrix basename."""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import re
import statistics


def spice_cycle_seconds(row: dict[str, object]) -> float:
    if "spice_cycle_seconds" in row:
        return float(row["spice_cycle_seconds"])
    analysis = float(row["analysis_seconds"])
    factor = float(row.get("initial_factor_seconds", row["factor_seconds_avg"]))
    solve = float(row["solve_seconds_avg"])
    first = float(row.get("refactor_first_seconds", row["refactor_seconds_avg"]))
    steady = float(
        row.get("refactor_steady_seconds_avg", row["refactor_seconds_avg"])
    )
    return analysis + factor + first + 98.0 * steady + 100.0 * solve


def load_rows(path: pathlib.Path) -> dict[str, dict[str, object]]:
    rows: dict[str, dict[str, object]] = {}
    with path.open("r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, 1):
            if not line.strip():
                continue
            row = json.loads(line)
            matrix = pathlib.Path(str(row["matrix"])).name
            if matrix in rows:
                raise ValueError(f"{path}:{line_no}: duplicate matrix basename {matrix}")
            row["spice_cycle_seconds"] = spice_cycle_seconds(row)
            rows[matrix] = row
    return rows


def matrix_name_variants(name: str) -> set[str]:
    path_name = pathlib.Path(name).name
    variants = {path_name}
    stem = pathlib.Path(path_name).stem
    if stem:
        variants.add(stem)
    return variants


def load_manifest_names(path: pathlib.Path) -> set[str]:
    names: set[str] = set()
    with path.open("r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            token = stripped.split()[0]
            variants = matrix_name_variants(token)
            if not variants:
                raise ValueError(f"{path}:{line_no}: empty matrix name")
            names.update(variants)
    return names


def matches_manifest_name(name: str, manifest_names: set[str]) -> bool:
    return bool(matrix_name_variants(name) & manifest_names)


def filter_matrix_names(
    names: set[str],
    include_names: set[str],
    exclude_names: set[str],
) -> set[str]:
    if include_names:
        names = {name for name in names if matches_manifest_name(name, include_names)}
    if exclude_names:
        names = {name for name in names if not matches_manifest_name(name, exclude_names)}
    return names


def load_failures(path: pathlib.Path) -> dict[str, dict[str, object]]:
    failure_path = path.with_suffix(".failures")
    rows: dict[str, dict[str, object]] = {}
    if not failure_path.exists():
        return rows
    with failure_path.open("r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, 1):
            if not line.strip():
                continue
            row = json.loads(line)
            matrix = pathlib.Path(str(row["matrix"])).name
            if matrix in rows:
                raise ValueError(
                    f"{failure_path}:{line_no}: duplicate matrix basename {matrix}"
                )
            rows[matrix] = row
    return rows


def geometric_mean(values: list[float]) -> float:
    positives = [v for v in values if v > 0.0 and math.isfinite(v)]
    if not positives:
        return math.nan
    return math.exp(sum(math.log(v) for v in positives) / len(positives))


def print_failures(name: str, failures: dict[str, dict[str, object]]) -> None:
    if not failures:
        return
    print(f"\nFailures in {name}:")
    for matrix in sorted(failures):
        reason = failures[matrix].get("reason", "")
        print(f"  {matrix}: {reason}")


def parsed_timeout_seconds(row: dict[str, object]) -> float | None:
    candidates: list[str] = []
    reason = row.get("reason")
    if reason is not None:
        candidates.append(str(reason))
    sample_failures = row.get("sample_failures")
    if isinstance(sample_failures, list):
        candidates.extend(str(item) for item in sample_failures)
    pattern = re.compile(r"timeout after ([0-9]+(?:\.[0-9]+)?)s")
    for text in candidates:
        match = pattern.search(text)
        if match:
            return float(match.group(1))
    return None


def effective_cycle_seconds(
    name: str,
    rows: dict[str, dict[str, object]],
    failures: dict[str, dict[str, object]],
    fallback_seconds: float | None,
) -> tuple[float, str] | None:
    if name in rows:
        return float(rows[name]["spice_cycle_seconds"]), "ok"
    if name in failures:
        seconds = fallback_seconds
        if seconds is None:
            seconds = parsed_timeout_seconds(failures[name])
        if seconds is None:
            return None
        return seconds, "failed"
    if fallback_seconds is None:
        return None
    return fallback_seconds, "missing"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=pathlib.Path, required=True)
    parser.add_argument("--reference", type=pathlib.Path, required=True)
    parser.add_argument("--candidate-name", default="candidate")
    parser.add_argument("--reference-name", default="reference")
    parser.add_argument("--max-rows", type=int, default=12)
    parser.add_argument(
        "--include-manifest",
        action="append",
        type=pathlib.Path,
        default=[],
        help=(
            "only compare matrices listed in this manifest; may be repeated. "
            "Bare SuiteSparse names and .mtx basenames both match."
        ),
    )
    parser.add_argument(
        "--exclude-manifest",
        action="append",
        type=pathlib.Path,
        default=[],
        help=(
            "drop matrices listed in this manifest; may be repeated. "
            "Bare SuiteSparse names and .mtx basenames both match."
        ),
    )
    parser.add_argument(
        "--include-failures",
        action="store_true",
        help=(
            "Score failed or missing matrices instead of restricting the "
            "summary to successful rows common to both runs."
        ),
    )
    parser.add_argument(
        "--failure-seconds",
        type=float,
        default=None,
        help=(
            "Fallback cycle-time penalty for non-timeout failures or missing "
            "rows when --include-failures is used. Timeout failures infer "
            "their process cap from the .failures file if this is omitted; "
            "that cap is only a lower bound for SPICE-cycle comparisons."
        ),
    )
    args = parser.parse_args()

    candidate = load_rows(args.candidate)
    reference = load_rows(args.reference)
    candidate_failures = load_failures(args.candidate)
    reference_failures = load_failures(args.reference)
    include_names: set[str] = set()
    for manifest in args.include_manifest:
        include_names.update(load_manifest_names(manifest))
    exclude_names: set[str] = set()
    for manifest in args.exclude_manifest:
        exclude_names.update(load_manifest_names(manifest))
    if args.include_failures:
        common_set = (
            set(candidate)
            | set(reference)
            | set(candidate_failures)
            | set(reference_failures)
        )
    else:
        common_set = set(candidate) & set(reference)
    common = sorted(filter_matrix_names(common_set, include_names, exclude_names))
    if not common:
        print_failures(args.candidate_name, candidate_failures)
        print_failures(args.reference_name, reference_failures)
        raise SystemExit("no common matrix basenames")

    ratios: list[float] = []
    records: list[tuple[float, str, float, float, str, str]] = []
    for name in common:
        if args.include_failures:
            candidate_effective = effective_cycle_seconds(
                name, candidate, candidate_failures, args.failure_seconds
            )
            reference_effective = effective_cycle_seconds(
                name, reference, reference_failures, args.failure_seconds
            )
            if candidate_effective is None or reference_effective is None:
                continue
            c, c_status = candidate_effective
            r, r_status = reference_effective
        else:
            c = float(candidate[name]["spice_cycle_seconds"])
            r = float(reference[name]["spice_cycle_seconds"])
            c_status = "ok"
            r_status = "ok"
        if c <= 0.0 or r <= 0.0 or not math.isfinite(c) or not math.isfinite(r):
            continue
        ratio = c / r
        ratios.append(ratio)
        records.append((ratio, name, c, r, c_status, r_status))

    if not records:
        raise SystemExit("no finite positive common timings")

    wins = sum(1 for ratio, *_ in records if ratio < 0.98)
    ties = sum(1 for ratio, *_ in records if 0.98 <= ratio <= 1.02)
    losses = sum(1 for ratio, *_ in records if ratio > 1.02)
    speedups = [1.0 / ratio for ratio in ratios]
    summary = {
        "candidate": args.candidate_name,
        "reference": args.reference_name,
        "matrices_common": len(records),
        "include_failures": args.include_failures,
        "failure_seconds": args.failure_seconds,
        "candidate_failed_or_missing_scored": sum(
            1 for _, _, _, _, c_status, _ in records if c_status != "ok"
        ),
        "reference_failed_or_missing_scored": sum(
            1 for _, _, _, _, _, r_status in records if r_status != "ok"
        ),
        "candidate_geomean_seconds": geometric_mean(
            [c for _, _, c, _, _, _ in records]
        ),
        "reference_geomean_seconds": geometric_mean(
            [r for _, _, _, r, _, _ in records]
        ),
        "geomean_ratio_candidate_over_reference": geometric_mean(ratios),
        "geomean_speedup_reference_over_candidate": geometric_mean(speedups),
        "median_ratio_candidate_over_reference": statistics.median(ratios),
        "wins_over_2pct": wins,
        "ties_within_2pct": ties,
        "losses_over_2pct": losses,
    }
    print(json.dumps(summary, indent=2, sort_keys=True))

    print("\nLargest candidate wins:")
    for ratio, name, c, r, c_status, r_status in sorted(records)[: args.max_rows]:
        print(
            f"  {name}: ratio={ratio:.3f} "
            f"{args.candidate_name}={c:.6g}s/{c_status} "
            f"{args.reference_name}={r:.6g}s/{r_status}"
        )

    print("\nLargest candidate losses:")
    for ratio, name, c, r, c_status, r_status in sorted(records, reverse=True)[: args.max_rows]:
        print(
            f"  {name}: ratio={ratio:.3f} "
            f"{args.candidate_name}={c:.6g}s/{c_status} "
            f"{args.reference_name}={r:.6g}s/{r_status}"
        )

    visible_candidate = filter_matrix_names(
        set(candidate), include_names, exclude_names
    )
    visible_reference = filter_matrix_names(
        set(reference), include_names, exclude_names
    )
    missing_candidate = sorted(visible_reference - visible_candidate)
    missing_reference = sorted(visible_candidate - visible_reference)
    if missing_candidate:
        print(f"\nMissing from {args.candidate_name}: {', '.join(missing_candidate)}")
    if missing_reference:
        print(f"\nMissing from {args.reference_name}: {', '.join(missing_reference)}")
    print_failures(args.candidate_name, candidate_failures)
    print_failures(args.reference_name, reference_failures)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
