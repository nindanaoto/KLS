#!/usr/bin/env python3
"""Compare prior-work-style repeated factor and refactor averages.

The CKTSO and SubtreeLU released benchmarks reuse one unchanged values array
and report factorization and refactorization separately.  This scorer rejects
zero-repeat records so an initial/cold factor cannot be mislabeled as their
warm repeated-factor metric.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


METRICS = {
    "factor": ("factor_seconds_avg", "factor_repeat"),
    "refactor": ("refactor_seconds_avg", "refactor_repeat"),
}


def matrix_name(record: dict[str, object]) -> str:
    return Path(str(record.get("matrix", ""))).stem.lower()


def relative_residual(record: dict[str, object]) -> float | None:
    raw = record.get("relative_residual_l2", record.get("max_relative_residual"))
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def geometric_mean(values: list[float]) -> float:
    positive = [value for value in values if math.isfinite(value) and value > 0.0]
    if not positive:
        return math.nan
    return math.exp(sum(math.log(value) for value in positive) / len(positive))


def load(
    path: Path,
    residual_threshold: float,
    selection: str,
    minimum_repeats: dict[str, int],
) -> dict[str, dict[str, dict[str, float | int]]]:
    samples: dict[str, dict[str, list[dict[str, float | int]]]] = {
        metric: {} for metric in METRICS
    }
    with path.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"{path}:{line_number}: invalid JSON: {exc}") from exc
            if record.get("timeout_or_fail"):
                continue
            name = matrix_name(record)
            residual = relative_residual(record)
            if not name or residual is None or residual > residual_threshold:
                continue
            for metric, (timing_key, repeat_key) in METRICS.items():
                try:
                    timing = float(record[timing_key])
                    repeats = int(record[repeat_key])
                except (KeyError, TypeError, ValueError):
                    continue
                if (
                    repeats < minimum_repeats[metric]
                    or not math.isfinite(timing)
                    or timing <= 0.0
                ):
                    continue
                samples[metric].setdefault(name, []).append(
                    {"seconds": timing, "repeats": repeats}
                )

    selected: dict[str, dict[str, dict[str, float | int]]] = {}
    for metric, by_name in samples.items():
        selected[metric] = {}
        for name, values in by_name.items():
            ordered = sorted(values, key=lambda item: float(item["seconds"]))
            selected[metric][name] = (
                ordered[0] if selection == "min" else ordered[len(ordered) // 2]
            )
    return selected


def report(
    kls: dict[str, dict[str, dict[str, float | int]]],
    reference: dict[str, dict[str, dict[str, float | int]]],
    reference_name: str,
) -> None:
    for metric in METRICS:
        rows = []
        for name in kls[metric].keys() & reference[metric].keys():
            candidate = kls[metric][name]
            baseline = reference[metric][name]
            ratio = float(candidate["seconds"]) / float(baseline["seconds"])
            rows.append((ratio, name, candidate, baseline))
        rows.sort(reverse=True)
        print(
            f"\n== repeated {metric}: KLS vs {reference_name} "
            f"gm={geometric_mean([row[0] for row in rows]):.3f} "
            f"wins={sum(row[0] < 1.0 for row in rows)}/{len(rows)} =="
        )
        print(
            f"{'matrix':20} {'ratio':>7} {'KLS(s)':>11} "
            f"{reference_name + '(s)':>11} {'repeats K/R':>13}"
        )
        for ratio, name, candidate, baseline in rows:
            print(
                f"{name:20} {ratio:7.3f} "
                f"{float(candidate['seconds']):11.6g} "
                f"{float(baseline['seconds']):11.6g} "
                f"{int(candidate['repeats']):6d}/{int(baseline['repeats']):<6d}"
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kls", type=Path, help="KLS JSONL run")
    parser.add_argument("ck", type=Path, help="CKTSO JSONL run")
    parser.add_argument("st", type=Path, nargs="?", help="optional SubtreeLU JSONL run")
    parser.add_argument("--residual-threshold", type=float, default=1e-8)
    parser.add_argument("--selection", choices=("median", "min"), default="median")
    parser.add_argument(
        "--min-factor-repeat",
        type=int,
        default=1,
        help="reject records with fewer repeated factor calls (use 100 for paper parity)",
    )
    parser.add_argument(
        "--min-refactor-repeat",
        type=int,
        default=1,
        help="reject records with fewer refactor calls (use 100 for paper parity)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.min_factor_repeat < 1 or args.min_refactor_repeat < 1:
        raise SystemExit("minimum repeat counts must be positive")
    minimum_repeats = {
        "factor": args.min_factor_repeat,
        "refactor": args.min_refactor_repeat,
    }
    kls = load(args.kls, args.residual_threshold, args.selection, minimum_repeats)
    ck = load(args.ck, args.residual_threshold, args.selection, minimum_repeats)
    report(kls, ck, "CK")
    if args.st is not None:
        st = load(args.st, args.residual_threshold, args.selection, minimum_repeats)
        report(kls, st, "ST")


if __name__ == "__main__":
    main()
