#!/usr/bin/env python3
"""Compare two JSONL benchmark runs by matrix basename."""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import statistics


def spice_cycle_seconds(row: dict[str, object]) -> float:
    if "spice_cycle_seconds" in row:
        return float(row["spice_cycle_seconds"])
    analysis = float(row["analysis_seconds"])
    factor = float(row.get("initial_factor_seconds", row["factor_seconds_avg"]))
    refactor = float(row["refactor_seconds_avg"])
    solve = float(row["solve_seconds_avg"])
    return analysis + factor + solve + 99.0 * (refactor + solve)


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


def geometric_mean(values: list[float]) -> float:
    positives = [v for v in values if v > 0.0 and math.isfinite(v)]
    if not positives:
        return math.nan
    return math.exp(sum(math.log(v) for v in positives) / len(positives))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=pathlib.Path, required=True)
    parser.add_argument("--reference", type=pathlib.Path, required=True)
    parser.add_argument("--candidate-name", default="candidate")
    parser.add_argument("--reference-name", default="reference")
    parser.add_argument("--max-rows", type=int, default=12)
    args = parser.parse_args()

    candidate = load_rows(args.candidate)
    reference = load_rows(args.reference)
    common = sorted(set(candidate) & set(reference))
    if not common:
        raise SystemExit("no common matrix basenames")

    ratios: list[float] = []
    records: list[tuple[float, str, float, float]] = []
    for name in common:
        c = float(candidate[name]["spice_cycle_seconds"])
        r = float(reference[name]["spice_cycle_seconds"])
        if c <= 0.0 or r <= 0.0 or not math.isfinite(c) or not math.isfinite(r):
            continue
        ratio = c / r
        ratios.append(ratio)
        records.append((ratio, name, c, r))

    if not records:
        raise SystemExit("no finite positive common timings")

    wins = sum(1 for ratio, _, _, _ in records if ratio < 0.98)
    ties = sum(1 for ratio, _, _, _ in records if 0.98 <= ratio <= 1.02)
    losses = sum(1 for ratio, _, _, _ in records if ratio > 1.02)
    speedups = [1.0 / ratio for ratio in ratios]
    summary = {
        "candidate": args.candidate_name,
        "reference": args.reference_name,
        "matrices_common": len(records),
        "candidate_geomean_seconds": geometric_mean(
            [float(candidate[name]["spice_cycle_seconds"]) for _, name, _, _ in records]
        ),
        "reference_geomean_seconds": geometric_mean(
            [float(reference[name]["spice_cycle_seconds"]) for _, name, _, _ in records]
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
    for ratio, name, c, r in sorted(records)[: args.max_rows]:
        print(f"  {name}: ratio={ratio:.3f} {args.candidate_name}={c:.6g}s {args.reference_name}={r:.6g}s")

    print("\nLargest candidate losses:")
    for ratio, name, c, r in sorted(records, reverse=True)[: args.max_rows]:
        print(f"  {name}: ratio={ratio:.3f} {args.candidate_name}={c:.6g}s {args.reference_name}={r:.6g}s")

    missing_candidate = sorted(set(reference) - set(candidate))
    missing_reference = sorted(set(candidate) - set(reference))
    if missing_candidate:
        print(f"\nMissing from {args.candidate_name}: {', '.join(missing_candidate)}")
    if missing_reference:
        print(f"\nMissing from {args.reference_name}: {', '.join(missing_reference)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
