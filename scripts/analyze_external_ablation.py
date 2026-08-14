#!/usr/bin/env python3
"""Summarize post-reveal route ablations and AUTO regret."""

from __future__ import annotations

import argparse
import collections
import json
import pathlib
import statistics
import sys
from typing import Any

from analyze_external_validation import geometric_mean, horizon_seconds, validity
from run_external_ablation import ARMS


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign", type=pathlib.Path, required=True)
    parser.add_argument("--results", type=pathlib.Path, required=True)
    parser.add_argument("--horizon", type=int, default=100)
    parser.add_argument("--json", type=pathlib.Path)
    args = parser.parse_args()
    if args.horizon < 2:
        parser.error("--horizon must be at least two")
    try:
        campaign = json.loads(args.campaign.read_text(encoding="utf-8"))
        residual_limit = float(campaign["residual_limit"])
        grouped: dict[tuple[str, str], list[dict[str, Any]]] = collections.defaultdict(list)
        for line_number, raw in enumerate(
            args.results.read_text(encoding="utf-8").splitlines(), 1
        ):
            if not raw.strip():
                continue
            row = json.loads(raw)
            arm = str(row.get("validation_arm", ""))
            name = pathlib.Path(str(row.get("matrix", ""))).stem.casefold()
            if arm not in dict(ARMS) or not name:
                raise ValueError(f"line {line_number}: invalid arm or matrix")
            grouped[(name, arm)].append(row)
        medians: dict[tuple[str, str], float] = {}
        invalid: dict[str, list[str]] = collections.defaultdict(list)
        for key, rows in grouped.items():
            valid_times = [
                horizon_seconds(row, args.horizon)
                for row in rows if validity(row, residual_limit)[0]
            ]
            if len(valid_times) != len(rows):
                invalid[key[1]].append(key[0])
            elif valid_times:
                medians[key] = statistics.median(valid_times)
        names = sorted({name for name, _ in grouped})
        arms = [arm for arm, _ in ARMS]
        ratios: dict[str, list[float]] = collections.defaultdict(list)
        auto_regret: list[float] = []
        for name in names:
            auto = medians.get((name, "auto"))
            if auto is None:
                continue
            alternatives = [
                medians[(name, arm)] for arm in arms if (name, arm) in medians
            ]
            if alternatives:
                auto_regret.append(auto / min(alternatives))
            for arm in arms:
                timing = medians.get((name, arm))
                if timing is not None:
                    ratios[arm].append(timing / auto)
        summary = {
            "diagnostic_only": True,
            "horizon": args.horizon,
            "matrix_count": len(names),
            "auto_oracle_regret_geomean": geometric_mean(auto_regret),
            "auto_oracle_regret_max": max(auto_regret) if auto_regret else None,
            "arms": {
                arm: {
                    "valid_pairs_with_auto": len(ratios[arm]),
                    "arm_over_auto_geomean": geometric_mean(ratios[arm]),
                    "invalid_matrices": sorted(set(invalid[arm])),
                }
                for arm in arms
            },
        }
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(exc, file=sys.stderr)
        return 1
    rendered = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
