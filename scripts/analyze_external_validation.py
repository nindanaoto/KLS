#!/usr/bin/env python3
"""Strictly reduce a KLS external-validation campaign."""

from __future__ import annotations

import argparse
import collections
import json
import math
import pathlib
import re
import statistics
import sys
from typing import Any


SOLVERS = ("kls", "ck", "st", "klu")
CHALLENGERS = ("ck", "st", "klu")
HORIZONS = (10, 100, 1000)


def geometric_mean(values: list[float]) -> float | None:
    if not values or any(not math.isfinite(value) or value <= 0 for value in values):
        return None
    return math.exp(sum(math.log(value) for value in values) / len(values))


def manifest_names(path: pathlib.Path) -> list[str]:
    names: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        token = raw.split("#", 1)[0].strip()
        if token:
            names.append(pathlib.Path(token.split()[0]).stem.casefold())
    return names


def manifest_cohorts(path: pathlib.Path) -> dict[str, set[str]]:
    cohorts: dict[str, set[str]] = {
        "all": set(), "circuit_like": set(), "non_circuit": set()
    }
    for raw in path.read_text(encoding="utf-8").splitlines():
        token = raw.split("#", 1)[0].strip()
        if not token:
            continue
        name = pathlib.Path(token.split()[0]).stem.casefold()
        cohorts["all"].add(name)
        if re.search(r"\bcircuit_like=1\b", raw):
            cohorts["circuit_like"].add(name)
        else:
            cohorts["non_circuit"].add(name)
    return cohorts


def matrix_name(row: dict[str, Any]) -> str:
    return pathlib.Path(str(row.get("matrix", ""))).stem.casefold()


def finite(row: dict[str, Any], key: str) -> bool:
    try:
        return math.isfinite(float(row[key]))
    except (KeyError, TypeError, ValueError):
        return False


def validity(row: dict[str, Any], residual_limit: float) -> tuple[bool, str]:
    if row.get("status") not in (None, 0):
        return False, f"status={row.get('status')}"
    if row.get("verify_each_refactor") is not True:
        return False, "refactors-not-verified"
    for key in ("relative_residual_l2", "refactor_max_relative_residual"):
        if not finite(row, key):
            return False, f"missing-{key}"
        if float(row[key]) > residual_limit:
            return False, f"{key}>{residual_limit:g}"
    timing_keys = (
        "analysis_seconds", "solve_seconds_avg", "refactor_seconds_avg",
    )
    if any(not finite(row, key) or float(row[key]) < 0 for key in timing_keys):
        return False, "invalid-timing"
    return True, "valid"


def horizon_seconds(row: dict[str, Any], horizon: int) -> float:
    analysis = float(row["analysis_seconds"])
    factor = float(row.get("initial_factor_seconds", row["factor_seconds_avg"]))
    solve = float(row["solve_seconds_avg"])
    first_refactor = float(row.get("refactor_first_seconds", row["refactor_seconds_avg"]))
    steady_refactor = float(
        row.get("refactor_steady_seconds_avg", row["refactor_seconds_avg"])
    )
    first_paired_solve = float(row.get("refactor_solve_first_seconds", solve))
    steady_paired_solve = float(
        row.get("refactor_solve_steady_seconds_avg", solve)
    )
    return (
        analysis + factor + solve + first_refactor + first_paired_solve
        + max(horizon - 2, 0) * (steady_refactor + steady_paired_solve)
    )


def load_rows(path: pathlib.Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not path.exists():
        return rows
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip():
            continue
        try:
            row = json.loads(raw)
        except json.JSONDecodeError as exc:
            raise ValueError(f"{path}:{line_number}: {exc}") from exc
        if not isinstance(row, dict):
            raise ValueError(f"{path}:{line_number}: expected JSON object")
        rows.append(row)
    return rows


def reduce_solver(
    rows: list[dict[str, Any]], expected: list[str], passes: int, residual_limit: float
) -> tuple[dict[str, dict[str, Any]], dict[str, str]]:
    grouped: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for row in rows:
        grouped[matrix_name(row)].append(row)
    valid: dict[str, dict[str, Any]] = {}
    invalid: dict[str, str] = {}
    for name in expected:
        samples = grouped.get(name, [])
        if not samples:
            invalid[name] = "missing"
            continue
        verdicts = [validity(sample, residual_limit) for sample in samples]
        bad = [reason for ok, reason in verdicts if not ok]
        if bad:
            invalid[name] = bad[0]
            continue
        if len(samples) != passes:
            invalid[name] = f"incomplete-passes:{len(samples)}/{passes}"
            continue
        samples.sort(key=lambda row: horizon_seconds(row, 100))
        valid[name] = samples[len(samples) // 2]
    extras = sorted(set(grouped) - set(expected))
    if extras:
        raise ValueError("unexpected matrices in results: " + ", ".join(extras))
    return valid, invalid


def pairwise(
    kls: dict[str, dict[str, Any]], challenger: dict[str, dict[str, Any]], horizon: int,
    allowed: set[str] | None = None,
) -> dict[str, Any]:
    common = sorted(set(kls) & set(challenger))
    if allowed is not None:
        common = [name for name in common if name in allowed]
    ratios = [
        horizon_seconds(challenger[name], horizon) / horizon_seconds(kls[name], horizon)
        for name in common
    ]
    return {
        "valid_pairs": len(common),
        "kls_wins": sum(ratio > 1.02 for ratio in ratios),
        "ties_within_2pct": sum(0.98 <= ratio <= 1.02 for ratio in ratios),
        "kls_losses": sum(ratio < 0.98 for ratio in ratios),
        "challenger_over_kls_geomean": geometric_mean(ratios),
        "challenger_over_kls_median": statistics.median(ratios) if ratios else None,
        "minimum_ratio": min(ratios) if ratios else None,
        "minimum_matrix": common[ratios.index(min(ratios))] if ratios else None,
    }


def markdown(summary: dict[str, Any]) -> str:
    lines = [
        "# KLS external validation", "",
        f"Frozen commit: `{summary['campaign']['frozen_kls_commit']}`", "",
        "No selector or threshold changes are admissible after manifest reveal.", "",
    ]
    for ident, config in summary["configurations"].items():
        lines.extend([f"## {ident}", ""])
        lines.append(
            "| Solver | Valid | Invalid | H100 pairs | KLS wins | Ratio GM | Circuit GM |"
        )
        lines.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: |")
        kls_coverage = config["coverage"]["kls"]
        lines.append(
            f"| KLS | {kls_coverage['valid']} | {kls_coverage['invalid']} | — | — | — | — |"
        )
        for challenger in CHALLENGERS:
            coverage = config["coverage"][challenger]
            result = config["pairwise"][challenger]["100"]
            ratio = result["challenger_over_kls_geomean"]
            ratio_text = "—" if ratio is None else f"{ratio:.3f}"
            circuit = config["cohort_pairwise"]["circuit_like"][challenger]["100"]
            circuit_ratio = circuit["challenger_over_kls_geomean"]
            circuit_text = "—" if circuit_ratio is None else f"{circuit_ratio:.3f}"
            lines.append(
                f"| {challenger} | {coverage['valid']} | {coverage['invalid']} | "
                f"{result['valid_pairs']} | {result['kls_wins']} | {ratio_text} | "
                f"{circuit_text} |"
            )
        lines.append("")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=pathlib.Path, required=True)
    parser.add_argument("--json", type=pathlib.Path)
    parser.add_argument("--markdown", type=pathlib.Path)
    parser.add_argument("--allow-incomplete", action="store_true")
    args = parser.parse_args()
    results = args.results.resolve()
    try:
        campaign = json.loads((results / "campaign.json").read_text(encoding="utf-8"))
        manifest = pathlib.Path(campaign["manifest"])
        expected = manifest_names(manifest)
        cohorts = manifest_cohorts(manifest)
        configurations = json.loads(
            (results / "selected_configurations.json").read_text(encoding="utf-8")
        )
        summary: dict[str, Any] = {"campaign": campaign, "configurations": {}}
        incomplete: list[str] = []
        for ident in configurations:
            reduced: dict[str, dict[str, dict[str, Any]]] = {}
            invalid: dict[str, dict[str, str]] = {}
            for solver in SOLVERS:
                rows = load_rows(results / f"{ident}.{solver}.jsonl")
                reduced[solver], invalid[solver] = reduce_solver(
                    rows, expected, int(campaign["passes"]),
                    float(campaign["residual_limit"]),
                )
            coverage = {
                solver: {
                    "valid": len(reduced[solver]),
                    "invalid": len(invalid[solver]),
                    "invalid_matrices": invalid[solver],
                }
                for solver in SOLVERS
            }
            comparisons = {
                challenger: {
                    str(horizon): pairwise(
                        reduced["kls"], reduced[challenger], horizon
                    )
                    for horizon in HORIZONS
                }
                for challenger in CHALLENGERS
            }
            cohort_comparisons = {
                cohort: {
                    challenger: {
                        str(horizon): pairwise(
                            reduced["kls"], reduced[challenger], horizon, names
                        )
                        for horizon in HORIZONS
                    }
                    for challenger in CHALLENGERS
                }
                for cohort, names in cohorts.items()
            }
            summary["configurations"][ident] = {
                "coverage": coverage,
                "pairwise": comparisons,
                "cohort_pairwise": cohort_comparisons,
            }
            if not (results / f"{ident}.complete").exists():
                incomplete.append(ident)
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(exc, file=sys.stderr)
        return 1
    summary["incomplete_configurations"] = incomplete
    if args.json:
        args.json.write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    report = markdown(summary)
    if args.markdown:
        args.markdown.write_text(report, encoding="utf-8")
    else:
        print(report, end="")
    if incomplete and not args.allow_incomplete:
        print("incomplete configurations: " + ", ".join(incomplete), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
