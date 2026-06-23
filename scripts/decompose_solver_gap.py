#!/usr/bin/env python3
"""Decompose SPICE-cycle benchmark gaps by solver phase."""

from __future__ import annotations

import argparse
import json
import math
import pathlib


def field(row: dict[str, object], *names: str, default: float = 0.0) -> float:
    for name in names:
        if name in row:
            return float(row[name])
    return default


def components(row: dict[str, object]) -> dict[str, float]:
    analysis = field(row, "analysis_seconds")
    initial_factor = field(row, "initial_factor_seconds", "factor_seconds_avg")
    refactor = field(row, "refactor_seconds_avg", "refactor_seconds")
    solve = field(row, "solve_seconds_avg", "solve_seconds")
    return {
        "analysis": analysis,
        "initial_factor": initial_factor,
        "refactor_99": 99.0 * refactor,
        "solve_100": 100.0 * solve,
    }


def cycle_seconds(row: dict[str, object]) -> float:
    if "spice_cycle_seconds" in row:
        return float(row["spice_cycle_seconds"])
    return sum(components(row).values())


def load_rows(path: pathlib.Path) -> dict[str, dict[str, object]]:
    rows: dict[str, dict[str, object]] = {}
    with path.open("r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, 1):
            if not line.strip():
                continue
            row = json.loads(line)
            name = pathlib.Path(str(row["matrix"])).name
            if name in rows:
                raise ValueError(f"{path}:{line_no}: duplicate matrix basename {name}")
            row["spice_cycle_seconds"] = cycle_seconds(row)
            rows[name] = row
    return rows


def ratio(numerator: float, denominator: float) -> float:
    if denominator <= 0.0 or not math.isfinite(denominator):
        return math.nan
    return numerator / denominator


def fmt_ratio(value: float) -> str:
    if math.isnan(value):
        return "n/a"
    return f"{value:.2f}x"


def int_value(row: dict[str, object], name: str) -> int:
    if name not in row:
        return 0
    return int(row[name])


def float_value(row: dict[str, object], name: str) -> float:
    if name not in row:
        return 0.0
    return float(row[name])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=pathlib.Path, required=True)
    parser.add_argument("--reference", type=pathlib.Path, required=True)
    parser.add_argument("--candidate-name", default="candidate")
    parser.add_argument("--reference-name", default="reference")
    parser.add_argument("--max-rows", type=int, default=20)
    args = parser.parse_args()

    candidate = load_rows(args.candidate)
    reference = load_rows(args.reference)
    common = sorted(set(candidate) & set(reference))
    if not common:
        raise SystemExit("no common matrix basenames")

    rows: list[
        tuple[
            float,
            str,
            dict[str, float],
            dict[str, float],
            dict[str, object],
        ]
    ] = []
    for name in common:
        cand_cycle = cycle_seconds(candidate[name])
        ref_cycle = cycle_seconds(reference[name])
        if cand_cycle <= 0.0 or ref_cycle <= 0.0:
            continue
        rows.append(
            (
                cand_cycle / ref_cycle,
                name,
                components(candidate[name]),
                components(reference[name]),
                candidate[name],
            )
        )
    rows.sort(reverse=True)

    header = (
        "matrix,candidate_cycle,reference_cycle,cycle_ratio,"
        "dominant_candidate_phase,dominant_candidate_share,"
        "refactor_ratio,solve_ratio,initial_factor_ratio,analysis_ratio,"
        "n,nblocks,max_block,scale,offdiag_pivots,"
        "refactor_dependency_levels,refactor_dependency_max_width,"
        "refactor_dependency_edges,refactor_dependency_root_columns,"
        "refactor_dependency_leaf_columns,refactor_dependency_max_fanout,"
        "refactor_dependency_max_column_work,"
        "refactor_dependency_pipeline_max_column_work,"
        "refactor_supernode_candidate_count,"
        "refactor_supernode_candidate_rows,"
        "refactor_supernode_candidate_max_width,"
        "refactor_supernode_candidate_dense_entries,"
        "refactor_supernode_candidate_trailing_entries,"
        "row_refactor_group_count,row_refactor_group_level_count,"
        "row_refactor_group_level_max_width,"
        "row_refactor_segment_count,row_refactor_segment_rows,"
        "row_refactor_segment_max_width,"
        "row_refactor_segment_dense_entries,"
        "row_refactor_segment_trailing_entries,"
        "refactor_dependency_cluster_levels,"
        "refactor_dependency_pipeline_columns,refactor_dependency_work,"
        "refactor_dependency_pipeline_work"
    )
    print(header)
    for cycle_ratio, name, cand, ref, cand_row in rows[: args.max_rows]:
        cand_cycle = sum(cand.values())
        ref_cycle = sum(ref.values())
        dominant_phase = max(cand, key=cand.get)
        dominant_share = cand[dominant_phase] / cand_cycle if cand_cycle > 0.0 else math.nan
        print(
            f"{name},{cand_cycle:.6g},{ref_cycle:.6g},{cycle_ratio:.3f},"
            f"{dominant_phase},{dominant_share:.1%},"
            f"{fmt_ratio(ratio(cand['refactor_99'], ref['refactor_99']))},"
            f"{fmt_ratio(ratio(cand['solve_100'], ref['solve_100']))},"
            f"{fmt_ratio(ratio(cand['initial_factor'], ref['initial_factor']))},"
            f"{fmt_ratio(ratio(cand['analysis'], ref['analysis']))},"
            f"{int_value(cand_row, 'n')},"
            f"{int_value(cand_row, 'nblocks')},"
            f"{int_value(cand_row, 'max_block')},"
            f"{int_value(cand_row, 'scale')},"
            f"{int_value(cand_row, 'offdiag_pivots')},"
            f"{int_value(cand_row, 'refactor_dependency_levels')},"
            f"{int_value(cand_row, 'refactor_dependency_max_width')},"
            f"{int_value(cand_row, 'refactor_dependency_edges')},"
            f"{int_value(cand_row, 'refactor_dependency_root_columns')},"
            f"{int_value(cand_row, 'refactor_dependency_leaf_columns')},"
            f"{int_value(cand_row, 'refactor_dependency_max_fanout')},"
            f"{float_value(cand_row, 'refactor_dependency_max_column_work'):.6g},"
            f"{float_value(cand_row, 'refactor_dependency_pipeline_max_column_work'):.6g},"
            f"{int_value(cand_row, 'refactor_supernode_candidate_count')},"
            f"{int_value(cand_row, 'refactor_supernode_candidate_rows')},"
            f"{int_value(cand_row, 'refactor_supernode_candidate_max_width')},"
            f"{float_value(cand_row, 'refactor_supernode_candidate_dense_entries'):.6g},"
            f"{float_value(cand_row, 'refactor_supernode_candidate_trailing_entries'):.6g},"
            f"{int_value(cand_row, 'row_refactor_group_count')},"
            f"{int_value(cand_row, 'row_refactor_group_level_count')},"
            f"{int_value(cand_row, 'row_refactor_group_level_max_width')},"
            f"{int_value(cand_row, 'row_refactor_segment_count')},"
            f"{int_value(cand_row, 'row_refactor_segment_rows')},"
            f"{int_value(cand_row, 'row_refactor_segment_max_width')},"
            f"{float_value(cand_row, 'row_refactor_segment_dense_entries'):.6g},"
            f"{float_value(cand_row, 'row_refactor_segment_trailing_entries'):.6g},"
            f"{int_value(cand_row, 'refactor_dependency_cluster_levels')},"
            f"{int_value(cand_row, 'refactor_dependency_pipeline_columns')},"
            f"{float_value(cand_row, 'refactor_dependency_work'):.6g},"
            f"{float_value(cand_row, 'refactor_dependency_pipeline_work'):.6g}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
