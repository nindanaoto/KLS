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


def str_value(row: dict[str, object], name: str) -> str:
    value = row.get(name, "")
    return str(value) if value is not None else ""


def share(numerator: float, denominator: float) -> float:
    if denominator <= 0.0 or not math.isfinite(denominator):
        return math.nan
    return numerator / denominator


def fmt_share(value: float) -> str:
    if math.isnan(value):
        return "n/a"
    return f"{value:.3f}"


def row_refactor_panel_overstage_entries(row: dict[str, object]) -> float:
    return max(
        float_value(row, "row_refactor_compact_supernode_update_entries"),
        float_value(row, "row_refactor_native_row_panel_blocked_entries"),
        float_value(row, "row_refactor_compact_supernode_batch_entries"),
    )


def egraph_scalar_numeric_owner_missing(row: dict[str, object]) -> bool:
    egraph_work = float_value(row, "refactor_dependency_work")
    pipeline_work = float_value(row, "refactor_dependency_pipeline_work")
    if egraph_work < 5.0e7 or share(pipeline_work, egraph_work) < 0.50:
        return False
    active_grouped_updates = (
        int_value(row, "refactor_last_supernode_update_runs")
        + int_value(row, "refactor_last_supernode_cached_probe_applied")
        + int_value(row, "refactor_last_u_supernode_l_update_runs")
        + int_value(row, "refactor_last_supernode_algorithm5_payoff_prefix_prep_runs")
        + int_value(row, "refactor_last_supernode_algorithm5_payoff_current_state_seed_runs")
        + int_value(row, "refactor_last_supernode_algorithm5_payoff_suffix_advance_deps")
        + int_value(row, "refactor_last_btf_scalar_run_exec_runs")
    )
    return active_grouped_updates == 0


def exact_shape_bounded_owner_payoff_tiny(row: dict[str, object]) -> bool:
    exact_rows = int_value(
        row, "refactor_supernode_consumer_plan_group_l_batch_candidate_run_rows"
    )
    payoff_rows = int_value(
        row,
        "refactor_supernode_consumer_plan_shape_bounded_advance_payoff_run_rows",
    )
    bounded_rows = int_value(
        row, "refactor_supernode_consumer_plan_shape_bounded_advance_run_rows"
    )
    if exact_rows < 100_000 or bounded_rows == 0:
        return False
    return share(float(payoff_rows), float(exact_rows)) < 0.01


def paper_gap_signal(
    cand: dict[str, float],
    cand_row: dict[str, object],
) -> str:
    dominant_phase = max(cand, key=cand.get)
    initial_path = str_value(cand_row, "initial_factor_path")
    last_path = str_value(cand_row, "last_factor_path")
    last_refactor_path = str_value(cand_row, "last_refactor_path")
    row_groups = int_value(cand_row, "row_refactor_group_count")
    row_run = int_value(cand_row, "row_refactor_last_run")
    first_skip_scaled_single = int_value(
        cand_row, "kls_first_auto_skipped_scaled_single_block_count"
    )
    egraph_work = float_value(cand_row, "refactor_dependency_work")
    row_work = float_value(cand_row, "row_refactor_total_group_work")
    compact_work = float_value(cand_row, "row_refactor_compact_dense_panel_update_work")
    native_panel = int_value(cand_row, "row_refactor_last_native_row_panel")
    native_fallbacks = int_value(
        cand_row, "row_refactor_native_row_panel_fallback_count"
    )
    native_rejects = int_value(
        cand_row, "row_refactor_native_row_panel_checked_reject_count"
    )
    auto_lower_rejected = int_value(
        cand_row, "row_refactor_auto_lower_bound_rejected"
    )
    panel_entries = row_refactor_panel_overstage_entries(cand_row)
    separator_flop_queue = int_value(
        cand_row, "row_refactor_last_separator_flop_queue"
    )
    separator_private_groups = int_value(
        cand_row, "row_refactor_last_separator_flop_private_groups"
    )
    separator_pipeline_groups = int_value(
        cand_row, "row_refactor_last_separator_flop_pipeline_groups"
    )
    compact_update_rows = int_value(
        cand_row, "row_refactor_compact_supernode_update_rows"
    )
    compact_eligible_rows = int_value(
        cand_row, "row_refactor_compact_dense_panel_eligible_rows"
    )
    btf_scalar_run_rows = int_value(
        cand_row, "refactor_last_btf_scalar_run_rows"
    )

    if first_skip_scaled_single:
        return "missing_parallel_rowup_first_factor"
    if dominant_phase == "refactor_99":
        if last_refactor_path == "row_refactor" or row_run:
            if native_rejects:
                return "native_row_panel_checked_reject"
            if native_fallbacks:
                return "native_row_panel_fallback"
            if (
                panel_entries > 0.0
                and row_work > 0.0
                and panel_entries > 1.5 * row_work
            ):
                return "row_refactor_panel_overstaged"
            if separator_flop_queue and separator_private_groups > 0:
                if compact_update_rows == 0 and compact_eligible_rows > 0:
                    return "row_refactor_private_scalar_scaffold"
                if separator_pipeline_groups > 0:
                    return "row_refactor_many_private_groups"
            if native_panel:
                return "native_row_panel_active"
            if egraph_work > 0.0 and row_work > egraph_work:
                return "row_kernel_more_work_than_egraph"
            if compact_work > 0.0:
                return "row_panel_kernel_active"
            return "row_kernel_active"
        if last_refactor_path == "egraph":
            if egraph_scalar_numeric_owner_missing(cand_row):
                if exact_shape_bounded_owner_payoff_tiny(cand_row):
                    return "egraph_scalar_bounded_owner_payoff_tiny"
                if int_value(cand_row, "refactor_last_btf_scalar_run_exec_runs") > 0:
                    return "egraph_btf_scalar_run_executor_active"
                if auto_lower_rejected:
                    if int_value(
                        cand_row,
                        "refactor_btf_scalar_run_group_reused_entries",
                    ) > 0:
                        return "egraph_scalar_tail_grouped_producer_runs_unowned"
                    if btf_scalar_run_rows > 0:
                        return "egraph_scalar_tail_producer_runs_unowned"
                    return "egraph_scalar_tail_numeric_owner_missing"
                return "egraph_scalar_pipeline_numeric_owner_missing"
            if auto_lower_rejected:
                return "row_refactor_lower_bound_rejected"
            if egraph_work > 0.0:
                return "column_egraph_refactor_missing_row_engine"
        if last_path == "kls_fast_refactor":
            if row_run:
                if native_rejects:
                    return "native_row_panel_checked_reject"
                if native_fallbacks:
                    return "native_row_panel_fallback"
                if native_panel:
                    return "native_row_panel_active"
                if egraph_work > 0.0 and row_work > egraph_work:
                    return "row_kernel_more_work_than_egraph"
                if compact_work > 0.0:
                    return "row_panel_kernel_active"
                return "row_kernel_active"
            if egraph_work > 0.0:
                return "column_egraph_refactor_missing_row_engine"
            return "fast_refactor_without_row_diagnostics"
        return "refactor_dominant_non_kls_fast_path"
    if dominant_phase == "initial_factor" and "klu" in initial_path:
        return "klu_first_factor_missing_row_engine"
    if row_groups > 0 and not row_run:
        return "row_metadata_built_but_not_selected"
    return "mixed"


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
        "candidate_initial_factor_path,candidate_last_factor_path,"
        "candidate_last_refactor_path,"
        "kls_tail_last_mapped_columns,kls_tail_mapped_column_count,"
        "kls_first_last_row_uplooking_columns,"
        "kls_first_row_uplooking_column_count,"
        "kls_first_last_row_refactor_seeded_rows,"
        "kls_first_row_refactor_seeded_row_count,"
        "kls_first_last_row_pipeline,"
        "kls_first_row_pipeline_run_count,"
        "kls_first_last_row_pipeline_rows,"
        "kls_first_last_row_pipeline_threads,"
        "kls_first_last_row_pipeline_partial,"
        "kls_first_row_pipeline_partial_run_count,"
        "kls_first_last_row_pipeline_partial_rows,"
        "kls_first_last_row_pipeline_partial_threads,"
        "kls_first_last_row_pipeline_pivot_tail,"
        "kls_first_row_pipeline_pivot_tail_run_count,"
        "kls_first_last_row_pipeline_pivot_tail_rows,"
        "kls_first_last_row_pipeline_pivot_restarts,"
        "kls_first_row_pipeline_pivot_restart_count,"
        "kls_first_last_row_pipeline_pivot_serial_rows,"
        "kls_first_last_row_pipeline_prefix_panel_rebuild,"
        "kls_first_row_pipeline_prefix_panel_rebuild_count,"
        "kls_first_last_row_pipeline_prefix_panel_rebuild_rows,"
        "kls_first_last_row_supernode_update,"
        "kls_first_row_supernode_update_run_count,"
        "kls_first_last_row_supernode_update_groups,"
        "kls_first_last_row_supernode_update_rows,"
        "kls_first_last_row_supernode_panel_update,"
        "kls_first_row_supernode_panel_update_run_count,"
        "kls_first_last_row_supernode_panel_update_groups,"
        "kls_first_last_row_supernode_panel_update_rows,"
        "kls_first_last_dynamic_column_pivots,"
        "kls_first_dynamic_column_pivot_count,"
        "kls_first_last_separator_dynamic_column_pivots,"
        "kls_first_separator_dynamic_column_pivot_count,"
        "kls_first_last_separator_extent_dynamic_column_pivots,"
        "kls_first_separator_extent_dynamic_column_pivot_count,"
        "kls_first_last_separator_dynamic_column_fallbacks,"
        "kls_first_separator_dynamic_column_fallback_count,"
        "kls_first_last_separator_dynamic_column_rejects,"
        "kls_first_separator_dynamic_column_reject_count,"
        "kls_first_last_separator_queue,"
        "kls_first_separator_queue_run_count,"
        "kls_first_last_separator_queue_private_components,"
        "kls_first_last_separator_queue_pipeline_components,"
        "kls_first_last_separator_queue_private_rows,"
        "kls_first_last_separator_queue_pipeline_rows,"
        "kls_first_last_separator_queue_nonempty_threads,"
        "kls_first_last_separator_queue_max_thread_rows,"
        "kls_first_last_separator_queue_partitioned,"
        "kls_first_separator_queue_partitioned_count,"
        "kls_first_last_separator_queue_split_components,"
        "kls_first_last_separator_queue_executed,"
        "kls_first_separator_queue_executed_run_count,"
        "kls_first_last_separator_queue_executed_private_rows,"
        "kls_first_last_separator_queue_executed_pipeline_rows,"
        "kls_first_last_separator_queue_parallel_private,"
        "kls_first_separator_queue_parallel_private_run_count,"
        "kls_first_last_separator_queue_parallel_private_rows,"
        "kls_first_last_separator_queue_parallel_private_threads,"
        "kls_first_last_separator_queue_parallel_pipeline,"
        "kls_first_separator_queue_parallel_pipeline_run_count,"
        "kls_first_last_separator_queue_parallel_pipeline_rows,"
        "kls_first_last_separator_queue_parallel_pipeline_threads,"
        "kls_first_last_separator_queue_pipeline_partial,"
        "kls_first_separator_queue_pipeline_partial_run_count,"
        "kls_first_last_separator_queue_pipeline_partial_rows,"
        "kls_first_last_separator_queue_pipeline_partial_threads,"
        "kls_first_last_separator_queue_pipeline_wait_partial,"
        "kls_first_separator_queue_pipeline_wait_partial_run_count,"
        "kls_first_last_separator_queue_pipeline_wait_partial_rows,"
        "kls_first_last_separator_queue_pipeline_wait_partial_deps,"
        "kls_first_last_separator_queue_pipeline_supernode_update,"
        "kls_first_separator_queue_pipeline_supernode_update_run_count,"
        "kls_first_last_separator_queue_pipeline_supernode_update_groups,"
        "kls_first_last_separator_queue_pipeline_supernode_update_rows,"
        "kls_first_last_separator_queue_pipeline_supernode_panel_update,"
        "kls_first_separator_queue_pipeline_supernode_panel_update_run_count,"
        "kls_first_last_separator_queue_pipeline_supernode_panel_update_groups,"
        "kls_first_last_separator_queue_pipeline_supernode_panel_update_rows,"
        "kls_first_last_separator_queue_pipeline_pivot_tail,"
        "kls_first_separator_queue_pipeline_pivot_tail_run_count,"
        "kls_first_last_separator_queue_pipeline_pivot_tail_rows,"
        "kls_first_last_separator_queue_pipeline_prefix_panel_rebuild,"
        "kls_first_separator_queue_pipeline_prefix_panel_rebuild_count,"
        "kls_first_last_separator_queue_pipeline_prefix_panel_rebuild_rows,"
        "kls_first_auto_skipped_scaled_single_block,"
        "kls_first_auto_skipped_scaled_single_block_count,"
        "paper_gap_signal,refactor_dependency_pipeline_share,"
        "parallel_model_r1,parallel_model_r2,"
        "parallel_model_recommends_parallel,"
        "parallel_task_flow_threads,parallel_task_flow_dependencies,"
        "parallel_task_flow_work,parallel_task_flow_finish_time,"
        "parallel_task_flow_speedup,"
        "parallel_task_flow_recommends_parallel,"
        "row_refactor_group_work_ratio,"
        "row_refactor_compact_panel_work_share,"
        "row_refactor_panel_overstage_ratio,"
        "row_refactor_private_group_share,"
        "refactor_ratio,solve_ratio,initial_factor_ratio,analysis_ratio,"
        "n,nblocks,max_block,scale,offdiag_pivots,"
        "fast_repaired_tail_restart_overcompute_columns,"
        "fast_repaired_tail_restart_overcompute_work,"
        "fast_repaired_tail_restart_exact_mask,"
        "fast_repaired_tail_restart_etree_mask,"
        "fast_kls_rebuild_restarts,"
        "fast_kls_block_restart_last_row_pipeline,"
        "fast_kls_block_restart_row_pipeline_count,"
        "fast_kls_block_restart_last_row_pipeline_rows,"
        "fast_kls_block_restart_last_row_pipeline_threads,"
        "fast_kls_block_restart_last_row_pipeline_prefix_rows,"
        "fast_kls_block_restart_last_row_pipeline_suffix_rows,"
        "fast_kls_block_restart_last_row_pipeline_gap_rows,"
        "fast_kls_block_restart_last_row_pipeline_etree_tail,"
        "fast_kls_block_restart_row_pipeline_etree_tail_count,"
        "fast_kls_block_restart_last_row_pipeline_etree_tail_rows,"
        "fast_kls_block_restart_last_row_pipeline_etree_tail_gap_rows,"
        "fast_kls_block_restart_last_row_pipeline_etree_tail_exact_mask,"
        "fast_kls_block_restart_last_row_pipeline_separator_tail_scope,"
        "fast_kls_block_restart_row_pipeline_separator_tail_scope_count,"
        "fast_kls_block_restart_last_row_pipeline_separator_tail_scope_rows,"
        "fast_kls_block_restart_last_row_pipeline_separator_queue,"
        "fast_kls_block_restart_row_pipeline_separator_queue_count,"
        "fast_kls_block_restart_last_row_pipeline_separator_private_rows,"
        "fast_kls_block_restart_last_row_pipeline_separator_pipeline_rows,"
        "fast_kls_block_restart_last_row_pipeline_separator_private_threads,"
        "fast_kls_block_restart_last_row_pipeline_separator_partitioned,"
        "fast_kls_block_restart_last_row_pipeline_separator_split_components,"
        "fast_kls_block_restart_last_row_pipeline_pivot_tail_rows,"
        "fast_kls_block_restart_last_row_pipeline_pivot_restarts,"
        "fast_kls_block_restart_last_row_pipeline_supernode_update_groups,"
        "fast_kls_block_restart_last_row_pipeline_supernode_update_rows,"
        "fast_kls_block_restart_last_row_pipeline_supernode_panel_update_groups,"
        "fast_kls_block_restart_last_row_pipeline_supernode_panel_update_rows,"
        "fast_rejected_pivoting_tail_contiguous,"
        "fast_rejected_pivoting_tail_suffix_exact,"
        "fast_rejected_pivoting_tail_gap_columns,"
        "fast_rejected_pivoting_tail_block_seed_columns,"
        "fast_rejected_pivoting_tail_suffix_overcompute_columns,"
        "fast_rejected_pivoting_tail_suffix_overcompute_work,"
        "fast_rejected_pivoting_tail_etree_edges,"
        "fast_rejected_pivoting_tail_etree_roots,"
        "fast_rejected_pivoting_tail_etree_leaves,"
        "fast_rejected_pivoting_tail_etree_max_fanout,"
        "fast_rejected_pivoting_tail_etree_levels,"
        "fast_rejected_pivoting_tail_etree_max_width,"
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
        "refactor_supernode_panel_count,"
        "refactor_supernode_panel_used_count,"
        "refactor_last_supernode_pipeline_tasks,"
        "refactor_last_supernode_pipeline_columns,"
        "refactor_supernode_pipeline_task_count,"
        "refactor_supernode_pipeline_column_count,"
        "refactor_last_supernode_update_runs,"
        "refactor_last_supernode_update_rows,"
        "refactor_last_supernode_update_entries,"
        "refactor_supernode_update_run_count,"
        "refactor_supernode_update_rows,"
        "refactor_supernode_update_entries,"
        "refactor_last_supernode_cblas_update_runs,"
        "refactor_last_supernode_cblas_update_rows,"
        "refactor_last_supernode_cblas_update_entries,"
        "refactor_supernode_cblas_update_run_count,"
        "refactor_supernode_cblas_update_rows,"
        "refactor_supernode_cblas_update_entries,"
        "refactor_last_supernode_blocked_update_runs,"
        "refactor_last_supernode_blocked_update_rows,"
        "refactor_last_supernode_blocked_update_entries,"
        "refactor_supernode_blocked_update_run_count,"
        "refactor_supernode_blocked_update_rows,"
        "refactor_supernode_blocked_update_entries,"
        "refactor_supernode_cached_probe_attempt_count,"
        "refactor_supernode_cached_probe_panel_hits,"
        "refactor_supernode_cached_probe_contiguous,"
        "refactor_supernode_cached_probe_allowed,"
        "refactor_supernode_cached_probe_allowed_rows,"
        "refactor_supernode_cached_probe_applied,"
        "refactor_supernode_cached_probe_applied_rows,"
        "refactor_supernode_cached_probe_shape_rejects,"
        "refactor_supernode_cached_probe_shape_reject_rows,"
        "refactor_supernode_cached_probe_stream_rejects,"
        "refactor_supernode_cached_probe_stream_reject_rows,"
        "refactor_supernode_cached_probe_work_rejects,"
        "refactor_supernode_cached_probe_work_reject_rows,"
        "refactor_supernode_cached_probe_workspace_rejects,"
        "refactor_supernode_cached_probe_workspace_reject_rows,"
        "refactor_last_btf_scalar_run_candidates,"
        "refactor_last_btf_scalar_run_rows,"
        "refactor_last_btf_scalar_run_entries,"
        "refactor_last_btf_scalar_run_max_rows,"
        "refactor_btf_scalar_run_candidate_count,"
        "refactor_btf_scalar_run_rows,"
        "refactor_btf_scalar_run_entries,"
        "refactor_btf_scalar_run_max_rows,"
        "refactor_last_btf_scalar_run_exec_runs,"
        "refactor_last_btf_scalar_run_exec_rows,"
        "refactor_last_btf_scalar_run_exec_entries,"
        "refactor_last_btf_scalar_run_exec_max_rows,"
        "refactor_btf_scalar_run_exec_count,"
        "refactor_btf_scalar_run_exec_rows,"
        "refactor_btf_scalar_run_exec_entries,"
        "refactor_btf_scalar_run_exec_max_rows,"
        "refactor_btf_scalar_run_group_built,"
        "refactor_btf_scalar_run_group_count,"
        "refactor_btf_scalar_run_group_current_total,"
        "refactor_btf_scalar_run_group_multi_count,"
        "refactor_btf_scalar_run_group_multi_current_total,"
        "refactor_btf_scalar_run_group_rows,"
        "refactor_btf_scalar_run_group_reused_rows,"
        "refactor_btf_scalar_run_group_entries,"
        "refactor_btf_scalar_run_group_reused_entries,"
        "refactor_btf_scalar_run_group_max_currents,"
        "refactor_btf_scalar_run_group_max_rows,"
        "refactor_last_btf_scalar_run_group_waits,"
        "refactor_last_btf_scalar_run_group_overlap_waits,"
        "refactor_last_btf_scalar_run_group_max_live,"
        "refactor_last_ready_queue_columns,"
        "refactor_ready_queue_run_count,"
        "row_refactor_group_count,"
        "row_refactor_group_single_count,"
        "row_refactor_group_batch_count,"
        "row_refactor_group_batch_rows,"
        "row_refactor_group_scalar_candidate_count,"
        "row_refactor_group_scalar_candidate_rows,"
        "row_refactor_group_scalar_short_count,"
        "row_refactor_group_scalar_short_rows,"
        "row_refactor_group_scalar_stop_level_mismatch_count,"
        "row_refactor_group_scalar_stop_internal_dep_count,"
        "row_refactor_group_scalar_stop_next_segment_count,"
        "row_refactor_group_scalar_stop_max_width_count,"
        "row_refactor_group_scalar_stop_matrix_end_count,"
        "row_refactor_group_generic_count,"
        "row_refactor_group_dense_count,"
        "row_refactor_group_level_count,"
        "row_refactor_group_level_max_width,"
        "row_refactor_group_cluster_levels,"
        "row_refactor_group_pipeline_groups,"
        "row_refactor_group_pipeline_rows,"
        "row_refactor_group_pipeline_work,"
        "row_refactor_total_group_work,"
        "row_refactor_auto_enabled,"
        "row_refactor_auto_values_ready,"
        "row_refactor_auto_work_allowed,"
        "row_refactor_auto_should_run,"
        "row_refactor_auto_lower_bound_work,"
        "row_refactor_auto_lower_bound_rejected,"
        "row_refactor_auto_pattern_build_failed,"
        "row_refactor_auto_value_copy_failed,"
        "row_refactor_auto_model_recommended,"
        "row_refactor_auto_model_attempted,"
        "row_refactor_auto_model_accepted,"
        "row_solve_parallel_run_count,row_solve_parallel_l_slice_runs,"
        "row_solve_parallel_u_slice_runs,"
        "row_solve_parallel_l_sparse_level_runs,"
        "row_solve_parallel_u_sparse_level_runs,"
        "row_solve_thread_count,row_solve_l_thread_max_rect_entries,"
        "row_solve_u_thread_max_rect_entries,"
        "row_solve_partition_ready,row_solve_partition_slices,"
        "row_solve_l_sparse_level_count,row_solve_l_sparse_cluster_levels,"
        "row_solve_l_sparse_level_max_width,"
        "row_solve_l_dense_tail_start,row_solve_l_dense_tail_rows,"
        "row_solve_l_dense_tail_entries,row_solve_l_slice_max_entries,"
        "row_solve_l_segmented_rows,row_solve_l_rect_entries,"
        "row_solve_l_tri_entries,"
        "row_solve_u_sparse_level_count,row_solve_u_sparse_cluster_levels,"
        "row_solve_u_sparse_level_max_width,"
        "row_solve_u_dense_tail_start,"
        "row_solve_u_dense_tail_rows,row_solve_u_dense_tail_entries,"
        "row_solve_u_slice_max_entries,"
        "row_solve_u_segmented_rows,row_solve_u_rect_entries,"
        "row_solve_u_tri_entries,"
        "row_refactor_last_run,"
        "row_refactor_last_checked,"
        "row_refactor_last_parallel,"
        "row_refactor_run_count,"
        "row_refactor_checked_run_count,"
        "row_refactor_parallel_run_count,"
        "row_refactor_segment_count,row_refactor_segment_rows,"
        "row_refactor_segment_max_width,"
        "row_refactor_segment_dense_entries,"
        "row_refactor_segment_trailing_entries,"
        "row_refactor_dense_segment_count,row_refactor_dense_segment_rows,"
        "row_refactor_dense_segment_max_width,"
        "row_refactor_dense_segment_dense_entries,"
        "row_refactor_dense_segment_trailing_entries,"
        "row_refactor_compact_dense_panel_eligible_count,"
        "row_refactor_compact_dense_panel_eligible_rows,"
        "row_refactor_compact_dense_panel_update_work,"
        "row_refactor_compact_dense_panel_entries,"
        "row_refactor_compact_dense_panel_persistent_groups,"
        "row_refactor_compact_dense_panel_persistent_entries,"
        "row_refactor_last_compact_dense_panel_persistent,"
        "row_refactor_compact_dense_panel_persistent_run_count,"
        "row_refactor_last_compact_dense_panel,"
        "row_refactor_compact_dense_panel_count,"
        "row_refactor_last_compact_dense_panel_blocked,"
        "row_refactor_compact_dense_panel_blocked_run_count,"
        "row_refactor_compact_dense_panel_blocked_rows,"
        "row_refactor_compact_dense_panel_blocked_entries,"
        "row_refactor_native_row_panel_enabled,"
        "row_refactor_last_native_row_panel,"
        "row_refactor_native_row_panel_count,"
        "row_refactor_native_row_panel_rows,"
        "row_refactor_native_row_panel_entries,"
        "row_refactor_native_row_panel_blocked_count,"
        "row_refactor_native_row_panel_blocked_rows,"
        "row_refactor_native_row_panel_blocked_entries,"
        "row_refactor_native_row_panel_fallback_count,"
        "row_refactor_native_row_panel_checked_reject_count,"
        "row_refactor_last_compact_supernode_partial_update,"
        "row_refactor_compact_supernode_partial_update_count,"
        "row_refactor_compact_supernode_partial_update_rows,"
        "row_refactor_compact_supernode_partial_update_entries,"
        "row_refactor_last_local_ready_groups,"
        "row_refactor_local_ready_group_count,"
        "row_refactor_last_private_ready_groups,"
        "row_refactor_private_ready_group_count,"
        "row_refactor_last_separator_flop_queue,"
        "row_refactor_separator_flop_queue_run_count,"
        "row_refactor_last_separator_flop_ordered_private,"
        "row_refactor_separator_flop_ordered_private_run_count,"
        "row_refactor_last_separator_flop_components,"
        "row_refactor_separator_flop_component_count,"
        "row_refactor_last_separator_flop_private_groups,"
        "row_refactor_last_separator_flop_pipeline_groups,"
        "row_refactor_last_separator_flop_closure_groups,"
        "row_refactor_last_separator_flop_private_threads,"
        "row_refactor_last_separator_flop_private_min_groups,"
        "row_refactor_last_separator_flop_private_max_groups,"
        "row_refactor_last_separator_flop_private_min_work,"
        "row_refactor_last_separator_flop_private_max_work,"
        "row_refactor_separator_flop_private_group_count,"
        "row_refactor_separator_flop_pipeline_group_count,"
        "row_refactor_separator_flop_closure_group_count,"
        "refactor_dependency_cluster_levels,"
        "refactor_dependency_pipeline_columns,refactor_dependency_work,"
        "refactor_dependency_pipeline_work,"
        "refactor_last_btf_scalar_run_candidates,"
        "refactor_last_btf_scalar_run_rows,"
        "refactor_last_btf_scalar_run_entries,"
        "refactor_last_btf_scalar_run_max_rows,"
        "refactor_last_btf_scalar_run_exec_runs,"
        "refactor_last_btf_scalar_run_exec_rows,"
        "refactor_last_btf_scalar_run_exec_entries,"
        "refactor_last_btf_scalar_run_exec_max_rows,"
        "refactor_btf_scalar_run_group_multi_count,"
        "refactor_btf_scalar_run_group_multi_current_total,"
        "refactor_btf_scalar_run_group_reused_entries,"
        "refactor_btf_scalar_run_group_max_currents,"
        "refactor_last_btf_scalar_run_group_waits,"
        "refactor_last_btf_scalar_run_group_overlap_waits,"
        "refactor_last_btf_scalar_run_group_max_live"
    )
    print(header)
    for cycle_ratio, name, cand, ref, cand_row in rows[: args.max_rows]:
        cand_cycle = sum(cand.values())
        ref_cycle = sum(ref.values())
        dominant_phase = max(cand, key=cand.get)
        dominant_share = cand[dominant_phase] / cand_cycle if cand_cycle > 0.0 else math.nan
        egraph_work = float_value(cand_row, "refactor_dependency_work")
        row_work = float_value(cand_row, "row_refactor_total_group_work")
        compact_work = float_value(
            cand_row, "row_refactor_compact_dense_panel_update_work"
        )
        row_groups = float_value(cand_row, "row_refactor_group_count")
        private_groups = float_value(
            cand_row, "row_refactor_last_separator_flop_private_groups"
        )
        panel_entries = row_refactor_panel_overstage_entries(cand_row)
        print(
            f"{name},{cand_cycle:.6g},{ref_cycle:.6g},{cycle_ratio:.3f},"
            f"{dominant_phase},{dominant_share:.1%},"
            f"{str_value(cand_row, 'initial_factor_path')},"
            f"{str_value(cand_row, 'last_factor_path')},"
            f"{str_value(cand_row, 'last_refactor_path')},"
            f"{int_value(cand_row, 'kls_tail_last_mapped_columns')},"
            f"{int_value(cand_row, 'kls_tail_mapped_column_count')},"
            f"{int_value(cand_row, 'kls_first_last_row_uplooking_columns')},"
            f"{int_value(cand_row, 'kls_first_row_uplooking_column_count')},"
            f"{int_value(cand_row, 'kls_first_last_row_refactor_seeded_rows')},"
            f"{int_value(cand_row, 'kls_first_row_refactor_seeded_row_count')},"
            f"{int_value(cand_row, 'kls_first_last_row_pipeline')},"
            f"{int_value(cand_row, 'kls_first_row_pipeline_run_count')},"
            f"{int_value(cand_row, 'kls_first_last_row_pipeline_rows')},"
            f"{int_value(cand_row, 'kls_first_last_row_pipeline_threads')},"
            f"{int_value(cand_row, 'kls_first_last_row_pipeline_partial')},"
            f"{int_value(cand_row, 'kls_first_row_pipeline_partial_run_count')},"
            f"{int_value(cand_row, 'kls_first_last_row_pipeline_partial_rows')},"
            f"{int_value(cand_row, 'kls_first_last_row_pipeline_partial_threads')},"
            f"{int_value(cand_row, 'kls_first_last_row_pipeline_pivot_tail')},"
            f"{int_value(cand_row, 'kls_first_row_pipeline_pivot_tail_run_count')},"
            f"{int_value(cand_row, 'kls_first_last_row_pipeline_pivot_tail_rows')},"
            f"{int_value(cand_row, 'kls_first_last_row_pipeline_pivot_restarts')},"
            f"{int_value(cand_row, 'kls_first_row_pipeline_pivot_restart_count')},"
            f"{int_value(cand_row, 'kls_first_last_row_pipeline_pivot_serial_rows')},"
            f"{int_value(cand_row, 'kls_first_last_row_pipeline_prefix_panel_rebuild')},"
            f"{int_value(cand_row, 'kls_first_row_pipeline_prefix_panel_rebuild_count')},"
            f"{int_value(cand_row, 'kls_first_last_row_pipeline_prefix_panel_rebuild_rows')},"
            f"{int_value(cand_row, 'kls_first_last_row_supernode_update')},"
            f"{int_value(cand_row, 'kls_first_row_supernode_update_run_count')},"
            f"{int_value(cand_row, 'kls_first_last_row_supernode_update_groups')},"
            f"{int_value(cand_row, 'kls_first_last_row_supernode_update_rows')},"
            f"{int_value(cand_row, 'kls_first_last_row_supernode_panel_update')},"
            f"{int_value(cand_row, 'kls_first_row_supernode_panel_update_run_count')},"
            f"{int_value(cand_row, 'kls_first_last_row_supernode_panel_update_groups')},"
            f"{int_value(cand_row, 'kls_first_last_row_supernode_panel_update_rows')},"
            f"{int_value(cand_row, 'kls_first_last_dynamic_column_pivots')},"
            f"{int_value(cand_row, 'kls_first_dynamic_column_pivot_count')},"
            f"{int_value(cand_row, 'kls_first_last_separator_dynamic_column_pivots')},"
            f"{int_value(cand_row, 'kls_first_separator_dynamic_column_pivot_count')},"
            f"{int_value(cand_row, 'kls_first_last_separator_extent_dynamic_column_pivots')},"
            f"{int_value(cand_row, 'kls_first_separator_extent_dynamic_column_pivot_count')},"
            f"{int_value(cand_row, 'kls_first_last_separator_dynamic_column_fallbacks')},"
            f"{int_value(cand_row, 'kls_first_separator_dynamic_column_fallback_count')},"
            f"{int_value(cand_row, 'kls_first_last_separator_dynamic_column_rejects')},"
            f"{int_value(cand_row, 'kls_first_separator_dynamic_column_reject_count')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue')},"
            f"{int_value(cand_row, 'kls_first_separator_queue_run_count')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_private_components')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_pipeline_components')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_private_rows')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_pipeline_rows')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_nonempty_threads')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_max_thread_rows')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_partitioned')},"
            f"{int_value(cand_row, 'kls_first_separator_queue_partitioned_count')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_split_components')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_executed')},"
            f"{int_value(cand_row, 'kls_first_separator_queue_executed_run_count')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_executed_private_rows')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_executed_pipeline_rows')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_parallel_private')},"
            f"{int_value(cand_row, 'kls_first_separator_queue_parallel_private_run_count')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_parallel_private_rows')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_parallel_private_threads')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_parallel_pipeline')},"
            f"{int_value(cand_row, 'kls_first_separator_queue_parallel_pipeline_run_count')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_parallel_pipeline_rows')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_parallel_pipeline_threads')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_pipeline_partial')},"
            f"{int_value(cand_row, 'kls_first_separator_queue_pipeline_partial_run_count')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_pipeline_partial_rows')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_pipeline_partial_threads')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_pipeline_wait_partial')},"
            f"{int_value(cand_row, 'kls_first_separator_queue_pipeline_wait_partial_run_count')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_pipeline_wait_partial_rows')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_pipeline_wait_partial_deps')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_pipeline_supernode_update')},"
            f"{int_value(cand_row, 'kls_first_separator_queue_pipeline_supernode_update_run_count')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_pipeline_supernode_update_groups')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_pipeline_supernode_update_rows')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_pipeline_supernode_panel_update')},"
            f"{int_value(cand_row, 'kls_first_separator_queue_pipeline_supernode_panel_update_run_count')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_pipeline_supernode_panel_update_groups')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_pipeline_supernode_panel_update_rows')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_pipeline_pivot_tail')},"
            f"{int_value(cand_row, 'kls_first_separator_queue_pipeline_pivot_tail_run_count')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_pipeline_pivot_tail_rows')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_pipeline_prefix_panel_rebuild')},"
            f"{int_value(cand_row, 'kls_first_separator_queue_pipeline_prefix_panel_rebuild_count')},"
            f"{int_value(cand_row, 'kls_first_last_separator_queue_pipeline_prefix_panel_rebuild_rows')},"
            f"{int_value(cand_row, 'kls_first_auto_skipped_scaled_single_block')},"
            f"{int_value(cand_row, 'kls_first_auto_skipped_scaled_single_block_count')},"
            f"{paper_gap_signal(cand, cand_row)},"
            f"{fmt_share(share(float_value(cand_row, 'refactor_dependency_pipeline_work'), egraph_work))},"
            f"{float_value(cand_row, 'parallel_model_r1'):.6g},"
            f"{float_value(cand_row, 'parallel_model_r2'):.6g},"
            f"{int_value(cand_row, 'parallel_model_recommends_parallel')},"
            f"{int_value(cand_row, 'parallel_task_flow_threads')},"
            f"{int_value(cand_row, 'parallel_task_flow_dependencies')},"
            f"{float_value(cand_row, 'parallel_task_flow_work'):.6g},"
            f"{float_value(cand_row, 'parallel_task_flow_finish_time'):.6g},"
            f"{float_value(cand_row, 'parallel_task_flow_speedup'):.6g},"
            f"{int_value(cand_row, 'parallel_task_flow_recommends_parallel')},"
            f"{fmt_share(share(row_work, egraph_work))},"
            f"{fmt_share(share(compact_work, egraph_work))},"
            f"{fmt_share(share(panel_entries, row_work))},"
            f"{fmt_share(share(private_groups, row_groups))},"
            f"{fmt_ratio(ratio(cand['refactor_99'], ref['refactor_99']))},"
            f"{fmt_ratio(ratio(cand['solve_100'], ref['solve_100']))},"
            f"{fmt_ratio(ratio(cand['initial_factor'], ref['initial_factor']))},"
            f"{fmt_ratio(ratio(cand['analysis'], ref['analysis']))},"
            f"{int_value(cand_row, 'n')},"
            f"{int_value(cand_row, 'nblocks')},"
            f"{int_value(cand_row, 'max_block')},"
            f"{int_value(cand_row, 'scale')},"
            f"{int_value(cand_row, 'offdiag_pivots')},"
            f"{int_value(cand_row, 'fast_repaired_tail_restart_overcompute_columns')},"
            f"{float_value(cand_row, 'fast_repaired_tail_restart_overcompute_work'):.6g},"
            f"{int_value(cand_row, 'fast_repaired_tail_restart_exact_mask')},"
            f"{int_value(cand_row, 'fast_repaired_tail_restart_etree_mask')},"
            f"{int_value(cand_row, 'fast_kls_rebuild_restarts')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_row_pipeline_count')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline_rows')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline_threads')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline_prefix_rows')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline_suffix_rows')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline_gap_rows')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline_etree_tail')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_row_pipeline_etree_tail_count')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline_etree_tail_rows')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline_etree_tail_gap_rows')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline_etree_tail_exact_mask')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline_separator_tail_scope')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_row_pipeline_separator_tail_scope_count')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline_separator_tail_scope_rows')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline_separator_queue')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_row_pipeline_separator_queue_count')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline_separator_private_rows')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline_separator_pipeline_rows')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline_separator_private_threads')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline_separator_partitioned')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline_separator_split_components')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline_pivot_tail_rows')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline_pivot_restarts')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline_supernode_update_groups')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline_supernode_update_rows')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline_supernode_panel_update_groups')},"
            f"{int_value(cand_row, 'fast_kls_block_restart_last_row_pipeline_supernode_panel_update_rows')},"
            f"{int_value(cand_row, 'fast_rejected_pivoting_tail_contiguous')},"
            f"{int_value(cand_row, 'fast_rejected_pivoting_tail_suffix_exact')},"
            f"{int_value(cand_row, 'fast_rejected_pivoting_tail_gap_columns')},"
            f"{int_value(cand_row, 'fast_rejected_pivoting_tail_block_seed_columns')},"
            f"{int_value(cand_row, 'fast_rejected_pivoting_tail_suffix_overcompute_columns')},"
            f"{float_value(cand_row, 'fast_rejected_pivoting_tail_suffix_overcompute_work'):.6g},"
            f"{int_value(cand_row, 'fast_rejected_pivoting_tail_etree_edges')},"
            f"{int_value(cand_row, 'fast_rejected_pivoting_tail_etree_roots')},"
            f"{int_value(cand_row, 'fast_rejected_pivoting_tail_etree_leaves')},"
            f"{int_value(cand_row, 'fast_rejected_pivoting_tail_etree_max_fanout')},"
            f"{int_value(cand_row, 'fast_rejected_pivoting_tail_etree_levels')},"
            f"{int_value(cand_row, 'fast_rejected_pivoting_tail_etree_max_width')},"
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
            f"{int_value(cand_row, 'refactor_supernode_panel_count')},"
            f"{int_value(cand_row, 'refactor_supernode_panel_used_count')},"
            f"{int_value(cand_row, 'refactor_last_supernode_pipeline_tasks')},"
            f"{int_value(cand_row, 'refactor_last_supernode_pipeline_columns')},"
            f"{int_value(cand_row, 'refactor_supernode_pipeline_task_count')},"
            f"{int_value(cand_row, 'refactor_supernode_pipeline_column_count')},"
            f"{int_value(cand_row, 'refactor_last_supernode_update_runs')},"
            f"{int_value(cand_row, 'refactor_last_supernode_update_rows')},"
            f"{int_value(cand_row, 'refactor_last_supernode_update_entries')},"
            f"{int_value(cand_row, 'refactor_supernode_update_run_count')},"
            f"{int_value(cand_row, 'refactor_supernode_update_rows')},"
            f"{int_value(cand_row, 'refactor_supernode_update_entries')},"
            f"{int_value(cand_row, 'refactor_last_supernode_cblas_update_runs')},"
            f"{int_value(cand_row, 'refactor_last_supernode_cblas_update_rows')},"
            f"{int_value(cand_row, 'refactor_last_supernode_cblas_update_entries')},"
            f"{int_value(cand_row, 'refactor_supernode_cblas_update_run_count')},"
            f"{int_value(cand_row, 'refactor_supernode_cblas_update_rows')},"
            f"{int_value(cand_row, 'refactor_supernode_cblas_update_entries')},"
            f"{int_value(cand_row, 'refactor_last_supernode_blocked_update_runs')},"
            f"{int_value(cand_row, 'refactor_last_supernode_blocked_update_rows')},"
            f"{int_value(cand_row, 'refactor_last_supernode_blocked_update_entries')},"
            f"{int_value(cand_row, 'refactor_supernode_blocked_update_run_count')},"
            f"{int_value(cand_row, 'refactor_supernode_blocked_update_rows')},"
            f"{int_value(cand_row, 'refactor_supernode_blocked_update_entries')},"
            f"{int_value(cand_row, 'refactor_supernode_cached_probe_attempt_count')},"
            f"{int_value(cand_row, 'refactor_supernode_cached_probe_panel_hits')},"
            f"{int_value(cand_row, 'refactor_supernode_cached_probe_contiguous')},"
            f"{int_value(cand_row, 'refactor_supernode_cached_probe_allowed')},"
            f"{int_value(cand_row, 'refactor_supernode_cached_probe_allowed_rows')},"
            f"{int_value(cand_row, 'refactor_supernode_cached_probe_applied')},"
            f"{int_value(cand_row, 'refactor_supernode_cached_probe_applied_rows')},"
            f"{int_value(cand_row, 'refactor_supernode_cached_probe_shape_rejects')},"
            f"{int_value(cand_row, 'refactor_supernode_cached_probe_shape_reject_rows')},"
            f"{int_value(cand_row, 'refactor_supernode_cached_probe_stream_rejects')},"
            f"{int_value(cand_row, 'refactor_supernode_cached_probe_stream_reject_rows')},"
            f"{int_value(cand_row, 'refactor_supernode_cached_probe_work_rejects')},"
            f"{int_value(cand_row, 'refactor_supernode_cached_probe_work_reject_rows')},"
            f"{int_value(cand_row, 'refactor_supernode_cached_probe_workspace_rejects')},"
            f"{int_value(cand_row, 'refactor_supernode_cached_probe_workspace_reject_rows')},"
            f"{int_value(cand_row, 'refactor_last_btf_scalar_run_candidates')},"
            f"{int_value(cand_row, 'refactor_last_btf_scalar_run_rows')},"
            f"{int_value(cand_row, 'refactor_last_btf_scalar_run_entries')},"
            f"{int_value(cand_row, 'refactor_last_btf_scalar_run_max_rows')},"
            f"{int_value(cand_row, 'refactor_btf_scalar_run_candidate_count')},"
            f"{int_value(cand_row, 'refactor_btf_scalar_run_rows')},"
            f"{int_value(cand_row, 'refactor_btf_scalar_run_entries')},"
            f"{int_value(cand_row, 'refactor_btf_scalar_run_max_rows')},"
            f"{int_value(cand_row, 'refactor_last_btf_scalar_run_exec_runs')},"
            f"{int_value(cand_row, 'refactor_last_btf_scalar_run_exec_rows')},"
            f"{int_value(cand_row, 'refactor_last_btf_scalar_run_exec_entries')},"
            f"{int_value(cand_row, 'refactor_last_btf_scalar_run_exec_max_rows')},"
            f"{int_value(cand_row, 'refactor_btf_scalar_run_exec_count')},"
            f"{int_value(cand_row, 'refactor_btf_scalar_run_exec_rows')},"
            f"{int_value(cand_row, 'refactor_btf_scalar_run_exec_entries')},"
            f"{int_value(cand_row, 'refactor_btf_scalar_run_exec_max_rows')},"
            f"{int_value(cand_row, 'refactor_btf_scalar_run_group_built')},"
            f"{int_value(cand_row, 'refactor_btf_scalar_run_group_count')},"
            f"{int_value(cand_row, 'refactor_btf_scalar_run_group_current_total')},"
            f"{int_value(cand_row, 'refactor_btf_scalar_run_group_multi_count')},"
            f"{int_value(cand_row, 'refactor_btf_scalar_run_group_multi_current_total')},"
            f"{int_value(cand_row, 'refactor_btf_scalar_run_group_rows')},"
            f"{int_value(cand_row, 'refactor_btf_scalar_run_group_reused_rows')},"
            f"{int_value(cand_row, 'refactor_btf_scalar_run_group_entries')},"
            f"{int_value(cand_row, 'refactor_btf_scalar_run_group_reused_entries')},"
            f"{int_value(cand_row, 'refactor_btf_scalar_run_group_max_currents')},"
            f"{int_value(cand_row, 'refactor_btf_scalar_run_group_max_rows')},"
            f"{int_value(cand_row, 'refactor_last_btf_scalar_run_group_waits')},"
            f"{int_value(cand_row, 'refactor_last_btf_scalar_run_group_overlap_waits')},"
            f"{int_value(cand_row, 'refactor_last_btf_scalar_run_group_max_live')},"
            f"{int_value(cand_row, 'refactor_last_ready_queue_columns')},"
            f"{int_value(cand_row, 'refactor_ready_queue_run_count')},"
            f"{int_value(cand_row, 'row_refactor_group_count')},"
            f"{int_value(cand_row, 'row_refactor_group_single_count')},"
            f"{int_value(cand_row, 'row_refactor_group_batch_count')},"
            f"{int_value(cand_row, 'row_refactor_group_batch_rows')},"
            f"{int_value(cand_row, 'row_refactor_group_scalar_candidate_count')},"
            f"{int_value(cand_row, 'row_refactor_group_scalar_candidate_rows')},"
            f"{int_value(cand_row, 'row_refactor_group_scalar_short_count')},"
            f"{int_value(cand_row, 'row_refactor_group_scalar_short_rows')},"
            f"{int_value(cand_row, 'row_refactor_group_scalar_stop_level_mismatch_count')},"
            f"{int_value(cand_row, 'row_refactor_group_scalar_stop_internal_dep_count')},"
            f"{int_value(cand_row, 'row_refactor_group_scalar_stop_next_segment_count')},"
            f"{int_value(cand_row, 'row_refactor_group_scalar_stop_max_width_count')},"
            f"{int_value(cand_row, 'row_refactor_group_scalar_stop_matrix_end_count')},"
            f"{int_value(cand_row, 'row_refactor_group_generic_count')},"
            f"{int_value(cand_row, 'row_refactor_group_dense_count')},"
            f"{int_value(cand_row, 'row_refactor_group_level_count')},"
            f"{int_value(cand_row, 'row_refactor_group_level_max_width')},"
            f"{int_value(cand_row, 'row_refactor_group_cluster_levels')},"
            f"{int_value(cand_row, 'row_refactor_group_pipeline_groups')},"
            f"{int_value(cand_row, 'row_refactor_group_pipeline_rows')},"
            f"{float_value(cand_row, 'row_refactor_group_pipeline_work'):.6g},"
            f"{float_value(cand_row, 'row_refactor_total_group_work'):.6g},"
            f"{int_value(cand_row, 'row_refactor_auto_enabled')},"
            f"{int_value(cand_row, 'row_refactor_auto_values_ready')},"
            f"{int_value(cand_row, 'row_refactor_auto_work_allowed')},"
            f"{int_value(cand_row, 'row_refactor_auto_should_run')},"
            f"{float_value(cand_row, 'row_refactor_auto_lower_bound_work'):.6g},"
            f"{int_value(cand_row, 'row_refactor_auto_lower_bound_rejected')},"
            f"{int_value(cand_row, 'row_refactor_auto_pattern_build_failed')},"
            f"{int_value(cand_row, 'row_refactor_auto_value_copy_failed')},"
            f"{int_value(cand_row, 'row_refactor_auto_model_recommended')},"
            f"{int_value(cand_row, 'row_refactor_auto_model_attempted')},"
            f"{int_value(cand_row, 'row_refactor_auto_model_accepted')},"
            f"{int_value(cand_row, 'row_solve_parallel_run_count')},"
            f"{int_value(cand_row, 'row_solve_parallel_l_slice_runs')},"
            f"{int_value(cand_row, 'row_solve_parallel_u_slice_runs')},"
            f"{int_value(cand_row, 'row_solve_parallel_l_sparse_level_runs')},"
            f"{int_value(cand_row, 'row_solve_parallel_u_sparse_level_runs')},"
            f"{int_value(cand_row, 'row_solve_thread_count')},"
            f"{int_value(cand_row, 'row_solve_l_thread_max_rect_entries')},"
            f"{int_value(cand_row, 'row_solve_u_thread_max_rect_entries')},"
            f"{int_value(cand_row, 'row_solve_partition_ready')},"
            f"{int_value(cand_row, 'row_solve_partition_slices')},"
            f"{int_value(cand_row, 'row_solve_l_sparse_level_count')},"
            f"{int_value(cand_row, 'row_solve_l_sparse_cluster_levels')},"
            f"{int_value(cand_row, 'row_solve_l_sparse_level_max_width')},"
            f"{int_value(cand_row, 'row_solve_l_dense_tail_start')},"
            f"{int_value(cand_row, 'row_solve_l_dense_tail_rows')},"
            f"{int_value(cand_row, 'row_solve_l_dense_tail_entries')},"
            f"{int_value(cand_row, 'row_solve_l_slice_max_entries')},"
            f"{int_value(cand_row, 'row_solve_l_segmented_rows')},"
            f"{int_value(cand_row, 'row_solve_l_rect_entries')},"
            f"{int_value(cand_row, 'row_solve_l_tri_entries')},"
            f"{int_value(cand_row, 'row_solve_u_sparse_level_count')},"
            f"{int_value(cand_row, 'row_solve_u_sparse_cluster_levels')},"
            f"{int_value(cand_row, 'row_solve_u_sparse_level_max_width')},"
            f"{int_value(cand_row, 'row_solve_u_dense_tail_start')},"
            f"{int_value(cand_row, 'row_solve_u_dense_tail_rows')},"
            f"{int_value(cand_row, 'row_solve_u_dense_tail_entries')},"
            f"{int_value(cand_row, 'row_solve_u_slice_max_entries')},"
            f"{int_value(cand_row, 'row_solve_u_segmented_rows')},"
            f"{int_value(cand_row, 'row_solve_u_rect_entries')},"
            f"{int_value(cand_row, 'row_solve_u_tri_entries')},"
            f"{int_value(cand_row, 'row_refactor_last_run')},"
            f"{int_value(cand_row, 'row_refactor_last_checked')},"
            f"{int_value(cand_row, 'row_refactor_last_parallel')},"
            f"{int_value(cand_row, 'row_refactor_run_count')},"
            f"{int_value(cand_row, 'row_refactor_checked_run_count')},"
            f"{int_value(cand_row, 'row_refactor_parallel_run_count')},"
            f"{int_value(cand_row, 'row_refactor_segment_count')},"
            f"{int_value(cand_row, 'row_refactor_segment_rows')},"
            f"{int_value(cand_row, 'row_refactor_segment_max_width')},"
            f"{float_value(cand_row, 'row_refactor_segment_dense_entries'):.6g},"
            f"{float_value(cand_row, 'row_refactor_segment_trailing_entries'):.6g},"
            f"{int_value(cand_row, 'row_refactor_dense_segment_count')},"
            f"{int_value(cand_row, 'row_refactor_dense_segment_rows')},"
            f"{int_value(cand_row, 'row_refactor_dense_segment_max_width')},"
            f"{float_value(cand_row, 'row_refactor_dense_segment_dense_entries'):.6g},"
            f"{float_value(cand_row, 'row_refactor_dense_segment_trailing_entries'):.6g},"
            f"{int_value(cand_row, 'row_refactor_compact_dense_panel_eligible_count')},"
            f"{int_value(cand_row, 'row_refactor_compact_dense_panel_eligible_rows')},"
            f"{float_value(cand_row, 'row_refactor_compact_dense_panel_update_work'):.6g},"
            f"{float_value(cand_row, 'row_refactor_compact_dense_panel_entries'):.6g},"
            f"{int_value(cand_row, 'row_refactor_compact_dense_panel_persistent_groups')},"
            f"{int_value(cand_row, 'row_refactor_compact_dense_panel_persistent_entries')},"
            f"{int_value(cand_row, 'row_refactor_last_compact_dense_panel_persistent')},"
            f"{int_value(cand_row, 'row_refactor_compact_dense_panel_persistent_run_count')},"
            f"{int_value(cand_row, 'row_refactor_last_compact_dense_panel')},"
            f"{int_value(cand_row, 'row_refactor_compact_dense_panel_count')},"
            f"{int_value(cand_row, 'row_refactor_last_compact_dense_panel_blocked')},"
            f"{int_value(cand_row, 'row_refactor_compact_dense_panel_blocked_run_count')},"
            f"{int_value(cand_row, 'row_refactor_compact_dense_panel_blocked_rows')},"
            f"{int_value(cand_row, 'row_refactor_compact_dense_panel_blocked_entries')},"
            f"{int_value(cand_row, 'row_refactor_native_row_panel_enabled')},"
            f"{int_value(cand_row, 'row_refactor_last_native_row_panel')},"
            f"{int_value(cand_row, 'row_refactor_native_row_panel_count')},"
            f"{int_value(cand_row, 'row_refactor_native_row_panel_rows')},"
            f"{int_value(cand_row, 'row_refactor_native_row_panel_entries')},"
            f"{int_value(cand_row, 'row_refactor_native_row_panel_blocked_count')},"
            f"{int_value(cand_row, 'row_refactor_native_row_panel_blocked_rows')},"
            f"{int_value(cand_row, 'row_refactor_native_row_panel_blocked_entries')},"
            f"{int_value(cand_row, 'row_refactor_native_row_panel_fallback_count')},"
            f"{int_value(cand_row, 'row_refactor_native_row_panel_checked_reject_count')},"
            f"{int_value(cand_row, 'row_refactor_last_compact_supernode_partial_update')},"
            f"{int_value(cand_row, 'row_refactor_compact_supernode_partial_update_count')},"
            f"{int_value(cand_row, 'row_refactor_compact_supernode_partial_update_rows')},"
            f"{int_value(cand_row, 'row_refactor_compact_supernode_partial_update_entries')},"
            f"{int_value(cand_row, 'row_refactor_last_local_ready_groups')},"
            f"{int_value(cand_row, 'row_refactor_local_ready_group_count')},"
            f"{int_value(cand_row, 'row_refactor_last_private_ready_groups')},"
            f"{int_value(cand_row, 'row_refactor_private_ready_group_count')},"
            f"{int_value(cand_row, 'row_refactor_last_separator_flop_queue')},"
            f"{int_value(cand_row, 'row_refactor_separator_flop_queue_run_count')},"
            f"{int_value(cand_row, 'row_refactor_last_separator_flop_ordered_private')},"
            f"{int_value(cand_row, 'row_refactor_separator_flop_ordered_private_run_count')},"
            f"{int_value(cand_row, 'row_refactor_last_separator_flop_components')},"
            f"{int_value(cand_row, 'row_refactor_separator_flop_component_count')},"
            f"{int_value(cand_row, 'row_refactor_last_separator_flop_private_groups')},"
            f"{int_value(cand_row, 'row_refactor_last_separator_flop_pipeline_groups')},"
            f"{int_value(cand_row, 'row_refactor_last_separator_flop_closure_groups')},"
            f"{int_value(cand_row, 'row_refactor_last_separator_flop_private_threads')},"
            f"{int_value(cand_row, 'row_refactor_last_separator_flop_private_min_groups')},"
            f"{int_value(cand_row, 'row_refactor_last_separator_flop_private_max_groups')},"
            f"{float_value(cand_row, 'row_refactor_last_separator_flop_private_min_work'):.6g},"
            f"{float_value(cand_row, 'row_refactor_last_separator_flop_private_max_work'):.6g},"
            f"{int_value(cand_row, 'row_refactor_separator_flop_private_group_count')},"
            f"{int_value(cand_row, 'row_refactor_separator_flop_pipeline_group_count')},"
            f"{int_value(cand_row, 'row_refactor_separator_flop_closure_group_count')},"
            f"{int_value(cand_row, 'refactor_dependency_cluster_levels')},"
            f"{int_value(cand_row, 'refactor_dependency_pipeline_columns')},"
            f"{float_value(cand_row, 'refactor_dependency_work'):.6g},"
            f"{float_value(cand_row, 'refactor_dependency_pipeline_work'):.6g},"
            f"{int_value(cand_row, 'refactor_last_btf_scalar_run_candidates')},"
            f"{int_value(cand_row, 'refactor_last_btf_scalar_run_rows')},"
            f"{int_value(cand_row, 'refactor_last_btf_scalar_run_entries')},"
            f"{int_value(cand_row, 'refactor_last_btf_scalar_run_max_rows')},"
            f"{int_value(cand_row, 'refactor_last_btf_scalar_run_exec_runs')},"
            f"{int_value(cand_row, 'refactor_last_btf_scalar_run_exec_rows')},"
            f"{int_value(cand_row, 'refactor_last_btf_scalar_run_exec_entries')},"
            f"{int_value(cand_row, 'refactor_last_btf_scalar_run_exec_max_rows')},"
            f"{int_value(cand_row, 'refactor_btf_scalar_run_group_multi_count')},"
            f"{int_value(cand_row, 'refactor_btf_scalar_run_group_multi_current_total')},"
            f"{int_value(cand_row, 'refactor_btf_scalar_run_group_reused_entries')},"
            f"{int_value(cand_row, 'refactor_btf_scalar_run_group_max_currents')},"
            f"{int_value(cand_row, 'refactor_last_btf_scalar_run_group_waits')},"
            f"{int_value(cand_row, 'refactor_last_btf_scalar_run_group_overlap_waits')},"
            f"{int_value(cand_row, 'refactor_last_btf_scalar_run_group_max_live')}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
