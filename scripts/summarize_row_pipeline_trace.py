#!/usr/bin/env python3
"""Summarize KLS row-pipeline trace stderr files.

The row-pipeline trace is intentionally line-oriented so long timeout runs can
be inspected after a process cap. This helper extracts the counters most useful
for comparing first-factor tail experiments.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
from typing import Iterable


KEY_VALUE_RE = re.compile(r"([A-Za-z0-9_]+)=([^ \n]+)")

PRODUCER_COUNTERS = [
    "producer_batches",
    "producer_targets",
    "producer_candidate_targets",
    "producer_candidate_target_u_entries",
    "producer_underfilled_target_u_entries",
    "producer_low_saved_stream_target_u_entries",
    "producer_active_catchup_attempts",
    "producer_active_catchup_deps",
    "producer_active_catchup_targets",
    "compact_window_fills",
    "compact_window_evictions",
    "compact_window_overflows",
    "compact_window_claim_attempts",
    "compact_window_claims",
    "compact_window_claim_misses",
    "compact_window_claim_stale_clears",
    "compact_window_claim_group_scatters",
    "compact_window_claim_group_detaches",
    "compact_window_claim_group_run_probes",
    "compact_window_claim_group_run_states",
    "compact_window_claim_group_run_consecutive",
    "compact_window_claim_group_run_near64",
    "compact_window_claim_group_run_near512",
    "compact_window_claim_group_run_max_consecutive",
    "compact_window_claim_group_run_stealable",
    "compact_window_claim_group_run_max_stealable",
    "compact_window_claim_group_run_unreserved_near64",
    "compact_window_claim_run_reservations",
    "compact_window_claim_run_rows",
    "compact_window_claim_run_recomputes",
    "compact_window_claim_run_updates",
    "compact_window_claim_run_update_targets",
    "compact_window_span_owner_reservations",
    "compact_window_span_owner_rows",
    "compact_window_span_owner_owned_rows",
    "compact_window_span_owner_hole_rows",
    "compact_window_span_owner_deps",
    "compact_window_span_owner_unique_deps",
    "compact_window_span_owner_scan_entries",
    "compact_window_span_owner_cols",
    "compact_window_span_owner_slots",
    "compact_window_span_owner_entries",
    "compact_window_span_owner_hole_skips",
    "compact_window_span_owner_hole_skip_rows",
    "compact_window_span_owner_prefix_shrinks",
    "compact_window_span_owner_prefix_shrink_rows",
    "compact_window_span_owner_oversize_skips",
    "compact_window_span_owner_oversize_slots",
    "compact_window_span_owner_link_skips",
    "compact_window_span_owner_link_skip_entries",
    "compact_window_span_owner_scan_skips",
    "compact_window_span_owner_scan_skip_entries",
    "compact_window_span_owner_payoff_skips",
    "compact_window_span_owner_payoff_skip_rows",
    "compact_window_span_owner_payoff_skip_entries",
    "compact_window_span_owner_payoff_skip_scan_entries",
    "compact_window_span_shared_owner_reservations",
    "compact_window_span_shared_owner_rows",
    "compact_window_span_shared_owner_applies",
    "compact_window_span_shared_owner_entries",
    "compact_window_claim_run_output_surfaces",
    "compact_window_claim_run_output_states",
    "compact_window_claim_run_output_deps",
    "compact_window_claim_run_output_unique_deps",
    "compact_window_claim_run_output_duplicate_deps",
    "compact_window_claim_run_output_scan_entries",
    "compact_window_claim_run_output_group_scan_entries",
    "compact_window_claim_span_output_surfaces",
    "compact_window_claim_span_output_states",
    "compact_window_claim_span_output_deps",
    "compact_window_claim_span_output_unique_deps",
    "compact_window_claim_span_output_duplicate_deps",
    "compact_window_claim_span_output_scan_entries",
    "compact_window_claim_span_output_group_scan_entries",
    "compact_window_claim_span_output_value_entries",
    "compact_window_claim_span_output_state_col_slots",
    "compact_window_claim_span_output_unique_cols",
    "compact_window_claim_span_output_panel_slots",
    "compact_window_claim_span_output_dense_slots",
    "compact_window_claim_span_output_max_states",
    "compact_window_claim_span_output_max_unique_cols",
    "compact_window_delayed_output_skips",
    "compact_window_delayed_output_replays",
    "compact_window_delayed_output_deps",
    "compact_window_delayed_output_entries",
    "compact_window_delayed_output_scan_entries",
    "compact_window_delayed_output_seek_skips",
    "compact_window_deferred_output_cache_stores",
    "compact_window_deferred_output_cache_applies",
    "compact_window_deferred_output_cache_l_entries",
    "compact_window_deferred_output_cache_entries",
    "compact_window_deferred_output_cache_disables",
    "compact_window_deferred_output_cache_overflows",
    "compact_window_delayed_group_replay_surfaces",
    "compact_window_delayed_group_replay_states",
    "compact_window_delayed_group_replay_deps",
    "compact_window_delayed_group_replay_unique_deps",
    "compact_window_delayed_group_replay_duplicate_deps",
    "compact_window_delayed_group_replay_scan_entries",
    "compact_window_delayed_group_replay_group_scan_entries",
    "compact_window_claim_group_output_surfaces",
    "compact_window_claim_group_output_states",
    "compact_window_claim_group_output_deps",
    "compact_window_claim_group_output_unique_deps",
    "compact_window_claim_group_output_duplicate_deps",
    "compact_window_claim_group_output_scan_entries",
    "compact_window_claim_group_output_group_scan_entries",
    "compact_window_probes",
    "compact_window_batches",
    "compact_window_stream_u_entries",
    "compact_window_targets",
    "compact_window_target_u_entries",
    "compact_window_state_rows",
    "compact_window_group_steps",
    "compact_window_group_targets",
    "compact_window_group_cols",
    "compact_window_group_values",
    "compact_window_group_merges",
    "compact_window_group_merge_targets",
    "compact_window_group_merge_cols",
    "compact_window_group_merge_values",
    "owner_surface_probes",
    "owner_surface_probe_u_entries",
    "owner_surface_scanned_rows",
    "owner_surface_targets",
    "owner_surface_target_u_entries",
    "owner_surface_scanned_input_entries",
    "owner_surface_target_input_entries",
    "owner_surface_max_targets",
    "owner_surface_max_target_input_entries",
    "producer_probe_workers",
    "producer_probe_lookahead",
    "producer_ready_roots",
    "producer_underfilled",
    "producer_underfilled_targets",
    "producer_low_saved_stream",
    "producer_reject_bad_state",
    "producer_reject_epoch",
    "producer_reject_dep_absent",
    "producer_reject_not_root",
    "producer_reject_not_ready",
    "producer_reject_supernode",
    "producer_reject_cached_panel",
]


def parse_int(value: str | None, default: int = 0) -> int:
    if value is None:
        return default
    try:
        return int(value)
    except ValueError:
        return default


def split_completed(value: str | None) -> tuple[int, int]:
    if value is None or "/" not in value:
        return 0, 0
    done, total = value.split("/", 1)
    return parse_int(done), parse_int(total)


def ratio(numerator: int, denominator: int) -> float | None:
    if denominator <= 0:
        return None
    return float(numerator) / float(denominator)


def fields(line: str) -> dict[str, str]:
    return {match.group(1): match.group(2) for match in KEY_VALUE_RE.finditer(line)}


def summarize(path: pathlib.Path) -> dict[str, object]:
    last_trace: dict[str, str] = {}
    trace_event_count = 0
    long_rows = 0
    long_max_row = 0
    long_pivot_rows = 0
    long_scalar = 0
    long_internal = 0
    long_output = 0
    long_producer_targets = 0
    long_producer_target_u = 0
    long_producer_state_rows = 0
    long_producer_unique_state_rows = 0
    long_candidate_targets = 0
    long_candidate_target_u = 0
    long_underfilled_target_u = 0
    long_low_saved_stream_target_u = 0
    long_active_catchup_attempts = 0
    long_active_catchup_deps = 0
    long_active_catchup_targets = 0
    long_compact_window_claim_attempts = 0
    long_compact_window_claims = 0
    long_compact_window_claim_misses = 0
    long_compact_window_claim_stale_clears = 0
    long_compact_window_claim_group_scatters = 0
    long_compact_window_claim_group_detaches = 0
    long_compact_window_claim_group_run_probes = 0
    long_compact_window_claim_group_run_states = 0
    long_compact_window_claim_group_run_consecutive = 0
    long_compact_window_claim_group_run_near64 = 0
    long_compact_window_claim_group_run_near512 = 0
    long_compact_window_claim_group_run_max_consecutive = 0
    long_compact_window_claim_group_run_stealable = 0
    long_compact_window_claim_group_run_max_stealable = 0
    long_compact_window_claim_group_run_unreserved_near64 = 0
    long_compact_window_claim_run_reservations = 0
    long_compact_window_claim_run_rows = 0
    long_compact_window_claim_run_recomputes = 0
    long_compact_window_claim_run_updates = 0
    long_compact_window_claim_run_update_targets = 0
    long_compact_window_claim_run_output_surfaces = 0
    long_compact_window_claim_run_output_states = 0
    long_compact_window_claim_run_output_deps = 0
    long_compact_window_claim_run_output_unique_deps = 0
    long_compact_window_claim_run_output_duplicate_deps = 0
    long_compact_window_claim_run_output_scan_entries = 0
    long_compact_window_claim_run_output_group_scan_entries = 0
    long_compact_window_claim_span_output_surfaces = 0
    long_compact_window_claim_span_output_states = 0
    long_compact_window_claim_span_output_deps = 0
    long_compact_window_claim_span_output_unique_deps = 0
    long_compact_window_claim_span_output_duplicate_deps = 0
    long_compact_window_claim_span_output_scan_entries = 0
    long_compact_window_claim_span_output_group_scan_entries = 0
    long_compact_window_claim_span_output_value_entries = 0
    long_compact_window_claim_span_output_state_col_slots = 0
    long_compact_window_claim_span_output_unique_cols = 0
    long_compact_window_claim_span_output_panel_slots = 0
    long_compact_window_claim_span_output_dense_slots = 0
    long_compact_window_claim_span_output_max_states = 0
    long_compact_window_claim_span_output_max_unique_cols = 0
    long_compact_window_delayed_output_skips = 0
    long_compact_window_delayed_output_replays = 0
    long_compact_window_delayed_output_deps = 0
    long_compact_window_delayed_output_entries = 0
    long_compact_window_delayed_output_scan_entries = 0
    long_compact_window_delayed_output_seek_skips = 0
    long_compact_window_delayed_group_replay_surfaces = 0
    long_compact_window_delayed_group_replay_states = 0
    long_compact_window_delayed_group_replay_deps = 0
    long_compact_window_delayed_group_replay_unique_deps = 0
    long_compact_window_delayed_group_replay_duplicate_deps = 0
    long_compact_window_delayed_group_replay_scan_entries = 0
    long_compact_window_delayed_group_replay_group_scan_entries = 0
    long_compact_window_claim_group_output_surfaces = 0
    long_compact_window_claim_group_output_states = 0
    long_compact_window_claim_group_output_deps = 0
    long_compact_window_claim_group_output_unique_deps = 0
    long_compact_window_claim_group_output_duplicate_deps = 0
    long_compact_window_claim_group_output_scan_entries = 0
    long_compact_window_claim_group_output_group_scan_entries = 0
    long_compact_window_targets = 0
    long_compact_window_target_u = 0
    long_compact_window_batches = 0
    long_compact_window_stream_u = 0
    long_compact_window_state_rows = 0
    long_compact_window_unique_state_rows = 0
    long_compact_union_batches = 0
    long_compact_union_targets = 0
    long_compact_union_cols = 0
    long_compact_union_values = 0
    long_compact_union_skips = 0
    long_compact_union_skip_targets = 0
    long_compact_union_skip_cols = 0
    long_compact_union_skip_values = 0
    long_compact_union_skip_sparse_values = 0
    long_compact_group_steps = 0
    long_compact_group_targets = 0
    long_compact_group_cols = 0
    long_compact_group_values = 0
    long_compact_group_merges = 0
    long_compact_group_merge_targets = 0
    long_compact_group_merge_cols = 0
    long_compact_group_merge_values = 0
    long_owner_surface_probes = 0
    long_owner_surface_probe_u = 0
    long_owner_surface_scanned = 0
    long_owner_surface_targets = 0
    long_owner_surface_target_u = 0
    long_owner_surface_scanned_input = 0
    long_owner_surface_target_input = 0
    long_owner_surface_max_targets = 0
    long_owner_surface_max_target_input = 0
    long_panel_rows = 0
    long_panel_entries = 0
    long_pivot_scalar = 0
    long_pivot_output = 0
    long_pivot_producer_rows = 0
    long_pivot_producer_targets = 0
    long_pivot_producer_target_u = 0
    long_pivot_probe_lookahead = 0
    long_pivot_ready_roots = 0
    long_pivot_underfilled = 0
    long_pivot_reject_bad_state = 0
    long_pivot_reject_epoch = 0

    with path.open("r", encoding="utf-8", errors="replace") as trace:
        for line in trace:
            if "KLS row-pipeline trace:" in line:
                last_trace = fields(line)
                trace_event_count += 1
                continue
            if "KLS row-pipeline long-row:" not in line:
                continue
            row = fields(line)
            long_rows += 1
            long_max_row = max(long_max_row, parse_int(row.get("row")))
            pivoted = parse_int(row.get("pivoted"))
            long_pivot_rows += 1 if pivoted else 0
            long_scalar += parse_int(row.get("scalar_u_entries"))
            long_internal += parse_int(row.get("scalar_u_internal"))
            long_output += parse_int(row.get("scalar_u_output"))
            long_producer_targets += parse_int(row.get("producer_targets"))
            long_producer_target_u += parse_int(row.get("producer_target_u_entries"))
            long_producer_state_rows += parse_int(row.get("producer_state_rows"))
            long_producer_unique_state_rows += parse_int(
                row.get("producer_unique_state_rows")
            )
            long_candidate_targets += parse_int(row.get("producer_candidate_targets"))
            long_candidate_target_u += parse_int(
                row.get("producer_candidate_target_u_entries")
            )
            long_underfilled_target_u += parse_int(
                row.get("producer_underfilled_target_u_entries")
            )
            long_low_saved_stream_target_u += parse_int(
                row.get("producer_low_saved_stream_target_u_entries")
            )
            long_active_catchup_attempts += parse_int(
                row.get("producer_active_catchup_attempts")
            )
            long_active_catchup_deps += parse_int(
                row.get("producer_active_catchup_deps")
            )
            long_active_catchup_targets += parse_int(
                row.get("producer_active_catchup_targets")
            )
            long_compact_window_claim_attempts += parse_int(
                row.get("compact_window_claim_attempts")
            )
            long_compact_window_claims += parse_int(
                row.get("compact_window_claims")
            )
            long_compact_window_claim_misses += parse_int(
                row.get("compact_window_claim_misses")
            )
            long_compact_window_claim_stale_clears += parse_int(
                row.get("compact_window_claim_stale_clears")
            )
            long_compact_window_claim_group_scatters += parse_int(
                row.get("compact_window_claim_group_scatters")
            )
            long_compact_window_claim_group_detaches += parse_int(
                row.get("compact_window_claim_group_detaches")
            )
            long_compact_window_claim_group_run_probes += parse_int(
                row.get("compact_window_claim_group_run_probes")
            )
            long_compact_window_claim_group_run_states += parse_int(
                row.get("compact_window_claim_group_run_states")
            )
            long_compact_window_claim_group_run_consecutive += parse_int(
                row.get("compact_window_claim_group_run_consecutive")
            )
            long_compact_window_claim_group_run_near64 += parse_int(
                row.get("compact_window_claim_group_run_near64")
            )
            long_compact_window_claim_group_run_near512 += parse_int(
                row.get("compact_window_claim_group_run_near512")
            )
            long_compact_window_claim_group_run_max_consecutive = max(
                long_compact_window_claim_group_run_max_consecutive,
                parse_int(row.get("compact_window_claim_group_run_max_consecutive")),
            )
            long_compact_window_claim_group_run_stealable += parse_int(
                row.get("compact_window_claim_group_run_stealable")
            )
            long_compact_window_claim_group_run_max_stealable = max(
                long_compact_window_claim_group_run_max_stealable,
                parse_int(row.get("compact_window_claim_group_run_max_stealable")),
            )
            long_compact_window_claim_group_run_unreserved_near64 += parse_int(
                row.get("compact_window_claim_group_run_unreserved_near64")
            )
            long_compact_window_claim_run_reservations += parse_int(
                row.get("compact_window_claim_run_reservations")
            )
            long_compact_window_claim_run_rows += parse_int(
                row.get("compact_window_claim_run_rows")
            )
            long_compact_window_claim_run_recomputes += parse_int(
                row.get("compact_window_claim_run_recomputes")
            )
            long_compact_window_claim_run_updates += parse_int(
                row.get("compact_window_claim_run_updates")
            )
            long_compact_window_claim_run_update_targets += parse_int(
                row.get("compact_window_claim_run_update_targets")
            )
            long_compact_window_claim_run_output_surfaces += parse_int(
                row.get("compact_window_claim_run_output_surfaces")
            )
            long_compact_window_claim_run_output_states += parse_int(
                row.get("compact_window_claim_run_output_states")
            )
            long_compact_window_claim_run_output_deps += parse_int(
                row.get("compact_window_claim_run_output_deps")
            )
            long_compact_window_claim_run_output_unique_deps += parse_int(
                row.get("compact_window_claim_run_output_unique_deps")
            )
            long_compact_window_claim_run_output_duplicate_deps += parse_int(
                row.get("compact_window_claim_run_output_duplicate_deps")
            )
            long_compact_window_claim_run_output_scan_entries += parse_int(
                row.get("compact_window_claim_run_output_scan_entries")
            )
            long_compact_window_claim_run_output_group_scan_entries += parse_int(
                row.get("compact_window_claim_run_output_group_scan_entries")
            )
            long_compact_window_claim_span_output_surfaces += parse_int(
                row.get("compact_window_claim_span_output_surfaces")
            )
            long_compact_window_claim_span_output_states += parse_int(
                row.get("compact_window_claim_span_output_states")
            )
            long_compact_window_claim_span_output_deps += parse_int(
                row.get("compact_window_claim_span_output_deps")
            )
            long_compact_window_claim_span_output_unique_deps += parse_int(
                row.get("compact_window_claim_span_output_unique_deps")
            )
            long_compact_window_claim_span_output_duplicate_deps += parse_int(
                row.get("compact_window_claim_span_output_duplicate_deps")
            )
            long_compact_window_claim_span_output_scan_entries += parse_int(
                row.get("compact_window_claim_span_output_scan_entries")
            )
            long_compact_window_claim_span_output_group_scan_entries += parse_int(
                row.get("compact_window_claim_span_output_group_scan_entries")
            )
            long_compact_window_claim_span_output_value_entries += parse_int(
                row.get("compact_window_claim_span_output_value_entries")
            )
            long_compact_window_claim_span_output_state_col_slots += parse_int(
                row.get("compact_window_claim_span_output_state_col_slots")
            )
            long_compact_window_claim_span_output_unique_cols += parse_int(
                row.get("compact_window_claim_span_output_unique_cols")
            )
            long_compact_window_claim_span_output_panel_slots += parse_int(
                row.get("compact_window_claim_span_output_panel_slots")
            )
            long_compact_window_claim_span_output_dense_slots += parse_int(
                row.get("compact_window_claim_span_output_dense_slots")
            )
            long_compact_window_claim_span_output_max_states = max(
                long_compact_window_claim_span_output_max_states,
                parse_int(row.get("compact_window_claim_span_output_max_states")),
            )
            long_compact_window_claim_span_output_max_unique_cols = max(
                long_compact_window_claim_span_output_max_unique_cols,
                parse_int(
                    row.get("compact_window_claim_span_output_max_unique_cols")
                ),
            )
            long_compact_window_delayed_output_skips += parse_int(
                row.get("compact_window_delayed_output_skips")
            )
            long_compact_window_delayed_output_replays += parse_int(
                row.get("compact_window_delayed_output_replays")
            )
            long_compact_window_delayed_output_deps += parse_int(
                row.get("compact_window_delayed_output_deps")
            )
            long_compact_window_delayed_output_entries += parse_int(
                row.get("compact_window_delayed_output_entries")
            )
            long_compact_window_delayed_output_scan_entries += parse_int(
                row.get("compact_window_delayed_output_scan_entries")
            )
            long_compact_window_delayed_output_seek_skips += parse_int(
                row.get("compact_window_delayed_output_seek_skips")
            )
            long_compact_window_delayed_group_replay_surfaces += parse_int(
                row.get("compact_window_delayed_group_replay_surfaces")
            )
            long_compact_window_delayed_group_replay_states += parse_int(
                row.get("compact_window_delayed_group_replay_states")
            )
            long_compact_window_delayed_group_replay_deps += parse_int(
                row.get("compact_window_delayed_group_replay_deps")
            )
            long_compact_window_delayed_group_replay_unique_deps += parse_int(
                row.get("compact_window_delayed_group_replay_unique_deps")
            )
            long_compact_window_delayed_group_replay_duplicate_deps += parse_int(
                row.get("compact_window_delayed_group_replay_duplicate_deps")
            )
            long_compact_window_delayed_group_replay_scan_entries += parse_int(
                row.get("compact_window_delayed_group_replay_scan_entries")
            )
            long_compact_window_delayed_group_replay_group_scan_entries += parse_int(
                row.get("compact_window_delayed_group_replay_group_scan_entries")
            )
            long_compact_window_claim_group_output_surfaces += parse_int(
                row.get("compact_window_claim_group_output_surfaces")
            )
            long_compact_window_claim_group_output_states += parse_int(
                row.get("compact_window_claim_group_output_states")
            )
            long_compact_window_claim_group_output_deps += parse_int(
                row.get("compact_window_claim_group_output_deps")
            )
            long_compact_window_claim_group_output_unique_deps += parse_int(
                row.get("compact_window_claim_group_output_unique_deps")
            )
            long_compact_window_claim_group_output_duplicate_deps += parse_int(
                row.get("compact_window_claim_group_output_duplicate_deps")
            )
            long_compact_window_claim_group_output_scan_entries += parse_int(
                row.get("compact_window_claim_group_output_scan_entries")
            )
            long_compact_window_claim_group_output_group_scan_entries += parse_int(
                row.get("compact_window_claim_group_output_group_scan_entries")
            )
            long_compact_window_targets += parse_int(
                row.get("compact_window_targets")
            )
            long_compact_window_target_u += parse_int(
                row.get("compact_window_target_u_entries")
            )
            long_compact_window_batches += parse_int(
                row.get("compact_window_batches")
            )
            long_compact_window_stream_u += parse_int(
                row.get("compact_window_stream_u_entries")
            )
            long_compact_window_state_rows += parse_int(
                row.get("compact_window_state_rows")
            )
            long_compact_window_unique_state_rows += parse_int(
                row.get("compact_window_unique_state_rows")
            )
            long_compact_union_batches += parse_int(
                row.get("compact_window_union_batches")
            )
            long_compact_union_targets += parse_int(
                row.get("compact_window_union_targets")
            )
            long_compact_union_cols += parse_int(
                row.get("compact_window_union_cols")
            )
            long_compact_union_values += parse_int(
                row.get("compact_window_union_values")
            )
            long_compact_union_skips += parse_int(
                row.get("compact_window_union_skips")
            )
            long_compact_union_skip_targets += parse_int(
                row.get("compact_window_union_skip_targets")
            )
            long_compact_union_skip_cols += parse_int(
                row.get("compact_window_union_skip_cols")
            )
            long_compact_union_skip_values += parse_int(
                row.get("compact_window_union_skip_values")
            )
            long_compact_union_skip_sparse_values += parse_int(
                row.get("compact_window_union_skip_sparse_values")
            )
            long_compact_group_steps += parse_int(
                row.get("compact_window_group_steps")
            )
            long_compact_group_targets += parse_int(
                row.get("compact_window_group_targets")
            )
            long_compact_group_cols += parse_int(
                row.get("compact_window_group_cols")
            )
            long_compact_group_values += parse_int(
                row.get("compact_window_group_values")
            )
            long_compact_group_merges += parse_int(
                row.get("compact_window_group_merges")
            )
            long_compact_group_merge_targets += parse_int(
                row.get("compact_window_group_merge_targets")
            )
            long_compact_group_merge_cols += parse_int(
                row.get("compact_window_group_merge_cols")
            )
            long_compact_group_merge_values += parse_int(
                row.get("compact_window_group_merge_values")
            )
            long_owner_surface_probes += parse_int(
                row.get("owner_surface_probes")
            )
            long_owner_surface_probe_u += parse_int(
                row.get("owner_surface_probe_u_entries")
            )
            long_owner_surface_scanned += parse_int(
                row.get("owner_surface_scanned_rows")
            )
            long_owner_surface_targets += parse_int(
                row.get("owner_surface_targets")
            )
            long_owner_surface_target_u += parse_int(
                row.get("owner_surface_target_u_entries")
            )
            long_owner_surface_scanned_input += parse_int(
                row.get("owner_surface_scanned_input_entries")
            )
            long_owner_surface_target_input += parse_int(
                row.get("owner_surface_target_input_entries")
            )
            long_owner_surface_max_targets = max(
                long_owner_surface_max_targets,
                parse_int(row.get("owner_surface_max_targets")),
            )
            long_owner_surface_max_target_input = max(
                long_owner_surface_max_target_input,
                parse_int(row.get("owner_surface_max_target_input_entries")),
            )
            long_panel_rows += parse_int(row.get("panel_update_rows"))
            long_panel_entries += parse_int(row.get("panel_update_entries"))
            if pivoted:
                producer_batches = parse_int(row.get("producer_batches"))
                producer_targets = parse_int(row.get("producer_targets"))
                producer_target_u = parse_int(row.get("producer_target_u_entries"))
                long_pivot_scalar += parse_int(row.get("scalar_u_entries"))
                long_pivot_output += parse_int(row.get("scalar_u_output"))
                long_pivot_producer_rows += 1 if producer_batches > 0 else 0
                long_pivot_producer_targets += producer_targets
                long_pivot_producer_target_u += producer_target_u
                long_pivot_probe_lookahead += parse_int(
                    row.get("producer_probe_lookahead")
                )
                long_pivot_ready_roots += parse_int(row.get("producer_ready_roots"))
                long_pivot_underfilled += parse_int(row.get("producer_underfilled"))
                long_pivot_reject_bad_state += parse_int(
                    row.get("producer_reject_bad_state")
                )
                long_pivot_reject_epoch += parse_int(
                    row.get("producer_reject_epoch")
                )

    completed, total = split_completed(last_trace.get("completed"))
    scalar = parse_int(last_trace.get("scalar_u_entries"))
    output = parse_int(last_trace.get("scalar_u_output"))
    producer_target_u = parse_int(last_trace.get("producer_target_u_entries"))
    producer_state_rows = parse_int(last_trace.get("producer_state_rows"))
    producer_unique_state_rows = parse_int(
        last_trace.get("producer_unique_state_rows")
    )
    candidate_targets = parse_int(last_trace.get("producer_candidate_targets"))
    candidate_target_u = parse_int(
        last_trace.get("producer_candidate_target_u_entries")
    )
    underfilled_target_u = parse_int(
        last_trace.get("producer_underfilled_target_u_entries")
    )
    low_saved_stream_target_u = parse_int(
        last_trace.get("producer_low_saved_stream_target_u_entries")
    )
    active_catchup_attempts = parse_int(
        last_trace.get("producer_active_catchup_attempts")
    )
    active_catchup_deps = parse_int(last_trace.get("producer_active_catchup_deps"))
    active_catchup_targets = parse_int(
        last_trace.get("producer_active_catchup_targets")
    )
    compact_window_targets = parse_int(last_trace.get("compact_window_targets"))
    compact_window_target_u = parse_int(
        last_trace.get("compact_window_target_u_entries")
    )
    compact_window_claim_attempts = parse_int(
        last_trace.get("compact_window_claim_attempts")
    )
    compact_window_claims = parse_int(last_trace.get("compact_window_claims"))
    compact_window_claim_misses = parse_int(
        last_trace.get("compact_window_claim_misses")
    )
    compact_window_claim_stale_clears = parse_int(
        last_trace.get("compact_window_claim_stale_clears")
    )
    compact_window_claim_group_scatters = parse_int(
        last_trace.get("compact_window_claim_group_scatters")
    )
    compact_window_claim_group_detaches = parse_int(
        last_trace.get("compact_window_claim_group_detaches")
    )
    compact_window_claim_group_run_probes = parse_int(
        last_trace.get("compact_window_claim_group_run_probes")
    )
    compact_window_claim_group_run_states = parse_int(
        last_trace.get("compact_window_claim_group_run_states")
    )
    compact_window_claim_group_run_consecutive = parse_int(
        last_trace.get("compact_window_claim_group_run_consecutive")
    )
    compact_window_claim_group_run_near64 = parse_int(
        last_trace.get("compact_window_claim_group_run_near64")
    )
    compact_window_claim_group_run_near512 = parse_int(
        last_trace.get("compact_window_claim_group_run_near512")
    )
    compact_window_claim_group_run_max_consecutive = parse_int(
        last_trace.get("compact_window_claim_group_run_max_consecutive")
    )
    compact_window_claim_group_run_stealable = parse_int(
        last_trace.get("compact_window_claim_group_run_stealable")
    )
    compact_window_claim_group_run_max_stealable = parse_int(
        last_trace.get("compact_window_claim_group_run_max_stealable")
    )
    compact_window_claim_group_run_unreserved_near64 = parse_int(
        last_trace.get("compact_window_claim_group_run_unreserved_near64")
    )
    compact_window_claim_run_reservations = parse_int(
        last_trace.get("compact_window_claim_run_reservations")
    )
    compact_window_claim_run_rows = parse_int(
        last_trace.get("compact_window_claim_run_rows")
    )
    compact_window_claim_run_recomputes = parse_int(
        last_trace.get("compact_window_claim_run_recomputes")
    )
    compact_window_claim_run_updates = parse_int(
        last_trace.get("compact_window_claim_run_updates")
    )
    compact_window_claim_run_update_targets = parse_int(
        last_trace.get("compact_window_claim_run_update_targets")
    )
    compact_window_span_owner_reservations = parse_int(
        last_trace.get("compact_window_span_owner_reservations")
    )
    compact_window_span_owner_rows = parse_int(
        last_trace.get("compact_window_span_owner_rows")
    )
    compact_window_span_owner_owned_rows = parse_int(
        last_trace.get("compact_window_span_owner_owned_rows")
    )
    compact_window_span_owner_hole_rows = parse_int(
        last_trace.get("compact_window_span_owner_hole_rows")
    )
    compact_window_span_owner_deps = parse_int(
        last_trace.get("compact_window_span_owner_deps")
    )
    compact_window_span_owner_unique_deps = parse_int(
        last_trace.get("compact_window_span_owner_unique_deps")
    )
    compact_window_span_owner_scan_entries = parse_int(
        last_trace.get("compact_window_span_owner_scan_entries")
    )
    compact_window_span_owner_cols = parse_int(
        last_trace.get("compact_window_span_owner_cols")
    )
    compact_window_span_owner_slots = parse_int(
        last_trace.get("compact_window_span_owner_slots")
    )
    compact_window_span_owner_entries = parse_int(
        last_trace.get("compact_window_span_owner_entries")
    )
    compact_window_span_owner_hole_skips = parse_int(
        last_trace.get("compact_window_span_owner_hole_skips")
    )
    compact_window_span_owner_hole_skip_rows = parse_int(
        last_trace.get("compact_window_span_owner_hole_skip_rows")
    )
    compact_window_span_owner_prefix_shrinks = parse_int(
        last_trace.get("compact_window_span_owner_prefix_shrinks")
    )
    compact_window_span_owner_prefix_shrink_rows = parse_int(
        last_trace.get("compact_window_span_owner_prefix_shrink_rows")
    )
    compact_window_span_owner_oversize_skips = parse_int(
        last_trace.get("compact_window_span_owner_oversize_skips")
    )
    compact_window_span_owner_oversize_slots = parse_int(
        last_trace.get("compact_window_span_owner_oversize_slots")
    )
    compact_window_span_owner_link_skips = parse_int(
        last_trace.get("compact_window_span_owner_link_skips")
    )
    compact_window_span_owner_link_skip_entries = parse_int(
        last_trace.get("compact_window_span_owner_link_skip_entries")
    )
    compact_window_span_owner_scan_skips = parse_int(
        last_trace.get("compact_window_span_owner_scan_skips")
    )
    compact_window_span_owner_scan_skip_entries = parse_int(
        last_trace.get("compact_window_span_owner_scan_skip_entries")
    )
    compact_window_span_owner_payoff_skips = parse_int(
        last_trace.get("compact_window_span_owner_payoff_skips")
    )
    compact_window_span_owner_payoff_skip_rows = parse_int(
        last_trace.get("compact_window_span_owner_payoff_skip_rows")
    )
    compact_window_span_owner_payoff_skip_entries = parse_int(
        last_trace.get("compact_window_span_owner_payoff_skip_entries")
    )
    compact_window_span_owner_payoff_skip_scan_entries = parse_int(
        last_trace.get("compact_window_span_owner_payoff_skip_scan_entries")
    )
    compact_window_span_shared_owner_reservations = parse_int(
        last_trace.get("compact_window_span_shared_owner_reservations")
    )
    compact_window_span_shared_owner_rows = parse_int(
        last_trace.get("compact_window_span_shared_owner_rows")
    )
    compact_window_span_shared_owner_applies = parse_int(
        last_trace.get("compact_window_span_shared_owner_applies")
    )
    compact_window_span_shared_owner_entries = parse_int(
        last_trace.get("compact_window_span_shared_owner_entries")
    )
    compact_window_claim_run_output_surfaces = parse_int(
        last_trace.get("compact_window_claim_run_output_surfaces")
    )
    compact_window_claim_run_output_states = parse_int(
        last_trace.get("compact_window_claim_run_output_states")
    )
    compact_window_claim_run_output_deps = parse_int(
        last_trace.get("compact_window_claim_run_output_deps")
    )
    compact_window_claim_run_output_unique_deps = parse_int(
        last_trace.get("compact_window_claim_run_output_unique_deps")
    )
    compact_window_claim_run_output_duplicate_deps = parse_int(
        last_trace.get("compact_window_claim_run_output_duplicate_deps")
    )
    compact_window_claim_run_output_scan_entries = parse_int(
        last_trace.get("compact_window_claim_run_output_scan_entries")
    )
    compact_window_claim_run_output_group_scan_entries = parse_int(
        last_trace.get("compact_window_claim_run_output_group_scan_entries")
    )
    compact_window_claim_span_output_surfaces = parse_int(
        last_trace.get("compact_window_claim_span_output_surfaces")
    )
    compact_window_claim_span_output_states = parse_int(
        last_trace.get("compact_window_claim_span_output_states")
    )
    compact_window_claim_span_output_deps = parse_int(
        last_trace.get("compact_window_claim_span_output_deps")
    )
    compact_window_claim_span_output_unique_deps = parse_int(
        last_trace.get("compact_window_claim_span_output_unique_deps")
    )
    compact_window_claim_span_output_duplicate_deps = parse_int(
        last_trace.get("compact_window_claim_span_output_duplicate_deps")
    )
    compact_window_claim_span_output_scan_entries = parse_int(
        last_trace.get("compact_window_claim_span_output_scan_entries")
    )
    compact_window_claim_span_output_group_scan_entries = parse_int(
        last_trace.get("compact_window_claim_span_output_group_scan_entries")
    )
    compact_window_claim_span_output_value_entries = parse_int(
        last_trace.get("compact_window_claim_span_output_value_entries")
    )
    compact_window_claim_span_output_state_col_slots = parse_int(
        last_trace.get("compact_window_claim_span_output_state_col_slots")
    )
    compact_window_claim_span_output_unique_cols = parse_int(
        last_trace.get("compact_window_claim_span_output_unique_cols")
    )
    compact_window_claim_span_output_panel_slots = parse_int(
        last_trace.get("compact_window_claim_span_output_panel_slots")
    )
    compact_window_claim_span_output_dense_slots = parse_int(
        last_trace.get("compact_window_claim_span_output_dense_slots")
    )
    compact_window_claim_span_output_max_states = parse_int(
        last_trace.get("compact_window_claim_span_output_max_states")
    )
    compact_window_claim_span_output_max_unique_cols = parse_int(
        last_trace.get("compact_window_claim_span_output_max_unique_cols")
    )
    compact_window_delayed_output_skips = parse_int(
        last_trace.get("compact_window_delayed_output_skips")
    )
    compact_window_delayed_output_replays = parse_int(
        last_trace.get("compact_window_delayed_output_replays")
    )
    compact_window_delayed_output_deps = parse_int(
        last_trace.get("compact_window_delayed_output_deps")
    )
    compact_window_delayed_output_entries = parse_int(
        last_trace.get("compact_window_delayed_output_entries")
    )
    compact_window_delayed_output_scan_entries = parse_int(
        last_trace.get("compact_window_delayed_output_scan_entries")
    )
    compact_window_delayed_output_seek_skips = parse_int(
        last_trace.get("compact_window_delayed_output_seek_skips")
    )
    compact_window_deferred_output_cache_stores = parse_int(
        last_trace.get("compact_window_deferred_output_cache_stores")
    )
    compact_window_deferred_output_cache_applies = parse_int(
        last_trace.get("compact_window_deferred_output_cache_applies")
    )
    compact_window_deferred_output_cache_l_entries = parse_int(
        last_trace.get("compact_window_deferred_output_cache_l_entries")
    )
    compact_window_deferred_output_cache_entries = parse_int(
        last_trace.get("compact_window_deferred_output_cache_entries")
    )
    compact_window_deferred_output_cache_disables = parse_int(
        last_trace.get("compact_window_deferred_output_cache_disables")
    )
    compact_window_deferred_output_cache_overflows = parse_int(
        last_trace.get("compact_window_deferred_output_cache_overflows")
    )
    compact_window_delayed_group_replay_surfaces = parse_int(
        last_trace.get("compact_window_delayed_group_replay_surfaces")
    )
    compact_window_delayed_group_replay_states = parse_int(
        last_trace.get("compact_window_delayed_group_replay_states")
    )
    compact_window_delayed_group_replay_deps = parse_int(
        last_trace.get("compact_window_delayed_group_replay_deps")
    )
    compact_window_delayed_group_replay_unique_deps = parse_int(
        last_trace.get("compact_window_delayed_group_replay_unique_deps")
    )
    compact_window_delayed_group_replay_duplicate_deps = parse_int(
        last_trace.get("compact_window_delayed_group_replay_duplicate_deps")
    )
    compact_window_delayed_group_replay_scan_entries = parse_int(
        last_trace.get("compact_window_delayed_group_replay_scan_entries")
    )
    compact_window_delayed_group_replay_group_scan_entries = parse_int(
        last_trace.get("compact_window_delayed_group_replay_group_scan_entries")
    )
    compact_window_claim_group_output_surfaces = parse_int(
        last_trace.get("compact_window_claim_group_output_surfaces")
    )
    compact_window_claim_group_output_states = parse_int(
        last_trace.get("compact_window_claim_group_output_states")
    )
    compact_window_claim_group_output_deps = parse_int(
        last_trace.get("compact_window_claim_group_output_deps")
    )
    compact_window_claim_group_output_unique_deps = parse_int(
        last_trace.get("compact_window_claim_group_output_unique_deps")
    )
    compact_window_claim_group_output_duplicate_deps = parse_int(
        last_trace.get("compact_window_claim_group_output_duplicate_deps")
    )
    compact_window_claim_group_output_scan_entries = parse_int(
        last_trace.get("compact_window_claim_group_output_scan_entries")
    )
    compact_window_claim_group_output_group_scan_entries = parse_int(
        last_trace.get("compact_window_claim_group_output_group_scan_entries")
    )
    compact_window_batches = parse_int(last_trace.get("compact_window_batches"))
    compact_window_stream_u = parse_int(
        last_trace.get("compact_window_stream_u_entries")
    )
    compact_window_state_rows = parse_int(
        last_trace.get("compact_window_state_rows")
    )
    compact_window_unique_state_rows = parse_int(
        last_trace.get("compact_window_unique_state_rows")
    )
    compact_union_batches = parse_int(
        last_trace.get("compact_window_union_batches")
    )
    compact_union_targets = parse_int(
        last_trace.get("compact_window_union_targets")
    )
    compact_union_cols = parse_int(last_trace.get("compact_window_union_cols"))
    compact_union_values = parse_int(
        last_trace.get("compact_window_union_values")
    )
    compact_union_skips = parse_int(last_trace.get("compact_window_union_skips"))
    compact_union_skip_targets = parse_int(
        last_trace.get("compact_window_union_skip_targets")
    )
    compact_union_skip_cols = parse_int(
        last_trace.get("compact_window_union_skip_cols")
    )
    compact_union_skip_values = parse_int(
        last_trace.get("compact_window_union_skip_values")
    )
    compact_union_skip_sparse_values = parse_int(
        last_trace.get("compact_window_union_skip_sparse_values")
    )
    compact_group_steps = parse_int(last_trace.get("compact_window_group_steps"))
    compact_group_targets = parse_int(
        last_trace.get("compact_window_group_targets")
    )
    compact_group_cols = parse_int(last_trace.get("compact_window_group_cols"))
    compact_group_values = parse_int(
        last_trace.get("compact_window_group_values")
    )
    compact_group_merges = parse_int(
        last_trace.get("compact_window_group_merges")
    )
    compact_group_merge_targets = parse_int(
        last_trace.get("compact_window_group_merge_targets")
    )
    compact_group_merge_cols = parse_int(
        last_trace.get("compact_window_group_merge_cols")
    )
    compact_group_merge_values = parse_int(
        last_trace.get("compact_window_group_merge_values")
    )
    owner_surface_probes = parse_int(last_trace.get("owner_surface_probes"))
    owner_surface_probe_u = parse_int(
        last_trace.get("owner_surface_probe_u_entries")
    )
    owner_surface_scanned = parse_int(
        last_trace.get("owner_surface_scanned_rows")
    )
    owner_surface_targets = parse_int(last_trace.get("owner_surface_targets"))
    owner_surface_target_u = parse_int(
        last_trace.get("owner_surface_target_u_entries")
    )
    owner_surface_scanned_input = parse_int(
        last_trace.get("owner_surface_scanned_input_entries")
    )
    owner_surface_target_input = parse_int(
        last_trace.get("owner_surface_target_input_entries")
    )
    owner_surface_max_targets = parse_int(
        last_trace.get("owner_surface_max_targets")
    )
    owner_surface_max_target_input = parse_int(
        last_trace.get("owner_surface_max_target_input_entries")
    )
    panel_rows = parse_int(last_trace.get("panel_update_rows"))
    panel_entries = parse_int(last_trace.get("panel_update_entries"))
    summary: dict[str, object] = {
        "path": str(path),
        "trace_events": trace_event_count,
        "last_event": last_trace.get("event", ""),
        "completed": completed,
        "total": total,
        "scalar_u_entries": scalar,
        "scalar_u_output": output,
        "producer_target_u_entries": producer_target_u,
        "producer_state_rows": producer_state_rows,
        "producer_unique_state_rows": producer_unique_state_rows,
        "producer_candidate_targets": candidate_targets,
        "producer_candidate_target_u_entries": candidate_target_u,
        "producer_underfilled_target_u_entries": underfilled_target_u,
        "producer_low_saved_stream_target_u_entries": low_saved_stream_target_u,
        "producer_active_catchup_attempts": active_catchup_attempts,
        "producer_active_catchup_deps": active_catchup_deps,
        "producer_active_catchup_targets": active_catchup_targets,
        "compact_window_claim_attempts": compact_window_claim_attempts,
        "compact_window_claims": compact_window_claims,
        "compact_window_claim_misses": compact_window_claim_misses,
        "compact_window_claim_stale_clears": compact_window_claim_stale_clears,
        "compact_window_claim_group_scatters":
            compact_window_claim_group_scatters,
        "compact_window_claim_group_detaches":
            compact_window_claim_group_detaches,
        "compact_window_claim_group_run_probes":
            compact_window_claim_group_run_probes,
        "compact_window_claim_group_run_states":
            compact_window_claim_group_run_states,
        "compact_window_claim_group_run_consecutive":
            compact_window_claim_group_run_consecutive,
        "compact_window_claim_group_run_near64":
            compact_window_claim_group_run_near64,
        "compact_window_claim_group_run_near512":
            compact_window_claim_group_run_near512,
        "compact_window_claim_group_run_max_consecutive":
            compact_window_claim_group_run_max_consecutive,
        "compact_window_claim_group_run_stealable":
            compact_window_claim_group_run_stealable,
        "compact_window_claim_group_run_max_stealable":
            compact_window_claim_group_run_max_stealable,
        "compact_window_claim_group_run_unreserved_near64":
            compact_window_claim_group_run_unreserved_near64,
        "compact_window_claim_run_reservations":
            compact_window_claim_run_reservations,
        "compact_window_claim_run_rows":
            compact_window_claim_run_rows,
        "compact_window_claim_run_recomputes":
            compact_window_claim_run_recomputes,
        "compact_window_claim_run_updates":
            compact_window_claim_run_updates,
        "compact_window_claim_run_update_targets":
            compact_window_claim_run_update_targets,
        "compact_window_span_owner_reservations":
            compact_window_span_owner_reservations,
        "compact_window_span_owner_rows":
            compact_window_span_owner_rows,
        "compact_window_span_owner_owned_rows":
            compact_window_span_owner_owned_rows,
        "compact_window_span_owner_hole_rows":
            compact_window_span_owner_hole_rows,
        "compact_window_span_owner_deps":
            compact_window_span_owner_deps,
        "compact_window_span_owner_unique_deps":
            compact_window_span_owner_unique_deps,
        "compact_window_span_owner_scan_entries":
            compact_window_span_owner_scan_entries,
        "compact_window_span_owner_cols":
            compact_window_span_owner_cols,
        "compact_window_span_owner_slots":
            compact_window_span_owner_slots,
        "compact_window_span_owner_entries":
            compact_window_span_owner_entries,
        "compact_window_span_owner_hole_skips":
            compact_window_span_owner_hole_skips,
        "compact_window_span_owner_hole_skip_rows":
            compact_window_span_owner_hole_skip_rows,
        "compact_window_span_owner_prefix_shrinks":
            compact_window_span_owner_prefix_shrinks,
        "compact_window_span_owner_prefix_shrink_rows":
            compact_window_span_owner_prefix_shrink_rows,
        "compact_window_span_owner_oversize_skips":
            compact_window_span_owner_oversize_skips,
        "compact_window_span_owner_oversize_slots":
            compact_window_span_owner_oversize_slots,
        "compact_window_span_owner_link_skips":
            compact_window_span_owner_link_skips,
        "compact_window_span_owner_link_skip_entries":
            compact_window_span_owner_link_skip_entries,
        "compact_window_span_owner_scan_skips":
            compact_window_span_owner_scan_skips,
        "compact_window_span_owner_scan_skip_entries":
            compact_window_span_owner_scan_skip_entries,
        "compact_window_span_owner_payoff_skips":
            compact_window_span_owner_payoff_skips,
        "compact_window_span_owner_payoff_skip_rows":
            compact_window_span_owner_payoff_skip_rows,
        "compact_window_span_owner_payoff_skip_entries":
            compact_window_span_owner_payoff_skip_entries,
        "compact_window_span_owner_payoff_skip_scan_entries":
            compact_window_span_owner_payoff_skip_scan_entries,
        "compact_window_span_shared_owner_reservations":
            compact_window_span_shared_owner_reservations,
        "compact_window_span_shared_owner_rows":
            compact_window_span_shared_owner_rows,
        "compact_window_span_shared_owner_applies":
            compact_window_span_shared_owner_applies,
        "compact_window_span_shared_owner_entries":
            compact_window_span_shared_owner_entries,
        "compact_window_claim_run_output_surfaces":
            compact_window_claim_run_output_surfaces,
        "compact_window_claim_run_output_states":
            compact_window_claim_run_output_states,
        "compact_window_claim_run_output_deps":
            compact_window_claim_run_output_deps,
        "compact_window_claim_run_output_unique_deps":
            compact_window_claim_run_output_unique_deps,
        "compact_window_claim_run_output_duplicate_deps":
            compact_window_claim_run_output_duplicate_deps,
        "compact_window_claim_run_output_scan_entries":
            compact_window_claim_run_output_scan_entries,
        "compact_window_claim_run_output_group_scan_entries":
            compact_window_claim_run_output_group_scan_entries,
        "compact_window_claim_span_output_surfaces":
            compact_window_claim_span_output_surfaces,
        "compact_window_claim_span_output_states":
            compact_window_claim_span_output_states,
        "compact_window_claim_span_output_deps":
            compact_window_claim_span_output_deps,
        "compact_window_claim_span_output_unique_deps":
            compact_window_claim_span_output_unique_deps,
        "compact_window_claim_span_output_duplicate_deps":
            compact_window_claim_span_output_duplicate_deps,
        "compact_window_claim_span_output_scan_entries":
            compact_window_claim_span_output_scan_entries,
        "compact_window_claim_span_output_group_scan_entries":
            compact_window_claim_span_output_group_scan_entries,
        "compact_window_claim_span_output_value_entries":
            compact_window_claim_span_output_value_entries,
        "compact_window_claim_span_output_state_col_slots":
            compact_window_claim_span_output_state_col_slots,
        "compact_window_claim_span_output_unique_cols":
            compact_window_claim_span_output_unique_cols,
        "compact_window_claim_span_output_panel_slots":
            compact_window_claim_span_output_panel_slots,
        "compact_window_claim_span_output_dense_slots":
            compact_window_claim_span_output_dense_slots,
        "compact_window_claim_span_output_max_states":
            compact_window_claim_span_output_max_states,
        "compact_window_claim_span_output_max_unique_cols":
            compact_window_claim_span_output_max_unique_cols,
        "compact_window_delayed_output_skips":
            compact_window_delayed_output_skips,
        "compact_window_delayed_output_replays":
            compact_window_delayed_output_replays,
        "compact_window_delayed_output_deps": compact_window_delayed_output_deps,
        "compact_window_delayed_output_entries":
            compact_window_delayed_output_entries,
        "compact_window_delayed_output_scan_entries":
            compact_window_delayed_output_scan_entries,
        "compact_window_delayed_output_seek_skips":
            compact_window_delayed_output_seek_skips,
        "compact_window_deferred_output_cache_stores":
            compact_window_deferred_output_cache_stores,
        "compact_window_deferred_output_cache_applies":
            compact_window_deferred_output_cache_applies,
        "compact_window_deferred_output_cache_l_entries":
            compact_window_deferred_output_cache_l_entries,
        "compact_window_deferred_output_cache_entries":
            compact_window_deferred_output_cache_entries,
        "compact_window_deferred_output_cache_disables":
            compact_window_deferred_output_cache_disables,
        "compact_window_deferred_output_cache_overflows":
            compact_window_deferred_output_cache_overflows,
        "compact_window_delayed_group_replay_surfaces":
            compact_window_delayed_group_replay_surfaces,
        "compact_window_delayed_group_replay_states":
            compact_window_delayed_group_replay_states,
        "compact_window_delayed_group_replay_deps":
            compact_window_delayed_group_replay_deps,
        "compact_window_delayed_group_replay_unique_deps":
            compact_window_delayed_group_replay_unique_deps,
        "compact_window_delayed_group_replay_duplicate_deps":
            compact_window_delayed_group_replay_duplicate_deps,
        "compact_window_delayed_group_replay_scan_entries":
            compact_window_delayed_group_replay_scan_entries,
        "compact_window_delayed_group_replay_group_scan_entries":
            compact_window_delayed_group_replay_group_scan_entries,
        "compact_window_claim_group_output_surfaces":
            compact_window_claim_group_output_surfaces,
        "compact_window_claim_group_output_states":
            compact_window_claim_group_output_states,
        "compact_window_claim_group_output_deps":
            compact_window_claim_group_output_deps,
        "compact_window_claim_group_output_unique_deps":
            compact_window_claim_group_output_unique_deps,
        "compact_window_claim_group_output_duplicate_deps":
            compact_window_claim_group_output_duplicate_deps,
        "compact_window_claim_group_output_scan_entries":
            compact_window_claim_group_output_scan_entries,
        "compact_window_claim_group_output_group_scan_entries":
            compact_window_claim_group_output_group_scan_entries,
        "compact_window_targets": compact_window_targets,
        "compact_window_target_u_entries": compact_window_target_u,
        "compact_window_batches": compact_window_batches,
        "compact_window_stream_u_entries": compact_window_stream_u,
        "compact_window_state_rows": compact_window_state_rows,
        "compact_window_unique_state_rows": compact_window_unique_state_rows,
        "compact_window_union_batches": compact_union_batches,
        "compact_window_union_targets": compact_union_targets,
        "compact_window_union_cols": compact_union_cols,
        "compact_window_union_values": compact_union_values,
        "compact_window_union_skips": compact_union_skips,
        "compact_window_union_skip_targets": compact_union_skip_targets,
        "compact_window_union_skip_cols": compact_union_skip_cols,
        "compact_window_union_skip_values": compact_union_skip_values,
        "compact_window_union_skip_sparse_values":
            compact_union_skip_sparse_values,
        "compact_window_group_steps": compact_group_steps,
        "compact_window_group_targets": compact_group_targets,
        "compact_window_group_cols": compact_group_cols,
        "compact_window_group_values": compact_group_values,
        "compact_window_group_merges": compact_group_merges,
        "compact_window_group_merge_targets": compact_group_merge_targets,
        "compact_window_group_merge_cols": compact_group_merge_cols,
        "compact_window_group_merge_values": compact_group_merge_values,
        "owner_surface_probes": owner_surface_probes,
        "owner_surface_probe_u_entries": owner_surface_probe_u,
        "owner_surface_scanned_rows": owner_surface_scanned,
        "owner_surface_targets": owner_surface_targets,
        "owner_surface_target_u_entries": owner_surface_target_u,
        "owner_surface_scanned_input_entries": owner_surface_scanned_input,
        "owner_surface_target_input_entries": owner_surface_target_input,
        "owner_surface_max_targets": owner_surface_max_targets,
        "owner_surface_max_target_input_entries":
            owner_surface_max_target_input,
        "panel_update_rows": panel_rows,
        "panel_update_entries": panel_entries,
        "long_rows": long_rows,
        "long_max_row": long_max_row,
        "long_pivot_rows": long_pivot_rows,
        "long_scalar_u_entries": long_scalar,
        "long_scalar_u_internal": long_internal,
        "long_scalar_u_output": long_output,
        "long_producer_targets": long_producer_targets,
        "long_producer_target_u_entries": long_producer_target_u,
        "long_producer_state_rows": long_producer_state_rows,
        "long_producer_unique_state_rows": long_producer_unique_state_rows,
        "long_producer_candidate_targets": long_candidate_targets,
        "long_producer_candidate_target_u_entries": long_candidate_target_u,
        "long_producer_underfilled_target_u_entries": long_underfilled_target_u,
        "long_producer_low_saved_stream_target_u_entries":
            long_low_saved_stream_target_u,
        "long_producer_active_catchup_attempts": long_active_catchup_attempts,
        "long_producer_active_catchup_deps": long_active_catchup_deps,
        "long_producer_active_catchup_targets": long_active_catchup_targets,
        "long_compact_window_claim_attempts":
            long_compact_window_claim_attempts,
        "long_compact_window_claims": long_compact_window_claims,
        "long_compact_window_claim_misses": long_compact_window_claim_misses,
        "long_compact_window_claim_stale_clears":
            long_compact_window_claim_stale_clears,
        "long_compact_window_claim_group_scatters":
            long_compact_window_claim_group_scatters,
        "long_compact_window_claim_group_detaches":
            long_compact_window_claim_group_detaches,
        "long_compact_window_claim_group_run_probes":
            long_compact_window_claim_group_run_probes,
        "long_compact_window_claim_group_run_states":
            long_compact_window_claim_group_run_states,
        "long_compact_window_claim_group_run_consecutive":
            long_compact_window_claim_group_run_consecutive,
        "long_compact_window_claim_group_run_near64":
            long_compact_window_claim_group_run_near64,
        "long_compact_window_claim_group_run_near512":
            long_compact_window_claim_group_run_near512,
        "long_compact_window_claim_group_run_max_consecutive":
            long_compact_window_claim_group_run_max_consecutive,
        "long_compact_window_claim_group_run_stealable":
            long_compact_window_claim_group_run_stealable,
        "long_compact_window_claim_group_run_max_stealable":
            long_compact_window_claim_group_run_max_stealable,
        "long_compact_window_claim_group_run_unreserved_near64":
            long_compact_window_claim_group_run_unreserved_near64,
        "long_compact_window_claim_run_reservations":
            long_compact_window_claim_run_reservations,
        "long_compact_window_claim_run_rows":
            long_compact_window_claim_run_rows,
        "long_compact_window_claim_run_recomputes":
            long_compact_window_claim_run_recomputes,
        "long_compact_window_claim_run_updates":
            long_compact_window_claim_run_updates,
        "long_compact_window_claim_run_update_targets":
            long_compact_window_claim_run_update_targets,
        "long_compact_window_claim_run_output_surfaces":
            long_compact_window_claim_run_output_surfaces,
        "long_compact_window_claim_run_output_states":
            long_compact_window_claim_run_output_states,
        "long_compact_window_claim_run_output_deps":
            long_compact_window_claim_run_output_deps,
        "long_compact_window_claim_run_output_unique_deps":
            long_compact_window_claim_run_output_unique_deps,
        "long_compact_window_claim_run_output_duplicate_deps":
            long_compact_window_claim_run_output_duplicate_deps,
        "long_compact_window_claim_run_output_scan_entries":
            long_compact_window_claim_run_output_scan_entries,
        "long_compact_window_claim_run_output_group_scan_entries":
            long_compact_window_claim_run_output_group_scan_entries,
        "long_compact_window_claim_span_output_surfaces":
            long_compact_window_claim_span_output_surfaces,
        "long_compact_window_claim_span_output_states":
            long_compact_window_claim_span_output_states,
        "long_compact_window_claim_span_output_deps":
            long_compact_window_claim_span_output_deps,
        "long_compact_window_claim_span_output_unique_deps":
            long_compact_window_claim_span_output_unique_deps,
        "long_compact_window_claim_span_output_duplicate_deps":
            long_compact_window_claim_span_output_duplicate_deps,
        "long_compact_window_claim_span_output_scan_entries":
            long_compact_window_claim_span_output_scan_entries,
        "long_compact_window_claim_span_output_group_scan_entries":
            long_compact_window_claim_span_output_group_scan_entries,
        "long_compact_window_claim_span_output_value_entries":
            long_compact_window_claim_span_output_value_entries,
        "long_compact_window_claim_span_output_state_col_slots":
            long_compact_window_claim_span_output_state_col_slots,
        "long_compact_window_claim_span_output_unique_cols":
            long_compact_window_claim_span_output_unique_cols,
        "long_compact_window_claim_span_output_panel_slots":
            long_compact_window_claim_span_output_panel_slots,
        "long_compact_window_claim_span_output_dense_slots":
            long_compact_window_claim_span_output_dense_slots,
        "long_compact_window_claim_span_output_max_states":
            long_compact_window_claim_span_output_max_states,
        "long_compact_window_claim_span_output_max_unique_cols":
            long_compact_window_claim_span_output_max_unique_cols,
        "long_compact_window_delayed_output_skips":
            long_compact_window_delayed_output_skips,
        "long_compact_window_delayed_output_replays":
            long_compact_window_delayed_output_replays,
        "long_compact_window_delayed_output_deps":
            long_compact_window_delayed_output_deps,
        "long_compact_window_delayed_output_entries":
            long_compact_window_delayed_output_entries,
        "long_compact_window_delayed_output_scan_entries":
            long_compact_window_delayed_output_scan_entries,
        "long_compact_window_delayed_output_seek_skips":
            long_compact_window_delayed_output_seek_skips,
        "long_compact_window_delayed_group_replay_surfaces":
            long_compact_window_delayed_group_replay_surfaces,
        "long_compact_window_delayed_group_replay_states":
            long_compact_window_delayed_group_replay_states,
        "long_compact_window_delayed_group_replay_deps":
            long_compact_window_delayed_group_replay_deps,
        "long_compact_window_delayed_group_replay_unique_deps":
            long_compact_window_delayed_group_replay_unique_deps,
        "long_compact_window_delayed_group_replay_duplicate_deps":
            long_compact_window_delayed_group_replay_duplicate_deps,
        "long_compact_window_delayed_group_replay_scan_entries":
            long_compact_window_delayed_group_replay_scan_entries,
        "long_compact_window_delayed_group_replay_group_scan_entries":
            long_compact_window_delayed_group_replay_group_scan_entries,
        "long_compact_window_claim_group_output_surfaces":
            long_compact_window_claim_group_output_surfaces,
        "long_compact_window_claim_group_output_states":
            long_compact_window_claim_group_output_states,
        "long_compact_window_claim_group_output_deps":
            long_compact_window_claim_group_output_deps,
        "long_compact_window_claim_group_output_unique_deps":
            long_compact_window_claim_group_output_unique_deps,
        "long_compact_window_claim_group_output_duplicate_deps":
            long_compact_window_claim_group_output_duplicate_deps,
        "long_compact_window_claim_group_output_scan_entries":
            long_compact_window_claim_group_output_scan_entries,
        "long_compact_window_claim_group_output_group_scan_entries":
            long_compact_window_claim_group_output_group_scan_entries,
        "long_compact_window_targets": long_compact_window_targets,
        "long_compact_window_target_u_entries": long_compact_window_target_u,
        "long_compact_window_batches": long_compact_window_batches,
        "long_compact_window_stream_u_entries": long_compact_window_stream_u,
        "long_compact_window_state_rows": long_compact_window_state_rows,
        "long_compact_window_unique_state_rows":
            long_compact_window_unique_state_rows,
        "long_compact_window_union_batches": long_compact_union_batches,
        "long_compact_window_union_targets": long_compact_union_targets,
        "long_compact_window_union_cols": long_compact_union_cols,
        "long_compact_window_union_values": long_compact_union_values,
        "long_compact_window_union_skips": long_compact_union_skips,
        "long_compact_window_union_skip_targets":
            long_compact_union_skip_targets,
        "long_compact_window_union_skip_cols": long_compact_union_skip_cols,
        "long_compact_window_union_skip_values": long_compact_union_skip_values,
        "long_compact_window_union_skip_sparse_values":
            long_compact_union_skip_sparse_values,
        "long_compact_window_group_steps": long_compact_group_steps,
        "long_compact_window_group_targets": long_compact_group_targets,
        "long_compact_window_group_cols": long_compact_group_cols,
        "long_compact_window_group_values": long_compact_group_values,
        "long_compact_window_group_merges": long_compact_group_merges,
        "long_compact_window_group_merge_targets":
            long_compact_group_merge_targets,
        "long_compact_window_group_merge_cols": long_compact_group_merge_cols,
        "long_compact_window_group_merge_values":
            long_compact_group_merge_values,
        "long_owner_surface_probes": long_owner_surface_probes,
        "long_owner_surface_probe_u_entries": long_owner_surface_probe_u,
        "long_owner_surface_scanned_rows": long_owner_surface_scanned,
        "long_owner_surface_targets": long_owner_surface_targets,
        "long_owner_surface_target_u_entries": long_owner_surface_target_u,
        "long_owner_surface_scanned_input_entries":
            long_owner_surface_scanned_input,
        "long_owner_surface_target_input_entries":
            long_owner_surface_target_input,
        "long_owner_surface_max_targets": long_owner_surface_max_targets,
        "long_owner_surface_max_target_input_entries":
            long_owner_surface_max_target_input,
        "long_panel_update_rows": long_panel_rows,
        "long_panel_update_entries": long_panel_entries,
        "long_pivot_scalar_u_entries": long_pivot_scalar,
        "long_pivot_scalar_u_output": long_pivot_output,
        "long_pivot_producer_rows": long_pivot_producer_rows,
        "long_pivot_producer_targets": long_pivot_producer_targets,
        "long_pivot_producer_target_u_entries": long_pivot_producer_target_u,
        "long_pivot_probe_lookahead": long_pivot_probe_lookahead,
        "long_pivot_ready_roots": long_pivot_ready_roots,
        "long_pivot_underfilled": long_pivot_underfilled,
        "long_pivot_reject_bad_state": long_pivot_reject_bad_state,
        "long_pivot_reject_epoch": long_pivot_reject_epoch,
        "completion_fraction": ratio(completed, total),
        "scalar_u_output_share": ratio(output, scalar),
        "producer_target_u_per_scalar_u": ratio(producer_target_u, scalar),
        "producer_candidate_target_u_per_scalar_u": ratio(candidate_target_u, scalar),
        "scalar_u_output_per_producer_candidate_target_u": ratio(
            output, candidate_target_u
        ),
        "producer_accepted_target_u_share_of_candidates": ratio(
            producer_target_u, candidate_target_u
        ),
        "compact_window_target_u_per_scalar_u": ratio(
            compact_window_target_u, scalar
        ),
        "compact_window_claim_rate": ratio(
            compact_window_claims, compact_window_claim_attempts
        ),
        "compact_window_claim_miss_rate": ratio(
            compact_window_claim_misses, compact_window_claim_attempts
        ),
        "compact_window_claim_group_detach_rate": ratio(
            compact_window_claim_group_detaches,
            compact_window_claim_group_scatters
            + compact_window_claim_group_detaches,
        ),
        "compact_window_claim_group_run_states_per_probe": ratio(
            compact_window_claim_group_run_states,
            compact_window_claim_group_run_probes,
        ),
        "compact_window_claim_group_run_consecutive_per_probe": ratio(
            compact_window_claim_group_run_consecutive,
            compact_window_claim_group_run_probes,
        ),
        "compact_window_claim_group_run_consecutive_share": ratio(
            compact_window_claim_group_run_consecutive,
            compact_window_claim_group_run_states,
        ),
        "compact_window_claim_group_run_near64_per_probe": ratio(
            compact_window_claim_group_run_near64,
            compact_window_claim_group_run_probes,
        ),
        "compact_window_claim_group_run_near512_per_probe": ratio(
            compact_window_claim_group_run_near512,
            compact_window_claim_group_run_probes,
        ),
        "compact_window_claim_group_run_stealable_per_probe": ratio(
            compact_window_claim_group_run_stealable,
            compact_window_claim_group_run_probes,
        ),
        "compact_window_claim_group_run_stealable_share": ratio(
            compact_window_claim_group_run_stealable,
            compact_window_claim_group_run_consecutive,
        ),
        "compact_window_claim_group_run_unreserved_near64_share": ratio(
            compact_window_claim_group_run_unreserved_near64,
            compact_window_claim_group_run_near64,
        ),
        "compact_window_claim_run_rows_per_reservation": ratio(
            compact_window_claim_run_rows,
            compact_window_claim_run_reservations,
        ),
        "compact_window_claim_run_recompute_share": ratio(
            compact_window_claim_run_recomputes,
            compact_window_claim_run_rows,
        ),
        "compact_window_claim_run_updates_per_reservation": ratio(
            compact_window_claim_run_updates,
            compact_window_claim_run_reservations,
        ),
        "compact_window_claim_run_update_targets_per_update": ratio(
            compact_window_claim_run_update_targets,
            compact_window_claim_run_updates,
        ),
        "compact_window_span_owner_rows_per_reservation": ratio(
            compact_window_span_owner_rows,
            compact_window_span_owner_reservations,
        ),
        "compact_window_span_owner_owned_share": ratio(
            compact_window_span_owner_owned_rows,
            compact_window_span_owner_rows,
        ),
        "compact_window_span_owner_holes_per_reservation": ratio(
            compact_window_span_owner_hole_rows,
            compact_window_span_owner_reservations,
        ),
        "compact_window_span_owner_unique_dep_share": ratio(
            compact_window_span_owner_unique_deps,
            compact_window_span_owner_deps,
        ),
        "compact_window_span_owner_cols_per_reservation": ratio(
            compact_window_span_owner_cols,
            compact_window_span_owner_reservations,
        ),
        "compact_window_span_owner_entries_per_scan": ratio(
            compact_window_span_owner_entries,
            compact_window_span_owner_scan_entries,
        ),
        "compact_window_span_owner_entries_per_slot": ratio(
            compact_window_span_owner_entries,
            compact_window_span_owner_slots,
        ),
        "compact_window_span_owner_hole_skip_rows_per_skip": ratio(
            compact_window_span_owner_hole_skip_rows,
            compact_window_span_owner_hole_skips,
        ),
        "compact_window_span_owner_prefix_shrink_rows_per_shrink": ratio(
            compact_window_span_owner_prefix_shrink_rows,
            compact_window_span_owner_prefix_shrinks,
        ),
        "compact_window_span_owner_oversize_slots_per_skip": ratio(
            compact_window_span_owner_oversize_slots,
            compact_window_span_owner_oversize_skips,
        ),
        "compact_window_span_owner_link_skip_entries_per_skip": ratio(
            compact_window_span_owner_link_skip_entries,
            compact_window_span_owner_link_skips,
        ),
        "compact_window_span_owner_scan_skip_entries_per_skip": ratio(
            compact_window_span_owner_scan_skip_entries,
            compact_window_span_owner_scan_skips,
        ),
        "compact_window_span_owner_payoff_skip_rows_per_skip": ratio(
            compact_window_span_owner_payoff_skip_rows,
            compact_window_span_owner_payoff_skips,
        ),
        "compact_window_span_owner_payoff_skip_entries_per_scan": ratio(
            compact_window_span_owner_payoff_skip_entries,
            compact_window_span_owner_payoff_skip_scan_entries,
        ),
        "compact_window_span_shared_owner_rows_per_reservation": ratio(
            compact_window_span_shared_owner_rows,
            compact_window_span_shared_owner_reservations,
        ),
        "compact_window_span_shared_owner_entries_per_apply": ratio(
            compact_window_span_shared_owner_entries,
            compact_window_span_shared_owner_applies,
        ),
        "compact_window_claim_run_output_states_per_surface": ratio(
            compact_window_claim_run_output_states,
            compact_window_claim_run_output_surfaces,
        ),
        "compact_window_claim_run_output_unique_dep_share": ratio(
            compact_window_claim_run_output_unique_deps,
            compact_window_claim_run_output_deps,
        ),
        "compact_window_claim_run_output_group_scan_ratio": ratio(
            compact_window_claim_run_output_group_scan_entries,
            compact_window_claim_run_output_scan_entries,
        ),
        "compact_window_claim_span_output_states_per_surface": ratio(
            compact_window_claim_span_output_states,
            compact_window_claim_span_output_surfaces,
        ),
        "compact_window_claim_span_output_unique_dep_share": ratio(
            compact_window_claim_span_output_unique_deps,
            compact_window_claim_span_output_deps,
        ),
        "compact_window_claim_span_output_group_scan_ratio": ratio(
            compact_window_claim_span_output_group_scan_entries,
            compact_window_claim_span_output_scan_entries,
        ),
        "compact_window_claim_span_output_values_per_surface": ratio(
            compact_window_claim_span_output_value_entries,
            compact_window_claim_span_output_surfaces,
        ),
        "compact_window_claim_span_output_state_col_slots_per_surface": ratio(
            compact_window_claim_span_output_state_col_slots,
            compact_window_claim_span_output_surfaces,
        ),
        "compact_window_claim_span_output_unique_cols_per_surface": ratio(
            compact_window_claim_span_output_unique_cols,
            compact_window_claim_span_output_surfaces,
        ),
        "compact_window_claim_span_output_values_per_state_col_slot": ratio(
            compact_window_claim_span_output_value_entries,
            compact_window_claim_span_output_state_col_slots,
        ),
        "compact_window_claim_span_output_panel_slots_per_state_col_slot": ratio(
            compact_window_claim_span_output_panel_slots,
            compact_window_claim_span_output_state_col_slots,
        ),
        "compact_window_claim_span_output_dense_slots_per_state_col_slot": ratio(
            compact_window_claim_span_output_dense_slots,
            compact_window_claim_span_output_state_col_slots,
        ),
        "compact_window_claim_span_output_panel_slots_per_value": ratio(
            compact_window_claim_span_output_panel_slots,
            compact_window_claim_span_output_value_entries,
        ),
        "compact_window_claim_span_output_dense_slots_per_value": ratio(
            compact_window_claim_span_output_dense_slots,
            compact_window_claim_span_output_value_entries,
        ),
        "compact_window_delayed_output_replay_entries_per_skip": ratio(
            compact_window_delayed_output_entries,
            compact_window_delayed_output_skips,
        ),
        "compact_window_delayed_output_deps_per_replay": ratio(
            compact_window_delayed_output_deps,
            compact_window_delayed_output_replays,
        ),
        "compact_window_delayed_output_entries_per_scan": ratio(
            compact_window_delayed_output_entries,
            compact_window_delayed_output_scan_entries,
        ),
        "compact_window_delayed_output_seek_skip_share": ratio(
            compact_window_delayed_output_seek_skips,
            compact_window_delayed_output_scan_entries
            + compact_window_delayed_output_seek_skips,
        ),
        "compact_window_deferred_output_cache_entries_per_apply": ratio(
            compact_window_deferred_output_cache_entries,
            compact_window_deferred_output_cache_applies,
        ),
        "compact_window_deferred_output_cache_store_share": ratio(
            compact_window_deferred_output_cache_stores,
            compact_window_delayed_output_skips,
        ),
        "compact_window_delayed_group_replay_states_per_surface": ratio(
            compact_window_delayed_group_replay_states,
            compact_window_delayed_group_replay_surfaces,
        ),
        "compact_window_delayed_group_replay_unique_dep_share": ratio(
            compact_window_delayed_group_replay_unique_deps,
            compact_window_delayed_group_replay_deps,
        ),
        "compact_window_delayed_group_replay_group_scan_ratio": ratio(
            compact_window_delayed_group_replay_group_scan_entries,
            compact_window_delayed_group_replay_scan_entries,
        ),
        "compact_window_claim_group_output_states_per_surface": ratio(
            compact_window_claim_group_output_states,
            compact_window_claim_group_output_surfaces,
        ),
        "compact_window_claim_group_output_unique_dep_share": ratio(
            compact_window_claim_group_output_unique_deps,
            compact_window_claim_group_output_deps,
        ),
        "compact_window_claim_group_output_group_scan_ratio": ratio(
            compact_window_claim_group_output_group_scan_entries,
            compact_window_claim_group_output_scan_entries,
        ),
        "compact_window_target_u_per_stream_u": ratio(
            compact_window_target_u, compact_window_stream_u
        ),
        "scalar_u_output_per_compact_window_target_u": ratio(
            output, compact_window_target_u
        ),
        "compact_window_state_rows_per_unique": ratio(
            compact_window_state_rows, compact_window_unique_state_rows
        ),
        "compact_window_union_values_per_col": ratio(
            compact_union_values, compact_union_cols
        ),
        "compact_window_union_targets_per_batch": ratio(
            compact_union_targets, compact_union_batches
        ),
        "compact_window_union_skip_targets_per_skip": ratio(
            compact_union_skip_targets, compact_union_skips
        ),
        "compact_window_union_skip_values_per_sparse_value": ratio(
            compact_union_skip_values, compact_union_skip_sparse_values
        ),
        "compact_window_group_values_per_col": ratio(
            compact_group_values, compact_group_cols
        ),
        "compact_window_group_targets_per_step": ratio(
            compact_group_targets, compact_group_steps
        ),
        "compact_window_group_merge_targets_per_merge": ratio(
            compact_group_merge_targets, compact_group_merges
        ),
        "compact_window_group_merge_values_per_col": ratio(
            compact_group_merge_values, compact_group_merge_cols
        ),
        "owner_surface_targets_per_probe": ratio(
            owner_surface_targets, owner_surface_probes
        ),
        "owner_surface_targets_per_scanned_row": ratio(
            owner_surface_targets, owner_surface_scanned
        ),
        "owner_surface_target_u_per_probe_u": ratio(
            owner_surface_target_u, owner_surface_probe_u
        ),
        "owner_surface_target_input_per_scanned_input": ratio(
            owner_surface_target_input, owner_surface_scanned_input
        ),
        "owner_surface_target_input_per_target": ratio(
            owner_surface_target_input, owner_surface_targets
        ),
        "panel_update_entries_per_scalar_u": ratio(panel_entries, scalar),
        "scalar_u_output_per_panel_update_entries": ratio(output, panel_entries),
        "scalar_u_output_per_producer_target_u": ratio(
            output, producer_target_u
        ),
        "producer_state_rows_per_unique": ratio(
            producer_state_rows, producer_unique_state_rows
        ),
        "long_pivot_share": ratio(long_pivot_rows, long_rows),
        "long_scalar_u_output_share": ratio(long_output, long_scalar),
        "long_producer_target_u_per_scalar_u": ratio(
            long_producer_target_u, long_scalar
        ),
        "long_producer_candidate_target_u_per_scalar_u": ratio(
            long_candidate_target_u, long_scalar
        ),
        "long_compact_window_target_u_per_scalar_u": ratio(
            long_compact_window_target_u, long_scalar
        ),
        "long_compact_window_claim_rate": ratio(
            long_compact_window_claims, long_compact_window_claim_attempts
        ),
        "long_compact_window_claim_miss_rate": ratio(
            long_compact_window_claim_misses,
            long_compact_window_claim_attempts,
        ),
        "long_compact_window_claim_group_detach_rate": ratio(
            long_compact_window_claim_group_detaches,
            long_compact_window_claim_group_scatters
            + long_compact_window_claim_group_detaches,
        ),
        "long_compact_window_claim_group_run_states_per_probe": ratio(
            long_compact_window_claim_group_run_states,
            long_compact_window_claim_group_run_probes,
        ),
        "long_compact_window_claim_group_run_consecutive_per_probe": ratio(
            long_compact_window_claim_group_run_consecutive,
            long_compact_window_claim_group_run_probes,
        ),
        "long_compact_window_claim_group_run_consecutive_share": ratio(
            long_compact_window_claim_group_run_consecutive,
            long_compact_window_claim_group_run_states,
        ),
        "long_compact_window_claim_group_run_near64_per_probe": ratio(
            long_compact_window_claim_group_run_near64,
            long_compact_window_claim_group_run_probes,
        ),
        "long_compact_window_claim_group_run_near512_per_probe": ratio(
            long_compact_window_claim_group_run_near512,
            long_compact_window_claim_group_run_probes,
        ),
        "long_compact_window_claim_group_run_stealable_per_probe": ratio(
            long_compact_window_claim_group_run_stealable,
            long_compact_window_claim_group_run_probes,
        ),
        "long_compact_window_claim_group_run_stealable_share": ratio(
            long_compact_window_claim_group_run_stealable,
            long_compact_window_claim_group_run_consecutive,
        ),
        "long_compact_window_claim_group_run_unreserved_near64_share": ratio(
            long_compact_window_claim_group_run_unreserved_near64,
            long_compact_window_claim_group_run_near64,
        ),
        "long_compact_window_claim_run_rows_per_reservation": ratio(
            long_compact_window_claim_run_rows,
            long_compact_window_claim_run_reservations,
        ),
        "long_compact_window_claim_run_recompute_share": ratio(
            long_compact_window_claim_run_recomputes,
            long_compact_window_claim_run_rows,
        ),
        "long_compact_window_claim_run_updates_per_reservation": ratio(
            long_compact_window_claim_run_updates,
            long_compact_window_claim_run_reservations,
        ),
        "long_compact_window_claim_run_update_targets_per_update": ratio(
            long_compact_window_claim_run_update_targets,
            long_compact_window_claim_run_updates,
        ),
        "long_compact_window_delayed_output_replay_entries_per_skip": ratio(
            long_compact_window_delayed_output_entries,
            long_compact_window_delayed_output_skips,
        ),
        "long_compact_window_delayed_output_deps_per_replay": ratio(
            long_compact_window_delayed_output_deps,
            long_compact_window_delayed_output_replays,
        ),
        "long_compact_window_delayed_output_entries_per_scan": ratio(
            long_compact_window_delayed_output_entries,
            long_compact_window_delayed_output_scan_entries,
        ),
        "long_compact_window_delayed_output_seek_skip_share": ratio(
            long_compact_window_delayed_output_seek_skips,
            long_compact_window_delayed_output_scan_entries
            + long_compact_window_delayed_output_seek_skips,
        ),
        "long_compact_window_delayed_group_replay_states_per_surface": ratio(
            long_compact_window_delayed_group_replay_states,
            long_compact_window_delayed_group_replay_surfaces,
        ),
        "long_compact_window_delayed_group_replay_unique_dep_share": ratio(
            long_compact_window_delayed_group_replay_unique_deps,
            long_compact_window_delayed_group_replay_deps,
        ),
        "long_compact_window_delayed_group_replay_group_scan_ratio": ratio(
            long_compact_window_delayed_group_replay_group_scan_entries,
            long_compact_window_delayed_group_replay_scan_entries,
        ),
        "long_compact_window_claim_group_output_states_per_surface": ratio(
            long_compact_window_claim_group_output_states,
            long_compact_window_claim_group_output_surfaces,
        ),
        "long_compact_window_claim_group_output_unique_dep_share": ratio(
            long_compact_window_claim_group_output_unique_deps,
            long_compact_window_claim_group_output_deps,
        ),
        "long_compact_window_claim_group_output_group_scan_ratio": ratio(
            long_compact_window_claim_group_output_group_scan_entries,
            long_compact_window_claim_group_output_scan_entries,
        ),
        "long_compact_window_target_u_per_stream_u": ratio(
            long_compact_window_target_u, long_compact_window_stream_u
        ),
        "long_compact_window_state_rows_per_unique": ratio(
            long_compact_window_state_rows,
            long_compact_window_unique_state_rows,
        ),
        "long_compact_window_union_values_per_col": ratio(
            long_compact_union_values, long_compact_union_cols
        ),
        "long_compact_window_union_targets_per_batch": ratio(
            long_compact_union_targets, long_compact_union_batches
        ),
        "long_compact_window_union_skip_targets_per_skip": ratio(
            long_compact_union_skip_targets, long_compact_union_skips
        ),
        "long_compact_window_union_skip_values_per_sparse_value": ratio(
            long_compact_union_skip_values,
            long_compact_union_skip_sparse_values,
        ),
        "long_compact_window_group_values_per_col": ratio(
            long_compact_group_values, long_compact_group_cols
        ),
        "long_compact_window_group_targets_per_step": ratio(
            long_compact_group_targets, long_compact_group_steps
        ),
        "long_compact_window_group_merge_targets_per_merge": ratio(
            long_compact_group_merge_targets, long_compact_group_merges
        ),
        "long_compact_window_group_merge_values_per_col": ratio(
            long_compact_group_merge_values, long_compact_group_merge_cols
        ),
        "long_owner_surface_targets_per_probe": ratio(
            long_owner_surface_targets, long_owner_surface_probes
        ),
        "long_owner_surface_targets_per_scanned_row": ratio(
            long_owner_surface_targets, long_owner_surface_scanned
        ),
        "long_owner_surface_target_u_per_probe_u": ratio(
            long_owner_surface_target_u, long_owner_surface_probe_u
        ),
        "long_owner_surface_target_input_per_scanned_input": ratio(
            long_owner_surface_target_input,
            long_owner_surface_scanned_input,
        ),
        "long_owner_surface_target_input_per_target": ratio(
            long_owner_surface_target_input, long_owner_surface_targets
        ),
        "long_scalar_u_output_per_producer_target_u": ratio(
            long_output, long_producer_target_u
        ),
        "long_panel_update_entries_per_scalar_u": ratio(
            long_panel_entries, long_scalar
        ),
        "long_scalar_u_output_per_panel_update_entries": ratio(
            long_output, long_panel_entries
        ),
        "long_producer_state_rows_per_unique": ratio(
            long_producer_state_rows, long_producer_unique_state_rows
        ),
        "long_pivot_producer_row_share": ratio(
            long_pivot_producer_rows, long_pivot_rows
        ),
        "long_pivot_producer_target_u_per_scalar_u": ratio(
            long_pivot_producer_target_u, long_pivot_scalar
        ),
        "long_pivot_scalar_u_output_per_producer_target_u": ratio(
            long_pivot_output, long_pivot_producer_target_u
        ),
    }
    for name in PRODUCER_COUNTERS:
        summary[name] = parse_int(last_trace.get(name))
    return summary


def print_table(rows: Iterable[dict[str, object]]) -> None:
    columns = [
        "path",
        "last_event",
        "completed",
        "total",
        "scalar_u_entries",
        "scalar_u_output",
        "producer_batches",
        "producer_targets",
        "producer_target_u_entries",
        "producer_state_rows",
        "producer_unique_state_rows",
        "producer_state_rows_per_unique",
        "producer_candidate_targets",
        "producer_candidate_target_u_entries",
        "producer_candidate_target_u_per_scalar_u",
        "scalar_u_output_per_producer_candidate_target_u",
        "producer_accepted_target_u_share_of_candidates",
        "producer_underfilled_target_u_entries",
        "producer_low_saved_stream_target_u_entries",
        "producer_active_catchup_attempts",
        "producer_active_catchup_deps",
        "producer_active_catchup_targets",
        "compact_window_fills",
        "compact_window_evictions",
        "compact_window_overflows",
        "compact_window_claim_attempts",
        "compact_window_claims",
        "compact_window_claim_misses",
        "compact_window_claim_stale_clears",
        "compact_window_claim_group_scatters",
        "compact_window_claim_rate",
        "compact_window_claim_miss_rate",
        "compact_window_claim_group_detaches",
        "compact_window_claim_group_detach_rate",
        "compact_window_claim_group_run_probes",
        "compact_window_claim_group_run_states",
        "compact_window_claim_group_run_states_per_probe",
        "compact_window_claim_group_run_consecutive",
        "compact_window_claim_group_run_consecutive_per_probe",
        "compact_window_claim_group_run_consecutive_share",
        "compact_window_claim_group_run_near64",
        "compact_window_claim_group_run_near64_per_probe",
        "compact_window_claim_group_run_near512",
        "compact_window_claim_group_run_near512_per_probe",
        "compact_window_claim_group_run_max_consecutive",
        "compact_window_delayed_output_skips",
        "compact_window_delayed_output_replays",
        "compact_window_delayed_output_deps",
        "compact_window_delayed_output_entries",
        "compact_window_delayed_output_scan_entries",
        "compact_window_delayed_output_seek_skips",
        "compact_window_delayed_output_replay_entries_per_skip",
        "compact_window_delayed_output_deps_per_replay",
        "compact_window_delayed_output_entries_per_scan",
        "compact_window_delayed_output_seek_skip_share",
        "compact_window_deferred_output_cache_stores",
        "compact_window_deferred_output_cache_applies",
        "compact_window_deferred_output_cache_l_entries",
        "compact_window_deferred_output_cache_entries",
        "compact_window_deferred_output_cache_disables",
        "compact_window_deferred_output_cache_overflows",
        "compact_window_deferred_output_cache_entries_per_apply",
        "compact_window_deferred_output_cache_store_share",
        "compact_window_delayed_group_replay_surfaces",
        "compact_window_delayed_group_replay_states",
        "compact_window_delayed_group_replay_states_per_surface",
        "compact_window_delayed_group_replay_deps",
        "compact_window_delayed_group_replay_unique_deps",
        "compact_window_delayed_group_replay_duplicate_deps",
        "compact_window_delayed_group_replay_unique_dep_share",
        "compact_window_delayed_group_replay_scan_entries",
        "compact_window_delayed_group_replay_group_scan_entries",
        "compact_window_delayed_group_replay_group_scan_ratio",
        "compact_window_claim_span_output_surfaces",
        "compact_window_claim_span_output_states",
        "compact_window_claim_span_output_states_per_surface",
        "compact_window_claim_span_output_deps",
        "compact_window_claim_span_output_unique_deps",
        "compact_window_claim_span_output_unique_dep_share",
        "compact_window_claim_span_output_scan_entries",
        "compact_window_claim_span_output_group_scan_entries",
        "compact_window_claim_span_output_group_scan_ratio",
        "compact_window_claim_span_output_value_entries",
        "compact_window_claim_span_output_state_col_slots",
        "compact_window_claim_span_output_unique_cols",
        "compact_window_claim_span_output_panel_slots",
        "compact_window_claim_span_output_dense_slots",
        "compact_window_claim_span_output_max_states",
        "compact_window_claim_span_output_max_unique_cols",
        "compact_window_claim_span_output_values_per_surface",
        "compact_window_claim_span_output_state_col_slots_per_surface",
        "compact_window_claim_span_output_unique_cols_per_surface",
        "compact_window_claim_span_output_values_per_state_col_slot",
        "compact_window_claim_span_output_panel_slots_per_state_col_slot",
        "compact_window_claim_span_output_dense_slots_per_state_col_slot",
        "compact_window_claim_span_output_panel_slots_per_value",
        "compact_window_claim_span_output_dense_slots_per_value",
        "compact_window_span_owner_reservations",
        "compact_window_span_owner_rows",
        "compact_window_span_owner_rows_per_reservation",
        "compact_window_span_owner_owned_rows",
        "compact_window_span_owner_hole_rows",
        "compact_window_span_owner_owned_share",
        "compact_window_span_owner_holes_per_reservation",
        "compact_window_span_owner_unique_dep_share",
        "compact_window_span_owner_scan_entries",
        "compact_window_span_owner_cols",
        "compact_window_span_owner_cols_per_reservation",
        "compact_window_span_owner_slots",
        "compact_window_span_owner_entries",
        "compact_window_span_owner_entries_per_scan",
        "compact_window_span_owner_entries_per_slot",
        "compact_window_span_owner_hole_skips",
        "compact_window_span_owner_hole_skip_rows_per_skip",
        "compact_window_span_owner_prefix_shrinks",
        "compact_window_span_owner_prefix_shrink_rows_per_shrink",
        "compact_window_span_owner_oversize_skips",
        "compact_window_span_owner_oversize_slots_per_skip",
        "compact_window_span_owner_link_skips",
        "compact_window_span_owner_link_skip_entries_per_skip",
        "compact_window_span_owner_scan_skips",
        "compact_window_span_owner_scan_skip_entries_per_skip",
        "compact_window_span_owner_payoff_skips",
        "compact_window_span_owner_payoff_skip_rows_per_skip",
        "compact_window_span_owner_payoff_skip_entries_per_scan",
        "compact_window_span_shared_owner_reservations",
        "compact_window_span_shared_owner_rows_per_reservation",
        "compact_window_span_shared_owner_applies",
        "compact_window_span_shared_owner_entries_per_apply",
        "compact_window_batches",
        "compact_window_stream_u_entries",
        "compact_window_targets",
        "compact_window_target_u_entries",
        "compact_window_target_u_per_stream_u",
        "compact_window_target_u_per_scalar_u",
        "scalar_u_output_per_compact_window_target_u",
        "compact_window_state_rows",
        "compact_window_unique_state_rows",
        "compact_window_state_rows_per_unique",
        "compact_window_union_batches",
        "compact_window_union_targets",
        "compact_window_union_cols",
        "compact_window_union_values",
        "compact_window_union_values_per_col",
        "compact_window_union_targets_per_batch",
        "compact_window_union_skips",
        "compact_window_union_skip_targets",
        "compact_window_union_skip_cols",
        "compact_window_union_skip_values",
        "compact_window_union_skip_sparse_values",
        "compact_window_union_skip_targets_per_skip",
        "compact_window_union_skip_values_per_sparse_value",
        "compact_window_group_steps",
        "compact_window_group_targets",
        "compact_window_group_cols",
        "compact_window_group_values",
        "compact_window_group_values_per_col",
        "compact_window_group_targets_per_step",
        "compact_window_group_merges",
        "compact_window_group_merge_targets",
        "compact_window_group_merge_cols",
        "compact_window_group_merge_values",
        "compact_window_group_merge_targets_per_merge",
        "compact_window_group_merge_values_per_col",
        "owner_surface_probes",
        "owner_surface_probe_u_entries",
        "owner_surface_scanned_rows",
        "owner_surface_targets",
        "owner_surface_target_u_entries",
        "owner_surface_scanned_input_entries",
        "owner_surface_target_input_entries",
        "owner_surface_targets_per_probe",
        "owner_surface_targets_per_scanned_row",
        "owner_surface_target_u_per_probe_u",
        "owner_surface_target_input_per_scanned_input",
        "owner_surface_target_input_per_target",
        "owner_surface_max_targets",
        "owner_surface_max_target_input_entries",
        "producer_underfilled",
        "producer_low_saved_stream",
        "producer_reject_bad_state",
        "producer_reject_epoch",
        "producer_reject_dep_absent",
        "producer_reject_not_root",
        "panel_update_rows",
        "panel_update_entries",
        "panel_update_entries_per_scalar_u",
        "scalar_u_output_per_panel_update_entries",
        "long_rows",
        "long_max_row",
        "long_pivot_rows",
        "long_pivot_producer_rows",
        "long_scalar_u_entries",
        "long_scalar_u_output",
        "long_producer_target_u_entries",
        "long_producer_state_rows",
        "long_producer_unique_state_rows",
        "long_producer_state_rows_per_unique",
        "long_producer_candidate_targets",
        "long_producer_candidate_target_u_entries",
        "long_producer_candidate_target_u_per_scalar_u",
        "long_producer_underfilled_target_u_entries",
        "long_producer_low_saved_stream_target_u_entries",
        "long_producer_active_catchup_attempts",
        "long_producer_active_catchup_deps",
        "long_producer_active_catchup_targets",
        "long_compact_window_claim_attempts",
        "long_compact_window_claims",
        "long_compact_window_claim_misses",
        "long_compact_window_claim_stale_clears",
        "long_compact_window_claim_group_scatters",
        "long_compact_window_claim_rate",
        "long_compact_window_claim_miss_rate",
        "long_compact_window_claim_group_detaches",
        "long_compact_window_claim_group_detach_rate",
        "long_compact_window_claim_group_run_probes",
        "long_compact_window_claim_group_run_states",
        "long_compact_window_claim_group_run_states_per_probe",
        "long_compact_window_claim_group_run_consecutive",
        "long_compact_window_claim_group_run_consecutive_per_probe",
        "long_compact_window_claim_group_run_consecutive_share",
        "long_compact_window_claim_group_run_near64",
        "long_compact_window_claim_group_run_near64_per_probe",
        "long_compact_window_claim_group_run_near512",
        "long_compact_window_claim_group_run_near512_per_probe",
        "long_compact_window_claim_group_run_max_consecutive",
        "long_compact_window_delayed_output_skips",
        "long_compact_window_delayed_output_replays",
        "long_compact_window_delayed_output_deps",
        "long_compact_window_delayed_output_entries",
        "long_compact_window_delayed_output_scan_entries",
        "long_compact_window_delayed_output_seek_skips",
        "long_compact_window_delayed_output_replay_entries_per_skip",
        "long_compact_window_delayed_output_deps_per_replay",
        "long_compact_window_delayed_output_entries_per_scan",
        "long_compact_window_delayed_output_seek_skip_share",
        "long_compact_window_delayed_group_replay_surfaces",
        "long_compact_window_delayed_group_replay_states",
        "long_compact_window_delayed_group_replay_states_per_surface",
        "long_compact_window_delayed_group_replay_deps",
        "long_compact_window_delayed_group_replay_unique_deps",
        "long_compact_window_delayed_group_replay_duplicate_deps",
        "long_compact_window_delayed_group_replay_unique_dep_share",
        "long_compact_window_delayed_group_replay_scan_entries",
        "long_compact_window_delayed_group_replay_group_scan_entries",
        "long_compact_window_delayed_group_replay_group_scan_ratio",
        "long_compact_window_targets",
        "long_compact_window_target_u_entries",
        "long_compact_window_batches",
        "long_compact_window_stream_u_entries",
        "long_compact_window_target_u_per_stream_u",
        "long_compact_window_target_u_per_scalar_u",
        "long_compact_window_state_rows",
        "long_compact_window_unique_state_rows",
        "long_compact_window_state_rows_per_unique",
        "long_compact_window_union_batches",
        "long_compact_window_union_targets",
        "long_compact_window_union_cols",
        "long_compact_window_union_values",
        "long_compact_window_union_values_per_col",
        "long_compact_window_union_targets_per_batch",
        "long_compact_window_union_skips",
        "long_compact_window_union_skip_targets",
        "long_compact_window_union_skip_cols",
        "long_compact_window_union_skip_values",
        "long_compact_window_union_skip_sparse_values",
        "long_compact_window_union_skip_targets_per_skip",
        "long_compact_window_union_skip_values_per_sparse_value",
        "long_compact_window_group_steps",
        "long_compact_window_group_targets",
        "long_compact_window_group_cols",
        "long_compact_window_group_values",
        "long_compact_window_group_values_per_col",
        "long_compact_window_group_targets_per_step",
        "long_compact_window_group_merges",
        "long_compact_window_group_merge_targets",
        "long_compact_window_group_merge_cols",
        "long_compact_window_group_merge_values",
        "long_compact_window_group_merge_targets_per_merge",
        "long_compact_window_group_merge_values_per_col",
        "long_owner_surface_probes",
        "long_owner_surface_probe_u_entries",
        "long_owner_surface_scanned_rows",
        "long_owner_surface_targets",
        "long_owner_surface_target_u_entries",
        "long_owner_surface_scanned_input_entries",
        "long_owner_surface_target_input_entries",
        "long_owner_surface_targets_per_probe",
        "long_owner_surface_targets_per_scanned_row",
        "long_owner_surface_target_u_per_probe_u",
        "long_owner_surface_target_input_per_scanned_input",
        "long_owner_surface_target_input_per_target",
        "long_owner_surface_max_targets",
        "long_owner_surface_max_target_input_entries",
        "long_panel_update_rows",
        "long_panel_update_entries",
        "long_panel_update_entries_per_scalar_u",
        "long_scalar_u_output_per_panel_update_entries",
        "long_pivot_scalar_u_entries",
        "long_pivot_scalar_u_output",
        "long_pivot_producer_target_u_entries",
        "long_pivot_probe_lookahead",
        "long_pivot_ready_roots",
        "long_pivot_underfilled",
        "long_pivot_reject_bad_state",
        "long_pivot_reject_epoch",
        "completion_fraction",
        "scalar_u_output_share",
        "producer_target_u_per_scalar_u",
        "scalar_u_output_per_producer_target_u",
        "long_pivot_share",
        "long_pivot_producer_row_share",
        "long_scalar_u_output_share",
        "long_producer_target_u_per_scalar_u",
        "long_scalar_u_output_per_producer_target_u",
        "long_pivot_producer_target_u_per_scalar_u",
        "long_pivot_scalar_u_output_per_producer_target_u",
    ]
    print(",".join(columns))
    for row in rows:
        print(",".join(format_value(row[column]) for column in columns))


def format_value(value: object) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("traces", type=pathlib.Path, nargs="+")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    rows = [summarize(path) for path in args.traces]
    if args.json:
        print(json.dumps(rows, indent=2, sort_keys=True))
    else:
        print_table(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
