#!/usr/bin/env python3
"""Summarize KLS fast-reject tail-restart opportunities from JSONL rows."""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys
from collections import Counter
from collections.abc import Iterable


def matrix_name(row: dict[str, object]) -> str:
    return pathlib.Path(str(row.get("matrix", "<unknown>"))).name


def as_float(row: dict[str, object], key: str, default: float = 0.0) -> float:
    value = row.get(key, default)
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    return number if math.isfinite(number) else default


def as_int(row: dict[str, object], key: str, default: int = 0) -> int:
    value = row.get(key, default)
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def read_jsonl(path_arg: str) -> Iterable[dict[str, object]]:
    if path_arg == "-":
        for line_no, line in enumerate(sys.stdin, 1):
            if not line.strip():
                continue
            row = json.loads(line)
            row["_source"] = "<stdin>"
            row["_line_no"] = line_no
            yield row
        return

    path = pathlib.Path(path_arg)
    with path.open("r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, 1):
            if not line.strip():
                continue
            row = json.loads(line)
            row["_source"] = str(path)
            row["_line_no"] = line_no
            yield row


def classify(row: dict[str, object]) -> str:
    if as_int(row, "fast_rejected_pivot", -1) < 0:
        return "no_reject"
    if as_int(row, "fast_block_restarts", 0) <= 0:
        return "no_block_repair"
    if as_int(row, "fast_tail_restarts", 0) > 0:
        return "tail_restart_executed"
    if as_int(row, "fast_repaired_tail_restart_ready", 0) == 1:
        return "strict_ready"
    refresh_state = as_int(row, "fast_rejected_refresh_state", 0)
    if refresh_state == 0:
        return "unknown_refresh_state"
    if refresh_state == 2:
        return "all_current_reject_state"
    if refresh_state != 1:
        return "non_prefix_reject_state"
    if (
        as_int(row, "fast_rejected_tail_candidate_row", -1) < 0
        or as_int(row, "fast_rejected_tail_candidate_count", 0) <= 0
    ):
        return "no_row_tail_candidate_diagnostic"
    if as_int(row, "fast_rejected_tail_repair_ready", 0) != 1:
        return "tail_candidate_not_repair_ready"
    if as_int(row, "fast_repaired_pivot_matches_tail_candidate", 0) != 1:
        return "repair_picked_different_row"
    if as_int(row, "fast_repaired_prefix_changed_pivots", 0) > 0:
        return "prefix_pivots_changed"
    if as_int(row, "fast_rejected_pivot", -1) == as_int(
        row, "fast_rejected_block_start", -2
    ):
        return "root_reject_no_serial_prefix"
    tail_columns = as_int(
        row,
        "fast_rejected_pivoting_tail_columns",
        as_int(row, "fast_rejected_etree_columns", 0),
    )
    tail_work = as_float(
        row,
        "fast_rejected_pivoting_tail_work",
        as_float(row, "fast_rejected_etree_work"),
    )
    if tail_columns <= 0 or tail_work <= 0.0:
        return "missing_pivoting_tail"
    if (
        "fast_rejected_pivoting_tail_contains_reject" in row
        and as_int(row, "fast_rejected_pivoting_tail_contains_reject", 0) != 1
    ):
        return "pivoting_tail_missing_reject"
    if (
        "fast_rejected_pivoting_tail_topological" in row
        and as_int(row, "fast_rejected_pivoting_tail_topological", 0) != 1
    ):
        return "pivoting_tail_not_topological"
    if as_int(row, "fast_repaired_suffix_changed_pivots", 0) <= 0:
        return "no_suffix_pivot_change"
    return "other_not_ready"


def compact_record(row: dict[str, object], reason: str) -> dict[str, object]:
    suffix_tail_work = as_float(row, "fast_repaired_tail_restart_work")
    pivoting_tail_work = as_float(
        row,
        "fast_rejected_pivoting_tail_work",
        as_float(row, "fast_rejected_etree_work"),
    )
    suffix_tail_columns = as_int(row, "fast_repaired_tail_restart_columns", 0)
    pivoting_tail_columns = as_int(
        row,
        "fast_rejected_pivoting_tail_columns",
        as_int(row, "fast_rejected_etree_columns", 0),
    )
    return {
        "matrix": matrix_name(row),
        "reason": reason,
        "fast_rejected_pivot": as_int(row, "fast_rejected_pivot", -1),
        "fast_rejected_block_start": as_int(
            row, "fast_rejected_block_start", -1
        ),
        "fast_rejected_tail_candidate_row": as_int(
            row, "fast_rejected_tail_candidate_row", -1
        ),
        "fast_repaired_pivot_row": as_int(row, "fast_repaired_pivot_row", -1),
        "fast_repaired_prefix_changed_pivots": as_int(
            row, "fast_repaired_prefix_changed_pivots", 0
        ),
        "fast_repaired_suffix_changed_pivots": as_int(
            row, "fast_repaired_suffix_changed_pivots", 0
        ),
        "fast_repaired_block_work": as_float(row, "fast_repaired_block_work"),
        "fast_repaired_tail_restart_columns": suffix_tail_columns,
        "fast_repaired_tail_restart_work": suffix_tail_work,
        "fast_repaired_tail_restart_saved_work": as_float(
            row, "fast_repaired_tail_restart_saved_work"
        ),
        "fast_tail_restarts": as_int(row, "fast_tail_restarts", 0),
        "fast_rejected_etree_columns": as_int(row, "fast_rejected_etree_columns", 0),
        "fast_rejected_etree_work": as_float(row, "fast_rejected_etree_work"),
        "fast_rejected_pivoting_tail_columns": pivoting_tail_columns,
        "fast_rejected_pivoting_tail_work": pivoting_tail_work,
        "suffix_tail_overcompute_columns": max(
            0, suffix_tail_columns - pivoting_tail_columns
        ),
        "suffix_tail_overcompute_work": max(
            0.0, suffix_tail_work - pivoting_tail_work
        ),
        "fast_rejected_pivoting_tail_first": as_int(
            row, "fast_rejected_pivoting_tail_first", -1
        ),
        "fast_rejected_pivoting_tail_last": as_int(
            row, "fast_rejected_pivoting_tail_last", -1
        ),
        "fast_rejected_pivoting_tail_contains_reject": as_int(
            row, "fast_rejected_pivoting_tail_contains_reject", 0
        ),
        "fast_rejected_pivoting_tail_topological": as_int(
            row, "fast_rejected_pivoting_tail_topological", 0
        ),
        "fast_rejected_row_tail_columns": as_int(
            row, "fast_rejected_row_tail_columns", 0
        ),
        "fast_rejected_row_tail_work": as_float(
            row, "fast_rejected_row_tail_work"
        ),
    }


def print_records(title: str, records: list[dict[str, object]], limit: int) -> None:
    if not records or limit <= 0:
        return
    print(f"\n{title}:")
    for record in records[:limit]:
        print(
            "  {matrix}: reason={reason} block_work={block:.6g} "
            "tail_cols={tail_cols} tail_work={tail:.6g} "
            "row_tail_work={row_tail:.6g} "
            "pivoting_tail_cols={pivoting_tail_cols} "
            "pivoting_tail_work={pivoting_tail:.6g} "
            "overcompute_work={overcompute:.6g} "
            "saved_work={saved:.6g} "
            "block_start={block_start} "
            "tail_first={tail_first} tail_last={tail_last} "
            "tail_topo={tail_topo} "
            "tail_row={tail_row} repair_row={repair_row} "
            "prefix_changes={prefix} suffix_changes={suffix}".format(
                matrix=record["matrix"],
                reason=record["reason"],
                block=float(record["fast_repaired_block_work"]),
                tail_cols=record["fast_repaired_tail_restart_columns"],
                tail=float(record["fast_repaired_tail_restart_work"]),
                row_tail=float(record["fast_rejected_row_tail_work"]),
                pivoting_tail_cols=record[
                    "fast_rejected_pivoting_tail_columns"
                ],
                pivoting_tail=float(
                    record["fast_rejected_pivoting_tail_work"]
                ),
                overcompute=float(record["suffix_tail_overcompute_work"]),
                saved=float(record["fast_repaired_tail_restart_saved_work"]),
                block_start=record["fast_rejected_block_start"],
                tail_first=record["fast_rejected_pivoting_tail_first"],
                tail_last=record["fast_rejected_pivoting_tail_last"],
                tail_topo=record[
                    "fast_rejected_pivoting_tail_topological"
                ],
                tail_row=record["fast_rejected_tail_candidate_row"],
                repair_row=record["fast_repaired_pivot_row"],
                prefix=record["fast_repaired_prefix_changed_pivots"],
                suffix=record["fast_repaired_suffix_changed_pivots"],
            )
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--jsonl",
        nargs="+",
        required=True,
        help="One or more KLS benchmark JSONL files, or '-' for stdin.",
    )
    parser.add_argument("--max-rows", type=int, default=12)
    parser.add_argument(
        "--json-only",
        action="store_true",
        help="Only print the machine-readable summary.",
    )
    args = parser.parse_args()

    rows: list[dict[str, object]] = []
    for path_arg in args.jsonl:
        rows.extend(read_jsonl(path_arg))
    if not rows:
        raise SystemExit("no JSON rows found")

    reason_counts: Counter[str] = Counter()
    records: list[dict[str, object]] = []
    rows_with_new_fields = 0
    for row in rows:
        if "fast_repaired_tail_restart_ready" in row:
            rows_with_new_fields += 1
        reason = classify(row)
        reason_counts[reason] += 1
        records.append(compact_record(row, reason))

    ready = [record for record in records if record["reason"] == "strict_ready"]
    executed = [
        record for record in records if record["reason"] == "tail_restart_executed"
    ]
    rejected = [record for record in records if record["reason"] != "no_reject"]
    repaired = [
        record for record in rejected if float(record["fast_repaired_block_work"]) > 0.0
    ]
    blocked = [
        record for record in repaired
        if record["reason"] not in {"strict_ready", "tail_restart_executed"}
    ]

    ready_block_work = sum(float(r["fast_repaired_block_work"]) for r in ready)
    ready_tail_work = sum(float(r["fast_repaired_tail_restart_work"]) for r in ready)
    ready_saved_work = sum(
        float(r["fast_repaired_tail_restart_saved_work"]) for r in ready
    )
    executed_block_work = sum(float(r["fast_repaired_block_work"]) for r in executed)
    executed_tail_work = sum(
        float(r["fast_repaired_tail_restart_work"]) for r in executed
    )
    executed_saved_work = sum(
        float(r["fast_repaired_tail_restart_saved_work"]) for r in executed
    )
    executed_pivoting_tail_work = sum(
        float(r["fast_rejected_pivoting_tail_work"]) for r in executed
    )
    executed_overcompute_work = sum(
        float(r["suffix_tail_overcompute_work"]) for r in executed
    )
    executed_overcompute_columns = sum(
        int(r["suffix_tail_overcompute_columns"]) for r in executed
    )
    repaired_block_work = sum(float(r["fast_repaired_block_work"]) for r in repaired)
    blocked_block_work = sum(float(r["fast_repaired_block_work"]) for r in blocked)
    row_tail = [
        record for record in rejected
        if int(record["fast_rejected_row_tail_columns"]) > 0
    ]
    row_tail_work = sum(float(r["fast_rejected_row_tail_work"]) for r in row_tail)
    pivoting_tail = [
        record for record in rejected
        if int(record["fast_rejected_pivoting_tail_columns"]) > 0
    ]
    pivoting_tail_work = sum(
        float(r["fast_rejected_pivoting_tail_work"]) for r in pivoting_tail
    )

    summary = {
        "rows_total": len(rows),
        "rows_with_tail_restart_fields": rows_with_new_fields,
        "rejected_rows": len(rejected),
        "repaired_rows_with_block_work": len(repaired),
        "strict_ready_rows": len(ready),
        "tail_restart_executed_rows": len(executed),
        "reason_counts": dict(sorted(reason_counts.items())),
        "strict_ready_block_work_total": ready_block_work,
        "strict_ready_tail_work_total": ready_tail_work,
        "strict_ready_saved_work_total": ready_saved_work,
        "strict_ready_saved_work_fraction": (
            ready_saved_work / ready_block_work if ready_block_work > 0.0 else 0.0
        ),
        "tail_restart_executed_block_work_total": executed_block_work,
        "tail_restart_executed_tail_work_total": executed_tail_work,
        "tail_restart_executed_pivoting_tail_work_total": (
            executed_pivoting_tail_work
        ),
        "tail_restart_executed_overcompute_work_total": (
            executed_overcompute_work
        ),
        "tail_restart_executed_overcompute_columns_total": (
            executed_overcompute_columns
        ),
        "tail_restart_executed_overcompute_work_fraction": (
            executed_overcompute_work / executed_tail_work
            if executed_tail_work > 0.0 else 0.0
        ),
        "tail_restart_executed_saved_work_total": executed_saved_work,
        "tail_restart_executed_saved_work_fraction": (
            executed_saved_work / executed_block_work
            if executed_block_work > 0.0 else 0.0
        ),
        "blocked_repaired_block_work_total": blocked_block_work,
        "repaired_block_work_total": repaired_block_work,
        "rows_with_row_tail_scope": len(row_tail),
        "row_tail_work_total": row_tail_work,
        "rows_with_pivoting_tail_scope": len(pivoting_tail),
        "pivoting_tail_work_total": pivoting_tail_work,
    }
    print(json.dumps(summary, indent=2, sort_keys=True))

    if args.json_only:
        return 0

    ready_by_saved = sorted(
        ready,
        key=lambda record: float(record["fast_repaired_tail_restart_saved_work"]),
        reverse=True,
    )
    blocked_by_work = sorted(
        blocked,
        key=lambda record: float(record["fast_repaired_block_work"]),
        reverse=True,
    )
    executed_by_saved = sorted(
        executed,
        key=lambda record: float(record["fast_repaired_tail_restart_saved_work"]),
        reverse=True,
    )
    print_records("Largest strict-ready saved work", ready_by_saved, args.max_rows)
    print_records("Largest executed tail restarts", executed_by_saved, args.max_rows)
    print_records("Largest repaired blockers", blocked_by_work, args.max_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
