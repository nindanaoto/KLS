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
    if as_int(row, "fast_repaired_tail_restart_ready", 0) == 1:
        return "strict_ready"
    if as_int(row, "fast_rejected_tail_repair_ready", 0) != 1:
        return "no_tolerance_valid_tail_candidate"
    if as_int(row, "fast_repaired_pivot_matches_tail_candidate", 0) != 1:
        return "repair_picked_different_row"
    if as_int(row, "fast_repaired_prefix_changed_pivots", 0) > 0:
        return "prefix_pivots_changed"
    if (
        as_int(row, "fast_rejected_etree_columns", 0) <= 0
        or as_float(row, "fast_rejected_etree_work") <= 0.0
    ):
        return "missing_etree_tail"
    if as_int(row, "fast_repaired_suffix_changed_pivots", 0) <= 0:
        return "no_suffix_pivot_change"
    return "other_not_ready"


def compact_record(row: dict[str, object], reason: str) -> dict[str, object]:
    return {
        "matrix": matrix_name(row),
        "reason": reason,
        "fast_rejected_pivot": as_int(row, "fast_rejected_pivot", -1),
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
        "fast_repaired_tail_restart_columns": as_int(
            row, "fast_repaired_tail_restart_columns", 0
        ),
        "fast_repaired_tail_restart_work": as_float(
            row, "fast_repaired_tail_restart_work"
        ),
        "fast_repaired_tail_restart_saved_work": as_float(
            row, "fast_repaired_tail_restart_saved_work"
        ),
        "fast_rejected_etree_columns": as_int(row, "fast_rejected_etree_columns", 0),
        "fast_rejected_etree_work": as_float(row, "fast_rejected_etree_work"),
    }


def print_records(title: str, records: list[dict[str, object]], limit: int) -> None:
    if not records or limit <= 0:
        return
    print(f"\n{title}:")
    for record in records[:limit]:
        print(
            "  {matrix}: reason={reason} block_work={block:.6g} "
            "tail_work={tail:.6g} saved_work={saved:.6g} "
            "prefix_changes={prefix} suffix_changes={suffix}".format(
                matrix=record["matrix"],
                reason=record["reason"],
                block=float(record["fast_repaired_block_work"]),
                tail=float(record["fast_repaired_tail_restart_work"]),
                saved=float(record["fast_repaired_tail_restart_saved_work"]),
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
    rejected = [record for record in records if record["reason"] != "no_reject"]
    repaired = [
        record for record in rejected if float(record["fast_repaired_block_work"]) > 0.0
    ]
    blocked = [record for record in repaired if record["reason"] != "strict_ready"]

    ready_block_work = sum(float(r["fast_repaired_block_work"]) for r in ready)
    ready_tail_work = sum(float(r["fast_repaired_tail_restart_work"]) for r in ready)
    ready_saved_work = sum(
        float(r["fast_repaired_tail_restart_saved_work"]) for r in ready
    )
    repaired_block_work = sum(float(r["fast_repaired_block_work"]) for r in repaired)
    blocked_block_work = sum(float(r["fast_repaired_block_work"]) for r in blocked)

    summary = {
        "rows_total": len(rows),
        "rows_with_tail_restart_fields": rows_with_new_fields,
        "rejected_rows": len(rejected),
        "repaired_rows_with_block_work": len(repaired),
        "strict_ready_rows": len(ready),
        "reason_counts": dict(sorted(reason_counts.items())),
        "strict_ready_block_work_total": ready_block_work,
        "strict_ready_tail_work_total": ready_tail_work,
        "strict_ready_saved_work_total": ready_saved_work,
        "strict_ready_saved_work_fraction": (
            ready_saved_work / ready_block_work if ready_block_work > 0.0 else 0.0
        ),
        "blocked_repaired_block_work_total": blocked_block_work,
        "repaired_block_work_total": repaired_block_work,
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
    print_records("Largest strict-ready saved work", ready_by_saved, args.max_rows)
    print_records("Largest repaired blockers", blocked_by_work, args.max_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
