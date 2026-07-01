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
    long_panel_rows = 0

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
            long_panel_rows += parse_int(row.get("panel_update_rows"))

    completed, total = split_completed(last_trace.get("completed"))
    scalar = parse_int(last_trace.get("scalar_u_entries"))
    output = parse_int(last_trace.get("scalar_u_output"))
    producer_target_u = parse_int(last_trace.get("producer_target_u_entries"))
    panel_rows = parse_int(last_trace.get("panel_update_rows"))
    summary: dict[str, object] = {
        "path": str(path),
        "trace_events": trace_event_count,
        "last_event": last_trace.get("event", ""),
        "completed": completed,
        "total": total,
        "scalar_u_entries": scalar,
        "scalar_u_output": output,
        "producer_target_u_entries": producer_target_u,
        "panel_update_rows": panel_rows,
        "long_rows": long_rows,
        "long_max_row": long_max_row,
        "long_pivot_rows": long_pivot_rows,
        "long_scalar_u_entries": long_scalar,
        "long_scalar_u_internal": long_internal,
        "long_scalar_u_output": long_output,
        "long_producer_targets": long_producer_targets,
        "long_producer_target_u_entries": long_producer_target_u,
        "long_panel_update_rows": long_panel_rows,
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
        "producer_underfilled",
        "producer_low_saved_stream",
        "producer_reject_bad_state",
        "producer_reject_epoch",
        "producer_reject_dep_absent",
        "producer_reject_not_root",
        "panel_update_rows",
        "long_rows",
        "long_max_row",
        "long_pivot_rows",
        "long_scalar_u_entries",
        "long_scalar_u_output",
        "long_producer_target_u_entries",
    ]
    print(",".join(columns))
    for row in rows:
        print(",".join(str(row[column]) for column in columns))


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
