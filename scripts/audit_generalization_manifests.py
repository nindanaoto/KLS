#!/usr/bin/env python3
"""Validate the committed SuiteSparse generalization split without network I/O."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import pathlib
import re
import sys


ENTRY_RE = re.compile(
    r"^(?P<name>\S+)\s+#\s+group=(?P<group>\S+)\s+"
    r"rows=(?P<rows>\d+)\s+nnz=(?P<nnz>\d+)\s+"
    r"stratum=(?P<stratum>\S+)\s+kind="
)


@dataclass(frozen=True)
class Manifest:
    label: str
    index_digest: str
    seed: str
    excluded_groups: frozenset[str]
    names: tuple[str, ...]
    groups: tuple[str, ...]


def header_value(headers: dict[str, str], key: str, path: pathlib.Path) -> str:
    value = headers.get(key, "")
    if not value:
        raise ValueError(f"{path}: missing '# {key}:' header")
    return value


def load_manifest(path: pathlib.Path) -> Manifest:
    headers: dict[str, str] = {}
    names: list[str] = []
    groups: list[str] = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            header = line[1:].strip()
            key, separator, value = header.partition(":")
            if separator:
                headers[key.strip()] = value.strip()
            continue
        match = ENTRY_RE.match(line)
        if match is None:
            raise ValueError(f"{path}:{line_number}: malformed manifest entry")
        if int(match.group("rows")) <= 0 or int(match.group("nnz")) <= 0:
            raise ValueError(f"{path}:{line_number}: non-positive matrix size")
        if len(match.group("stratum").split("/")) != 3:
            raise ValueError(f"{path}:{line_number}: malformed stratum")
        names.append(match.group("name").casefold())
        groups.append(match.group("group").casefold())

    excluded = frozenset(
        group.strip().casefold()
        for group in header_value(headers, "Excluded paper groups", path).split(",")
        if group.strip()
    )
    return Manifest(
        label=header_value(headers, "Set", path),
        index_digest=header_value(headers, "SuiteSparse index SHA256", path),
        seed=header_value(headers, "Seed", path),
        excluded_groups=excluded,
        names=tuple(names),
        groups=tuple(groups),
    )


def duplicates(values: tuple[str, ...]) -> list[str]:
    seen: set[str] = set()
    repeated: set[str] = set()
    for value in values:
        if value in seen:
            repeated.add(value)
        seen.add(value)
    return sorted(repeated)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--development", type=pathlib.Path, required=True)
    parser.add_argument("--holdout", type=pathlib.Path, required=True)
    parser.add_argument("--expected-count", type=int, default=24)
    args = parser.parse_args()

    try:
        development = load_manifest(args.development)
        holdout = load_manifest(args.holdout)
        if development.label != "development" or holdout.label != "holdout":
            raise ValueError("manifest Set headers do not match their roles")
        for manifest in (development, holdout):
            if len(manifest.names) != args.expected_count:
                raise ValueError(
                    f"{manifest.label}: expected {args.expected_count} entries, "
                    f"found {len(manifest.names)}"
                )
            repeated_names = duplicates(manifest.names)
            repeated_groups = duplicates(manifest.groups)
            if repeated_names or repeated_groups:
                raise ValueError(
                    f"{manifest.label}: duplicate names={repeated_names} "
                    f"groups={repeated_groups}"
                )
            leaked = set(manifest.groups) & manifest.excluded_groups
            if leaked:
                raise ValueError(
                    f"{manifest.label}: paper groups leaked into split: "
                    + ", ".join(sorted(leaked))
                )
        if development.index_digest != holdout.index_digest:
            raise ValueError("development and holdout use different index digests")
        if development.seed != holdout.seed:
            raise ValueError("development and holdout use different seeds")
        if development.excluded_groups != holdout.excluded_groups:
            raise ValueError("development and holdout exclude different paper groups")
        overlap = set(development.groups) & set(holdout.groups)
        if overlap:
            raise ValueError(
                "development/holdout group overlap: " + ", ".join(sorted(overlap))
            )
        name_overlap = set(development.names) & set(holdout.names)
        if name_overlap:
            raise ValueError(
                "development/holdout name overlap: " + ", ".join(sorted(name_overlap))
            )
    except (OSError, ValueError) as exc:
        print(exc, file=sys.stderr)
        return 1

    print(
        f"generalization manifests: development={len(development.names)} "
        f"holdout={len(holdout.names)} paper_groups="
        f"{len(development.excluded_groups)} group_overlap=0"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
