#!/usr/bin/env python3
"""Build SuiteSparse-group-disjoint development and holdout manifests.

The paper suite is useful for comparison with prior work, but it is a poor
place to tune and validate the same policy.  This script excludes every
SuiteSparse group represented by that suite, partitions the remaining groups
before selecting matrices, and chooses at most one matrix from each group.
"""

from __future__ import annotations

import argparse
import collections
import csv
from dataclasses import dataclass
import hashlib
import io
import pathlib
import sys
import urllib.request


SSTATS_URL = "https://sparse.tamu.edu/files/ssstats.csv"
DEFAULT_SEED = "kls-suitesparse-generalization-v1"


@dataclass(frozen=True)
class MatrixInfo:
    group: str
    name: str
    nrows: int
    ncols: int
    nnz: int
    is_real: bool
    pattern_symmetry: float
    kind: str

    @property
    def canonical_name(self) -> str:
        return f"{self.group}/{self.name}"


def read_text(path: pathlib.Path | None) -> str:
    if path is not None:
        return path.read_text(encoding="utf-8")
    with urllib.request.urlopen(SSTATS_URL, timeout=60) as response:
        return response.read().decode("utf-8", errors="replace")


def load_index(path: pathlib.Path | None) -> tuple[str, str, list[MatrixInfo]]:
    text = read_text(path)
    lines = text.splitlines()
    if len(lines) < 3:
        raise ValueError("SuiteSparse index has no matrix rows")
    timestamp = lines[1].strip()
    matrices: list[MatrixInfo] = []
    for row in csv.reader(io.StringIO("\n".join(lines[2:]))):
        if len(row) < 12:
            continue
        try:
            nrows = int(float(row[2]))
            ncols = int(float(row[3]))
            nnz = int(float(row[4]))
            is_real = int(float(row[5])) == 1
            pattern_symmetry = float(row[9])
        except ValueError:
            continue
        group = row[0].strip()
        name = row[1].strip()
        if not group or not name:
            continue
        matrices.append(
            MatrixInfo(
                group=group,
                name=name,
                nrows=nrows,
                ncols=ncols,
                nnz=nnz,
                is_real=is_real,
                pattern_symmetry=pattern_symmetry,
                kind=row[11].strip(),
            )
        )
    if not matrices:
        raise ValueError("SuiteSparse index contained no usable matrix rows")
    digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
    return timestamp, digest, matrices


def read_manifest_names(path: pathlib.Path) -> list[str]:
    names: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        name = raw.split("#", 1)[0].strip()
        if name:
            names.append(name)
    return names


def resolve_paper_groups(matrices: list[MatrixInfo], names: list[str]) -> set[str]:
    by_name: dict[str, list[MatrixInfo]] = collections.defaultdict(list)
    by_canonical: dict[str, MatrixInfo] = {}
    for matrix in matrices:
        by_name[matrix.name.casefold()].append(matrix)
        by_canonical[matrix.canonical_name.casefold()] = matrix

    missing: list[str] = []
    ambiguous: list[str] = []
    resolved: list[MatrixInfo] = []
    for requested in names:
        canonical = by_canonical.get(requested.casefold())
        matches = (
            [canonical]
            if canonical is not None
            else by_name.get(requested.casefold(), [])
        )
        if not matches:
            missing.append(requested)
        elif len(matches) != 1:
            ambiguous.append(requested)
        else:
            resolved.append(matches[0])
    if missing or ambiguous:
        details: list[str] = []
        if missing:
            details.append("missing=" + ", ".join(sorted(missing)))
        if ambiguous:
            details.append("ambiguous=" + ", ".join(sorted(ambiguous)))
        raise ValueError(
            "paper manifest did not resolve uniquely: " + "; ".join(details)
        )
    return {matrix.group.casefold() for matrix in resolved}


def stable_score(seed: str, purpose: str, value: str) -> bytes:
    payload = f"{seed}\0{purpose}\0{value.casefold()}".encode("utf-8")
    return hashlib.sha256(payload).digest()


def stratum(matrix: MatrixInfo) -> tuple[str, str, str]:
    if matrix.nrows < 10_000:
        size = "small"
    elif matrix.nrows < 100_000:
        size = "medium"
    else:
        size = "large"
    degree = matrix.nnz / max(matrix.nrows, 1)
    if degree < 5.0:
        density = "sparse"
    elif degree < 20.0:
        density = "moderate"
    else:
        density = "dense"
    symmetry = "symmetric" if matrix.pattern_symmetry >= 0.9 else "unsymmetric"
    return size, density, symmetry


def representative_by_group(matrices: list[MatrixInfo], seed: str) -> list[MatrixInfo]:
    grouped: dict[str, list[MatrixInfo]] = collections.defaultdict(list)
    for matrix in matrices:
        grouped[matrix.group.casefold()].append(matrix)
    representatives: list[MatrixInfo] = []
    for candidates in grouped.values():
        representatives.append(
            min(
                candidates,
                key=lambda matrix: stable_score(seed, "matrix", matrix.canonical_name),
            )
        )
    return representatives


def balanced_sample(
    matrices: list[MatrixInfo], count: int, seed: str, label: str
) -> list[MatrixInfo]:
    buckets: dict[tuple[str, str, str], list[MatrixInfo]] = collections.defaultdict(
        list
    )
    for matrix in matrices:
        buckets[stratum(matrix)].append(matrix)
    for bucket in buckets.values():
        bucket.sort(
            key=lambda matrix: stable_score(
                seed, f"{label}-sample", matrix.canonical_name
            ),
            reverse=True,
        )

    selected: list[MatrixInfo] = []
    while len(selected) < count:
        progressed = False
        for key in sorted(buckets):
            if buckets[key] and len(selected) < count:
                selected.append(buckets[key].pop())
                progressed = True
        if not progressed:
            raise ValueError(
                f"only {len(selected)} eligible {label} groups for requested {count}"
            )
    return selected


def sanitize_comment(text: str) -> str:
    return " ".join(text.replace("#", " ").split())


def render_manifest(
    label: str,
    selected: list[MatrixInfo],
    timestamp: str,
    index_digest: str,
    paper_groups: set[str],
    args: argparse.Namespace,
) -> str:
    excluded = ",".join(sorted(paper_groups))
    lines = [
        "# Generated by scripts/build_generalization_manifests.py; do not hand-pick rows.",
        f"# Set: {label}",
        f"# SuiteSparse index timestamp: {timestamp}",
        f"# SuiteSparse index SHA256: {index_digest}",
        f"# Seed: {args.seed}",
        (
            "# Filters: real square matrices; "
            f"rows={args.min_rows}..{args.max_rows}; "
            f"nnz={args.min_nnz}..{args.max_nnz}; one matrix per group."
        ),
        "# All groups in the paper benchmark are excluded before partitioning.",
        f"# Excluded paper groups: {excluded}",
        "# Entry comments preserve the SuiteSparse family and selection stratum.",
        "",
    ]
    for matrix in selected:
        size, density, symmetry = stratum(matrix)
        lines.append(
            f"{matrix.name}  # group={matrix.group} rows={matrix.nrows} "
            f"nnz={matrix.nnz} stratum={size}/{density}/{symmetry} "
            f"kind={sanitize_comment(matrix.kind)}"
        )
    return "\n".join(lines) + "\n"


def print_summary(label: str, selected: list[MatrixInfo], eligible_groups: int) -> None:
    strata = collections.Counter(stratum(matrix) for matrix in selected)
    print(f"{label}: selected={len(selected)} eligible_groups={eligible_groups}")
    for key, count in sorted(strata.items()):
        print(f"  {'/'.join(key)}: {count}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--index",
        type=pathlib.Path,
        help=f"local sstats.csv snapshot; default downloads {SSTATS_URL}",
    )
    parser.add_argument(
        "--paper-manifest",
        type=pathlib.Path,
        default=pathlib.Path("bench/suitesparse_paper_manifest.txt"),
    )
    parser.add_argument(
        "--development-out",
        type=pathlib.Path,
        default=pathlib.Path("bench/suitesparse_generalization_dev_manifest.txt"),
    )
    parser.add_argument(
        "--holdout-out",
        type=pathlib.Path,
        default=pathlib.Path("bench/suitesparse_generalization_holdout_manifest.txt"),
    )
    parser.add_argument("--development-count", type=int, default=24)
    parser.add_argument("--holdout-count", type=int, default=24)
    parser.add_argument("--min-rows", type=int, default=1_000)
    parser.add_argument("--max-rows", type=int, default=300_000)
    parser.add_argument("--min-nnz", type=int, default=3_000)
    parser.add_argument("--max-nnz", type=int, default=1_000_000)
    parser.add_argument("--seed", default=DEFAULT_SEED)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail instead of writing if committed manifests are stale",
    )
    args = parser.parse_args()

    if args.dry_run and args.check:
        parser.error("--dry-run and --check are mutually exclusive")
    if args.development_count <= 0 or args.holdout_count <= 0:
        parser.error("manifest counts must be positive")
    if not (0 < args.min_rows <= args.max_rows):
        parser.error("row bounds are invalid")
    if not (0 < args.min_nnz <= args.max_nnz):
        parser.error("nnz bounds are invalid")

    try:
        timestamp, index_digest, matrices = load_index(args.index)
        paper_groups = resolve_paper_groups(
            matrices, read_manifest_names(args.paper_manifest)
        )
    except (OSError, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        return 1

    name_counts = collections.Counter(matrix.name.casefold() for matrix in matrices)
    candidates = [
        matrix
        for matrix in matrices
        if matrix.is_real
        and matrix.nrows == matrix.ncols
        and args.min_rows <= matrix.nrows <= args.max_rows
        and args.min_nnz <= matrix.nnz <= args.max_nnz
        and matrix.group.casefold() not in paper_groups
        # Existing suite runners intentionally key files by basename.  Keep
        # manifests unambiguous until those interfaces accept group/name.
        and name_counts[matrix.name.casefold()] == 1
    ]
    representatives = representative_by_group(candidates, args.seed)
    development_pool = [
        matrix
        for matrix in representatives
        if (stable_score(args.seed, "partition", matrix.group)[0] & 1) == 0
    ]
    holdout_pool = [
        matrix
        for matrix in representatives
        if (stable_score(args.seed, "partition", matrix.group)[0] & 1) == 1
    ]
    try:
        development = balanced_sample(
            development_pool, args.development_count, args.seed, "development"
        )
        holdout = balanced_sample(
            holdout_pool, args.holdout_count, args.seed, "holdout"
        )
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    development_groups = {matrix.group.casefold() for matrix in development}
    holdout_groups = {matrix.group.casefold() for matrix in holdout}
    if (
        development_groups & holdout_groups
        or (development_groups | holdout_groups) & paper_groups
    ):
        print("internal error: generalization groups are not disjoint", file=sys.stderr)
        return 1

    development_text = render_manifest(
        "development", development, timestamp, index_digest, paper_groups, args
    )
    holdout_text = render_manifest(
        "holdout", holdout, timestamp, index_digest, paper_groups, args
    )
    if args.dry_run:
        print(development_text, end="")
        print(holdout_text, end="")
    elif args.check:
        stale: list[pathlib.Path] = []
        for path, expected in (
            (args.development_out, development_text),
            (args.holdout_out, holdout_text),
        ):
            try:
                actual = path.read_text(encoding="utf-8")
            except OSError:
                actual = ""
            if actual != expected:
                stale.append(path)
        if stale:
            print(
                "stale generalization manifest(s): "
                + ", ".join(str(path) for path in stale),
                file=sys.stderr,
            )
            return 1
    else:
        args.development_out.parent.mkdir(parents=True, exist_ok=True)
        args.holdout_out.parent.mkdir(parents=True, exist_ok=True)
        args.development_out.write_text(development_text, encoding="utf-8")
        args.holdout_out.write_text(holdout_text, encoding="utf-8")

    print_summary("development", development, len(development_pool))
    print_summary("holdout", holdout, len(holdout_pool))
    print(f"paper_groups_excluded={len(paper_groups)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
