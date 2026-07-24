#!/usr/bin/env python3
"""Create a deterministic simultaneous row/column permutation holdout."""

from __future__ import annotations

import argparse
import pathlib
import random


def read_coordinate_matrix(
    path: pathlib.Path,
) -> tuple[list[str], str, str, int, list[tuple[int, int, tuple[str, ...]]]]:
    with path.open("r", encoding="utf-8") as source:
        header = source.readline().strip()
        fields = header.split()
        if len(fields) != 5 or fields[:3] != [
            "%%MatrixMarket",
            "matrix",
            "coordinate",
        ]:
            raise ValueError("only MatrixMarket coordinate matrices are supported")
        value_kind = fields[3].lower()
        symmetry = fields[4].lower()
        if symmetry not in {"general", "symmetric"}:
            raise ValueError(f"unsupported MatrixMarket symmetry: {symmetry}")

        comments: list[str] = []
        for line in source:
            if line.startswith("%"):
                comments.append(line.rstrip("\n"))
                continue
            dimensions = line.split()
            break
        else:
            raise ValueError("missing MatrixMarket dimensions")
        if len(dimensions) != 3:
            raise ValueError("invalid MatrixMarket dimensions")
        rows, columns, declared_nnz = map(int, dimensions)
        if rows != columns:
            raise ValueError("simultaneous permutation requires a square matrix")

        entries: list[tuple[int, int, tuple[str, ...]]] = []
        for line in source:
            if not line.strip() or line.startswith("%"):
                continue
            parts = line.split()
            required = 2 if value_kind == "pattern" else 3
            if len(parts) < required:
                raise ValueError("invalid MatrixMarket coordinate entry")
            row = int(parts[0]) - 1
            column = int(parts[1]) - 1
            if row < 0 or row >= rows or column < 0 or column >= columns:
                raise ValueError("MatrixMarket coordinate is out of bounds")
            entries.append((row, column, tuple(parts[2:])))
        if len(entries) != declared_nnz:
            raise ValueError(
                f"declared {declared_nnz} entries but read {len(entries)}"
            )
    return comments, value_kind, symmetry, rows, entries


def adjacent_permutation(order: int, swaps: int, seed: int) -> list[int]:
    if swaps < 0:
        raise ValueError("adjacent swap count must be nonnegative")
    permutation = list(range(order))
    if order < 2 and swaps:
        raise ValueError("cannot swap a matrix with fewer than two rows")
    rng = random.Random(seed)
    for _ in range(swaps):
        left = rng.randrange(order - 1)
        permutation[left], permutation[left + 1] = (
            permutation[left + 1],
            permutation[left],
        )
    old_to_new = [0] * order
    for new, old in enumerate(permutation):
        old_to_new[old] = new
    return old_to_new


def permute_entries(
    entries: list[tuple[int, int, tuple[str, ...]]],
    old_to_new: list[int],
    symmetry: str,
) -> list[tuple[int, int, tuple[str, ...]]]:
    transformed: list[tuple[int, int, tuple[str, ...]]] = []
    for row, column, value in entries:
        new_row = old_to_new[row]
        new_column = old_to_new[column]
        if symmetry == "symmetric" and new_row < new_column:
            new_row, new_column = new_column, new_row
        transformed.append((new_row, new_column, value))
    transformed.sort(key=lambda entry: (entry[1], entry[0]))
    return transformed


def write_coordinate_matrix(
    path: pathlib.Path,
    comments: list[str],
    value_kind: str,
    symmetry: str,
    order: int,
    entries: list[tuple[int, int, tuple[str, ...]]],
    swaps: int,
    seed: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as output:
        output.write(
            f"%%MatrixMarket matrix coordinate {value_kind} {symmetry}\n"
        )
        output.write(
            f"% metamorphic simultaneous permutation: adjacent_swaps={swaps} seed={seed}\n"
        )
        for comment in comments:
            output.write(f"{comment}\n")
        output.write(f"{order} {order} {len(entries)}\n")
        for row, column, value in entries:
            payload = "" if not value else " " + " ".join(value)
            output.write(f"{row + 1} {column + 1}{payload}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=pathlib.Path)
    parser.add_argument("output", type=pathlib.Path)
    parser.add_argument("--adjacent-swaps", type=int, required=True)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    if args.input.resolve() == args.output.resolve():
        parser.error("input and output must differ")
    try:
        comments, value_kind, symmetry, order, entries = read_coordinate_matrix(
            args.input
        )
        permutation = adjacent_permutation(
            order, args.adjacent_swaps, args.seed
        )
        transformed = permute_entries(entries, permutation, symmetry)
        write_coordinate_matrix(
            args.output,
            comments,
            value_kind,
            symmetry,
            order,
            transformed,
            args.adjacent_swaps,
            args.seed,
        )
    except (OSError, ValueError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
