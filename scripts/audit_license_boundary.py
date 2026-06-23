#!/usr/bin/env python3
"""Audit KLS's MC64/HSL source boundary.

KLS may use MC64-style matching/scaling ideas and redistributable source, but
the LGPL distribution must not accidentally absorb HSL MC64 or solver-tree
copies that retain HSL-style redistribution limits.
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys


SCAN_ROOTS = (
    "src",
    "include",
    "bench",
    "tests",
    "third_party/suitesparse",
    "third_party/metis",
    "third_party/gklib",
    "third_party/scotch",
)

TEXT_SUFFIXES = {
    ".c",
    ".cc",
    ".cpp",
    ".cxx",
    ".f",
    ".f90",
    ".f95",
    ".h",
    ".hh",
    ".hpp",
    ".hxx",
    ".in",
    ".md",
    ".txt",
}

SPRAL_BSD_MARKERS = (
    "Redistribution and use in source and binary forms",
    "with or without",
    "modification, are permitted provided that",
    "Neither the name of the STFC nor the names of its contributors",
)

NOTICE_MARKERS = (
    "Existing MC64-style source can be used only when",
    "compatible with LGPL-2.1-or-later KLS",
    "HSL MC64",
)

SPRAL_SCALING_MARKERS = (
    "Hungarian code derives from HSL MC64 code",
    "This code is adapted from HSL_MC64",
    "This code is adapted from MC64",
)

REQUIRED_SPRAL_CMAKE_FRAGMENTS = (
    "${KLS_SPRAL_ROOT}/src/matrix_util.f90",
    "${KLS_SPRAL_ROOT}/src/scaling.f90",
    "${KLS_SPRAL_ROOT}/interfaces/C/scaling.f90",
)

DISALLOWED_SPRAL_CMAKE_FRAGMENTS = (
    "${KLS_SPRAL_ROOT}/src/match_order.f90",
    "${KLS_SPRAL_ROOT}/src/rutherford_boeing.f90",
    "${KLS_SPRAL_ROOT}/src/core_analyse.f90",
    "${KLS_SPRAL_ROOT}/src/ssids/",
    "${KLS_SPRAL_ROOT}/driver/",
)

RESTRICTED_MARKERS = (
    "HSL_MC64",
    "HSL MC64",
    "This code is adapted from MC64",
    "Subroutine MC64",
    "subroutine mc64",
    "MC64AD",
    "MC64BD",
    "MC64CD",
    "MC64DD",
    "MC64ED",
    "MC64FD",
    "redistribution is not permitted",
)


def read_text(path: pathlib.Path) -> str:
    return path.read_text(encoding="utf-8", errors="ignore")


def relpath(path: pathlib.Path, root: pathlib.Path) -> str:
    return path.relative_to(root).as_posix()


def tracked_files(root: pathlib.Path) -> list[pathlib.Path]:
    try:
        proc = subprocess.run(
            ["git", "-C", str(root), "ls-files", "--recurse-submodules"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        files: list[pathlib.Path] = []
        for scan_root in SCAN_ROOTS:
            directory = root / scan_root
            if directory.exists():
                files.extend(path for path in directory.rglob("*") if path.is_file())
        return sorted(files)

    return sorted(
        root / line.strip()
        for line in proc.stdout.splitlines()
        if line.strip()
    )


def under_scan_roots(path: pathlib.Path, root: pathlib.Path) -> bool:
    rel = relpath(path, root)
    return any(rel == scan_root or rel.startswith(f"{scan_root}/")
               for scan_root in SCAN_ROOTS)


def is_text_candidate(path: pathlib.Path) -> bool:
    return path.suffix.lower() in TEXT_SUFFIXES or path.name in {
        "CMakeLists.txt",
        "Makefile",
    }


def check_project_license(root: pathlib.Path, errors: list[str]) -> None:
    license_text = read_text(root / "LICENSE")
    if "GNU LESSER GENERAL PUBLIC LICENSE" not in license_text:
        errors.append("LICENSE is not the expected GNU LGPL text")

    notice_text = read_text(root / "THIRD_PARTY_NOTICES.md")
    for marker in NOTICE_MARKERS:
        if marker not in notice_text:
            errors.append(f"THIRD_PARTY_NOTICES.md is missing marker: {marker}")


def check_spral_boundary(root: pathlib.Path, errors: list[str]) -> None:
    spral_license = root / "third_party/spral/LICENCE"
    spral_scaling = root / "third_party/spral/src/scaling.f90"
    cmake = root / "CMakeLists.txt"

    if not spral_license.exists():
        errors.append("third_party/spral/LICENCE is missing")
    else:
        license_text = read_text(spral_license)
        for marker in SPRAL_BSD_MARKERS:
            if marker not in license_text:
                errors.append(f"SPRAL LICENCE is missing BSD marker: {marker}")

    if not spral_scaling.exists():
        errors.append("third_party/spral/src/scaling.f90 is missing")
    else:
        scaling_text = read_text(spral_scaling)
        for marker in SPRAL_SCALING_MARKERS:
            if marker not in scaling_text:
                errors.append(f"SPRAL scaling source is missing marker: {marker}")

    cmake_text = read_text(cmake)
    for fragment in REQUIRED_SPRAL_CMAKE_FRAGMENTS:
        if fragment not in cmake_text:
            errors.append(f"CMake SPRAL scaling build is missing: {fragment}")
    for fragment in DISALLOWED_SPRAL_CMAKE_FRAGMENTS:
        if fragment in cmake_text:
            errors.append(f"CMake SPRAL scaling build includes non-scaling code: {fragment}")


def check_restricted_markers(root: pathlib.Path, errors: list[str]) -> int:
    scanned = 0
    for path in tracked_files(root):
        if not path.exists() or not under_scan_roots(path, root) or not is_text_candidate(path):
            continue
        scanned += 1
        text = read_text(path)
        rel = relpath(path, root)
        lowered = text.lower()
        for marker in RESTRICTED_MARKERS:
            if marker.lower() in lowered:
                errors.append(
                    f"{rel}: contains restricted MC64/HSL marker {marker!r}; "
                    "use KLS code, pinned BSD SPRAL scaling, or a verified "
                    "LGPL-compatible source with notices"
                )
    return scanned


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Audit KLS LGPL/MC64 source boundary"
    )
    parser.add_argument(
        "--root",
        type=pathlib.Path,
        default=pathlib.Path(__file__).resolve().parents[1],
        help="KLS repository root",
    )
    args = parser.parse_args()

    root = args.root.resolve()
    errors: list[str] = []

    check_project_license(root, errors)
    check_spral_boundary(root, errors)
    scanned = check_restricted_markers(root, errors)

    if errors:
        for error in errors:
            print(f"license-boundary audit: {error}", file=sys.stderr)
        return 1

    print(
        "MC64/LGPL boundary ok: verified LGPL project notice, "
        "BSD SPRAL scaling subset, and scanned "
        f"{scanned} first-party/vendor source files"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
