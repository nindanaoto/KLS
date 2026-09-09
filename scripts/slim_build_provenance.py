"""Fail-closed provenance checks for paired Unix Makefiles benchmark builds.

Inspect linked targets, not just kls.dir: SPRAL's Fortran flags affect matching
even when every KLS C option agrees. Frozen binaries alone are not build proof.
"""
from __future__ import annotations

import difflib
import hashlib
import json
import pathlib
import re
import shlex
import subprocess


def digest(path: pathlib.Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def output(command: list[str], cwd: pathlib.Path | None = None) -> str:
    return subprocess.check_output(command, cwd=cwd, text=True).strip()


def normalize(text: str, source: pathlib.Path, build: pathlib.Path) -> str:
    # Replace build first because it normally lives beneath source. Do not
    # sort/deduplicate options: the last optimization flag wins.
    return text.replace(str(build), "<BUILD>").replace(str(source), "<SOURCE>")


def collect_build_provenance(
        binary: pathlib.Path, *, dependency_patch_sha256: str | None = None) -> dict:
    """Collect evidence, optionally recording an explicitly reviewed dependency patch.

    The opt-in hash covers output()'s stripped UTF-8 git binary diff. It does
    not make changed dependencies equivalent: the returned proof records it,
    and require_matching_builds continues to reject dirty dependencies.
    """
    binary = binary.resolve(strict=True)
    build = binary.parent
    cache = (build / "CMakeCache.txt").read_text()
    settings = dict(re.findall(r"^([^/#\n][^:\n]*):[^=\n]+=(.*)$", cache, re.M))
    if settings.get("CMAKE_GENERATOR") != "Unix Makefiles":
        raise ValueError("build provenance requires CMake Unix Makefiles metadata")
    source = pathlib.Path(settings["CMAKE_HOME_DIRECTORY"]).resolve(strict=True)
    norm = lambda text: normalize(text, source, build)

    # Resolve archive outputs relative to each target's actual link directory,
    # including nested SCOTCH targets. Only reachable targets enter the proof.
    targets = {}
    for link in build.rglob("link.txt"):
        lines = link.read_text().splitlines()
        if not lines:
            continue
        command = shlex.split(lines[0])
        if not command:
            continue
        if "-o" in command:
            artifact = command[command.index("-o") + 1]
        elif pathlib.Path(command[0]).name.endswith("ar"):
            artifact = command[2]
        else:
            continue
        link_cwd = link.parent.parent.parent
        targets[(link_cwd / artifact).resolve()] = (link, link_cwd, command)
    if binary not in targets:
        raise ValueError(f"missing link provenance for {binary}")

    compiled, external, compilers = {}, {}, {}

    def visit(artifact: pathlib.Path) -> None:
        if artifact not in targets:
            checksum = digest(artifact)
            if artifact.name in external and external[artifact.name] != checksum:
                raise ValueError(f"ambiguous external dependency name: {artifact.name}")
            external[artifact.name] = checksum
            return
        link, link_cwd, command = targets[artifact]
        key = str(link.parent.relative_to(build))
        if key in compiled:
            return
        flags_path = link.parent / "flags.make"
        flags = flags_path.read_text()
        if artifact.stat().st_mtime_ns < max(
                flags_path.stat().st_mtime_ns, link.stat().st_mtime_ns):
            raise ValueError(f"stale artifact; rebuild {artifact}")
        options = dict(re.findall(r"^(\w+_(?:FLAGS|DEFINES|INCLUDES)) = (.*)$",
                                  flags, re.M))
        languages = re.findall(r"^# compile (\w+) with (.+)$", flags, re.M)
        if not languages or not options:
            raise ValueError(f"missing compiler/flags evidence in {flags_path}")
        for language, compiler in languages:
            if f"{language}_FLAGS" not in options:
                raise ValueError(f"missing {language} flags in {flags_path}")
            identity = {"command": compiler,
                        "version": output(shlex.split(compiler) + ["--version"])}
            if language in compilers and compilers[language] != identity:
                raise ValueError(f"inconsistent {language} compilers within build")
            compilers[language] = identity
        compiled[key] = {"options": {k: norm(v) for k, v in options.items()},
                         "link": norm(link.read_text())}
        for token in command[1:]:
            if token.startswith("-") or not re.search(r"\.(?:a|o|so(?:\.\d+)*)$", token):
                continue
            dependency = (link_cwd / token).resolve(strict=True)
            if dependency == artifact:
                continue
            if binary.stat().st_mtime_ns < dependency.stat().st_mtime_ns:
                raise ValueError(f"stale benchmark; relink after {dependency}")
            # Object files compiled by this target are covered by its flags;
            # externally supplied objects (MT-METIS) require content hashes.
            if dependency.suffix == ".o" and dependency.is_relative_to(link.parent):
                continue
            visit(dependency)

    visit(binary)
    # -l dependencies must also be identical, including runtime BLAS/Fortran.
    loaded = {}
    for line in output(["ldd", str(binary)]).splitlines():
        if "not found" in line:
            raise ValueError(f"unresolved runtime dependency: {line.strip()}")
        match = re.search(r"(?:=>\s*)?(/\S+)\s+\(", line)
        if match:
            library = pathlib.Path(match[1]).resolve(strict=True)
            loaded[library.name] = digest(library)
    if not loaded:
        raise ValueError("missing runtime dependency provenance")
    tree = output(["git", "ls-tree", "HEAD", "third_party"], source)
    if not tree:
        raise ValueError("missing third_party revision provenance")
    entries = output(["git", "ls-tree", "-r", "HEAD", "third_party"], source)
    submodules = {}
    for line in entries.splitlines():
        mode, kind, rest = line.split(None, 2)
        revision, name = rest.split("\t", 1)
        if mode != "160000":
            continue
        dependency = (source / name).resolve(strict=True)
        actual = output(["git", "rev-parse", "HEAD"], dependency)
        if actual != revision or output(
                ["git", "diff", "--no-ext-diff", "HEAD", "--"], dependency):
            raise ValueError(f"dependency revision/contents differ from gitlink: {name}")
        submodules[name] = actual
    changed = output(["git", "diff", "--name-only", "HEAD", "--", "third_party"], source)
    # Historical worktrees may symlink their gitlinks to the identical clean
    # dependency checkout. Check those actual revisions above, not link text.
    dependency_patch = {}
    if set(changed.splitlines()) - submodules.keys():
        if dependency_patch_sha256 is None:
            raise ValueError("dirty dependency sources; cannot certify paired builds")
        patch = output(["git", "diff", "--no-ext-diff", "--no-textconv",
                        "--binary", "HEAD", "--", "third_party"], source)
        actual = hashlib.sha256(patch.encode("utf-8")).hexdigest()
        if not patch or actual != dependency_patch_sha256:
            raise ValueError("dependency patch does not match reviewed SHA-256")
        dependency_patch = {"dependency_patch_sha256": actual}
    elif dependency_patch_sha256 is not None:
        raise ValueError("expected dependency patch is missing")
    return {**dependency_patch, "schema": 1, "generator": settings["CMAKE_GENERATOR"],
            "build_type": settings.get("CMAKE_BUILD_TYPE", ""),
            "compilers": compilers, "targets": compiled,
            "dependency_tree": tree, "submodules": submodules, "external_objects": external,
            "runtime_libraries": loaded}


def require_matching_builds(binaries: dict[str, pathlib.Path]) -> dict:
    try:
        proof = {side: collect_build_provenance(path) for side, path in binaries.items()}
        before = json.dumps(proof["before"], indent=2, sort_keys=True).splitlines()
        after = json.dumps(proof["after"], indent=2, sort_keys=True).splitlines()
        if before != after:
            difference = "\n".join(difflib.unified_diff(
                before, after, fromfile="before build", tofile="after build"))
            raise ValueError("paired build provenance differs:\n" + difference)
        return proof
    except (OSError, KeyError, subprocess.CalledProcessError, ValueError) as error:
        raise ValueError(f"benchmark preflight failed: {error}") from error
