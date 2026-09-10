#!/usr/bin/env python3
"""Fast end-to-end check for generated, application-managed tuning profiles."""

import argparse
import json
import os
import pathlib
import subprocess
import tempfile


def run(command: list[str], ok: bool = True) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True)
    if (result.returncode == 0) != ok:
        raise RuntimeError(
            f"unexpected exit {result.returncode}: {' '.join(command)}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--autotune", type=pathlib.Path, required=True)
    parser.add_argument("--bench", type=pathlib.Path, required=True)
    parser.add_argument("--source", type=pathlib.Path, required=True)
    args = parser.parse_args()

    tuner_text = (args.source / "bench" / "kls_autotune.c").read_text()
    forbidden = ("suitesparse_paper", "paper_manifest", "matrix-dir", "manifest")
    found = [word for word in forbidden if word in tuner_text.lower()]
    if found:
        raise RuntimeError(f"autotuner must not consume benchmark corpora: {found}")

    with tempfile.TemporaryDirectory(prefix="kls-autotune-test-") as directory:
        root = pathlib.Path(directory)
        profile = root / "profile.conf"
        matrix = root / "tiny.mtx"
        matrix.write_text(
            "%%MatrixMarket matrix coordinate real general\n"
            "3 3 7\n"
            "1 1 4\n1 2 -1\n2 1 -1\n2 2 4\n2 3 -1\n3 2 -1\n3 3 4\n"
        )
        cpu = min(os.sched_getaffinity(0))
        affinity = ["taskset", "-c", str(cpu)]
        run([*affinity, str(args.autotune), "--threads", "1", "--cpus",
             str(cpu), "--output", str(profile), "--self-test"])
        result = run([*affinity, str(args.bench), str(matrix), "--threads", "1",
                      "--tuning-profile", str(profile), "--analyze-only", "--json"])
        record = json.loads(result.stdout)
        if not record["tuning_profile_active"] or not record["tuning_profile_id"]:
            raise RuntimeError(f"profile provenance missing: {record}")

        damaged = root / "damaged.conf"
        text = profile.read_text()
        damaged.write_text(text.replace("row_batch_min_rows=8", "row_batch_min_rows=7"))
        run([*affinity, str(args.bench), str(matrix), "--threads", "1",
             "--tuning-profile", str(damaged), "--analyze-only"], ok=False)


if __name__ == "__main__":
    main()
