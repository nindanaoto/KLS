"""Smoke-check the versioned public-call lifecycle timing output."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bench', type=Path, required=True)
    args = parser.parse_args()
    env = dict(os.environ, KLS_BENCH_VERIFY_EACH_REFACTOR='1', OMP_NUM_THREADS='1')
    with tempfile.TemporaryDirectory(prefix='kls-bench-timing-') as directory:
        matrix = Path(directory) / 'diagonal.mtx'
        matrix.write_text('%%MatrixMarket matrix coordinate real general\n'
                          '3 3 3\n1 1 2\n2 2 3\n3 3 4\n')
        for systems in (1, 4, 0):
            command = [str(args.bench.resolve()), str(matrix), '--threads', '1', '--json']
            command += (['--lifecycle-systems', str(systems), '--refactor-values', 'entrywise']
                        if systems else ['--repeat', '1', '--refactor-repeat', '1'])
            started = time.monotonic()
            result = subprocess.run(command, env=env, capture_output=True, text=True,
                                    check=True, timeout=120)
            elapsed = time.monotonic() - started
            row = json.loads(result.stdout)
            assert row['status'] == 0 and row['lifecycle_timing_schema'] == 2
            assert row['relative_residual_l2'] <= 1e-8
            public, internal = row['measured_lifecycle_seconds'], row['internal_lifecycle_seconds']
            if systems:
                assert row['lifecycle_systems'] == systems and row['refactor_repeat'] == systems - 1
                assert math.isfinite(public) and 0 < public <= elapsed
                assert math.isfinite(internal) and internal >= 0
            else:
                assert public == internal == -1
    print('lifecycle timing schema checks passed')


if __name__ == '__main__':
    main()
