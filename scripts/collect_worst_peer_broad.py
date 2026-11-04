"""Collect fresh peer evidence for every original standalone regression.

This collector makes no acceptance decision. Its raw ledger must be checked
independently before the peer-aware broad guard may admit any exception.
"""
import argparse
import fcntl
import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

from worst_peer_broad_evidence import bind_original, read, sha


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('campaign', 'peer-build', 'license', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    campaign, build, out = (p.resolve() for p in
                            (args.campaign, args.peer_build, args.output))
    repo = Path(__file__).resolve().parents[1]
    original = read(campaign / 'metadata.json')
    workspace = Path(original['workspace'])
    lock = (workspace / '.spice-timing.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    # Prevent a build or solver outside this lock from contaminating timings.
    busy = []
    for process in Path('/proc').iterdir():
        if not process.name.isdigit() or int(process.name) == os.getpid():
            continue
        try:
            name = (process / 'comm').read_text().strip()
            if name in ('cc1', 'cc1plus', 'f951', 'ctest', 'cmake', 'Xyce',
                        'kls_bench', 'baseline-bench', 'candidate-bench',
                        'ck-bench', 'st-bench'):
                busy.append((process.name, name))
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            pass
    assert not busy, ('timing requires idle builds/tests/solvers', busy)
    out.mkdir(exist_ok=False)
    # Re-run the original independent auditor; a stale receipt is insufficient.
    subprocess.run([sys.executable, str(repo / 'scripts/audit_auto_broad.py'),
                    '--campaign', str(campaign), '--output', str(out / 'original-audit')],
                   check=True)
    bound = bind_original(campaign, out / 'original-audit/audit.json')
    peers = read(build / 'build.json')
    assert peers['passed'] is True
    assert peers['candidate_core'] == original['candidate_core']
    assert peers['commands_sha256'] == sha(build / 'commands.json')
    commands = read(build / 'commands.json')
    assert [c['peer'] for c in commands] == ['ck', 'st']
    assert all(c['exit_code'] == 0 and c['closure']['all_session_processes_absent']
               for c in commands)
    value_header = Path(original['candidate_core']) / 'source/bench/bench_value_sequence.h'
    assert peers['value_sequence_sha256'] == sha(value_header)
    pins = dict(bound['input_sha256'])
    for path, digest in peers['input_sha256'].items():
        assert path not in pins or pins[path] == digest
        pins[path] = digest
    paper_path = Path(original['manifest']).parent / 'paper_runner.py'
    cleanup_path = repo / 'bench/process_session_cleanup.py'
    image = workspace / 'KLS/build/kls-xyce-private.sif'
    for path in (Path(__file__).resolve(), build / 'build.json', build / 'commands.json',
                 paper_path, cleanup_path, image):
        pins[str(path)] = sha(path)
    assert all(sha(p) == h for p, h in pins.items())
    paper = load_module('paper', paper_path)
    cleanup = load_module('cleanup', cleanup_path)
    binaries = dict(candidate=bound['candidate'], **peers['binaries'])
    assert set(binaries) == {'candidate', 'ck', 'st'}
    assert all(pins.get(path) == sha(path) for path in binaries.values())
    env = {k: v for k, v in paper.clean_environment({}).items()
           if not k.startswith(('APPTAINERENV_', 'LD_', 'ASAN_', 'UBSAN_',
                                'TSAN_', 'GOMP_'))}
    license_dir = args.license.resolve().parent
    assert args.license.name == 'cktso.lic' and args.license.is_file()
    prefix = ['apptainer', 'exec', '--bind', f'{workspace}:/work', '--bind',
              f'{workspace}:{workspace}', '--bind', '/tmp:/tmp', '--bind',
              f'{license_dir}:/peer-license:ro', '--pwd', '/peer-license', str(image)]
    targets = [(g['config'], g['matrix']) for g in bound['targets']]
    metadata = dict(original_campaign=str(campaign), peer_build=str(build),
                    binaries=binaries, targets=targets, expected=len(targets) * 12,
                    input_sha256=pins, configs=original['configs'],
                    matrix_paths=original['matrix_paths'], prefix=prefix,
                    environment_overrides={k: env[k] for k in (
                        'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
                        'KLS_BENCH_VERIFY_EACH_REFACTOR', 'BENCH_VERIFY_EACH_REFACTOR')},
                    timeout_seconds=900, lifecycle_systems=100,
                    residual_limit=1e-8, acceptance_decision=False)
    save(out / 'metadata.json', metadata)
    configs = {c['id']: c for c in original['configs']}
    count = 0
    for index, (config_id, name) in enumerate(targets):
        config = configs[config_id]
        matrix = Path(original['matrix_paths'][name])
        for trial in range(4):
            # Forward/reverse pairs balance every candidate/peer comparison.
            order = ['candidate', 'ck', 'st']
            if (index + trial) % 2:
                order.reverse()
            for position, side in enumerate(order):
                dest = out / config_id / name / str(trial) / side
                dest.mkdir(parents=True)
                command = ['taskset', '-c', config['cpus'], *prefix,
                           *paper.command_for('kls' if side == 'candidate' else side,
                                              Path(binaries[side]), matrix, config)]
                start, timed_out = time.monotonic(), False
                with (dest / 'stdout.log').open('x') as stdout, (dest / 'stderr.log').open('x') as stderr:
                    child = subprocess.Popen(command, env=env, stdout=stdout,
                        stderr=stderr, start_new_session=True)
                    save(dest / 'started.json', dict(pid=child.pid, command=command))
                    try:
                        code = child.wait(timeout=900)
                    except subprocess.TimeoutExpired:
                        timed_out = True
                        os.killpg(child.pid, signal.SIGKILL)
                        child.wait()
                        code = None
                save(dest / 'numeric-exit.json', dict(returncode=code, timed_out=timed_out,
                    wall_seconds=time.monotonic() - start))
                closure = cleanup.drain_exited_apptainer_session(child.pid, image)
                assert closure['all_session_processes_absent']
                save(dest / 'closure.json', closure)
                result, valid = None, False
                reason = 'timeout' if timed_out else 'process-failed'
                if code == 0:
                    try:
                        result = paper.parse_json_output((dest / 'stdout.log').read_text())
                        valid, reason = paper.result_valid(result, 1e-8, 100)
                    except (ValueError, TypeError, KeyError) as error:
                        reason = 'parse-error:' + repr(error)
                row = dict(config=config_id, matrix=name, trial=trial, side=side,
                    position=position, directory=str(dest), command=command,
                    returncode=code, timed_out=timed_out, valid=valid, reason=reason,
                    result=result, files_sha256={n: sha(dest / n) for n in
                        ('stdout.log', 'stderr.log', 'started.json', 'numeric-exit.json', 'closure.json')})
                with (out / 'observations.jsonl').open('a') as stream:
                    stream.write(json.dumps(row) + '\n')
                count += 1
                print(count, '/', metadata['expected'], config_id, name, trial,
                      side, valid, reason, flush=True)
    (out / 'observations.jsonl').touch(exist_ok=True)
    assert count == metadata['expected']
    assert all(sha(p) == h for p, h in pins.items())
    save(out / 'completed.json', dict(observations=count, acceptance_decision=False,
         metadata_sha256=sha(out / 'metadata.json'),
         ledger_sha256=sha(out / 'observations.jsonl')))


if __name__ == '__main__':
    main()
