"""Run the original 110-matrix regression protocol with a matched AUTO driver.

Build one benchmark object and link it to each qualified core using identical
remaining dependencies. Retain all 1,760 observations, including failures.
A complete, immutable matrix prefix may be continued with the identical
binaries after an independently retained teardown infrastructure failure.
"""
import argparse
import fcntl
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shlex
import signal
import statistics
import subprocess
import time


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')



def observation_id(row):
    return tuple(row[k] for k in ('config_id','matrix_id','trial','side','position'))


def load_completed_prefix(root, expected, required, failure_path):
    """Keep every recorded row; accept only a complete matrix boundary.

    This is continuation of identical frozen numerical work, not selection of
    successful observations. The abandoned unrecorded attempt must be the
    first baseline arm of the next matrix, with its failure retained outside
    the immutable campaign. No numeric result is inferred from that attempt.
    """
    root=root.resolve();prior=read(root/'metadata.json');failure=read(failure_path)
    assert not (root/'summary.json').exists(), 'completed campaigns are immutable'
    assert all(prior[k]==value for k,value in required.items()), 'numeric protocol changed'
    assert all(sha(p)==digest for p,digest in prior['input_sha256'].items()), 'frozen input changed'
    assert sha(root/'build-commands.json')==prior['build_commands_sha256']
    data=(root/'observations.jsonl').read_bytes()
    rows=[json.loads(line) for line in data.splitlines()]
    count=len(rows)
    assert count>0 and count<len(expected) and count%8==0, 'prefix must end at a full matrix'
    assert [observation_id(r) for r in rows]==expected[:count], 'prefix is incomplete or reordered'
    assert failure['campaign']==str(root) and failure['terminal_exit_code']!=0
    assert failure['completed_observations']==count and failure['ledger_sha256']==sha(root/'observations.jsonl')
    assert failure['failure']=='RuntimeError: unrecognized process in completed child session'
    stage=failure['stage'];next_id=expected[count]
    assert (stage['config'],stage['matrix'],stage['trial'],stage['side'],0)==next_id
    assert stage['side']=='baseline' and stage['trial']==0
    assert stage['numeric_outcome']=='not persisted before cleanup raised'
    assert failure['later_strict_cleanup']['session_id']==stage['session_id']
    assert failure['later_strict_cleanup']['all_session_processes_absent']
    # No timing is live during preparation; the caller also takes the timing lock.
    assert not Path('/proc',str(stage['session_id'])).exists()
    for row in rows:
        directory=Path(row['directory'])
        assert row['closure']['all_session_processes_absent']
        assert sha(directory/'stdout.log')==row['stdout_sha256']
        assert sha(directory/'stderr.log')==row['stderr_sha256']
        assert row['command'].count(prior['binaries'][row['side']])==1
    return prior,rows,data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('workspace', 'baseline-core', 'candidate-core', 'protocol', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--dependency-build', type=Path,
        help='Common dependency archives for both arms; each arm retains its own modified ordering archive')
    parser.add_argument('--continue-from', type=Path)
    parser.add_argument('--continuation-failure', type=Path)
    args = parser.parse_args()
    assert bool(args.continue_from)==bool(args.continuation_failure)
    workspace, out = args.workspace.resolve(), args.output.resolve()
    repo = Path(__file__).resolve().parents[1]
    lock = (workspace / '.spice-timing.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    busy = []
    for process in Path('/proc').iterdir():
        if not process.name.isdigit() or int(process.name) == os.getpid():
            continue
        try:
            name = (process / 'comm').read_text().strip()
            if name in ('cc1', 'cc1plus', 'f951', 'ctest', 'cmake', 'Xyce', 'kls_bench', 'baseline-bench', 'candidate-bench'):
                busy.append((int(process.name), name))
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
    assert not busy, ('timing requires idle build/test/solver processes', busy)
    protocol = args.protocol.resolve()
    image = workspace / 'KLS/build/kls-xyce-private.sif'
    cores = {'baseline': args.baseline_core.resolve(), 'candidate': args.candidate_core.resolve()}
    pins = {str(p): sha(p) for p in [Path(__file__), image, protocol / 'paper_runner.py',
        protocol / 'manifest.txt', protocol / 'configs.json', repo / 'bench/process_session_cleanup.py',
        repo / 'bench/process_session_settle.py', repo / 'scripts/audit_auto_broad.py',
        repo / 'tests/spice_session_settle.py', repo / 'tests/auto_broad_continuation.py']}
    for core in cores.values():
        completed, source = read(core / 'completed.json'), read(core / 'source.json')
        assert completed['passed'] and 'native' in completed['profiles']
        assert completed['source_sha256'] == sha(core / 'source.json')
        assert completed['commands_sha256'] == sha(core / 'commands.json')
        assert completed['libraries']['native'] == sha(core / 'native/libkls.a')
        pins.update(source['input_sha256'])
        pins.update({str(p): sha(p) for p in [core / 'completed.json', core / 'source.json', core / 'native/libkls.a']})
    candidate, baseline = cores['candidate'], cores['baseline']
    dependency_build = args.dependency_build.resolve() if args.dependency_build else baseline / 'native'
    assert sha(candidate / 'source/bench/kls_bench.c') == sha(baseline / 'source/bench/kls_bench.c')
    headers = lambda core: {str(p.relative_to(core / 'source/include')): sha(p)
        for p in (core / 'source/include').rglob('*.h')}
    assert headers(candidate) == headers(baseline)
    revisions = lambda core: {name: pin['revision'] for name, pin in read(core / 'source.json')['dependencies'].items()}
    assert revisions(candidate) == revisions(baseline)
    out.mkdir(exist_ok=False)
    spec = importlib.util.spec_from_file_location('paper', protocol / 'paper_runner.py')
    paper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(paper)
    spec = importlib.util.spec_from_file_location('cleanup', repo / 'bench/process_session_settle.py')
    cleanup = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cleanup)
    names = paper.manifest_names(protocol / 'manifest.txt')
    matrices = paper.resolve_matrices(workspace / 'KLS/data/suitesparse', names)
    configs = read(protocol / 'configs.json')
    assert len(names) == len(set(names)) == 110 and len(configs) == 2
    assert {(c['id'], c['cpus']) for c in configs} == {('headline_llc32_t8', '8-15'), ('headline_llc96_t8', '0-7')}
    assert all(c['threads'] == 8 and c['lifecycle_systems'] == 100 and
        c['values'] == 'entrywise' and c['amplitude'] == .001 and c['env'] == {} for c in configs)
    pins.update({str(p): sha(p) for p in matrices.values()})
    prefix = ['apptainer', 'exec', '--bind', f'{workspace}:/work', '--bind',
        f'{workspace}:{workspace}', '--bind', '/tmp:/tmp', '--pwd', '/work', str(image)]
    env = {k: v for k, v in paper.clean_environment({}).items() if not k.startswith(
        ('APPTAINERENV_', 'LD_', 'ASAN_', 'UBSAN_', 'TSAN_', 'GOMP_'))}
    commands = []
    build = candidate / 'native'
    flags_path = build / 'CMakeFiles/kls_bench.dir/flags.make'
    link_path = build / 'CMakeFiles/kls_bench.dir/link.txt'
    pins.update({str(p): sha(p) for p in [flags_path, link_path]})
    flags = []
    for line in flags_path.read_text().splitlines():
        if line.startswith(('C_DEFINES =', 'C_INCLUDES =', 'C_FLAGS =')):
            flags.extend(shlex.split(line.split(' = ', 1)[1]))

    def prepare(stage, command):
        assert all(sha(p) == h for p, h in pins.items())
        save(out / 'preparation-inputs.json', pins)
        with (out / (stage + '.log')).open('x') as log:
            child = subprocess.Popen(prefix + list(map(str, command)), cwd=build, env=env,
                stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            try:
                code = child.wait(timeout=1800)
            except BaseException:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
                cleanup.drain_exited_apptainer_session(child.pid, image)
                raise
        closure = cleanup.drain_exited_apptainer_session(child.pid, image)
        commands.append(dict(stage=stage, command=prefix + list(map(str, command)), exit_code=code, closure=closure))
        save(out / 'build-commands.json', commands)
        assert code == 0 and closure['all_session_processes_absent'], stage

    expected=[]
    for config in configs:
        for index,name in enumerate(names):
            for trial in range(4):
                order=['baseline','candidate'] if (index+trial)%2==0 else ['candidate','baseline']
                expected.extend((config['id'],name,trial,side,position) for position,side in enumerate(order))
    inherited_rows=[];continuation=None
    if args.continue_from:
        required=dict(expected=1760,passes=4,matrix_count=110,configs=configs,
            workspace=str(workspace),baseline_core=str(baseline),candidate_core=str(candidate),
            dependency_build=str(dependency_build),lifecycle_systems=100,timeout_seconds=900,residual_limit=1e-8,
            matrix_paths={name:str(path) for name,path in matrices.items()},
            manifest=str(protocol/'manifest.txt'),configs_path=str(protocol/'configs.json'),
            driver_sha256=sha(candidate/'source/bench/kls_bench.c'),matched_benchmark_object=True)
        prior,inherited_rows,data=load_completed_prefix(args.continue_from,expected,required,args.continuation_failure)
        # Reuse the actual object and binaries, rather than rebuilding presumed
        # equivalent files. Original rows and all their failed outcomes survive.
        obj=Path(prior['benchmark_object'])
        binaries={side:Path(path) for side,path in prior['binaries'].items()}
        commands=read(args.continue_from/'build-commands.json')
        (out/'build-commands.json').write_bytes((args.continue_from/'build-commands.json').read_bytes())
        (out/'observations.jsonl').write_bytes(data)
        pins.update(prior['input_sha256'])
        origin_files=[args.continue_from/'metadata.json',args.continue_from/'observations.jsonl',
            args.continue_from/'build-commands.json',args.continuation_failure]
        pins.update({str(p.resolve()):sha(p) for p in origin_files})
        continuation=dict(source=str(args.continue_from.resolve()),observations=len(inherited_rows),
            prefix_sha256=hashlib.sha256(data).hexdigest(),failure_receipt=str(args.continuation_failure.resolve()),
            identical_binaries=True,complete_matrix_boundary=True,all_recorded_failures_retained=True,
            abandoned_unrecorded_attempt_included_in_performance=False)
        save(out/'preparation-inputs.json',pins)
        print('Retained complete frozen prefix:',len(inherited_rows),'observations',flush=True)
    else:
        obj = out / 'bench.o'
        prepare('compile', ['gcc', *flags, '-c', candidate / 'source/bench/kls_bench.c', '-o', obj])
        pins[str(obj)] = sha(obj)
        binaries = {}
        for side, core in cores.items():
            binary = out / (side + '-bench')
            link = shlex.split(link_path.read_text())
            common = ('libkls_metis.a', 'libkls_gklib.a', 'lib/libscotch.a', 'lib/libscotcherr.a', 'libkls_spral_scaling.a')
            link = [str(obj) if p == 'CMakeFiles/kls_bench.dir/bench/kls_bench.c.o' else
                str(core / 'native' / p) if p in ('libkls.a', 'libkls_suitesparse.a') else
                str(dependency_build / p) if p in common else str(binary) if p == 'kls_bench' else
                '-Wl,--dependency-file=' + str(out / (side + '.link.d')) if p.startswith('-Wl,--dependency-file=') else p for p in link]
            for item in link[1:]:
                path = Path(item) if item.startswith('/') else build / item
                if path.is_file() and path not in (obj, binary):
                    pins[str(path)] = sha(path)
            prepare(side + '-link', link)
            binaries[side] = binary
            pins[str(binary)] = sha(binary)
    metadata = dict(matrix_count=110, passes=4, arms=list(cores), configs=configs,
        workspace=str(workspace),
        baseline_core=str(baseline), candidate_core=str(candidate), dependency_build=str(dependency_build),
        ordering_archives={side: str(core / 'native/libkls_suitesparse.a') for side, core in cores.items()},
        lifecycle_systems=100, expected=1760, timeout_seconds=900, residual_limit=1e-8,
        binaries={side: str(path) for side, path in binaries.items()},
        matrix_paths={name: str(path) for name, path in matrices.items()},
        manifest=str(protocol / 'manifest.txt'), configs_path=str(protocol / 'configs.json'),
        driver_sha256=sha(candidate / 'source/bench/kls_bench.c'), matched_benchmark_object=True,
        benchmark_object=str(obj), build_commands_sha256=sha(out / 'build-commands.json'),
        input_sha256=pins, scope='Qualified main versus structural AUTO; complete original matrix manifest; shared benchmark object and common dependencies; modified ordering is part of each core; no successful-subset aggregate.')
    if continuation:metadata['continuation']=continuation
    save(out / 'metadata.json', metadata)
    rows = list(inherited_rows)
    prefix_count=len(inherited_rows);visited=0
    for config in configs:
        for index, name in enumerate(names):
            for trial in range(4):
                order = ['baseline', 'candidate'] if (index + trial) % 2 == 0 else ['candidate', 'baseline']
                for position, side in enumerate(order):
                    visited+=1
                    if visited<=prefix_count:continue
                    dest = out / config['id'] / name / ('pass' + str(trial)) / side
                    dest.mkdir(parents=True)
                    matrix = Path('/work') / matrices[name].relative_to(workspace)
                    command = ['taskset', '-c', config['cpus'], *prefix,
                        *paper.command_for('kls', binaries[side], matrix, config)]
                    timed_out, start = False, time.monotonic()
                    with (dest / 'stdout.log').open('x') as stdout, (dest / 'stderr.log').open('x') as stderr:
                        child = subprocess.Popen(command, env=env, stdout=stdout, stderr=stderr, start_new_session=True)
                        save(dest/'started.json',dict(pid=child.pid,command=command,timeout_seconds=900))
                        try:
                            code = child.wait(timeout=900)
                        except subprocess.TimeoutExpired:
                            timed_out = True
                            os.killpg(child.pid, signal.SIGKILL)
                            child.wait()
                            code = None
                    save(dest/'numeric-exit.json',dict(pid=child.pid,returncode=code,timed_out=timed_out,
                        numeric_process_wall_seconds=time.monotonic()-start))
                    try:
                        closure = cleanup.drain_exited_apptainer_session(child.pid, image)
                    except BaseException as error:
                        save(dest/'teardown-failure.json',dict(pid=child.pid,error=repr(error)))
                        raise
                    save(dest/'closure.json',closure)
                    result, valid, reason = None, False, 'timeout' if timed_out else 'process-failed'
                    if code == 0:
                        try:
                            result = paper.parse_json_output((dest / 'stdout.log').read_text())
                            valid, reason = paper.result_valid(result, 1e-8, 100)
                        except (ValueError, TypeError, KeyError) as error:
                            reason = 'parse-error:' + repr(error)
                    assert closure['all_session_processes_absent']
                    row = dict(config_id=config['id'], matrix_id=name, trial=trial, side=side,
                        position=position, directory=str(dest), command=command, returncode=code,
                        timed_out=timed_out, wall_seconds=time.monotonic() - start, valid=valid,
                        reason=reason, result=result, closure=closure,
                        stdout_sha256=sha(dest / 'stdout.log'), stderr_sha256=sha(dest / 'stderr.log'),
                        lifecycle_receipt_sha256={name:sha(dest/name) for name in
                            ('started.json','numeric-exit.json','closure.json')})
                    rows.append(row)
                    with (out / 'observations.jsonl').open('a') as stream:
                        stream.write(json.dumps(row) + '\n')
                    print(len(rows), '/1760', config['id'], name, trial, side, valid, reason, flush=True)
    assert all(sha(p) == h for p, h in pins.items())
    indexed = {(r['config_id'], r['matrix_id'], r['trial'], r['side']): r for r in rows}
    assert len(indexed) == len(rows) == 1760
    groups, regressions = [], []
    for config in configs:
        for name in names:
            ratios = []
            for trial in range(4):
                a, b = (indexed[config['id'], name, trial, side] for side in ('candidate', 'baseline'))
                if b['valid'] and not a['valid']:
                    regressions.append(dict(config=config['id'], matrix=name, trial=trial, reason=a['reason']))
                if a['valid'] and b['valid']:
                    ratios.append(a['result']['measured_lifecycle_seconds'] / b['result']['measured_lifecycle_seconds'])
            groups.append(dict(config=config['id'], matrix=name, valid_pairs=len(ratios),
                ratios=ratios, median=statistics.median(ratios) if len(ratios) == 4 else None))
    complete = all(g['median'] is not None for g in groups)
    save(out / 'summary.json', dict(observations=len(rows), expected=1760, changed_inputs=[],
        coverage_regressions=regressions, groups=groups, all_matrices_qualified=complete,
        candidate_over_baseline_geomean=math.exp(statistics.mean(math.log(g['median']) for g in groups)) if complete else None,
        invalid=[r for r in rows if not r['valid']], metadata_sha256=sha(out / 'metadata.json')))


if __name__ == '__main__':
    main()
