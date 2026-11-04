"""Build pinned standalone peer drivers from a qualified candidate snapshot."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('workspace', 'candidate-core', 'image', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    workspace, core, out = (p.resolve() for p in
                            (args.workspace, args.candidate_core, args.output))
    image = args.image.resolve()
    done = json.loads((core / 'completed.json').read_text())
    source = json.loads((core / 'source.json').read_text())
    assert done['passed'] and 'native' in done['profiles']
    assert done['source_sha256'] == sha(core / 'source.json')
    assert done['commands_sha256'] == sha(core / 'commands.json')
    assert done['libraries']['native'] == sha(core / 'native/libkls.a')
    assert source['image_sha256'] == sha(image)
    pins = dict(source['input_sha256'])
    for path in (Path(__file__).resolve(), image, core / 'completed.json',
                 core / 'source.json', core / 'commands.json', core / 'native/libkls.a'):
        pins[str(path)] = sha(path)
    repo = Path(__file__).resolve().parents[1]
    cleanup_path = repo / 'bench/process_session_cleanup.py'
    pins[str(cleanup_path)] = sha(cleanup_path)
    spec = importlib.util.spec_from_file_location('cleanup', cleanup_path)
    cleanup = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cleanup)
    out.mkdir(exist_ok=False)
    prefix = ['apptainer', 'exec', '--bind', f'{workspace}:{workspace}',
              '--bind', '/tmp:/tmp', str(image)]
    env = {k: v for k, v in os.environ.items() if not k.startswith(
        ('LD_', 'KLS_', 'OMP_', 'GOMP_', 'APPTAINERENV_', 'ASAN_', 'UBSAN_', 'TSAN_'))}
    commands, binaries, libraries = [], {}, {}
    for peer, root, library, driver in (
        ('ck', workspace / 'cktso', 'rocky8_x64_gcc850/libcktso.so', 'cktso_compare.cpp'),
        ('st', workspace / 'SubtreeLU', 'lib/libsubtree_lu.so', 'subtreelu_compare.cpp')):
        lib = root / library
        header = root / 'include' / ('cktso.h' if peer == 'ck' else 'subtree_lu.h')
        src = core / 'source/bench' / driver
        for path in (lib, header, src, core / 'source/bench/bench_value_sequence.h'):
            pins[str(path)] = sha(path)
        binary = out / (peer + '-bench')
        command = prefix + ['g++', '-O3', '-march=native', '-std=c++17',
            '-I' + str(root / 'include'), str(src), str(lib), '-lm',
            '-Wl,-rpath,' + str(lib.parent), '-o', str(binary)]
        with (out / (peer + '-build.log')).open('x') as log:
            child = subprocess.Popen(command, env=env, stdout=log,
                stderr=subprocess.STDOUT, start_new_session=True)
            try:
                code = child.wait(timeout=300)
            except BaseException:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
                cleanup.drain_exited_apptainer_session(child.pid, image)
                raise
        closure = cleanup.drain_exited_apptainer_session(child.pid, image)
        commands.append(dict(peer=peer, command=command, exit_code=code, closure=closure))
        (out / 'commands.json').write_text(json.dumps(commands, indent=2) + '\n')
        assert code == 0 and closure['all_session_processes_absent']
        binaries[peer], libraries[peer] = str(binary), str(lib)
        pins[str(binary)] = sha(binary)
    assert all(sha(path) == digest for path, digest in pins.items())
    receipt = dict(passed=True, candidate_core=str(core), binaries=binaries,
        libraries=libraries, input_sha256=pins,
        commands_sha256=sha(out / 'commands.json'),
        value_sequence_sha256=sha(core / 'source/bench/bench_value_sequence.h'))
    (out / 'build.json').write_text(json.dumps(receipt, indent=2) + '\n')


if __name__ == '__main__':
    main()
