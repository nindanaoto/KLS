#!/usr/bin/env python3
"""Relink pinned Xyce artifacts with fresh adapters and a selected KLS build."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess

from run_spice_campaign import sha, write_json


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('workspace', 'image', 'artifacts', 'provenance', 'generated',
                 'trilinos', 'xyce-source', 'kls-source', 'kls-build', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--common-archive', type=Path,
                   help='Reuse an already rebuilt common archive; headers must match its provenance')
    p.add_argument('--dependency-build', type=Path,
                   help='Use this common KLS dependency build for matched links')
    p.add_argument('--coarse-profile', action='store_true',
                   help='Diagnostic only: resolve the coarse profiler allocation wrappers')
    args = p.parse_args()
    builder_hash = sha(__file__)
    args.output.mkdir(parents=True)
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location('cleanup', root / 'bench/process_session_cleanup.py')
    cleanup = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cleanup)
    records = []
    inputs = [Path(__file__), root / 'scripts/run_spice_campaign.py', root / 'bench/spice_link.json',
              args.image, args.kls_source / 'include/kls/kls.h']
    inputs += list(args.kls_build.glob('*.a')) + list((args.kls_build / 'lib').glob('*.a'))
    env = {k: v for k, v in os.environ.items() if not k.startswith(('KLS_', 'LD_', 'APPTAINERENV_'))}

    def execute(stage, command):
        cmd = ['apptainer', 'exec', '--bind', f'{args.workspace}:/work',
               '--bind', f'{args.workspace}:{args.workspace}', '--bind', '/tmp:/tmp',
               str(args.image), *map(str, command)]
        with (args.output / (stage + '.log')).open('x') as log:
            child = subprocess.Popen(cmd, env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            code = child.wait()
        closure = cleanup.drain_exited_apptainer_session(child.pid, args.image)
        records.append(dict(stage=stage, command=cmd, pid=child.pid, exit_code=code, closure=closure))
        write_json(args.output / 'commands.json', records)
        if code:
            raise RuntimeError(f'{stage} failed: see {args.output / (stage + ".log")}')

    if args.common_archive:
        archive = args.common_archive.resolve()
        common = json.loads((archive.parent / 'build.json').read_text())
        header = str(args.kls_source / 'include/kls/kls.h')
        if sha(header) != common['public_header_sha256']:
            raise ValueError('common Xyce archive has a different KLS public header')
        if sha(archive) != common['output_sha256'][str(archive)]:
            raise ValueError('common archive changed after qualification')
        inputs += [archive, archive.parent / 'build.json']
    else:
        archive = args.output.resolve() / 'libxyce.a'
        original = args.artifacts / 'src/libxyce.a'
        inputs.append(original)
        shutil.copy2(original, archive)
        flags_file = args.provenance / 'src/CMakeFiles/XyceLib.dir/flags.make'
        cache_file = args.provenance / 'CMakeCache.txt'
        flags_text, cache = flags_file.read_text(), cache_file.read_text()
        inputs += [flags_file, cache_file]
        old_generated = next(x.split('=', 1)[1] for x in cache.splitlines() if x.startswith('CMAKE_CACHEFILE_DIR:'))
        flags = []
        for key in ('CXX_DEFINES', 'CXX_INCLUDES', 'CXX_FLAGS'):
            flags.extend(shlex.split(next(x.split(' = ', 1)[1] for x in flags_text.splitlines() if x.startswith(key + ' = '))))
        replaced = []
        for flag in flags:
            flag = flag.replace('/work/Xyce', str(args.xyce_source)).replace(old_generated, str(args.generated))
            flag = flag.replace('/tmp/kls-trilinos-install', str(args.trilinos))
            if flag.startswith('-I') and flag.endswith('/include') and 'kls-' in flag.lower() and 'trilinos' not in flag.lower():
                flag = '-I' + str(args.kls_source / 'include')
            replaced.append(flag)
        for name in ('N_LAS_KLSSolver', 'N_LAS_ThirdPartySolver'):
            source = args.xyce_source / 'src/LinearAlgebraServicesPKG' / (name + '.C')
            obj = args.output.resolve() / (name + '.C.o')
            inputs.append(source)
            execute(name, ['g++', *replaced, '-c', source, '-o', obj])
            execute(name + '-archive', ['ar', 'r', archive, obj])
    values = dict(artifacts=args.artifacts, xyce_archive=archive, kls_build=args.kls_build,
                  trilinos=args.trilinos, workspace=args.workspace, output=args.output.resolve() / 'Xyce')
    template = json.loads((root / 'bench/spice_link.json').read_text())['command']
    if args.dependency_build:
        template = [x.replace('{kls_build}', str(args.dependency_build))
                    if x.startswith('{kls_build}') and not x.endswith('/libkls.a') else x for x in template]
    command = [x.format(**values) for x in template]
    if args.coarse_profile:
        command += ['-Wl,--wrap=' + name for name in
                    ('malloc', 'calloc', 'realloc', 'free', 'aligned_alloc', 'posix_memalign')]
    # Include every explicitly linked file in the input manifest.
    inputs += [Path(x) for x in command if x.startswith('/') and Path(x).is_file() and x != str(values['output'])]
    hashes = {str(x.resolve()): sha(x) for x in inputs}
    execute('link', command)
    if sha(__file__) != builder_hash or any(sha(path) != h for path, h in hashes.items()):
        raise RuntimeError('input changed during link')
    write_json(args.output / 'build.json', dict(input_sha256=hashes,
        public_header_sha256=sha(args.kls_source / 'include/kls/kls.h'),
        output_sha256={str(values['output']): sha(values['output']), str(archive): sha(archive)}, commands=records))


if __name__ == '__main__':
    main()
