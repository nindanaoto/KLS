#!/usr/bin/env python3
"""Freeze and qualify a KLS worktree with pinned local submodule sources."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess

from run_spice_campaign import sha, write_json


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('workspace', 'source', 'dependencies', 'image', 'output'):
        p.add_argument('--' + name, required=True, type=Path)
    p.add_argument('--profiles', nargs='+', choices=['native', 'portable', 'asan', 'tsan'], default=['native','portable','asan'])
    p.add_argument('--jobs', type=int, default=8)
    args = p.parse_args()
    args.output.mkdir(parents=True)
    source = args.output / 'source'
    source.mkdir()
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location('cleanup', root / 'bench/process_session_cleanup.py')
    cleanup = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cleanup)
    def git(*cmd):
        return subprocess.check_output(['git','-C',str(args.source),*cmd])
    # Copy the actual tracked/untracked worktree state, including candidate
    # files; never use a checkout's HEAD as a substitute for its dirty bytes.
    names = git('ls-files','--cached','--others','--exclude-standard','-z').decode().split('\0')
    names = [n for n in names if n and (args.source/n).is_file()]
    for name in names:
        dest=source/name;dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(args.source/name,dest)
    pins = {}
    for name in ['metis','gklib','scotch','spral']:
        dep=args.dependencies/'third_party'/name
        revision=subprocess.check_output(['git','-C',str(dep),'rev-parse','HEAD'],text=True).strip()
        expected=git('ls-tree','HEAD','third_party/'+name).decode().split()[2]
        if revision!=expected:
            raise ValueError(f'dependency {name}: expected {expected}, got {revision}')
        state=subprocess.check_output(['git','-C',str(dep),'status','--porcelain'],text=True)
        if state:
            raise ValueError(f'dirty dependency {dep}')
        target=source/'third_party'/name
        target.symlink_to(dep.resolve(),target_is_directory=True)
        pins[name]=dict(revision=revision,path=str(dep.resolve()))
    files={str(source/name):sha(source/name) for name in names}
    write_json(args.output/'source.json',dict(head=git('rev-parse','HEAD').decode().strip(),
        source=str(args.source),input_sha256=files,dependencies=pins,image_sha256=sha(args.image)))
    (args.output/'source.patch').write_bytes(git('diff','HEAD','--'))
    commands=[]
    env={k:v for k,v in os.environ.items() if not k.startswith(('LD_','KLS_','OMP_','GOMP_','OPENBLAS_','APPTAINERENV_','ASAN_','UBSAN_','TSAN_'))}
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',
               ASAN_OPTIONS='detect_leaks=1:halt_on_error=1',UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1',
               TSAN_OPTIONS='halt_on_error=1')
    def run(stage,cmd,expected=0):
        full=['apptainer','exec','--bind',f'{args.workspace}:/work','--bind',f'{args.workspace}:{args.workspace}',
              '--bind','/tmp:/tmp',str(args.image),*map(str,cmd)]
        with (args.output/(stage+'.log')).open('x') as f:
            child=subprocess.Popen(full,env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
            try:
                code=child.wait(timeout=1800)
            except BaseException:
                os.killpg(child.pid,signal.SIGKILL);child.wait()
                cleanup.drain_exited_apptainer_session(child.pid,args.image)
                raise
        commands.append(dict(stage=stage,command=full,pid=child.pid,exit_code=code,
            closure=cleanup.drain_exited_apptainer_session(child.pid,args.image)))
        write_json(args.output/'commands.json',commands)
        print(stage,code,flush=True)
        if code != expected:
            raise RuntimeError(f'{stage} failed; retained at {args.output}')
    for profile in args.profiles:
        build=args.output/profile
        flags=[]
        if profile=='portable':
            flags=['-DKLS_NATIVE_ARCH=OFF','-DKLS_ENABLE_HOSTED_SIMD=OFF']
        if profile in ('asan','tsan'):
            sanitize='address,undefined' if profile=='asan' else 'thread'
            flags=[f'-DCMAKE_C_FLAGS=-fsanitize={sanitize} -fno-omit-frame-pointer -fno-pie',
                   f'-DCMAKE_EXE_LINKER_FLAGS=-fsanitize={sanitize} -no-pie','-DCMAKE_C_FLAGS_RELEASE=-O2 -DNDEBUG']
        if profile=='tsan':
            # GCC's target_clones resolver runs before TSan's stack-hook
            # initialization. Keep memory instrumentation, disable only the
            # function entry/exit hooks, and prove race detection still works.
            flags[0] += ' --param=tsan-instrument-func-entry-exit=0'
            control=root/'tests/spice_tsan_control.c'
            for label,race in [('clean',False),('race',True)]:
                binary=args.output/('tsan-'+label)
                run('tsan-'+label+'-compile',['gcc','-O2','-g','-fsanitize=thread',
                    '--param=tsan-instrument-func-entry-exit=0','-fno-pie','-no-pie',
                    *(['-DRACE'] if race else []),control,'-pthread','-o',binary])
                run('tsan-'+label+'-run',[binary],66 if race else 0)
                log=(args.output/('tsan-'+label+'-run.log')).read_text()
                if 'CONTROL_MAIN_ENTERED' not in log or (race and 'WARNING: ThreadSanitizer: data race' not in log):
                    raise RuntimeError('TSan runtime control did not execute correctly')
        run(profile+'-configure',['cmake','-S',source,'-B',build,'-DCMAKE_BUILD_TYPE=Release',*flags])
        # Direct OpenMP reduction testing with uninstrumented libgomp reports
        # races even on the unchanged baseline. This scope tests the owned
        # pthread team and public solver ownership, without suppressing races.
        targets=['kls_accuracy_owner','kls_componentwise_team','kls_klu_allocator_ownership'] if profile=='tsan' else []
        run(profile+'-build',['cmake','--build',build,'-j'+str(args.jobs),*(['--target',*targets] if targets else [])])
        run(profile+'-test',['ctest','--test-dir',build,'--output-on-failure','-j1',*(['-R','^('+'|'.join(targets)+')$'] if targets else [])])
    if any(sha(path)!=h for path,h in files.items()):
        raise RuntimeError('frozen source changed during qualification')
    write_json(args.output/'completed.json',dict(passed=True,profiles=args.profiles,source_sha256=sha(args.output/'source.json'),
        commands_sha256=sha(args.output/'commands.json'),libraries={p:sha(args.output/p/'libkls.a') for p in args.profiles if (args.output/p/'libkls.a').exists()}))


if __name__=='__main__':
    main()
