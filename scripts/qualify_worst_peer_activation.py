"""Qualify observed-parallel and forest execution on frozen circuit builds.

All times collected here are diagnostic. Ordinary performance campaigns use
candidate-link/Xyce only, never activation/Xyce.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import sys


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def save(path,value):
    Path(path).write_text(json.dumps(value,indent=2)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('workspace','candidate-core','baseline-core','output'):
        p.add_argument('--'+name,required=True,type=Path)
    a=p.parse_args();workspace=a.workspace.resolve();out=a.output.resolve()
    cores={'candidate':a.candidate_core.resolve(),'baseline':a.baseline_core.resolve()}
    repo=Path(__file__).resolve().parents[1];image=workspace/'KLS/build/kls-xyce-private.sif'
    pins={str(path):sha(path) for path in (Path(__file__),image,
        repo/'scripts/build_spice_links.py',repo/'scripts/run_spice_campaign.py',
        repo/'scripts/audit_spice_campaign.py',repo/'scripts/auto_btf_merge_gate.py')}
    for side,core in cores.items():
        done=read(core/'completed.json');source=read(core/'source.json')
        assert done['passed'] and 'native' in done['profiles']
        if side=='candidate':assert done['profiles']==['native','portable','asan','tsan']
        assert done['source_sha256']==sha(core/'source.json')
        assert done['commands_sha256']==sha(core/'commands.json')
        assert done['libraries']['native']==sha(core/'native/libkls.a')
        pins.update(source['input_sha256'])
        for path in (core/'completed.json',core/'source.json',core/'commands.json'):
            pins[str(path)]=sha(path)
        for path in list((core/'native').glob('*.a'))+list((core/'native/lib').glob('*.a')):
            pins[str(path)]=sha(path)
    headers=lambda core:{str(p.relative_to(core/'source/include')):sha(p) for p in (core/'source/include').rglob('*.h')}
    assert headers(cores['candidate'])==headers(cores['baseline'])
    assert read(cores['candidate']/'source.json')['dependencies']==read(cores['baseline']/'source.json')['dependencies']
    out.mkdir(exist_ok=False);save(out/'inputs.json',dict(diagnostic_only=True,input_sha256=pins))
    spec=importlib.util.spec_from_file_location('cleanup',repo/'bench/process_session_cleanup.py')
    cleanup=importlib.util.module_from_spec(spec);spec.loader.exec_module(cleanup)
    env={k:v for k,v in os.environ.items() if not k.startswith(('LD_','KLS_','OMP_','GOMP_',
        'OPENBLAS_','APPTAINERENV_','ASAN_','UBSAN_','TSAN_'))}
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    prefix=['apptainer','exec','--bind',f'{workspace}:/work','--bind',f'{workspace}:{workspace}',
            '--bind','/tmp:/tmp',str(image)]
    commands=[]

    def verify():
        assert all(Path(path).exists() and sha(path)==digest for path,digest in pins.items())

    def run(stage,command,timeout=7200):
        verify();command=list(map(str,command));print('START',stage,flush=True)
        with (out/(stage+'.log')).open('x') as log:
            child=subprocess.Popen(command,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            try:code=child.wait(timeout=timeout)
            except BaseException:
                os.killpg(child.pid,signal.SIGKILL);child.wait()
                cleanup.drain_exited_apptainer_session(child.pid,image);raise
        closure=cleanup.drain_exited_apptainer_session(child.pid,image)
        commands.append(dict(stage=stage,command=command,exit_code=code,closure=closure))
        save(out/'commands.json',commands);print('END',stage,code,flush=True)
        assert code==0 and closure['all_session_processes_absent'],stage
        verify()

    for side,core in cores.items():
        dest=out/(side+'-link')
        run(side+'-link',[sys.executable,repo/'scripts/build_spice_links.py','--workspace',workspace,
            '--image',image,'--artifacts',workspace/'KLS/build/xyce-testing-bump-20261023/artifacts',
            '--provenance',workspace/'KLS/build/xyce-testing-bump-20261023/build-provenance',
            '--generated','/tmp/kls-xyce-testing-20261023','--trilinos','/tmp/kls-trilinos-install',
            '--xyce-source',workspace/'Xyce-spice-fast','--kls-source',core/'source',
            '--kls-build',core/'native','--dependency-build',core/'native',
            '--common-archive','/tmp/kls-spice-20261027/control-link-v3/libxyce.a','--output',dest])
        build=read(dest/'build.json');pins.update(build['input_sha256'])
        pins.update(build['output_sha256'])
        for name in ('build.json','commands.json'):pins[str(dest/name)]=sha(dest/name)
        save(out/'inputs.json',dict(diagnostic_only=True,input_sha256=pins))

    def campaign(label,binary,provenance):
        dest=out/label;audit=out/(label+'-independent')
        run(label,[sys.executable,repo/'scripts/run_spice_campaign.py','--workspace',workspace,
            '--regression','/tmp/kls-xyce-regression','--image',image,'--license','/tmp/kls-goal-cktso/cktso.lic',
            '--baseline',out/'baseline-link/Xyce','--candidate',binary,'--provenance',provenance,
            '--diagnostic','--passes','1','--output',dest])
        run(label+'-review',[sys.executable,repo/'scripts/run_spice_campaign.py','--review',dest])
        run(label+'-audit',[sys.executable,repo/'scripts/audit_spice_campaign.py','--campaign',dest,'--output',audit])
        r=read(dest/'review.json');q=read(audit/'audit.json')
        assert r['all_qualified'] and q['all_qualified'] and r['observations']==q['observations']==70
        key=lambda c:(c['domain'],c['case'],c['threads'],c['peer'])
        assert len(r['cells'])==len(q['cells'])==28
        assert {key(c):(c['median'],c['wall_median']) for c in r['cells']}=={key(c):(c['median'],c['wall_median']) for c in q['cells']}
        save(out/(label+'-agreement.json'),dict(all70_qualified=True,independent_exact_agreement=True))
        return dest

    campaign('api',out/'candidate-link/Xyce',out/'candidate-link/build.json')
    core=cores['candidate'];diagnostic=out/'activation-source'
    shutil.copytree(core/'source/src',diagnostic);insertions={}

    def insert(name,needle,added):
        path=diagnostic/name;before=path.read_text();assert before.count(needle)==1,(name,needle)
        path.write_text(before.replace(needle,needle+added,1));insertions.setdefault(name,[]).append(added)

    insert('kls.c','  unsigned observed_parallel_epoch;', '\n  uint64_t diagnostic_observed_accepted,diagnostic_forest_accepted,diagnostic_plan,diagnostic_stage,diagnostic_register_accepted;\n')
    insert('kls.c','#include <stdatomic.h>', '\nstatic _Atomic uint64_t diagnostic_register_calls;\n#define KLS_SNB_REGISTER_TEST_EVENT(width) atomic_fetch_add_explicit(&diagnostic_register_calls,1,memory_order_relaxed)\n#define KLS_FOREST_REMAINING_TEST_EVENT() (++s->diagnostic_plan)\n'
        '#define KLS_FOREST_STAGE_TEST_EVENT() (++s->diagnostic_stage)\n')
    insert('kls.c','  ++s->forest_combined_depth;', '\n  const uint64_t diagnostic_leaves=s->forest_leaf_numerics,diagnostic_solves=s->forest_leaf_solves,diagnostic_plan=s->diagnostic_plan,diagnostic_stage=s->diagnostic_stage;\n'
        '  const int diagnostic_observed=s->observed_parallel.choice==2;\n'
        '  const uint64_t diagnostic_register=atomic_load_explicit(&diagnostic_register_calls,memory_order_relaxed);\n')
    insert('kls.c','  --s->forest_combined_depth;', '''
  if(status==KLS_OK && atomic_load_explicit(&diagnostic_register_calls,memory_order_relaxed)>diagnostic_register)
    ++s->diagnostic_register_accepted;
  if(status==KLS_OK && diagnostic_observed && s->observed_parallel.choice==2 &&
     s->stats.last_refactor_path==KLS_REFACTOR_PATH_EGRAPH)++s->diagnostic_observed_accepted;
  if(status==KLS_OK && s->forest_unscale_attempted>0 &&
     s->forest_leaf_numerics>diagnostic_leaves && s->forest_leaf_solves>diagnostic_solves &&
     s->diagnostic_plan>diagnostic_plan && s->diagnostic_stage>diagnostic_stage)
    ++s->diagnostic_forest_accepted;
''')
    insert('kls.c','void kls_destroy(kls_solver *solver) {', '''
  if(solver)fprintf(stderr,"WORST_PEER_ACTIVATION threads=%d observed=%llu forest=%llu plan=%llu stage=%llu registers=%llu\\n",
    solver->options.threads,(unsigned long long)solver->diagnostic_observed_accepted,
    (unsigned long long)solver->diagnostic_forest_accepted,
    (unsigned long long)solver->diagnostic_plan,(unsigned long long)solver->diagnostic_stage,
    (unsigned long long)solver->diagnostic_register_accepted);
''')
    for name,adds in insertions.items():
        text=(diagnostic/name).read_text()
        for added in reversed(adds):assert text.count(added)==1;text=text.replace(added,'',1)
        assert text==(core/'source/src'/name).read_text()
    save(out/'activation-insertions.json',insertions)
    flags=[]
    for line in (core/'native/CMakeFiles/kls.dir/flags.make').read_text().splitlines():
        if line.startswith(('C_DEFINES =','C_INCLUDES =','C_FLAGS =')):flags+=shlex.split(line.split(' = ',1)[1])
    activation=out/'activation';activation.mkdir();obj=activation/'kls.c.o'
    run('activation-compile',prefix+['gcc',*flags,'-c',diagnostic/'kls.c','-o',obj])
    shutil.copyfile(core/'native/libkls.a',activation/'libkls.a')
    run('activation-archive',prefix+['ar','r',activation/'libkls.a',obj])
    run('activation-index',prefix+['ranlib',activation/'libkls.a'])
    ordinary=out/'candidate-link';build=read(ordinary/'build.json')
    command=build['commands'][-1]['command']
    assert command.count(str(core/'native/libkls.a'))==1 and command.count(str(ordinary/'Xyce'))==1
    command=[str(activation/'libkls.a') if x==str(core/'native/libkls.a') else str(activation/'Xyce') if x==str(ordinary/'Xyce') else x for x in command]
    run('activation-link',command)
    activation_pins={str(p):sha(p) for p in diagnostic.rglob('*') if p.is_file()}
    activation_pins.update({str(p):sha(p) for p in (activation/'Xyce',activation/'libkls.a',out/'activation-insertions.json')})
    save(activation/'build.json',dict(instrumented=True,input_sha256={**pins,**activation_pins}))
    pins.update(activation_pins);save(out/'inputs.json',dict(diagnostic_only=True,input_sha256=pins))
    dest=campaign('activation-api',activation/'Xyce',activation/'build.json')
    rows=[json.loads(s) for s in (dest/'observations.jsonl').read_text().splitlines()];events=[]
    for row in rows:
        if row['side']!='candidate':continue
        log=(Path(row['directory'])/'stdout.log').read_text()
        entries=[dict(zip(('threads','observed','forest','plan','stage','registers'),map(int,m)))
                 for m in re.findall(r'WORST_PEER_ACTIVATION threads=(\d+) observed=(\d+) forest=(\d+) plan=(\d+) stage=(\d+) registers=(\d+)',log)]
        assert entries, (row['case'], 'missing activation records')
        if row['case']=='NEURON/rallpack3.cir':
            assert any(e['forest']>0 and e['plan']>=e['forest'] and e['stage']>=e['forest'] for e in entries)
        events.append(dict(domain=row['domain'],case=row['case'],threads=row['threads'],owners=entries))
    assert len(events)==28
    assert any(e['observed']>0 for row in events for e in row['owners'])
    save(out/'activation-events.json',events);verify()
    save(out/'completed.json',dict(passed=True,api70_qualified=True,activation70_qualified=True,
        retained_forest_all_four_neuron_cells=True,observed_parallel_execution_seen=True,independent_exact_agreement=True,
        performance_qualified=False,source_pins_unchanged=True,
        ordinary_candidate=str(out/'candidate-link/Xyce'),ordinary_baseline=str(out/'baseline-link/Xyce'),
        input_manifest_sha256=sha(out/'inputs.json'),commands_sha256=sha(out/'commands.json')))
    print('API and retained AUTO execution qualified; performance remains separate.',flush=True)


if __name__=='__main__':main()
