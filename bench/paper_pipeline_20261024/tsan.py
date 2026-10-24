from common import *
import sys,shutil
role,version=sys.argv[1:];B=build(role,version,'tsan');S=source(role,version);runner=Runner(role+'-'+version+'-tsan','16-23');W=Path('/tmp/kls-paper-'+role+'-'+version+'-tsan-controls-20261024');W.mkdir()
control=R/'build/native-btf-v8-tsan-20261021/control.c';flags=['-O1','-g','-fsanitize=thread','--param=tsan-instrument-func-entry-exit=0','-fno-omit-frame-pointer'];ENV.update(TSAN_OPTIONS='halt_on_error=1:exitcode=66:second_deadlock_stack=1')
for name,race in [('clean',False),('race',True)]:
 runner.run(name+'-compile',['gcc','-std=c11',*flags,*(['-DRACE'] if race else []),control,'-pthread','-o',W/name]);row=runner.run(name+'-run',[W/name],required=False);assert row['exit_code']==(66 if race else 0)
 if race:assert 'WARNING: ThreadSanitizer: data race' in Path(row['log']).read_text()
targets=['kls_subtree_domains','kls_subtree_first','kls_accuracy_owner','kls_certificate_pool'] if role=='subtree' else ['kls_scaled_recovery','kls_accuracy_owner','kls_accuracy_epochs','kls_certificate_pool']
targets+=['kls_klu_allocator_ownership']
runner.run('configure',['cmake','-S',S,'-B',B,'-DCMAKE_BUILD_TYPE=Release','-DCMAKE_C_FLAGS='+' '.join(flags),'-DCMAKE_C_FLAGS_RELEASE=-O1 -g -DNDEBUG','-DCMAKE_EXE_LINKER_FLAGS=-fsanitize=thread'])
runner.run('build',['cmake','--build',B,'-j8','--target',*targets])
for name in targets:runner.run(name,[B/name])
assert '__tsan_read' in subprocess.check_output(['nm','-u',str(B/'CMakeFiles/kls.dir/src/kls.c.o')],text=True)
(O/(runner.name+'-completed.json')).write_text(json.dumps(dict(passed=True,scope=targets,processes=runner.rows,controls='clean passes; positive race detected; no suppressions; GCC IFUNC entry/exit workaround only'),indent=2))
