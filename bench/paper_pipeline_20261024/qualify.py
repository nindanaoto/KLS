from common import *
import sys
role,version,mode=sys.argv[1:4];full='--full' in sys.argv[4:];S=source(role,version);B=build(role,version,mode);runner=Runner(role+'-'+version+'-'+mode+('-full' if full else ''),('16-23' if mode=='asan' else '24-31') if mode!='native' else ('0-7' if role=='subtree' else '8-15'));flags=[]
if mode=='asan':flags=['-DCMAKE_C_FLAGS=-fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie','-DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address,undefined -no-pie','-DCMAKE_C_FLAGS_RELEASE=-O1 -DNDEBUG']
if mode=='portable':flags=['-DKLS_NATIVE_ARCH=OFF','-DKLS_ENABLE_HOSTED_SIMD=OFF']
runner.run('configure',['cmake','-S',S,'-B',B,'-DCMAKE_BUILD_TYPE=Release',*flags])
targets=['kls_subtree_domains','kls_subtree_first','kls_accuracy_owner'] if role=='subtree' else ['kls_scaled_recovery','kls_componentwise_caller','kls_accuracy_epochs']
runner.run('build',['cmake','--build',B,'-j8',*([] if full else ['--target',*targets])])
runner.run('tests',['ctest','--test-dir',B,'--output-on-failure','-j1',*([] if full else ['-R','^('+'|'.join(targets)+')$'])])
(O/(runner.name+'-completed.json')).write_text(json.dumps(dict(passed=True,library=str(B/'libkls.a'),library_sha256=sha(B/'libkls.a'),processes=runner.rows),indent=2))
