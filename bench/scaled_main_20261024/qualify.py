from common import *
import sys
role,mode=sys.argv[1:];assert role in ['control','candidate'] and mode in ['native','portable','asan'];S=source(role);B=build(role,mode);run=Runner(role+'-'+mode,{'control':'8-15','candidate':'0-7'}[role] if mode=='native' else '16-23' if mode=='asan' else '24-31');flags=[]
if mode=='portable':flags=['-DKLS_NATIVE_ARCH=OFF','-DKLS_ENABLE_HOSTED_SIMD=OFF']
if mode=='asan':flags=['-DCMAKE_C_FLAGS=-fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie','-DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address,undefined -no-pie','-DCMAKE_C_FLAGS_RELEASE=-O1 -DNDEBUG']
run.run('configure',['cmake','-S',S,'-B',B,'-DCMAKE_BUILD_TYPE=Release',*flags])
selection=['kls_scaled_recovery','kls_componentwise_caller','kls_accuracy_epochs'] if mode=='asan' else []
run.run('build',['cmake','--build',B,'-j8',*(['--target',*selection] if selection else [])])
run.run('tests',['ctest','--test-dir',B,'--output-on-failure','-j1',*(['-R','^('+ '|'.join(selection)+')$'] if selection else [])])
(O/(run.name+'-completed.json')).write_text(json.dumps(dict(passed=True,source=str(S),library=str(B/'libkls.a'),library_sha256=sha(B/'libkls.a'),processes=run.rows),indent=2))
