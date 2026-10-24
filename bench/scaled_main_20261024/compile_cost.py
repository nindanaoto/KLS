from common import *
import shlex
W=Path('/tmp/kls-main-recovery-cost-20261024');W.mkdir(exist_ok=True);runner=Runner('cost-build','16-23')
for role in ['control','candidate']:
 B=build(role);flags=[];text=(B/'CMakeFiles/kls.dir/flags.make').read_text()
 for key in ['C_DEFINES','C_INCLUDES','C_FLAGS']:flags+=shlex.split(next(x for x in text.splitlines() if x.startswith(key+' = ')).split(' = ',1)[1])
 obj=W/(role+'.o');runner.run(role+'-compile',['gcc',*flags,'-c',Path(__file__).with_name('cost.c'),'-o',obj])
 cmd=shlex.split((B/'CMakeFiles/kls_smoke.dir/link.txt').read_text());cmd=[str(obj) if x=='CMakeFiles/kls_smoke.dir/tests/kls_smoke.c.o' else x for x in cmd if not x.startswith('-Wl,--dependency-file=')];cmd[cmd.index('-o')+1]=str(W/role);runner.run(role+'-link',cmd,cwd=B)
(O/'cost-build-completed.json').write_text(json.dumps(dict(passed=True,processes=runner.rows,binaries={role:dict(path=str(W/role),sha256=sha(W/role)) for role in ['control','candidate']}),indent=2))
