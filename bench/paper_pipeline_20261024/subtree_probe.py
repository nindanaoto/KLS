from common import *
import sys,shlex
version=sys.argv[1];B=build('subtree',sys.argv[2] if len(sys.argv)>2 else version);S=source('subtree',version);W=Path('/tmp/kls-paper-subtree-probe-'+version+'-20261024');W.mkdir();runner=Runner('subtree-probe-'+version,'16-23');f=(B/'CMakeFiles/kls.dir/flags.make').read_text();flags=[]
for key in ['C_DEFINES','C_INCLUDES','C_FLAGS']:flags+=shlex.split(next(x for x in f.splitlines() if x.startswith(key+' = ')).split(' = ',1)[1])
runner.run('compile',['gcc',*flags,'-g','-DKLS_CORE_SOURCE="'+str(S/'src/kls.c')+'"','-c',Path(__file__).with_name('subtree_probe.c'),'-o',W/'probe.o'])
cmd=shlex.split((B/'CMakeFiles/kls_subtree_first.dir/link.txt').read_text());cmd=[str(W/'probe.o') if x=='CMakeFiles/kls_subtree_first.dir/tests/kls_subtree_first.c.o' else x for x in cmd if not x.startswith('-Wl,--dependency-file=')];cmd[cmd.index('-o')+1]=str(W/'probe');runner.run('link',cmd,cwd=B)
frame=R/'build/subtree-investigation-20261024/followup/ram2k/initial-frame.bin'
for matched in [0,1]:
 for crew in [2,8]:runner.run(('matched' if matched else 'original')+'-t'+str(crew),[W/'probe',frame,str(crew),str(matched)],required=False)
(O/(runner.name+'-completed.json')).write_text(json.dumps(dict(diagnostic_only=True,frame_sha256=sha(frame),processes=runner.rows),indent=2))
