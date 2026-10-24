from pathlib import Path
import subprocess,os,signal,json,importlib.util,hashlib,resource
WORK=Path(__file__).resolve().parents[3];R=WORK/'KLS';O=R/'build/paper-pipeline-20261024';IMAGE=R/'build/kls-xyce-private.sif'
def worktree(role):
 role='subtree' if role=='subtree' else 'recovery';durable=WORK/('KLS-paper-'+role)
 return durable if durable.exists() else Path('/tmp/kls-paper-'+role+'-work-20261024')
def build(role,version='v1',mode='native'):return Path('/tmp/kls-paper-'+role+'-'+version+'-'+mode+'-20261024')
def source(role,version='v1'):return Path('/tmp/kls-paper-'+role+'-'+version+'-source-20261024')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
spec=importlib.util.spec_from_file_location('cleanup',R/'bench/process_session_cleanup.py');cleanup=importlib.util.module_from_spec(spec);spec.loader.exec_module(cleanup)
ENV={k:v for k,v in os.environ.items() if not k.startswith(('LD_','KLS_','OMP_','GOMP_','OPENBLAS_','APPTAINERENV_','ASAN_','UBSAN_','TSAN_'))};ENV.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',ASAN_OPTIONS='detect_leaks=1:halt_on_error=1',UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1',XYCE_NO_TRACKING='1')
class Runner:
 def __init__(self,name,cpus):
  self.name=name;self.cpus=cpus;self.rows=[]
  lo,hi=map(int,cpus.split('-'));assert set(range(lo,hi+1))<=os.sched_getaffinity(0),cpus
 def run(self,stage,cmd,cwd=None,timeout=1800,required=True,memory_gb=None):
  args=['taskset','-c',self.cpus,'apptainer','exec','--bind',str(WORK)+':/work','--bind','/tmp:/tmp']
  if cwd:args+=['--pwd',str(cwd).replace(str(WORK)+'/', '/work/')]
  args+=[str(IMAGE),*[str(x).replace(str(WORK)+'/', '/work/') for x in cmd]]
  def cap():resource.setrlimit(resource.RLIMIT_AS,(memory_gb*1024**3,)*2)
  log=O/(self.name+'-'+stage+'.log')
  with log.open('x') as f:
   p=subprocess.Popen(args,env=ENV,stdout=f,stderr=subprocess.STDOUT,start_new_session=True,preexec_fn=cap if memory_gb else None)
   timed=False
   try:rc=p.wait(timeout=timeout)
   except subprocess.TimeoutExpired:timed=True;os.killpg(p.pid,signal.SIGKILL);rc=p.wait()
   except BaseException:os.killpg(p.pid,signal.SIGKILL);p.wait();cleanup.drain_exited_apptainer_session(p.pid,IMAGE);raise
  row=dict(stage=stage,command=args,pid=p.pid,exit_code=rc,timed_out=timed,closure=cleanup.drain_exited_apptainer_session(p.pid,IMAGE),log=str(log));self.rows.append(row);(O/(self.name+'-processes.json')).write_text(json.dumps(self.rows,indent=2));print(self.name,stage,rc,flush=True)
  assert row['closure']['all_session_processes_absent']
  if required:assert rc==0 and not timed,row
  return row
