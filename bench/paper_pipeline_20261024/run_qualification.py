import subprocess,concurrent.futures
from pathlib import Path
from common import R as r,O as o
p=Path(__file__).resolve().parent
def series(mode):
 for role,version in [('recovery','v6'),('subtree','compact-v11')]:
  args=['python3',str(p/('tsan.py' if mode=='tsan' else 'qualify.py')),role,version]
  if mode!='tsan':args+=[mode]+(['--full'] if mode in ['native','portable'] else [])
  name=role+'-'+version+'-'+mode+('-full' if mode in ['native','portable'] else '')
  with (o/(name+'-progress.log')).open('x') as f:
   result=subprocess.run(args,cwd=r,stdout=f,stderr=subprocess.STDOUT)
  print(name,result.returncode,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as e:list(e.map(series,['native','portable','asan','tsan']))
