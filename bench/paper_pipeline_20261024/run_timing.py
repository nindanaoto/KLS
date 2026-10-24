import subprocess,time,json
from pathlib import Path
from common import R as r,O as o
p=Path(__file__).resolve().parent
required=[role+'-'+v+'-'+m+('-full' if m in ['native','portable'] else '')+'-completed.json' for role,v in [('recovery','v6'),('subtree','compact-v11')] for m in ['native','portable','asan','tsan']]+['recovery-v6-xyce-completed.json','recovery-v6-corpus-t1-completed.json','recovery-v6-corpus-t8-completed.json']
# Bounded orchestration wait: an absent receipt stops this campaign rather than timing an unqualified source.
deadline=time.monotonic()+3600
while not all((o/x).exists() for x in required):
 if time.monotonic()>deadline:raise RuntimeError('Qualification receipts incomplete; timing not started')
 time.sleep(10)
for name in required:
 d=json.loads((o/name).read_text());assert d.get('passed',d.get('all_qualified',d.get('all_preserved',False))),name
for campaign in ['initial','confirmation']:
 for script,stem in [('recovery_epochs.py','recovery-epochs'),('healthy.py','healthy')]:
  name=stem+'-v6-'+campaign
  with (o/(name+'-progress.log')).open('x') as f:result=subprocess.run(['python3',str(p/script),'v6',campaign],cwd=r,stdout=f,stderr=subprocess.STDOUT)
  print(name,result.returncode,flush=True)
  if result.returncode:raise RuntimeError(name)
