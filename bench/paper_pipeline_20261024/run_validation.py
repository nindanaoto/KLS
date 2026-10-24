import subprocess,time,concurrent.futures
from pathlib import Path
from common import R as r,O as o
p=Path(__file__).resolve().parent
while not (o/'recovery-v6-native-full-completed.json').exists():time.sleep(10)
def run(item):
 script,args,name=item
 with (o/(name+'-progress.log')).open('x') as f:result=subprocess.run(['python3',str(p/script),*args],cwd=r,stdout=f,stderr=subprocess.STDOUT)
 print(name,result.returncode,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as e:list(e.map(run,[('corpus.py',['recovery','v6','1'],'recovery-v6-corpus-t1'),('corpus.py',['recovery','v6','8'],'recovery-v6-corpus-t8'),('applications.py',['v6'],'recovery-v6-xyce')]))
