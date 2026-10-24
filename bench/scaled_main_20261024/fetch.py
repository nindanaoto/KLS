from pathlib import Path
import concurrent.futures,urllib.request,tarfile,io,json,re,hashlib
R=Path('/home/work-bleaker/kmatsuoka/sources/kls-workspace/KLS');O=R/'build/scaled-main-20261024';D=Path('/tmp/kls-main-recovery-corpus-20261024');D.mkdir();rows=[]
manifest=R/'bench/suitesparse_generalization_dev_manifest.txt'
def fetch(line):
 name=line.split()[0];group=re.search(r'group=(\S+)',line)[1];url=f'https://sparse.tamu.edu/MM/{group}/{name}.tar.gz';row=dict(name=name,group=group,url=url)
 try:
  with urllib.request.urlopen(url,timeout=30) as response:payload=response.read()
  with tarfile.open(fileobj=io.BytesIO(payload),mode='r:gz') as archive:
   content=archive.extractfile(name+'/'+name+'.mtx').read()
  p=D/(name+'.mtx');p.write_bytes(content);row.update(path=str(p),sha256=hashlib.sha256(content).hexdigest(),status='downloaded')
 except Exception as e:row.update(status='unavailable',error=str(e))
 print(name,row['status'],flush=True);return row
lines=[x for x in manifest.read_text().splitlines() if x and not x.startswith('#')]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:rows=list(pool.map(fetch,lines))
(O/'generalization-downloads.json').write_text(json.dumps(dict(manifest=str(manifest),manifest_sha256=hashlib.sha256(manifest.read_bytes()).hexdigest(),rows=rows),indent=2))
