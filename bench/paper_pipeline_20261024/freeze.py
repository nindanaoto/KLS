from common import *
import sys,shutil
role,version=sys.argv[1:];S=source(role,version);S.mkdir();W=worktree(role)
files=subprocess.check_output(['git','-C',str(W),'ls-files','--cached','--others','--exclude-standard','-z']).decode().split('\0')
inputs={}
for name in files:
 if not name:continue
 p=W/name
 if p.is_file():
  target=S/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target);inputs[name]=sha(p)
for module in ['gklib','metis','scotch','spral']:
 p=S/'third_party'/module;p.symlink_to('/work/KLS-scaled-main/third_party/'+module,target_is_directory=True)
(O/(role+'-'+version+'-source.json')).write_text(json.dumps(dict(worktree=str(W),head=subprocess.check_output(['git','-C',str(W),'rev-parse','HEAD'],text=True).strip(),source=str(S),sha256=inputs),indent=2));print('FROZEN',role,version,len(inputs),flush=True)
