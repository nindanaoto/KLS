#define _POSIX_C_SOURCE 200809L
#include "kls/kls.h"
#include "../src/kls_tuning_profile.h"
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Recompute checksums so malformed-section tests exercise the parser, not
   merely the checksum guard. All writes target this test's own mkstemp. */
static int rewrite(const char *path,const char *omit,const char *extra) {
  char text[16384]={0},line[512];size_t used=0;
  FILE *f=fopen(path,"r");if(!f)return 0;
  while(fgets(line,sizeof(line),f)) {
    if(!strncmp(line,"checksum=",9)||(omit&&!strncmp(line,omit,strlen(omit))))continue;
    size_t n=strlen(line);if(used+n>=sizeof(text)){fclose(f);return 0;}
    memcpy(text+used,line,n);used+=n;
  }
  fclose(f);
  if(extra){size_t n=strlen(extra);if(used+n>=sizeof(text))return 0;memcpy(text+used,extra,n);used+=n;}
  uint64_t hash=UINT64_C(1469598103934665603);
  for(size_t i=0;i<used;++i){hash^=(unsigned char)text[i];hash*=UINT64_C(1099511628211);}
  f=fopen(path,"w");if(!f)return 0;
  int ok=fwrite(text,1,used,f)==used && fprintf(f,"checksum=%016" PRIx64 "\n",hash)>0;
  return fclose(f)==0&&ok;
}

int main(void) {
  char path[]="/tmp/kls-snb-profile-XXXXXX";int fd=mkstemp(path);
  if(fd<0)return 1;close(fd);
  kls_tuning_values values,loaded;kls_tuning_defaults(&values);
  kls_tuning_metadata meta={0};kls_snb_cost_model model={0};
  int failed=0;
  if(kls_tuning_write_host_profile(path,1,&values,0,1,NULL)||
     kls_tuning_load_host_profile(path,1,&loaded,&meta)||meta.snb_cost.enabled)failed=2;
  model.enabled=1;
  for(int i=0;i<3;++i)model.narrow[i]=(i+1)*1e-9;
  for(int i=0;i<6;++i)model.dense[i]=(i+1)*1e-9;
  for(int i=0;i<4;++i)model.panel[i]=(i+1)*1e-9;
  if(kls_tuning_write_host_profile_model(path,1,&values,&model,0,1,NULL)||
     kls_tuning_load_host_profile(path,1,&loaded,&meta)||
     !meta.snb_cost.enabled||memcmp(model.narrow,meta.snb_cost.narrow,sizeof(model.narrow))||
     memcmp(model.dense,meta.snb_cost.dense,sizeof(model.dense))||
     memcmp(model.panel,meta.snb_cost.panel,sizeof(model.panel)))failed=3;
  const char *omits[]={"snb_panel_entries=",NULL,"snb_cost_version=", "snb_narrow_edge="};
  const char *extras[]={NULL,"snb_panel_entries=1e-9\n",NULL,"snb_narrow_edge=nan\n"};
  for(int i=0;i<4;++i) {
    if(kls_tuning_write_host_profile_model(path,1,&values,&model,0,1,NULL)||
       !rewrite(path,omits[i],extras[i])||
       kls_tuning_load_host_profile(path,1,&loaded,&meta)!=KLS_ERR_TUNING_PROFILE)failed=4;
  }
  if(kls_tuning_write_host_profile_model(path,8,&values,&model,0,1,NULL)!=KLS_ERR_INVALID_ARGUMENT)failed=5;
  model.dense[0]=NAN;
  if(kls_tuning_write_host_profile_model(path,1,&values,&model,0,1,NULL)!=KLS_ERR_INVALID_ARGUMENT)failed=6;
  unlink(path);
  if(failed)fprintf(stderr,"SNB profile test failed at %d\n",failed);
  return failed;
}
