#ifndef KLS_SNB_AUTO_POLICY_H
#define KLS_SNB_AUTO_POLICY_H
#include <math.h>
#include <string.h>
#include "kls_tuning.inc"
/* One consultation: two ABBA blocks of certified production calls. */
typedef struct {
  int choice,phase,count[2]; /* 0 untried, 1 prepared, 2 sampling, 3 retained, -1 declined */
  double seconds[2][4],factor_credit,limit,spent;
} kls_snb_auto_policy;
static inline int kls_snb_auto_arm(const kls_snb_auto_policy *p) {
  static const unsigned char arms[8]={0,1,1,0,0,1,1,0};
  return p->choice==3?1:(p->phase>=0 && p->phase<8?arms[p->phase]:0);
}
static inline double kls_snb_auto_median(const double *a) {
  double v[4];memcpy(v,a,sizeof v);
  for(int i=1;i<4;++i){double x=v[i];int j=i;while(j && v[j-1]>x){v[j]=v[j-1];--j;}v[j]=x;}
  return .5*v[1]+.5*v[2];
}
static inline void kls_snb_auto_sample(kls_snb_auto_policy *p,int changed,
    int clean,double seconds) {
  if(p->choice<1 || p->choice>2 || !changed)return;
  if(!clean || !(seconds>0) || !isfinite(seconds) || p->phase<0 || p->phase>=8) {
    p->choice=-1;return;
  }
  if(!p->phase){
    p->limit=p->factor_credit+KLS_SNB_AUTO_EARLY_CALL_ALLOWANCE*seconds;
    if(!(p->limit>0) || !isfinite(p->limit) || p->spent>p->limit){p->choice=-1;return;}
  }
  const int arm=kls_snb_auto_arm(p);
  if(p->count[arm]>=4 || !isfinite(p->spent) || p->spent>p->limit){p->choice=-1;return;}
  p->seconds[arm][p->count[arm]++]=seconds;++p->phase;p->choice=2;
  if(p->phase==8){
    const double incumbent=kls_snb_auto_median(p->seconds[0]);
    const double candidate=kls_snb_auto_median(p->seconds[1]);
    p->choice=candidate<=(1.0-KLS_COMPACT_CYCLE_MIN_RELATIVE_SAVING)*incumbent &&
      incumbent-candidate>=KLS_COMPACT_CYCLE_MIN_SECONDS_SAVED?3:-1;
  }
}
#endif
