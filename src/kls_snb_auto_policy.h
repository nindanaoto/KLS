#ifndef KLS_SNB_AUTO_POLICY_H
#define KLS_SNB_AUTO_POLICY_H
#include <math.h>
#include <string.h>
#include "kls_tuning.inc"
/* Initial ABBA consultation; at most one late ABBA recheck, without
 * increasing setup credit or performing duplicate numeric calls. */
enum {KLS_SNB_AUTO_MONITOR_SAMPLES=16};
typedef struct {
  int choice,phase,count[2]; /* 0 untried, 1 prepared, 2 sampling, 3 retained, -1 declined */
  double seconds[2][4],factor_credit,limit,spent;
  double retained[KLS_SNB_AUTO_MONITOR_SAMPLES];
  unsigned retained_count,retained_cursor;
  int rechecked;
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
    if(!p->rechecked)p->limit=p->factor_credit+KLS_SNB_AUTO_EARLY_CALL_ALLOWANCE*seconds;
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
static inline double kls_snb_auto_retained_median(const double *a) {
  double v[KLS_SNB_AUTO_MONITOR_SAMPLES];memcpy(v,a,sizeof(v));
  for(unsigned i=1;i<KLS_SNB_AUTO_MONITOR_SAMPLES;++i) {
    double x=v[i];unsigned j=i;while(j&&v[j-1]>x){v[j]=v[j-1];--j;}v[j]=x;
  }
  return .5*v[KLS_SNB_AUTO_MONITOR_SAMPLES/2-1]+.5*v[KLS_SNB_AUTO_MONITOR_SAMPLES/2];
}
/* Called only for fresh, certified, changed-operator retained calls. A sixteen-call
 * median filters short scheduling and exceptional-repair bursts. A recheck uses the existing plan
 * and the original setup budget; later slower calls cannot mint new credit. */
static inline int kls_snb_auto_recheck(kls_snb_auto_policy *p,double seconds) {
  if(p->choice!=3||p->rechecked||!(seconds>0)||!isfinite(seconds))return 0;
  p->retained[p->retained_cursor]=seconds;
  p->retained_cursor=(p->retained_cursor+1u)%KLS_SNB_AUTO_MONITOR_SAMPLES;
  if(p->retained_count<KLS_SNB_AUTO_MONITOR_SAMPLES)++p->retained_count;
  if(p->retained_count<KLS_SNB_AUTO_MONITOR_SAMPLES)return 0;
  const double incumbent=kls_snb_auto_median(p->seconds[0]);
  const double current=kls_snb_auto_retained_median(p->retained);
  if(!(incumbent>0)||!isfinite(incumbent)||
     current<=(1.+KLS_COMPACT_CYCLE_MIN_RELATIVE_SAVING)*incumbent||
     current-incumbent<KLS_COMPACT_CYCLE_MIN_SECONDS_SAVED)return 0;
  p->rechecked=1;p->choice=1;p->phase=0;p->count[0]=p->count[1]=0;
  memset(p->seconds,0,sizeof(p->seconds));return 1;
}
#endif
