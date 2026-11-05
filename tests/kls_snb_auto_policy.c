#include "../src/kls_snb_auto_policy.h"
#undef NDEBUG
#include <assert.h>
static kls_snb_auto_policy retained_policy(void) {
 kls_snb_auto_policy p={.choice=3,.phase=8,.count={4,4},
   .factor_credit=1000.,.limit=.02,.spent=.001};
 for(int i=0;i<4;++i){p.seconds[0][i]=.002;p.seconds[1][i]=.001;}
 return p;
}
static void late_recheck_checks(void) {
 kls_snb_auto_policy p=retained_policy();
 assert(!kls_snb_auto_recheck(&p,NAN)&&!p.retained_count);
 for(int i=0;i<KLS_SNB_AUTO_MONITOR_SAMPLES-1;++i)assert(!kls_snb_auto_recheck(&p,.001));
 assert(!kls_snb_auto_recheck(&p,.1)&&p.choice==3); /* One outlier. */
 p=retained_policy();for(int i=0;i<KLS_SNB_AUTO_MONITOR_SAMPLES-1;++i)assert(!kls_snb_auto_recheck(&p,.004));
 assert(kls_snb_auto_recheck(&p,.004)&&p.choice==1&&p.rechecked&&p.phase==0);
 double limit=p.limit,spent=p.spent;
 for(int i=0;i<8;++i) {
  int arm=kls_snb_auto_arm(&p);kls_snb_auto_sample(&p,1,1,arm?.004:.002);
  assert(p.limit==limit&&p.spent==spent);
 }
 assert(p.choice==-1);
 p=retained_policy();for(int i=0;i<KLS_SNB_AUTO_MONITOR_SAMPLES;++i)(void)kls_snb_auto_recheck(&p,.004);
 for(int i=0;i<8;++i){int arm=kls_snb_auto_arm(&p);kls_snb_auto_sample(&p,1,1,arm?.002:.004);}
 assert(p.choice==3&&p.rechecked&&p.limit==limit);
 for(int i=0;i<100;++i)assert(!kls_snb_auto_recheck(&p,1.)&&p.choice==3);
 p=retained_policy();for(int i=0;i<KLS_SNB_AUTO_MONITOR_SAMPLES;++i)(void)kls_snb_auto_recheck(&p,.004);
 p.spent=p.limit+.001;kls_snb_auto_sample(&p,1,1,.004);
 assert(p.choice==-1&&p.limit==limit);
}
static void column_confidence_checks(void) {
 for(int column=0;column<2;++column)for(int strong=0;strong<2;++strong) {
  kls_snb_auto_policy p={.choice=1,.factor_credit=.01,.spent=.001,.column_only=column};
  for(int i=0;i<8;++i) {
   const int arm=kls_snb_auto_arm(&p);
   kls_snb_auto_sample(&p,1,1,arm?(strong?.0075:.0085):.01);
  }
  assert(p.choice==((!column||strong)?3:-1));
  assert(p.phase==8&&p.count[0]==4&&p.count[1]==4);
  assert(fabs(p.limit-.09)<1e-15&&p.spent==.001);
 }
}
int main(void) {
  late_recheck_checks();
  column_confidence_checks();
  const int arms[8]={0,1,1,0,0,1,1,0};
  for(int verdict=0;verdict<4;++verdict) {
    kls_snb_auto_policy p={.choice=1,.factor_credit=.01,.spent=.001};
    kls_snb_auto_sample(&p,0,1,.01);assert(!p.phase);
    for(int i=0;i<8;++i) {
      assert(kls_snb_auto_arm(&p)==arms[i]);
      double t=arms[i]?(verdict==1?.0099:verdict==2?NAN:.009):.01;
      kls_snb_auto_sample(&p,1,verdict!=3,t);
      if(p.choice<0)break;
    }
    if(!verdict)assert(p.choice==3 && p.count[0]==4 && p.count[1]==4);
    else assert(p.choice==-1);
  }
  kls_snb_auto_policy p={.choice=1,.factor_credit=.01,.spent=.1};
  kls_snb_auto_sample(&p,1,1,.01);assert(p.choice==-1);
  p=(kls_snb_auto_policy){.choice=1,.factor_credit=.01};
  kls_snb_auto_sample(&p,1,1,.01);assert(p.phase==1 && fabs(p.limit-.09)<1e-15);
  p.spent=p.limit+.001;kls_snb_auto_sample(&p,1,1,.009);assert(p.choice==-1);
  return 0;
}
