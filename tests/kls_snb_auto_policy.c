#include "../src/kls_snb_auto_policy.h"
#undef NDEBUG
#include <assert.h>
int main(void) {
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
