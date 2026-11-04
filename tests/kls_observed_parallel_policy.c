#include "../src/kls_observed_parallel_policy.h"
#include <assert.h>
#undef NDEBUG
#include <assert.h>

int main(void) {
  kls_observed_parallel_state g={0};
  for(int i=0;i<500;++i) {
    kls_observed_parallel_observe(&g,0,1,1);
    kls_observed_parallel_observe(&g,1,0,1);
    kls_observed_parallel_observe(&g,1,1,NAN);
  }
  assert(g.cycles==0 && g.total==0 && !kls_observed_parallel_start(&g));
  for(int i=0;i<199;++i) {
    kls_observed_parallel_observe(&g,1,1,1);
    assert(!kls_observed_parallel_start(&g));
  }
  kls_observed_parallel_observe(&g,1,1,1);
  assert(kls_observed_parallel_start(&g));
  assert(g.limit==4 && g.cycle==1);
  assert(kls_observed_parallel_preparation_budget(&g)==2);
  kls_observed_parallel_observe(&g,1,1,1000);
  assert(g.limit==4 && g.total==200 && g.cycles==200);
  g.spent=1.5;
  assert(!kls_observed_parallel_record(&g,1,0.94,1));
  assert(g.samples==1 && g.choice==1);
  assert(kls_observed_parallel_preparation_budget(&g)==1.5);
  g.spent=2.5;
  assert(kls_observed_parallel_record(&g,1,0.9,1));
  assert(g.choice==2 && g.samples==2);
  for(int reason=0;reason<6;++reason) {
    kls_observed_parallel_state bad=g;bad.choice=1;bad.samples=0;
    if(reason==0) bad.spent=bad.limit+0.01;
    assert(!kls_observed_parallel_record(&bad,reason!=1,
        reason==2?0.951:reason==3?NAN:reason==4?0:0.9,
        reason==5?INFINITY:1));
    assert(bad.choice==-1 && !kls_observed_parallel_start(&bad));
  }
  /* A timed-out joined job spends its budget permanently. Only later
     changed inputs can fund one retry, with two new scored comparisons. */
  kls_observed_parallel_state retry={.choice=1,.cycles=200,.total=200,
    .cycle=1,.limit=4,.spent=3.5,.samples=1};
  for(int i=0;i<8;++i)retry.history[i]=1;
  assert(kls_observed_parallel_defer_deadline(&retry));
  assert(retry.choice==0 && retry.samples==0 && retry.retry_pending && retry.deadline_retries==1);
  assert(retry.spent==3.5 && retry.total==200 && !kls_observed_parallel_start(&retry));
  for(int i=0;i<100;++i)kls_observed_parallel_observe(&retry,0,1,100);
  assert(retry.total==200 && !kls_observed_parallel_start(&retry));
  for(int i=0;i<74;++i) {
    kls_observed_parallel_observe(&retry,1,1,1);assert(!kls_observed_parallel_start(&retry));
  }
  kls_observed_parallel_observe(&retry,1,1,1);assert(kls_observed_parallel_start(&retry));
  assert(retry.limit==5.5 && retry.spent==3.5);
  assert(!kls_observed_parallel_defer_deadline(&retry)); /* No second retry. */
  retry.spent+=.8;assert(!kls_observed_parallel_record(&retry,1,.8,1));
  retry.spent+=.8;assert(kls_observed_parallel_record(&retry,1,.8,1));
  assert(retry.choice==2 && !retry.retry_pending && retry.spent<retry.limit);
  g=(kls_observed_parallel_state){0};
  for(int i=0;i<32;++i) kls_observed_parallel_observe(&g,1,1,i==0?10000:1);
  assert(kls_observed_parallel_start(&g) && g.cycle==1);
  return 0;
}
