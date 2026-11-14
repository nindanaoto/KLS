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
  for(int i=0;i<39;++i) {
    kls_observed_parallel_observe(&g,1,1,1);
    assert(!kls_observed_parallel_start(&g));
  }
  kls_observed_parallel_observe(&g,1,1,1);
  assert(kls_observed_parallel_start(&g));
  assert(g.limit==4 && g.cycle==1);
  assert(kls_observed_parallel_preparation_budget(&g)==2);
  kls_observed_parallel_observe(&g,1,1,1000);
  assert(g.limit==4 && g.total==40 && g.cycles==40);
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
    if(reason==2) {
      /* One cold loss is a warm-up: priced, unscored, and never repeated. */
      assert(bad.choice==1 && bad.warmed && bad.samples==0);
      /* One warm loss is tolerated; the second declines. */
      assert(!kls_observed_parallel_record(&bad,1,0.951,1));
      assert(bad.choice==1 && bad.lost && bad.samples==0);
      assert(!kls_observed_parallel_record(&bad,1,0.951,1));
    }
    assert(bad.choice==-1 && !kls_observed_parallel_start(&bad));
  }
  {
    /* After a warm-up the engine still needs two fresh wins; the verdict is
       two wins before two warm losses, in any order. */
    kls_observed_parallel_state w=g;w.choice=1;w.samples=0;w.warmed=0;
    assert(!kls_observed_parallel_record(&w,1,2.0,1) && w.choice==1 && w.warmed);
    assert(!kls_observed_parallel_record(&w,1,0.9,1) && w.samples==1);
    kls_observed_parallel_state lose=w;
    assert(!kls_observed_parallel_record(&lose,1,0.96,1) && lose.choice==1 &&
      lose.lost && lose.samples==1);
    kls_observed_parallel_state recover=lose;
    assert(kls_observed_parallel_record(&recover,1,0.9,1) && recover.choice==2);
    assert(!kls_observed_parallel_record(&lose,1,0.96,1) && lose.choice==-1);
    /* A warm loss never counts as the cold warm-up, and the reverse. */
    kls_observed_parallel_state first=g;first.choice=1;first.samples=0;
    first.warmed=0;first.lost=1;
    assert(!kls_observed_parallel_record(&first,1,2.0,1) && first.choice==-1);
    assert(kls_observed_parallel_record(&w,1,0.9,1) && w.choice==2);
    /* A warm-up never rescues a trial that overran its allowance. */
    kls_observed_parallel_state over=g;over.choice=1;over.samples=0;over.warmed=0;
    over.spent=over.limit+0.01;
    assert(!kls_observed_parallel_record(&over,1,2.0,1) && over.choice==-1);
  }
  /* A timed-out joined job spends its budget permanently. Only later
     changed inputs can fund one retry, with two new scored comparisons. */
  kls_observed_parallel_state retry={.choice=1,.cycles=40,.total=40,
    .cycle=1,.limit=4,.spent=3.5,.samples=1};
  for(int i=0;i<8;++i)retry.history[i]=1;
  assert(kls_observed_parallel_defer_deadline(&retry));
  assert(retry.choice==0 && retry.samples==0 && retry.retry_pending && retry.deadline_retries==1);
  assert(retry.spent==3.5 && retry.total==40 && !kls_observed_parallel_start(&retry));
  for(int i=0;i<100;++i)kls_observed_parallel_observe(&retry,0,1,100);
  assert(retry.total==40 && !kls_observed_parallel_start(&retry));
  for(int i=0;i<14;++i) {
    kls_observed_parallel_observe(&retry,1,1,1);assert(!kls_observed_parallel_start(&retry));
  }
  kls_observed_parallel_observe(&retry,1,1,1);assert(kls_observed_parallel_start(&retry));
  assert(retry.limit==5.5 && retry.spent==3.5);
  assert(!kls_observed_parallel_defer_deadline(&retry)); /* No second retry. */
  retry.spent+=.8;assert(!kls_observed_parallel_record(&retry,1,.8,1));
  retry.spent+=.8;assert(kls_observed_parallel_record(&retry,1,.8,1));
  assert(retry.choice==2 && !retry.retry_pending && retry.spent<retry.limit);
  {
    /* The first comparison never waits; later ones wait for the allowance
       to cover both the two-cycle estimate and the observed price. */
    kls_observed_parallel_state f={.choice=1,.cycle=1,.limit=10,.spent=9.9};
    assert(kls_observed_parallel_funded(&f,0));
    f.samples=1;f.spent=5;f.price=3;
    assert(!kls_observed_parallel_funded(&f,1));       /* estimate not covered */
    assert(kls_observed_parallel_funded(&f,1.5));      /* 5 left >= 3.75 */
    f.spent=6.5;assert(!kls_observed_parallel_funded(&f,1.5)); /* 3.5 left */
    f.samples=0;f.warmed=1;assert(!kls_observed_parallel_funded(&f,1.5));
    f.spent=6.25;assert(kls_observed_parallel_funded(&f,1.5));
    /* So does the retry of a candidate that hit its deadline: it executed
       and set a price, although it left neither a sample nor a warm-up. */
    f.warmed=0;f.spent=6.5;assert(!kls_observed_parallel_funded(&f,1.5));
    f.total=100;kls_observed_parallel_fund(&f,1,10);
    assert(f.total==110 && f.limit==11 && kls_observed_parallel_funded(&f,1.5));
    f.lost=1;f.spent=10.9;assert(!kls_observed_parallel_funded(&f,1.5));
  }
  g=(kls_observed_parallel_state){0};
  for(int i=0;i<32;++i) kls_observed_parallel_observe(&g,1,1,i==0?10000:1);
  assert(kls_observed_parallel_start(&g) && g.cycle==1);
  return 0;
}
