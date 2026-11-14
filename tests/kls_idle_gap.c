/* Gap-predicted worker polling: extended polling is granted only after a
   full history of short gaps, is bounded by that history, and is withdrawn
   by one long gap. */
#include "../src/kls.c"
#undef NDEBUG
#include <assert.h>

static void feed(kls_idle_gap *g,uint64_t gap_ns) {
  g->last_finish_ns=kls_idle_gap_now()-gap_ns;
  kls_idle_gap_dispatch(g);
}
int main(void) {
  kls_idle_gap g={0};
  uint64_t start=0;unsigned tick=0;
  kls_idle_gap_dispatch(&g);                    /* no finish seen yet */
  assert(!atomic_load(&g.budget_ns) && !g.count);
  assert(!kls_idle_gap_poll(&g,&start,&tick));
  for(int i=0;i<KLS_IDLE_GAP_HISTORY-1;++i) {
    feed(&g,1000000);
    assert(!atomic_load(&g.budget_ns));         /* history not yet full */
  }
  feed(&g,1000000);
  const uint64_t budget=atomic_load(&g.budget_ns);
  assert(budget>=1500000 && budget<=3000000);
  /* A worker polls for about the budget and then gives up. */
  const uint64_t begin=kls_idle_gap_now();
  unsigned long polls=0;
  while(kls_idle_gap_poll(&g,&start,&tick))++polls;
  const uint64_t waited=kls_idle_gap_now()-begin;
  assert(polls>0 && waited>=budget/2 && waited<=budget+50000000);
  /* A longest gap that fits under the cap sets the budget. */
  feed(&g,4000000);
  assert(atomic_load(&g.budget_ns)>=6000000);
  for(int i=0;i<KLS_IDLE_GAP_HISTORY;++i) feed(&g,1000000);
  /* Otherwise the upper-quartile gap does: up to a quarter of the recent
     gaps may be long pauses without changing it. */
  for(int i=0;i<KLS_IDLE_GAP_HISTORY/4;++i) {
    feed(&g,10*KLS_IDLE_GAP_SPIN_CAP_NS);
    const uint64_t b=atomic_load(&g.budget_ns);
    assert(b>=1500000 && b<=3000000);
  }
  /* One more long pause makes long gaps typical: polling is withdrawn. */
  feed(&g,10*KLS_IDLE_GAP_SPIN_CAP_NS);
  assert(!atomic_load(&g.budget_ns));
  start=0;tick=0;
  assert(!kls_idle_gap_poll(&g,&start,&tick));
  /* It returns once a long pause has aged out of the history, leaving short
     gaps as three quarters of it again. */
  int restored=0;
  for(int i=0;i<KLS_IDLE_GAP_HISTORY;++i) {
    feed(&g,1000000);
    if(atomic_load(&g.budget_ns)){restored=i+1;break;}
  }
  assert(restored==KLS_IDLE_GAP_HISTORY-KLS_IDLE_GAP_HISTORY/4);
  /* The typical gap, not the shortest, sets the size. */
  for(int i=0;i<KLS_IDLE_GAP_HISTORY;++i) feed(&g,4000000);
  assert(atomic_load(&g.budget_ns)>=6000000);
  /* Typical gaps beyond the cap never get extended polling. */
  for(int i=0;i<KLS_IDLE_GAP_HISTORY;++i) feed(&g,KLS_IDLE_GAP_SPIN_CAP_NS);
  assert(!atomic_load(&g.budget_ns));
  return 0;
}
