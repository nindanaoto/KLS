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
  /* The longest recent gap sets the budget. */
  feed(&g,4000000);
  assert(atomic_load(&g.budget_ns)>=6000000);
  /* One gap beyond the cap withdraws extended polling until it ages out. */
  feed(&g,KLS_IDLE_GAP_SPIN_CAP_NS);
  assert(!atomic_load(&g.budget_ns));
  start=0;tick=0;
  assert(!kls_idle_gap_poll(&g,&start,&tick));
  for(int i=0;i<KLS_IDLE_GAP_HISTORY-1;++i) {
    feed(&g,1000000);
    assert(!atomic_load(&g.budget_ns));
  }
  feed(&g,1000000);
  assert(atomic_load(&g.budget_ns)>=1500000);
  return 0;
}
