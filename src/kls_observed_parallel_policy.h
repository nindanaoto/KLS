/* Private bounded consultation; one refunded retry after worker deadline. */
#ifndef KLS_OBSERVED_PARALLEL_POLICY_H
#define KLS_OBSERVED_PARALLEL_POLICY_H
#include <math.h>
#include <limits.h>
#include <string.h>

/* Share of the observed changed-call time that exploration may spend. The
 * consultation costs a fixed handful of cycles whenever it runs, so its
 * relative weight is that handful over the calls seen so far: a tenth lets it
 * start after a few dozen calls, and the same absolute cost keeps shrinking
 * as a share of every longer run. */
#define KLS_OBSERVED_PARALLEL_ALLOWANCE 0.10

typedef struct {
  /* choice: 0 collecting, 1 sampling, 2 retained, -1 declined.
   * active: 1 speculative candidate, 2 mandatory incumbent restoration. */
  int choice,active,samples;
  unsigned deadline_retries;
  int retry_pending;
  /* The first candidate ran cold and lost; it was a warm-up, not a sample. */
  int warmed;
  unsigned cycles;
  double history[8],total,cycle,limit,spent,deadline;
  double candidate[2],incumbent[2];
  /* Largest whole extra cost of one executed comparison so far. */
  double price;
} kls_observed_parallel_state;

/* A later comparison waits until the allowance covers what a comparison has
 * actually cost on this factor (with a quarter of headroom), not only the
 * two-cycle estimate: starting one that then overruns would turn a funding
 * shortfall of a few microseconds into a permanent decline. */
static int kls_observed_parallel_funded(const kls_observed_parallel_state *g,
    double budget) {
  if(g->samples<=0 && !g->warmed) return 1;
  return budget>g->cycle && g->limit-g->spent>=1.25*g->price;
}

static void kls_observed_parallel_observe(kls_observed_parallel_state *g,
    int changed,int successful,double seconds) {
  if(g->choice || !changed || !successful || !(seconds>0) ||
      !isfinite(seconds)) return;
  g->history[g->cycles%8]=seconds;
  if(g->cycles<UINT_MAX) ++g->cycles;
  if(isfinite(g->total+seconds)) g->total+=seconds;
}

static int kls_observed_parallel_start(kls_observed_parallel_state *g) {
  if(g->choice || g->cycles<32) return g->choice==1;
  _Alignas(64) double v[8];memcpy(v,g->history,sizeof(v));
  for(int i=1;i<8;++i) {
    double x=v[i];int j=i;
    while(j && v[j-1]>x) {v[j]=v[j-1];--j;}v[j]=x;
  }
  const double cycle=0.5*v[3]+0.5*v[4],
    allowance=KLS_OBSERVED_PARALLEL_ALLOWANCE*g->total;
  const double reserve=(g->retry_pending?2.0:4.0)*cycle;
  if(!(cycle>0) || !isfinite(cycle) || !isfinite(allowance) ||
      !isfinite(g->spent) || g->spent<0 || allowance-g->spent<reserve) return 0;
  g->cycle=cycle;g->limit=allowance;g->choice=1;
  return 1;
}

/* A first paired sample was taken but the allowance cannot cover the second.
 * Declining here would let a one-time validation cost decide the engine, so
 * later changed incumbent calls extend the same allowance instead. The
 * confirming sample is taken only once the allowance again covers it. The
 * samples after a cold warm-up wait for funding in the same way. */
static void kls_observed_parallel_fund(kls_observed_parallel_state *g,
    int successful,double seconds) {
  if(g->choice!=1 || (g->samples<=0 && !g->warmed) || !successful ||
      !(seconds>0) || !isfinite(seconds) || !isfinite(g->total+seconds)) return;
  g->total+=seconds;g->limit=KLS_OBSERVED_PARALLEL_ALLOWANCE*g->total;
}

static double kls_observed_parallel_preparation_budget(
    const kls_observed_parallel_state *g) {
  if(g->choice!=1 || g->samples<0 || g->samples>=2) return 0;
  const double left=g->limit-g->spent-(2-g->samples)*g->cycle;
  return left>0 && isfinite(left)?left:0;
}

/* Only the caller can identify a joined worker deadline. Preserve its full
 * cost and collect later changed-input funding before two fresh comparisons.
 * A numeric/certification failure or measured loss never requests this retry. */
static int kls_observed_parallel_defer_deadline(kls_observed_parallel_state *g) {
  if(g->choice!=1 || g->deadline_retries || !isfinite(g->spent) || g->spent<0)return 0;
  g->deadline_retries=1;g->retry_pending=1;g->choice=0;g->samples=0;
  memset(g->candidate,0,sizeof g->candidate);memset(g->incumbent,0,sizeof g->incumbent);
  return 1;
}

static int kls_observed_parallel_record(kls_observed_parallel_state *g,
    int successful,double candidate,double incumbent) {
  if(g->choice!=1) return 0;
  if(successful && g->samples>=0 && g->samples<2 && candidate>0 &&
      incumbent>0 && isfinite(candidate) && isfinite(incumbent)) {
    g->candidate[g->samples]=candidate;g->incumbent[g->samples]=incumbent;
  }
  if(!successful || g->samples<0 || g->samples>=2 ||
      !(candidate>0) || !(incumbent>0) || !isfinite(candidate) ||
      !isfinite(incumbent) || !isfinite(g->spent) ||
      g->spent>g->limit) {
    g->choice=-1;return 0;
  }
  if(candidate>0.95*incumbent) {
    /* The first candidate execution creates its crew and faults in its
       representation, costs a retained engine never pays again. Losing that
       one comparison proves nothing about the warm engine: keep its price in
       the budget, discard its time, and require two fresh wins. A second
       loss, or any loss after a win, declines. */
    if(g->samples==0 && !g->warmed) {g->warmed=1;return 0;}
    g->choice=-1;return 0;
  }
  if(++g->samples==2) {g->choice=2;g->retry_pending=0;}
  return g->choice==2;
}
#endif
