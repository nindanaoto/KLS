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

/* What a crew of `crew` threads could plausibly take off a serial cycle, as
 * a share of that cycle. An ideal crew leaves cycle/crew; every dispatch
 * also pays a fixed wake-and-join cost, and real factor dependencies deliver
 * about half of the ideal. A cycle shorter than its own dispatch cost has
 * nothing to gain. */
#define KLS_OBSERVED_PARALLEL_REALIZED_SHARE 0.5
static double kls_observed_parallel_plausible_gain(double cycle,int crew,
    double dispatch) {
  if(!(cycle>0) || !isfinite(cycle) || crew<2 || !(dispatch>0) ||
      !isfinite(dispatch)) return 0;
  const double saved=cycle*(1.0-1.0/crew)-dispatch;
  return saved>0?KLS_OBSERVED_PARALLEL_REALIZED_SHARE*saved/cycle:0;
}

typedef struct {
  /* choice: 0 collecting, 1 sampling, 2 retained, -1 declined.
   * active: 1 speculative candidate, 2 mandatory incumbent restoration. */
  int choice,active,samples;
  unsigned deadline_retries;
  int retry_pending;
  /* The first candidate ran cold and lost; it was a warm-up, not a sample. */
  int warmed;
  /* One warm comparison has been lost; a second loss declines. */
  int lost;
  unsigned cycles;
  double history[8],total,cycle,limit,spent,deadline;
  /* How long a candidate may run once its crew has been dispatched. The
     deadline above is the budget's hard stop and is on the clock from the
     start of the call; this bound starts at dispatch, so the time taken to
     create the crew, which the start-up reserve already paid for, is not
     taken out of the candidate's own running time. */
  double bound;
  double candidate[2],incumbent[2];
  /* Time the caller had already spent on this matrix before the first
     observed cycle (its analysis and first factorization). It is observed
     time like the cycles are, and it grows with the same matrix, so the
     allowance counts it: a run pays for exploration out of everything it
     has spent so far, not only out of its refactorizations. */
  double prior;
  /* Largest whole extra cost of one executed comparison so far. */
  double price;
  /* What starting a crew for the first candidate is expected to cost. */
  double startup;
  /* Retained-engine calls seen in the current audit window. */
  unsigned audited;
  /* Time the retained engine has taken off the incumbent's cycle so far,
     summed over every audited call (negative when it has cost time). */
  double saved;
  /* The caller's cost model for a retained crew: its size and the fixed
     cost of one dispatch. Zero means no model, and the base allowance. */
  int crew;
  double dispatch;
  /* Share of observed time exploration may spend (see _share below). */
  double gain;
} kls_observed_parallel_state;

/* The base allowance bounds what a consultation can cost a run it does not
 * help, but it also makes a run it would help wait ten times the
 * consultation's cost on the serial engine first. When the cost model says a
 * crew could take a larger share than that off every cycle, exploration that
 * costs no more than that share of the time already observed repays itself
 * if the run merely lasts as long again, so the larger share is allowed. */
static double kls_observed_parallel_share(const kls_observed_parallel_state *g) {
  return g->gain>KLS_OBSERVED_PARALLEL_ALLOWANCE && g->gain<1?g->gain:
    KLS_OBSERVED_PARALLEL_ALLOWANCE;
}

static double kls_observed_parallel_seen(const kls_observed_parallel_state *g) {
  return g->prior>0 && isfinite(g->total+g->prior)?g->total+g->prior:g->total;
}

/* A later comparison waits until the allowance covers what a comparison has
 * actually cost on this factor (with a quarter of headroom), not only the
 * two-cycle estimate: starting one that then overruns would turn a funding
 * shortfall of a few microseconds into a permanent decline. "Later" is any
 * comparison after one has executed and set a price, which includes the
 * retry of a candidate that ran into its deadline. */
static int kls_observed_parallel_priced(const kls_observed_parallel_state *g) {
  return g->samples>0 || g->warmed || g->lost || g->price>0;
}

static int kls_observed_parallel_funded(const kls_observed_parallel_state *g,
    double budget) {
  if(!kls_observed_parallel_priced(g)) return 1;
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
  /* The median needs its eight samples; without a cost model that promises
     a gain the consultation also waits for a longer record. */
  if(g->choice || g->cycles<8) return g->choice==1;
  _Alignas(64) double v[8];memcpy(v,g->history,sizeof(v));
  for(int i=1;i<8;++i) {
    double x=v[i];int j=i;
    while(j && v[j-1]>x) {v[j]=v[j-1];--j;}v[j]=x;
  }
  const double cycle=0.5*v[3]+0.5*v[4];
  g->gain=kls_observed_parallel_plausible_gain(cycle,g->crew,g->dispatch);
  if(g->cycles<32 && !(g->gain>KLS_OBSERVED_PARALLEL_ALLOWANCE)) return 0;
  const double allowance=kls_observed_parallel_share(g)*
    kls_observed_parallel_seen(g);
  /* The first candidate also starts the crew, which costs the same on a
     ten-microsecond cycle as on a ten-millisecond one: the allowance must
     cover it too, so a run too short to repay a crew never starts one. */
  const double reserve=(g->retry_pending?2.0:4.0)*cycle+
    (kls_observed_parallel_priced(g) || g->retry_pending?0.0:g->startup);
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
  if(g->choice!=1 || !kls_observed_parallel_priced(g) || !successful ||
      !(seconds>0) || !isfinite(seconds) || !isfinite(g->total+seconds)) return;
  g->total+=seconds;
  g->limit=kls_observed_parallel_share(g)*kls_observed_parallel_seen(g);
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
       the budget, discard its time, and require two fresh wins. Among warm
       comparisons one loss is tolerated: a single sample is one solver call
       long, so interference on either arm can flip it, and the verdict is
       two wins before two losses. */
    if(g->samples==0 && !g->warmed && !g->lost) {g->warmed=1;return 0;}
    if(!g->lost) {g->lost=1;return 0;}
    g->choice=-1;return 0;
  }
  if(++g->samples==2) {g->choice=2;g->retry_pending=0;}
  return g->choice==2;
}
/* Two paired samples retain the engine, but each is one call long and taken
 * with the crew already awake. What the caller then pays is the engine at
 * the application's own call rhythm, which the samples cannot see (workers
 * that fell asleep between calls, a crew shared with other work). So every
 * retained call is observed like the incumbent calls were, and whenever the
 * median of eight of them is no better than the incumbent's median cycle the
 * retention is withdrawn, unless the engine is still ahead on balance. A
 * window can be slow because the engine does not suit the rhythm, or because
 * a crew thread lost its processor for a few milliseconds: the first shows
 * as a balance that never turned positive, the second as a dent in savings
 * already made. Withdrawal is permanent, so it waits until the engine has
 * given back everything it saved; a run can then lose at most what the
 * engine had earned it. Returns 1 when the engine must be dropped. */
static int kls_observed_parallel_audit(kls_observed_parallel_state *g,
    int successful,double seconds) {
  if(g->choice!=2 || !(g->cycle>0)) return 0;
  if(!successful || !(seconds>0) || !isfinite(seconds)) return 0;
  g->history[g->audited%8]=seconds;
  if(isfinite(g->saved+g->cycle-seconds)) g->saved+=g->cycle-seconds;
  if(++g->audited%8) return 0;
  _Alignas(64) double v[8];memcpy(v,g->history,sizeof(v));
  for(int i=1;i<8;++i) {
    double x=v[i];int j=i;
    while(j && v[j-1]>x) {v[j]=v[j-1];--j;}v[j]=x;
  }
  if(0.5*v[3]+0.5*v[4]<g->cycle || g->saved>0) return 0;
  g->choice=-1;g->retry_pending=0;return 1;
}
#endif
