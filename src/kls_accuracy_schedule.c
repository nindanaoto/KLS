/* Sweep-level cold recovery. The caller owns one absolute deadline across
 * internal/caller frames and RHSs. A completed sweep publishes its best answer;
 * widening therefore preserves progress instead of spending half the budget
 * at the narrow width. No two-sweep rejection rule is imposed. */
#define _POSIX_C_SOURCE 200809L
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include <time.h>

extern int kls_accuracy_lattice_sweep(int,const int64_t *,const int64_t *,
    const double *,const double *,double *,int,int,double,double,int *,int *);

static double schedule_now(void) {
  struct timespec ts;
  if(clock_gettime(CLOCK_MONOTONIC,&ts)) return INFINITY;
  return (double)ts.tv_sec+1e-9*(double)ts.tv_nsec;
}
int kls_accuracy_lattice_schedule(int n,const int64_t *p,const int64_t *rows,
    const double *a,const double *b,double *x,int first_width,int last_width,
    int sweeps,int row_limit,double regularizer,double deadline,int *proposals_out) {
  if(proposals_out) *proposals_out=0;
  if(n<=0 || !p || !b || !x || first_width<=0 || last_width>64 ||
     first_width>last_width || last_width%first_width || sweeps<=0 ||
     row_limit<=0 || !isfinite(regularizer) || regularizer<=0 ||
     !isfinite(deadline)) return -1;
  const unsigned widths=(unsigned)(last_width/first_width);
  if(widths&(widths-1)) return -1;
  int proposals=0;
  for(int sweep=0;sweep<sweeps;++sweep) {
    int all_stalled=1;
    for(int width=first_width;width<=last_width;width*=2) {
      const double remaining=deadline-schedule_now();
      if(!(remaining>0)) return 0;
      int count=0,no_progress=0;
      int status=kls_accuracy_lattice_sweep(n,p,rows,a,b,x,width,row_limit,
                                           regularizer,remaining,&count,&no_progress);
      if(!no_progress) all_stalled=0;
      /* Each completed sweep starts from the best answer returned by its
       * predecessor. The fixed-width implementation never discards that
       * incumbent on failed setup or an inferior proposal. */
      if(count>0) proposals=count>INT_MAX-proposals?INT_MAX:proposals+count;
      if(proposals_out) *proposals_out=proposals;
      if(status!=0) return status;
    }
    /* Every width saw the same incumbent and completed without an internal
     * change. Repeating this deterministic cycle cannot advance the search.
     * This is search exhaustion, not certification or proof of infeasibility. */
    if(all_stalled) return 0;
  }
  return 0;
}
