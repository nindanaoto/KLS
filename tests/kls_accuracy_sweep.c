#define _POSIX_C_SOURCE 200809L
#include <time.h>
#include <stddef.h>
#include <stdint.h>
#undef NDEBUG
#include <assert.h>
#include <stdio.h>

static int clock_calls,expire,clock_failure,proposal_kind,certified;
static int fake_clock(clockid_t id,struct timespec *ts) {
  assert(id==CLOCK_MONOTONIC);
  if(clock_failure) return -1;
  ts->tv_sec=(expire && clock_calls++)?2:0;ts->tv_nsec=0;
  return 0;
}
#define clock_gettime fake_clock
#include "../src/kls_accuracy_recovery.c"
#undef clock_gettime

int kls_accuracy_lattice_propose(int m,int n,const double *a,
    const double *b,double *delta,int steps) {
  assert(m==2 && n==1 && steps>0);
  if(proposal_kind==3) return -1;
  if(!proposal_kind) return 0;
  delta[0]=proposal_kind==1?0:b[0]/a[0];
  return 1;
}
int kls_accuracy_certify(size_t n,const int64_t *p,const int64_t *rows,
    const double *a,const double *b,const double *x,int transpose,double *ratio) {
  assert(n==1 && p && rows && a && b && x && !transpose);
  /* Keep the published incumbent unchanged even if internal work improves. */
  *ratio=2;
  return certified;
}
int main(void) {
  int64_t p[]={0,1},rows[]={0};double a[]={1},b[]={1},x[]={0};
  int count=-1,stalled=-1;
  assert(!kls_accuracy_lattice_sweep(1,p,rows,a,b,x,1,1,1e-7,1,&count,&stalled));
  assert(stalled==1 && count==0 && x[0]==0);
  proposal_kind=1;stalled=-1;
  assert(!kls_accuracy_lattice_sweep(1,p,rows,a,b,x,1,1,1e-7,1,&count,&stalled));
  assert(stalled==1 && count==1 && x[0]==0);
  proposal_kind=2;stalled=-1;
  assert(!kls_accuracy_lattice_sweep(1,p,rows,a,b,x,1,1,1e-7,1,&count,&stalled));
  assert(stalled==0 && count==1 && x[0]==0);
  /* A transient proposal workspace failure is not a deterministic stop. */
  proposal_kind=3;stalled=1;
  assert(!kls_accuracy_lattice_sweep(1,p,rows,a,b,x,1,1,1e-7,1,&count,&stalled));
  assert(stalled==0 && count==0 && x[0]==0);
  /* Deadline truncation must never claim the configured sweep completed. */
  proposal_kind=0;expire=1;clock_calls=0;stalled=1;
  assert(!kls_accuracy_lattice_sweep(1,p,rows,a,b,x,1,1,1e-7,1,&count,&stalled));
  assert(stalled==0 && count==0 && x[0]==0);
  expire=0;certified=1;stalled=1;
  assert(kls_accuracy_lattice_sweep(1,p,rows,a,b,x,1,1,1e-7,1,&count,&stalled)==1);
  assert(stalled==0 && count==0);
  certified=0;clock_failure=1;stalled=1;
  assert(kls_accuracy_lattice_sweep(1,p,rows,a,b,x,1,1,1e-7,1,&count,&stalled)==-1);
  assert(stalled==0);
  clock_failure=0;stalled=1;
  assert(kls_accuracy_lattice_sweep(0,p,rows,a,b,x,1,1,1e-7,1,&count,&stalled)==-1);
  assert(stalled==0);
  puts("PASS sweep progress: complete, truncated, unpublished progress, success, failure");
  return 0;
}
