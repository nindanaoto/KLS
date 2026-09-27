#define _POSIX_C_SOURCE 200809L
#include <time.h>
#include <stdint.h>
#include <stddef.h>
#undef NDEBUG
#include <assert.h>
#include <stdio.h>

static double elapsed;
static int clock_failure,calls,success_at,failure_at,expected_first=32,expected_last=64;
static int stagnation_mode;
static double published;
static int fake_clock(clockid_t id,struct timespec *ts) {
  assert(id==CLOCK_MONOTONIC);
  if(clock_failure) return -1;
  ts->tv_sec=(time_t)elapsed;ts->tv_nsec=(long)((elapsed-ts->tv_sec)*1e9);return 0;
}
#define clock_gettime fake_clock
#include "../src/kls_accuracy_schedule.c"
#undef clock_gettime

int kls_accuracy_lattice_sweep(int n,const int64_t *p,const int64_t *rows,
    const double *a,const double *b,double *x,int width,int limit,
    double regularizer,double remaining,int *proposals,int *no_progress) {
  assert(n==1 && p && rows && a && b && x);
  assert(width==(calls%2?expected_last:expected_first));
  assert(limit==1024 && regularizer==1e-7 && remaining>0);
  *no_progress=0;
  assert(x[0]==published); /* Every call receives all retained progress. */
  ++calls;elapsed+=1;
  if(calls==failure_at) { *proposals=0;return -1; }
  if(stagnation_mode==1 || (stagnation_mode==2 && width==expected_first) ||
     (stagnation_mode==4 && width==expected_last)) {
    *no_progress=1;*proposals=0;
  } else if(stagnation_mode==3) {
    /* Unchanged best but no completed/no-internal-progress proof. */
    *proposals=3;
  } else {
    x[0]=calls;published=x[0];*proposals=3;
  }
  return calls==success_at?1:0;
}
static void reset(void) {
  elapsed=0;clock_failure=calls=success_at=failure_at=0;
  stagnation_mode=0;published=0;
  expected_first=32;expected_last=64;
}
int main(void) {
  int64_t p[]={0,1},rows[]={0};double a[]={1},b[]={1},x[]={0};int count;
  /* Success after the initial pair: continuation cannot become early rejection. */
  reset();success_at=5;
  assert(kls_accuracy_lattice_schedule(1,p,rows,a,b,x,32,64,32,1024,1e-7,100,&count)==1);
  assert(calls==5 && x[0]==5 && count==15);
  /* Both widths retain their complete configured sweep allowance. */
  reset();x[0]=0;
  assert(!kls_accuracy_lattice_schedule(1,p,rows,a,b,x,32,64,32,1024,1e-7,100,&count));
  assert(calls==64 && x[0]==64 && count==192);
  /* One absolute deadline, also when another RHS/frame re-enters afterward.
   * An in-flight atomic sweep may finish; a further one must not start. */
  reset();x[0]=0;
  assert(!kls_accuracy_lattice_schedule(1,p,rows,a,b,x,32,64,32,1024,1e-7,2.5,&count));
  assert(calls==3 && x[0]==3 && count==9);
  assert(!kls_accuracy_lattice_schedule(1,p,rows,a,b,x,32,64,32,1024,1e-7,2.5,&count));
  assert(calls==3 && x[0]==3 && count==0);
  /* Failed workspace setup retains the last published best candidate. */
  reset();x[0]=0;failure_at=3;
  assert(kls_accuracy_lattice_schedule(1,p,rows,a,b,x,32,64,32,1024,1e-7,100,&count)==-1);
  assert(calls==3 && x[0]==2 && count==6);
  reset();x[0]=0;clock_failure=1;
  assert(!kls_accuracy_lattice_schedule(1,p,rows,a,b,x,32,64,32,1024,1e-7,100,&count));
  assert(!calls && x[0]==0);
  reset();x[0]=0;expected_last=32;
  assert(!kls_accuracy_lattice_schedule(1,p,rows,a,b,x,32,32,3,1024,1e-7,100,&count));
  assert(calls==3 && x[0]==3);
  /* A complete fixed-point cycle stops, including the single-width policy. */
  reset();x[0]=0;stagnation_mode=1;
  assert(!kls_accuracy_lattice_schedule(1,p,rows,a,b,x,32,64,32,1024,1e-7,100,&count));
  assert(calls==2 && x[0]==0 && count==0);
  reset();x[0]=0;stagnation_mode=1;expected_last=32;
  assert(!kls_accuracy_lattice_schedule(1,p,rows,a,b,x,32,32,32,1024,1e-7,100,&count));
  assert(calls==1 && x[0]==0);
  /* Either width can advance the incumbent; the other must be retried. */
  reset();x[0]=0;stagnation_mode=2;success_at=6;
  assert(kls_accuracy_lattice_schedule(1,p,rows,a,b,x,32,64,32,1024,1e-7,100,&count)==1);
  assert(calls==6 && x[0]==6 && count==9);
  reset();x[0]=0;stagnation_mode=4;success_at=5;
  assert(kls_accuracy_lattice_schedule(1,p,rows,a,b,x,32,64,32,1024,1e-7,100,&count)==1);
  assert(calls==5 && x[0]==5 && count==9);
  /* Equality of the published vector is not a fixed-point certificate. */
  reset();x[0]=0;stagnation_mode=3;
  assert(!kls_accuracy_lattice_schedule(1,p,rows,a,b,x,32,64,32,1024,1e-7,100,&count));
  assert(calls==64 && x[0]==0 && count==192);
  reset();x[0]=0;
  assert(kls_accuracy_lattice_schedule(1,p,rows,a,b,x,24,64,32,1024,1e-7,100,&count)==-1);
  assert(kls_accuracy_lattice_schedule(1,p,rows,a,b,x,0,64,32,1024,1e-7,100,&count)==-1);
  assert(!calls && x[0]==0);
  puts("PASS sweep scheduler: continuation, fixed points, effort, retained best, shared deadline, setup failure");
  return 0;
}
