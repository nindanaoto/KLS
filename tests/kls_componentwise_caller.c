#define _GNU_SOURCE
#define _POSIX_C_SOURCE 200809L
#include <stdlib.h>
static int fail_malloc,fail_calloc;
static void *caller_malloc(size_t n) {
  if(fail_malloc) {fail_malloc=0;return NULL;}return malloc(n);
}
static void *caller_calloc(size_t n,size_t size) {
  if(fail_calloc) {fail_calloc=0;return NULL;}return calloc(n,size);
}
#define malloc caller_malloc
#define calloc caller_calloc
#define kls_componentwise_certify caller_certificate
#define kls_accuracy_certify caller_strict_certificate
#define kls_accuracy_lattice_recover caller_forbidden_lattice
#define kls_accuracy_lattice_schedule caller_forbidden_schedule
#include "../src/kls.c"
#undef malloc
#undef calloc
#undef kls_componentwise_certify
#undef kls_accuracy_certify
#undef kls_accuracy_lattice_recover
#undef kls_accuracy_lattice_schedule
#undef NDEBUG
#include <assert.h>
extern int kls_componentwise_certify(size_t,const int64_t*,const int64_t*,const double*,const double*,const double*,int,double*,long double*);
extern int kls_accuracy_certify(size_t,const int64_t*,const int64_t*,const double*,const double*,const double*,int,double*);
static int calls,residual_calls,strict_calls,reject_all,oom_call;
/* Trap both external cold-search entrances, including recovery after forced
 * exhaustion. No production instrumentation or acceptance shortcut is added. */
int caller_forbidden_lattice(int n,const int64_t *p,const int64_t *r,const double *a,
    const double *b,double *x,int width,int sweeps,int limit,double reg,double deadline,int *count) {
  (void)n;(void)p;(void)r;(void)a;(void)b;(void)x;(void)width;(void)sweeps;
  (void)limit;(void)reg;(void)deadline;(void)count;abort();
}
int caller_forbidden_schedule(int n,const int64_t *p,const int64_t *r,const double *a,
    const double *b,double *x,int first,int last,int sweeps,int limit,double reg,double deadline,int *count) {
  (void)n;(void)p;(void)r;(void)a;(void)b;(void)x;(void)first;(void)last;
  (void)sweeps;(void)limit;(void)reg;(void)deadline;(void)count;abort();
}
int caller_certificate(size_t n,const int64_t *p,const int64_t *r,const double *a,
    const double *b,const double *x,int transpose,double *res,long double *upper) {
  ++calls;residual_calls+=res!=NULL;
  if(calls==oom_call) return -1;
  if(reject_all) {
    *upper=1;
    if(res) memset(res,0,n*sizeof(*res));
    return 0;
  }
  return kls_componentwise_certify(n,p,r,a,b,x,transpose,res,upper);
}
int caller_strict_certificate(size_t n,const int64_t *p,const int64_t *r,const double *a,
    const double *b,const double *x,int transpose,double *upper) {
  ++strict_calls;return kls_accuracy_certify(n,p,r,a,b,x,transpose,upper);
}
int main(void) {
  int64_t p[]={0,1,2},r[]={0,1};double a[]={2,4};
  for(int orientation=0;orientation<3;++orientation) for(int api=0;api<3;++api) {
    kls_solver *s=NULL;kls_options o;kls_default_options(&o);o.orientation=orientation;o.threads=1;
    assert(kls_create(&s)==KLS_OK);
    assert(kls_set_accuracy_policy(s,KLS_ACCURACY_COMPONENTWISE_BACKWARD_ERROR)==KLS_OK);
    assert(kls_analyze_csc(s,KLS_INDEX_INT64,2,p,r,0,&o)==KLS_OK);
    assert(kls_factor(s,a)==KLS_OK);
    for(int fault=0;fault<6;++fault) {
      double x[]={2,4,77,4,8,88},before[6];memcpy(before,x,sizeof x);
      calls=residual_calls=strict_calls=0;
      fail_malloc=fault==0;fail_calloc=fault==1;
      oom_call=fault==2?1:(fault==3?2:0);reject_all=fault==4;
      int status=api==2?kls_refactor_solve(s,a,2,x,3,x,3):
          api?kls_solve_transpose(s,2,x,3,x,3):kls_solve(s,2,x,3,x,3);
      assert(strict_calls==0 && s->private_cold_deadline==0);
      assert(!fail_malloc && !fail_calloc);
      if(fault<4) assert(status==KLS_ERR_OUT_OF_MEMORY);
      else if(fault==4) {
        assert(status==KLS_ERR_SOLVE_FAILED);
        assert(residual_calls==KLS_COMPONENTWISE_MAX_CORRECTIONS);
        assert(calls==1+2*KLS_COMPONENTWISE_MAX_CORRECTIONS);
      } else {
        assert(status==KLS_OK && calls==2 && residual_calls==0);
        assert(x[0]==1 && x[1]==1 && x[2]==77 && x[3]==2 && x[4]==2 && x[5]==88);
      }
      if(fault<5) assert(!memcmp(x,before,sizeof x));
    }
    calls=0;double b[]={2,4},x[]={1,1};
    assert(private_componentwise_finish(s,0,1,b,x,0,KLS_ERR_SOLVE_FAILED)==KLS_ERR_SOLVE_FAILED);
    assert(calls==0); /* A failed kernel must never enter acceptance. */
    kls_destroy(s);
  }
  puts("PASS componentwise caller: atomic faults, second RHS, correction budget, no strict gate");
  return 0;
}
