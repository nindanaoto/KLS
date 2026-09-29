/* Test-only allocator/worker-environment injection and plan lifecycle checks. */
#include <stdlib.h>
#include <stdio.h>
#include <stdatomic.h>
#include <fenv.h>
#include <pthread.h>
#ifdef _OPENMP
#include <omp.h>
#endif
static atomic_int allocations,live;
static int fail_at,bad_worker;
static void *fault_malloc(size_t size) {
  if(atomic_fetch_add(&allocations,1)+1==fail_at) return NULL;
  void *p=malloc(size);if(p) atomic_fetch_add(&live,1);return p;
}
static void *fault_calloc(size_t n,size_t size) {
  if(atomic_fetch_add(&allocations,1)+1==fail_at) return NULL;
  void *p=calloc(n,size);if(p) atomic_fetch_add(&live,1);return p;
}
static void fault_free(void *p) {if(p) atomic_fetch_sub(&live,1);free(p);}
static int fault_round(void) {
#ifdef _OPENMP
  if(bad_worker && omp_in_parallel() && omp_get_thread_num()==1) return FE_UPWARD;
#endif
  return fegetround();
}
#define malloc fault_malloc
#define calloc fault_calloc
#define free fault_free
#define fegetround fault_round
#include "../src/kls_componentwise.c"
#undef malloc
#undef calloc
#undef free
#undef fegetround
#define CHECK(c) do {if(!(c)) {fprintf(stderr,"line %d: %s\n",__LINE__,#c);return 1;}} while(0)

typedef struct {kls_componentwise_plan *plan;int failed;} context;
typedef struct {int calls,mode;} executor_context;
static int test_executor(void *opaque,void (*work)(void *,int),void *job,int threads) {
  executor_context *c=opaque;++c->calls;
  if(c->mode==1) return 0;
  for(int tid=threads-1;tid>=0;--tid) {
    if(c->mode==2 && tid==1) fesetround(FE_UPWARD);
    work(job,tid);
    if(c->mode==2 && tid==1) fesetround(FE_TONEAREST);
  }
  return 1;
}
static void *independent(void *opaque) {
  context *c=opaque;
  int64_t p[]={0,2,4},r[]={0,1,0,1};double a[]={2,3,1,4},x[]={1,2};
  for(int i=0;i<30;++i) {
    double b[]={4,11},bt[]={8,9};long double upper;
    if(kls_componentwise_certify_planned(&c->plan,2,2,p,r,a,i%2?bt:b,x,i%2,&upper)!=1)
      {c->failed=1;break;}
  }
  return NULL;
}

int main(void) {
#if KLS_COMPONENTWISE_HAS_ACCUMULATED && defined(_OPENMP)
  omp_set_dynamic(0);
  int64_t p[]={0,2,4},r[]={0,1,0,1};double a[]={2,3,1,4},x[]={1,2},b[]={4,11},bt[]={8,9};
  long double upper;
  for(fail_at=1;fail_at<=5;++fail_at) {
    allocations=0;kls_componentwise_plan *plan=NULL;
    CHECK(kls_componentwise_certify_planned(&plan,2,2,p,r,a,b,x,0,&upper)==1);
    kls_componentwise_plan_destroy(plan);CHECK(live==0);
  }
  fail_at=0;
  kls_componentwise_plan *external=NULL;
  executor_context executor={0};
  for(int d=0;d<2;++d) for(int i=0;i<5;++i)
    CHECK(kls_componentwise_certify_executor(&external,4,2,p,r,a,d?bt:b,x,d,&upper,
              test_executor,&executor)==1);
  CHECK(executor.calls>=2*KLS_CERTIFICATE_PROBE_PAIRS);
  for(int mode=1;mode<=2;++mode) {
    external->direction[0].choice=1;executor.mode=mode;
    CHECK(kls_componentwise_certify_executor(&external,4,2,p,r,a,b,x,0,&upper,
              test_executor,&executor)==1); /* unchanged safe fallback */
    CHECK(external->eligible==-1 && fegetround()==FE_TONEAREST);
    kls_componentwise_plan_reset(external);
  }
  kls_componentwise_plan_destroy(external);CHECK(live==0);
  /* Serial row requests must not launch even trial workers. The injected
   * unsupported worker would disable a plan if any parallel trial ran. */
  kls_componentwise_plan *serial_only=NULL;
  bad_worker=1;
  for(int d=0;d<2;++d) for(int i=0;i<5;++i)
    CHECK(kls_componentwise_certify_planned(&serial_only,1,2,p,r,a,d?bt:b,x,d,&upper)==1);
  CHECK(serial_only && serial_only->eligible==1 && serial_only->threads==0);
  CHECK(!serial_only->direction[0].samples && !serial_only->direction[1].samples);
  bad_worker=0;kls_componentwise_plan_destroy(serial_only);CHECK(live==0);
  for(int threads=2;threads<=8;threads*=2) {
    kls_componentwise_plan *plan=NULL;
    for(int d=0;d<2;++d) for(int i=0;i<KLS_CERTIFICATE_PROBE_PAIRS+2;++i)
      CHECK(kls_componentwise_certify_planned(&plan,threads,2,p,r,a,d?bt:b,x,d,&upper)==1);
    CHECK(plan && plan->eligible==1 && plan->threads==threads);
    for(int d=0;d<2;++d) CHECK(plan->direction[d].choice && plan->direction[d].samples==KLS_CERTIFICATE_PROBE_PAIRS);
    size_t *rows=plan->rowptr;
    for(int d=0;d<2;++d) {
      long double serial,parallel;
      CHECK(componentwise_plan_run(plan,1,p,r,a,d?bt:b,x,d,&serial)==1);
      CHECK(componentwise_plan_run(plan,threads,p,r,a,d?bt:b,x,d,&parallel)==1);
      CHECK(serial==parallel);
    }
    /* Changing values must not reuse an old residual or restart selection. */
    double updated[]={4,6,2,8},updated_b[]={8,22};
    CHECK(kls_componentwise_certify_planned(&plan,threads,2,p,r,updated,updated_b,x,0,&upper)==1);
    CHECK(plan->rowptr==rows && plan->direction[0].samples==KLS_CERTIFICATE_PROBE_PAIRS);
    plan->direction[0].choice=1;bad_worker=1;
    CHECK(kls_componentwise_certify_planned(&plan,threads,2,p,r,a,b,x,0,&upper)==1);
    CHECK(plan->eligible==-1);bad_worker=0;
    kls_componentwise_plan_reset(plan);CHECK(plan->eligible==0 && !plan->rowptr);
    int64_t changed[]={1,0,0,1};double reordered[]={3,2,1,4};
    CHECK(kls_componentwise_certify_planned(&plan,threads,2,p,changed,reordered,b,x,0,&upper)==1);
    CHECK(plan->eligible==1);
    /* A configured-thread change redoes only the bounded timing selection. */
    CHECK(kls_componentwise_certify_planned(&plan,threads==2?4:2,2,p,changed,reordered,b,x,0,&upper)==1);
    CHECK(plan->direction[0].samples==1 && plan->direction[1].samples==0);
    kls_componentwise_plan_destroy(plan);CHECK(live==0);
  }
  kls_componentwise_plan *duplicate=NULL;
  int64_t dp[]={0,3},dr[]={0,0,0};double da[]={0x1p54,1,-0x1p54},db[]={1},dx[]={2};
  CHECK(kls_componentwise_certify_planned(&duplicate,2,1,dp,dr,da,db,dx,0,&upper)==0);
  CHECK(duplicate && duplicate->eligible==-1);
  kls_componentwise_plan_destroy(duplicate);CHECK(live==0);
  context c[2]={0};pthread_t workers[2];
  for(int i=0;i<2;++i) CHECK(pthread_create(workers+i,NULL,independent,c+i)==0);
  for(int i=0;i<2;++i) CHECK(pthread_join(workers[i],NULL)==0 && !c[i].failed);
  CHECK(c[0].plan!=c[1].plan && c[0].plan->rowptr!=c[1].plan->rowptr);
  for(int i=0;i<2;++i) kls_componentwise_plan_destroy(c[i].plan);
  CHECK(live==0);
#endif
  puts("PASS certificate row-plan faults, selection, lifecycle and independent concurrency");
  return 0;
}
