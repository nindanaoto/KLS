#define _GNU_SOURCE
#undef NDEBUG
#include "../src/kls_policy.h"
#include "kls/kls.h"
#include <assert.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

static void *reserve_worker(void *arg) {
  kls_policy_context *p=arg;
  for(int i=0;i<100;++i)if(kls_policy_reserve(p,.001))
    kls_policy_charge(p,.001,.001);
  return NULL;
}
static void budget_checks(void) {
  kls_policy_context p;kls_policy_init(&p,NULL);
  assert(!kls_policy_model_only(&p));
  p.cost.enabled=1;assert(kls_policy_model_only(&p));
  p.kernel_overrides=1;assert(!kls_policy_model_only(&p));
  kls_policy_init(&p,NULL);
  assert(!kls_policy_reserve(&p,1));
  kls_policy_fund(&p,NAN);assert(!p.funded);
  kls_policy_fund(&p,1);double limit=p.limit;
  kls_policy_fund(&p,1e12);assert(p.limit==limit);
  assert(!kls_policy_reserve(&p,NAN));assert(!kls_policy_reserve(&p,-1));
  pthread_t workers[8];for(int i=0;i<8;++i)assert(!pthread_create(&workers[i],NULL,reserve_worker,&p));
  for(int i=0;i<8;++i)assert(!pthread_join(workers[i],NULL));
  assert(p.reserved==0 && p.spent<=p.limit && p.admitted>0);
  assert(!kls_policy_reserve(&p,.001));
  kls_policy_init(&p,NULL);kls_policy_fund(&p,1);
  assert(kls_policy_reserve(&p,.001));kls_policy_charge(&p,.001,.1);
  assert(p.spent>p.limit && !kls_policy_reserve(&p,.001));
  kls_policy_init(&p,NULL);kls_policy_fund(&p,1);
  kls_policy_charge(&p,0,NAN);assert(!kls_policy_reserve(&p,.001));

}
static void estimate_checks(void) {
  kls_policy_context p;kls_policy_init(&p,NULL);
  assert(isnan(kls_policy_numeric_seconds(&p,10,20,-1,1,1,1)));
  assert(isnan(kls_policy_numeric_seconds(&p,10,DBL_MAX,0,1,1,1)));
  assert(isnan(kls_policy_numeric_seconds(&p,10,20,INFINITY,1,1,1)));
  assert(isnan(kls_policy_graph_seconds(&p,10,20,NAN)));
  double cold=kls_policy_numeric_seconds(&p,10,20,100,1,0,0);
  double factor=kls_policy_numeric_seconds(&p,10,20,100,1,1,0);
  double solve=kls_policy_numeric_seconds(&p,10,20,100,1,0,1);
  assert(cold>0 && factor>cold && solve>cold && factor!=solve);
  assert(kls_policy_graph_seconds(&p,10,20,100)>
         kls_policy_graph_seconds(&p,10,20,10));
  p.cost.seconds[KLS_POLICY_FIRST_REFACTOR_ENTRY]=1e-3;
  double t0=kls_policy_numeric_seconds(&p,10,20,100,1,0,0);
  double t1=kls_policy_numeric_seconds(&p,10,20,100,1,1,0);
  double t2=kls_policy_numeric_seconds(&p,10,20,100,1,2,0);
  double t3=kls_policy_numeric_seconds(&p,10,20,100,1,3,0);
  assert(t1-t0>t2-t1);
  assert(fabs((t3-t2)-(t2-t1))<1e-12);
}
static void profile_checks(void) {
  char path[]="/tmp/kls-policy-test-XXXXXX";int fd=mkstemp(path);assert(fd>=0);close(fd);unlink(path);
  kls_tuning_values v;kls_tuning_defaults(&v);kls_tuning_metadata m={0};
  assert(kls_tuning_write_host_profile(path,1,&v,17,0,NULL)==KLS_OK);
  assert(kls_tuning_load_host_profile(path,1,&v,&m)==KLS_OK && !m.policy_cost.enabled);
  unlink(path);
  kls_policy_context p;kls_policy_init(&p,NULL);p.cost.enabled=1;
  p.cost.cache_bytes=kls_tuning_policy_cache_bytes();
  if(!p.cost.cache_bytes)return; /* Portable library needs no Linux topology. */
  assert(kls_tuning_write_policy_profile(path,1,&p.cost,17,0)==KLS_OK);
  assert(kls_tuning_load_host_profile(path,1,&v,&m)==KLS_OK);
  assert(m.policy_cost.enabled && !m.kernel_overrides);
  for(int i=0;i<KLS_POLICY_COST_COUNT;++i)assert(m.policy_cost.seconds[i]==p.cost.seconds[i]);
  assert(kls_tuning_load_host_profile(path,2,&v,&m)==KLS_ERR_TUNING_PROFILE);
  assert(kls_tuning_write_policy_profile(path,1,&p.cost,17,0)==KLS_ERR_TUNING_PROFILE);
  p.cost.seconds[0]=NAN;
  assert(kls_tuning_write_policy_profile(path,1,&p.cost,17,1)==KLS_ERR_INVALID_ARGUMENT);
  assert(kls_tuning_load_host_profile(path,1,&v,&m)==KLS_OK);
  FILE *f=fopen(path,"a");assert(f);fputs("policy_cost_version=1\n",f);fclose(f);
  assert(kls_tuning_load_host_profile(path,1,&v,&m)==KLS_ERR_TUNING_PROFILE);
  unlink(path);
}
int main(void){budget_checks();estimate_checks();profile_checks();puts("policy budget, concurrency, estimates and profile checks passed");return 0;}
