#include "../src/kls.c"
#undef NDEBUG
#include <assert.h>
static void reset_owned_cache(void) {
  enum{N=32};int64_t ptr[N+1],rows[3*N];double values[3*N];int at=0;
  for(int c=0;c<N;++c){ptr[c]=at;for(int r=c-1;r<=c+1;++r)if(r>=0&&r<N){rows[at]=r;values[at++]=r==c?4:-.25;}}ptr[N]=at;
  kls_options options;kls_default_options(&options);options.use_btf=0;options.ordering=KLS_ORDERING_NATURAL;options.orientation=KLS_ORIENTATION_NORMAL;options.scale=0;
  kls_solver *s=NULL;assert(kls_create(&s)==KLS_OK);assert(kls_analyze_csc(s,KLS_INDEX_INT64,N,ptr,rows,0,&options)==KLS_OK);assert(kls_factor(s,values)==KLS_OK);
  s->native_i16=kls_native_i16_build(s,0,SIZE_MAX);assert(s->native_i16);
  kls_observed_parallel_state old;memcpy(&old,&s->observed_parallel,sizeof(old));
  unsigned epoch=s->native_i16_epoch;s->native_i16_combined_depth=3;
  s->native_i16_policy.choice=2;s->native_i16_policy.spent=.001;
  s->native_i16_serial_calls=17;s->native_i16_runs=5;s->native_i16_last_native=1;
  /* Even a solver without a forest owns the separate compact cache. */
  kls_free_forest_map(s);
  assert(!s->native_i16&&s->native_i16_policy.choice==0&&s->native_i16_epoch==epoch+1);
  assert(s->native_i16_combined_depth==3&&s->native_i16_serial_calls==17&&s->native_i16_runs==5&&!s->native_i16_last_native);
  assert(!memcmp(&old,&s->observed_parallel,sizeof(old)));
  s->native_i16_combined_depth=0;kls_destroy(s);
}
static void stale_consultation_does_not_poison_new_layout(void) {
  kls_solver *s=NULL;assert(kls_create(&s)==KLS_OK);
  kls_native_i16_call call={.outer=1,.epoch=s->native_i16_epoch,.sampled=1,.candidate_ok=1,.candidate_seconds=.001};
  call.incumbent_start=kls_now_seconds();
  kls_native_i16_invalidate(s);unsigned epoch=s->native_i16_epoch;
  kls_native_i16_after(s,call,KLS_OK,0);
  assert(s->native_i16_epoch==epoch&&!s->native_i16_policy.choice&&!s->native_i16_policy.active&&!s->native_i16_policy.spent);
  kls_destroy(s);
}
static void fresh_public_normal_calls_remain_unfunded(void) {
  enum{N=96,NZ=3*N};int64_t ptr[N+1],rows[NZ];double a[NZ],rhs[N],out[N];int at=0;
  for(int c=0;c<N;++c){ptr[c]=at;for(int r=c-1;r<=c+1;++r)if(r>=0&&r<N){rows[at]=r;a[at++]=r==c?4:-.25;}}ptr[N]=at;
  kls_options options;kls_default_options(&options);options.threads=1;options.use_btf=0;options.orientation=KLS_ORIENTATION_NORMAL;options.ordering=KLS_ORDERING_NATURAL;options.scale=0;
  options.expected_refactorizations=options.expected_solves=0;
  kls_solver *s=NULL;assert(kls_create(&s)==KLS_OK);assert(kls_set_accuracy_policy(s,KLS_ACCURACY_COMPONENTWISE_BACKWARD_ERROR)==KLS_OK);
  assert(kls_analyze_csc(s,KLS_INDEX_INT64,N,ptr,rows,0,&options)==KLS_OK);assert(kls_factor(s,a)==KLS_OK);
  for(int epoch=0;epoch<600;++epoch){memset(rhs,0,sizeof(rhs));
    for(int c=0;c<N;++c)for(int q=(int)ptr[c];q<(int)ptr[c+1];++q){a[q]=rows[q]==c?4+.03125*(epoch%5):-.25*(1+.015625*(epoch%3));rhs[rows[q]]+=a[q]*(1+.125*(c%7));}
    assert(kls_refactor_solve(s,a,1,rhs,0,out,0)==KLS_OK);
    for(int c=0;c<N;++c)assert(fabs(out[c]-(1+.125*(c%7)))<1e-12);
    assert(s->native_i16_combined_depth==0&&!s->native_i16);
  }
  assert(s->native_i16_policy.cycles==0&&s->native_i16_policy.spent==0&&s->native_i16_runs==0);
  kls_destroy(s);
}
int main(void){reset_owned_cache();stale_consultation_does_not_poison_new_layout();fresh_public_normal_calls_remain_unfunded();puts("PASS native mixed-width policy ownership, independent epochs, stale consultation and unfunded public normal calls");return 0;}
