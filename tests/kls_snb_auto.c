#define KLS_SNB_AUTO_TEST_SECONDS(arm,actual) ((arm)?.001:.002)
#include "../src/kls.c"
#undef NDEBUG
#include <assert.h>
static void exercise(int transpose) {
 enum{N=129};int64_t ap[N+1],ai[N*N];double a[N*N],b[N],x[N];
 for(int c=0;c<N;++c){ap[c]=c*N;for(int r=0;r<N;++r){ai[c*N+r]=r;a[c*N+r]=r==c?N+1.:.03125;}}ap[N]=N*N;
 kls_options o;kls_default_options(&o);o.threads=2;o.use_btf=0;o.scale=0;
 o.ordering=KLS_ORDERING_NATURAL;o.orientation=transpose?KLS_ORIENTATION_TRANSPOSE:KLS_ORIENTATION_NORMAL;
 kls_solver *s=NULL;assert(kls_create(&s)==KLS_OK);
 assert(kls_analyze_csc(s,KLS_INDEX_INT64,N,ap,ai,0,&o)==KLS_OK);
 assert(kls_factor(s,a)==KLS_OK);assert(s->snb_auto.choice==1);
 const kls_snb_auto_policy prepared=s->snb_auto;
 for(int other=0;other<3;++other) {
  int *active=other==0?&s->transpose_cycle.active:
    other==1?&s->observed_parallel.active:&s->native_i16_policy.active;
  assert(!*active);*active=1;s->snb_auto.choice=0;
  double elapsed=0;assert(!kls_snb_auto_factor_prepare(s,&elapsed));
  assert(s->snb_auto.choice==0);*active=0;s->snb_auto=prepared;
 }
 assert(!kls_compact_transpose_cycle_capable(s));
 kls_observed_parallel_state observed=s->observed_parallel;
 (void)kls_observed_parallel_before(s,NULL,1,NULL,NULL);
 assert(!memcmp(&observed,&s->observed_parallel,sizeof observed));

 for(int epoch=1;epoch<=12;++epoch) {
  memset(b,0,sizeof b);
  for(int c=0;c<N;++c)for(int r=0;r<N;++r){
   a[c*N+r]=r==c?N+1.+.125*epoch:.03125*(1+.0625*epoch);
   b[r]+=a[c*N+r]*(1.+.125*((c+epoch)%7));
  }
  assert(kls_refactor_solve(s,a,1,b,0,x,0)==KLS_OK);
  for(int c=0;c<N;++c)assert(fabs(x[c]-(1.+.125*((c+epoch)%7)))<1e-11);
  if(epoch<=8)assert(s->snb_auto.phase==epoch);
 }
 assert(s->snb_auto.choice==3 && s->snb_auto_candidate);
 assert(s->snb_auto_candidate->stream_jobs>=4);
 if(!transpose)assert(s->native_i16_normal_runs>=4);
 /* An independent schedule must survive retirement of the candidate. */
 kls_snb_state *ordinary=kls_snb_build(s,NULL,0,SIZE_MAX);assert(ordinary);
 assert(ordinary!=s->snb_auto_candidate);s->snb=ordinary;
 s->private_caller_depth=1;
 (void)kls_snb_auto_before(s,a,1,kls_now_seconds());
 s->private_caller_depth=0;
 assert(s->snb==ordinary && !s->snb_auto_candidate && s->snb_auto.choice<0);
 kls_destroy(s);puts("PASS AUTO stream balanced certified calls and current RHS");
}

int main(void){exercise(0);exercise(1);return 0;}
