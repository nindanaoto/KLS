#define _GNU_SOURCE
#define _POSIX_C_SOURCE 200809L
#undef NDEBUG
#include <stdlib.h>
#include <assert.h>
static size_t fail_size;
static int fault_seen;
static void *test_malloc(size_t n) {
  if(fail_size && n==fail_size) { fault_seen=1; return NULL; }
  return malloc(n);
}
#define malloc test_malloc
#include "../src/kls.c"
#undef malloc
#undef NDEBUG
#include <assert.h>

int main(void) {
  UF_long ptr[]={0,1,2},rows[]={0,1};
  double a[]={1,1},b[]={1,1},good[]={1-8e-9,1},bad[]={2,2},best[]={0,0};
  kls_solver plain={0}; plain.n=2;plain.col_ptr=ptr;plain.row_idx=rows;
  for(int transpose=0;transpose<2;++transpose) {
    long double upper=INFINITY;
    int certified=0;
    private_keep_best(&plain,a,b,good,transpose,best,&upper,&certified);
    assert(isfinite(upper) && certified);
    long double saved=upper;
    /* Re-ranking the identical incumbent must neither allocate nor discard
       its bound. The final acceptance certificate is deliberately separate. */
    fault_seen=0;fail_size=4*sizeof(long double);
    private_keep_best(&plain,a,b,good,transpose,best,&upper,&certified);
    assert(!fault_seen && upper==saved && certified);
    assert(!interval_csc_certificate(2,ptr,rows,a,b,best,transpose,NULL,NULL));
    assert(fault_seen);
    fail_size=0;
    private_keep_best(&plain,a,b,bad,transpose,best,&upper,&certified);
    assert(upper==saved && certified && !memcmp(best,good,sizeof(best)));
    assert(interval_csc_certificate(2,ptr,rows,a,b,best,transpose,NULL,NULL));
    /* An unavailable first bound is not a reusable ranking result. */
    upper=INFINITY;certified=0;fault_seen=0;fail_size=4*sizeof(long double);
    private_keep_best(&plain,a,b,good,transpose,best,&upper,&certified);
    assert(fault_seen && !isfinite(upper) && !certified);
    fail_size=0;
    private_keep_best(&plain,a,b,good,transpose,best,&upper,&certified);
    assert(isfinite(upper) && certified);
    /* A finite ranking bound is not by itself an acceptance proof. */
    upper=INFINITY;certified=0;
    private_keep_best(&plain,a,b,bad,transpose,best,&upper,&certified);
    assert(isfinite(upper) && !certified);
    private_keep_best(&plain,a,b,good,transpose,best,&upper,&certified);
    assert(certified && !memcmp(best,good,sizeof(best)));
  }
  enum {N=17};
  int64_t cp[N+1],ri[N*N];double av[N*N],rhs[N],x[N];
  memset(rhs,0,sizeof(rhs));
  for(int j=0;j<N;++j) {
    cp[j]=N*j;
    for(int i=0;i<N;++i) {
      int p=j*N+i;ri[p]=i;av[p]=i==j?5.:1./(i+j+3);
      rhs[i]+=av[p];
    }
  }
  cp[N]=N*N;
  kls_solver *s=NULL;kls_options o;kls_default_options(&o);
  o.threads=1;o.scale=0;o.orientation=KLS_ORIENTATION_NORMAL;
  assert(kls_create(&s)==KLS_OK);
  assert(kls_set_accuracy_policy(s,KLS_ACCURACY_STRICT_RHS_L2)==KLS_OK);
  assert(kls_analyze_csc(s,KLS_INDEX_INT64,N,cp,ri,0,&o)==KLS_OK);
  assert(kls_factor(s,av)==KLS_OK);
  /* Isolated cold proposals never install their factors in the parent. */
  void *original_numeric=s->numeric;
  double cold_deadline=kls_now_seconds()+100;
  s->private_cold_deadline=cold_deadline;
  for(int j=0;j<N;++j) x[j]=-7;
  double unchanged[N];memcpy(unchanged,x,sizeof x);
  assert(!private_fresh_proposal(s,N,cp,ri,av,rhs,x,kls_now_seconds()-1));
  assert(!memcmp(x,unchanged,sizeof x));
  fault_seen=0;fail_size=sizeof x;
  assert(!private_fresh_proposal(s,N,cp,ri,av,rhs,x,cold_deadline));
  fail_size=0;
  assert(fault_seen && !memcmp(x,unchanged,sizeof x));
  s->private_cold_trial=1;
  assert(!private_fresh_proposal(s,N,cp,ri,av,rhs,x,cold_deadline));
  s->private_cold_trial=0;
  assert(!memcmp(x,unchanged,sizeof x));
  assert(private_fresh_proposal(s,N,cp,ri,av,rhs,x,cold_deadline));
  double cold_bound;
  assert(kls_accuracy_certify(N,cp,ri,av,rhs,x,0,&cold_bound));
  assert(s->numeric==original_numeric && s->private_cold_deadline==cold_deadline);
  /* Inconsistent singular trial: failure must retain the entire seed. */
  int64_t singular_p[]={0,1,2},singular_rows[]={0,1};
  double singular_a[]={0,0},singular_b[]={1,1},seed[]={7,9};
  assert(!private_fresh_proposal(s,2,singular_p,singular_rows,
                                singular_a,singular_b,seed,cold_deadline));
  assert(seed[0]==7 && seed[1]==9 && s->numeric==original_numeric);
  /* Exercise caller-frame transpose assembly with a nonsymmetric operator. */
  int64_t transpose_p[]={0,1,3},transpose_rows[]={0,0,1};
  double transpose_a[]={2,1,3},transpose_b[]={2,4},transpose_x[]={0,0};
  private_original_operator transpose_operator={0};
  transpose_operator.n=2;transpose_operator.nnz=3;
  transpose_operator.p=transpose_p;transpose_operator.rows=transpose_rows;
  transpose_operator.values=transpose_a;
  assert(private_caller_recover(s,&transpose_operator,transpose_b,transpose_x,1,cold_deadline));
  assert(kls_accuracy_certify(2,transpose_p,transpose_rows,transpose_a,
                             transpose_b,transpose_x,1,&cold_bound));
  assert(s->numeric==original_numeric && s->private_cold_deadline==cold_deadline);
  /* Fail the factor-side snapshot allocation itself, then verify the public
     solve cannot fall back to stale analyze-time values. */
  free(s->solve_refine_values);s->solve_refine_values=NULL;
  s->solve_contract_probe=0;
  fault_seen=0;fail_size=sizeof(av);
  kls_factor_solve_contract_classify(s,av);
  fail_size=0;
  assert(fault_seen && s->solve_contract_probe==2 && !s->solve_refine_values);
  for(int transpose=0;transpose<2;++transpose) {
    int status=transpose?kls_solve_transpose(s,1,rhs,N,x,N)
                        :kls_solve(s,1,rhs,N,x,N);
    assert(status==KLS_ERR_OUT_OF_MEMORY);
  }
  kls_factor_solve_contract_classify(s,av);
  assert(s->solve_refine_values);
  /* Aliased RHS preservation is mandatory when verification is armed. */
  for(int transpose=0;transpose<2;++transpose) {
    memcpy(x,rhs,sizeof x);fault_seen=0;fail_size=sizeof x;
    int status=transpose?kls_solve_transpose(s,1,x,N,x,N)
                        :kls_solve(s,1,x,N,x,N);
    fail_size=0;
    assert(fault_seen && status==KLS_ERR_OUT_OF_MEMORY);
    assert(!memcmp(x,rhs,sizeof x));
  }
  s->row_solve_self_check=1;
  free(s->solve_refine_values);
  s->solve_refine_values=malloc(sizeof(av));
  assert(s->solve_refine_values);
  memcpy(s->solve_refine_values,av,sizeof(av));
  free(s->solve_refine_workspace);s->solve_refine_workspace=NULL;
  s->solve_contract_probe=1;
  s->solve_contract_verified=1;
  s->low_rcond_solve_contract_state=1;
  s->verified_rhs_valid=1;
  fault_seen=0;fail_size=4*N*sizeof(double);
  int status=kls_solve(s,1,rhs,N,x,N);
  fail_size=0;
  fprintf(stderr,"allocation fault_seen=%d status=%d expected=%d\n",fault_seen,status,KLS_ERR_OUT_OF_MEMORY);
  assert(fault_seen && status==KLS_ERR_OUT_OF_MEMORY);
  assert(s->solve_contract_probe==0 && !s->solve_contract_verified &&
         !s->low_rcond_solve_contract_state && !s->verified_rhs_valid);
  kls_destroy(s);
  puts("PASS best-bound preservation and required-refinement allocation failure");
  return 0;
}
