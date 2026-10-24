#define _GNU_SOURCE
#define _POSIX_C_SOURCE 200809L
#include "../src/kls.c"
#undef NDEBUG
#include <assert.h>
#include "kls_scaled_recovery_fixtures.h"

static void witness(size_t n,size_t nz,const int64_t *p,const int64_t *j,
    const double *a,const double *b) {
  double *x=malloc(2*(n+1)*sizeof(double));
  double *rhs=malloc(2*(n+1)*sizeof(double));
  double *updated=malloc(nz*sizeof(double));
  assert(x && rhs && updated);
  kls_options o;kls_default_options(&o);
  o.ordering=KLS_ORDERING_METIS;o.use_btf=0;o.threads=1;
  for(int mode=0;mode<3;++mode) {
    o.scale=mode==0?-1:mode==1?2:-2;
    kls_solver *s=NULL;assert(kls_create(&s)==KLS_OK);
    assert(kls_analyze_csr(s,KLS_INDEX_INT64,n,p,j,0,&o)==KLS_OK);
    assert(kls_factor(s,a)==KLS_OK);
    for(size_t i=0;i<n;++i) {rhs[i]=b[i];rhs[n+1+i]=2*b[i];}
    rhs[n]=rhs[2*n+1]=77;
    for(size_t i=0;i<2*(n+1);++i) x[i]=91;
    /* Stand in for a loaded numeric profile, after the incumbent factor. */
    if(mode==2) ++s->tuning.snode_min_batch;
    trilinos_klu_l_numeric *incumbent=s->numeric;
    int status=kls_solve(s,2,rhs,n+1,x,n+1);
    assert(s->numeric==incumbent);
    if(mode==0) {
      assert(status==KLS_ERR_SOLVE_FAILED && !s->private_scaled_recovery);
      for(size_t i=0;i<2*(n+1);++i) assert(x[i]==91);
      kls_destroy(s);continue;
    }
    assert(status==KLS_OK && x[n]==91 && x[2*n+1]==91);
    for(int k=0;k<2;++k) {
      long double bound;
      assert(kls_componentwise_certify(n,p,j,a,rhs+k*(n+1),x+k*(n+1),
          1,NULL,&bound)==1);
    }
    if(mode==1) {assert(!s->private_scaled_recovery);kls_destroy(s);continue;}
    kls_solver *retained=s->private_scaled_recovery;
    assert(s->private_scaled_recovery_builds==1);
    assert(retained && retained->common.scale==2 && retained->private_cold_trial);
    assert(retained->tuning.snode_min_batch==s->tuning.snode_min_batch);
    assert(!retained->private_scaled_recovery && s->options.scale==-2);
    /* Reuse only a factor, never a previous answer or RHS certificate. */
    assert(kls_solve(s,2,rhs,n+1,rhs,n+1)==KLS_OK);
    assert(s->private_scaled_recovery==retained && rhs[n]==77);
    const private_original_operator *original=s->private_original;
    memset(rhs,0,n*sizeof(double));
    for(size_t c=0;c<n;++c)
      for(int64_t k=original->p[c];k<original->p[c+1];++k)
        rhs[c]+=original->values[k]; /* A^T * ones */
    assert(kls_solve_transpose(s,1,rhs,0,x,0)==KLS_OK);
    long double bound;
    assert(kls_componentwise_certify(n,original->p,original->rows,
        original->values,rhs,x,1,NULL,&bound)==1);
    /* Failure in a later RHS cannot publish the earlier RHS or replace the
     * retained factor. Nonfinite input remains a failure under DEFAULT. */
    memcpy(rhs,b,n*sizeof(double));memcpy(rhs+n+1,b,n*sizeof(double));
    rhs[n+1]=NAN;
    for(size_t i=0;i<2*(n+1);++i) x[i]=91;
    assert(kls_solve(s,2,rhs,n+1,x,n+1)!=KLS_OK);
    for(size_t i=0;i<2*(n+1);++i) assert(x[i]==91);
    assert(s->private_scaled_recovery==retained);
    /* Public updates invalidate numeric validity but retain compatible analysis. */
    for(size_t i=0;i<nz;++i) updated[i]=2*a[i];
    assert(kls_refactor(s,updated)==KLS_OK && s->private_scaled_recovery==retained);
    assert(s->private_scaled_recovery_pending);
    for(size_t i=0;i<n;++i) rhs[i]=2*b[i];
    assert(kls_solve(s,1,rhs,0,x,0)==KLS_OK);
    assert(kls_componentwise_certify(n,p,j,updated,rhs,x,1,NULL,&bound)==1);
    assert(s->private_scaled_recovery==retained && !s->private_scaled_recovery_pending);
    assert(s->private_scaled_recovery_updates==1 && s->private_scaled_recovery_accepts==1);
    assert(s->private_scaled_recovery_builds==1);
    assert(kls_factor(s,a)==KLS_OK && s->private_scaled_recovery==retained);
    assert(s->private_scaled_recovery_pending);
    assert(kls_solve(s,1,b,0,x,0)==KLS_OK && s->private_scaled_recovery);
    assert(kls_analyze_csr(s,KLS_INDEX_INT64,n,p,j,0,&o)==KLS_OK);
    assert(!s->private_scaled_recovery);
    assert(kls_factor(s,a)==KLS_OK);
    assert(kls_refactor_solve(s,a,1,b,0,x,0)==KLS_OK);
    assert(kls_componentwise_certify(n,p,j,a,b,x,1,NULL,&bound)==1);
    kls_destroy(s);
  }
  free(x);free(rhs);free(updated);
}

int main(void) {
  witness(34,147,fixture_p0,fixture_j0,fixture_a0,fixture_b0);
  witness(58,349,fixture_p1,fixture_j1,fixture_a1,fixture_b1);
  puts("PASS scaled recovery: witnesses, multi-RHS, aliasing, transpose, epochs, rollback");
  return 0;
}
