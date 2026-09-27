#include "../src/kls.c"
#undef NDEBUG
#include <assert.h>

int main(void) {
  int64_t p[]={0,2,4,6},rows[]={0,2,0,1,1,2};
  double a[]={4,1,2,8,1,16};
  int count=0;
  for(int orientation=0;orientation<3;++orientation)
  for(int scaling=-2;scaling<=2;scaling+=2)
  for(int threads=1;threads<=8;threads*=8) {
    kls_solver *s=NULL;kls_options o;kls_default_options(&o);
    o.orientation=orientation;o.scale=scaling;o.threads=threads;
    assert(kls_create(&s)==KLS_OK);
    assert(kls_set_accuracy_policy(s,KLS_ACCURACY_STRICT_RHS_L2)==KLS_OK);
    assert(kls_analyze_csc(s,KLS_INDEX_INT64,3,p,rows,0,&o)==KLS_OK);
    assert(kls_factor(s,a)==KLS_OK);
    for(int transpose=0;transpose<2;++transpose) {
      double b[3]={0,0,0},x[]={1.0001,.9999,1.0002};
      for(int j=0;j<3;++j)for(int64_t k=p[j];k<p[j+1];++k)
        b[transpose?j:rows[k]]+=a[k];
      double ratio=INFINITY;
      assert(!kls_accuracy_certify(3,p,rows,a,b,x,transpose,&ratio));
      double *saved=malloc(sizeof b);assert(saved);memcpy(saved,b,sizeof b);
      assert(private_caller_finish(s,transpose,1,saved,x,0,KLS_ERR_SOLVE_FAILED)==KLS_OK);
      assert(kls_accuracy_certify(3,p,rows,a,b,x,transpose,&ratio));
      ++count;
    }
    /* An exhausted allowance is shared across bridges and cannot renew at
     * final certification. An already-correct answer still needs no search. */
    double b[]={6,9,17},bad[]={1.0001,.9999,1.0002},before[3];
    memcpy(before,bad,sizeof bad);
    s->private_caller_depth=1;
    s->private_cold_deadline=kls_now_seconds()-1.;
    const double expired=s->private_cold_deadline;
    assert(private_recovery_deadline(s)==expired);
    assert(!private_lattice_recover(s,a,b,bad,0));
    assert(!memcmp(before,bad,sizeof bad));
    double *saved=malloc(sizeof b);assert(saved);memcpy(saved,b,sizeof b);
    assert(private_caller_finish(s,0,1,saved,bad,0,KLS_ERR_SOLVE_FAILED)==KLS_ERR_SOLVE_FAILED);
    assert(s->private_cold_deadline==expired && !memcmp(before,bad,sizeof bad));
    double good[]={1,1,1};saved=malloc(sizeof b);assert(saved);memcpy(saved,b,sizeof b);
    assert(private_caller_finish(s,0,1,saved,good,0,KLS_OK)==KLS_OK);
    assert(s->private_cold_deadline==expired);
    s->private_caller_depth=0;
    assert(kls_solve(s,1,b,0,good,0)==KLS_OK);
    assert(s->private_caller_depth==0 && s->private_cold_deadline!=expired);
    s->private_cold_deadline=expired;
    assert(kls_refactor_solve(s,a,1,b,0,good,0)==KLS_OK);
    assert(s->private_caller_depth==0 && s->private_cold_deadline!=expired);
    kls_destroy(s);
  }
  printf("PASS caller recovery: %d oriented/scaled/thread configurations\n",count);
  return 0;
}
