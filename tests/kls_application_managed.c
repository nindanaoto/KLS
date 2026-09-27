#include "../src/kls.c"
#undef NDEBUG
#include <assert.h>
int main(void) {
  int64_t p[]={0,1,2},r[]={0,1};double a[]={2,4};
  for(int orientation=0;orientation<3;++orientation) for(int api=0;api<3;++api) {
    kls_solver *s=NULL;kls_options o;kls_default_options(&o);o.orientation=orientation;
    assert(kls_create(&s)==KLS_OK);
    kls_accuracy_policy policy;
    assert(kls_get_accuracy_policy(s,&policy)==KLS_OK && policy==KLS_ACCURACY_COMPONENTWISE_BACKWARD_ERROR);
    assert(kls_set_accuracy_policy(s,KLS_ACCURACY_APPLICATION_MANAGED)==KLS_OK);
    assert(kls_analyze_csc(s,KLS_INDEX_INT64,2,p,r,0,&o)==KLS_OK);
    assert(kls_set_accuracy_policy(s,KLS_ACCURACY_COMPONENTWISE_BACKWARD_ERROR)==KLS_ERR_INVALID_ARGUMENT);
    assert(kls_factor(s,a)==KLS_OK);
    double x[]={2,4,77,4,8,88};
    int status=api==2?kls_refactor_solve(s,a,2,x,3,x,3):
      api?kls_solve_transpose(s,2,x,3,x,3):kls_solve(s,2,x,3,x,3);
    assert(status==KLS_OK && x[0]==1 && x[1]==1 && x[2]==77 && x[3]==2 && x[4]==2 && x[5]==88);
    double bad[]={2,4,0,NAN,8,0},before[6];memcpy(before,x,sizeof x);
    status=api==2?kls_refactor_solve(s,a,2,bad,3,x,3):
      api?kls_solve_transpose(s,2,bad,3,x,3):kls_solve(s,2,bad,3,x,3);
    assert(status!=KLS_OK && !memcmp(x,before,sizeof x));
    /* Application-managed acceptance is NOT a residual certificate. */
    for(int kind=0;kind<4;++kind) {
      double *saved=malloc(2*sizeof(*saved));assert(saved);saved[0]=2;saved[1]=4;
      double proposed[]={0,kind==1?NAN:kind==2?INFINITY:0};
      int input_status=kind==3?KLS_ERR_SOLVE_FAILED:KLS_OK;
      status=private_caller_finish(s,0,1,saved,proposed,0,input_status);
      assert(status==(kind==0?KLS_OK:KLS_ERR_SOLVE_FAILED));
      assert(s->private_cold_deadline==0);
    }
    double invalid[]={INFINITY,4};
    assert(kls_refactor(s,invalid)!=KLS_OK);
    kls_destroy(s);
  }
  return 0;
}
