#define _GNU_SOURCE
#define _POSIX_C_SOURCE 200809L
#include <stdlib.h>
static int fail_next;
static int fail_candidate;
static void *gate_malloc(size_t n) {
  if (fail_next) { fail_next = 0; return NULL; }
  return malloc(n);
}
static void *gate_calloc(size_t n,size_t size) {
  if (fail_candidate) { fail_candidate=0;return NULL; }
  return calloc(n,size);
}
#define malloc gate_malloc
#define calloc gate_calloc
#include "../src/kls.c"
#undef malloc
#undef calloc
#undef NDEBUG
#include <assert.h>

int main(void) {
  int64_t p[] = {0,1,2}, rows[] = {0,1};
  double values[] = {2,4};
  for (int orientation = 0; orientation < 3; ++orientation) {
    kls_options o; kls_default_options(&o);
    o.orientation = orientation; o.threads = 1;
    kls_solver *s = NULL;
    assert(kls_create(&s) == KLS_OK);
    assert(kls_set_accuracy_policy(s,KLS_ACCURACY_STRICT_RHS_L2)==KLS_OK);
    assert(kls_analyze_csc(s,KLS_INDEX_INT64,2,p,rows,0,&o) == KLS_OK);
    assert(kls_factor(s,values) == KLS_OK);
    for (int api = 0; api < 3; ++api) {
      double x[] = {2,4,77,4,8,88};
      double before[6]; memcpy(before,x,sizeof x);
      fail_next = 1;
      int status = api == 2 ? kls_refactor_solve(s,values,2,x,3,x,3) :
        api ? kls_solve_transpose(s,2,x,3,x,3) : kls_solve(s,2,x,3,x,3);
      assert(status == KLS_ERR_OUT_OF_MEMORY && !fail_next);
      assert(!memcmp(x,before,sizeof x));
      fail_candidate=1;
      status = api == 2 ? kls_refactor_solve(s,values,2,x,3,x,3) :
        api ? kls_solve_transpose(s,2,x,3,x,3) : kls_solve(s,2,x,3,x,3);
      assert(status == KLS_ERR_OUT_OF_MEMORY && !fail_candidate);
      assert(!memcmp(x,before,sizeof x));
      status = api == 2 ? kls_refactor_solve(s,values,2,x,3,x,3) :
        api ? kls_solve_transpose(s,2,x,3,x,3) : kls_solve(s,2,x,3,x,3);
      assert(status == KLS_OK);
      assert(x[0] == 1 && x[1] == 1 && x[2] == 77);
      assert(x[3] == 2 && x[4] == 2 && x[5] == 88);
    }
    double b[] = {2,4}, x[] = {0,0}, *saved = NULL;
    assert(private_caller_rhs(s,INT64_MAX,b,2,x,2,&saved) == KLS_ERR_INVALID_ARGUMENT);
    assert(private_caller_rhs(s,2,b,INT64_MAX,x,2,&saved) == KLS_ERR_INVALID_ARGUMENT);
    assert(private_caller_rhs(s,2,b,2,x,INT64_MAX,&saved) == KLS_ERR_INVALID_ARGUMENT);
    assert(private_caller_rhs(s,1,b,-1,x,2,&saved) == KLS_ERR_INVALID_ARGUMENT);
    assert(private_caller_rhs(s,1,b,1,x,2,&saved) == KLS_ERR_INVALID_ARGUMENT);
    assert(!saved);
    /* A numerically successful internal solve must not override disagreement
     * with the original operator. Inject disagreement only into this fixture. */
    s->private_original->values[0] *= 2;
    for(int transpose=0;transpose<2;++transpose) {
      int status=transpose?kls_solve_transpose(s,1,b,0,x,0):kls_solve(s,1,b,0,x,0);
      double ratio;
      assert(status==KLS_OK &&
        kls_accuracy_certify(2,s->private_original->p,s->private_original->rows,
                     s->private_original->values,b,x,transpose,&ratio));
    }
    /* An explicit public update replaces the deliberately corrupted snapshot,
     * including identical values. Then exercise the combined fast return. */
    assert(kls_refactor(s,values) == KLS_OK);
    assert(kls_refactor_solve(s,values,1,b,0,x,0) == KLS_OK);
    assert(x[0] == 1 && x[1] == 1);
    assert(private_caller_rhs(s,1,b,0,x,0,&saved) == KLS_OK);
    x[0] = NAN;
    /* With cold recovery suppressed, the gate must reject the bad answer.
       With isolated fresh repair enabled, this nonsingular system is solvable;
       success must describe the repaired vector, not the rejected seed. */
    s->private_cold_trial=1;
    assert(private_caller_finish(s,0,1,saved,x,0,KLS_OK) == KLS_ERR_SOLVE_FAILED);
    s->private_cold_trial=0;
    assert(private_caller_rhs(s,1,b,0,x,0,&saved) == KLS_OK);
    assert(private_caller_finish(s,0,1,saved,x,0,KLS_OK) == KLS_OK);
    assert(x[0] == 1 && x[1] == 1);
    kls_destroy(s);
  }
  puts("PASS caller gate: aliasing, strides, overflow, OOM, disagreement, renewal, nonfinite");
  return 0;
}
