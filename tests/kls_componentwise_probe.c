/* Test-only access: no diagnostic entry points in the production library. */
#include <stdio.h>
#include "../src/kls_componentwise.c"

int kls_test_fast_available(void) { return KLS_COMPONENTWISE_HAS_ACCUMULATED; }

int kls_test_fast(size_t n, const int64_t *p, const int64_t *rows,
    const double *a, const double *b, const double *x, int transpose,
    double *residual, long double *upper) {
  (void)residual;
  *upper=INFINITY;
  /* Fixtures supply valid CSC and the default arithmetic environment. */
  return componentwise_accumulated_certificate(n,p,rows,a,b,x,transpose,upper);
}

int kls_test_bound(const long double *upper, char *out, size_t capacity) {
  /* Hexadecimal text retains every significand bit, unlike ctypes.value. */
  return snprintf(out,capacity,"%.*La",(LDBL_MANT_DIG+3)/4,*upper);
}

int kls_test_parallel_available(void) {
#if KLS_COMPONENTWISE_HAS_ACCUMULATED && defined(_OPENMP)
  return 1;
#else
  return 0;
#endif
}

int kls_test_parallel(size_t n,const int64_t *p,const int64_t *rows,
    const double *a,const double *b,const double *x,int transpose,
    double *residual,long double *upper,int threads) {
  (void)residual;*upper=INFINITY;
#if KLS_COMPONENTWISE_HAS_ACCUMULATED && defined(_OPENMP)
  kls_componentwise_plan plan={0};
  if(!n) {*upper=0;return 1;}
  componentwise_plan_build(&plan,n,p,rows);
  int result=plan.eligible>0?componentwise_plan_run(&plan,threads,p,rows,a,b,x,transpose,upper):0;
  kls_componentwise_plan_reset(&plan);
  return result;
#else
  (void)n;(void)p;(void)rows;(void)a;(void)b;(void)x;(void)transpose;(void)threads;
  return 0;
#endif
}
