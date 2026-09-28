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
