/* Conservative hardware interval certificate shared by recovery and the
 * original-system caller gate.
 * No ambient rounding-mode changes. Each elementary result is enclosed by
 * conservative outward long-double bounds. A false result means NOT CERTIFIED,
 * never proof that the answer is inaccurate. Uncertain caller certificates
 * fall back to the unchanged binary128 implementation. */
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <fenv.h>
#if defined(__SSE__)
#include <xmmintrin.h>
#endif

/* Called only under the round-to-nearest/environment checks below.
 * In radix two, EPSILON is an exact power of two. Unless it underflows,
 * |value|*EPSILON is exact and at least the larger adjacent spacing at value.
 * The normal floor covers underflow and zero. Nearest-rounded addition of
 * this delta therefore reaches or exceeds nextafter in the requested
 * direction. The enclosure can be wider, never narrower; uncertain cases
 * still fall back to the caller's binary128 certificate. */
static long double interval_outward(long double value,long double direction) {
#if FLT_RADIX == 2
  if(isfinite(value)) {
    long double delta=fabsl(value)*LDBL_EPSILON;
    if(delta<LDBL_MIN) delta=LDBL_MIN;
    return direction>0?value+delta:value-delta;
  }
#endif
  return nextafterl(value,direction);
}

static int interval_csc_certificate(size_t n, const int64_t *ptr,
    const int64_t *rows, const double *a, const double *b, const double *x,
    int transpose, long double *upper_squared_out,
    long double *budget_squared_out) {
#if defined(__FAST_MATH__) || \
    (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__ > 0) || \
    LDBL_MANT_DIG < DBL_MANT_DIG || LDBL_HAS_SUBNORM != 1
  (void)n;(void)ptr;(void)rows;(void)a;(void)b;(void)x;(void)transpose;
  (void)upper_squared_out;(void)budget_squared_out;
  return 0;
#else
  if (fegetround()!=FE_TONEAREST || !ptr || !b || !x ||
      n>SIZE_MAX/(2*sizeof(long double))) return 0;
#if defined(__SSE__)
#if defined(__x86_64__) && LDBL_MANT_DIG == 64 && LDBL_MAX_EXP == 16384
  /* x87 extended arithmetic is independent of SSE FTZ/DAZ. Require full
   * precision and nearest rounding; volatile extended inputs below prevent
   * double comparisons from silently treating subnormal inputs as zero. */
  unsigned short control;
  __asm__ volatile ("fnstcw %0" : "=m" (control));
  if ((control & 0x0f00u) != 0x0300u) return 0;
#else
  if (_mm_getcsr() & ((1u<<15)|(1u<<6))) return 0;
#endif
#endif
  if (ptr[0]!=0 || ptr[n]<0 || (ptr[n] && (!a || !rows))) return 0;
  long double *lo=malloc((n ? 2*n : 2)*sizeof(*lo));
  if (!lo) return 0;
  long double *hi=lo+n;
  int certified=0, zero_rhs=1;
  long double bn_lower=0, rn_upper=0;
  for(size_t i=0;i<n;++i) {
    volatile long double bi=(long double)b[i], xi=(long double)x[i];
    if(!isfinite(bi) || !isfinite(xi)) goto done;
    lo[i]=hi[i]=bi;
    if(bi!=0) zero_rhs=0;
    long double sq=interval_outward(bi*bi,-INFINITY);
    if(!isfinite(sq)) goto done;
    if(sq<0) sq=0;
    bn_lower=interval_outward(bn_lower+sq,-INFINITY);
    if(bn_lower<0) bn_lower=0;
  }
  for(size_t col=0;col<n;++col) {
    if(ptr[col]<0 || ptr[col]>ptr[col+1] || ptr[col+1]>ptr[n]) goto done;
    for(size_t p=(size_t)ptr[col];p<(size_t)ptr[col+1];++p) {
      if(rows[p]<0 || (uintmax_t)rows[p]>=n) goto done;
      const size_t row=(size_t)rows[p];
      size_t dest=transpose?col:row;
      size_t source=transpose?row:col;
      volatile long double ai=(long double)a[p], xi=(long double)x[source];
      /* Every x component was validated above, including empty columns.
       * Keep its extended load for FTZ/DAZ independence, not another check
       * for each nonzero that references the same immutable component. */
      if(!isfinite(ai)) goto done;
      const long double product=ai*xi;
      if(!isfinite(product)) goto done;
      const long double lower=interval_outward(product,-INFINITY);
      const long double upper=interval_outward(product,INFINITY);
      lo[dest]=interval_outward(lo[dest]-upper,-INFINITY);
      hi[dest]=interval_outward(hi[dest]-lower,INFINITY);
      if(!isfinite(lo[dest]) || !isfinite(hi[dest])) goto done;
    }
  }
  for(size_t i=0;i<n;++i) {
    const long double r=fmaxl(fabsl(lo[i]),fabsl(hi[i]));
    const long double sq=interval_outward(r*r,INFINITY);
    rn_upper=interval_outward(rn_upper+sq,INFINITY);
    if(!isfinite(rn_upper)) goto done;
  }
  /* Literal rounding is also enclosed: do not silently raise 1e-8. */
  const long double tol_lower=interval_outward(1.0e-8L,-INFINITY);
  const long double tol2_lower=interval_outward(tol_lower*tol_lower,-INFINITY);
  const long double denominator=zero_rhs?1:bn_lower;
  if(!(denominator>0) || !isfinite(denominator)) goto done;
  const long double budget=interval_outward(tol2_lower*denominator,-INFINITY);
  if(upper_squared_out) *upper_squared_out=rn_upper;
  if(budget_squared_out) *budget_squared_out=budget;
  certified=budget>=0 && rn_upper<=budget;
done:
  free(lo);
  return certified;
#endif
}
