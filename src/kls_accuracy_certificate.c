/* Binary128 sufficient certificate. Basic arithmetic uses compiler
 * runtime support, not an external solver or quadmath library. */
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <fenv.h>
#include <math.h>
#include <float.h>
#if defined(__FAST_MATH__) || \
    (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__ > 0)
/* Finite-only can be enabled without FAST_MATH and deletes explicit NaN/Inf
 * checks, including for solution components in empty matrix columns. */
#error KLS accuracy certificate requires strict finite/nonfinite semantics
#endif
#if !defined(__SIZEOF_FLOAT128__) || \
    __SIZEOF_FLOAT128__ != 16 || !defined(__SIZEOF_INT128__) || __SIZEOF_INT128__ != 16
#error IEEE binary128 and 128-bit integers with strict arithmetic are required
#endif
/* Clang exposes __float128 and its size but not GCC's FLT128 format macros.
 * Check advertised format properties where available; precision checks
 * below also cover that Clang frontend. */
#if (defined(__FLT128_MANT_DIG__) && __FLT128_MANT_DIG__ != 113) || \
    (defined(__FLT128_MAX_EXP__) && __FLT128_MAX_EXP__ != 16384)
#error Unsupported 128-bit floating format
#endif
typedef __float128 quad;
#include "kls_strict_interval.inc"
#include "kls_gamma_certificate.inc"
__extension__ typedef unsigned __int128 bits;
_Static_assert(sizeof(quad)==16 && sizeof(bits)==16,"binary128 bit layout required");
_Static_assert(sizeof(double)==8 && DBL_MANT_DIG==53 && DBL_MAX_EXP==1024,
               "IEEE binary64 inputs required");
/* Exact bit conversion avoids any dependence on the caller's DAZ mode or
 * on the long-double format. Normal and subnormal binary64 values fit exactly
 * in binary128, so no floating conversion or rounding is needed. */
static quad from_double(double value) {
  uint64_t raw;memcpy(&raw,&value,sizeof raw);
  unsigned exponent=(unsigned)((raw>>52)&2047);
  uint64_t fraction=raw&UINT64_C(0x000fffffffffffff);
  bits encoded=(bits)(raw>>63)<<127;
  if(exponent==2047) encoded|=((bits)32767<<112)|((bits)fraction<<60);
  else if(exponent) encoded|=((bits)(exponent+15360)<<112)|((bits)fraction<<60);
  else if(fraction) {
    int leading=63-__builtin_clzll(fraction);
    encoded|=((bits)(leading-1074+16383)<<112)|
             ((bits)(fraction-(UINT64_C(1)<<leading))<<(112-leading));
  }
  quad result;memcpy(&result,&encoded,sizeof result);return result;
}
static quad adjacent(quad x,int up) {
  bits u;memcpy(&u,&x,sizeof u);
  if(x==0) u=up?1:(((bits)1<<127)|1);
  else if((x>0)==up) ++u;
  else --u;
  memcpy(&x,&u,sizeof x);return x;
}
int kls_accuracy_certify(size_t n,const int64_t *p,const int64_t *rows,const double *a,
                 const double *b,const double *x,int transpose,double *scaled_upper) {
  /* An inexpensive sufficient proof may finish the caller gate, but a miss
   * must not change acceptance or recovery policy. Keep the binary128 path
   * below as the authoritative fallback, including unusual FP modes. */
  if (!scaled_upper) return 0;
  long double fast_upper, fast_budget;
  if (gamma_csc_certificate(n,p,rows,a,b,x,transpose,&fast_upper,&fast_budget) ||
      interval_csc_certificate(n,p,rows,a,b,x,transpose,&fast_upper,&fast_budget)) {
    *scaled_upper=nextafter((double)sqrtl(fast_upper/fast_budget),INFINITY);
    return 1;
  }
  /* Fail closed if the compiler's floating representation does not match
   * the bit-exact conversion below. Constant-folded on supported targets. */
  const quad one=1;bits one_bits;memcpy(&one_bits,&one,sizeof one);
  if(one_bits != ((bits)16383<<112) ||
     (quad)0x1p112 + 1 == (quad)0x1p112 ||
     (quad)0x1p113 + 1 != (quad)0x1p113) return 0;
  if(fegetround()!=FE_TONEAREST || !p || !b || !x || !scaled_upper ||
     n>SIZE_MAX/sizeof(quad)/2 || p[0]!=0 || p[n]<0 || (p[n]&&(!rows||!a))) return 0;
  quad *lo=malloc((n?2*n:2)*sizeof(*lo));if(!lo) return 0;
  quad *hi=lo+n,bn=0;int ok=0,nonzero=0;
  for(size_t j=0;j<n;++j) {
    quad bi=from_double(b[j]),xi=from_double(x[j]);
    if(!__builtin_isfinite(bi)||!__builtin_isfinite(xi)) goto done;
    lo[j]=hi[j]=bi;
    nonzero|=bi!=0;
    /* A product of two finite binary64 values is exact in binary128. */
    bn=adjacent(bn+bi*bi,0);
  }
  for(size_t j=0;j<n;++j) {
    if(p[j]<0 || p[j]>p[j+1] || p[j+1]>p[n]) goto done;
    for(int64_t k=p[j];k<p[j+1];++k) {
      if(rows[k]<0 || (uint64_t)rows[k]>=n) goto done;
      size_t dest=transpose?j:(size_t)rows[k],src=transpose?(size_t)rows[k]:j;
      quad av=from_double(a[k]);if(!__builtin_isfinite(av)) goto done;
      quad product=av*from_double(x[src]);
      lo[dest]=adjacent(lo[dest]-product,0);
      hi[dest]=adjacent(hi[dest]-product,1);
    }
  }
  quad upper=0;
  for(size_t j=0;j<n;++j) {
    quad lower=lo[j]<0?-lo[j]:lo[j],higher=hi[j]<0?-hi[j]:hi[j];
    quad value=lower>higher?lower:higher;
    upper=adjacent(upper+adjacent(value*value,1),1);
  }
  quad tol=adjacent((quad)1/(quad)100000000,0);
  quad budget=adjacent(adjacent(tol*tol,0)*(nonzero?bn:1),0);
  if(!(budget>0) || !__builtin_isfinite(upper)) goto done;
  *scaled_upper=sqrt((double)(upper/budget));
  ok=upper<=budget;
done:
  free(lo);return ok;
}
