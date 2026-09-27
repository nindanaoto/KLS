/* Conservative componentwise certificate in caller coordinates. Factors and
 * solution updates remain binary64; no binary128 arithmetic is used here. */
#include "kls_componentwise.h"
#include "kls_strict_interval.inc"
#include <limits.h>

int kls_componentwise_certify(size_t n, const int64_t *p,
    const int64_t *rows, const double *a, const double *b, const double *x,
    int transpose, double *residual, long double *upper) {
  if (upper) *upper=INFINITY;
#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__ > 0)
#error Componentwise certification requires strict floating point semantics
#endif
  if (!p || !b || !x || !upper || (transpose!=0 && transpose!=1) ||
      n>SIZE_MAX/(5*sizeof(long double)) || n>SIZE_MAX/(2*sizeof(size_t)) ||
      p[0]!=0 || p[n]<0 || (p[n] && (!a || !rows)) ||
      fegetround()!=FE_TONEAREST || LDBL_MANT_DIG<DBL_MANT_DIG ||
      LDBL_HAS_SUBNORM!=1) return 0;
#if defined(__SSE__)
#if defined(__x86_64__) && LDBL_MANT_DIG == 64 && LDBL_MAX_EXP == 16384
  unsigned short control;
  __asm__ volatile ("fnstcw %0" : "=m" (control));
  if ((control&0x0f00u)!=0x0300u) return 0;
#else
  if (_mm_getcsr()&((1u<<15)|(1u<<6))) return 0;
#endif
#endif
  long double *w=malloc((n?5*n:5)*sizeof(*w));
  size_t *indices=malloc((n?2*n:2)*sizeof(*indices));
  if (!w || !indices) { free(w);free(indices);return -1; }
  long double *lo=w,*hi=lo+n,*den=hi+n,*al=den+n,*ah=al+n;
  size_t *mark=indices,*touched=mark+n;
  int accepted=0;
  for(size_t i=0;i<n;++i) {
    volatile long double bi=(long double)b[i],xi=(long double)x[i];
    if(!isfinite(bi)||!isfinite(xi)) goto done;
    lo[i]=hi[i]=bi;den[i]=fabsl(bi);mark[i]=SIZE_MAX;
  }
  for(size_t j=0;j<n;++j) {
    if(p[j]<0 || p[j]>p[j+1] || p[j+1]>p[n]) goto done;
    size_t count=0;
    /* Coalesce signed duplicates before forming |A|. Taking absolute values
     * of separate duplicate entries could falsely inflate the denominator. */
    for(int64_t k=p[j];k<p[j+1];++k) {
      if(rows[k]<0 || (uintmax_t)rows[k]>=n) goto done;
      const size_t i=(size_t)rows[k];
      volatile long double av=(long double)a[k];
      if(!isfinite(av)) goto done;
      if(mark[i]!=j) { mark[i]=j;touched[count++]=i;al[i]=ah[i]=av; }
      else if(av!=0) {
        al[i]=interval_outward(al[i]+av,-INFINITY);
        ah[i]=interval_outward(ah[i]+av,INFINITY);
      }
    }
    for(size_t k=0;k<count;++k) {
      const size_t i=touched[k],dest=transpose?j:i,src=transpose?i:j;
      volatile long double xv=(long double)x[src];
      if(xv==0 || (al[i]==0 && ah[i]==0)) continue;
      long double pl=interval_outward((xv>0?al[i]:ah[i])*xv,-INFINITY);
      long double ph=interval_outward((xv>0?ah[i]:al[i])*xv,INFINITY);
      if(!isfinite(pl)||!isfinite(ph)) goto done;
      lo[dest]=interval_outward(lo[dest]-ph,-INFINITY);
      hi[dest]=interval_outward(hi[dest]-pl,INFINITY);
      const long double lower=pl>0?pl:(ph<0?-ph:0);
      if(lower>0) {
        const long double next=interval_outward(den[dest]+lower,-INFINITY);
        /* This lower bound is nonnegative. A comparison avoids an out-of-line
         * fmaxl call for every nonzero; NaN still maps to zero as in fmaxl. */
        den[dest]=next>0?next:0;
      }
    }
  }
  {
    const long double tolerance=interval_outward(
        1.0L/KLS_COMPONENTWISE_TOLERANCE_DENOMINATOR,-INFINITY);
    long double maximum=0;accepted=1;
    for(size_t i=0;i<n;++i) {
      if(!isfinite(lo[i])||!isfinite(hi[i])||!isfinite(den[i])) {accepted=0;goto done;}
      /* Both endpoints were checked finite above: no NaN-selection semantics
       * are needed for this maximum. Keep the same outward-rounded bounds. */
      const long double left=fabsl(lo[i]),right=fabsl(hi[i]);
      const long double error=left>right?left:right;
      if(residual) {
        long double midpoint=lo[i]/2+hi[i]/2;
        if(!isfinite(midpoint) || fabsl(midpoint)>DBL_MAX) {accepted=0;goto done;}
        residual[i]=(double)midpoint;
      }
      if(error==0) continue;
      if(den[i]==0) {maximum=INFINITY;accepted=0;continue;}
      const long double ratio=interval_outward(error/den[i],INFINITY);
      if(ratio>maximum) maximum=ratio;
      if(error>interval_outward(tolerance*den[i],-INFINITY)) accepted=0;
    }
    *upper=maximum;
  }
done:
  free(w);free(indices);return accepted;
}
