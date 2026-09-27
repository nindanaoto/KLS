/* Independent, exactly specified cases for the sufficient certificate. */
#undef NDEBUG
#include <assert.h>
#include <fenv.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#if defined(__SSE__)
#include <xmmintrin.h>
#endif
extern int kls_accuracy_certify(size_t,const int64_t*,const int64_t*,
    const double*,const double*,const double*,int,double*);

int main(void) {
  const int rounding=fegetround();assert(fesetround(FE_TONEAREST)==0);
  int64_t p[]={0,1},rows[]={0};
  double a[]={1},b[]={1},x[]={1},bound=INFINITY;
  assert(kls_accuracy_certify(1,p,rows,a,b,x,0,&bound));
  x[0]=1-5e-9;assert(kls_accuracy_certify(1,p,rows,a,b,x,0,&bound));
  x[0]=1-2e-8;assert(!kls_accuracy_certify(1,p,rows,a,b,x,0,&bound));
  b[0]=0;x[0]=5e-9;assert(kls_accuracy_certify(1,p,rows,a,b,x,0,&bound));
  x[0]=2e-8;assert(!kls_accuracy_certify(1,p,rows,a,b,x,0,&bound));
  a[0]=DBL_MAX;b[0]=DBL_MAX;x[0]=1;
  assert(kls_accuracy_certify(1,p,rows,a,b,x,0,&bound));
  x[0]=DBL_MAX;assert(!kls_accuracy_certify(1,p,rows,a,b,x,0,&bound));
  a[0]=DBL_TRUE_MIN;b[0]=DBL_TRUE_MIN;x[0]=1;
  assert(kls_accuracy_certify(1,p,rows,a,b,x,0,&bound));
#if defined(__SSE__)
  unsigned csr=_mm_getcsr();_mm_setcsr(csr|(1u<<15)|(1u<<6));
  assert(kls_accuracy_certify(1,p,rows,a,b,x,0,&bound));
  x[0]=0;assert(!kls_accuracy_certify(1,p,rows,a,b,x,0,&bound));
  _mm_setcsr(csr);
#endif
  a[0]=1;b[0]=1;x[0]=1;
  assert(fesetround(FE_UPWARD)==0);
  assert(!kls_accuracy_certify(1,p,rows,a,b,x,0,&bound));
  assert(fesetround(FE_TONEAREST)==0);
  x[0]=NAN;assert(!kls_accuracy_certify(1,p,rows,a,b,x,0,&bound));
  /* Nonfinite components in empty columns must be rejected explicitly:
   * there is no matrix product to propagate their NaN/Inf into the residual. */
  int64_t empty[]={0,0};
  b[0]=0;
  assert(!kls_accuracy_certify(1,empty,NULL,NULL,b,x,0,&bound));
  x[0]=INFINITY;
  assert(!kls_accuracy_certify(1,empty,NULL,NULL,b,x,0,&bound));
  x[0]=-INFINITY;
  assert(!kls_accuracy_certify(1,empty,NULL,NULL,b,x,1,&bound));
  b[0]=1;
  x[0]=1;a[0]=INFINITY;
  assert(!kls_accuracy_certify(1,p,rows,a,b,x,0,&bound));
  a[0]=1;rows[0]=1;
  assert(!kls_accuracy_certify(1,p,rows,a,b,x,0,&bound));
  rows[0]=0;p[0]=1;
  assert(!kls_accuracy_certify(1,p,rows,a,b,x,0,&bound));
  /* Exact first-row residual is 1, even though a low-precision dot product
   * can lose that unit next to 2^100. Both normal and transpose are checked. */
  int64_t cp[]={0,1,3},ri[]={0,0,1};
  double av[]={1,0x1p100,1},rhs[]={1,1},answer[]={-0x1p100,1};
  assert(!kls_accuracy_certify(2,cp,ri,av,rhs,answer,0,&bound));
  answer[0]=1;answer[1]=-0x1p100;
  assert(!kls_accuracy_certify(2,cp,ri,av,rhs,answer,1,&bound));
  rhs[0]=0x1p100;rhs[1]=1;answer[0]=0;answer[1]=1;
  assert(kls_accuracy_certify(2,cp,ri,av,rhs,answer,0,&bound));
  assert(fesetround(rounding)==0);
  puts("PASS accuracy certificate: cancellation, boundaries, formats, modes");
  return 0;
}
