/* Check outward enclosure at binade boundaries and representational extremes. */
#undef NDEBUG
#include <assert.h>
#include <stdio.h>
#include "../src/kls_strict_interval.inc"
#include "../src/kls_gamma_certificate.inc"

static void check(long double x) {
  if(!isfinite(x)) return;
  assert(interval_outward(x,INFINITY)>=nextafterl(x,INFINITY));
  assert(interval_outward(x,-INFINITY)<=nextafterl(x,-INFINITY));
}
int main(void) {
  int rounding=fegetround();assert(fesetround(FE_TONEAREST)==0);
  check(0);check(-0.L);check(LDBL_MAX);check(-LDBL_MAX);
  check(nextafterl(0,1));check(nextafterl(0,-1));
  for(int e=LDBL_MIN_EXP-LDBL_MANT_DIG;e<LDBL_MAX_EXP;++e) {
    long double x=ldexpl(1.L,e);
    for(int sign=-1;sign<=1;sign+=2) {
      long double v=sign*x;
      check(v);check(nextafterl(v,INFINITY));check(nextafterl(v,-INFINITY));
      check(v*1.5L);check(v*1.999L);
    }
  }
  uint64_t state=19;
  for(int i=0;i<100000;++i) {
    state=state*UINT64_C(6364136223846793005)+1;
    int e=(int)(state%32000)-16000;
    check(ldexpl((long double)(state|1)/18446744073709551616.L,e));
  }
  int64_t p[]={0,1},rows[]={0};double a[]={1},b[]={1},x[]={1};
  assert(interval_csc_certificate(1,p,rows,a,b,x,0,NULL,NULL));
  x[0]=1-2e-8;
  assert(!interval_csc_certificate(1,p,rows,a,b,x,0,NULL,NULL));
  assert(fesetround(FE_UPWARD)==0);
  x[0]=1;assert(!interval_csc_certificate(1,p,rows,a,b,x,0,NULL,NULL));
  assert(fesetround(rounding)==0);
  /* Both orientations must reject nonfinite x even in empty columns. */
  assert(fesetround(FE_TONEAREST)==0);
  int64_t ep[]={0,1,1}, er[]={0};
  double ea[]={1}, eb[]={1,0}, ex[]={1,0};
  for(int t=0;t<2;++t) {
    for(int j=0;j<2;++j) {
      for(int k=0;k<3;++k) {
        ex[0]=1;ex[1]=0;ex[j]=k==0?NAN:k==1?INFINITY:-INFINITY;
        assert(!interval_csc_certificate(2,ep,er,ea,eb,ex,t,NULL,NULL));
        assert(!gamma_csc_certificate(2,ep,er,ea,eb,ex,t,NULL,NULL));
      }
    }
  }
  assert(fesetround(rounding)==0);
  puts("PASS outward interval bounds and rounding-mode guard");
}
