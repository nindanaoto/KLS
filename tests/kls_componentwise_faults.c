/* Isolate certificate allocation/environment faults without production hooks. */
#include <stdlib.h>
#include <stdio.h>
#include <fenv.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#if defined(__SSE__)
#include <xmmintrin.h>
#endif
static int allocation, fail_at, live;
static void *fault_malloc(size_t size) {
  if (++allocation==fail_at) return NULL;
  void *p=malloc(size);if(p) ++live;return p;
}
static void fault_free(void *p) { if(p) --live;free(p); }
#define malloc fault_malloc
#define free fault_free
#include "../src/kls_componentwise.c"
#undef malloc
#undef free
#define CHECK(c) do { if(!(c)) {fprintf(stderr,"line %d: %s\n",__LINE__,#c);return 1;} } while(0)
int main(void) {
  int64_t p[]={0,1},rows[]={0};double a[]={1},b[]={1},x[]={1};
  long double upper;double residual=41;
  for(fail_at=1;fail_at<=2;++fail_at) {
    allocation=0;
    CHECK(kls_componentwise_certify(1,p,rows,a,b,x,0,&residual,&upper)==-1);
    CHECK(live==0 && residual==41 && isinf(upper));
  }
  fail_at=0;
  CHECK(kls_componentwise_certify(1,p,rows,a,b,x,0,NULL,&upper)==1);
  CHECK(live==0);
  /* Exercise allocation faults in the fast proof as well as the residual
   * path above, then in the interval fallback after a duplicate-coordinate
   * fast miss. Every allocated scratch block must be released. */
  for(fail_at=1;fail_at<=2;++fail_at) {
    allocation=0;
    CHECK(kls_componentwise_certify(1,p,rows,a,b,x,0,NULL,&upper)==-1);
    CHECK(live==0);
  }
#if KLS_COMPONENTWISE_HAS_ACCUMULATED
  int64_t dp[]={0,2},dr[]={0,0};double da[]={0.5,0.5};
  for(fail_at=3;fail_at<=4;++fail_at) {
    allocation=0;
    CHECK(kls_componentwise_certify(1,dp,dr,da,b,x,0,NULL,&upper)==-1);
    CHECK(live==0);
  }
#endif
  fail_at=0;
  int rounding=fegetround();
  int modes[]={FE_UPWARD,FE_DOWNWARD,FE_TOWARDZERO};
  for(size_t i=0;i<sizeof(modes)/sizeof(modes[0]);++i) {
    CHECK(fesetround(modes[i])==0);
    int result=kls_componentwise_certify(1,p,rows,a,b,x,0,&residual,&upper);
    CHECK(fesetround(rounding)==0);
    CHECK(result==0 && residual==41 && live==0);
  }
#if defined(__SSE__)
  unsigned csr=_mm_getcsr();
  a[0]=b[0]=0x0.0000000000001p-1022;
  _mm_setcsr(csr|(1u<<15)|(1u<<6));
  int result=kls_componentwise_certify(1,p,rows,a,b,x,0,NULL,&upper);
  _mm_setcsr(csr);
#if defined(__x86_64__) && LDBL_MANT_DIG==64 && LDBL_MAX_EXP==16384
  /* x87 input loads and arithmetic intentionally support SSE FTZ/DAZ. */
  CHECK(result==1);
  unsigned short control,reduced;
  __asm__ volatile("fnstcw %0":"=m"(control));
  reduced=(unsigned short)((control&~0x0300u)|0x0200u);
  __asm__ volatile("fldcw %0"::"m"(reduced));
  result=kls_componentwise_certify(1,p,rows,a,b,x,0,NULL,&upper);
  __asm__ volatile("fldcw %0"::"m"(control));
  CHECK(result==0);
#else
  CHECK(result==0);
#endif
#endif
  CHECK(live==0);return 0;
}
