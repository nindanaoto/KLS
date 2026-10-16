/* Capture the current correction only after a fast proof misses. */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
static int allocations,live,fail_at;
static void *fault_malloc(size_t n) {
  if(++allocations==fail_at) return NULL;
  void *p=malloc(n);if(p) ++live;return p;
}
static void *fault_calloc(size_t n,size_t size) {
  if(++allocations==fail_at) return NULL;
  void *p=calloc(n,size);if(p) ++live;return p;
}
static void fault_free(void *p) {if(p) --live;free(p);}
#define malloc fault_malloc
#define calloc fault_calloc
#define free fault_free
#include "../src/kls_componentwise.c"
#undef malloc
#undef calloc
#undef free
#define CHECK(c) do {if(!(c)) {fprintf(stderr,"line %d: %s\n",__LINE__,#c);return 1;}} while(0)

int main(void) {
  int64_t p[]={0,2,4},rows[]={0,1,0,1};
  double a[]={2,3,1,4},x[]={1,2},b[]={4,11},bt[]={8,9};
  for(int planned=0;planned<2;++planned) {
    kls_componentwise_plan *plan=NULL;double *work=NULL;long double upper;
    for(int transpose=0;transpose<2;++transpose) {
      double rhs[2],reference[2];
      memcpy(rhs,transpose?bt:b,sizeof(rhs));
      CHECK(kls_componentwise_certify_correction(planned?&plan:NULL,2,
          2,p,rows,a,rhs,x,transpose,&upper,NULL,NULL,&work)==1);
#if KLS_COMPONENTWISE_HAS_ACCUMULATED
      if(!transpose) CHECK(!work); /* Fast success allocates no correction. */
#endif
      for(int generation=1;generation<=4;++generation) {
        rhs[transpose]=generation+(transpose?bt[transpose]:b[transpose]);
        long double reference_upper;
        CHECK(kls_componentwise_certify(2,p,rows,a,rhs,x,transpose,
            reference,&reference_upper)==0);
        CHECK(kls_componentwise_certify_correction(planned?&plan:NULL,2,
            2,p,rows,a,rhs,x,transpose,&upper,NULL,NULL,&work)==0);
        CHECK(work && isfinite(upper));
        CHECK(upper==reference_upper && !memcmp(work,reference,sizeof(reference)));
      }
      /* Retained scratch cannot certify a different answer or numeric. */
      x[0]=3;rhs[0]=transpose?12:8;rhs[1]=transpose?11:17;
      CHECK(kls_componentwise_certify_correction(planned?&plan:NULL,2,
          2,p,rows,a,rhs,x,transpose,&upper,NULL,NULL,&work)==1);
      x[0]=1;
    }
    /* An unsupported environment must invalidate the bound even with a
     * retained finite correction from the preceding numeric/RHS. */
    {
      double before[2];memcpy(before,work,sizeof(before));
      const int rounding=fegetround(),before_allocations=allocations;
      CHECK(fesetround(FE_UPWARD)==0);
      int result=kls_componentwise_certify_correction(planned?&plan:NULL,2,
          2,p,rows,a,b,x,0,&upper,NULL,NULL,&work);
      CHECK(fesetround(rounding)==0);
      CHECK(result==0 && !isfinite(upper));
      CHECK(allocations==before_allocations && !memcmp(before,work,sizeof(before)));
    }
    fault_free(work);kls_componentwise_plan_destroy(plan);CHECK(live==0);
  }
  for(fail_at=1;fail_at<=16;++fail_at) {
    allocations=0;kls_componentwise_plan *plan=NULL;double *work=NULL;
    double bad[]={5,11};long double upper;
    int result=kls_componentwise_certify_correction(&plan,2,
        2,p,rows,a,bad,x,0,&upper,NULL,NULL,&work);
    CHECK(result==0 || result==-1);
    if(result==0 && isfinite(upper)) CHECK(work);
    fault_free(work);kls_componentwise_plan_destroy(plan);CHECK(live==0);
  }
  fail_at=0;
  /* Streaming rows must decline overlaps with any caller input. The CSC
   * fallback reads the whole operator before overwriting the correction. */
  for(int input=0;input<3;++input) {
    int64_t ap[]={0,1,2},ar[]={0,1};
    double aa[]={2,3},ab[]={3,7},ax[]={1,2},reference[2];
    double *work=fault_malloc(4*sizeof(*work));CHECK(work);
    const double *av=aa,*bv=ab,*xv=ax;
    if(input==0) {memcpy(work,aa,sizeof(aa));av=work;}
    if(input==1) {memcpy(work,ab,sizeof(ab));bv=work;}
    if(input==2) {memcpy(work,ax,sizeof(ax));xv=work;}
    long double upper,reference_upper;
    CHECK(kls_componentwise_certify(2,ap,ar,av,bv,xv,0,
        reference,&reference_upper)==0);
    kls_componentwise_plan *plan=NULL;
    CHECK(kls_componentwise_certify_correction(&plan,1,2,ap,ar,av,bv,xv,0,
        &upper,NULL,NULL,&work)==0);
    CHECK(upper==reference_upper && !memcmp(work,reference,sizeof(reference)));
    fault_free(work);kls_componentwise_plan_destroy(plan);CHECK(live==0);
  }
#if KLS_COMPONENTWISE_HAS_ACCUMULATED
  {
    int64_t sp[]={0,1,2},sr[]={0,1};
    double sa[]={2,3},sx[]={0x1p-1074,0x1p-1074},sb[]={0,0};
    const unsigned saved_csr=_mm_getcsr();
    _mm_setcsr(saved_csr|(1u<<15)|(1u<<6));
    for(int transpose=0;transpose<2;++transpose) {
      double reference[2],*work=NULL;long double upper,reference_upper;
      kls_componentwise_plan *plan=NULL;
      CHECK(kls_componentwise_certify(2,sp,sr,sa,sb,sx,transpose,
          reference,&reference_upper)==0);
      CHECK(kls_componentwise_certify_correction(&plan,1,2,sp,sr,sa,sb,sx,
          transpose,&upper,NULL,NULL,&work)==0);
      CHECK(upper==reference_upper && work && !memcmp(work,reference,sizeof(reference)));
      CHECK(_mm_getcsr()==(saved_csr|(1u<<15)|(1u<<6)));
      fault_free(work);kls_componentwise_plan_destroy(plan);CHECK(live==0);
    }
    _mm_setcsr(saved_csr);
  }
#endif
  {
    int64_t op[]={0,1,2},orr[]={0,1};
    double oa[]={DBL_MAX,DBL_MAX},ox[]={DBL_MAX,DBL_MAX},ob[]={0,0},*work=NULL;
    kls_componentwise_plan *plan=NULL;long double upper;
    CHECK(kls_componentwise_certify_correction(&plan,1,2,op,orr,oa,ob,ox,0,
        &upper,NULL,NULL,&work)==0);
    CHECK(!isfinite(upper)); /* No binary64 correction can represent this RHS. */
    ox[0]=NAN;
    CHECK(kls_componentwise_certify_correction(&plan,1,2,op,orr,oa,ob,ox,0,
        &upper,NULL,NULL,&work)==0 && !isfinite(upper));
    fault_free(work);kls_componentwise_plan_destroy(plan);CHECK(live==0);
  }
  /* Duplicate coordinates force the interval path. Very large cancelling
   * products may certify even though a correction cannot fit in binary64. */
  {
    int64_t dp[]={0,2,3},dr[]={0,0,0};
    double da[]={0x1p1022,0x1p1022,-nextafter(0x1p1023,0)};
    double dx[]={0x1p1023,0x1p1023},db[]={0,0},*work=NULL;
    long double upper;
    CHECK(kls_componentwise_certify(2,dp,dr,da,db,dx,0,NULL,&upper)==1);
    CHECK(kls_componentwise_certify_correction(NULL,1,
        2,dp,dr,da,db,dx,0,&upper,NULL,NULL,&work)==1);
    da[2]=-0x1p1022;
    CHECK(kls_componentwise_certify_correction(NULL,1,
        2,dp,dr,da,db,dx,0,&upper,NULL,NULL,&work)==0);
    CHECK(!isfinite(upper)); /* No usable correction may be consumed. */
    fault_free(work);CHECK(live==0);
  }
  puts("PASS correction capture: fresh RHS, transpose, reuse, allocation faults, exponent extremes");
  return 0;
}
