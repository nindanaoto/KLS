/* Standalone cache contract: compare every hit/update with a fresh rebuild. */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <fenv.h>
#undef NDEBUG
#include <assert.h>
typedef int64_t UF_long;
enum {KLS_OK=0,KLS_ERR_INVALID_ARGUMENT=-1,KLS_ERR_OUT_OF_MEMORY=-2,
      KLS_ORIENTATION_NORMAL=1,KLS_ORIENTATION_TRANSPOSE=2,
      KLS_ACCURACY_APPLICATION_MANAGED=3};
typedef struct {
  UF_long n,nnz,*col_ptr,*row_idx,*row_perm,*user_col_perm,*input_to_csc;
  int orientation,accuracy_policy;
  void *componentwise_plan;
} kls_solver;
static int resets,allocation,fail;
static void kls_componentwise_plan_reset(void *p) {(void)p;++resets;}
static void *test_malloc(size_t n) {return ++allocation==fail?NULL:malloc(n);}
static void *test_calloc(size_t n,size_t s) {return ++allocation==fail?NULL:calloc(n,s);}
#define malloc test_malloc
#define calloc test_calloc
#include "../src/kls_original_operator.inc"
#undef malloc
#undef calloc
static void compare(kls_solver *s,double *values,private_original_operator *out) {
  private_original_operator fresh={0};
  assert(private_original_operator_build(s,values,out)==KLS_OK);
  assert(private_original_operator_build(s,values,&fresh)==KLS_OK);
  assert(!memcmp(out->p,fresh.p,4*sizeof(int64_t)));
  assert(!memcmp(out->rows,fresh.rows,5*sizeof(int64_t)));
  assert(!memcmp(out->values,fresh.values,5*sizeof(double)));
  private_original_operator_free(&fresh);
}
static void representations(void) {
  const uint64_t patterns[]={0,UINT64_C(0x8000000000000000),1,
    UINT64_C(0x8000000000000001),UINT64_C(0x7ff0000000000000),
    UINT64_C(0xfff0000000000000),UINT64_C(0x7ff8000000000123),
    UINT64_C(0x7ff0000000000123),UINT64_C(0xffffffffffffffff)};
  const int rounds[]={FE_TONEAREST,FE_DOWNWARD,FE_UPWARD,FE_TOWARDZERO};
  for(int orientation=1;orientation<=2;++orientation)for(int flags=0;flags<8;++flags) {
    UF_long p[]={0,2,3,5},r[]={0,2,1,0,2},rp[]={2,0,1},cp[]={1,2,0},map[]={2,0,4,1,3};
    double a[5];
    kls_solver solver={3,5,p,r,flags&1?rp:NULL,flags&2?cp:NULL,
      flags&4?map:NULL,orientation,0,NULL};
    private_original_operator out={0};
    for(size_t offset=0;offset<sizeof(patterns)/sizeof(*patterns);++offset) {
      for(size_t k=0;k<5;++k)memcpy(a+k,patterns+(offset+k)%9,sizeof(*a));
      for(size_t mode=0;mode<4;++mode) {
        assert(!fesetround(rounds[mode]));compare(&solver,a,&out);
        compare(&solver,a,&out);assert(fegetround()==rounds[mode]);
        /* A numeric update may read the snapshot that it replaces. */
        double saved[5];memcpy(saved,out.values,sizeof(saved));
        private_original_operator expected={0};
        assert(private_original_operator_build(&solver,saved,&expected)==KLS_OK);
        assert(private_original_operator_build(&solver,out.values,&out)==KLS_OK);
        assert(!memcmp(out.values,expected.values,sizeof(saved)));
        private_original_operator_free(&expected);
      }
    }
    private_original_operator_free(&out);
  }
  assert(!fesetround(FE_TONEAREST));
}

int main(void) {
  for(int orientation=1;orientation<=2;++orientation) for(int flags=0;flags<8;++flags) {
    UF_long p[]={0,2,3,5},r[]={0,2,1,0,2},rp[]={2,0,1},cp[]={1,2,0},map[]={2,0,4,1,3};
    double a[]={1,-0.0,3,4,5};
    kls_solver s={3,5,p,r,flags&1?rp:NULL,flags&2?cp:NULL,flags&4?map:NULL,orientation,0,NULL};
    private_original_operator out={0};compare(&s,a,&out);
    UF_long *layout=out.layout;int64_t *sources=out.value_source;
    for(int epoch=0;epoch<4;++epoch) {
      a[0]+=1;a[4]*=2;compare(&s,a,&out);
      assert(out.layout==layout && out.value_source==sources);
    }
    double old[5];memcpy(old,out.values,sizeof old);
    allocation=0;fail=1;
    assert(private_original_operator_build(&s,a,&out)==KLS_ERR_OUT_OF_MEMORY);
    assert(!memcmp(old,out.values,sizeof old));fail=0;
    s.accuracy_policy=KLS_ACCURACY_APPLICATION_MANAGED;a[1]=NAN;
    assert(private_original_operator_build(&s,a,&out)==KLS_ERR_INVALID_ARGUMENT);
    assert(!memcmp(old,out.values,sizeof old));a[1]=2;
    r[0]=1;compare(&s,a,&out);
    p[1]=1;compare(&s,a,&out);
    UF_long t=rp[0];rp[0]=rp[1];rp[1]=t;compare(&s,a,&out);
    t=cp[0];cp[0]=cp[1];cp[1]=t;compare(&s,a,&out);
    t=map[0];map[0]=map[1];map[1]=t;compare(&s,a,&out);
    s.row_perm=s.row_perm?NULL:rp;compare(&s,a,&out);
    s.orientation=3-s.orientation;compare(&s,a,&out);
    private_original_operator_free(&out);
    for(int at=1;at<=8;++at) {
      allocation=0;fail=at;
      int status=private_original_operator_build(&s,a,&out);
      assert(status==(at<=6?KLS_ERR_OUT_OF_MEMORY:KLS_OK));
      if(at>6) assert(!out.layout && !out.value_source);
      private_original_operator_free(&out);
    }
    fail=0;
  }
  representations();
  puts("PASS original-operator cache: changes, permutations, orientations, atomic failures");
}
