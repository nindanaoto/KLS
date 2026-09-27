#define _GNU_SOURCE
#define _POSIX_C_SOURCE 200809L
#include <stdlib.h>
static int allocation_countdown=-1;
static int fail_allocation(void) {
  if (allocation_countdown<0) return 0;
  if (allocation_countdown--==0) return 1;
  return 0;
}
static void *owner_malloc(size_t n) { return fail_allocation() ? NULL : malloc(n); }
static void *owner_calloc(size_t n,size_t size) { return fail_allocation() ? NULL : calloc(n,size); }
#define malloc owner_malloc
#define calloc owner_calloc
#include "../src/kls.c"
#undef malloc
#undef calloc
#undef NDEBUG
#include <assert.h>

static double entry(const kls_solver *s,int row,int col) {
  assert(s->private_original_valid && s->private_original);
  const private_original_operator *a=s->private_original;
  for (int64_t k=a->p[col];k<a->p[col+1];++k)
    if (a->rows[k]==row) return a->values[k];
  return 0;
}

int main(void) {
  int configurations=0;
  int64_t p[]={0,1,2},rows[]={0,1};
  for (int orientation=0;orientation<3;++orientation)
  for (int threads=1;threads<=8;threads+=7)
  for (int scale=-2;scale<=-1;++scale) {
    kls_solver *s=NULL;kls_options o;kls_default_options(&o);
    o.orientation=orientation;o.threads=threads;o.scale=scale;
    o.expected_refactorizations=99;o.expected_solves=100;
    assert(kls_create(&s)==KLS_OK);
    assert(kls_set_accuracy_policy(s,KLS_ACCURACY_STRICT_RHS_L2)==KLS_OK);
    assert(kls_analyze_csc(s,KLS_INDEX_INT64,2,p,rows,0,&o)==KLS_OK);
    double a[]={2,4},b[]={2,4},x[2];
    assert(kls_factor(s,a)==KLS_OK);
    assert(entry(s,0,0)==2 && entry(s,1,1)==4);
    a[0]=200;a[1]=400;
    assert(entry(s,0,0)==2 && entry(s,1,1)==4);
    assert(kls_solve(s,1,b,0,x,0)==KLS_OK);
    assert(x[0]==1 && x[1]==1);
    a[0]=3;a[1]=6;
    assert(kls_refactor(s,a)==KLS_OK);
    assert(entry(s,0,0)==3 && entry(s,1,1)==6);
    assert(kls_refactor(s,a)==KLS_OK); /* Identical-value shortcut. */
    assert(entry(s,0,0)==3 && entry(s,1,1)==6);
    private_original_operator *saved=s->private_original;
    double *saved_values=saved->values;
    /* Model an internal recovery factor: its working operator is NOT A. */
    s->private_solve_depth=1;
    double repaired[]={5,10};
    assert(kls_factor(s,repaired)==KLS_OK);
    s->private_solve_depth=0;
    assert(s->private_original==saved && saved->values==saved_values);
    assert(entry(s,0,0)==3 && entry(s,1,1)==6);
    assert(kls_factor(s,a)==KLS_OK);
    saved=s->private_original;saved_values=saved->values;
    /* Inject each allocation failure at the post-numeric publication boundary. */
    for (int failure=0;failure<6;++failure) {
      s->private_original_valid=0;++s->private_numeric_depth;
      allocation_countdown=failure;
      int status=private_numeric_snapshot_finish(s,a,KLS_OK,1);
      allocation_countdown=-1;
      assert(status==KLS_ERR_OUT_OF_MEMORY && !s->private_original_valid);
      assert(s->private_original==saved && saved->values==saved_values);
      assert(s->private_numeric_depth==0);
      assert(kls_solve(s,1,b,0,x,0)==KLS_ERR_SOLVE_FAILED);
    }
    assert(kls_factor(s,a)==KLS_OK);
    double singular[]={0,0};
    assert(kls_refactor(s,singular)!=KLS_OK);
    assert(!s->private_original_valid && s->private_numeric_depth==0);
    assert(kls_solve(s,1,b,0,x,0)==KLS_ERR_SOLVE_FAILED);
    assert(kls_factor(s,a)==KLS_OK);
    assert(kls_analyze_csc(s,KLS_INDEX_INT64,2,p,rows,0,&o)==KLS_OK);
    assert(!s->private_original && !s->private_original_valid);
    assert(kls_factor(s,a)==KLS_OK);
    kls_destroy(s);++configurations;
  }
  printf("PASS snapshot ownership configurations=%d, six allocation failures each\n",configurations);
  return configurations==12 ? 0 : 1;
}
