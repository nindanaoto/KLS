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
  const int csr=private_original_operator_input_order(a) &&
    a->orientation==KLS_ORIENTATION_TRANSPOSE;
  const int outer=csr?row:col,inner=csr?col:row;
  for (int64_t k=a->p[outer];k<a->p[outer+1];++k)
    if (a->rows[k]==inner) return a->values[k];
  return 0;
}

static void fill_caller_values(int csr,const double *matrix,double *values) {
  static const int64_t ptr[2][4]={{0,3,5,7},{0,2,4,7}};
  static const int64_t idx[2][7]={{0,1,2,1,2,0,2},{0,2,0,1,0,1,2}};
  for(int j=0;j<3;++j) for(int64_t k=ptr[csr][j];k<ptr[csr][j+1];++k)
    values[k]=matrix[csr?3*j+idx[csr][k]:3*idx[csr][k]+j];
}

static void caller_rhs(const double *matrix,int transpose,double *b) {
  const double solution[]={1,-2,3};
  for(int rhs=0;rhs<2;++rhs) {
    for(int row=0;row<3;++row) {
      double value=0;
      for(int col=0;col<3;++col)
        value+=matrix[transpose?3*col+row:3*row+col]*solution[col];
      b[5*rhs+row]=(rhs+1)*value;
    }
    b[5*rhs+3]=77;b[5*rhs+4]=88;
  }
}

static void check_caller_solution(const double *x) {
  const double expected[]={1,-2,3};
  for(int rhs=0;rhs<2;++rhs) {
    for(int row=0;row<3;++row)
      assert(fabs(x[5*rhs+row]-(rhs+1)*expected[row])<2e-12);
    assert(x[5*rhs+3]==77 && x[5*rhs+4]==88);
  }
}

static void caller_pattern_frames(void) {
  static const int64_t ptr[2][4]={{0,3,5,7},{0,2,4,7}};
  static const int64_t idx[2][7]={{0,1,2,1,2,0,2},{0,2,0,1,0,1,2}};
  static const double matrices[2][9]={{6,0,1,2,7,0,3,4,8},
                                           {7,0,1,2,9,0,3,4,11}};
  int configurations=0;
  for(int csr=0;csr<2;++csr) for(int wide=0;wide<2;++wide)
  for(int base=0;base<2;++base) for(int orientation=0;orientation<3;++orientation)
  for(int policy=0;policy<3;++policy) for(int threads=1;threads<=8;threads+=7)
  for(int scale=-2;scale<=-1;++scale) {
    const kls_index_type type=wide?KLS_INDEX_INT64:KLS_INDEX_INT32;
    void *p=malloc(4*(size_t)type),*rows=malloc(7*(size_t)type);
    assert(p && rows);
    for(int j=0;j<4;++j) {
      if(wide) ((int64_t *)p)[j]=ptr[csr][j]+base;
      else ((int32_t *)p)[j]=(int32_t)(ptr[csr][j]+base);
    }
    for(int k=0;k<7;++k) {
      if(wide) ((int64_t *)rows)[k]=idx[csr][k]+base;
      else ((int32_t *)rows)[k]=(int32_t)(idx[csr][k]+base);
    }
    kls_solver *s=NULL;kls_options o;kls_default_options(&o);
    o.orientation=orientation;o.threads=threads;o.scale=scale;
    o.expected_refactorizations=99;o.expected_solves=100;
    assert(kls_create(&s)==KLS_OK);
    assert(kls_set_accuracy_policy(s,(kls_accuracy_policy)policy)==KLS_OK);
    int status=csr?kls_analyze_csr(s,type,3,p,rows,base,&o)
                  :kls_analyze_csc(s,type,3,p,rows,base,&o);
    assert(status==KLS_OK && !s->private_original_valid);
    assert(private_original_operator_input_order(s->private_original));
    assert(!s->private_original->values);
    /* Public analysis owns the input structure, not a borrowed pointer. */
    memset(p,0xff,4*(size_t)type);memset(rows,0xff,7*(size_t)type);
    free(p);free(rows);
    double values[7],b[10];
    for(int update=0;update<2;++update) {
      fill_caller_values(csr,matrices[update],values);
      assert((update?kls_refactor(s,values):kls_factor(s,values))==KLS_OK);
      for(int row=0;row<3;++row) for(int col=0;col<3;++col)
        assert(entry(s,row,col)==matrices[update][3*row+col]);
      for(int k=0;k<7;++k) values[k]=777;
      for(int transpose=0;transpose<2;++transpose) {
        caller_rhs(matrices[update],transpose,b);
        assert((transpose?kls_solve_transpose(s,2,b,5,b,5)
                         :kls_solve(s,2,b,5,b,5))==KLS_OK);
        check_caller_solution(b);
      }
    }
    fill_caller_values(csr,matrices[1],values);
    caller_rhs(matrices[1],0,b);
    assert(kls_refactor_solve(s,values,2,b,5,b,5)==KLS_OK);
    check_caller_solution(b);
    if(csr && policy==KLS_ACCURACY_STRICT_RHS_L2) {
      /* Cold repair must solve the snapshot's operator, not its storage
       * transpose. Change the whole operator independently of the factor. */
      for(size_t k=0;k<s->private_original->nnz;++k)
        s->private_original->values[k]*=2;
      double doubled[9];for(int k=0;k<9;++k) doubled[k]=2*matrices[1][k];
      for(int transpose=0;transpose<2;++transpose) {
        caller_rhs(doubled,transpose,b);
        assert((transpose?kls_solve_transpose(s,2,b,5,b,5)
                         :kls_solve(s,2,b,5,b,5))==KLS_OK);
        check_caller_solution(b);
      }
    }
    kls_destroy(s);++configurations;
  }
  assert(configurations==288);
  printf("PASS caller input frames: %d configurations, both operations, epochs and CSR cold repair\n",configurations);
}

static void caller_pattern_faults(void) {
  int64_t p[]={0,1,2},rows[]={0,1};
  for(int failure=0;failure<3;++failure) {
    allocation_countdown=failure;
    assert(!private_original_input_pattern(KLS_INDEX_INT64,2,p,rows,0,1));
    allocation_countdown=-1;
  }
  int64_t invalid[]={1,1,2};
  assert(!private_original_input_pattern(KLS_INDEX_INT64,2,invalid,rows,0,1));
  kls_solver s={0};s.n=2;
  s.accuracy_policy=KLS_ACCURACY_COMPONENTWISE_BACKWARD_ERROR;
  s.private_original=private_original_input_pattern(KLS_INDEX_INT64,2,p,rows,0,1);
  assert(s.private_original);
  uint64_t bits[]={UINT64_C(0x7ff0000000000017),UINT64_C(0x8000000000000000)};
  double values[2];memcpy(values,bits,sizeof values);
  fenv_t saved_environment;assert(!fegetenv(&saved_environment));
  assert(!feclearexcept(FE_ALL_EXCEPT));assert(!feraiseexcept(FE_DIVBYZERO));
  s.private_numeric_depth=1;
  assert(private_numeric_snapshot_finish(&s,values,KLS_OK,1)==KLS_OK);
  assert(s.private_original_valid && !memcmp(s.private_original->values,bits,sizeof bits));
  assert(fetestexcept(FE_ALL_EXCEPT)==FE_DIVBYZERO);
  double *old=s.private_original->values;
  s.private_original_valid=0;s.private_numeric_depth=1;
  assert(private_numeric_snapshot_finish(&s,old,KLS_OK,1)==KLS_OK);
  assert(!memcmp(s.private_original->values,bits,sizeof bits));
  assert(fetestexcept(FE_ALL_EXCEPT)==FE_DIVBYZERO);
  old=s.private_original->values;
  s.private_original_valid=0;s.private_numeric_depth=1;allocation_countdown=0;
  assert(private_numeric_snapshot_finish(&s,values,KLS_OK,1)==KLS_ERR_OUT_OF_MEMORY);
  allocation_countdown=-1;
  assert(!s.private_original_valid && s.private_numeric_depth==0);
  assert(s.private_original->values==old && !memcmp(old,bits,sizeof bits));
  assert(fetestexcept(FE_ALL_EXCEPT)==FE_DIVBYZERO);
  assert(!fesetenv(&saved_environment));
  s.accuracy_policy=KLS_ACCURACY_APPLICATION_MANAGED;
  double finite[]={2,4};s.private_numeric_depth=1;
  assert(private_numeric_snapshot_finish(&s,finite,KLS_OK,1)==KLS_OK);
  old=s.private_original->values;
  double bad[]={2,INFINITY};s.private_original_valid=0;s.private_numeric_depth=1;
  assert(private_numeric_snapshot_finish(&s,bad,KLS_OK,1)==KLS_ERR_INVALID_ARGUMENT);
  assert(!s.private_original_valid && s.private_numeric_depth==0);
  assert(s.private_original->values==old && !memcmp(old,finite,sizeof finite));
  private_original_operator_free(s.private_original);free(s.private_original);
  puts("PASS caller input capture: optional faults, raw bits/fenv, alias, OOM and finite precommit");
}

int main(void) {
  caller_pattern_faults();caller_pattern_frames();
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
    /* Repeated small-factor preparation must publish one solve cache per
     * numeric. ASan/TSan cover leaked allocations and overlapping builders. */
    for(int epoch=0;epoch<8;++epoch) {
      assert(kls_factor(s,a)==KLS_OK);
      assert(kls_solve(s,1,b,0,x,0)==KLS_OK);
      assert(fabs(x[0]-2.0/3.0)<1e-14 && fabs(x[1]-2.0/3.0)<1e-14);
    }
    saved=s->private_original;saved_values=saved->values;
    /* A cache hit allocates only fresh values. Failure still invalidates the
     * numeric epoch without changing the previously published snapshot. */
    assert(private_original_operator_input_order(saved) ||
           (saved->layout && saved->value_source));
    s->private_original_valid=0;++s->private_numeric_depth;
    allocation_countdown=0;
    int cached_status=private_numeric_snapshot_finish(s,a,KLS_OK,1);
    allocation_countdown=-1;
    assert(cached_status==KLS_ERR_OUT_OF_MEMORY && !s->private_original_valid);
    assert(s->private_original==saved && saved->values==saved_values);
    assert(s->private_numeric_depth==0);
    assert(kls_solve(s,1,b,0,x,0)==KLS_ERR_SOLVE_FAILED);
    /* Disable optional metadata to retain coverage of every required rebuild
     * allocation, rather than assuming that every update rebuilds structure. */
    free(saved->layout);saved->layout=NULL;
    free(saved->value_source);saved->value_source=NULL;
    saved->map_flags&=~PRIVATE_ORIGINAL_INPUT_ORDER;
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
    assert(!s->private_original_valid &&
           (!s->private_original || !s->private_original->values));
    assert(kls_factor(s,a)==KLS_OK);
    kls_destroy(s);++configurations;
  }
  printf("PASS snapshot ownership configurations=%d, cached and six rebuild allocation failures each\n",configurations);
  return configurations==12 ? 0 : 1;
}
