static void native_i16_event(int phase);
#define KLS_NATIVE_I16_TEST_EVENT(phase) native_i16_event(phase)
#include "../src/kls.c"
#include <fenv.h>
#undef NDEBUG
#include <assert.h>
static _Atomic int test_phase,test_rows,test_cancel;
static void native_i16_event(int phase) {
  if(phase==atomic_load(&test_phase)&&atomic_fetch_add(&test_rows,1)>=20)atomic_store(&test_cancel,1);
}
static void check_equal(const double *actual,const double *expected,UF_long n) {
  for(UF_long k=0;k<n;++k) {
    if(isnan(expected[k]))assert(isnan(actual[k]));
    else assert(!memcmp(&actual[k],&expected[k],sizeof(double)));
  }
}
static void fp_checks(kls_solver *s,const kls_native_i16_plan *p,
    const double *rhs,double *actual,double *expected) {
#if defined(__x86_64__) && FLT_EVAL_METHOD==0
  fenv_t original;assert(!fegetenv(&original));unsigned saved_mx;
  __asm__ volatile("stmxcsr %0":"=m"(saved_mx)::"memory");
  const int rounds[]={FE_TONEAREST,FE_UPWARD,FE_DOWNWARD,FE_TOWARDZERO};
  for(int round=0;round<4;++round)for(unsigned mode=0;mode<4;++mode) {
    assert(!fesetround(rounds[round]));unsigned mx;
    __asm__ volatile("stmxcsr %0":"=m"(mx)::"memory");
    mx=(mx|0x1f80u)&~(0x3fu|0x8040u);mx|=(mode&1?1u<<15:0)|(mode&2?1u<<6:0);
    assert(!feclearexcept(FE_ALL_EXCEPT));__asm__ volatile("ldmxcsr %0"::"m"(mx):"memory");
    memcpy(expected,rhs,(size_t)s->n*sizeof(double));assert(kls_native_transpose_solve(s,expected)==1);
    unsigned expected_mx;unsigned short expected_sw;
    __asm__ volatile("stmxcsr %0":"=m"(expected_mx)::"memory");__asm__ volatile("fnstsw %0":"=am"(expected_sw)::"memory");
    assert(!feclearexcept(FE_ALL_EXCEPT));__asm__ volatile("ldmxcsr %0"::"m"(mx):"memory");
    memcpy(actual,rhs,(size_t)s->n*sizeof(double));assert(kls_native_i16_solve(s,p,actual,0,NULL)==1);
    unsigned actual_mx;unsigned short actual_sw;
    __asm__ volatile("stmxcsr %0":"=m"(actual_mx)::"memory");__asm__ volatile("fnstsw %0":"=am"(actual_sw)::"memory");
    assert(actual_mx==expected_mx);assert((actual_sw&0x3fu)==(expected_sw&0x3fu));check_equal(actual,expected,s->n);
  }
  assert(!fesetenv(&original));
  (void)saved_mx;
#else
  (void)s;(void)p;(void)rhs;(void)actual;(void)expected;
#endif
}
static void exercise(int amd,int scale) {
  enum{N=384,NZ=N*N};int64_t *ptr=calloc(N+1,sizeof(*ptr)),*rows=calloc(NZ,sizeof(*rows));
  double *a=calloc(NZ,sizeof(*a)),*truth=calloc(N,sizeof(*truth)),*rhs=calloc(N,sizeof(*rhs));
  double *actual=calloc(N,sizeof(*actual)),*expected=calloc(N,sizeof(*expected));assert(ptr&&rows&&a&&truth&&rhs&&actual&&expected);
  for(int c=0;c<N;++c){ptr[c]=c*N;for(int r=0;r<N;++r){rows[c*N+r]=r;a[c*N+r]=r==c?64:-.0005*(1+(r*13+c*7)%11);}}ptr[N]=NZ;
  kls_options options;kls_default_options(&options);options.threads=1;options.backend=KLS_BACKEND_SERIAL;options.use_btf=0;
  options.ordering=amd?KLS_ORDERING_AMD:KLS_ORDERING_NATURAL;options.orientation=KLS_ORIENTATION_NORMAL;options.scale=scale;
  options.expected_refactorizations=options.expected_solves=0;kls_solver *s=NULL;assert(kls_create(&s)==KLS_OK);
  assert(kls_set_accuracy_policy(s,KLS_ACCURACY_COMPONENTWISE_BACKWARD_ERROR)==KLS_OK);
  assert(kls_analyze_csc(s,KLS_INDEX_INT64,N,ptr,rows,0,&options)==KLS_OK);assert(kls_factor(s,a)==KLS_OK);
  uint64_t total=0;for(int k=0;k<N;++k)total+=s->numeric->Llen[k]+s->numeric->Ulen[k];
  assert(total>2u*UINT16_MAX);
  size_t bytes=sizeof(kls_native_i16_plan)+2u*N*sizeof(kls_native_i16_column)+(size_t)total*sizeof(uint16_t);
  int old_i32=s->i32solve_state;kls_observed_parallel_state old_policy;memcpy(&old_policy,&s->observed_parallel,sizeof(old_policy));
  assert(!kls_native_i16_build(s,0,bytes-1));assert(!kls_native_i16_build(s,1,bytes));
  kls_native_i16_plan *p=kls_native_i16_build(s,0,bytes);assert(p&&p->allocation_bytes==bytes&&p->count==total);
  assert(p->columns[N].offset>UINT16_MAX);assert(s->i32solve_state==old_i32&&!memcmp(&old_policy,&s->observed_parallel,sizeof(old_policy)));
  UF_long old_n=s->numeric->n;s->numeric->n=old_n-1;assert(!kls_native_i16_build(s,0,bytes));s->numeric->n=old_n;
  old_n=s->n;s->n=(UF_long)UINT16_MAX+1;assert(!kls_native_i16_build(s,0,SIZE_MAX));s->n=old_n;
  UF_long *bad=(UF_long *)((double *)s->numeric->LUbx[0]+s->numeric->Uip[N-1]);UF_long saved_index=bad[0];bad[0]=N;
  assert(!kls_native_i16_build(s,0,bytes));bad[0]=saved_index;
  for(int epoch=0;epoch<6;++epoch) {
    memset(rhs,0,N*sizeof(double));
    for(int c=0;c<N;++c){truth[c]=1+.125*(c%7-3)+.015625*(epoch%3);for(int r=0;r<N;++r)a[c*N+r]=r==c?64+.03125*epoch:epoch==0?-DBL_MIN:-.0005*(1+(r*13+c*7+epoch)%11);}
    for(int c=0;c<N;++c)for(int r=0;r<N;++r)rhs[c]+=a[c*N+r]*truth[r];
    double *values=NULL;assert(prepare_numeric_values(s,a,&values)==KLS_OK);assert(kls_serial_klu_refactor_dense_tail(s,values,0)==1);
    memcpy(expected,rhs,N*sizeof(double));assert(kls_native_transpose_solve(s,expected)==1);memcpy(actual,rhs,N*sizeof(double));
    _Atomic int cancel=1;assert(kls_native_i16_solve(s,p,actual,0,&cancel)==-1);assert(!memcmp(actual,rhs,N*sizeof(double)));
    assert(kls_native_i16_solve(s,p,actual,1,NULL)==-1);assert(!memcmp(actual,rhs,N*sizeof(double)));
    for(int phase=1;phase<=2;++phase){atomic_store(&test_phase,phase);atomic_store(&test_rows,0);atomic_store(&test_cancel,0);
      assert(kls_native_i16_solve(s,p,actual,0,&test_cancel)==-1);assert(!memcmp(actual,rhs,N*sizeof(double)));assert(atomic_load(&test_rows)>20);
      atomic_store(&test_phase,0);
    }
    uint64_t now;assert(kls_native_i16_clock(&now));assert(kls_native_i16_solve(s,p,actual,now+10000000000u,NULL)==1);
    check_equal(actual,expected,N);for(int k=0;k<N;++k)assert(fabs(actual[k]-truth[k])<1e-12);
    if(epoch==0)fp_checks(s,p,rhs,actual,expected);
    for(int k=0;k<N;++k)rhs[k]*=2;
    memcpy(expected,rhs,N*sizeof(double));assert(kls_native_transpose_solve(s,expected)==1);memcpy(actual,rhs,N*sizeof(double));
    assert(kls_native_i16_solve(s,p,actual,0,NULL)==1);check_equal(actual,expected,N);
    double normal_rhs[N];memset(normal_rhs,0,sizeof normal_rhs);
    for(int c=0;c<N;++c)for(int row=0;row<N;++row)
      normal_rhs[row]+=a[c*N+row]*truth[c];
    memcpy(expected,normal_rhs,sizeof normal_rhs);
    assert(trilinos_klu_l_solve(s->symbolic,s->numeric,N,1,expected,&s->common));
    memcpy(actual,normal_rhs,sizeof normal_rhs);
    assert(kls_native_i16_normal_solve(s,p,actual,actual,0,NULL)==1);
    check_equal(actual,expected,N);
    memcpy(actual,normal_rhs,sizeof normal_rhs);
    assert(kls_native_i16_normal_solve(s,p,actual,actual,1,NULL)==-1);
    assert(!memcmp(actual,normal_rhs,sizeof normal_rhs));
    _Atomic int normal_cancel=1;
    assert(kls_native_i16_normal_solve(s,p,actual,actual,0,&normal_cancel)==-1);
    assert(!memcmp(actual,normal_rhs,sizeof normal_rhs));
    if(epoch==2){unsigned owner_epoch=p->owner_epoch;assert(kls_parallel_lu_sort(s));assert(s->native_i16_epoch!=owner_epoch);
      memcpy(actual,rhs,N*sizeof(double));assert(kls_native_i16_solve(s,p,actual,0,NULL)==-1);assert(!memcmp(actual,rhs,N*sizeof(double)));
      kls_native_i16_destroy(p);p=kls_native_i16_build(s,0,bytes);assert(p);
    }
  }
  const void *owner=p->numeric_owner;p->numeric_owner=NULL;assert(kls_native_i16_solve(s,p,actual,0,NULL)==-1);p->numeric_owner=owner;
  owner=p->lu_owner;p->lu_owner=NULL;assert(kls_native_i16_solve(s,p,actual,0,NULL)==-1);p->lu_owner=owner;
  double saved_diagonal=((double *)s->numeric->Udiag)[N-1];const uint64_t special[]={UINT64_C(0x7ff0000000000001),UINT64_C(0x7ff0000000000000),UINT64_C(0x7fefffffffffffff)};
  for(int i=0;i<3;++i){memcpy(&((double *)s->numeric->Udiag)[N-1],&special[i],sizeof(double));fp_checks(s,p,rhs,actual,expected);}
  ((double *)s->numeric->Udiag)[N-1]=saved_diagonal;
#if defined(__x86_64__) && FLT_EVAL_METHOD==0
  {unsigned original;assert(!feclearexcept(FE_ALL_EXCEPT));
    __asm__ volatile("stmxcsr %0":"=m"(original)::"memory");
    unsigned unmasked=original&~(1u<<7);__asm__ volatile("ldmxcsr %0"::"m"(unmasked):"memory");
    memcpy(actual,rhs,N*sizeof(double));assert(kls_native_i16_solve(s,p,actual,0,NULL)==-1);
    __asm__ volatile("ldmxcsr %0"::"m"(original):"memory");
    assert(!memcmp(actual,rhs,N*sizeof(double)));
  }
#endif
  free_numeric(s);assert(kls_factor(s,a)==KLS_OK);memcpy(actual,rhs,N*sizeof(double));
  assert(kls_native_i16_solve(s,p,actual,0,NULL)==-1);assert(!memcmp(actual,rhs,N*sizeof(double)));
  fprintf(stderr,"NATIVE_I16 amd=%d scale=%d entries=%llu bytes=%zu exact=1\n",amd,scale,(unsigned long long)total,bytes);
  kls_native_i16_destroy(p);kls_destroy(s);free(ptr);free(rows);free(a);free(truth);free(rhs);free(actual);free(expected);
}
int main(void){exercise(0,0);exercise(1,0);exercise(0,2);exercise(1,2);puts("PASS native uint16 rows/uint32 offsets: fresh matrices/RHS, exact dots, FP, alias/cancel, sort/NUM ownership, budgets");return 0;}
