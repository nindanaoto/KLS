#define _GNU_SOURCE
#define _POSIX_C_SOURCE 200809L
#include <stdatomic.h>
static atomic_uint activations;
#define KLS_CACHED_PLAIN_TEST_EVENT() ((void)atomic_fetch_add(&activations,1))
#include "../src/kls.c"
#undef NDEBUG
#include <assert.h>
#include "kls_cached_plain_reference.inc"
static size_t capture_plain(kls_solver *s,double *out) {
  size_t at=0;trilinos_klu_l_numeric *v=s->numeric;
  memcpy(out+at,v->Udiag,(size_t)s->n*sizeof(double));at+=s->n;
  for(UF_long b=0;b<s->symbolic->nblocks;++b) {
    memcpy(out+at,v->LUbx[b],(size_t)v->LUsize[b]*sizeof(double));at+=v->LUsize[b];
  }
  memcpy(out+at,v->Offx,(size_t)v->Offp[s->n]*sizeof(double));return at+v->Offp[s->n];
}
static unsigned plain_flags(void) {
#if defined(__SSE2__)
  return _mm_getcsr();
#else
  return (unsigned)fetestexcept(FE_ALL_EXCEPT);
#endif
}
static void plain_clear_flags(unsigned control) {
#if defined(__SSE2__)
  _mm_setcsr(control&~63u);
#else
  (void)control;assert(!feclearexcept(FE_ALL_EXCEPT));
#endif
}
int main(void) {
  enum{N=68};int64_t ap[N+1],ai[N*N];double a[N*N];int nz=0;
  for(int c=0;c<N;++c) {
    ap[c]=nz;
    for(int r=0;r<N;++r)if(r/17==c/17 || r==c-17) {ai[nz]=r;a[nz++]=r==c?100.:.001;}
  }ap[N]=nz;
  kls_options o;kls_default_options(&o);o.threads=1;o.scale=0;o.ordering=KLS_ORDERING_NATURAL;o.orientation=KLS_ORIENTATION_NORMAL;
  kls_solver *s=NULL;assert(kls_create(&s)==KLS_OK);
  assert(kls_analyze_csc(s,KLS_INDEX_INT64,N,ap,ai,0,&o)==KLS_OK);
  assert(kls_factor(s,a)==KLS_OK);assert(s->symbolic->nblocks==4);
  assert(kls_build_refactor_map(s));assert(kls_build_refactor_lu_pointer_cache(s));
  assert(!s->refactor_l_packed_values && !s->refactor_u_packed_values);
  assert(kls_build_refactor_user_input_pos32(s));
  size_t capacity=N+(size_t)s->numeric->Offp[N];
  for(UF_long b=0;b<s->symbolic->nblocks;++b)capacity+=s->numeric->LUsize[b];
  double *expected=malloc(capacity*sizeof(double)),*actual=malloc(capacity*sizeof(double));assert(expected&&actual);
  double x[N];kls_parallel_refactor_shared h={0};kls_parallel_refactor_worker w={0};w.shared=&h;w.x=x;
  h.solver=s;h.numeric=s->numeric;h.symbolic=s->symbolic;h.n=s->n;
  h.map_col_ptr=s->refactor_col_ptr;h.map_row_idx=s->refactor_row_idx;h.map_input_pos=s->refactor_input_pos;h.map_block_start=s->refactor_block_start;
  h.col_ptr=s->col_ptr;h.row_idx=s->row_idx;h.nnz=s->nnz;h.scale=0;
  fenv_t saved;assert(!feholdexcept(&saved));unsigned control=plain_flags();
  for(int native=0;native<2;++native)for(int direct=0;direct<2;++direct)for(int epoch=0;epoch<8;++epoch) {
    for(int c=0;c<N;++c)for(int p=ap[c];p<ap[c+1];++p)a[p]=ai[p]==c?100.+epoch:.001*(1+(p+epoch)%5);
    double *values=NULL;assert(prepare_numeric_values(s,a,&values)==KLS_OK);
    h.values=direct?a:values;h.direct_user_values=direct;h.native_short_l=native;
    memset(x,0,sizeof x);w.singular=w.invalid=0;plain_clear_flags(control);
    for(UF_long b=0;b<4;++b)reference_plain_block(&w,b);
    assert(!w.singular&&!w.invalid);unsigned flags=plain_flags();assert(capture_plain(s,expected)==capacity);
    memset(x,0,sizeof x);w.singular=w.invalid=0;plain_clear_flags(control);unsigned before=activations;
    for(UF_long b=0;b<4;++b)kls_refactor_block_dispatch(&w,b);
    assert(activations==before+4 && !w.singular&&!w.invalid);assert(flags==plain_flags());
    assert(capture_plain(s,actual)==capacity && !memcmp(expected,actual,capacity*sizeof(double)));
    s->options.threads=8;assert(!kls_cached_plain_refactor_block(&w,0));s->options.threads=1;
    h.check_pivots=1;assert(!kls_cached_plain_refactor_block(&w,0));h.check_pivots=0;
    h.snode_run_end=(const UF_long *)1;assert(!kls_cached_plain_refactor_block(&w,0));h.snode_run_end=NULL;
  }
  /* Exercise the supported multiworker fallback and its join. The serial
     fast path must remain disabled and the numeric must match the oracle. */
  double *pool_values=NULL;assert(prepare_numeric_values(s,a,&pool_values)==KLS_OK);
  int invalid=0,rejected=0,singular=0,prefix=0;UF_long rank=0,column=0;
  unsigned pool_before=activations;
  s->options.threads=2;
  assert(run_refactor_pool(s,pool_values,2,0,0,&invalid,&rejected,
    NULL,NULL,NULL,NULL,NULL,NULL,&singular,&rank,&column,&prefix));
  assert(!invalid&&!rejected&&!singular && activations==pool_before);
  s->options.threads=1;
  assert(capture_plain(s,actual)==capacity && !memcmp(expected,actual,capacity*sizeof(double)));
  /* Singular exits preserve the generic body's partial-factor semantics. */
  h.direct_user_values=0;h.native_short_l=1;
  for(int halt=0;halt<2;++halt) {
    memset(a,0,(size_t)nz*sizeof(double));double *values=NULL;
    assert(prepare_numeric_values(s,a,&values)==KLS_OK);h.values=values;h.halt_if_singular=halt;
    memset(x,0,sizeof x);w.singular=w.invalid=0;plain_clear_flags(control);
    for(UF_long b=0;b<4;++b)reference_plain_block(&w,b);
    assert(w.singular&&!w.invalid);UF_long expected_rank=w.numerical_rank;
    unsigned flags=plain_flags();assert(capture_plain(s,expected)==capacity);
    memset(x,0,sizeof x);w.singular=w.invalid=0;plain_clear_flags(control);
    for(UF_long b=0;b<4;++b)kls_refactor_block_dispatch(&w,b);
    assert(w.singular&&!w.invalid && w.numerical_rank==expected_rank);
    assert(flags==plain_flags());assert(capture_plain(s,actual)==capacity);
    assert(!memcmp(expected,actual,capacity*sizeof(double)));
  }
  assert(!fesetenv(&saved));free(expected);free(actual);kls_destroy(s);
  fprintf(stderr,"CACHED_PLAIN_EXACT_CASES 34 activation=%u pool_join=1\n",activations);return 0;
}
