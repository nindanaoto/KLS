/* Private differential oracle; no test switches enter the installed API. */
#include "../src/kls.c"
#undef NDEBUG
#include <assert.h>
#include <fenv.h>

static void equal_numeric(const kls_solver *r, const kls_solver *f) {
  assert(!memcmp(r->numeric->Udiag, f->numeric->Udiag,
                  (size_t)r->n * sizeof(double)));
  if (r->numeric->Rs)
    assert(!memcmp(r->numeric->Rs, f->numeric->Rs,
                    (size_t)r->n * sizeof(double)));
  for (UF_long k = 0; k < r->n; ++k) for (int upper = 0; upper < 2; ++upper) {
    UF_long *ri, *fi, rn, fn;
    double *rv, *fv;
    kls_klu_get_pointer(r->numeric->LUbx[0],
      upper ? r->numeric->Uip : r->numeric->Lip,
      upper ? r->numeric->Ulen : r->numeric->Llen,k,&ri,&rv,&rn);
    kls_klu_get_pointer(f->numeric->LUbx[0],
      upper ? f->numeric->Uip : f->numeric->Lip,
      upper ? f->numeric->Ulen : f->numeric->Llen,k,&fi,&fv,&fn);
    assert(rn == fn);
    if (rn) {
      assert(!memcmp(ri,fi,(size_t)rn*sizeof(*ri)));
      assert(!memcmp(rv,fv,(size_t)rn*sizeof(*rv)));
    }
  }
}

static kls_solver *make(int n, const int64_t *p, const int64_t *rows,
                        const double *a, int scale, int reorder) {
  kls_options o; kls_default_options(&o);
  o.backend=KLS_BACKEND_SERIAL; o.threads=1; o.use_btf=0;
  o.ordering=reorder?KLS_ORDERING_AMD:KLS_ORDERING_NATURAL;
  o.orientation=KLS_ORIENTATION_NORMAL;
  o.static_pivoting=0; o.scale=scale;
  kls_solver *s=NULL;
  assert(kls_create(&s)==KLS_OK);
  assert(kls_analyze_csc(s,KLS_INDEX_INT64,n,p,rows,0,&o)==KLS_OK);
  assert(kls_factor(s,a)==KLS_OK);
  assert(kls_build_refactor_map(s));
  s->fused_refactor_solve_work=calloc((size_t)n*2,sizeof(double));
  s->fused_refactor_solve_work_n=n;
  assert(s->fused_refactor_solve_work);
  return s;
}

static void exercise(int n, int scale, int star, int rounding, int reorder) {
  int *parent=malloc((size_t)n*sizeof(*parent));
  int64_t *p=malloc(((size_t)n+1)*sizeof(*p));
  int64_t *rows=malloc((size_t)(3*n-2)*sizeof(*rows));
  double *a=malloc((size_t)(3*n-2)*sizeof(*a));
  double *initial=malloc((size_t)(3*n-2)*sizeof(*initial));
  double *rhs=malloc((size_t)n*sizeof(*rhs));
  assert(parent&&p&&rows&&a&&initial&&rhs);
  unsigned seed=15213;
  for(int j=0;j<n-1;++j) {
    seed=1664525u*seed+1013904223u;
    parent[j]=star?n-1:j+1+(int)(seed%(unsigned)(n-j-1));
  }
  parent[n-1]=-1;
  int count=0;
  for(int j=0;j<n;++j) {
    p[j]=count;
    const int oldcol=reorder?(j*17+13)%n:j;
    for(int i=0;i<n;++i) {
      const int oldrow=reorder?(i*17+13)%n:i;
      if(oldrow==oldcol||parent[oldcol]==oldrow||parent[oldrow]==oldcol) {
      rows[count]=i;
      initial[count]=a[count]=i==j?8.0+(oldrow%7):-.125*(1+(oldrow+2*oldcol)%3);
      ++count;
      }
    }
  }
  p[n]=count;assert(count==3*n-2);
  kls_solver *r=make(n,p,rows,a,scale,reorder),*f=make(n,p,rows,a,scale,reorder);
  r->forest_map_state=-1; /* Force the unchanged generic loop as oracle. */
  assert(fesetround(rounding)==0);
  for(int epoch=0;epoch<6;++epoch) {
    for(int t=0;t<count;++t)a[t]=initial[t]*(1+.000125*(epoch+1)*(1+t%5));
    for(int i=0;i<n;++i)rhs[i]=.125*(i%11-5)+epoch*.03125;
    r->fused_refactor_solve_rhs=f->fused_refactor_solve_rhs=rhs;
    feclearexcept(FE_ALL_EXCEPT);
    const int ro=kls_lean_single_block_map32_refactor(r,a,1);
    const int rf=fetestexcept(FE_ALL_EXCEPT);
    feclearexcept(FE_ALL_EXCEPT);
    const int fo=kls_lean_single_block_map32_refactor(f,a,1);
    const int ff=fetestexcept(FE_ALL_EXCEPT);
    assert(ro==1&&fo==ro&&rf==ff);
    assert(f->forest_map_state==1&&f->forest_refactor_count==(UF_long)epoch+1u);
    equal_numeric(r,f);
    assert(!memcmp(r->fused_refactor_solve_work,f->fused_refactor_solve_work,
      (size_t)n*(scale>0?2u:1u)*sizeof(double)));
  }
  /* Early/late singularity and exceptional inputs must preserve the incumbent
     arithmetic, partial factor publication and exception flags as well. */
  for(int halt=0;halt<2;++halt) for(int edge=0;edge<5;++edge) {
    r->common.halt_if_singular=f->common.halt_if_singular=halt;
    memcpy(a,initial,(size_t)count*sizeof(*a));
    const UF_long k=edge==0?0u:(UF_long)n-1u;
    const UF_long col=r->symbolic->Q[k];
    for(int64_t t=p[col];t<p[col+1u];++t) {
      if(edge<2)a[t]=0;
      else if((UF_long)rows[t]==r->numeric->Pnum[k])
        a[t]=edge==2?NAN:edge==3?INFINITY:DBL_TRUE_MIN;
    }
    feclearexcept(FE_ALL_EXCEPT);
    const int ro=kls_lean_single_block_map32_refactor(r,a,1);
    const int rf=fetestexcept(FE_ALL_EXCEPT);
    feclearexcept(FE_ALL_EXCEPT);
    const int fo=kls_lean_single_block_map32_refactor(f,a,1);
    const int ff=fetestexcept(FE_ALL_EXCEPT);
    assert(ro==fo&&rf==ff&&r->common.status==f->common.status);
    assert(r->common.numerical_rank==f->common.numerical_rank);
    assert(r->common.singular_col==f->common.singular_col);
    equal_numeric(r,f);
    assert(!memcmp(r->fused_refactor_solve_work,f->fused_refactor_solve_work,
      (size_t)n*(scale>0?2u:1u)*sizeof(double)));
  }
  /* Map and numeric-owner invalidation must discard the structural plan. */
  kls_invalidate_i32_solve(f);
  assert(f->forest_map_state==0&&!f->forest_input_slot);
  assert(kls_prepare_forest_map(f));
  kls_free_forest_map(f);
  const int32_t saved=f->refactor_row_idx32[0];
  f->refactor_row_idx32[0]=n;
  assert(!kls_prepare_forest_map(f));
  assert(!f->forest_input_slot&&f->forest_map_state<0);
  f->refactor_row_idx32[0]=saved;kls_free_forest_map(f);
  assert(kls_prepare_forest_map(f));
  free_refactor_map(f);
  assert(!f->forest_input_slot&&f->forest_map_state==0);
  assert(fesetround(FE_TONEAREST)==0);
  kls_destroy(r);kls_destroy(f);
  free(parent);free(p);free(rows);free(a);free(initial);free(rhs);
}

static void decline(void) {
  const int64_t p[]={0,3,6,9},rows[]={0,1,2,0,1,2,0,1,2};
  const double a[]={8,1,1,1,8,1,1,1,8};
  kls_solver *s=make(3,p,rows,a,-1,0);
  assert(!kls_prepare_forest_map(s));
  assert(s->forest_map_state<0&&!s->forest_input_slot);
  kls_destroy(s);
}

int main(void) {
  if(!kls_forest_masked_environment()) return 0;
#if defined(__x86_64__)
  const unsigned saved_csr=_mm_getcsr();
  _mm_setcsr((saved_csr & ~0x3fu) & ~(1u<<7));
  assert(!kls_forest_masked_environment());
  _mm_setcsr(saved_csr);
#endif
  for(int scale=-1;scale<=2;++scale) for(int star=0;star<2;++star)
    for(int direction=0;direction<2;++direction)
      for(int reorder=0;reorder<2;++reorder)
        exercise(star?129:97,scale,star,direction?FE_UPWARD:FE_TONEAREST,reorder);
  decline();
  return 0;
}
