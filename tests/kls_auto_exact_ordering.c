/* Public current-value solves and symbolic rollback ownership. */
#include "../src/kls.c"
#undef NDEBUG
#include <assert.h>
static void auto_exact_public_checks(void) {
  enum {N=129};
  UF_long ap[N+1],ai[N*N];double a[N*N],b[N],x[N];
  for(int c=0;c<N;++c){ap[c]=c*N;for(int r=0;r<N;++r)ai[c*N+r]=r;}ap[N]=N*N;
  kls_options o;kls_default_options(&o);o.threads=2;o.use_btf=0;o.scale=0;
  o.orientation=KLS_ORIENTATION_TRANSPOSE;o.ordering=KLS_ORDERING_AUTO;
  kls_solver *solver=NULL;assert(kls_create(&solver)==KLS_OK);
  assert(kls_analyze_csc(solver,KLS_INDEX_INT64,N,ap,ai,0,&o)==KLS_OK);
  assert(solver->auto_amf_analysis_seconds>0);
  assert(!kls_amf3_exact_clique_bound && !kls_amf3_exact_selected);
  for(int epoch=0;epoch<3;++epoch){
    memset(b,0,sizeof b);
    for(int c=0;c<N;++c)for(int r=0;r<N;++r){
      a[c*N+r]=r==c?N+1.+epoch:.003*(1+(r+2*c+epoch)%5);
      b[r]+=a[c*N+r]*(1.+.125*((c+epoch)%7));
    }
    if(!epoch)assert(kls_factor(solver,a)==KLS_OK);
    assert(kls_refactor_solve(solver,a,1,b,0,x,0)==KLS_OK);
    for(int c=0;c<N;++c)assert(fabs(x[c]-(1.+.125*((c+epoch)%7)))<1e-11);
  }
  o.ordering=KLS_ORDERING_NATURAL;
  assert(kls_analyze_csc(solver,KLS_INDEX_INT64,N,ap,ai,0,&o)==KLS_OK);
  assert(!solver->auto_amf_exact && !solver->auto_amf_analysis_seconds);
  /* Explicit backend requests must not enter the extra AUTO ordering arm. */
  o.ordering=KLS_ORDERING_AUTO;
  const kls_backend backends[]={KLS_BACKEND_SERIAL,KLS_BACKEND_KLS};
  for(size_t i=0;i<sizeof backends/sizeof backends[0];++i) {
    o.backend=backends[i];o.threads=1;
    assert(kls_analyze_csc(solver,KLS_INDEX_INT64,N,ap,ai,0,&o)==KLS_OK);
    assert(!solver->auto_amf_analysis_seconds);
  }
  kls_destroy(solver);
}

static void auto_exact_fallback_checks(void) {
  enum {N=17};UF_long ap[N+1],ai[N*N];double a[N*N],b[N],x[N];
  for(int c=0;c<N;++c){ap[c]=c*N;for(int r=0;r<N;++r)ai[c*N+r]=r;}ap[N]=N*N;
  for(int scenario=0;scenario<3;++scenario) {
    kls_options o;kls_default_options(&o);o.threads=1;o.use_btf=0;o.scale=0;
    o.orientation=KLS_ORIENTATION_NORMAL;o.ordering=KLS_ORDERING_NATURAL;
    o.fast_factor=0;o.static_pivoting=0;o.halt_if_singular=1;
    kls_solver *s=NULL,*exact=NULL;assert(kls_create(&s)==KLS_OK);
    assert(kls_create(&exact)==KLS_OK);
    assert(kls_analyze_csc(s,KLS_INDEX_INT64,N,ap,ai,0,&o)==KLS_OK);
    trilinos_klu_l_symbolic *ordinary=s->symbolic;
    kls_options eo=o;eo.ordering=KLS_ORDERING_AMF3;kls_amf3_exact_clique_bound=1;
    assert(kls_analyze_csc(exact,KLS_INDEX_INT64,N,ap,ai,0,&eo)==KLS_OK);
    kls_amf3_exact_clique_bound=0;
    kls_pattern_candidate c={0};c.n=s->n;c.nnz=s->nnz;c.col_ptr=s->col_ptr;
    c.row_idx=s->row_idx;c.input_to_csc=s->input_to_csc;c.values=s->values;
    c.orientation=s->orientation;c.symbolic=exact->symbolic;c.common=exact->common;
    c.selected_ordering=KLS_ORDERING_AMF3;c.auto_amf_exact=1;exact->symbolic=NULL;
    kls_auto_amf_capture(&c,&s->symbolic,&s->common,KLS_ORDERING_NATURAL,
      &s->separator,c.symbolic);
    adopt_candidate(s,&c);free_candidate(&c);kls_destroy(exact);
    assert(s->auto_amf_fallback==ordinary && s->auto_amf_owner==s->symbolic);
    assert(kls_factor(s,NULL)==KLS_ERR_INVALID_ARGUMENT);
    assert(s->auto_amf_fallback==ordinary);
    if(scenario==0) {
      memset(a,0,sizeof a);assert(kls_factor(s,a)!=KLS_OK);
      assert(s->symbolic==ordinary && !s->auto_amf_fallback && !s->auto_amf_exact);
    } else if(scenario==1) {
      trilinos_klu_l_symbolic *selected=s->symbolic;
      s->symbolic=ordinary;kls_invalidate_factor_etree_stats(s);
      assert(!s->auto_amf_owner && !s->auto_amf_exact);
      s->symbolic=selected;assert(!kls_auto_amf_restore(s));
      assert(s->auto_amf_fallback==ordinary);
      s->auto_amf_owner=s->symbolic;assert(kls_auto_amf_restore(s));
      assert(s->symbolic==ordinary && !kls_auto_amf_restore(s));
    } else {
      assert(kls_analyze_csc(s,KLS_INDEX_INT64,N,ap,ai,0,&o)==KLS_OK);
      assert(!s->auto_amf_fallback && !s->auto_amf_owner && !s->auto_amf_exact);
    }
    memset(b,0,sizeof b);
    for(int c0=0;c0<N;++c0)for(int r=0;r<N;++r){
      a[c0*N+r]=r==c0?N+1.:.03125;
      b[r]+=a[c0*N+r]*(1.+.125*(c0%7));
    }
    assert(kls_factor(s,a)==KLS_OK);assert(kls_solve(s,1,b,0,x,0)==KLS_OK);
    for(int c0=0;c0<N;++c0)assert(fabs(x[c0]-(1.+.125*(c0%7)))<1e-11);
    kls_destroy(s);
  }
}


int main(void) {
  assert(kls_auto_exact_improves(89,110,0,100,100,0));
  assert(kls_auto_exact_improves(100,79,0,100,100,0));
  assert(kls_auto_exact_improves(101,102,89,100,100,100));
  assert(!kls_auto_exact_improves(101,103,89,100,100,100));
  assert(!kls_auto_exact_improves(100,100,0,100,100,0));
  assert(!kls_auto_exact_improves(NAN,50,1,100,100,100));
  auto_exact_public_checks();
  auto_exact_fallback_checks();
  return 0;
}
