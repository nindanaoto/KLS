#include "data/small_btf_inverter.inc"

/* Deterministic transversal produces a well-conditioned value deck over the
 * unchanged captured pattern. Selection must depend on structure, never A. */
static int small_btf_match(int col,const int64_t *p,const int64_t *rows,
    int *row_match,unsigned char *seen) {
  for(int64_t z=p[col];z<p[col+1];++z) {
    const int row=(int)rows[z];if(seen[row])continue;seen[row]=1;
    if(row_match[row]<0 || small_btf_match(row_match[row],p,rows,row_match,seen)) {
      row_match[row]=col;return 1;
    }
  }return 0;
}
static void small_btf_rhs(const double *a,const double *truth,double *b,int transpose) {
  const int n=SMALL_BTF_FIXTURE_N;memset(b,0,(size_t)n*sizeof(*b));
  for(int col=0;col<n;++col)for(int64_t z=small_btf_fixture_p[col];z<small_btf_fixture_p[col+1];++z) {
    const int row=(int)small_btf_fixture_rows[z];
    const int r=SMALL_BTF_FIXTURE_CSR?col:row,c=SMALL_BTF_FIXTURE_CSR?row:col;
    b[transpose?c:r]+=a[z]*truth[transpose?r:c];
  }
}
static _Atomic int small_btf_fail_numeric_calls;
static void *small_btf_fail_numeric(size_t bytes) {
  (void)bytes;++small_btf_fail_numeric_calls;return NULL;
}
static int small_btf_selection_checks(void) {
  enum{N=SMALL_BTF_FIXTURE_N,NZ=SMALL_BTF_FIXTURE_NNZ};
  double a[NZ],base[NZ],truth[N],b[N],x[N];int matching[N];unsigned char seen[N];
  kls_solver *s=NULL;int result=1;
#define SELECT_REQUIRE(c) do{if(!(c)){fprintf(stderr,"small BTF selection failed line %d\n",__LINE__);goto cleanup;}}while(0)
  SELECT_REQUIRE(kls_small_btf_metis_improves(1410,2150,988,-1));
  SELECT_REQUIRE(!kls_small_btf_metis_improves(974,1887,1002,-1));
  SELECT_REQUIRE(!kls_small_btf_metis_improves(1410,2150,988,NAN));
  SELECT_REQUIRE(!kls_small_btf_metis_improves(1410,2150,988,-2));
  SELECT_REQUIRE(!kls_small_btf_metis_improves(1410,2150,988,2151));
  for(int i=0;i<N;++i)matching[i]=-1;
  for(int col=0;col<N;++col){memset(seen,0,sizeof seen);SELECT_REQUIRE(small_btf_match(col,small_btf_fixture_p,small_btf_fixture_rows,matching,seen));}
  for(int col=0;col<N;++col)for(int64_t z=small_btf_fixture_p[col];z<small_btf_fixture_p[col+1];++z)
    base[z]=matching[small_btf_fixture_rows[z]]==col?8.0:.00390625;
  for(int threads=1;threads<=8;threads+=7)for(int scenario=0;scenario<4;++scenario) {
    kls_options o;kls_default_options(&o);o.threads=threads;
    if(scenario==1){o.fast_factor=0;o.orientation=KLS_ORIENTATION_TRANSPOSE;}
    SELECT_REQUIRE(kls_create(&s)==KLS_OK);
    SELECT_REQUIRE(kls_analyze_csr(s,KLS_INDEX_INT64,N,small_btf_fixture_p,small_btf_fixture_rows,0,&o)==KLS_OK);
    SELECT_REQUIRE(s->options.ordering==KLS_ORDERING_AUTO && s->stats.selected_ordering==KLS_ORDERING_METIS);
    SELECT_REQUIRE(kls_small_btf_owner_matches(s));
    const trilinos_klu_l_symbolic *ordinary=s->small_btf_fallback.symbolic;
    SELECT_REQUIRE(kls_factor(s,NULL)==KLS_ERR_INVALID_ARGUMENT && kls_small_btf_owner_matches(s));
    memcpy(a,base,sizeof a);
    if(scenario==1) {
      small_btf_fail_numeric_calls=0;s->common.malloc_memory=small_btf_fail_numeric;
      SELECT_REQUIRE(kls_factor(s,a)==KLS_OK);
      SELECT_REQUIRE(small_btf_fail_numeric_calls>0 && s->symbolic==ordinary && !s->small_btf_fallback.symbolic);
      SELECT_REQUIRE(s->stats.selected_ordering==KLS_ORDERING_AMD &&
        s->stats.nblocks==ordinary->nblocks && s->stats.max_block==ordinary->maxblock &&
        s->stats.structural_rank==ordinary->structural_rank &&
        s->stats.estimated_flops==ordinary->est_flops);
    } else {
      if(scenario==3) {
        const UF_long old=s->row_idx[0];s->row_idx[0]=(old+1)%N;
        SELECT_REQUIRE(!kls_small_btf_restore(s));s->row_idx[0]=old;
        SELECT_REQUIRE(!s->small_btf_fallback.symbolic);
      }
      SELECT_REQUIRE(kls_factor(s,a)==KLS_OK);
    }
    for(int epoch=1;epoch<=3;++epoch) {
      for(int z=0;z<NZ;++z)a[z]=base[z]*(1.+.03125*epoch);
      for(int k=0;k<N;++k)truth[k]=.125*((k+epoch)%11)-.5;
      small_btf_rhs(a,truth,b,0);
      if(epoch==2){memcpy(x,b,sizeof b);SELECT_REQUIRE(kls_refactor_solve(s,a,1,x,0,x,0)==KLS_OK);}
      else SELECT_REQUIRE(kls_refactor_solve(s,a,1,b,0,x,0)==KLS_OK);
      for(int k=0;k<N;++k)SELECT_REQUIRE(fabs(x[k]-truth[k])<1e-10);
    }
    if(scenario==0 || scenario==2)SELECT_REQUIRE(kls_small_btf_owner_matches(s));
    for(int transpose=0;transpose<2;++transpose) {
      for(int k=0;k<N;++k)truth[k]=.0625*((k+transpose)%13)-.25;
      small_btf_rhs(a,truth,b,transpose);
      SELECT_REQUIRE((transpose?kls_solve_transpose:kls_solve)(s,1,b,0,x,0)==KLS_OK);
      for(int k=0;k<N;++k)SELECT_REQUIRE(fabs(x[k]-truth[k])<1e-10);
    }
    if(scenario==2) {
      memset(a,0,sizeof a);for(int k=0;k<N;++k)x[k]=123.;
      SELECT_REQUIRE(kls_refactor_solve(s,a,1,b,0,x,0)!=KLS_OK);
      SELECT_REQUIRE(s->symbolic==ordinary && !s->small_btf_fallback.symbolic);
      SELECT_REQUIRE(s->stats.selected_ordering==KLS_ORDERING_AMD &&
        s->stats.nblocks==ordinary->nblocks && s->stats.max_block==ordinary->maxblock &&
        s->stats.structural_rank==ordinary->structural_rank &&
        s->stats.estimated_flops==ordinary->est_flops);
      for(int k=0;k<N;++k)SELECT_REQUIRE(x[k]==123.);
    }
    o.ordering=KLS_ORDERING_AMD;
    SELECT_REQUIRE(kls_analyze_csr(s,KLS_INDEX_INT64,N,small_btf_fixture_p,small_btf_fixture_rows,0,&o)==KLS_OK);
    SELECT_REQUIRE(!s->small_btf_fallback.symbolic && s->stats.selected_ordering==KLS_ORDERING_AMD);
    kls_destroy(s);s=NULL;
  }
  /* An optional snapshot allocation failure must retain the original winner. */
  kls_options o;kls_default_options(&o);small_btf_fail_pattern=1;
  SELECT_REQUIRE(kls_create(&s)==KLS_OK);
  SELECT_REQUIRE(kls_analyze_csr(s,KLS_INDEX_INT64,N,small_btf_fixture_p,small_btf_fixture_rows,0,&o)==KLS_OK);
  SELECT_REQUIRE(!small_btf_fail_pattern && !s->small_btf_fallback.symbolic && s->stats.selected_ordering==KLS_ORDERING_AMD);
  kls_destroy(s);s=NULL;
  /* Explicit choices and non-default accuracy must never gain a proposal. */
  const kls_ordering explicit_orders[]={KLS_ORDERING_AMD,KLS_ORDERING_METIS,KLS_ORDERING_NATURAL,KLS_ORDERING_COLAMD};
  for(size_t k=0;k<sizeof(explicit_orders)/sizeof(*explicit_orders);++k) {
    kls_default_options(&o);o.ordering=explicit_orders[k];
    SELECT_REQUIRE(kls_create(&s)==KLS_OK);
    SELECT_REQUIRE(kls_analyze_csr(s,KLS_INDEX_INT64,N,small_btf_fixture_p,small_btf_fixture_rows,0,&o)==KLS_OK);
    SELECT_REQUIRE(!s->small_btf_fallback.symbolic && s->stats.selected_ordering==o.ordering);
    kls_destroy(s);s=NULL;
  }
  const kls_accuracy_policy policies[]={KLS_ACCURACY_STRICT_RHS_L2,KLS_ACCURACY_APPLICATION_MANAGED};
  for(size_t k=0;k<sizeof(policies)/sizeof(*policies);++k) {
    kls_default_options(&o);SELECT_REQUIRE(kls_create(&s)==KLS_OK);
    SELECT_REQUIRE(kls_set_accuracy_policy(s,policies[k])==KLS_OK);
    SELECT_REQUIRE(kls_analyze_csr(s,KLS_INDEX_INT64,N,small_btf_fixture_p,small_btf_fixture_rows,0,&o)==KLS_OK);
    SELECT_REQUIRE(!s->small_btf_fallback.symbolic);kls_destroy(s);s=NULL;
  }
  kls_default_options(&o);
  /* A real grid is a rejection fixture, not an artificial positive case. */
  {
    int64_t gp[14],gi[64];int nz=0;
    for(int col=0;col<13;++col){gp[col]=nz;for(int row=0;row<13;++row)
      if(row==col || (row<9 && col<9 && ((row/3==col/3 && abs(row-col)==1) || abs(row-col)==3)))gi[nz++]=row;
    }gp[13]=nz;
    kls_default_options(&o);o.orientation=KLS_ORIENTATION_TRANSPOSE;
    const unsigned long long before=small_btf_proposals;
    SELECT_REQUIRE(kls_create(&s)==KLS_OK);
    SELECT_REQUIRE(kls_analyze_csc(s,KLS_INDEX_INT64,13,gp,gi,0,&o)==KLS_OK);
    SELECT_REQUIRE(small_btf_proposals>before && !s->small_btf_fallback.symbolic && s->stats.selected_ordering==KLS_ORDERING_AMD);
    kls_destroy(s);s=NULL;kls_default_options(&o);
  }
  /* Actual captured numerical values must also factor and solve successfully. */
  SELECT_REQUIRE(kls_create(&s)==KLS_OK);
  SELECT_REQUIRE(kls_analyze_csr(s,KLS_INDEX_INT64,N,small_btf_fixture_p,small_btf_fixture_rows,0,&o)==KLS_OK);
  SELECT_REQUIRE(kls_factor(s,small_btf_fixture_a)==KLS_OK);
  memset(b,0,sizeof b);SELECT_REQUIRE(kls_solve(s,1,b,0,x,0)==KLS_OK);
  for(int k=0;k<N;++k)SELECT_REQUIRE(x[k]==0.0);
  fprintf(stderr,"SMALL_BTF_PUBLIC_SELECTION_CHECKS 8\n");result=0;
cleanup:kls_destroy(s);
#undef SELECT_REQUIRE
  return result;
}
