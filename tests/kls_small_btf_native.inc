/* Current input/scaling, off-block values, and public-factor publication. */
static int native_btf_checks(void) {
  enum {N=40};double dense[N*N],a[N*N],b[N],x[N],truth[N];
  int64_t ptr[N+1],idx[N*N];kls_solver *s=NULL,*reference=NULL;int result=1,calls=0;
#define BTF_REQUIRE(condition) do {if(!(condition)){fprintf(stderr,"native BTF check failed line %d\n",__LINE__);goto cleanup;}}while(0)
  BTF_REQUIRE(getenv("KLS_EXPERIMENT_NATIVE_BTF_TRIAL")==NULL);
  for(int scaling=0;scaling<=0;++scaling)for(int csr=0;csr<2;++csr)
    for(int threads=1;threads<=8;threads+=7) {
      memset(dense,0,sizeof dense);
      for(int block=0;block<5;++block) {
        for(int r=0;r<7;++r)for(int c=0;c<7;++c) {
          const int row=(7*(8*block+r)+3)%N,col=(7*(8*block+c)+3)%N;
          dense[row*N+col]=r==c?(r==0?.000001:5.+r):.25*(1+r+c);
        }
        const int singleton=(7*(8*block+7)+3)%N;
        dense[singleton*N+singleton]=7.5+block;
        dense[((7*(8*block+1)+3)%N)*N+singleton]=-.375;
        if(block<4)dense[singleton*N+(7*(8*(block+1)+2)+3)%N]=.1875;
      }
      int nnz=0;
      for(int outer=0;outer<N;++outer) {
        ptr[outer]=nnz;
        for(int inner=0;inner<N;++inner) {
          const int row=csr?outer:inner,col=csr?inner:outer;
          if(dense[row*N+col]!=0){idx[nnz]=inner;a[nnz++]=dense[row*N+col];}
        }
      }
      ptr[N]=nnz;
      kls_options o;kls_default_options(&o);o.threads=threads;o.scale=scaling;
      o.ordering=KLS_ORDERING_NATURAL;o.static_pivoting=0;
      o.orientation=csr?KLS_ORIENTATION_TRANSPOSE:KLS_ORIENTATION_NORMAL;
      const int compact_fixture=csr && kls_small_btf_environment();
      BTF_REQUIRE(kls_create(&s)==KLS_OK);
      BTF_REQUIRE((csr?kls_analyze_csr:kls_analyze_csc)(s,KLS_INDEX_INT64,N,ptr,idx,0,&o)==KLS_OK);
      BTF_REQUIRE(kls_factor(s,a)==KLS_OK);
      BTF_REQUIRE(s->symbolic->nblocks>1 && !s->numeric->Rs);
      int permuted=0;for(int k=0;k<N;++k)permuted+=s->numeric->Pnum[k]!=k;
      BTF_REQUIRE(permuted>0);
      if(compact_fixture) {
        BTF_REQUIRE(kls_create(&reference)==KLS_OK);
        BTF_REQUIRE(kls_analyze_csr(reference,KLS_INDEX_INT64,N,ptr,idx,0,&o)==KLS_OK);
        BTF_REQUIRE(kls_factor(reference,a)==KLS_OK);
        BTF_REQUIRE(kls_build_refactor_map(reference) && kls_build_refactor_map(s));
        BTF_REQUIRE(kls_i32_solve_ready(s) && kls_small_btf_cache_ready(s));
        BTF_REQUIRE(btf_descriptor_checks(s,reference,a));
      }
      const unsigned long long factors_before=compact_btf_factor_calls,solves_before=compact_btf_solve_calls;
      int mapped=0;
      for(int epoch=1;epoch<=3;++epoch) {
        for(int outer=0;outer<N;++outer)for(int64_t p=ptr[outer];p<ptr[outer+1];++p) {
          const int row=csr?outer:(int)idx[p],col=csr?(int)idx[p]:outer;
          a[p]=dense[row*N+col]*(.5+.125*((row+epoch)%7))*(1.+.0625*epoch);
        }
        for(int c=0;c<N;++c)truth[c]=.25*((c+epoch)%9)-1.;
        memset(b,0,sizeof b);
        for(int outer=0;outer<N;++outer)for(int64_t p=ptr[outer];p<ptr[outer+1];++p) {
          const int row=csr?outer:(int)idx[p],col=csr?(int)idx[p]:outer;
          b[row]+=a[p]*truth[col];
        }
        if(compact_fixture) {
          memcpy(reference->numeric->Xwork,s->numeric->Xwork,
            (size_t)s->symbolic->maxblock*sizeof(double));
          BTF_REQUIRE(kls_native_btf_map32_refactor(reference,a)==1);
          BTF_REQUIRE(kls_compact_btf_map32_refactor(s,a)==1);
          BTF_REQUIRE(compact_btf_factors_equal(s,reference));
        }
        const unsigned long long public_solves_before=compact_btf_solve_calls;
        const unsigned long long public_direct_before=direct_btf_solve_calls;
        const unsigned long long public_factors_before=compact_btf_factor_calls;
        const unsigned long long descriptors_before=compact_btf_descriptor_calls;
        BTF_REQUIRE(kls_refactor_solve(s,a,1,b,0,x,0)==KLS_OK);
        if(compact_fixture) {
          BTF_REQUIRE(compact_btf_solve_calls>public_solves_before);
          BTF_REQUIRE(direct_btf_solve_calls>public_direct_before);
          BTF_REQUIRE(compact_btf_factor_calls>public_factors_before);
          BTF_REQUIRE(compact_btf_descriptor_calls>descriptors_before);
        } else {
          BTF_REQUIRE(compact_btf_solve_calls==public_solves_before);
          BTF_REQUIRE(direct_btf_solve_calls==public_direct_before);
          BTF_REQUIRE(compact_btf_factor_calls==public_factors_before);
        }
        if(compact_fixture)BTF_REQUIRE(compact_btf_solve_oracle(s));

        mapped+=s->stats.last_refactor_path==KLS_REFACTOR_PATH_MAPPED;
        if(s->stats.last_refactor_path!=KLS_REFACTOR_PATH_MAPPED)
          fprintf(stderr,"native BTF warmup epoch=%d path=%d unchanged=%d direct=%d\n",epoch,
            (int)s->stats.last_refactor_path,s->unchanged_refactor_state,s->direct_klu_choice);
        if(compact_fixture && epoch>1)BTF_REQUIRE(s->stats.last_refactor_path==KLS_REFACTOR_PATH_MAPPED);
        for(int k=0;k<N;++k)BTF_REQUIRE(fabs(x[k]-truth[k])<1e-10);
        ++calls;
      }
      if(compact_fixture)BTF_REQUIRE(mapped>=2);
      if(compact_fixture) {
        BTF_REQUIRE(compact_btf_factor_calls>=factors_before+2);
        BTF_REQUIRE(compact_btf_solve_calls>solves_before+3);
      } else {
        BTF_REQUIRE(compact_btf_factor_calls==factors_before);
        BTF_REQUIRE(compact_btf_solve_calls==solves_before);
      }

      /* Exercise the prepared mapped entry itself, including factors below
         the incumbent's high-work admission floor. */
#ifdef KLS_SCALED_BTF_MAPPED_TEST
      BTF_REQUIRE(kls_mapped_refactor(s,a,0)==1);
#endif
      for(int transpose=0;transpose<2;++transpose) {
        for(int c=0;c<N;++c)truth[c]=.125*((3*c+transpose)%11)-.5;
        memset(b,0,sizeof b);
        for(int outer=0;outer<N;++outer)for(int64_t p=ptr[outer];p<ptr[outer+1];++p) {
          const int row=csr?outer:(int)idx[p],col=csr?(int)idx[p]:outer;
          b[transpose?col:row]+=a[p]*truth[transpose?row:col];
        }
        BTF_REQUIRE((transpose?kls_solve_transpose:kls_solve)(s,1,b,0,x,0)==KLS_OK);
        for(int k=0;k<N;++k)BTF_REQUIRE(fabs(x[k]-truth[k])<1e-10);
      }
      memset(a,0,(size_t)nnz*sizeof(*a));for(int k=0;k<N;++k)x[k]=123.;
      BTF_REQUIRE(kls_refactor_solve(s,a,1,b,0,x,0)!=KLS_OK);
      for(int k=0;k<N;++k)BTF_REQUIRE(x[k]==123.);
      if(compact_fixture) {
        const int64_t small_ptr[4]={0,1,2,3},small_rows[3]={0,1,2};
        const double small_a[3]={2,3,4},small_b[3]={2,-6,12};double small_x[3];
        const unsigned long long old_factors=compact_btf_factor_calls,old_solves=compact_btf_solve_calls;
        BTF_REQUIRE(kls_analyze_csr(s,KLS_INDEX_INT64,3,small_ptr,small_rows,0,&o)==KLS_OK);
        BTF_REQUIRE(atomic_load_explicit(&s->compact_btf_policy,memory_order_relaxed)==0);
        BTF_REQUIRE(kls_factor(s,small_a)==KLS_OK);
        BTF_REQUIRE(!kls_small_btf_cache_capable(s));
        BTF_REQUIRE(kls_refactor_solve(s,small_a,1,small_b,0,small_x,0)==KLS_OK);
        BTF_REQUIRE(small_x[0]==1 && small_x[1]==-2 && small_x[2]==3);
        BTF_REQUIRE(compact_btf_factor_calls==old_factors && compact_btf_solve_calls==old_solves);
      }
      kls_destroy(s);s=NULL;kls_destroy(reference);reference=NULL;
    }
  fprintf(stderr,"DIRECT_BTF_SOLVE_CALLS %llu\n",direct_btf_solve_calls);
  fprintf(stderr,"COMPACT_BTF_ORACLE_CALLS factor=%llu solve=%llu\n",compact_btf_factor_calls,compact_btf_solve_calls);
  fprintf(stderr,"NATIVE_BTF_CERTIFIED_CALLS %d\n",calls);result=0;
cleanup:
  kls_destroy(reference);kls_destroy(s);unsetenv("KLS_EXPERIMENT_NATIVE_BTF_TRIAL");
#undef BTF_REQUIRE
  return result;
}
