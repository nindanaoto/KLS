static void forest_unscaled_canonical(kls_solver *s,const double *a) {
  const size_t n=(size_t)s->n;
  double *l=calloc(n*n,sizeof(*l)),*u=calloc(n*n,sizeof(*u)),*original=calloc(n*n,sizeof(*original));
  assert(l&&u&&original&&!s->numeric->Rs&&s->symbolic->nblocks==1);
  assert(kls_direct_klu_public_frame_capable(s));
  for(size_t col=0;col<n;++col) {
    l[col*n+col]=1;u[col*n+col]=((double *)s->numeric->Udiag)[col];
    for(UF_long z=s->col_ptr[col];z<s->col_ptr[col+1];++z)original[s->row_idx[z]*n+col]=a[z];
    for(int upper=0;upper<2;++upper) {
      UF_long *rows,length;double *values;
      kls_klu_get_pointer(s->numeric->LUbx[0],upper?s->numeric->Uip:s->numeric->Lip,
        upper?s->numeric->Ulen:s->numeric->Llen,col,&rows,&values,&length);
      for(UF_long z=0;z<length;++z)(upper?u:l)[rows[z]*n+col]=values[z];
    }
  }
  for(size_t row=0;row<n;++row)for(size_t col=0;col<n;++col) {
    long double sum=0;
    for(size_t k=0;k<n;++k)sum+=(long double)l[row*n+k]*u[k*n+col];
    const double expected=original[s->numeric->Pnum[row]*n+s->symbolic->Q[col]];
    assert(fabsl(sum-expected)<=1e-12L*(1+fabs(expected)));
  }
  free(l);free(u);free(original);
}

static void forest_unscaled_public_checks(int n,const int64_t *ptr,
    const int64_t *rows,const double *input,int nnz) {
  double *a=malloc((size_t)nnz*sizeof(*a)),*b=malloc((size_t)n*sizeof(*b));
  double *truth=malloc((size_t)n*sizeof(*truth)),*x=malloc((size_t)n*sizeof(*x));
  UF_long *pivots=malloc((size_t)n*sizeof(*pivots));
  assert(a&&b&&truth&&x&&pivots);
  int calls=0;
  int epochs=3;
#ifdef KLS_FOREST_RESUME_TEST
  epochs=12;
#endif
  for(int threads=1;threads<=8;threads+=7) {
    memcpy(a,input,(size_t)nnz*sizeof(*a));
    kls_solver *s=create(n,(int64_t *)ptr,(int64_t *)rows,a,1,1,2);
    assert(s->numeric->Rs&&s->i32solve_state>0);
#ifdef KLS_FOREST_RESUME_TEST
    /* create() deliberately uses the serial KLU constructor. Model the
       identical native KLU-first ownership tag, as exercise() already does
       for the fused-cycle fixture; real AUTO ownership is checked in Xyce. */
    kls_set_last_factor_path(s,KLS_FACTOR_PATH_KLU_FIRST);
#endif
    memcpy(pivots,s->numeric->Pnum,(size_t)n*sizeof(*pivots));
    s->options.backend=KLS_BACKEND_AUTO;s->options.scale=KLS_SCALE_AUTO;s->options.threads=threads;
    s->observed_parallel.choice=-1;
    for(int epoch=0;epoch<epochs;++epoch) {
      for(int z=0;z<nnz;++z)a[z]=input[z]*(1+.0009765625*(epoch+1)*(1+z%3));
      for(int i=0;i<n;++i)truth[i]=1+.015625*((i+epoch)%9);
      rhs(n,ptr,rows,a,1,0,truth,b);
      double *output=epoch==2?b:x;
      assert(kls_refactor_solve(s,a,1,b,0,output,0)==KLS_OK);
      assert(s->forest_unscale_attempted==1&&s->common.scale==-1&&!s->numeric->Rs);
      assert(!memcmp(pivots,s->numeric->Pnum,(size_t)n*sizeof(*pivots)));
      ++calls;answer(n,output,truth);forest_unscaled_canonical(s,a);
      for(int i=0;i<n;++i)truth[i]=-.5+.03125*((2*i+epoch)%7);
      rhs(n,ptr,rows,a,1,0,truth,b);
      assert(kls_solve(s,1,b,0,x,0)==KLS_OK);answer(n,x,truth);
    }
#ifdef KLS_FOREST_RESUME_TEST
    assert(s->transpose_cycle.cycles>0);
#endif
    memset(a,0,(size_t)nnz*sizeof(*a));for(int i=0;i<n;++i)x[i]=123.;
    assert(kls_refactor_solve(s,a,1,b,0,x,0)!=KLS_OK);
    for(int i=0;i<n;++i)assert(x[i]==123.);
    assert(kls_solve(s,1,b,0,x,0)!=KLS_OK);
    kls_options o=s->options;o.scale=2;
    assert(kls_analyze_csr(s,KLS_INDEX_INT64,n,ptr,rows,0,&o)==KLS_OK);
    assert(!s->forest_unscale_attempted);
    memcpy(a,input,(size_t)nnz*sizeof(*a));
    assert(kls_factor(s,a)==KLS_OK);
    rhs(n,ptr,rows,a,1,0,truth,b);
    assert(kls_refactor_solve(s,a,1,b,0,x,0)==KLS_OK);
    assert(s->common.scale==2&&s->numeric->Rs&&!s->forest_unscale_attempted);
    answer(n,x,truth);kls_destroy(s);
  }
  free(a);free(b);free(truth);free(x);free(pivots);
  fprintf(stderr,"FOREST_UNSCALED_PUBLIC_CALLS %d\n",calls);
}
