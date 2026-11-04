/* Public-call lifetime, aliasing, boundary and failed-output publication. */
static int small_frame_checks(void) {
  fenv_t environment;if(fegetenv(&environment))return 0;
  for(int n=512;n<=513;++n) {
    int *ptr=malloc(((size_t)n+1)*sizeof(int));
    int *idx=malloc(2*(size_t)n*sizeof(int));
    double *a=malloc(2*(size_t)n*sizeof(double));
    double *b=malloc((size_t)n*sizeof(double)),*x=malloc((size_t)n*sizeof(double));
    kls_solver *s=NULL;int ok=ptr&&idx&&a&&b&&x;
    if(!ok){free(ptr);free(idx);free(a);free(b);free(x);(void)fesetenv(&environment);return 0;}
    int nz=0;
    for(int c=0;c<n;++c) {
      ptr[c]=nz;idx[nz]=c;a[nz++]=4.;
      if((c^1)<n){idx[nz]=c^1;a[nz++]=.125;}
    }
    ptr[n]=nz;
    kls_options o;kls_default_options(&o);o.threads=1;o.backend=KLS_BACKEND_SERIAL;
    o.ordering=KLS_ORDERING_NATURAL;o.orientation=KLS_ORIENTATION_TRANSPOSE;
    ok=kls_create(&s)==KLS_OK &&
      kls_analyze_csr(s,KLS_INDEX_INT32,n,ptr,idx,0,&o)==KLS_OK &&
      kls_factor(s,a)==KLS_OK && s->symbolic->nblocks>1;
    for(int epoch=0;epoch<3&&ok;++epoch) {
      for(int c=0;c<n;++c) {
        a[ptr[c]]=4.+.25*epoch;
        b[c]=(4.+.25*epoch+.125*((c^1)<n))*(1.+epoch);
        x[c]=b[c];
      }
      const unsigned long long local=caller_local_frames,heap=caller_heap_frames;
      ok=kls_refactor_solve(s,a,1,x,0,x,0)==KLS_OK;
      ok=ok && (n==512?caller_local_frames==local+1:caller_heap_frames==heap+1);
      for(int c=0;c<n&&ok;++c)ok=fabs(x[c]-(1.+epoch))<1e-12;
    }
    if(ok) {
      for(int c=0;c<n;++c)x[c]=917.;
      ok=kls_refactor_solve(s,a,1,b,n-1,x,n)==KLS_ERR_INVALID_ARGUMENT;
      for(int c=0;c<n&&ok;++c)ok=x[c]==917.;
      a[ptr[0]]=NAN;
      ok=ok && kls_refactor_solve(s,a,1,b,0,x,0)!=KLS_OK;
      for(int c=0;c<n&&ok;++c)ok=x[c]==917.;
    }
    kls_destroy(s);free(ptr);free(idx);free(a);free(b);free(x);
    if(!ok){fprintf(stderr,"small frame check failed n=%d\n",n);(void)fesetenv(&environment);return 0;}
  }
  if(fesetenv(&environment))return 0;
  puts("SMALL_FRAME_ORACLES boundary=512/513 inplace=1 fresh_rhs=1 failure_staging=1");
  return 1;
}
