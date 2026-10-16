/* Full interval proof for a validated, duplicate-free row plan. Stable
 * scatter retains each row's CSC encounter order. Every elementary operation
 * below matches componentwise_certify_body, including its correction RHS.
 * No per-entry mark/coalescing arrays or residual scatter are needed. */
static int componentwise_workspace_overlaps(const double *work,size_t n,
    const void *input,size_t bytes) {
  if(!bytes) return 0;
  const uintptr_t w=(uintptr_t)work,i=(uintptr_t)input;
  const size_t size=2*n*sizeof(*work);
  if(size>UINTPTR_MAX-w || bytes>UINTPTR_MAX-i) return 1;
  return w<i+bytes && i<w+size;
}

/* -2 declines this optional traversal; the original full pass remains the
 * fallback. A finite upper on return zero still promises all n corrections. */
static int componentwise_interval_rows(const kls_componentwise_plan *plan,
    const int64_t *p,const int64_t *rows,const double *a,const double *b,
    const double *x,int transpose,long double *upper,double **correction) {
  *upper=INFINITY;
  if(!correction || !componentwise_plan_environment()) return -2;
  const size_t n=plan->n,nnz=(size_t)plan->nnz;
  /* The CSC pass reads all inputs before publishing corrections. Streaming
   * rows may publish sooner, so retain that pass for overlapping workspaces. */
  if(*correction &&
     (componentwise_workspace_overlaps(*correction,n,b,n*sizeof(*b)) ||
      componentwise_workspace_overlaps(*correction,n,x,n*sizeof(*x)) ||
      componentwise_workspace_overlaps(*correction,n,a,nnz*sizeof(*a)) ||
      componentwise_workspace_overlaps(*correction,n,p,(n+1)*sizeof(*p)) ||
      componentwise_workspace_overlaps(*correction,n,rows,nnz*sizeof(*rows)))) return -2;
  if(!*correction) {
    *correction=malloc(2*n*sizeof(**correction));
    if(!*correction) return -1;
  }
  const long double tolerance=interval_outward(
      1.0L/KLS_COMPONENTWISE_TOLERANCE_DENOMINATOR,-INFINITY);
  long double maximum=0;
  int accepted=1,correction_valid=1;
  for(size_t row=0;row<n;++row) {
    volatile long double bi=(long double)b[row],xi=(long double)x[row];
    if(!isfinite(bi) || !isfinite(xi)) return 0;
    long double lo=bi,hi=bi,den=fabsl(bi);
    const size_t first=transpose?(size_t)p[row]:plan->rowptr[row];
    const size_t last=transpose?(size_t)p[row+1]:plan->rowptr[row+1];
    for(size_t k=first;k<last;++k) {
      const size_t at=transpose?k:(size_t)plan->positions[k];
      const size_t col=transpose?(size_t)rows[k]:(size_t)plan->columns[k];
      volatile long double av=(long double)a[at],xv=(long double)x[col];
      if(!isfinite(av)) return 0;
      if(xv==0 || av==0) continue;
      const long double pl=interval_outward(av*xv,-INFINITY);
      const long double ph=interval_outward(av*xv,INFINITY);
      if(!isfinite(pl) || !isfinite(ph)) return 0;
      lo=interval_outward(lo-ph,-INFINITY);
      hi=interval_outward(hi-pl,INFINITY);
      const long double lower=pl>0?pl:(ph<0?-ph:0);
      if(lower>0) {
        const long double next=interval_outward(den+lower,-INFINITY);
        den=next>0?next:0;
      }
    }
    if(!isfinite(lo) || !isfinite(hi) || !isfinite(den)) return 0;
    const long double left=fabsl(lo),right=fabsl(hi);
    const long double error=left>right?left:right;
    const long double midpoint=lo/2+hi/2;
    if(!isfinite(midpoint) || fabsl(midpoint)>DBL_MAX) correction_valid=0;
    else (*correction)[row]=(double)midpoint;
    if(error==0) continue;
    if(den==0) {maximum=INFINITY;accepted=0;continue;}
    const long double ratio=interval_outward(error/den,INFINITY);
    if(ratio>maximum) maximum=ratio;
    if(error>interval_outward(tolerance*den,-INFINITY)) accepted=0;
  }
  *upper=(!accepted && !correction_valid)?INFINITY:maximum;
  return accepted;
}
