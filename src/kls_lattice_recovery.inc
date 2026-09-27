/* Private cold recovery bridge. Fixed bounded search policy is experimental;
 * only the unchanged original-system certificate authorizes success. */
extern int kls_accuracy_lattice_recover(int,const int64_t*,const int64_t*,const double*,
                           const double*,double*,int,int,int,double,double,int*);
extern int kls_accuracy_lattice_schedule(int,const int64_t*,const int64_t*,const double*,
                           const double*,double*,int,int,int,int,double,double,int*);
_Static_assert(KLS_COLD_LATTICE_FIRST_WIDTH > 0 &&
               KLS_COLD_LATTICE_FIRST_WIDTH <= KLS_COLD_LATTICE_LAST_WIDTH &&
               KLS_COLD_LATTICE_LAST_WIDTH <= 64,
               "cold lattice widths exceed supported workspace");
/* One lazy cold-search allowance per outer caller operation. Nested correction
 * solves, old-frame proposals and caller-frame proposals cannot restart it. */
static double private_recovery_deadline(kls_solver *s) {
  if (!(s->private_cold_deadline > 0))
    s->private_cold_deadline = kls_now_seconds() + KLS_COLD_RECOVERY_SECONDS;
  return s->private_cold_deadline;
}
extern int kls_accuracy_certify(size_t,const int64_t*,const int64_t*,
    const double*,const double*,const double*,int,double*);
/* A rejected refactor answer need not imply a bad original system. Try new
 * pivots without mutating the live factors or sacrificing the search seed.
 * The matrix here already represents the requested (possibly transposed)
 * operator. Each indivisible API stage may finish, but no subsequent stage
 * starts after the parent's absolute deadline. Nested cold search is disabled.
 */
static int private_fresh_proposal(kls_solver *parent,size_t n,
    const int64_t *p,const int64_t *rows,const double *a,const double *b,
    double *x,double deadline) {
  if(parent->private_cold_trial || !n || n>SIZE_MAX/sizeof(double) ||
     !(deadline>kls_now_seconds())) return 0;
  double *candidate=malloc(n*sizeof(*candidate));
  kls_solver *trial=NULL;
  int accepted=0;
  if(!candidate || kls_create(&trial)!=KLS_OK) goto done;
  /* The public default must not change an explicit strict recovery search. */
  if(kls_set_accuracy_policy(trial,parent->accuracy_policy)!=KLS_OK) goto done;
  trial->private_cold_trial=1;
  kls_options options=parent->options;
  /* Analyze owns profile loading; the parent's stored options contain no
   * borrowed path. Copy the already loaded numeric tuning after analysis. */
  options.tuning_profile_path=NULL;
  if(!(deadline>kls_now_seconds()) ||
     kls_analyze_csc(trial,KLS_INDEX_INT64,(int64_t)n,p,rows,0,&options)!=KLS_OK)
    goto done;
  trial->tuning=parent->tuning;
  if(!(deadline>kls_now_seconds()) || kls_factor(trial,a)!=KLS_OK) goto done;
  if(!(deadline>kls_now_seconds()) ||
     kls_solve(trial,1,b,0,candidate,0)!=KLS_OK) goto done;
  double bound;
  if(kls_accuracy_certify(n,p,rows,a,b,candidate,0,&bound)) {
    memcpy(x,candidate,n*sizeof(*x));accepted=1;
  }
done:
  kls_destroy(trial);free(candidate);return accepted;
}
static int private_lattice_recover(kls_solver *s,const double *a,const double *b,
                                   double *x,int transpose) {
  if(!s || s->private_cold_trial || s->n<=0 || s->n>INT_MAX-64 || s->nnz<0 || !a || !b || !x ||
     s->row_perm || s->user_col_perm || s->row_scale || s->col_scale ||
     (uintmax_t)s->n>SIZE_MAX/sizeof(int64_t)-1 ||
     (uintmax_t)s->nnz>SIZE_MAX/sizeof(int64_t)) return 0;
  const double deadline = private_recovery_deadline(s);
  if (!(deadline > kls_now_seconds())) return 0;
  const size_t n=(size_t)s->n,nnz=(size_t)s->nnz;
  int64_t *p=calloc(n+1,sizeof(*p)),*next=malloc(n*sizeof(*next));
  int64_t *rows=malloc((nnz?nnz:1)*sizeof(*rows));
  double *values=malloc((nnz?nnz:1)*sizeof(*values));int accepted=0;
  if(!p||!next||!rows||!values) goto done;
  if(s->col_ptr[0]!=0 || s->col_ptr[n]!=s->nnz) goto done;
  for(size_t j=0;j<n;++j) {
    if(s->col_ptr[j]<0 || s->col_ptr[j]>s->col_ptr[j+1] || s->col_ptr[j+1]>s->nnz) goto done;
    for(UF_long k=s->col_ptr[j];k<s->col_ptr[j+1];++k) {
      if(s->row_idx[k]<0 || (uintmax_t)s->row_idx[k]>=n) goto done;
      ++p[(transpose?(size_t)s->row_idx[k]:j)+1];
    }
  }
  for(size_t j=0;j<n;++j) { p[j+1]+=p[j];next[j]=p[j]; }
  for(size_t j=0;j<n;++j) for(UF_long k=s->col_ptr[j];k<s->col_ptr[j+1];++k) {
    size_t col=transpose?(size_t)s->row_idx[k]:j;int64_t at=next[col]++;
    rows[at]=transpose?(int64_t)j:(int64_t)s->row_idx[k];values[at]=a[k];
  }
  {
    int proposals=0;
    accepted=private_fresh_proposal(s,n,p,rows,values,b,x,deadline);
    if(!accepted) accepted=kls_accuracy_lattice_schedule((int)n,p,rows,values,b,x,
      KLS_COLD_LATTICE_FIRST_WIDTH,KLS_COLD_LATTICE_LAST_WIDTH,
      KLS_COLD_LATTICE_SWEEPS,KLS_COLD_LATTICE_ROW_LIMIT,
      KLS_COLD_LATTICE_REGULARIZER,deadline,&proposals)==1;
  }
done:
  free(p);free(next);free(rows);free(values);return accepted;
}
