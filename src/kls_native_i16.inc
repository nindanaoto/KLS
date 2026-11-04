static int kls_native_i16_clock(uint64_t *ns) {
  struct timespec t;
  if(clock_gettime(CLOCK_MONOTONIC,&t)||t.tv_sec<0||
     (uint64_t)t.tv_sec>(UINT64_MAX-(uint64_t)t.tv_nsec)/1000000000u)return 0;
  *ns=(uint64_t)t.tv_sec*1000000000u+(uint64_t)t.tv_nsec;return 1;
}
static int kls_native_i16_stop(uint64_t deadline,const _Atomic int *cancel) {
  if(cancel&&atomic_load_explicit(cancel,memory_order_acquire))return 1;
  if(deadline){uint64_t now;if(!kls_native_i16_clock(&now)||now>=deadline)return 1;}
  return 0;
}
static kls_native_i16_plan *kls_native_i16_build(kls_solver *s,uint64_t deadline,size_t limit) {
  if(!s||!s->numeric||!s->symbolic||s->n<=0||s->n>UINT16_MAX||
     s->numeric->n!=s->n||s->symbolic->n!=s->n||s->symbolic->nblocks!=1||
     !s->numeric->LUbx||!s->numeric->LUbx[0]||s->symbolic->R[0]!=0||
     s->symbolic->R[1]!=s->n||kls_native_i16_stop(deadline,NULL))return NULL;
  const UF_long n=s->n;uint64_t count=0;
  for(UF_long k=0;k<n;++k) {
    if((k&63)==0&&kls_native_i16_stop(deadline,NULL))return NULL;
    const UF_long l=s->numeric->Llen[k],u=s->numeric->Ulen[k];
    if(l<0||l>n-k-1||u<0||u>k||count>UINT32_MAX-(uint64_t)l-(uint64_t)u)return NULL;
    count+=(uint64_t)l+(uint64_t)u;
  }
  const size_t cols_bytes=2u*(size_t)n*sizeof(kls_native_i16_column);
  const size_t row_count=count?(size_t)count:1;
  if(cols_bytes>SIZE_MAX-sizeof(kls_native_i16_plan)||
     row_count>(SIZE_MAX-sizeof(kls_native_i16_plan)-cols_bytes)/sizeof(uint16_t))return NULL;
  const size_t bytes=sizeof(kls_native_i16_plan)+cols_bytes+row_count*sizeof(uint16_t);
  if(bytes>limit)return NULL;
  kls_native_i16_plan *p=calloc(1,sizeof(*p));if(!p)return NULL;
  p->n=n;p->numeric_owner=s->numeric;p->lu_owner=s->numeric->LUbx[0];
  p->owner_epoch=s->native_i16_epoch;p->count=(uint32_t)count;p->allocation_bytes=bytes;
  p->columns=malloc(cols_bytes);p->rows=malloc(row_count*sizeof(*p->rows));
  if(!p->columns||!p->rows)goto fail;
  uint32_t cursor=0;const double *lu=s->numeric->LUbx[0];
  for(int upper=0;upper<2;++upper)for(UF_long k=0;k<n;++k) {
    if((k&63)==0&&kls_native_i16_stop(deadline,NULL))goto fail;
    const UF_long len=upper?s->numeric->Ulen[k]:s->numeric->Llen[k];
    const UF_long offset=upper?s->numeric->Uip[k]:s->numeric->Lip[k];
    const UF_long index_units=kls_klu_units_for_indices(len);
    const UF_long capacity=s->numeric->LUsize[0];
    if(offset<0||index_units>capacity||len>capacity-index_units||offset>capacity-index_units-len)goto fail;
    const UF_long *indices=(const UF_long *)(lu+offset);
    kls_native_i16_column *col=&p->columns[(upper?n:0)+k];
    col->values=lu+offset+index_units;col->offset=cursor;col->length=(uint32_t)len;
    for(UF_long j=0;j<len;++j) {
      if((cursor&4095u)==0&&kls_native_i16_stop(deadline,NULL))goto fail;
      const UF_long row=indices[j];
      if(row<0||row>=n||(upper?row>=k:row<=k))goto fail;
      p->rows[cursor++]=(uint16_t)row;
    }
  }
  if(cursor!=p->count||kls_native_i16_stop(deadline,NULL))goto fail;
  return p;
fail:kls_native_i16_destroy(p);return NULL;
}
static int kls_native_i16_environment(void) {
#if (defined(__x86_64__)||defined(__i386__)) && FLT_EVAL_METHOD==0
  unsigned mx;unsigned short cw;
  __asm__ volatile("stmxcsr %0":"=m"(mx)::"memory");
  __asm__ volatile("fnstcw %0":"=m"(cw)::"memory");
  return (mx&0x1f80u)==0x1f80u&&(cw&0x3fu)==0x3fu;
#else
  return 0;
#endif
}
/* Current Numeric values remain authoritative. Work is always seeded from the
 * current RHS. Stop/decline leaves out unchanged so native fallback may retry. */
static int kls_native_i16_solve(kls_solver *s,const kls_native_i16_plan *p,
    double *out,uint64_t deadline,const _Atomic int *cancel) {
  if(!s||!p||!out||!s->numeric||!s->symbolic||s->symbolic->nblocks!=1||
     p->numeric_owner!=s->numeric||!s->numeric->LUbx||p->lu_owner!=s->numeric->LUbx[0]||
     p->n!=s->n||s->numeric->n!=s->n||s->symbolic->n!=s->n||
     p->owner_epoch!=s->native_i16_epoch||!s->numeric->Xwork||
     s->refactor_l_packed_valid||s->refactor_u_packed_valid||
     s->row_refactor_values_ready||s->row_refactor_values_dirty||
     !kls_native_i16_environment()||kls_native_i16_stop(deadline,cancel))return -1;
  double *work=s->numeric->Xwork;const double *diag=s->numeric->Udiag;
  const kls_native_i16_column *lower=p->columns,*upper=lower+p->n;
  for(UF_long k=0;k<p->n;++k)work[k]=out[s->symbolic->Q[k]];
  for(UF_long k=0;k<p->n;++k) {
    if((k&63)==0&&kls_native_i16_stop(deadline,cancel))return -1;
    if(p->n-k>8) {
      const kls_native_i16_column *next=&upper[k+8];
      if(next->length)__builtin_prefetch(next->values,0,2);
    }
    const kls_native_i16_column *col=&upper[k];
    work[k]=kls_lean_i16_row_solve_dot(col->values,p->rows+col->offset,work,0,col->length,work[k])/diag[k];
#ifdef KLS_NATIVE_I16_TEST_EVENT
    KLS_NATIVE_I16_TEST_EVENT(1);
#endif
  }
  for(UF_long k=p->n;k-- >0;) {
    if((k&63)==0&&kls_native_i16_stop(deadline,cancel))return -1;
    if(k>=8) {
      const kls_native_i16_column *next=&lower[k-8];
      if(next->length)__builtin_prefetch(next->values,0,2);
    }
    const kls_native_i16_column *col=&lower[k];
    work[k]=kls_lean_i16_row_solve_dot(col->values,p->rows+col->offset,work,0,col->length,work[k]);
#ifdef KLS_NATIVE_I16_TEST_EVENT
    KLS_NATIVE_I16_TEST_EVENT(2);
#endif
  }
  if(kls_native_i16_stop(deadline,cancel))return -1;
  const double *rs=s->numeric->Rs;
  for(UF_long k=0;k<p->n;++k)out[s->numeric->Pnum[k]]=rs?work[k]/rs[k]:work[k];
  s->common.status=TRILINOS_KLU_OK;return 1;
}

static KLS_ALWAYS_INLINE void kls_native_i16_normal_scatter(
    double *restrict x,const uint16_t *restrict rows,
    const double *restrict values,UF_long len,double alpha) {
#if defined(__FMA__)
  kls_i16_solve_scatter(x,rows,values,len,alpha);
#else
  for(UF_long p=0;p<len;++p)x[rows[p]]-=values[p]*alpha;
#endif
}

static int kls_native_i16_normal_solve(kls_solver *s,const kls_native_i16_plan *p,
    const double *rhs,double *out,uint64_t deadline,const _Atomic int *cancel) {
  if(!s||!p||!rhs||!out||!s->numeric||!s->symbolic||s->symbolic->nblocks!=1||
     p->numeric_owner!=s->numeric||!s->numeric->LUbx||p->lu_owner!=s->numeric->LUbx[0]||
     p->n!=s->n||s->numeric->n!=s->n||s->symbolic->n!=s->n||
     p->owner_epoch!=s->native_i16_epoch||!s->numeric->Xwork||!p->columns||!p->rows||
     !s->numeric->Pnum||!s->symbolic->Q||!s->numeric->Udiag||
     s->refactor_l_packed_valid||s->refactor_u_packed_valid||
     s->row_refactor_values_ready||s->row_refactor_values_dirty||
     !kls_native_i16_environment()||kls_native_i16_stop(deadline,cancel))return -1;
  double *work=s->numeric->Xwork;
  const size_t bytes=(size_t)p->n*sizeof(double);
  const uintptr_t w=(uintptr_t)work,r=(uintptr_t)rhs,o=(uintptr_t)out;
  if(bytes>UINTPTR_MAX-w||bytes>UINTPTR_MAX-r||bytes>UINTPTR_MAX-o||
     (w<r+bytes && r<w+bytes)||(w<o+bytes && o<w+bytes))return -1;
  const double *rs=s->numeric->Rs,*diag=s->numeric->Udiag;
  for(UF_long k=0;k<p->n;++k)work[k]=rs?rhs[s->numeric->Pnum[k]]/rs[k]:rhs[s->numeric->Pnum[k]];
  for(UF_long k=0;k<p->n;++k) {
    if((k&63)==0&&kls_native_i16_stop(deadline,cancel))return -1;
    const kls_native_i16_column *col=p->columns+k;
    kls_native_i16_normal_scatter(work,p->rows+col->offset,col->values,col->length,work[k]);
#ifdef KLS_NATIVE_I16_TEST_EVENT
    KLS_NATIVE_I16_TEST_EVENT(3);
#endif
  }
  for(UF_long k=p->n;k-- >0;) {
    if((k&63)==0&&kls_native_i16_stop(deadline,cancel))return -1;
    const kls_native_i16_column *col=p->columns+p->n+k;
    const double value=work[k]/diag[k];work[k]=value;
    kls_native_i16_normal_scatter(work,p->rows+col->offset,col->values,col->length,value);
#ifdef KLS_NATIVE_I16_TEST_EVENT
    KLS_NATIVE_I16_TEST_EVENT(4);
#endif
  }
  if(kls_native_i16_stop(deadline,cancel))return -1;
  for(UF_long k=0;k<p->n;++k)out[s->symbolic->Q[k]]=work[k];
  s->common.status=TRILINOS_KLU_OK;return 1;
}


/* The public combined-call governor owns admission. A decline leaves RHS/out
 * intact, so the unchanged native walk supplies the current caller answer. */
static int kls_native_i16_try_solve(kls_solver *s,double *out) {
  kls_observed_parallel_state *g=&s->native_i16_policy;
  if(!s->native_i16 || s->native_i16_combined_depth!=1 ||
     !(g->active==1 || (g->choice==2 && g->active!=2)))return 0;
  uint64_t deadline=0;
  if(g->active==1) {
    const double left=g->deadline-kls_now_seconds();uint64_t now;
    if(left>0 && left<(double)(UINT64_MAX/1000000000u) && kls_native_i16_clock(&now)) {
      const uint64_t span=(uint64_t)(left*1e9);
      if(span && span<=UINT64_MAX-now)deadline=now+span;
    }
  }
  const int result=g->active==1&&!deadline?-1:kls_native_i16_solve(s,s->native_i16,out,deadline,NULL);
  if(result==1){++s->native_i16_runs;return 1;}
  g->choice=-1;
  const double cleanup_start=kls_now_seconds();kls_native_i16_drop(s);
  g->spent+=kls_now_seconds()-cleanup_start;
  return 0;
}
