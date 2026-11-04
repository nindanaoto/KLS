/* Hoist immutable mapped/unscaled modes and reuse the already-certified LU
 * pointer cache. The generic block remains the fallback for every other mode. */
#if defined(__ELF__) && (defined(__GNUC__) || defined(__clang__))
__attribute__((noinline,aligned(64),section(".text.kls_cached_plain")))
#elif defined(__GNUC__) || defined(__clang__)
__attribute__((noinline,aligned(64)))
#endif
static int kls_cached_plain_refactor_block(kls_parallel_refactor_worker *worker,
                                          UF_long block) {
  kls_parallel_refactor_shared *h=worker->shared;kls_solver *s=h->solver;
  if(s->options.threads!=1 || s->refactor_l_packed_values || s->refactor_u_packed_values ||
     h->scale>0 || h->check_pivots || h->snode_run_end || h->padded_src ||
     !h->map_col_ptr || !h->map_block_start || !s->refactor_row_idx32 ||
     s->refactor_lu_pointer_count!=s->n || !s->refactor_l_indices ||
     !s->refactor_l_values || !s->refactor_u_indices || !s->refactor_u_values ||
     (h->direct_user_values && s->prepared_value_scale))return 0;
  const int32_t *source=h->direct_user_values?s->refactor_input_user_pos32:s->refactor_input_pos32;
  if(!source)return 0;
  const UF_long k1=h->symbolic->R[block],k2=h->symbolic->R[block+1];
  if(k2-k1<2)return 0;
#ifdef KLS_CACHED_PLAIN_TEST_EVENT
  KLS_CACHED_PLAIN_TEST_EVENT();
#endif
  const UF_long *cp=h->map_col_ptr,*bs=h->map_block_start;
  const int32_t *ri=s->refactor_row_idx32;
  const UF_long *ll=h->numeric->Llen,*ul=h->numeric->Ulen;
  const UF_long *offp=h->numeric->Offp;
  double *offx=h->numeric->Offx,*diag=h->numeric->Udiag,*x=worker->x;
  UF_long **li=s->refactor_l_indices,**ui=s->refactor_u_indices;
  double **lv=s->refactor_l_values,**uv=s->refactor_u_values;
  const int native=h->native_short_l || !s->refactor_l_indices32;
  for(UF_long g=k1;g<k2;++g) {
    UF_long po=offp[g];
    for(UF_long p=cp[g];p<bs[g];++p)offx[po++]=h->values[(UF_long)source[p]];
    for(UF_long p=bs[g];p<cp[g+1];++p)x[(UF_long)ri[p]-k1]=h->values[(UF_long)source[p]];
    const UF_long *upper=ui[g];double *uvalues=uv[g];
    for(UF_long p=0;p<ul[g];++p) {
      UF_long j=upper[p];double u=x[j];x[j]=0.;uvalues[p]=u;
      if(u!=0.) {
        UF_long dep=k1+j,len=ll[dep];const UF_long *rows=li[dep];const double *values=lv[dep];
        if(len<=3) {
          if(len>0)x[rows[0]]-=values[0]*u;
          if(len>1)x[rows[1]]-=values[1]*u;
          if(len>2)x[rows[2]]-=values[2]*u;
        } else if(native)kls_scatter_subtract(x,rows,values,len,u);
        else kls_scatter_subtract_refactor_l(s,x,dep,rows,values,len,u);
      }
    }
    const UF_long k=g-k1;double pivot=x[k];x[k]=0.;
    if(pivot==0.) {
      kls_worker_record_singular(worker,g,h->symbolic->Q[g]);
      if(h->halt_if_singular)return 1;
    }
    diag[g]=pivot;
    const UF_long *rows=li[g];double *values=lv[g];const UF_long len=ll[g];
    if(len==1) {values[0]=x[rows[0]]/pivot;x[rows[0]]=0.;}
    else if(len>1) {
      const double inverse=1./pivot;
      for(UF_long p=0;p<len;++p) {UF_long r=rows[p];values[p]=x[r]*inverse;x[r]=0.;}
    }
  }
  return 1;
}

/* Keep the original generic numeric function unchanged. Caller-side dispatch
 * sends only serial, proven mapped modes into the optional compact body. */
static KLS_ALWAYS_INLINE void kls_refactor_block_dispatch(
    kls_parallel_refactor_worker *worker,UF_long block) {
  if(worker->shared->solver->options.threads==1 &&
     kls_cached_plain_refactor_block(worker,block))return;
  kls_parallel_refactor_block(worker,block);
}
