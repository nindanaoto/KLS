/* Independent leaves of a proved fixed-pattern forest have no U dependency.
 * Only structural slots persist; every value and RHS is recomputed per call. */
#include <fenv.h>

#include "kls_forest_remaining_plan.inc"

#include "kls_forest_staged_leaves.inc"

static int kls_prepare_forest_leaves(kls_solver *s) {
  if (s->forest_leaf_state) return s->forest_leaf_state > 0;
  if (s->forest_map_state <= 0 || !s->forest_input_slot || !s->numeric ||
      !s->symbolic || s->symbolic->nblocks != 1 || s->n < 64 ||
      s->n > INT32_MAX || s->nnz > INT32_MAX) return 0;
  s->forest_leaf_state = -1;
  const UF_long n = s->n;
  UF_long count = 0;
  for (UF_long k=0;k<n;++k)
    count += s->numeric->Ulen[k]==0 && s->numeric->Llen[k]==1;
  if (count < 32 || count < n-count) return 0;
#ifdef __linux__
  const uint64_t llc = kls_execution_domain_llc_bytes();
#else
  const uint64_t llc = 0;
#endif
  const uintmax_t index_count=(uintmax_t)n+5*(uintmax_t)count;
  const long double bytes = sizeof(kls_forest_leaf_plan) +
    (long double)index_count*sizeof(int32_t) +
    (long double)s->nnz*sizeof(int32_t); /* includes the retained forest proof */
  if (!llc || bytes > llc/32 ||
      index_count > SIZE_MAX/sizeof(int32_t)) return 0;
  kls_forest_leaf_plan *p=calloc(1,sizeof(*p));
  if (!p) return 0;
  p->n=n;p->count=count;p->remaining=n-count;
  p->indices=malloc((size_t)index_count*sizeof(*p->indices));
  if (!p->indices) goto decline;
  int32_t *ks=p->indices,*parents=ks+count,*ds=parents+count;
  int32_t *ls=ds+count,*qs=ls+count,*offset=qs+count,*rest=offset+count;
  UF_long at=0,other=0;
  for (UF_long k=0;k<n;++k) {
    if (s->numeric->Ulen[k]!=0 || s->numeric->Llen[k]!=1) {
      rest[other++]=(int32_t)k;continue;
    }
    UF_long *li,length;double *lx;
    kls_klu_get_pointer(s->numeric->LUbx[0],s->numeric->Lip,
      s->numeric->Llen,k,&li,&lx,&length);
    if (length!=1 || li[0]<=k || li[0]>=n ||
        s->numeric->Ulen[li[0]]==0 || s->symbolic->Q[k]<0 ||
        s->symbolic->Q[k]>=n) goto decline;
    int32_t diagonal=-1,lower=-1;
    for (UF_long j=s->refactor_col_ptr[k];j<s->refactor_col_ptr[k+1];++j) {
      const int32_t source=s->refactor_input_pos32[j];
      if (source<0 || (UF_long)source>=s->nnz) goto decline;
      if (s->forest_input_slot[j]==-1 && diagonal<0 &&
          s->refactor_row_idx32[j]==(int32_t)k) diagonal=source;
      else if (s->forest_input_slot[j]==-2 && lower<0 &&
          s->refactor_row_idx32[j]==(int32_t)li[0]) lower=source;
      else goto decline;
    }
    if (diagonal<0 || lower<0) goto decline;
    ks[at]=(int32_t)k;parents[at]=(int32_t)li[0];
    ds[at]=diagonal;ls[at]=lower;qs[at]=(int32_t)s->symbolic->Q[k];
    const UF_long units=kls_klu_units_for_indices(1);
    if(s->numeric->Lip[k]<0 || s->numeric->Lip[k]>INT32_MAX-units) goto decline;
    offset[at++]=(int32_t)(s->numeric->Lip[k]+units);
  }
  if (at!=count || other!=p->remaining) goto decline;
  p->numeric_owner=s->numeric;p->lu_owner=s->numeric->LUbx[0];
  kls_prepare_forest_remaining_plan(s,p,llc);
  s->forest_leaves=p;s->forest_leaf_state=1;
  return 1;
decline:
  free(p->indices);free(p);return 0;
}

/* Restrict only the optional accelerator. These exponent bounds prove a
 * normal finite reciprocal and finite lower/RHS products, so publication
 * needs neither speculative native writes nor a full temporary value deck. */
static int kls_forest_leaf_validate_scalar(const kls_forest_leaf_plan *p,
    const double *values,const double *rhs,const double *scale,UF_long first) {
  const UF_long count=p->count;
  const int32_t *ks=p->indices,*parents=ks+count,*ds=parents+count;
  const int32_t *ls=ds+count,*qs=ls+count;
  for(UF_long i=first;i<count;++i) {
    const double d=scale?values[ds[i]]*scale[ks[i]]:values[ds[i]];
    const double l=scale?values[ls[i]]*scale[parents[i]]:values[ls[i]];
    if(!(fabs(d)>=0x1p-256 && fabs(d)<=0x1p256 &&
         fabs(l)<=0x1p256 && fabs(rhs[qs[i]])<=0x1p256)) return 0;
  }
  return 1;
}
static void kls_forest_leaf_publish_scalar(kls_solver *s,
    const kls_forest_leaf_plan *p,const double *values,const double *rhs,
    const double *scale,double *forward,double *diag_recip,UF_long first) {
  const UF_long count=p->count;
  const int32_t *ks=p->indices,*parents=ks+count,*ds=parents+count;
  const int32_t *ls=ds+count,*qs=ls+count,*offset=qs+count;
  double *lu=(double *)s->numeric->LUbx[0],*diagonal=s->numeric->Udiag;
  for(UF_long i=first;i<count;++i) {
    const UF_long k=(UF_long)ks[i];
    const double d=scale?values[ds[i]]*scale[k]:values[ds[i]];
    const double l=scale?values[ls[i]]*scale[parents[i]]:values[ls[i]];
    const double inverse=1.0/d;
    diagonal[k]=d;lu[offset[i]]=l*inverse;forward[k]=rhs[qs[i]]*inverse;
    if(diag_recip) diag_recip[k]=inverse;
  }
}
#if KLS_HAVE_AVX512_KERNELS
__attribute__((target("avx512f"),noinline))
static int kls_forest_leaf_validate_avx512(const kls_forest_leaf_plan *p,
    const double *values,const double *rhs,const double *scale) {
  const UF_long count=p->count;
  const int32_t *ks=p->indices,*parents=ks+count,*ds=parents+count;
  const int32_t *ls=ds+count,*qs=ls+count;
  const __m512i magnitude=_mm512_set1_epi64(INT64_MAX);
  const __m512d small=_mm512_set1_pd(0x1p-256),large=_mm512_set1_pd(0x1p256);
  UF_long i=0;
  for(;i+8<=count;i+=8) {
    __m512d d=_mm512_i32gather_pd(_mm256_loadu_si256((const __m256i *)(ds+i)),values,8);
    __m512d l=_mm512_i32gather_pd(_mm256_loadu_si256((const __m256i *)(ls+i)),values,8);
    if(scale) {
      d=_mm512_mul_pd(d,_mm512_i32gather_pd(_mm256_loadu_si256((const __m256i *)(ks+i)),scale,8));
      l=_mm512_mul_pd(l,_mm512_i32gather_pd(_mm256_loadu_si256((const __m256i *)(parents+i)),scale,8));
    }
    const __m512d b=_mm512_i32gather_pd(_mm256_loadu_si256((const __m256i *)(qs+i)),rhs,8);
    d=_mm512_castsi512_pd(_mm512_and_si512(_mm512_castpd_si512(d),magnitude));
    l=_mm512_castsi512_pd(_mm512_and_si512(_mm512_castpd_si512(l),magnitude));
    const __m512d ab=_mm512_castsi512_pd(_mm512_and_si512(_mm512_castpd_si512(b),magnitude));
    if((_mm512_cmp_pd_mask(d,small,_CMP_GE_OQ)&_mm512_cmp_pd_mask(d,large,_CMP_LE_OQ)&
        _mm512_cmp_pd_mask(l,large,_CMP_LE_OQ)&_mm512_cmp_pd_mask(ab,large,_CMP_LE_OQ))!=0xff) return 0;
  }
  return kls_forest_leaf_validate_scalar(p,values,rhs,scale,i);
}
__attribute__((target("avx512f"),noinline))
static UF_long kls_forest_leaf_publish_avx512(kls_solver *s,
    const kls_forest_leaf_plan *p,const double *values,const double *rhs,
    const double *scale,double *forward,double *diag_recip) {
  const UF_long count=p->count;
  const int32_t *ks=p->indices,*parents=ks+count,*ds=parents+count;
  const int32_t *ls=ds+count,*qs=ls+count,*offset=qs+count;
  double *lu=(double *)s->numeric->LUbx[0],*diagonal=s->numeric->Udiag;
  UF_long i=0;
  for(;i+8<=count;i+=8) {
    const __m256i columns=_mm256_loadu_si256((const __m256i *)(ks+i));
    __m512d d=_mm512_i32gather_pd(_mm256_loadu_si256((const __m256i *)(ds+i)),values,8);
    __m512d l=_mm512_i32gather_pd(_mm256_loadu_si256((const __m256i *)(ls+i)),values,8);
    if(scale) {
      d=_mm512_mul_pd(d,_mm512_i32gather_pd(columns,scale,8));
      l=_mm512_mul_pd(l,_mm512_i32gather_pd(_mm256_loadu_si256((const __m256i *)(parents+i)),scale,8));
    }
    const __m512d inverse=_mm512_div_pd(_mm512_set1_pd(1.0),d);
    const __m512d b=_mm512_i32gather_pd(_mm256_loadu_si256((const __m256i *)(qs+i)),rhs,8);
    _mm512_i32scatter_pd(diagonal,columns,d,8);
    _mm512_i32scatter_pd(lu,_mm256_loadu_si256((const __m256i *)(offset+i)),_mm512_mul_pd(l,inverse),8);
    _mm512_i32scatter_pd(forward,columns,_mm512_mul_pd(b,inverse),8);
    if(diag_recip) _mm512_i32scatter_pd(diag_recip,columns,inverse,8);
  }
  return i;
}
#endif

/* Reordering independent divisions requires masked, nearest arithmetic. */
static int kls_forest_vector_environment(const kls_solver *s) {
#if defined(__x86_64__) && FLT_EVAL_METHOD==0 && KLS_HAVE_AVX512_KERNELS
  return s->n>=64 && s->options.backend==KLS_BACKEND_AUTO &&
    !s->stats.tuning_profile_active &&
    s->accuracy_policy==KLS_ACCURACY_COMPONENTWISE_BACKWARD_ERROR &&
    s->private_caller_depth==1 && !s->private_scaled_recovery &&
    !s->in_solve_refinement && s->solve_contract_probe==2 &&
    kls_forest_masked_environment() && fegetround()==FE_TONEAREST &&
    (_mm_getcsr()&0x6000u)==0 && kls_avx512_scatter_enabled();
#else
  (void)s;return 0;
#endif
}
#if KLS_HAVE_AVX512_KERNELS
__attribute__((target("avx512f"),noinline))
static void kls_forest_scale_reciprocals(UF_long n,const double *scale,double *inverse) {
  const __m512d one=_mm512_set1_pd(1.0);
  UF_long k=0;
  for(;k+8<=n;k+=8)
    _mm512_storeu_pd(inverse+k,_mm512_div_pd(one,_mm512_loadu_pd(scale+k)));
  for(;k<n;++k) inverse[k]=1.0/scale[k];
}
#endif

/* A declined prepass leaves native factors, forward RHS and FP status intact. */
static int kls_forest_leaf_begin(kls_solver *s,const double *values,
    const double *rhs,const double *scale,double *forward,double *diag_recip) {
#if defined(__x86_64__) && FLT_EVAL_METHOD==0 && KLS_HAVE_AVX512_KERNELS
  if(!rhs || !forward || !kls_forest_vector_environment(s) ||
     !kls_prepare_forest_leaves(s)) return 0;
  kls_forest_leaf_plan *p=s->forest_leaves;
  if(p->n!=s->n || p->numeric_owner!=s->numeric || p->lu_owner!=s->numeric->LUbx[0]) return 0;
  /* A read-only validator can be classified pure by IPA. An ordinary
     _mm_getcsr then sinks past that call; volatile assembly fixes the capture
     point and prevents input loads from moving before it. */
  unsigned saved;
  __asm__ volatile("stmxcsr %0":"=m"(saved)::"memory");
  const uintmax_t staged_entries=(uintmax_t)s->n+3*(uintmax_t)p->count;
  if(s->numeric->Work &&
     s->numeric->Xwork==s->numeric->Work && staged_entries<=SIZE_MAX/sizeof(double) &&
     s->numeric->worksize>=staged_entries*sizeof(double)) {
    double *stage=(double *)s->numeric->Xwork+s->n;
    if(!kls_forest_stage_leaves(p,values,rhs,scale,stage)) {
      __asm__ volatile("ldmxcsr %0"::"m"(saved):"memory");
      return 0;
    }
    kls_forest_publish_staged_leaves(s,p,stage,forward,diag_recip);
#ifdef KLS_FOREST_STAGE_TEST_EVENT
    KLS_FOREST_STAGE_TEST_EVENT();
#endif
#ifdef KLS_FOREST_STAGE_DIAGNOSTIC_EVENT
    KLS_FOREST_STAGE_DIAGNOSTIC_EVENT(s);
#endif
    return 1;
  }
  int vector=0,valid;
#if KLS_HAVE_AVX512_KERNELS
  vector=kls_avx512_scatter_enabled();
  valid=vector?kls_forest_leaf_validate_avx512(p,values,rhs,scale):
               kls_forest_leaf_validate_scalar(p,values,rhs,scale,0);
#else
  valid=kls_forest_leaf_validate_scalar(p,values,rhs,scale,0);
#endif
  if(!valid) {
    __asm__ volatile("ldmxcsr %0"::"m"(saved):"memory");
    return 0;
  }
  UF_long first=0;
#if KLS_HAVE_AVX512_KERNELS
  if(vector) first=kls_forest_leaf_publish_avx512(s,p,values,rhs,scale,forward,diag_recip);
#else
  (void)vector;
#endif
  kls_forest_leaf_publish_scalar(s,p,values,rhs,scale,forward,diag_recip,first);
  return 1;
#else
  (void)s;(void)values;(void)rhs;(void)scale;(void)forward;(void)diag_recip;return 0;
#endif
}

/* This stream is selected only after the complete exact forest proof and
 * independent-leaf publication. Every U producer has exactly one L entry,
 * targeting this column; no dense KLU workspace is read or written here. */
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline,aligned(64)))
#endif
static int kls_forest_leaf_remaining(kls_solver *s,const double *values,
    const double *rhs,const double *scale,double *y,double *diag_recip) {
  if(s->forest_leaves->remaining_plan)
    return kls_forest_remaining_planned(s,values,rhs,scale,y,diag_recip);
  const kls_forest_leaf_plan *p=s->forest_leaves;
  const int32_t *rest=p->indices+6*p->count;
  const trilinos_klu_l_numeric *numeric=s->numeric;
  double *lu=numeric->LUbx[0],*udiag=numeric->Udiag;
  const UF_long *lip=numeric->Lip,*uip=numeric->Uip,*ulen=numeric->Ulen;
  const UF_long *q=s->symbolic->Q;
  for(UF_long at=0;at<p->remaining;++at) {
    const UF_long k=(UF_long)rest[at],length=ulen[k];
    const UF_long *ui=length?(const UF_long *)(lu+uip[k]):NULL;
    double *ux=length?lu+uip[k]+kls_klu_units_for_indices(length):NULL;
    double diagonal=0,lower=0,yk=rhs[q[k]];
    for(UF_long j=s->refactor_col_ptr[k];j<s->refactor_col_ptr[k+1];++j) {
      double v=values[s->refactor_input_pos32[j]];
      if(scale) {
        const UF_long row=(UF_long)s->refactor_row_idx32[j];
        v=isfinite(scale[row])?v*scale[row]:v/((double *)numeric->Rs)[row];
      }
      const int32_t slot=s->forest_input_slot[j];
      if(slot>=0) ux[slot]=v;else if(slot==-1) diagonal=v;else lower=v;
    }
    for(UF_long up=0;up<length;++up) {
      const UF_long j=ui[up];const double u=ux[up];
#if defined(__FMA__)
      yk=fma(-u,y[j],yk);
#else
      yk-=u*y[j];
#endif
      if(u!=0) diagonal-=lu[lip[j]+kls_klu_units_for_indices(1)]*u;
    }
    if(diagonal==0) {
      s->common.status=TRILINOS_KLU_SINGULAR;
      if(s->common.numerical_rank==KLS_KLU_EMPTY) {
        s->common.numerical_rank=k;s->common.singular_col=q[k];
      }
      if(s->common.halt_if_singular) return 0;
    }
    udiag[k]=diagonal;const double inverse=1.0/diagonal;
    y[k]=isfinite(inverse)?yk*inverse:yk/diagonal;
    if(diag_recip) diag_recip[k]=inverse;
    if(numeric->Llen[k]) lu[lip[k]+kls_klu_units_for_indices(1)]=
      isfinite(inverse)?lower*inverse:lower/diagonal;
  }
  return 1;
}

/* Remaining columns cannot depend on a selected leaf in L^T: each selected
 * leaf's parent has a U child and is therefore in the remaining stream. */
static UF_long kls_forest_leaf_finish(kls_solver *s,double *work,double *out) {
  kls_forest_leaf_plan *p=s->forest_leaves;
  const UF_long count=p->count,n=p->n;
  const int32_t *ks=p->indices,*parents=ks+count,*rest=ks+6*count;
  const double *rs=s->numeric->Rs;
  const double *scale=rs?work+n:NULL;
  for (UF_long at=p->remaining;at-- >0;) {
    const UF_long k=(UF_long)rest[at];
    double value=work[k];
    if (s->numeric->Llen[k]) {
      UF_long *li,length;double *lx;
      kls_klu_get_pointer(s->numeric->LUbx[0],s->numeric->Lip,
        s->numeric->Llen,k,&li,&lx,&length);
      value-=lx[0]*work[li[0]];
    }
    work[k]=value;
    out[s->numeric->Pnum[k]]=rs==NULL?value:
      (isfinite(scale[k])?value*scale[k]:value/rs[k]);
  }
  const int32_t *offset=ks+5*count;
  const double *lu=(const double *)s->numeric->LUbx[0];
  for(UF_long i=0;i<count;++i) {
    const UF_long k=(UF_long)ks[i];
    const double value=work[k]-lu[offset[i]]*work[parents[i]];
    work[k]=value;
    out[s->numeric->Pnum[k]]=rs==NULL?value:
      (isfinite(scale[k])?value*scale[k]:value/rs[k]);
  }
  if (s->forest_leaf_solves<UINT64_MAX) ++s->forest_leaf_solves;
  return 1;
}
