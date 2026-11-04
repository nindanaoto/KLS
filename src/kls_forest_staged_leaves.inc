/* Current-call scratch only. Validate all leaves before publishing any native
   factor, then reuse the rounded values instead of gathering them twice. */
#if KLS_HAVE_AVX512_KERNELS
__attribute__((target("avx512f"),noinline))
static int kls_forest_stage_leaves(const kls_forest_leaf_plan *p,
    const double *values,const double *rhs,const double *scale,double *stage) {
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
    const __m512d ad=_mm512_castsi512_pd(_mm512_and_si512(_mm512_castpd_si512(d),magnitude));
    const __m512d al=_mm512_castsi512_pd(_mm512_and_si512(_mm512_castpd_si512(l),magnitude));
    const __m512d ab=_mm512_castsi512_pd(_mm512_and_si512(_mm512_castpd_si512(b),magnitude));
    if((_mm512_cmp_pd_mask(ad,small,_CMP_GE_OQ)&_mm512_cmp_pd_mask(ad,large,_CMP_LE_OQ)&
        _mm512_cmp_pd_mask(al,large,_CMP_LE_OQ)&_mm512_cmp_pd_mask(ab,large,_CMP_LE_OQ))!=0xff)return 0;
    _mm512_storeu_pd(stage+i,d);
    _mm512_storeu_pd(stage+count+i,l);
    _mm512_storeu_pd(stage+2*count+i,b);
  }
  for(;i<count;++i) {
    const double d=scale?values[ds[i]]*scale[ks[i]]:values[ds[i]];
    const double l=scale?values[ls[i]]*scale[parents[i]]:values[ls[i]];
    const double b=rhs[qs[i]];
    if(!(fabs(d)>=0x1p-256 && fabs(d)<=0x1p256 &&
         fabs(l)<=0x1p256 && fabs(b)<=0x1p256))return 0;
    stage[i]=d;stage[count+i]=l;stage[2*count+i]=b;
  }
  return 1;
}

__attribute__((target("avx512f"),noinline))
static void kls_forest_publish_staged_leaves(kls_solver *s,
    const kls_forest_leaf_plan *p,const double *stage,double *forward,double *diag_recip) {
  const UF_long count=p->count;
  const int32_t *ks=p->indices,*offset=ks+5*count;
  double *lu=s->numeric->LUbx[0],*diagonal=s->numeric->Udiag;
  UF_long i=0;
  for(;i+8<=count;i+=8) {
    const __m256i columns=_mm256_loadu_si256((const __m256i *)(ks+i));
    const __m512d d=_mm512_loadu_pd(stage+i),l=_mm512_loadu_pd(stage+count+i);
    const __m512d b=_mm512_loadu_pd(stage+2*count+i);
    const __m512d inverse=_mm512_div_pd(_mm512_set1_pd(1.0),d);
    _mm512_i32scatter_pd(diagonal,columns,d,8);
    _mm512_i32scatter_pd(lu,_mm256_loadu_si256((const __m256i *)(offset+i)),_mm512_mul_pd(l,inverse),8);
    _mm512_i32scatter_pd(forward,columns,_mm512_mul_pd(b,inverse),8);
    if(diag_recip)_mm512_i32scatter_pd(diag_recip,columns,inverse,8);
  }
  for(;i<count;++i) {
    const UF_long k=(UF_long)ks[i];
    const double d=stage[i],l=stage[count+i],inverse=1.0/d;
    diagonal[k]=d;lu[offset[i]]=l*inverse;forward[k]=stage[2*count+i]*inverse;
    if(diag_recip)diag_recip[k]=inverse;
  }
}
#endif
