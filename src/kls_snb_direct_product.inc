#if defined(__AVX512F__)
static KLS_ALWAYS_INLINE void kls_snb_direct_store8(double *W,int height,
    int c,int r,const int32_t *dst,__m512d value) {
  int contiguous=dst[r]>=0;
  for(int k=1;k<8;++k)contiguous&=(int64_t)dst[r+k]==(int64_t)dst[r]+k;
  double *target=W+(size_t)c*height;
  if(contiguous) {
    double *out=target+dst[r];
    _mm512_storeu_pd(out,_mm512_sub_pd(_mm512_loadu_pd(out),value));
  } else {
    _Alignas(64) double tmp[8];_mm512_storeu_pd(tmp,value);
    for(int k=0;k<8;++k)if(dst[r+k]>=0)target[dst[r+k]]-=tmp[k];
  }
}
#endif

/* Development direct-store product: the builder proves destinations are
   disjoint from the upper B block. Retain producer reduction order. */
static int kls_snb_direct_product(const double *L,int ld,const int32_t *sub,
    int sc,const double *B,int bld,int w,int rows,const int32_t *dst,
    double *W,int height) {
#if defined(__AVX512F__)
  for(int c=0;c<w;c+=4) {
    const int nc=w-c<4?w-c:4;
    int r=0;
    if(nc==4)for(;r+15<rows;r+=16) {
      const double *b0=B+(size_t)c*bld,*b1=b0+bld,*b2=b1+bld,*b3=b2+bld;
      __m512d l0=_mm512_loadu_pd(L+(size_t)sub[0]*ld+r);
      __m512d l1=_mm512_loadu_pd(L+(size_t)sub[0]*ld+r+8);
      __m512d t00=_mm512_mul_pd(l0,_mm512_set1_pd(b0[0]));
      __m512d t01=_mm512_mul_pd(l1,_mm512_set1_pd(b0[0]));
      __m512d t10=_mm512_mul_pd(l0,_mm512_set1_pd(b1[0]));
      __m512d t11=_mm512_mul_pd(l1,_mm512_set1_pd(b1[0]));
      __m512d t20=_mm512_mul_pd(l0,_mm512_set1_pd(b2[0]));
      __m512d t21=_mm512_mul_pd(l1,_mm512_set1_pd(b2[0]));
      __m512d t30=_mm512_mul_pd(l0,_mm512_set1_pd(b3[0]));
      __m512d t31=_mm512_mul_pd(l1,_mm512_set1_pd(b3[0]));
      for(int i=1;i<sc;++i) {
        l0=_mm512_loadu_pd(L+(size_t)sub[i]*ld+r);
        l1=_mm512_loadu_pd(L+(size_t)sub[i]*ld+r+8);
        t00=_mm512_add_pd(t00,_mm512_mul_pd(l0,_mm512_set1_pd(b0[i])));
        t01=_mm512_add_pd(t01,_mm512_mul_pd(l1,_mm512_set1_pd(b0[i])));
        t10=_mm512_add_pd(t10,_mm512_mul_pd(l0,_mm512_set1_pd(b1[i])));
        t11=_mm512_add_pd(t11,_mm512_mul_pd(l1,_mm512_set1_pd(b1[i])));
        t20=_mm512_add_pd(t20,_mm512_mul_pd(l0,_mm512_set1_pd(b2[i])));
        t21=_mm512_add_pd(t21,_mm512_mul_pd(l1,_mm512_set1_pd(b2[i])));
        t30=_mm512_add_pd(t30,_mm512_mul_pd(l0,_mm512_set1_pd(b3[i])));
        t31=_mm512_add_pd(t31,_mm512_mul_pd(l1,_mm512_set1_pd(b3[i])));
      }
      kls_snb_direct_store8(W,height,c+0,r+0,dst,t00);
      kls_snb_direct_store8(W,height,c+0,r+8,dst,t01);
      kls_snb_direct_store8(W,height,c+1,r+0,dst,t10);
      kls_snb_direct_store8(W,height,c+1,r+8,dst,t11);
      kls_snb_direct_store8(W,height,c+2,r+0,dst,t20);
      kls_snb_direct_store8(W,height,c+2,r+8,dst,t21);
      kls_snb_direct_store8(W,height,c+3,r+0,dst,t30);
      kls_snb_direct_store8(W,height,c+3,r+8,dst,t31);
    }
    for(;r<rows;r+=8) {
      const int nr=rows-r<8?rows-r:8;
      const __mmask8 mask=(__mmask8)((1u<<nr)-1u);
      __m512d sums[4];
      __m512d l=_mm512_maskz_loadu_pd(mask,L+(size_t)sub[0]*ld+r);
      if(nc==4) {
        const double *b0=B+(size_t)c*bld,*b1=b0+bld,*b2=b1+bld,*b3=b2+bld;
        __m512d t0=_mm512_maskz_mul_pd(mask,l,_mm512_set1_pd(b0[0]));
        __m512d t1=_mm512_maskz_mul_pd(mask,l,_mm512_set1_pd(b1[0]));
        __m512d t2=_mm512_maskz_mul_pd(mask,l,_mm512_set1_pd(b2[0]));
        __m512d t3=_mm512_maskz_mul_pd(mask,l,_mm512_set1_pd(b3[0]));
        for(int i=1;i<sc;++i) {
          l=_mm512_maskz_loadu_pd(mask,L+(size_t)sub[i]*ld+r);
          t0=_mm512_add_pd(t0,_mm512_maskz_mul_pd(mask,l,_mm512_set1_pd(b0[i])));
          t1=_mm512_add_pd(t1,_mm512_maskz_mul_pd(mask,l,_mm512_set1_pd(b1[i])));
          t2=_mm512_add_pd(t2,_mm512_maskz_mul_pd(mask,l,_mm512_set1_pd(b2[i])));
          t3=_mm512_add_pd(t3,_mm512_maskz_mul_pd(mask,l,_mm512_set1_pd(b3[i])));
        }
        sums[0]=t0;sums[1]=t1;sums[2]=t2;sums[3]=t3;
      } else {
        for(int j=0;j<nc;++j)
          sums[j]=_mm512_maskz_mul_pd(mask,l,_mm512_set1_pd(B[(size_t)(c+j)*bld]));
        for(int i=1;i<sc;++i) {
          l=_mm512_maskz_loadu_pd(mask,L+(size_t)sub[i]*ld+r);
          for(int j=0;j<nc;++j)
            sums[j]=_mm512_add_pd(sums[j],_mm512_maskz_mul_pd(mask,l,
              _mm512_set1_pd(B[(size_t)(c+j)*bld+i])));
        }
      }
      int contiguous=dst[r]>=0;
      for(int k=1;k<nr;++k)contiguous&=(int64_t)dst[r+k]==(int64_t)dst[r]+k;
      for(int j=0;j<nc;++j) {
        double *target=W+(size_t)(c+j)*height;
        if(contiguous) {
          double *out=target+dst[r];
          _mm512_mask_storeu_pd(out,mask,_mm512_mask_sub_pd(
            _mm512_setzero_pd(),mask,_mm512_maskz_loadu_pd(mask,out),sums[j]));
        } else {
          _Alignas(64) double tmp[8];
          _mm512_mask_storeu_pd(tmp,mask,sums[j]);
          for(int k=0;k<nr;++k)if(dst[r+k]>=0)target[dst[r+k]]-=tmp[k];
        }
      }
    }
  }
  return 1;
#else
  (void)L;(void)ld;(void)sub;(void)sc;(void)B;(void)bld;(void)w;
  (void)rows;(void)dst;(void)W;(void)height;
  return 0;
#endif
}

