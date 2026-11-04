/* Width-16 target updates in a temporary pair of row-major eight-column
 * tiles. Producer panels remain column-major. No factor is published here. */
#if defined(__AVX512F__) && defined(__AVX512DQ__) && FLT_EVAL_METHOD==0
static KLS_ALWAYS_INLINE void kls_snb_transpose8(__m512d a[8]) {
  __m512d t[8],s[8];
  for(int i=0;i<8;i+=2){t[i]=_mm512_unpacklo_pd(a[i],a[i+1]);t[i+1]=_mm512_unpackhi_pd(a[i],a[i+1]);}
  for(int g=0;g<8;g+=4) {
    s[g]=_mm512_shuffle_f64x2(t[g],t[g+2],0x44);
    s[g+1]=_mm512_shuffle_f64x2(t[g+1],t[g+3],0x44);
    s[g+2]=_mm512_shuffle_f64x2(t[g],t[g+2],0xee);
    s[g+3]=_mm512_shuffle_f64x2(t[g+1],t[g+3],0xee);
  }
  a[0]=_mm512_shuffle_f64x2(s[0],s[4],0x88);a[1]=_mm512_shuffle_f64x2(s[1],s[5],0x88);
  a[2]=_mm512_shuffle_f64x2(s[0],s[4],0xdd);a[3]=_mm512_shuffle_f64x2(s[1],s[5],0xdd);
  a[4]=_mm512_shuffle_f64x2(s[2],s[6],0x88);a[5]=_mm512_shuffle_f64x2(s[3],s[7],0x88);
  a[6]=_mm512_shuffle_f64x2(s[2],s[6],0xdd);a[7]=_mm512_shuffle_f64x2(s[3],s[7],0xdd);
}
static KLS_ALWAYS_INLINE int kls_snb_row8_scalar(double value) {
  uint64_t bits;memcpy(&bits,&value,sizeof bits);bits&=UINT64_C(0x7fffffffffffffff);
  return !bits || (bits>=UINT64_C(0x0010000000000000) && bits<UINT64_C(0x7ff0000000000000));
}
static KLS_ALWAYS_INLINE int kls_snb_row8_mask(__m512d value,__mmask8 *active) {
  __m512i bits=_mm512_and_si512(_mm512_castpd_si512(value),_mm512_set1_epi64(INT64_MAX));
  __m512i exp=_mm512_and_si512(bits,_mm512_set1_epi64(INT64_C(0x7ff0000000000000)));
  *active=_mm512_cmpneq_epi64_mask(bits,_mm512_setzero_si512());
  __mmask8 bad=_mm512_cmpeq_epi64_mask(exp,_mm512_set1_epi64(INT64_C(0x7ff0000000000000))) |
      (*active&_mm512_cmpeq_epi64_mask(exp,_mm512_setzero_si512()));
  return !bad;
}
static __attribute__((noinline)) int kls_snb_row8_edge(const kls_snb_block *blk,
    const kls_snb_snode *sn,const kls_snb_edge *edge,int narrow,double *W,double *T) {
  if(sn->width!=16)return 0;
  const kls_snb_snode *ps=&blk->snodes[edge->producer];
  const int H=sn->height,sc=edge->sub_count,ph=ps->below,ld=ps->height;
  if(sc<1 || sc>32)return 0;
  const int32_t *sub=blk->subcols+edge->sub_off,*dst=blk->edge_dsts+edge->dst_off;
  const double *L=blk->panels+ps->panel+ps->uw,*Lb=L+ps->width;
  /* Below destinations cannot overlap this edge's upper block. Finish its
     unit-lower solve first, keeping producer order within every result. */
  for(int group=0;group<2;++group)for(int i=0;i<sc;++i) {
    double *base=W+(size_t)group*8*H;
    __m512d u=_mm512_loadu_pd(base+(size_t)(edge->u_start+i)*8);__mmask8 active;
    if(!kls_snb_row8_mask(u,&active))return 0;
    if(!active)continue;
    for(int j=i+1;j<sc;++j) {
      double l=L[(size_t)sub[i]*ld+sub[j]];if(!kls_snb_row8_scalar(l))return 0;
      double *out=base+(size_t)(edge->u_start+j)*8;
      _mm512_storeu_pd(out,_mm512_mask3_fnmadd_pd(_mm512_set1_pd(l),u,_mm512_loadu_pd(out),active));
    }
  }
  if(!ph)return 1;
  if(sc<narrow) {
    for(int group=0;group<2;++group)for(int i=0;i<sc;++i) {
      double *base=W+(size_t)group*8*H;
      __m512d u=_mm512_loadu_pd(base+(size_t)(edge->u_start+i)*8);__mmask8 active;
      if(!kls_snb_row8_mask(u,&active))return 0;
      if(!active)continue;
      for(int r=0;r<ph;++r)if(dst[r]>=0) {
        double l=Lb[(size_t)sub[i]*ld+r];if(!kls_snb_row8_scalar(l))return 0;
        double *out=base+(size_t)dst[r]*8;
        _mm512_storeu_pd(out,_mm512_mask3_fnmadd_pd(_mm512_set1_pd(l),u,_mm512_loadu_pd(out),active));
      }
    }
  } else {
    double B[32*16];
    for(int c=0;c<16;++c)for(int i=0;i<sc;++i)
      B[(size_t)c*sc+i]=W[(size_t)(c/8)*8*H+(size_t)(edge->u_start+i)*8+c%8];
    kls_snb_dense_product(Lb,ld,sub,sc,B,sc,16,0,ph,T);
    for(int group=0;group<2;++group) {
      int r=0;double *base=W+(size_t)group*8*H;
      for(;r+7<ph;r+=8) {
        __m512d tile[8];
#pragma GCC unroll 8
        for(int c=0;c<8;++c)tile[c]=_mm512_loadu_pd(T+(size_t)(group*8+c)*ph+r);
        kls_snb_transpose8(tile);
#pragma GCC unroll 8
        for(int j=0;j<8;++j)if(dst[r+j]>=0) {
          double *out=base+(size_t)dst[r+j]*8;
          _mm512_storeu_pd(out,_mm512_sub_pd(_mm512_loadu_pd(out),tile[j]));
        }
      }
      for(;r<ph;++r)if(dst[r]>=0) {
        double values[8];for(int c=0;c<8;++c)values[c]=T[(size_t)(group*8+c)*ph+r];
        double *out=base+(size_t)dst[r]*8;
        _mm512_storeu_pd(out,_mm512_sub_pd(_mm512_loadu_pd(out),_mm512_loadu_pd(values)));
      }
    }
  }
  return 1;
}
#endif
