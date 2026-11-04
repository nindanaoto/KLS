#ifndef KLS_SNB_ROW_PANEL_INC
#define KLS_SNB_ROW_PANEL_INC
/* Owned transactional row-panel arithmetic; exceptional jobs are declined.
 * Width-16 targets use two eight-column row streams. Producer panels may use
 * this layout or the incumbent column layout; native LU remains unchanged.
 * A production caller must publish a row panel only after all getrf stores:
 * full-vector stores also rewrite earlier, already-finalized columns. */
#if defined(__AVX512F__) && defined(__AVX512DQ__) && FLT_EVAL_METHOD==0 && defined(__GNUC__) && !defined(__clang__) && defined(__OPTIMIZE__)
#define KLS_SNB_HAVE_ROW_PANEL 1
/* The incumbent masked row tail rounds multiplication before addition;
 * its complete eight-row vectors contract to FMA. Preserve that distinction
 * when vectorizing over columns instead of producer rows. */
static KLS_ALWAYS_INLINE __m512d kls_rpanel_mul_vec(__m512d a,__m512d b) {
  __asm__ volatile("vmulpd %1, %0, %0":"+v"(a):"v"(b));return a;
}
static KLS_ALWAYS_INLINE __m512d kls_rpanel_product_add(__m512d t,__m512d l,__m512d b,int separate) {
  if(!separate)return _mm512_fmadd_pd(l,b,t);
  l=kls_rpanel_mul_vec(l,b);
  __asm__ volatile("vaddpd %1, %0, %0":"+v"(t):"v"(l));return t;
}
static KLS_ALWAYS_INLINE size_t kls_rpanel_at(int H,int c,int r) {
  return (size_t)(c/8)*8*(size_t)H+(size_t)r*8+(unsigned)c%8;
}
static KLS_ALWAYS_INLINE const double *kls_rpanel_column(
    const double *panel,int H,int c,int row_layout) {
  return row_layout?panel+(size_t)(c/8)*8*H+(unsigned)c%8:panel+(size_t)c*H;
}
static KLS_ALWAYS_INLINE __mmask8 kls_rpanel_nonzero(__m512d value,__mmask8 valid) {
  if(!valid)return 0;
  __mmask8 mask;const __m512d zero=_mm512_setzero_pd();
  if(valid==255u) {
    __asm__ volatile("vcmppd $4, %2, %1, %0":"=k"(mask):"v"(value),"v"(zero):"cc");
  } else {
    __asm__ volatile("vcmppd $4, %2, %1, %0%{%3%}":"=k"(mask):"v"(value),"v"(zero),"Yk"(valid):"cc");
  }
  /* Preserve the scalar comparison's denormal flag even on compare forms
     whose packed implementation does not report it. One scalar comparison
     suffices for the sticky flag; its DAZ behavior is the caller's behavior. */
  const __m512i bits=_mm512_and_si512(_mm512_castpd_si512(value),_mm512_set1_epi64(INT64_MAX));
  __mmask8 denormal=valid & _mm512_cmpneq_epi64_mask(bits,_mm512_setzero_si512()) &
    _mm512_cmplt_epi64_mask(bits,_mm512_set1_epi64(INT64_C(0x0010000000000000)));
  if(denormal) {
    double lanes[8];_mm512_storeu_pd(lanes,value);
    const double one=lanes[__builtin_ctz((unsigned)denormal)],z=0.;
    __asm__ volatile("ucomisd %1, %0"::"x"(one),"x"(z):"cc");
  }
  return mask;
}

/* Each result retains the original producer order. A producer's below rows
 * cannot overlap the upper rows solved by this edge, so dense row products
 * may be subtracted immediately, without the old full-edge T buffer. */
static __attribute__((noinline)) int kls_rpanel_edge(
    const kls_snb_block *blk,const kls_snb_snode *sn,const kls_snb_edge *edge,
    int narrow,double *R,int producer_row_layout,int product_math) {
  const kls_snb_snode *ps=&blk->snodes[edge->producer];
  if(sn->width!=16 || edge->sub_count<1 || edge->sub_count>32 ||
     (producer_row_layout && ps->width!=16))return 0;
#ifdef KLS_RPANEL_TEST_EDGE
  KLS_RPANEL_TEST_EDGE();
#endif
  const int H=sn->height,sc=edge->sub_count,ph=ps->below,ld=ps->height;
  const int stride=producer_row_layout?8:1;
  const int32_t *sub=blk->subcols+edge->sub_off,*dst=blk->edge_dsts+edge->dst_off;
  const double *panel=blk->panels+ps->panel;
  double *R0=R,*R1=R+(size_t)8*H;
  for(int i=0;i<sc;++i) {
    const __m512d u0=_mm512_loadu_pd(R0+(size_t)(edge->u_start+i)*8);
    const __m512d u1=_mm512_loadu_pd(R1+(size_t)(edge->u_start+i)*8);
    const __mmask8 m0=kls_rpanel_nonzero(u0,255),m1=kls_rpanel_nonzero(u1,255);
    const double *col=kls_rpanel_column(panel,ld,sub[i],producer_row_layout);
    if(m0||m1) {
      for(int j=i+1;j<sc;++j) {
        const __m512d l=_mm512_set1_pd(col[(size_t)(ps->uw+sub[j])*stride]);
        double *out0=R0+(size_t)(edge->u_start+j)*8,*out1=R1+(size_t)(edge->u_start+j)*8;
        if(m0)_mm512_storeu_pd(out0,_mm512_mask3_fnmadd_pd(l,u0,_mm512_loadu_pd(out0),m0));
        if(m1)_mm512_storeu_pd(out1,_mm512_mask3_fnmadd_pd(l,u1,_mm512_loadu_pd(out1),m1));
      }
      if(sc<narrow)for(int r=0;r<ph;++r)if(dst[r]>=0) {
        const __m512d l=_mm512_set1_pd(col[(size_t)(ps->uw+ps->width+r)*stride]);
        double *out0=R0+(size_t)dst[r]*8,*out1=R1+(size_t)dst[r]*8;
        if(m0)_mm512_storeu_pd(out0,_mm512_mask3_fnmadd_pd(l,u0,_mm512_loadu_pd(out0),m0));
        if(m1)_mm512_storeu_pd(out1,_mm512_mask3_fnmadd_pd(l,u1,_mm512_loadu_pd(out1),m1));
      }
    }
  }
  if(sc<narrow || !ph)return 1;
  int r=0;
  for(;r+3<ph;r+=4) {
    const int category=r<(ph&~15)?0:r<(ph&~7)?1:2;
    const int separate=!(product_math&(1<<category));
    const double *lc=kls_rpanel_column(panel,ld,sub[0],producer_row_layout)+(size_t)(ps->uw+ps->width+r)*stride;
    __m512d b0=_mm512_loadu_pd(R0+(size_t)edge->u_start*8),b1=_mm512_loadu_pd(R1+(size_t)edge->u_start*8);
    __m512d l0=_mm512_set1_pd(lc[0]),l1=_mm512_set1_pd(lc[stride]);
    __m512d l2=_mm512_set1_pd(lc[2*stride]),l3=_mm512_set1_pd(lc[3*stride]);
    __m512d t00=kls_rpanel_mul_vec(l0,b0),t01=kls_rpanel_mul_vec(l0,b1);
    __m512d t10=kls_rpanel_mul_vec(l1,b0),t11=kls_rpanel_mul_vec(l1,b1);
    __m512d t20=kls_rpanel_mul_vec(l2,b0),t21=kls_rpanel_mul_vec(l2,b1);
    __m512d t30=kls_rpanel_mul_vec(l3,b0),t31=kls_rpanel_mul_vec(l3,b1);
    for(int i=1;i<sc;++i) {
      lc=kls_rpanel_column(panel,ld,sub[i],producer_row_layout)+(size_t)(ps->uw+ps->width+r)*stride;
      b0=_mm512_loadu_pd(R0+(size_t)(edge->u_start+i)*8);b1=_mm512_loadu_pd(R1+(size_t)(edge->u_start+i)*8);
      l0=_mm512_set1_pd(lc[0]);l1=_mm512_set1_pd(lc[stride]);l2=_mm512_set1_pd(lc[2*stride]);l3=_mm512_set1_pd(lc[3*stride]);
      t00=kls_rpanel_product_add(t00,l0,b0,separate);t01=kls_rpanel_product_add(t01,l0,b1,separate);
      t10=kls_rpanel_product_add(t10,l1,b0,separate);t11=kls_rpanel_product_add(t11,l1,b1,separate);
      t20=kls_rpanel_product_add(t20,l2,b0,separate);t21=kls_rpanel_product_add(t21,l2,b1,separate);
      t30=kls_rpanel_product_add(t30,l3,b0,separate);t31=kls_rpanel_product_add(t31,l3,b1,separate);
    }
#define KLS_RP_SUB(J,A,B) do {if(dst[r+(J)]>=0) {double *o0=R0+(size_t)dst[r+(J)]*8,*o1=R1+(size_t)dst[r+(J)]*8; _mm512_storeu_pd(o0,_mm512_sub_pd(_mm512_loadu_pd(o0),(A)));_mm512_storeu_pd(o1,_mm512_sub_pd(_mm512_loadu_pd(o1),(B)));}}while(0)
    KLS_RP_SUB(0,t00,t01);KLS_RP_SUB(1,t10,t11);KLS_RP_SUB(2,t20,t21);KLS_RP_SUB(3,t30,t31);
#undef KLS_RP_SUB
  }
  for(;r<ph;++r) {
    const int category=r<(ph&~15)?0:r<(ph&~7)?1:2;
    const int separate=!(product_math&(1<<category));
    const double *lc=kls_rpanel_column(panel,ld,sub[0],producer_row_layout);
    __m512d l=_mm512_set1_pd(lc[(size_t)(ps->uw+ps->width+r)*stride]);
    __m512d t0=kls_rpanel_mul_vec(l,_mm512_loadu_pd(R0+(size_t)edge->u_start*8));
    __m512d t1=kls_rpanel_mul_vec(l,_mm512_loadu_pd(R1+(size_t)edge->u_start*8));
    for(int i=1;i<sc;++i) {
      lc=kls_rpanel_column(panel,ld,sub[i],producer_row_layout);l=_mm512_set1_pd(lc[(size_t)(ps->uw+ps->width+r)*stride]);
      t0=kls_rpanel_product_add(t0,l,_mm512_loadu_pd(R0+(size_t)(edge->u_start+i)*8),separate);
      t1=kls_rpanel_product_add(t1,l,_mm512_loadu_pd(R1+(size_t)(edge->u_start+i)*8),separate);
    }
    if(dst[r]>=0) {
      double *o0=R0+(size_t)dst[r]*8,*o1=R1+(size_t)dst[r]*8;
      _mm512_storeu_pd(o0,_mm512_sub_pd(_mm512_loadu_pd(o0),t0));
      _mm512_storeu_pd(o1,_mm512_sub_pd(_mm512_loadu_pd(o1),t1));
    }
  }
  return 1;
}
static KLS_ALWAYS_INLINE void kls_rpanel_getrf_row(double *R,int H,int row,
    int pivot_col,double value,__m512d u0,__m512d u1,__mmask8 m0,__mmask8 m1) {
  double *p0=R+(size_t)row*8,*p1=R+(size_t)8*H+(size_t)row*8;
  __m512d v0=_mm512_loadu_pd(p0),v1=_mm512_loadu_pd(p1),l=_mm512_set1_pd(value);
  if(pivot_col<8)v0=_mm512_mask_mov_pd(v0,(__mmask8)(1u<<pivot_col),l);
  else v1=_mm512_mask_mov_pd(v1,(__mmask8)(1u<<(pivot_col-8)),l);
  v0=_mm512_mask3_fnmadd_pd(l,u0,v0,m0);v1=_mm512_mask3_fnmadd_pd(l,u1,v1,m1);
  _mm512_storeu_pd(p0,v0);_mm512_storeu_pd(p1,v1);
}
static KLS_ALWAYS_INLINE void kls_rpanel_transpose8(__m512d a[8]) {
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
/* Transpose each completed producer once, before publication. Bitwise moves
 * preserve every value and flag; only producers with mixed consumers pay. */
static void kls_rpanel_to_columns(const double *R,int H,double *C) {
  for(int group=0;group<2;++group) {
    int r=0;
    for(;r+7<H;r+=8) {
      __m512d tile[8];
      for(int i=0;i<8;++i)tile[i]=_mm512_loadu_pd(R+(size_t)group*8*H+(size_t)(r+i)*8);
      kls_rpanel_transpose8(tile);
      for(int c=0;c<8;++c)_mm512_storeu_pd(C+(size_t)(group*8+c)*H+r,tile[c]);
    }
    for(;r<H;++r)for(int c=0;c<8;++c)C[(size_t)(group*8+c)*H+r]=R[(size_t)group*8*H+(size_t)r*8+c];
  }
}
static void kls_rpanel_from_columns(const double *C,int H,double *R) {
  for(int group=0;group<2;++group) {
    int r=0;
    for(;r+7<H;r+=8) {
      __m512d tile[8];
      for(int c=0;c<8;++c)tile[c]=_mm512_loadu_pd(C+(size_t)(group*8+c)*H+r);
      kls_rpanel_transpose8(tile);
      for(int i=0;i<8;++i)_mm512_storeu_pd(R+(size_t)group*8*H+(size_t)(r+i)*8,tile[i]);
    }
    for(;r<H;++r)for(int c=0;c<8;++c)R[(size_t)group*8*H+(size_t)r*8+c]=C[(size_t)(group*8+c)*H+r];
  }
}
/* Returns the first zero pivot or -1. No publication occurs in this kernel. */
static __attribute__((noinline)) int kls_rpanel_getrf(double *R,int H,int uw,double *D) {
#ifdef KLS_RPANEL_TEST_GETRF
  KLS_RPANEL_TEST_GETRF();
#endif
  int singular=-1;
  for(int i=0;i<16;++i) {
    const int pivot_row=uw+i;const double pivot=R[kls_rpanel_at(H,i,pivot_row)];D[i]=pivot;
    if(pivot==0. && singular<0)singular=i;
    const __m512d u0=_mm512_loadu_pd(R+(size_t)pivot_row*8);
    const __m512d u1=_mm512_loadu_pd(R+(size_t)8*H+(size_t)pivot_row*8);
    const __mmask8 valid0=i<7?(__mmask8)(255u<<(i+1)):0;
    const __mmask8 valid1=i<8?255:(__mmask8)(255u<<(i-7));
    const __mmask8 m0=kls_rpanel_nonzero(u0,valid0),m1=kls_rpanel_nonzero(u1,valid1);
    const __m256d divisor=_mm256_set1_pd(pivot);
    int r=pivot_row+1;
    for(;r+3<H;r+=4) {
      const __m256d numerators=_mm256_setr_pd(R[kls_rpanel_at(H,i,r)],R[kls_rpanel_at(H,i,r+1)],R[kls_rpanel_at(H,i,r+2)],R[kls_rpanel_at(H,i,r+3)]);
      double q[4];_mm256_storeu_pd(q,_mm256_div_pd(numerators,divisor));
      kls_rpanel_getrf_row(R,H,r,i,q[0],u0,u1,m0,m1);kls_rpanel_getrf_row(R,H,r+1,i,q[1],u0,u1,m0,m1);
      kls_rpanel_getrf_row(R,H,r+2,i,q[2],u0,u1,m0,m1);kls_rpanel_getrf_row(R,H,r+3,i,q[3],u0,u1,m0,m1);
    }
    for(;r<H;++r) {
      const double q=R[kls_rpanel_at(H,i,r)]/pivot;
      kls_rpanel_getrf_row(R,H,r,i,q,u0,u1,m0,m1);
    }
  }
  return singular;
}
#endif

#endif
