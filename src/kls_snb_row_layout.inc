/* Change only private panel addressing. Native factors and all dependency
 * geometry remain unchanged; coefficients are freshly scattered each call. */
#ifdef KLS_SNB_HAVE_ROW_PANEL
/* Preserve the incumbent's actual contraction choices in this binary. GCC
 * instrumentation can change contraction around intrinsic loads. This is a
 * coefficient-independent capability probe, not a saved numerical answer.
 * The volatile function pointer prevents a constant-argument clone. */
static int kls_snb_row_product_math(void) {
  const unsigned csr=_mm_getcsr();
  if((csr&0x1f80u)!=0x1f80u||(csr&0x6000u))return -1;
  double L[32],B[8],T[64];const int32_t sub[2]={0,1};
  for(int r=0;r<16;++r){L[r]=-1.;L[16+r]=1.+0x1p-27;}
  for(int c=0;c<4;++c){B[2*c]=1.;B[2*c+1]=1.-0x1p-27;}
  void (*volatile product)(const double*,int32_t,const int32_t*,int32_t,
      const double*,int32_t,int32_t,int32_t,int32_t,double*)=kls_snb_dense_product;
  const int widths[5]={4,4,4,1,1},rows[5]={16,8,1,8,1};
  const double fused=-0x1p-54;uint64_t fused_bits;memcpy(&fused_bits,&fused,8);
  int result=0;
  for(int mode=0;mode<5;++mode) {
    product(L,16,sub,2,B,2,widths[mode],0,rows[mode],T);
    uint64_t bits;memcpy(&bits,T,8);
    if(bits==fused_bits)result|=1<<mode;
    else if(bits!=0){result=-1;break;}
  }
  _mm_setcsr(csr);return result;
}
#endif
/* A repaired row node uses the incumbent column kernel on its next visit.
 * A clean column visit permits another row trial; persistent underflow keeps
 * the column choice. Current coefficients are always factored afresh.
 * Atomic bytes let independent workers mark nodes sharing the same byte. */
static size_t kls_snb_row_preference_bytes(const kls_snb_state *st) {
  if(!st||st->nblocks!=1||!st->blocks||st->blocks[0].snode_count<1)return SIZE_MAX;
  if(st->row_column_preferred)return 0;
  const uint64_t slots=((uint64_t)st->blocks[0].snode_count+7)/8;
  if(slots>SIZE_MAX/sizeof(*st->row_column_preferred))return SIZE_MAX;
  return (size_t)slots*sizeof(*st->row_column_preferred);
}
static int kls_snb_prefers_column(const kls_snb_state *st,int64_t node) {
  return st->row_column_preferred&&
    (atomic_load_explicit(st->row_column_preferred+(size_t)node/8,memory_order_relaxed)&
      (unsigned char)(1u<<((unsigned)node&7u)));
}
static void kls_snb_prefer_column(kls_snb_state *st,int64_t node) {
  if(st->row_column_preferred)atomic_fetch_or_explicit(st->row_column_preferred+(size_t)node/8,
    (unsigned char)(1u<<((unsigned)node&7u)),memory_order_relaxed);
}
static void kls_snb_retry_row(kls_snb_state *st,int64_t node) {
  if(st->row_column_preferred)atomic_fetch_and_explicit(st->row_column_preferred+(size_t)node/8,
    (unsigned char)~(1u<<((unsigned)node&7u)),memory_order_relaxed);
}
/* Duplication must fit alongside the native factor, which every certified
 * solve still reads. A miss keeps the existing column candidate unchanged. */
static int kls_snb_row_cache_fits(const kls_snb_state *st,size_t native_bytes,
    size_t capacity) {
  if(!st||st->nblocks!=1||!st->blocks||native_bytes>capacity||
      st->allocation_bytes>capacity-native_bytes)return 0;
  size_t used=native_bytes+st->allocation_bytes;
  const size_t metadata=kls_snb_row_preference_bytes(st);
  if(metadata==SIZE_MAX||metadata>capacity-used)return 0;
  used+=metadata;
  const kls_snb_block *b=st->blocks;
  for(int64_t j=0;j<b->snode_count;++j)if(b->snodes[j].width==16) {
    const kls_snb_snode *sn=b->snodes+j;
    if(sn->height<16||(size_t)sn->height>SIZE_MAX/(16*sizeof(double)))return 0;
    if(sn->column_panel)continue;
    const size_t extra=(size_t)sn->height*16*sizeof(double);
    if(extra>capacity-used)return 0;
    used+=extra;
  }
  return 1;
}
#ifdef KLS_SNB_HAVE_ROW_PANEL
static int kls_snb_row_inplace_storage_fits(const kls_snb_state *st) {
  if(!st||st->nblocks!=1||!st->blocks||st->nworkers<1||st->nworkers>64||
      !st->gemm_w)return 0;
  for(int t=0;t<st->nworkers;++t)if(!st->gemm_w[t])return 0;
  const kls_snb_block *b=st->blocks;
  for(int64_t j=0;j<b->snode_count;++j)if(b->snodes[j].width==16) {
    const int32_t H=b->snodes[j].height;
    if(H<16||(size_t)H>SIZE_MAX/(16*sizeof(double))||
        (size_t)H>st->gemm_scratch_bytes/2)return 0;
  }
  return 1;
}
static int kls_snb_row_inplace_cache_fits(const kls_snb_state *st,
    size_t native_bytes,size_t capacity) {
  return kls_snb_row_inplace_storage_fits(st)&&native_bytes<=capacity&&
    st->allocation_bytes<=capacity-native_bytes&&
    kls_snb_row_preference_bytes(st)<=capacity-native_bytes-st->allocation_bytes;
}
#endif
#ifdef KLS_SNB_HAVE_ROW_PANEL
static int kls_snb_row_compact_images(kls_snb_state *st,double deadline,size_t limit) {
  kls_snb_block *b=st->blocks;
    /* Count avoidable packs in existing descriptor padding. Negative values
     * are temporary reuse counts; no coefficient or new allocation is used.
     * Only mixed-width dense edges need a packed L-below rectangle. */
    for(int64_t j=0;j<b->snode_count;++j) {
      if((j&63)==0&&deadline>0&&kls_now_seconds()>=deadline)return 0;
      if(b->snodes[j].width==16)continue;
      const kls_snb_snode *target=b->snodes+j;
      for(int32_t e=0;e<target->edge_count;++e) {
        const kls_snb_edge *edge=b->edges+target->edges_off+e;
        kls_snb_snode *producer=b->snodes+edge->producer;
        if(producer->width==16&&!producer->column_panel&&producer->below>=8&&
            edge->sub_count>=st->wnarrow&&producer->column_ld>INT32_MIN)
          --producer->column_ld;
      }
    }
    for(;;) {
      kls_snb_snode *best=NULL;uint64_t best_score=0;
      if(deadline>0&&kls_now_seconds()>=deadline)return 0;
      for(int64_t j=0;j<b->snode_count;++j) {
        kls_snb_snode *producer=b->snodes+j;
        if(producer->column_ld>=0)continue;
        const size_t ld=(size_t)16+producer->below;
        if(ld>SIZE_MAX/(16*sizeof(double)))continue;
        const size_t bytes=ld*16*sizeof(double);
        if(st->allocation_bytes>limit||bytes>limit-st->allocation_bytes)continue;
        const uint64_t uses=(uint64_t)(-(int64_t)producer->column_ld);
        const uint64_t score=uses*(uint64_t)producer->below/ld;
        if(!best||score>best_score){best=producer;best_score=score;}
      }
      if(!best)break;
      const int32_t ld=16+best->below;
      best->column_panel=kls_snb_allocate(st,(size_t)ld*16,sizeof(double),0,deadline,limit);
      if(!best->column_panel)break; /* Optional image; existing row plan survives. */
      best->column_ld=ld;
    }
    for(int64_t j=0;j<b->snode_count;++j)if(b->snodes[j].column_ld<0)b->snodes[j].column_ld=0;
  return 1;
}
#endif
static int kls_snb_enable_row_layout_mode(kls_solver *solver,kls_snb_state *st,double deadline,
    size_t limit,int column_mirrors) {
#ifdef KLS_SNB_HAVE_ROW_PANEL
  if(!st||st->nblocks!=1)return 0;
  if(st->row_layout)return 1;
  if(!column_mirrors&&!kls_snb_row_inplace_storage_fits(st))return 0;
  kls_snb_block *b=st->blocks;int found=0;
  for(int64_t j=0;j<b->snode_count;++j)if(b->snodes[j].width==16)found=1;
  if(!found)return 1;
  int product_math=kls_snb_row_product_math();if(product_math<0)return 0;
  if(!st->row_column_preferred) {
    const size_t bytes=kls_snb_row_preference_bytes(st);
    if(bytes==SIZE_MAX)return 0;
    st->row_column_preferred=kls_snb_allocate(st,bytes,1,1,deadline,limit);
    if(!st->row_column_preferred)return 0;
    for(size_t j=0;j<bytes/sizeof(*st->row_column_preferred);++j)atomic_init(st->row_column_preferred+j,0);
  }
  /* Existing admitted plans retain their column images. A budgeted
   * mirrorless plan reads row producers directly and repairs in place. */
  if(column_mirrors)for(int64_t j=0;j<b->snode_count;++j)if(b->snodes[j].width==16) {
    kls_snb_snode *producer=b->snodes+j;
    producer->column_panel=kls_snb_allocate(st,(size_t)producer->height*16,sizeof(double),0,deadline,limit);
    if(!producer->column_panel)return 0;
  }
  /* The first row group owns 8*uw dead upper doubles after writeback.
   * A 16-column L image fits there iff uw >= 2*(16+below), keeping
   * destination and both live row groups disjoint. No new allocation. */
  if(!column_mirrors)for(int64_t j=0;j<b->snode_count;++j) {
    kls_snb_snode *producer=b->snodes+j;
    if(producer->width==16&&producer->below>=0&&
        producer->uw>=2*((int64_t)16+producer->below))
      {producer->column_panel=b->panels+producer->panel;producer->column_ld=16+producer->below;}
  }
  if(!column_mirrors&&!kls_snb_row_compact_images(st,deadline,limit))return 0;
  for(int64_t j=0;j<b->snode_count;++j) {
    if(deadline>0&&kls_now_seconds()>=deadline)return 0;
    const kls_snb_snode *sn=b->snodes+j;if(sn->width!=16)continue;
    for(int c=0;c<16;++c) {
      const int64_t k=sn->start+c;const UF_long g=b->k1+(UF_long)k;
      if(st->scatter_dst16) {
        const UF_long oldcol=solver->symbolic->Q[g];
        for(UF_long p=solver->col_ptr[oldcol];p<solver->col_ptr[oldcol+1];++p) {
          const int64_t row=(int64_t)st->scatter_dst16[p]-(int64_t)c*sn->height;
          if(row<0||row>=sn->height)return 0;
          const size_t at=kls_rpanel_at(sn->height,c,(int)row);
          if(at>UINT16_MAX)return 0;
          st->scatter_dst16[p]=(uint16_t)at;
        }
      } else if(st->scatter_dst) {
        const UF_long oldcol=solver->symbolic->Q[g];
        for(UF_long p=solver->col_ptr[oldcol];p<solver->col_ptr[oldcol+1];++p) {
          int64_t r=(int64_t)st->scatter_dst[p]-sn->panel-(int64_t)c*sn->height;
          if(r<0||r>=sn->height)return 0;
          uint64_t at=(uint64_t)sn->panel+kls_rpanel_at(sn->height,c,(int)r);
          if(at>UINT32_MAX)return 0;
          st->scatter_dst[p]=(uint32_t)at;
        }
      } else for(int64_t p=st->col_blk_ptr[g];p<st->col_blk_ptr[g+1];++p) {
        int64_t r=st->blk[p].dst-sn->panel-(int64_t)c*sn->height;
        if(r<0||r>=sn->height)return 0;
        st->blk[p].dst=sn->panel+(int64_t)kls_rpanel_at(sn->height,c,(int)r);
      }
      for(int64_t p=b->wb_col_off[k];p<b->wb_col_off[k+1];++p) {
        int64_t r=b->writeback[p]-(int64_t)c*sn->height;
        if(r<0||r>=sn->height)return 0;
        b->writeback[p]=(int32_t)kls_rpanel_at(sn->height,c,(int)r);
      }
    }
  }
  st->row_product_math=product_math;st->row_layout=1;st->transactional=1;
#else
  (void)solver;(void)st;(void)deadline;(void)limit;(void)column_mirrors;
#endif
  return 1;
}
static int kls_snb_enable_row_layout(kls_solver *solver,kls_snb_state *st,double deadline,size_t limit) {
  return kls_snb_enable_row_layout_mode(solver,st,deadline,limit,1);
}
static int kls_snb_add_resident_bytes(size_t *bytes,uint64_t count,size_t unit) {
  if(unit&&count>(SIZE_MAX-*bytes)/unit)return 0;
  *bytes+=(size_t)count*unit;return 1;
}
static int kls_snb_try_row_layout(kls_solver *s,kls_snb_state *st,
    double deadline,size_t capacity) {
  if(!s||!s->numeric||!st)return 0;
  if(st->row_layout)return 1;
  const size_t per_row=4*sizeof(UF_long)+2*sizeof(double);
  if((uint64_t)s->n>SIZE_MAX/per_row)return 1;
  size_t native_bytes=(size_t)s->n*per_row;
  for(UF_long b=0;b<s->symbolic->nblocks;++b) {
    if((uint64_t)s->numeric->LUsize[b]>(SIZE_MAX-native_bytes)/sizeof(double))return 1;
    native_bytes+=(size_t)s->numeric->LUsize[b]*sizeof(double);
  }
  const size_t off_unit=sizeof(double)+sizeof(UF_long);
  if((uint64_t)s->numeric->Offp[s->n]>(SIZE_MAX-native_bytes)/off_unit)return 1;
  native_bytes+=(size_t)s->numeric->Offp[s->n]*off_unit;
  /* The public frame still streams CSC input, its certified original
   * operator, and retained refactor/solve caches. Count allocated streams;
   * aliases of compact index storage are charged only once. */
#define KLS_SNB_RESIDENT(pointer,count,type) do { \
  if((pointer)&&!kls_snb_add_resident_bytes(&native_bytes,(uint64_t)(count),sizeof(type)))return 1; \
} while(0)
  KLS_SNB_RESIDENT(s->col_ptr,s->n+1,UF_long);
  KLS_SNB_RESIDENT(s->row_idx,s->nnz,UF_long);
  KLS_SNB_RESIDENT(s->input_to_csc,s->nnz,UF_long);
  /* One current coefficient stream is present even in direct-user mode. */
  if(!kls_snb_add_resident_bytes(&native_bytes,(uint64_t)s->nnz,sizeof(double)))return 1;
  KLS_SNB_RESIDENT(s->prepared_value_scale,s->nnz,double);
  KLS_SNB_RESIDENT(s->prepared_value_input_pos,s->nnz,UF_long);
  if(s->refactor_input_snapshot!=s->values)KLS_SNB_RESIDENT(s->refactor_input_snapshot,s->nnz,double);
  if(s->solve_refine_values!=s->values&&s->solve_refine_values!=s->refactor_input_snapshot)
    KLS_SNB_RESIDENT(s->solve_refine_values,s->nnz,double);
  KLS_SNB_RESIDENT(s->refactor_col_ptr,s->n+1,UF_long);
  if(s->refactor_row_idx!=s->row_idx)KLS_SNB_RESIDENT(s->refactor_row_idx,s->nnz,UF_long);
  KLS_SNB_RESIDENT(s->refactor_input_pos,s->nnz,UF_long);
  KLS_SNB_RESIDENT(s->refactor_row_idx32,s->nnz,int32_t);
  KLS_SNB_RESIDENT(s->refactor_input_pos32,s->nnz,int32_t);
  KLS_SNB_RESIDENT(s->refactor_input_oldrow32,s->nnz,int32_t);
  KLS_SNB_RESIDENT(s->refactor_input_user_pos32,s->nnz,int32_t);
  KLS_SNB_RESIDENT(s->refactor_scale_row_ptr,s->n+1,UF_long);
  KLS_SNB_RESIDENT(s->refactor_scale_input_pos32,s->nnz,int32_t);
  KLS_SNB_RESIDENT(s->refactor_scale_permute_values,s->nnz,double);
  KLS_SNB_RESIDENT(s->refactor_l_indices,s->n,void *);
  KLS_SNB_RESIDENT(s->refactor_u_indices,s->n,void *);
  KLS_SNB_RESIDENT(s->refactor_l_values,s->n,void *);
  KLS_SNB_RESIDENT(s->refactor_u_values,s->n,void *);
  KLS_SNB_RESIDENT(s->refactor_l_indices32,s->n,void *);
  KLS_SNB_RESIDENT(s->refactor_u_indices32,s->n,void *);
  KLS_SNB_RESIDENT(s->refactor_l_indices32_storage,s->refactor_l_indices32_count,int32_t);
  KLS_SNB_RESIDENT(s->refactor_u_indices32_storage,s->refactor_u_indices32_count,int32_t);
  KLS_SNB_RESIDENT(s->refactor_l_sorted_indices32_storage,s->refactor_l_indices32_count,int32_t);
  KLS_SNB_RESIDENT(s->refactor_l_sorted_pos32_storage,s->refactor_l_indices32_count,int32_t);
  KLS_SNB_RESIDENT(s->refactor_l_sorted_values_storage,s->refactor_l_indices32_count,double);
  KLS_SNB_RESIDENT(s->refactor_l_packed_storage,s->refactor_l_indices32_count,double);
  KLS_SNB_RESIDENT(s->refactor_u_packed_storage,s->refactor_u_indices32_count,double);
  if(!s->i32solve_indices_alias_refactor) {
    KLS_SNB_RESIDENT(s->i32solve_l,s->numeric->lnz,int32_t);
    KLS_SNB_RESIDENT(s->i32solve_u,s->numeric->unz,int32_t);
  }
  KLS_SNB_RESIDENT(s->i16solve_l,s->numeric->lnz,uint16_t);
  KLS_SNB_RESIDENT(s->i16solve_u,s->numeric->unz,uint16_t);
  KLS_SNB_RESIDENT(s->i32solve_loff,s->n+1,int64_t);
  KLS_SNB_RESIDENT(s->i32solve_uoff,s->n+1,int64_t);
  KLS_SNB_RESIDENT(s->i32solve_pnum,s->n,uint32_t);
  KLS_SNB_RESIDENT(s->i32solve_rhs_perm32,s->n,uint32_t);
  KLS_SNB_RESIDENT(s->i32solve_q,s->n,uint32_t);
  KLS_SNB_RESIDENT(s->i32solve_llen,s->n,uint32_t);
  KLS_SNB_RESIDENT(s->i32solve_ulen,s->n,uint32_t);
  KLS_SNB_RESIDENT(s->i32solve_loff32,s->n,uint32_t);
  KLS_SNB_RESIDENT(s->i32solve_uoff32,s->n,uint32_t);
  if(s->native_i16&&!kls_snb_add_resident_bytes(&native_bytes,s->native_i16->allocation_bytes,1))return 1;
  const private_original_operator *original=s->private_original;
  if(original) {
    KLS_SNB_RESIDENT(original->p,original->n+1,int64_t);
    KLS_SNB_RESIDENT(original->rows,original->nnz,int64_t);
    KLS_SNB_RESIDENT(original->values,original->nnz,double);
    KLS_SNB_RESIDENT(original->value_source,original->nnz,int64_t);
    if(original->layout) {
      if(!kls_snb_add_resident_bytes(&native_bytes,original->n+1,sizeof(UF_long))||
         !kls_snb_add_resident_bytes(&native_bytes,original->n,2*sizeof(UF_long))||
         !kls_snb_add_resident_bytes(&native_bytes,original->nnz,2*sizeof(UF_long)))return 1;
    }
  }
#undef KLS_SNB_RESIDENT
#ifdef KLS_SNB_ROW_RESIDENT_BYTES
  KLS_SNB_ROW_RESIDENT_BYTES(s,st,native_bytes);
#endif
  if(kls_snb_row_cache_fits(st,native_bytes,capacity))
    return kls_snb_enable_row_layout(s,st,deadline,capacity);
#ifdef KLS_SNB_HAVE_ROW_PANEL
  /* AUTO builds its required solve cache before this resident-ledger check. */
  if(kls_snb_row_inplace_cache_fits(st,native_bytes,capacity))
    return kls_snb_enable_row_layout_mode(s,st,deadline,capacity-native_bytes,0);
#endif
  return 1;
}
#ifdef KLS_SNB_HAVE_ROW_PANEL
static KLS_ALWAYS_INLINE double kls_rpanel_mul(double a,double b) {
  __asm__ volatile("mulsd %1, %0":"+x"(a):"x"(b));return a;
}
static KLS_ALWAYS_INLINE double kls_rpanel_add(double a,double b) {
  __asm__ volatile("addsd %1, %0":"+x"(a):"x"(b));return a;
}
static KLS_ALWAYS_INLINE double kls_rpanel_sub(double a,double b) {
  __asm__ volatile("subsd %1, %0":"+x"(a):"x"(b));return a;
}
/* Rare mixed-width consumers retain the original per-output producer order.
 * GCC also leaves its final one-to-three output columns uncontracted.
 * The first dense product and final subtraction remain separately rounded. */
/* Move only the producer's L-below rectangle into an owned worker buffer.
 * Each8x8 transpose consists entirely of bit moves, including exceptional
 * payloads. The caller checks space for both this block and the product T. */
static void kls_rpanel_pack_below(const double *R,int H,int first,int rows,double *L) {
  for(int group=0;group<2;++group) {
    int r=0;
    for(;r+7<rows;r+=8) {
      __m512d tile[8];
      for(int i=0;i<8;++i)tile[i]=_mm512_loadu_pd(R+(size_t)group*8*H+(size_t)(first+r+i)*8);
      kls_rpanel_transpose8(tile);
      for(int c=0;c<8;++c)_mm512_storeu_pd(L+(size_t)(group*8+c)*rows+r,tile[c]);
    }
    for(;r<rows;++r)for(int c=0;c<8;++c)
      memcpy(L+(size_t)(group*8+c)*rows+r,R+(size_t)group*8*H+(size_t)(first+r)*8+c,sizeof(double));
  }
}
/* Native LU already owns U-upper when repairing an unpublished panel.
 * Stage only L-diagonal and L-below before writing their row representation;
 * consumers never read the dead upper region. Exceptional bits are moved
 * without arithmetic, and the worker buffer remains exclusively owned. */
static int kls_rpanel_repair_lower_inplace(double *R,int H,int uw,int below,
    void *work,size_t work_bytes) {
  if(!R||!work||H<16||uw<0||uw>H-16||below!=H-uw-16||
      (size_t)H>SIZE_MAX/(16*sizeof(double)))return 0;
  const size_t rows=(size_t)H-(size_t)uw,bytes=rows*16*sizeof(double);
  if(work_bytes<bytes)return 0;
  const uintptr_t first=(uintptr_t)R,scratch=(uintptr_t)work;
  const size_t panel_bytes=(size_t)H*16*sizeof(double);
  if(first>UINTPTR_MAX-panel_bytes||scratch>UINTPTR_MAX-bytes||
      (first<scratch+bytes&&scratch<first+panel_bytes))return 0;
  double *L=work;
  for(int c=0;c<16;++c)memcpy(L+(size_t)c*rows,R+(size_t)c*H+uw,rows*sizeof(double));
  for(int group=0;group<2;++group) {
    size_t r=0;
    for(;r+7<rows;r+=8) {
      __m512d tile[8];
      for(int c=0;c<8;++c)tile[c]=_mm512_loadu_pd(L+(size_t)(group*8+c)*rows+r);
      kls_rpanel_transpose8(tile);
      for(int i=0;i<8;++i)_mm512_storeu_pd(R+(size_t)group*8*H+((size_t)uw+r+i)*8,tile[i]);
    }
    for(;r<rows;++r)for(int c=0;c<8;++c)
      memcpy(R+(size_t)group*8*H+((size_t)uw+r)*8+c,L+(size_t)(group*8+c)*rows+r,sizeof(double));
  }
#ifdef KLS_SNB_LOWER_REPAIR_EVENT
  KLS_SNB_LOWER_REPAIR_EVENT();
#endif
  return 1;
}
static void kls_rpanel_small_target_edge_to(const kls_snb_block *blk,
    const kls_snb_snode *sn,const kls_snb_edge *edge,int narrow,int product_math,
    double *scratch,size_t scratch_doubles,double *W) {
#ifdef KLS_RPANEL_TEST_SMALL_EDGE
  KLS_RPANEL_TEST_SMALL_EDGE();
#endif
  const kls_snb_snode *ps=blk->snodes+edge->producer;
  const int sc=edge->sub_count,ph=ps->below,H=sn->height;
  if(ps->column_panel&&(sc<narrow||ph==0||
      (scratch&&(size_t)sn->width<=scratch_doubles/(size_t)ph))) {
    /* Reuse the incumbent vector kernels and the owned worker scratch. */
    kls_snb_update_edge_to(blk,sn,edge,narrow,scratch,W);return;
  }
  const int32_t *sub=blk->subcols+edge->sub_off,*dst=blk->edge_dsts+edge->dst_off;
  const double *panel=blk->panels+ps->panel;
  for(int c=0;c<sn->width;++c) {
    double *wc=W+(size_t)c*H;
    for(int i=0;i<sc;++i) {
      double u=wc[edge->u_start+i];if(u==0.)continue;
      const double *lc=kls_rpanel_column(panel,ps->height,sub[i],1);
      for(int j=i+1;j<sc;++j)
        wc[edge->u_start+j]=fma(-lc[(size_t)(ps->uw+sub[j])*8],u,wc[edge->u_start+j]);
      if(sc<narrow)for(int r=0;r<ph;++r)if(dst[r]>=0)
        wc[dst[r]]=fma(-lc[(size_t)(ps->uw+ps->width+r)*8],u,wc[dst[r]]);
    }
  }
  if(sc<narrow)return;
  if(ph>=8&&scratch&&(size_t)(16+sn->width)<=scratch_doubles/(size_t)ph) {
#ifdef KLS_SNB_PACKED_ROW_EDGE_EVENT
    KLS_SNB_PACKED_ROW_EDGE_EVENT();
#endif
    double *T=scratch+(size_t)16*ph;
    kls_rpanel_pack_below(panel,ps->height,ps->uw+ps->width,ph,scratch);
    kls_snb_dense_product(scratch,ph,sub,sc,W+edge->u_start,H,sn->width,0,ph,T);
    for(int c=0;c<sn->width;++c)for(int r=0;r<ph;++r)if(dst[r]>=0) {
      double *wc=W+(size_t)c*H;
      wc[dst[r]]=kls_rpanel_sub(wc[dst[r]],T[(size_t)c*ph+r]);
    }
    return;
  }
#ifdef KLS_SNB_SCALAR_ROW_EDGE_EVENT
  KLS_SNB_SCALAR_ROW_EDGE_EVENT();
#endif
  for(int c=0;c<sn->width;++c) {
    double *wc=W+(size_t)c*H;
    for(int r=0;r<ph;++r) {
      const double *lc=kls_rpanel_column(panel,ps->height,sub[0],1);
      double t=kls_rpanel_mul(lc[(size_t)(ps->uw+ps->width+r)*8],wc[edge->u_start]);
      for(int i=1;i<sc;++i) {
        lc=kls_rpanel_column(panel,ps->height,sub[i],1);
        double l=lc[(size_t)(ps->uw+ps->width+r)*8],u=wc[edge->u_start+i];
        int category=c>=(sn->width&~3)?(r<(ph&~7)?3:4):
            (r<(ph&~15)?0:r<(ph&~7)?1:2);
        t=(product_math&(1<<category))?fma(l,u,t):kls_rpanel_add(t,kls_rpanel_mul(l,u));
      }
      if(dst[r]>=0)wc[dst[r]]=kls_rpanel_sub(wc[dst[r]],t);
    }
  }
}
static void kls_rpanel_small_target_edge(const kls_snb_block *blk,
    const kls_snb_snode *sn,const kls_snb_edge *edge,int narrow,int product_math,
    double *scratch,size_t scratch_doubles) {
  kls_rpanel_small_target_edge_to(blk,sn,edge,narrow,product_math,scratch,scratch_doubles,
    blk->panels+sn->panel);
}
#endif
