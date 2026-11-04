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
/* Duplication must fit alongside the native factor, which every certified
 * solve still reads. A miss keeps the existing column candidate unchanged. */
static int kls_snb_row_cache_fits(const kls_snb_state *st,size_t native_bytes,
    size_t capacity) {
  if(!st||st->nblocks!=1||native_bytes>capacity||
      st->allocation_bytes>capacity-native_bytes)return 0;
  size_t used=native_bytes+st->allocation_bytes;
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
static int kls_snb_enable_row_layout(kls_snb_state *st,double deadline,size_t limit) {
#ifdef KLS_SNB_HAVE_ROW_PANEL
  if(!st||st->nblocks!=1)return 0;
  if(st->row_layout)return 1;
  kls_snb_block *b=st->blocks;int found=0;
  for(int64_t j=0;j<b->snode_count;++j)if(b->snodes[j].width==16)found=1;
  if(!found)return 1;
  int product_math=kls_snb_row_product_math();if(product_math<0)return 0;
  /* Every row producer has an ordinary column image, also used to repair
   * a rare underflowing target without discarding the whole owned job. */
  for(int64_t j=0;j<b->snode_count;++j)if(b->snodes[j].width==16) {
    kls_snb_snode *producer=b->snodes+j;
    producer->column_panel=kls_snb_allocate(st,(size_t)producer->height*16,sizeof(double),0,deadline,limit);
    if(!producer->column_panel)return 0;
  }
  for(int64_t j=0;j<b->snode_count;++j) {
    if(deadline>0&&kls_now_seconds()>=deadline)return 0;
    const kls_snb_snode *sn=b->snodes+j;if(sn->width!=16)continue;
    for(int c=0;c<16;++c) {
      const int64_t k=sn->start+c;const UF_long g=b->k1+(UF_long)k;
      for(int64_t p=st->col_blk_ptr[g];p<st->col_blk_ptr[g+1];++p) {
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
  (void)st;(void)deadline;(void)limit;
#endif
  return 1;
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
  if(!kls_snb_row_cache_fits(st,native_bytes,capacity))return 1;
  return kls_snb_enable_row_layout(st,deadline,capacity);
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
static void kls_rpanel_small_target_edge(const kls_snb_block *blk,
    const kls_snb_snode *sn,const kls_snb_edge *edge,int narrow,int product_math,double *scratch) {
#ifdef KLS_RPANEL_TEST_SMALL_EDGE
  KLS_RPANEL_TEST_SMALL_EDGE();
#endif
  const kls_snb_snode *ps=blk->snodes+edge->producer;
  if(ps->column_panel) {
    /* Reuse the incumbent vector kernels and the owned worker scratch. */
    kls_snb_update_edge(blk,sn,edge,narrow,scratch);return;
  }
  const int sc=edge->sub_count,ph=ps->below,H=sn->height;
  const int32_t *sub=blk->subcols+edge->sub_off,*dst=blk->edge_dsts+edge->dst_off;
  const double *panel=blk->panels+ps->panel;double *W=blk->panels+sn->panel;
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
    if(sc>=narrow)for(int r=0;r<ph;++r) {
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
#endif
