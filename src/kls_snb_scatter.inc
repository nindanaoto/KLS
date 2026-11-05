/* The single-block map borrows no storage: source ranges and scaling rows
 * are read from the solver's authoritative CSC geometry for this epoch. The
 * plan owns only its destination map. Multi-block plans keep the general
 * off-diagonal and reordered-column program. */
static KLS_ALWAYS_INLINE void kls_snb_scatter_input_column(kls_solver *solver,
    const kls_snb_state *st,const kls_snb_block *blk,const kls_snb_snode *sn,
    UF_long g,const double *az,const double *rs,int scale,double *offx,double *to) {
  if(st->scatter_dst) {
    const UF_long oldcol=solver->symbolic->Q[g];
    for(UF_long p=solver->col_ptr[oldcol];p<solver->col_ptr[oldcol+1];++p) {
      const double value=scale>0?az[p]/rs[solver->row_idx[p]]:az[p];
      const uint32_t dst=st->scatter_dst[p];
      if(to)to[kls_snb_row_to_column_index(sn->height,(int64_t)dst-sn->panel)]=value;
      else blk->panels[dst]=value;
    }
    return;
  }
  for(int64_t p=st->col_off_ptr[g];p<st->col_off_ptr[g+1];++p)
    offx[st->off[p].dst]=scale>0?az[st->off[p].src]/rs[st->off[p].oldrow]:az[st->off[p].src];
  for(int64_t p=st->col_blk_ptr[g];p<st->col_blk_ptr[g+1];++p) {
    const double value=scale>0?az[st->blk[p].src]/rs[st->blk[p].oldrow]:az[st->blk[p].src];
    if(to)to[kls_snb_row_to_column_index(sn->height,st->blk[p].dst-sn->panel)]=value;
    else blk->panels[st->blk[p].dst]=value;
  }
}
