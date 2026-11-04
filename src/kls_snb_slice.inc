/* Column slices retain producer order and four-column kernel tiles. */
static void kls_snb_process_target_slice(kls_solver *solver, kls_snb_state *st,
                                   kls_snb_block *blk, int64_t s,
                                   const double *az, const double *rs,
                                   int scale, double *gemm_scratch, int pipelined, int begin, int end) {
  const UF_long k1 = blk->k1;
  double *panels = blk->panels;
  double *offx = (double *)solver->numeric->Offx;
  const int wnarrow = st->wnarrow;
  const kls_snb_snode *original=&blk->snodes[s];
  if(begin<0 || begin>=end || end>original->width)return;
  kls_snb_snode slice=*original;
  slice.start+=begin;slice.panel+=(int64_t)begin*slice.height;slice.width=end-begin;
  const kls_snb_snode *sn=&slice;
  const int32_t w = sn->width;
  const int32_t H = sn->height;
  double *W = panels + sn->panel;
  memset(W, 0, (size_t)H * (size_t)w * sizeof(double));
  for (int32_t c = 0; c < w; ++c) {
    const UF_long gk = k1 + sn->start + c;
    for (int64_t p = st->col_off_ptr[gk]; p < st->col_off_ptr[gk + 1u];
         ++p) {
      offx[st->off[p].dst] =
        scale > 0 ? az[st->off[p].src] / rs[st->off[p].oldrow]
                  : az[st->off[p].src];
    }
    for (int64_t p = st->col_blk_ptr[gk]; p < st->col_blk_ptr[gk + 1u];
         ++p) {
      panels[st->blk[p].dst] =
        scale > 0 ? az[st->blk[p].src] / rs[st->blk[p].oldrow]
                  : az[st->blk[p].src];
    }
  }
  const kls_snb_edge *edges = blk->edges + sn->edges_off;
  for (int32_t e = 0; e < sn->edge_count; ++e) {
    if (pipelined) {
      const int32_t sc = edges[e].sub_count;
      const int32_t *sub = blk->subcols + edges[e].sub_off;
      /* the producer publishes each getrf column right after its
         division; this edge only reads columns <= sub[sc-1] */
      const int32_t need = sub[sc - 1] + 1;
      while (atomic_load_explicit(&blk->published[edges[e].producer],
                                  memory_order_acquire) < need) {
        if (atomic_load_explicit(&st->failed,memory_order_acquire)) return;
        kls_cpu_relax();
      }
    }
    kls_snb_update_edge(blk, sn, &edges[e], wnarrow,
                        gemm_scratch);
  }
}

static void kls_snb_finalize_target(kls_solver *solver,kls_snb_state *st,
    kls_snb_block *blk,int64_t s) {
  const UF_long k1=blk->k1;
  double *panels=blk->panels,*udiag=solver->numeric->Udiag;
  const kls_snb_snode *sn=&blk->snodes[s];
  const int32_t w=sn->width,H=sn->height,uw=sn->uw;
  double *W=panels+sn->panel;
  const int pipelined=0;
  /* fused no-pivot getrf/TRSM */
  {
    double *diag = W + uw;

    for (int32_t i = 0; i < w; ++i) {
      double *ci = diag + (size_t)i * H;
      const double pivot = ci[i];
      udiag[k1 + sn->start + i] = pivot;
      if (pivot == 0.0) {
        kls_snb_record_singular(st, solver, k1 + sn->start + i);
      }
      const int32_t hrem = H - uw - i - 1;
      double *lsub = ci + i + 1;
      for (int32_t r = 0; r < hrem; ++r) lsub[r] /= pivot;

      if (pipelined) {
        atomic_store_explicit(&blk->published[s], i + 1,
                              memory_order_release);
      }
      for (int32_t j = i + 1; j < w; ++j) {
        double *cj = diag + (size_t)j * H;
        const double u = cj[i];
        if (u == 0.0) continue;
        double *dstv = cj + i + 1;
        for (int32_t r = 0; r < hrem; ++r) {
          dstv[r] -= lsub[r] * u;
        }
      }
    }
  }
  /* immediate writeback of this snode's columns */
  {
    trilinos_klu_l_numeric *numeric = solver->numeric;
    const UF_long *lip = numeric->Lip + k1;
    const UF_long *llen = numeric->Llen + k1;
    const UF_long *uip = numeric->Uip + k1;
    const UF_long *ulen = numeric->Ulen + k1;
    double *lu = (double *)numeric->LUbx[blk - st->blocks];
    const double *base = panels + sn->panel;
    for (int64_t k = sn->start; k < sn->start + sn->width; ++k) {
      const int32_t *wb = blk->writeback + blk->wb_col_off[k];
      int64_t wc = 0;
      UF_long *ri; double *rv; UF_long rl;
      kls_klu_get_pointer(lu, lip, llen, (UF_long)k, &ri, &rv, &rl);
      for (UF_long p = 0; p < rl; ++p) {
        rv[p] = base[wb[wc++]];
      }
      kls_klu_get_pointer(lu, uip, ulen, (UF_long)k, &ri, &rv, &rl);
      for (UF_long p = 0; p < rl; ++p) {
        rv[p] = base[wb[wc++]];
      }
    }
  }
}

