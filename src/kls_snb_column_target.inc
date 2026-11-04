/* Keep the qualified column executor separate from optional row repair.
 * Mixing repair addressing into its writeback loops prevents the incumbent
 * compiler specialization and enlarges the hot instruction footprint.
 * The edge and target arithmetic below is copied from qualified v51. */
static KLS_ALWAYS_INLINE void kls_snb_column_update_edge(
  const kls_snb_block *blk, const kls_snb_snode *sn,
  const kls_snb_edge *edge, int wnarrow,
  double *gemm_scratch) {
  const kls_snb_snode *ps = &blk->snodes[edge->producer];
  const int32_t w = sn->width;
  const int32_t H = sn->height;
  double *W = blk->panels + sn->panel;
  const int32_t ph = ps->below;
  const int32_t sc = edge->sub_count;
  const int32_t *sub = blk->subcols + edge->sub_off;
  const int32_t *dsts = blk->edge_dsts + edge->dst_off;
  const int32_t ustart = edge->u_start;
  const int32_t pld = ps->height;
  const double *Ldiag = blk->panels + ps->panel + ps->uw;
  const double *Lbelow = Ldiag + ps->width;
  if (sc < wnarrow) {
    const int32_t first=kls_snb_narrow_columns(W,H,w,Ldiag,Lbelow,pld,sub,sc,dsts,ph,ustart);
    for (int32_t c = first; c < w; ++c) {
      double *Wc = W + (size_t)c * H;
      for (int32_t i = 0; i < sc; ++i) {
        const double u = Wc[ustart + i];
        if (u == 0.0) continue;
        {
          const double *lcol = Ldiag + (size_t)sub[i] * pld;
          for (int32_t i2 = i + 1; i2 < sc; ++i2) {
            Wc[ustart + i2] -= lcol[sub[i2]] * u;
          }
          const double *lb = Lbelow + (size_t)sub[i] * pld;
          for (int32_t r = 0; r < ph; ++r) {
            const int32_t dst = dsts[r];
            if (dst >= 0) {
              Wc[dst] -= lb[r] * u;
            }
          }
        }
      }
    }
    return;
  }
  /* This edge's producer below rows follow all its diagonal columns.
     Its scatter therefore cannot overwrite this upper block. Other edges
     may update it, so preserve edge order and finish the product first. */
  double *B = W + ustart;
  /* The builder emits strictly increasing producer-column subsets. For a
     contiguous subset, direct row loads avoid gathers without changing the
     triangular update order. Gapped subsets retain the indexed kernel. */
  const int contiguous = sub[sc-1]-sub[0] == sc-1;
  int32_t first_column=0;
#if defined(__AVX512F__) && defined(__AVX512VL__) && FLT_EVAL_METHOD == 0
  if(sc==16 && w>=4 && ph>0 && contiguous && (kls_snb_owned_masked_environment || fegetround()==FE_TONEAREST)) {
    unsigned mx=0x1f80u;unsigned short cw=63u;
    if(!kls_snb_owned_masked_environment) {
      __asm__ volatile("stmxcsr %0":"=m"(mx)::"memory");
      __asm__ volatile("fnstcw %0":"=m"(cw));
    }
    if((mx&0x7f80u)==0x1f80u && (cw&63u)==63u)
      for(;first_column+3<w;first_column+=4)
        kls_snb_register_trsm_16(B+(size_t)first_column*H,H,Ldiag,pld,sub[0]);
  }
#endif
  for (int32_t i = 0; i < sc; ++i) {
    {
      const double *lcol = Ldiag + (size_t)sub[i] * pld;
      for (int32_t c = first_column; c < w; ++c) {
        double *Bc = B + (size_t)c * H;
        const double u = Bc[i];
        if (u == 0.0) continue;
        if (contiguous) {
          const double *selected = lcol + sub[0];
#if defined(__AVX512F__)
          /* Vectorize rows only, preserving producer/target-column order.
             Explicit FMA gives full vectors and tails the same rounding.
             Mask arithmetic too: inactive lanes must not evaluate 0*Inf. */
          const __m512d vu = _mm512_set1_pd(u);
          int32_t i2 = i + 1;
          for (; i2 + 7 < sc; i2 += 8) {
            _mm512_storeu_pd(Bc+i2, _mm512_fnmadd_pd(
              _mm512_loadu_pd(selected+i2), vu, _mm512_loadu_pd(Bc+i2)));
          }
          if (i2 < sc) {
            const __mmask8 mask = (__mmask8)((1u << (sc-i2))-1u);
            _mm512_mask_storeu_pd(Bc+i2,mask,
              _mm512_mask3_fnmadd_pd(_mm512_maskz_loadu_pd(mask,selected+i2),
                vu, _mm512_maskz_loadu_pd(mask,Bc+i2), mask));
          }
#else
          for (int32_t i2 = i + 1; i2 < sc; ++i2)
            Bc[i2] -= selected[i2] * u;
#endif
        } else {
          for (int32_t i2 = i + 1; i2 < sc; ++i2)
            Bc[i2] -= lcol[sub[i2]] * u;
        }
      }
    }
  }
  if (ph == 0) return;
  double *T = gemm_scratch;
  kls_snb_dense_product(Lbelow, pld, sub, sc, B, H, w, 0, ph, T);
  /* Share destination decoding across the same four-column width used
     by the product. Each column retains producer-row subtraction order. */
  int32_t c = 0;
  for (; c + 3 < w; c += 4) {
    double *w0 = W + (size_t)c * H;
    const double *t0 = T + (size_t)c * ph;
    for (int32_t r = 0; r < ph; ++r) {
      const int32_t dst = dsts[r];
      if (dst >= 0) {
        w0[dst] -= t0[r];
        w0[(size_t)H + dst] -= t0[(size_t)ph + r];
        w0[2u * H + dst] -= t0[2u * ph + r];
        w0[3u * H + dst] -= t0[3u * ph + r];
      }
    }
  }
  for (; c < w; ++c)
    for (int32_t r = 0; r < ph; ++r)
      if (dsts[r] >= 0) W[(size_t)c * H + dsts[r]] -= T[(size_t)c * ph + r];
}

/* Process one supernode: panel compute (memset, A-scatter program, edge
 * updates, fused no-pivot getrf) followed by an immediate writeback of
 * its member columns into the packed numeric, so consumers that read
 * LUbx (including a future pipelined tail) always see finished columns.
 * Thread-safe across distinct snodes of one level. */
static __attribute__((noinline)) void kls_snb_column_process_target(kls_solver *solver, kls_snb_state *st,
                                   kls_snb_block *blk, int64_t s,
                                   const double *az, const double *rs,
                                   int scale, double *gemm_scratch, int pipelined) {
  const UF_long k1 = blk->k1;
  double *panels = blk->panels;
  double *offx = (double *)solver->numeric->Offx;
  double *udiag = (double *)solver->numeric->Udiag;
  const int wnarrow = st->wnarrow;
  const kls_snb_snode *sn = &blk->snodes[s];
  const int32_t w = sn->width;
  const int32_t H = sn->height;
  const int32_t uw = sn->uw;
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
    kls_snb_column_update_edge(blk, sn, &edges[e], wnarrow,
                        gemm_scratch);
  }
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

