/* Frozen ordinary edge executor before register TRSM integration. */
static __attribute__((noinline)) void reference_register_edge(
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
  for (int32_t i = 0; i < sc; ++i) {
    {
      const double *lcol = Ldiag + (size_t)sub[i] * pld;
      for (int32_t c = 0; c < w; ++c) {
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
