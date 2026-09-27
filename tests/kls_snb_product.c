/* Deterministic private-kernel coverage, independent of policy timing and
   circuit matrices. The test executable exports no additional library API. */
#include "../src/kls.c"

static KLS_ALWAYS_INLINE void reference_edge(
  const kls_snb_block *blk, const kls_snb_snode *sn,
  const kls_snb_edge *edge, int wnarrow,
  double *ub_scratch, double *gemm_scratch) {
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
    for (int32_t c = 0; c < w; ++c) {
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
  double *B = ub_scratch;
  for (int32_t c = 0; c < w; ++c) {
    const double *Wc = W + (size_t)c * H + ustart;
    double *Bc = B + (size_t)c * sc;
    for (int32_t i = 0; i < sc; ++i) Bc[i] = Wc[i];
  }
  for (int32_t i = 0; i < sc; ++i) {
    {
      const double *lcol = Ldiag + (size_t)sub[i] * pld;
      for (int32_t c = 0; c < w; ++c) {
        double *Bc = B + (size_t)c * sc;
        const double u = Bc[i];
        if (u == 0.0) continue;
        for (int32_t i2 = i + 1; i2 < sc; ++i2) {
          Bc[i2] -= lcol[sub[i2]] * u;
        }
      }
    }
  }
  for (int32_t c = 0; c < w; ++c) {
    double *Wc = W + (size_t)c * H + ustart;
    const double *Bc = B + (size_t)c * sc;
    for (int32_t i = 0; i < sc; ++i) Wc[i] = Bc[i];
  }
  if (ph == 0) return;
  double *T = gemm_scratch;
  kls_snb_dense_product(Lbelow, pld, sub, sc, B, sc, w, 0, ph, T);
  for (int32_t c = 0; c < w; ++c)
    for (int32_t r = 0; r < ph; ++r)
      if (dsts[r] >= 0) W[(size_t)c * H + dsts[r]] -= T[(size_t)c * ph + r];
}

static int check_edges(void) {
  const int widths[] = {1, 3, 4, 5, 7, 8, 9, 31, 32};
  const int heights[] = {0, 1, 2, 3, 4, 5, 6, 7, 17, 65};
  double panels[16384], expected[16384], ub[1024], tmp[4096];
  int32_t sub[32], dst[65];
  for (int stride = 1; stride <= 2; ++stride)
   for (int offset = 0; offset <= 1; ++offset)
    for (int sc = 1; sc <= 32; ++sc)
    for (unsigned wi = 0; wi < sizeof(widths)/sizeof(widths[0]); ++wi)
      for (unsigned hi = 0; hi < sizeof(heights)/sizeof(heights[0]); ++hi)
        for (int holes = 0; holes < 3; ++holes) {
          const int w = widths[wi], ph = heights[hi];
          kls_snb_snode nodes[2] = {{0}};
          nodes[0].width = 64; nodes[0].below = ph; nodes[0].height = 64 + ph;
          nodes[1].width = w; nodes[1].height = 129; nodes[1].panel = 9000;
          kls_snb_edge edge = {0}; edge.sub_count = sc; edge.u_start = 1;
          kls_snb_block blk = {0};
          blk.snodes = nodes; blk.subcols = sub; blk.edge_dsts = dst;
          for (int i = 0; i < sc; ++i) sub[i] = stride * i + offset;
          for (int row = 0; row < ph; ++row)
            dst[row] = (holes == 2 || (holes == 1 && row % 3 == 0)) ? -1 : 63 + row;
          for (int i = 0; i < 16384; ++i)
            panels[i] = expected[i] = ((i * 17) % 29 - 14) / 100.0;
          blk.panels = expected;
          reference_edge(&blk, &nodes[1], &edge, 8, ub, tmp);
          blk.panels = panels;
          kls_snb_update_edge(&blk, &nodes[1], &edge, 8, tmp);
          for (int i = 0; i < 16384; ++i)
            if (!isfinite(panels[i]) || fabs(panels[i] - expected[i]) > 1e-12) {
              fprintf(stderr, "edge mismatch sc=%d w=%d ph=%d holes=%d index=%d\n",
                      sc, w, ph, holes, i);
              return 1;
            }
        }
  return 0;
}

int main(void) {
  if (check_edges()) return 1;
  const int rows[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 16, 17, 63, 64, 65};
  const int widths[] = {1, 3, 4, 5, 7, 8, 9, 31, 32};
  const int counts[] = {1, 2, 7, 8, 17, 32};
  double L[73 * 37], B[66 * 32], scratch[65 * 32 + 2];
  int32_t sub[32];
  for (int p = 0; p < 73 * 37; ++p) L[p] = ((p * 17) % 29 - 14) / 11.0;
  for (int p = 0; p < 66 * 32; ++p) B[p] = ((p * 13) % 23 - 11) / 7.0;
  for (int i = 0; i < 32; ++i) sub[i] = (i * 7) % 37;
  for (unsigned a = 0; a < sizeof(rows) / sizeof(rows[0]); ++a)
    for (unsigned b = 0; b < sizeof(widths) / sizeof(widths[0]); ++b)
      for (unsigned d = 0; d < sizeof(counts) / sizeof(counts[0]); ++d)
       for (int pad = 0; pad <= 33; pad += 11) {
        const int nr = rows[a], w = widths[b], sc = counts[d];
        double *T = scratch + 1;
        scratch[0] = T[nr * w] = 12345.0;
        const int bld = sc + pad;
        kls_snb_dense_product(L, 73, sub, sc, B, bld, w, 3, 3 + nr, T);
        if (scratch[0] != 12345.0 || T[nr * w] != 12345.0) return 1;
        for (int c = 0; c < w; ++c)
          for (int r = 0; r < nr; ++r) {
            long double expected = 0.0;
            for (int i = 0; i < sc; ++i)
              expected += (long double)L[sub[i] * 73 + 3 + r] * B[c * bld + i];
            if (!isfinite(T[c * nr + r]) ||
                fabsl(T[c * nr + r] - expected) > 1e-12L) {
              fprintf(stderr, "SNB product mismatch rows=%d width=%d count=%d\n", nr, w, sc);
              return 1;
            }
          }
      }
  return 0;
}
