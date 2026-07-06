/* Type-parameterized supernodal refactor body (round 7: fp32 panels).
 * Included twice from supernodal_bench.c with:
 *   SNB_REAL    element type of the panels (double or float)
 *   SNB_NAME    name of the generated function
 *   SNB_PANELS  member of snb_block holding the panel storage
 * Scratch buffers (st->ub, st->gemm_temp) are allocated as doubles and
 * reused via cast; float needs strictly less space. */

#ifdef SNB_L_MIRROR
#define SNB_LPTR float
#else
#define SNB_LPTR SNB_REAL
#endif

#define SNB_GLUE2(a, b) a##b
#define SNB_GLUE(a, b) SNB_GLUE2(a, b)
#define SNB_TARGET_FN SNB_GLUE(SNB_NAME, _target)
#define SNB_WORKER_FN SNB_GLUE(SNB_NAME, _worker)

static void SNB_TARGET_FN(snb_state *st, snb_block *blk, int64_t s,
                          double *ub_scratch, double *gemm_scratch) {
  const double *rs = st->rs;
  const double *az = st->a->values;
  const scatter_program *sp = &st->sp;
  const int wnarrow = st->wnarrow;
  const UF_long k1 = blk->k1;
  SNB_REAL *panels = blk->SNB_PANELS;
  const int prof = st->profile && st->threads <= 1;

  {
    {
      const snb_snode *sn = &blk->snodes[s];
      const int32_t w = sn->width;
      const int32_t H = sn->height;
      const int32_t uw = sn->uw;
      SNB_REAL *W = panels + sn->panel;
      uint64_t t0 = prof ? snb_tsc() : 0;
      memset(W, 0, (size_t)H * (size_t)w * sizeof(SNB_REAL));
      /* scatter A columns (program dst = local panel index incl. col) */
      for (int32_t c = 0; c < w; ++c) {
        const UF_long gk = k1 + sn->start + c;
        for (int64_t p = sp->col_off_ptr[gk]; p < sp->col_off_ptr[gk + 1];
             ++p) {
          st->offx[sp->off[p].dst] =
            st->scale > 0 ? az[sp->off[p].src] / rs[sp->off[p].oldrow]
                          : az[sp->off[p].src];
        }
        for (int64_t p = sp->col_blk_ptr[gk]; p < sp->col_blk_ptr[gk + 1];
             ++p) {
          W[sp->blk[p].dst] = (SNB_REAL)(
            st->scale > 0 ? az[sp->blk[p].src] / rs[sp->blk[p].oldrow]
                          : az[sp->blk[p].src]);
        }
      }
      if (prof) {
        const uint64_t t1 = snb_tsc();
        st->phase_cycles[0] += t1 - t0;
        t0 = t1;
      }
      /* edges, ascending producer */
      const snb_edge *edges = blk->edges + sn->edges_off;
      for (int32_t e = 0; e < sn->edge_count; ++e) {
        const snb_snode *ps = &blk->snodes[edges[e].producer];
        const int32_t pw = ps->width;
        const int32_t ph = ps->below;
        const int32_t sc = edges[e].sub_count;
        const int32_t *sub = blk->subcols + edges[e].sub_off;
        const int32_t *dsts = blk->edge_dsts + edges[e].dst_off;
        const int32_t ustart = edges[e].u_start;
        const int32_t pld = ps->height;
        const int32_t puw = ps->uw;
#ifdef SNB_L_MIRROR
        const float *Lp = blk->lmirror + ps->panel;
        const float *Ldiag = Lp + puw;            /* pw x pw diag block */
        const float *Lbelow = Lp + puw + pw;      /* ph x pw below block */
#else
        const SNB_REAL *Lp = panels + ps->panel;  /* producer panel */
        const SNB_REAL *Ldiag = Lp + puw;         /* pw x pw diag block */
        const SNB_REAL *Lbelow = Lp + puw + pw;   /* ph x pw below block */
#endif
        if (sc < wnarrow) {
          /* hybrid scalar path: register u, immediate scatter, u==0
           * skip; precomputed dsts replace the map (round 3) */
          for (int32_t c = 0; c < w; ++c) {
            SNB_REAL *Wc = W + (size_t)c * H;
            for (int32_t i = 0; i < sc; ++i) {
              const SNB_REAL u = Wc[ustart + i];
              if (u == (SNB_REAL)0) continue;
              const SNB_LPTR *lcol = Ldiag + (size_t)sub[i] * pld;
              for (int32_t i2 = i + 1; i2 < sc; ++i2) {
                Wc[ustart + i2] -= (SNB_REAL)lcol[sub[i2]] * u;
              }
              const SNB_LPTR *lb = Lbelow + (size_t)sub[i] * pld;
              for (int32_t r = 0; r < ph; ++r) {
                const int32_t dst = dsts[r];
                if (dst >= 0) {
                  Wc[dst] -= (SNB_REAL)lb[r] * u;
                }
              }
            }
          }
          if (prof) {
            const uint64_t t1 = snb_tsc();
            st->phase_cycles[1] += t1 - t0;
            t0 = t1;
          }
          continue;
        }
        /* gather UB (sc x w), TRSM restricted to the touched producer
         * columns (untouched ones carry exact zeros), writeback */
        SNB_REAL *B = (SNB_REAL *)ub_scratch;  /* column-major sc x w */
        for (int32_t c = 0; c < w; ++c) {
          const SNB_REAL *Wc = W + (size_t)c * H + ustart;
          SNB_REAL *Bc = B + (size_t)c * sc;
          for (int32_t i = 0; i < sc; ++i) Bc[i] = Wc[i];
        }
        for (int32_t i = 0; i < sc; ++i) {
          const SNB_LPTR *lcol = Ldiag + (size_t)sub[i] * pld;
          for (int32_t c = 0; c < w; ++c) {
            SNB_REAL *Bc = B + (size_t)c * sc;
            const SNB_REAL u = Bc[i];
            if (u == (SNB_REAL)0) continue;
            for (int32_t i2 = i + 1; i2 < sc; ++i2) {
              Bc[i2] -= (SNB_REAL)lcol[sub[i2]] * u;
            }
          }
        }
        for (int32_t c = 0; c < w; ++c) {
          SNB_REAL *Wc = W + (size_t)c * H + ustart;
          const SNB_REAL *Bc = B + (size_t)c * sc;
          for (int32_t i = 0; i < sc; ++i) Wc[i] = Bc[i];
        }
        if (prof) {
          const uint64_t t1 = snb_tsc();
          st->phase_cycles[2] += t1 - t0;
          t0 = t1;
        }
        if (ph == 0) continue;
        /* block path: GEMM into row-major temp (long contiguous streams,
         * independent FMAs -- measured faster than a fused accumulate),
         * then row scatter through the precomputed dsts */
        SNB_REAL *Brow =
          (SNB_REAL *)ub_scratch + (size_t)sc * w; /* row-major mirror */
        for (int32_t c = 0; c < w; ++c) {
          const SNB_REAL *Bc = B + (size_t)c * sc;
          for (int32_t i = 0; i < sc; ++i) {
            Brow[(size_t)i * w + c] = Bc[i];
          }
        }
        SNB_REAL *T = (SNB_REAL *)gemm_scratch;
        {
          const SNB_LPTR *lb = Lbelow + (size_t)sub[0] * pld;
          const SNB_REAL *br = Brow;
          for (int32_t r = 0; r < ph; ++r) {
            const SNB_REAL l = (SNB_REAL)lb[r];
            SNB_REAL *tr = T + (size_t)r * w;
            for (int32_t c = 0; c < w; ++c) tr[c] = l * br[c];
          }
          for (int32_t i = 1; i < sc; ++i) {
            const SNB_LPTR *lbi = Lbelow + (size_t)sub[i] * pld;
            const SNB_REAL *bri = Brow + (size_t)i * w;
            for (int32_t r = 0; r < ph; ++r) {
              const SNB_REAL l = (SNB_REAL)lbi[r];
              SNB_REAL *tr = T + (size_t)r * w;
              for (int32_t c = 0; c < w; ++c) tr[c] += l * bri[c];
            }
          }
        }
        for (int32_t r = 0; r < ph; ++r) {
          const int32_t dst = dsts[r];
          const SNB_REAL *tr = T + (size_t)r * w;
          if (dst >= 0) {
            SNB_REAL *Wd = W + dst;
            for (int32_t c = 0; c < w; ++c) {
              Wd[(size_t)c * H] -= tr[c];
            }
          }
        }
        if (prof) {
          const uint64_t t1 = snb_tsc();
          st->phase_cycles[3] += t1 - t0;
          t0 = t1;
        }
      }
      /* fused no-pivot getrf/TRSM on rows [uw, H) */
      {
        SNB_REAL *diag = W + uw;
        for (int32_t i = 0; i < w; ++i) {
          SNB_REAL *ci = diag + (size_t)i * H;
          const SNB_REAL pivot = ci[i];
          st->udiag[k1 + sn->start + i] = (double)pivot;
          if (pivot == (SNB_REAL)0) st->status = 1;
          const int32_t hrem = H - uw - i - 1;
          SNB_REAL *lsub = ci + i + 1;
          /* division (not reciprocal) to match KLU's roundoff */
          for (int32_t r = 0; r < hrem; ++r) lsub[r] /= pivot;
          for (int32_t j = i + 1; j < w; ++j) {
            SNB_REAL *cj = diag + (size_t)j * H;
            const SNB_REAL u = cj[i];
            if (u == (SNB_REAL)0) continue;
            SNB_REAL *dstv = cj + i + 1;
            for (int32_t r = 0; r < hrem; ++r) {
              dstv[r] -= lsub[r] * u;
            }
          }
        }
      }
#ifdef SNB_L_MIRROR
      /* mirror the freshly factored L region to fp32 for consumers */
      {
        float *M = blk->lmirror + sn->panel;
        const SNB_REAL *Wd = W + uw;
        float *Md = M + uw;
        for (int32_t i = 0; i < w; ++i) {
          const SNB_REAL *ci = Wd + (size_t)i * H + i + 1;
          float *mi = Md + (size_t)i * H + i + 1;
          const int32_t hrem = H - uw - i - 1;
          for (int32_t r = 0; r < hrem; ++r) mi[r] = (float)ci[r];
        }
      }
#endif
      if (prof) {
        const uint64_t t1 = snb_tsc();
        st->phase_cycles[4] += t1 - t0;
      }
    }
  }
}

typedef struct SNB_GLUE(SNB_NAME, _worker_arg) {
  snb_state *st;
  int id;
} SNB_GLUE(SNB_NAME, _worker_arg);

static void *SNB_WORKER_FN(void *argp) {
  SNB_GLUE(SNB_NAME, _worker_arg) *warg =
    (SNB_GLUE(SNB_NAME, _worker_arg) *)argp;
  snb_state *st = warg->st;
  double *ub_scratch = st->ub_w[warg->id];
  double *gemm_scratch = st->gemm_w[warg->id];
  const uint32_t gen = st->generation;
  const int64_t count = st->task_count;
  for (;;) {
    const int64_t idx =
      atomic_fetch_add_explicit(&st->task_cursor, 1, memory_order_relaxed);
    if (idx >= count) break;
    snb_block *blk = &st->blocks[st->task_block[idx]];
    const int64_t s = st->task_snode[idx];
    const snb_snode *sn = &blk->snodes[s];
    const snb_edge *edges = blk->edges + sn->edges_off;
    for (int32_t e = 0; e < sn->edge_count; ++e) {
      while (atomic_load_explicit(&blk->done[edges[e].producer],
                                  memory_order_acquire) != gen) {
        __builtin_ia32_pause();
      }
    }
    SNB_TARGET_FN(st, blk, s, ub_scratch, gemm_scratch);
    atomic_store_explicit(&blk->done[s], gen, memory_order_release);
  }
  return NULL;
}

static int SNB_NAME(snb_state *st) {
  const matrix *a = st->a;
  st->status = 0;

  if (st->scale > 0) {
    double *rs = st->rs;
    const int64_t n = a->n;
    for (int64_t i = 0; i < n; ++i) rs[i] = 0.0;
    for (int64_t c = 0; c < n; ++c) {
      for (int64_t p = a->col_ptr[c]; p < a->col_ptr[c + 1]; ++p) {
        const double v = fabs(a->values[p]);
        if (st->scale == 1) {
          rs[a->row_idx[p]] += v;
        } else if (v > rs[a->row_idx[p]]) {
          rs[a->row_idx[p]] = v;
        }
      }
    }
    for (int64_t i = 0; i < n; ++i) {
      if (rs[i] == 0.0) rs[i] = 1.0;
    }
  }
  const double *rs = st->rs;
  const double *az = a->values;
  const scatter_program *sp = &st->sp;

  /* singleton blocks (trivial copies) run serially up front */
  for (UF_long block = 0; block < st->nblocks; ++block) {
    snb_block *blk = &st->blocks[block];
    const UF_long k1 = blk->k1;
    if (blk->nk != 1) continue;
    double sval = 0.0;
    for (int64_t p = sp->col_off_ptr[k1]; p < sp->col_off_ptr[k1 + 1];
         ++p) {
      st->offx[sp->off[p].dst] =
        st->scale > 0 ? az[sp->off[p].src] / rs[sp->off[p].oldrow]
                      : az[sp->off[p].src];
    }
    for (int64_t p = sp->col_blk_ptr[k1]; p < sp->col_blk_ptr[k1 + 1];
         ++p) {
      sval = st->scale > 0 ? az[sp->blk[p].src] / rs[sp->blk[p].oldrow]
                           : az[sp->blk[p].src];
    }
    st->udiag[k1] = sval;
  }

  if (st->threads <= 1) {
    for (UF_long block = 0; block < st->nblocks; ++block) {
      snb_block *blk = &st->blocks[block];
      if (blk->nk < 2) continue;
      for (int64_t s = 0; s < blk->snode_count; ++s) {
        SNB_TARGET_FN(st, blk, s, st->ub_w[0], st->gemm_w[0]);
      }
    }
    return 1;
  }

  st->generation++;
  atomic_store_explicit(&st->task_cursor, 0, memory_order_relaxed);
  const int nthreads = st->threads;
  pthread_t tids[64];
  SNB_GLUE(SNB_NAME, _worker_arg) wargs[64];
  for (int t = 1; t < nthreads; ++t) {
    wargs[t].st = st;
    wargs[t].id = t;
    pthread_create(&tids[t], NULL, SNB_WORKER_FN, &wargs[t]);
  }
  wargs[0].st = st;
  wargs[0].id = 0;
  SNB_WORKER_FN(&wargs[0]);
  for (int t = 1; t < nthreads; ++t) {
    pthread_join(tids[t], NULL);
  }
  return 1;
}

#undef SNB_LPTR
#undef SNB_GLUE2
#undef SNB_GLUE
#undef SNB_TARGET_FN
#undef SNB_WORKER_FN
