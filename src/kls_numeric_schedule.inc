/* Supernode preparation, sorting, and refactor schedule construction.
 * Included by kls_numeric_engines.inc in the original definition order.
 */

static void kls_snb_free(kls_solver *solver) {
  kls_snb_state *st = solver->snb;
  if (st == NULL) {
    return;
  }
  if (st->blocks != NULL) {
    for (UF_long b = 0; b < st->nblocks; ++b) {
      kls_snb_block *blk = &st->blocks[b];
      free(blk->snodes);
      free(blk->col2snode);
      free(blk->below_rows);
      free(blk->edges);
      free(blk->subcols);
      free(blk->edge_dsts);
      free(blk->panels);

      free(blk->writeback);
      free(blk->wb_col_off);
      free((void *)blk->published);
    }
    free(st->blocks);
  }
  free(st->off);
  free(st->blk);
  free(st->col_off_ptr);
  free(st->col_blk_ptr);
  for (int t = 0; t < st->nworkers; ++t) {
    if (st->ub_w != NULL) free(st->ub_w[t]);
    if (st->gemm_w != NULL) free(st->gemm_w[t]);
  }
  free(st->ub_w);
  free(st->gemm_w);
  free(st->rs_tmp);
  free(st->slevel_ptr);
  free(st->slevel_block);
  free(st->slevel_snode);
  free(st->seg_ptr);
  free(st->seg_coop);
  free(st);
  solver->snb = NULL;
  solver->snb_decision = 0;
}

static int kls_snb_trace_enabled(void) {
  const char *value = getenv("KLS_TRACE_SNB");
  return value != NULL && value[0] == '1';
}

/* Build the supernodal structures from the sorted numeric.  Returns 1 on
 * success; on failure or decline the solver is marked so we never retry. */
static int kls_snb_prepare(kls_solver *solver) {
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  if (solver->snb != NULL) {
    return 1;
  }
  if (solver->snb_declined) {
    return 0;
  }
  if (symbolic == NULL || numeric == NULL || numeric->LUbx == NULL ||
      !solver->snode_prepared) {
    /* the sorted-numeric prep has not run yet; retry on a later call */
    return 0;
  }
  if ((double)(numeric->lnz + numeric->unz) >
      (double)KLS_SNB_MAX_ENTRIES) {
    solver->snb_declined = 1;
    return 0;
  }
  kls_snb_state *st = (kls_snb_state *)calloc(1u, sizeof(*st));
  if (st == NULL) {
    solver->snb_declined = 1;
    return 0;
  }
  st->wnarrow = 8;
  {
    const char *wn = getenv("KLS_SNB_WNARROW");
    if (wn != NULL && wn[0] != '\0') {
      st->wnarrow = atoi(wn);
    }
  }

  st->wmax = KLS_SNB_WMAX;
  {
    const char *wm = getenv("KLS_SNB_WMAX");
    if (wm != NULL && wm[0] != '\0') {
      const int parsed = atoi(wm);
      if (parsed >= 4 && parsed <= 128) {
        st->wmax = parsed;
      }
    }
  }
  st->nblocks = symbolic->nblocks;
  st->blocks =
    (kls_snb_block *)calloc((size_t)st->nblocks, sizeof(*st->blocks));
  const UF_long maxblock = symbolic->maxblock;
  UF_long *union_rows =
    (UF_long *)malloc((size_t)maxblock * sizeof(*union_rows));
  UF_long *merge_tmp =
    (UF_long *)malloc((size_t)maxblock * sizeof(*merge_tmp));
  int64_t *starts = (int64_t *)malloc((size_t)maxblock * sizeof(*starts));
  int32_t *map = (int32_t *)malloc((size_t)maxblock * sizeof(*map));
  unsigned char *inset = (unsigned char *)calloc((size_t)maxblock, 1u);
  int64_t max_below = 0;
  int ok = st->blocks != NULL && union_rows != NULL && merge_tmp != NULL &&
           starts != NULL && map != NULL && inset != NULL;
  if (ok) {
    memset(map, 0xff, (size_t)maxblock * sizeof(*map));
  }

  for (UF_long block = 0; ok && block < st->nblocks; ++block) {
    const UF_long k1 = symbolic->R[block];
    const UF_long k2 = symbolic->R[block + 1u];
    const UF_long nk = k2 - k1;
    kls_snb_block *blk = &st->blocks[block];
    blk->k1 = k1;
    blk->nk = nk;
    if (nk < 2u || numeric->LUbx[block] == NULL) {
      continue;
    }
    const UF_long *lip = numeric->Lip + k1;
    const UF_long *llen = numeric->Llen + k1;
    const UF_long *uip = numeric->Uip + k1;
    const UF_long *ulen = numeric->Ulen + k1;
    double *lu = (double *)numeric->LUbx[block];

    /* greedy zeta amalgamation over sorted L patterns */
    double snb_zeta = (KLS_SNB_ZETA);
    {
      const char *z = getenv("KLS_SNB_ZETA");
      if (z != NULL && z[0] != '\0') {
        const double parsed = atof(z);
        if (parsed >= 0.0 && parsed <= 16.0) snb_zeta = parsed;
      }
    }
    blk->col2snode = (int32_t *)malloc((size_t)nk * sizeof(int32_t));
    if (blk->col2snode == NULL) { ok = 0; break; }
    int64_t ns = 0;
    {
      UF_long c = 0;
      while (c < nk) {
        UF_long *li; double *lx; UF_long len;
        kls_klu_get_pointer(lu, lip, llen, c, &li, &lx, &len);
        UF_long ulen_cur = len;
        memcpy(union_rows, li, (size_t)len * sizeof(UF_long));
        int64_t true_entries = len;
        UF_long w = 1;
        while (c + w < nk && w < (UF_long)st->wmax) {
          UF_long *li2; double *lx2; UF_long len2;
          kls_klu_get_pointer(lu, lip, llen, c + w, &li2, &lx2, &len2);
          UF_long i = 0, j = 0, m = 0;
          while (i < ulen_cur && j < len2) {
            const UF_long ri = union_rows[i];
            const UF_long rj = li2[j];
            const UF_long r = ri < rj ? ri : rj;
            if (ri == r) i++;
            if (rj == r) j++;
            merge_tmp[m++] = r;
          }
          while (i < ulen_cur) merge_tmp[m++] = union_rows[i++];
          while (j < len2) merge_tmp[m++] = li2[j++];
          UF_long below = 0;
          for (UF_long t = 0; t < m; ++t) {
            if (merge_tmp[t] > c + w) below++;
          }
          const UF_long nw = w + 1u;
          const double dense = (double)nw * (double)below +
                               (double)nw * (double)(nw - 1u) / 2.0;
          const double truth = (double)(true_entries + (int64_t)len2);
          if (dense > (1.0 + snb_zeta) * truth) {
            break;
          }
          memcpy(union_rows, merge_tmp, (size_t)m * sizeof(UF_long));
          ulen_cur = m;
          true_entries += len2;
          w++;
        }
        starts[ns] = (int64_t)c;
        for (UF_long i = 0; i < w; ++i) {
          blk->col2snode[c + i] = (int32_t)ns;
        }
        ns++;
        c += w;
      }
    }
    blk->snode_count = ns;
    blk->snodes =
      (kls_snb_snode *)calloc((size_t)ns, sizeof(*blk->snodes));
    if (blk->snodes == NULL) { ok = 0; break; }
    for (int64_t s = 0; s < ns; ++s) {
      blk->snodes[s].start = starts[s];
      const int64_t end = s + 1 < ns ? starts[s + 1] : (int64_t)nk;
      blk->snodes[s].width = (int32_t)(end - starts[s]);
    }

    /* below-row unions */
    int64_t rows_total = 0;
    for (int pass = 0; pass < 2 && ok; ++pass) {
      int64_t cursor = 0;
      for (int64_t s = 0; s < ns; ++s) {
        kls_snb_snode *sn = &blk->snodes[s];
        const int64_t end = sn->start + sn->width;
        UF_long ulen_cur = 0;
        for (int64_t c = sn->start; c < end; ++c) {
          UF_long *li; double *lx; UF_long len;
          kls_klu_get_pointer(lu, lip, llen, c, &li, &lx, &len);
          UF_long i = 0, j = 0, m = 0;
          while (i < ulen_cur && j < len) {
            const UF_long ri = union_rows[i];
            const UF_long rj = li[j];
            const UF_long r = ri < rj ? ri : rj;
            if (ri == r) i++;
            if (rj == r) j++;
            merge_tmp[m++] = r;
          }
          while (i < ulen_cur) merge_tmp[m++] = union_rows[i++];
          while (j < len) merge_tmp[m++] = li[j++];
          memcpy(union_rows, merge_tmp, (size_t)m * sizeof(UF_long));
          ulen_cur = m;
        }
        int64_t below = 0;
        for (UF_long t = 0; t < ulen_cur; ++t) {
          if (union_rows[t] >= (UF_long)end) {
            if (pass == 1) {
              blk->below_rows[cursor + below] = (int32_t)union_rows[t];
            }
            below++;
          }
        }
        if (pass == 0) {
          sn->below = (int32_t)below;
          sn->rows_off = rows_total;
          rows_total += below;
          if (below > max_below) max_below = below;
        }
        cursor += below;
      }
      if (pass == 0) {
        blk->below_rows = (int32_t *)malloc(
          (size_t)(rows_total > 0 ? rows_total : 1) * sizeof(int32_t));
        if (blk->below_rows == NULL) { ok = 0; }
      }
    }
    if (!ok) break;

    /* exact U-row unions -> edges with producer-column sublists */
    int64_t edges_total = 0;
    int64_t subs_total = 0;
    for (int pass = 0; pass < 2 && ok; ++pass) {
      int64_t sub_cursor = 0;
      int64_t edge_cursor = 0;
      for (int64_t s = 0; s < ns; ++s) {
        kls_snb_snode *sn = &blk->snodes[s];
        const int64_t end = sn->start + sn->width;
        int64_t m = 0;
        for (int64_t c = sn->start; c < end; ++c) {
          UF_long *ui; double *ux; UF_long len;
          kls_klu_get_pointer(lu, uip, ulen, c, &ui, &ux, &len);
          for (UF_long p = 0; p < len; ++p) {
            const UF_long j = ui[p];
            if (j >= (UF_long)sn->start) continue;
            if (!inset[j]) {
              inset[j] = 1;
              union_rows[m++] = j;
            }
          }
        }
        for (int64_t t = 0; t < m; ++t) inset[union_rows[t]] = 0;
        if (m > 1) {
          int sorted = 1;
          for (int64_t t = 1; t < m; ++t) {
            if (union_rows[t - 1] > union_rows[t]) { sorted = 0; break; }
          }
          if (!sorted) {
            for (int64_t i2 = 1; i2 < m; ++i2) {
              const UF_long v = union_rows[i2];
              int64_t j2 = i2 - 1;
              while (j2 >= 0 && union_rows[j2] > v) {
                union_rows[j2 + 1] = union_rows[j2];
                j2--;
              }
              union_rows[j2 + 1] = v;
            }
          }
        }
        int64_t cnt = 0;
        int32_t last_prod = -1;
        for (int64_t t = 0; t < m; ++t) {
          const int32_t prod = blk->col2snode[union_rows[t]];
          if (prod != last_prod) {
            if (pass == 1) {
              kls_snb_edge *e = &blk->edges[edge_cursor + cnt];
              e->producer = prod;
              e->u_start = (int32_t)t;
              e->sub_count = 0;
              e->sub_off = sub_cursor + t;
            }
            cnt++;
            last_prod = prod;
          }
          if (pass == 1) {
            blk->edges[edge_cursor + cnt - 1].sub_count++;
            blk->subcols[sub_cursor + t] =
              (int32_t)(union_rows[t] - blk->snodes[prod].start);
          }
        }
        if (pass == 0) {
          sn->edge_count = (int32_t)cnt;
          sn->edges_off = edges_total;
          edges_total += cnt;
          subs_total += m;
          sn->uw = (int32_t)m;
        }
        sub_cursor += m;
        edge_cursor += cnt;
      }
      if (pass == 0) {
        blk->edges = (kls_snb_edge *)malloc(
          (size_t)(edges_total > 0 ? edges_total : 1) * sizeof(kls_snb_edge));
        blk->subcols = (int32_t *)malloc(
          (size_t)(subs_total > 0 ? subs_total : 1) * sizeof(int32_t));
        if (blk->edges == NULL || blk->subcols == NULL) { ok = 0; }
      }
    }
    if (!ok) break;

    /* panel layout */
    int64_t panel_doubles = 0;
    for (int64_t s = 0; s < ns; ++s) {
      kls_snb_snode *sn = &blk->snodes[s];
      sn->height = sn->uw + sn->width + sn->below;
      sn->panel = panel_doubles;
      panel_doubles += (int64_t)sn->height * sn->width;
    }
    if (panel_doubles > KLS_SNB_MAX_ENTRIES) { ok = 0; break; }
    blk->panel_doubles = panel_doubles;
    /* panels materialize after the schedule/cost gate below */

    /* per-edge scatter destinations via a transient row map */
    int64_t dsts_total = 0;
    for (int64_t s = 0; s < ns; ++s) {
      kls_snb_snode *sn = &blk->snodes[s];
      kls_snb_edge *edges = blk->edges + sn->edges_off;
      for (int32_t e = 0; e < sn->edge_count; ++e) {
        edges[e].dst_off = dsts_total;
        dsts_total += blk->snodes[edges[e].producer].below;
      }
    }
    blk->edge_dsts = (int32_t *)malloc(
      (size_t)(dsts_total > 0 ? dsts_total : 1) * sizeof(int32_t));
    if (blk->edge_dsts == NULL) { ok = 0; break; }
    for (int64_t s = 0; s < ns; ++s) {
      kls_snb_snode *sn = &blk->snodes[s];
      const kls_snb_edge *edges = blk->edges + sn->edges_off;
      for (int32_t e = 0; e < sn->edge_count; ++e) {
        const kls_snb_snode *ps = &blk->snodes[edges[e].producer];
        const int32_t *sub = blk->subcols + edges[e].sub_off;
        for (int32_t i = 0; i < edges[e].sub_count; ++i) {
          map[ps->start + sub[i]] = edges[e].u_start + i;
        }
      }
      for (int32_t i = 0; i < sn->width; ++i) {
        map[sn->start + i] = sn->uw + i;
      }
      {
        const int32_t *rows = blk->below_rows + sn->rows_off;
        for (int32_t i = 0; i < sn->below; ++i) {
          map[rows[i]] = sn->uw + sn->width + i;
        }
      }
      for (int32_t e = 0; e < sn->edge_count; ++e) {
        const kls_snb_snode *ps = &blk->snodes[edges[e].producer];
        const int32_t *prows = blk->below_rows + ps->rows_off;
        int32_t *dsts = blk->edge_dsts + edges[e].dst_off;
        for (int32_t r = 0; r < ps->below; ++r) {
          dsts[r] = map[prows[r]];
        }
      }
      /* clear */
      for (int32_t e = 0; e < sn->edge_count; ++e) {
        const kls_snb_snode *ps = &blk->snodes[edges[e].producer];
        const int32_t *sub = blk->subcols + edges[e].sub_off;
        for (int32_t i = 0; i < edges[e].sub_count; ++i) {
          map[ps->start + sub[i]] = -1;
        }
      }
      for (int32_t i = 0; i < sn->width; ++i) map[sn->start + i] = -1;
      {
        const int32_t *rows = blk->below_rows + sn->rows_off;
        for (int32_t i = 0; i < sn->below; ++i) map[rows[i]] = -1;
      }
    }

    /* writeback lists: for each column, panel index per packed L then U
     * entry (block walk order matches the refactor's writeback loop) */
    {
      int64_t wb_total = 0;
      for (UF_long k = 0; k < nk; ++k) {
        wb_total += llen[k] + ulen[k];
      }
      blk->writeback = (int32_t *)malloc(
        (size_t)(wb_total > 0 ? wb_total : 1) * sizeof(int32_t));
      blk->wb_col_off =
        (int64_t *)malloc(((size_t)nk + 1u) * sizeof(int64_t));
      if (blk->writeback == NULL || blk->wb_col_off == NULL) {
        ok = 0;
        break;
      }
      int64_t wb = 0;
      for (UF_long k = 0; ok && k < nk; ++k) {
        blk->wb_col_off[k] = wb;
        const int32_t s = blk->col2snode[k];
        const kls_snb_snode *sn = &blk->snodes[s];
        const int32_t c = (int32_t)(k - sn->start);
        UF_long *ri; double *rv; UF_long rl;
        kls_klu_get_pointer(lu, lip, llen, k, &ri, &rv, &rl);
        for (UF_long p = 0; p < rl; ++p) {
          const UF_long row = ri[p];
          int32_t lr = -1;
          if (row < (UF_long)(sn->start + sn->width)) {
            lr = sn->uw + (int32_t)(row - sn->start);
          } else {
            const int32_t *rows = blk->below_rows + sn->rows_off;
            int32_t lo = 0, hi = sn->below - 1;
            while (lo <= hi) {
              const int32_t mid = (lo + hi) >> 1;
              if (rows[mid] < (int32_t)row) lo = mid + 1;
              else if (rows[mid] > (int32_t)row) hi = mid - 1;
              else { lr = sn->uw + sn->width + mid; break; }
            }
          }
          if (lr < 0) { ok = 0; break; }
          blk->writeback[wb++] = (int32_t)((int64_t)c * sn->height + lr);
        }
        kls_klu_get_pointer(lu, uip, ulen, k, &ri, &rv, &rl);
        for (UF_long p = 0; ok && p < rl; ++p) {
          const UF_long row = ri[p];
          int32_t lr = -1;
          if (row >= (UF_long)sn->start) {
            lr = sn->uw + (int32_t)(row - sn->start);
          } else {
            const int32_t prod = blk->col2snode[row];
            const kls_snb_edge *edges = blk->edges + sn->edges_off;
            for (int32_t e2 = 0; e2 < sn->edge_count; ++e2) {
              if (edges[e2].producer == prod) {
                const int32_t *sub = blk->subcols + edges[e2].sub_off;
                const int32_t want = (int32_t)(row - blk->snodes[prod].start);
                int32_t lo = 0, hi = edges[e2].sub_count - 1;
                while (lo <= hi) {
                  const int32_t mid = (lo + hi) >> 1;
                  if (sub[mid] < want) lo = mid + 1;
                  else if (sub[mid] > want) hi = mid - 1;
                  else { lr = edges[e2].u_start + mid; break; }
                }
                break;
              }
            }
          }
          if (lr < 0) { ok = 0; break; }
          blk->writeback[wb++] = (int32_t)((int64_t)c * sn->height + lr);
        }
      }
      if (ok) {
        blk->wb_col_off[nk] = wb;
      }
    }
  }

  /* snode-level schedule: level(s) = 1 + max(level of producers); blocks
   * are independent so their levels merge directly */
  if (ok) {
    int64_t total_snodes = 0;
    int64_t max_level = 0;
    for (UF_long b = 0; b < st->nblocks; ++b) {
      if (st->blocks[b].nk >= 2u) total_snodes += st->blocks[b].snode_count;
    }
    int32_t *levels = NULL;
    int64_t *counts = NULL;
    if (total_snodes > 0) {
      levels = (int32_t *)malloc((size_t)total_snodes * sizeof(int32_t));
    }
    if (levels != NULL) {
      int64_t base = 0;
      for (UF_long b = 0; b < st->nblocks; ++b) {
        const kls_snb_block *blk = &st->blocks[b];
        if (blk->nk < 2u) continue;
        for (int64_t s = 0; s < blk->snode_count; ++s) {
          const kls_snb_snode *sn = &blk->snodes[s];
          const kls_snb_edge *edges = blk->edges + sn->edges_off;
          int32_t lv = 0;
          for (int32_t e = 0; e < sn->edge_count; ++e) {
            const int32_t pl = levels[base + edges[e].producer];
            if (pl + 1 > lv) lv = pl + 1;
          }
          levels[base + s] = lv;
          if (lv > max_level) max_level = lv;
        }
        base += blk->snode_count;
      }
      st->slevel_count = max_level + 1;
      counts = (int64_t *)calloc((size_t)st->slevel_count + 1u,
                                 sizeof(int64_t));
      st->slevel_ptr = (int64_t *)calloc((size_t)st->slevel_count + 1u,
                                         sizeof(int64_t));
      st->slevel_block =
        (int32_t *)malloc((size_t)total_snodes * sizeof(int32_t));
      st->slevel_snode =
        (int32_t *)malloc((size_t)total_snodes * sizeof(int32_t));
      if (counts != NULL && st->slevel_ptr != NULL &&
          st->slevel_block != NULL && st->slevel_snode != NULL) {
        base = 0;
        for (UF_long b = 0; b < st->nblocks; ++b) {
          const kls_snb_block *blk = &st->blocks[b];
          if (blk->nk < 2u) continue;
          for (int64_t s = 0; s < blk->snode_count; ++s) {
            counts[levels[base + s] + 1]++;
          }
          base += blk->snode_count;
        }
        for (int64_t l = 0; l < st->slevel_count; ++l) {
          st->slevel_ptr[l + 1] = st->slevel_ptr[l] + counts[l + 1];
        }
        memcpy(counts, st->slevel_ptr,
               ((size_t)st->slevel_count + 1u) * sizeof(int64_t));
        base = 0;
        for (UF_long b = 0; b < st->nblocks; ++b) {
          const kls_snb_block *blk = &st->blocks[b];
          if (blk->nk < 2u) continue;
          for (int64_t s = 0; s < blk->snode_count; ++s) {
            const int64_t at = counts[levels[base + s]]++;
            st->slevel_block[at] = (int32_t)b;
            st->slevel_snode[at] = (int32_t)s;
          }
          base += blk->snode_count;
        }
      } else {
        free(st->slevel_ptr);
        free(st->slevel_block);
        free(st->slevel_snode);
        st->slevel_ptr = NULL;
        st->slevel_block = NULL;
        st->slevel_snode = NULL;
        st->slevel_count = 0;
      }
      free(counts);
      free(levels);
    }
    if (st->slevel_ptr != NULL && st->slevel_count > 0) {
      double coop = 0.0, wide = 0.0, total = 0.0;
      for (int64_t l = 0; l < st->slevel_count; ++l) {
        const int64_t begin = st->slevel_ptr[l];
        const int64_t end = st->slevel_ptr[l + 1];
        double level_work = 0.0;
        double level_max = 0.0;
        for (int64_t t = begin; t < end; ++t) {
          const kls_snb_snode *sn =
            &st->blocks[st->slevel_block[t]].snodes[st->slevel_snode[t]];
          const double wk = (double)sn->height * (double)sn->width;
          level_work += wk;
          if (wk > level_max) level_max = wk;
        }
        total += level_work;
        if (end - begin == 1 &&
            level_max >= (double)KLS_SNB_COOP_MIN_DOUBLES) {
          coop += level_work;
        } else if (end - begin >= 4) {
          wide += level_work;
        }
      }
      st->par_work_share = total > 0.0 ? (coop + wide) / total : 0.0;
      st->coop_work_share = total > 0.0 ? coop / total : 0.0;
    }
    /* pipelined segments: split the topological task list at cooperative
       levels; inside a segment cross-level dependencies are enforced by
       the per-edge published-column spins, so no barriers are needed */
    if (st->slevel_ptr != NULL && st->slevel_count > 0) {
      st->seg_ptr = (int64_t *)malloc(((size_t)st->slevel_count * 2u + 2u) *
                                      sizeof(int64_t));
      st->seg_coop = (unsigned char *)malloc((size_t)st->slevel_count + 1u);
      if (st->seg_ptr != NULL && st->seg_coop != NULL) {
        int64_t nseg = 0;
        int64_t open_start = -1;
        for (int64_t l = 0; l < st->slevel_count; ++l) {
          const int64_t begin = st->slevel_ptr[l];
          const int64_t end = st->slevel_ptr[l + 1];
          int is_coop = 0;
          if (end - begin == 1) {
            const kls_snb_snode *sn =
              &st->blocks[st->slevel_block[begin]]
                 .snodes[st->slevel_snode[begin]];
            is_coop = (int64_t)sn->height * sn->width >=
                      KLS_SNB_COOP_MIN_DOUBLES;
          }
          if (is_coop) {
            if (open_start >= 0) {
              st->seg_ptr[2 * nseg] = open_start;
              st->seg_ptr[2 * nseg + 1] = begin;
              st->seg_coop[nseg] = 0;
              nseg++;
              open_start = -1;
            }
            st->seg_ptr[2 * nseg] = begin;
            st->seg_ptr[2 * nseg + 1] = end;
            st->seg_coop[nseg] = 1;
            nseg++;
          } else if (open_start < 0) {
            open_start = begin;
          }
        }
        if (open_start >= 0) {
          st->seg_ptr[2 * nseg] = open_start;
          st->seg_ptr[2 * nseg + 1] =
            st->slevel_ptr[st->slevel_count];
          st->seg_coop[nseg] = 0;
          nseg++;
        }
        st->seg_count = nseg;
      } else {
        free(st->seg_ptr);
        free(st->seg_coop);
        st->seg_ptr = NULL;
        st->seg_coop = NULL;
        st->seg_count = 0;
      }
    }
    for (UF_long b = 0; ok && b < st->nblocks; ++b) {
      kls_snb_block *blk = &st->blocks[b];
      if (blk->nk < 2u || blk->snode_count <= 0) continue;
      blk->published = (_Atomic int32_t *)calloc(
        (size_t)blk->snode_count, sizeof(_Atomic int32_t));
      if (blk->published == NULL) ok = 0;
    }
  }

  /* cost gate: at t4 the trial+materialization cost must stay small
     relative to the factor already paid, and the schedule must have
     parallelizable shape; serial trials are always worth it (the mapped
     incumbents they race are far slower than any prep). */
  if (ok && solver->snb_trial_budget > 0.0) {
    /* t4: with the 0.6 decisive-adoption margin, the only measured
       adopter is pre2 (coop share 0.58); chain giants with tiny snodes
       (memchip 1.5-wide) and scattered shapes always reject, and their
       trial+panels would tax the once-charged factor term ~4s each.
       Gate on the cooperative work share computed from the schedule. */
    if (kls_snb_trace_enabled()) {
      fprintf(stderr, "KLS snb: shares coop %.2f par %.2f\n",
              st->coop_work_share, st->par_work_share);
    }
    if ((st->coop_work_share < 0.50) &&
        getenv("KLS_SNB_FORCE_TRIAL") == NULL) {
      if (kls_snb_trace_enabled()) {
        fprintf(stderr, "KLS snb: t4 gated off (coop %.2f)\n",
                st->coop_work_share);
      }
      ok = 0;
    }
  }
  if (ok) {
    for (UF_long b = 0; b < st->nblocks; ++b) {
      kls_snb_block *blk = &st->blocks[b];
      if (blk->nk < 2u || blk->panel_doubles <= 0) continue;
      blk->panels =
        (double *)malloc((size_t)blk->panel_doubles * sizeof(double));
      if (blk->panels == NULL) { ok = 0; break; }
      memset(blk->panels, 0, (size_t)blk->panel_doubles * sizeof(double));

    }
  }

  /* A-scatter program (klu_l_refactor walk order) */
  if (ok) {
    const UF_long n = symbolic->n;
    const UF_long nzoff = symbolic->nzoff;
    const UF_long *pinv = numeric->Pinv;
    int64_t off_count = 0;
    int64_t blk_count = 0;
    {
      int64_t poff = 0;
      for (UF_long block = 0; block < st->nblocks; ++block) {
        const UF_long k1 = symbolic->R[block];
        const UF_long k2 = symbolic->R[block + 1u];
        for (UF_long k = k1; k < k2; ++k) {
          const UF_long oldcol = symbolic->Q[k];
          for (UF_long p = solver->col_ptr[oldcol];
               p < solver->col_ptr[oldcol + 1u]; ++p) {
            const UF_long newrow = pinv[solver->row_idx[p]] - k1;
            if (newrow < 0 && poff < (int64_t)nzoff) {
              off_count++;
              poff++;
            } else {
              blk_count++;
            }
          }
        }
      }
    }
    st->off = (kls_snb_entry *)malloc(
      (size_t)(off_count > 0 ? off_count : 1) * sizeof(kls_snb_entry));
    st->blk = (kls_snb_entry *)malloc(
      (size_t)(blk_count > 0 ? blk_count : 1) * sizeof(kls_snb_entry));
    st->col_off_ptr = (int64_t *)calloc((size_t)n + 1u, sizeof(int64_t));
    st->col_blk_ptr = (int64_t *)calloc((size_t)n + 1u, sizeof(int64_t));
    ok = st->off != NULL && st->blk != NULL && st->col_off_ptr != NULL &&
         st->col_blk_ptr != NULL;
    if (ok) {
      int64_t poff = 0;
      int64_t pblk = 0;
      for (UF_long block = 0; ok && block < st->nblocks; ++block) {
        const UF_long k1 = symbolic->R[block];
        const UF_long k2 = symbolic->R[block + 1u];
        const kls_snb_block *blk = &st->blocks[block];
        for (UF_long k = k1; ok && k < k2; ++k) {
          st->col_off_ptr[k] = poff;
          st->col_blk_ptr[k] = pblk;
          const UF_long oldcol = symbolic->Q[k];
          for (UF_long p = solver->col_ptr[oldcol];
               p < solver->col_ptr[oldcol + 1u]; ++p) {
            const UF_long oldrow = solver->row_idx[p];
            const UF_long newrow = pinv[oldrow] - k1;
            if (newrow < 0 && poff < (int64_t)nzoff) {
              st->off[poff].dst = poff;
              st->off[poff].src = (int64_t)p;
              st->off[poff].oldrow = (int64_t)oldrow;
              poff++;
              continue;
            }
            int64_t dst = 0;
            if (blk->nk >= 2u) {
              const int32_t s = blk->col2snode[k - k1];
              const kls_snb_snode *sn = &blk->snodes[s];
              const int32_t c = (int32_t)(k - k1 - sn->start);
              int32_t lr = -1;
              if (newrow >= sn->start &&
                  newrow < sn->start + sn->width) {
                lr = sn->uw + (int32_t)(newrow - sn->start);
              } else if (newrow >= sn->start + sn->width) {
                const int32_t *rows = blk->below_rows + sn->rows_off;
                int32_t lo = 0, hi = sn->below - 1;
                while (lo <= hi) {
                  const int32_t mid = (lo + hi) >> 1;
                  if (rows[mid] < (int32_t)newrow) lo = mid + 1;
                  else if (rows[mid] > (int32_t)newrow) hi = mid - 1;
                  else { lr = sn->uw + sn->width + mid; break; }
                }
              } else {
                const int32_t prod = blk->col2snode[newrow];
                const kls_snb_edge *edges = blk->edges + sn->edges_off;
                for (int32_t e2 = 0; e2 < sn->edge_count; ++e2) {
                  if (edges[e2].producer == prod) {
                    const int32_t *sub = blk->subcols + edges[e2].sub_off;
                    const int32_t want =
                      (int32_t)(newrow - blk->snodes[prod].start);
                    int32_t lo = 0, hi = edges[e2].sub_count - 1;
                    while (lo <= hi) {
                      const int32_t mid = (lo + hi) >> 1;
                      if (sub[mid] < want) lo = mid + 1;
                      else if (sub[mid] > want) hi = mid - 1;
                      else { lr = edges[e2].u_start + mid; break; }
                    }
                    break;
                  }
                }
              }
              if (lr < 0) { ok = 0; break; }
              dst = sn->panel + (int64_t)c * sn->height + lr;
            }
            st->blk[pblk].dst = dst;
            st->blk[pblk].src = (int64_t)p;
            st->blk[pblk].oldrow = (int64_t)oldrow;
            pblk++;
          }
        }
      }
      st->col_off_ptr[n] = poff;
      st->col_blk_ptr[n] = pblk;
    }
  }

  if (ok) {
    const int64_t gmax_below = max_below > 0 ? max_below : 1;
    int nworkers = solver->options.threads;
    if (nworkers < 1) nworkers = 1;
    if (nworkers > 64) nworkers = 64;
    st->nworkers = nworkers;
    st->ub_w = (double **)calloc((size_t)nworkers, sizeof(double *));
    st->gemm_w = (double **)calloc((size_t)nworkers, sizeof(double *));
    st->rs_tmp = (double *)malloc((size_t)symbolic->n * sizeof(double));
    ok = st->ub_w != NULL && st->gemm_w != NULL && st->rs_tmp != NULL;
    for (int t = 0; ok && t < nworkers; ++t) {
      st->ub_w[t] = (double *)malloc(2u * (size_t)st->wmax *
                                     (size_t)st->wmax * sizeof(double));
      st->gemm_w[t] = (double *)malloc((size_t)gmax_below *
                                       (size_t)st->wmax *
                                       sizeof(double));
      ok = st->ub_w[t] != NULL && st->gemm_w[t] != NULL;
    }
    if (ok) {
      st->ub = st->ub_w[0];
      st->gemm_temp = st->gemm_w[0];
    }
  }

  free(union_rows);
  free(merge_tmp);
  free(starts);
  free(map);
  free(inset);
  if (!ok) {
    solver->snb = st;
    kls_snb_free(solver);
    solver->snb_declined = 1;
    return 0;
  }
  solver->snb = st;
  if (kls_snb_trace_enabled()) {
    int64_t snodes = 0, cols = 0, panels = 0;
    for (UF_long b = 0; b < st->nblocks; ++b) {
      if (st->blocks[b].nk < 2u) continue;
      snodes += st->blocks[b].snode_count;
      cols += st->blocks[b].nk;
      panels += st->blocks[b].panel_doubles;
    }
    fprintf(stderr,
            "KLS snb: prepared %ld snodes over %ld cols, %.1f MB panels, "
            "%ld snode levels, par share %.2f\n",
            (long)snodes, (long)cols, (double)panels * 8.0 / 1048576.0,
            (long)st->slevel_count, st->par_work_share);
  }
  return 1;
}

/* Record the first (lowest-index) singular pivot seen this refactor. */
static void kls_snb_record_singular(kls_snb_state *st, kls_solver *solver,
                                    UF_long gk) {
  int64_t seen = atomic_load_explicit(&st->singular_rank,
                                      memory_order_relaxed);
  while ((seen < 0 || (int64_t)gk < seen) &&
         !atomic_compare_exchange_weak_explicit(&st->singular_rank, &seen,
                                                (int64_t)gk,
                                                memory_order_relaxed,
                                                memory_order_relaxed)) {
  }
  if (solver->common.halt_if_singular) {
    atomic_store_explicit(&st->failed, 1, memory_order_release);
  }
  if (solver->zero_pivot_collect != NULL) {
    const long slot = atomic_fetch_add_explicit(
      &solver->zero_pivot_collect_count, 1, memory_order_relaxed);
    if (slot < solver->zero_pivot_collect_cap) {
      solver->zero_pivot_collect[slot] = gk;
    }
  }
}

/* Process one supernode: panel compute (memset, A-scatter program, edge
 * updates, fused no-pivot getrf) followed by an immediate writeback of
 * its member columns into the packed numeric, so consumers that read
 * LUbx (including a future pipelined tail) always see finished columns.
 * Thread-safe across distinct snodes of one level. */
static void kls_snb_process_target(kls_solver *solver, kls_snb_state *st,
                                   kls_snb_block *blk, int64_t s,
                                   const double *az, const double *rs,
                                   int scale, double *ub_scratch,
                                   double *gemm_scratch, int pipelined) {
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
    const kls_snb_snode *ps = &blk->snodes[edges[e].producer];
    const int32_t pw = ps->width;
    const int32_t ph = ps->below;
    const int32_t sc = edges[e].sub_count;
    const int32_t *sub = blk->subcols + edges[e].sub_off;
    const int32_t *dsts = blk->edge_dsts + edges[e].dst_off;
    const int32_t ustart = edges[e].u_start;
    const double *Lp = panels + ps->panel;
    const int32_t pld = ps->height;
    const double *Ldiag = Lp + ps->uw;
    const double *Lbelow = Lp + ps->uw + pw;

    if (pipelined) {
      /* the producer publishes each getrf column right after its
         division; this edge only reads columns <= sub[sc-1] */
      const int32_t need = sub[sc - 1] + 1;
      while (atomic_load_explicit(&blk->published[edges[e].producer],
                                  memory_order_acquire) < need) {
        kls_cpu_relax();
      }
    }
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
      continue;
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
    if (ph == 0) continue;
    double *Brow = ub_scratch + (size_t)sc * w;
    for (int32_t c = 0; c < w; ++c) {
      const double *Bc = B + (size_t)c * sc;
      for (int32_t i = 0; i < sc; ++i) {
        Brow[(size_t)i * w + c] = Bc[i];
      }
    }
    double *T = gemm_scratch;
    {
      const double *lb = Lbelow + (size_t)sub[0] * pld;
      const double *br = Brow;
      for (int32_t r = 0; r < ph; ++r) {
        const double l = lb[r];
        double *tr = T + (size_t)r * w;
        for (int32_t c = 0; c < w; ++c) tr[c] = l * br[c];
      }
      for (int32_t i = 1; i < sc; ++i) {
        const double *lbi = Lbelow + (size_t)sub[i] * pld;
        const double *bri = Brow + (size_t)i * w;
        for (int32_t r = 0; r < ph; ++r) {
          const double l = lbi[r];
          double *tr = T + (size_t)r * w;
          for (int32_t c = 0; c < w; ++c) tr[c] += l * bri[c];
        }
      }
    }
    for (int32_t r = 0; r < ph; ++r) {
      const int32_t dst = dsts[r];
      const double *tr = T + (size_t)r * w;
      if (dst >= 0) {
        double *Wd = W + dst;
        for (int32_t c = 0; c < w; ++c) {
          Wd[(size_t)c * H] -= tr[c];
        }
      }
    }
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

/* Shared prelude: scale recompute + singleton blocks + state reset. */
static int kls_snb_refactor_prelude(kls_solver *solver,
                                    const double *numeric_values) {
  kls_snb_state *st = solver->snb;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  trilinos_klu_l_common *common = &solver->common;
  const int scale = (int)common->scale;
  common->status = TRILINOS_KLU_OK;
  common->numerical_rank = KLS_KLU_EMPTY;
  common->singular_col = KLS_KLU_EMPTY;
  atomic_store_explicit(&st->singular_rank, -1, memory_order_relaxed);
  atomic_store_explicit(&st->failed, 0, memory_order_relaxed);
  if (scale > 0) {
    if (numeric->Rs == NULL ||
        !trilinos_klu_l_scale((UF_long)scale, solver->n, solver->col_ptr,
                              solver->row_idx, (double *)numeric_values,
                              numeric->Rs, NULL, common)) {
      return 0;
    }
  }
  const double *rs = numeric->Rs;
  const double *az = numeric_values;
  double *offx = (double *)numeric->Offx;
  double *udiag = (double *)numeric->Udiag;
  for (UF_long block = 0; block < st->nblocks; ++block) {
    kls_snb_block *blk = &st->blocks[block];
    if (blk->nk != 1u) continue;
    const UF_long k1 = blk->k1;
    double sv = 0.0;
    for (int64_t p = st->col_off_ptr[k1]; p < st->col_off_ptr[k1 + 1u];
         ++p) {
      offx[st->off[p].dst] =
        scale > 0 ? az[st->off[p].src] / rs[st->off[p].oldrow]
                  : az[st->off[p].src];
    }
    for (int64_t p = st->col_blk_ptr[k1]; p < st->col_blk_ptr[k1 + 1u];
         ++p) {
      sv = scale > 0 ? az[st->blk[p].src] / rs[st->blk[p].oldrow]
                     : az[st->blk[p].src];
    }
    udiag[k1] = sv;
    if (sv == 0.0) {
      kls_snb_record_singular(st, solver, k1);
    }
  }
  return 1;
}

/* Shared epilogue: singular reporting + the pivot-order Rs permute. */
static UF_long kls_snb_refactor_epilogue(kls_solver *solver) {
  kls_snb_state *st = solver->snb;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  trilinos_klu_l_common *common = &solver->common;
  const int64_t rank = atomic_load_explicit(&st->singular_rank,
                                            memory_order_relaxed);
  if (rank >= 0) {
    common->status = TRILINOS_KLU_SINGULAR;
    common->numerical_rank = (UF_long)rank;
    common->singular_col = solver->symbolic->Q[rank];
    if (common->halt_if_singular) {
      return 0;
    }
  }
  if ((int)common->scale > 0) {
    double *rs_tmp = st->rs_tmp;
    for (UF_long k = 0; k < solver->n; ++k) {
      rs_tmp[k] = numeric->Rs[numeric->Pnum[k]];
    }
    memcpy(numeric->Rs, rs_tmp, (size_t)solver->n * sizeof(double));
  }

  return 1;
}

static UF_long kls_snb_refactor(kls_solver *solver,
                                const double *numeric_values) {
  kls_snb_state *st = solver->snb;
  const int scale = (int)solver->common.scale;
  if (!kls_snb_refactor_prelude(solver, numeric_values)) {
    return 0;
  }
  const double *rs = solver->numeric->Rs;
  for (UF_long block = 0; block < st->nblocks; ++block) {
    kls_snb_block *blk = &st->blocks[block];
    if (blk->nk < 2u) continue;
    for (int64_t s = 0; s < blk->snode_count; ++s) {
      kls_snb_process_target(solver, st, blk, s, numeric_values, rs, scale,
                             st->ub_w[0], st->gemm_w[0], 0);
      if (atomic_load_explicit(&st->failed, memory_order_acquire)) {
        return 0;
      }
    }
  }
  return kls_snb_refactor_epilogue(solver);
}

/* Reusable spin barrier for the cooperative (within-snode) phases: the
 * participating crew is fixed per refactor and phases are microseconds
 * apart, so futex round trips would dominate a pthread barrier. */
typedef struct kls_snb_spin_barrier {
  _Atomic int arrived;
  _Atomic unsigned generation;
} kls_snb_spin_barrier;

static void kls_snb_spin_barrier_wait(kls_snb_spin_barrier *b, int crew) {
  const unsigned gen =
    atomic_load_explicit(&b->generation, memory_order_acquire);
  if (atomic_fetch_add_explicit(&b->arrived, 1, memory_order_acq_rel) ==
      crew - 1) {
    atomic_store_explicit(&b->arrived, 0, memory_order_relaxed);
    atomic_fetch_add_explicit(&b->generation, 1u, memory_order_release);
    return;
  }
  while (atomic_load_explicit(&b->generation, memory_order_acquire) == gen) {
    kls_cpu_relax();
  }
}

typedef struct kls_snb_par_ctx {
  kls_solver *solver;
  kls_snb_state *st;
  const double *values;
  const double *rs;
  int scale;
  int nthreads;
  pthread_barrier_t barrier;
  kls_snb_spin_barrier coop;
  _Atomic int64_t cursor;
  _Atomic int go; /* 0 = hold, 1 = run, 2 = abort (spawn failed) */
} kls_snb_par_ctx;

/* Cooperative processing of one large supernode by the whole crew.
 * Edge updates: tid 0 gathers/TRSMs/writes back the (small) UB block,
 * then every worker computes and scatters a disjoint row slice of the
 * GEMM (producer below-rows split across the crew; the shared B row
 * mirror is read-only).  Narrow edges are cheap and stay on tid 0.
 * getrf: i-steps run row-sliced with a spin barrier per step so row i
 * is fully updated before the crew reads it.  Writeback: column-sliced.
 * All slices are disjoint; the only shared writes are the barrier. */
static void kls_snb_process_target_coop(kls_snb_par_ctx *ctx,
                                        kls_snb_block *blk, int64_t s,
                                        int tid) {
  kls_solver *solver = ctx->solver;
  kls_snb_state *st = ctx->st;
  const int crew = ctx->nthreads;
  const double *az = ctx->values;
  const double *rs = ctx->rs;
  const int scale = ctx->scale;
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
  /* memset + A-scatter: column-sliced */
  for (int32_t c = tid; c < w; c += crew) {
    memset(W + (size_t)c * H, 0, (size_t)H * sizeof(double));
  }
  kls_snb_spin_barrier_wait(&ctx->coop, crew);
  for (int32_t c = tid; c < w; c += crew) {
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
  kls_snb_spin_barrier_wait(&ctx->coop, crew);
  const kls_snb_edge *edges = blk->edges + sn->edges_off;
  for (int32_t e = 0; e < sn->edge_count; ++e) {
    const kls_snb_snode *ps = &blk->snodes[edges[e].producer];
    const int32_t pw = ps->width;
    const int32_t ph = ps->below;
    const int32_t sc = edges[e].sub_count;
    const int32_t *sub = blk->subcols + edges[e].sub_off;
    const int32_t *dsts = blk->edge_dsts + edges[e].dst_off;
    const int32_t ustart = edges[e].u_start;
    const double *Lp = panels + ps->panel;
    const int32_t pld = ps->height;
    const double *Ldiag = Lp + ps->uw;
    const double *Lbelow = Lp + ps->uw + pw;

    if (sc < wnarrow || (int64_t)ph * w < KLS_SNB_COOP_MIN_DOUBLES) {
      /* small edge: tid 0 runs the ordinary serial edge body */
      if (tid == 0) {
        double *B = st->ub_w[0];
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
        } else {
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
          double *Brow = st->ub_w[0] + (size_t)sc * w;
          for (int32_t c = 0; c < w; ++c) {
            const double *Bc = B + (size_t)c * sc;
            for (int32_t i = 0; i < sc; ++i) {
              Brow[(size_t)i * w + c] = Bc[i];
            }
          }
          double *T = st->gemm_w[0];
          {
            const double *lb = Lbelow + (size_t)sub[0] * pld;
            const double *br = Brow;
            for (int32_t r = 0; r < ph; ++r) {
              const double l = lb[r];
              double *tr = T + (size_t)r * w;
              for (int32_t c = 0; c < w; ++c) tr[c] = l * br[c];
            }
            for (int32_t i = 1; i < sc; ++i) {
              const double *lbi = Lbelow + (size_t)sub[i] * pld;
              const double *bri = Brow + (size_t)i * w;
              for (int32_t r = 0; r < ph; ++r) {
                const double l = lbi[r];
                double *tr = T + (size_t)r * w;
                for (int32_t c = 0; c < w; ++c) tr[c] += l * bri[c];
              }
            }
          }
          for (int32_t r = 0; r < ph; ++r) {
            const int32_t dst = dsts[r];
            const double *tr = T + (size_t)r * w;
            if (dst >= 0) {
              double *Wd = W + dst;
              for (int32_t c = 0; c < w; ++c) {
                Wd[(size_t)c * H] -= tr[c];
              }
            }
          }
        }
      }
      kls_snb_spin_barrier_wait(&ctx->coop, crew);
      continue;
    }
    /* big edge: stage A on tid 0 (gather + TRSM + writeback + row
       mirror into ub_w[0]), stage B row-sliced across the crew */
    double *B = st->ub_w[0];
    double *Brow = st->ub_w[0] + (size_t)sc * w;
    if (tid == 0) {
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
        for (int32_t i = 0; i < sc; ++i) {
          Wc[i] = Bc[i];
          Brow[(size_t)i * w + c] = Bc[i];
        }
      }
    }
    kls_snb_spin_barrier_wait(&ctx->coop, crew);
    {
      const int32_t chunk = (ph + crew - 1) / crew;
      const int32_t r0 = tid * chunk;
      int32_t r1 = r0 + chunk;
      if (r1 > ph) r1 = ph;
      if (r0 < r1) {
        double *T = st->gemm_w[tid];
        {
          const double *lb = Lbelow + (size_t)sub[0] * pld;
          for (int32_t r = r0; r < r1; ++r) {
            const double l = lb[r];
            double *tr = T + (size_t)(r - r0) * w;
            for (int32_t c = 0; c < w; ++c) tr[c] = l * Brow[c];
          }
          for (int32_t i = 1; i < sc; ++i) {
            const double *lbi = Lbelow + (size_t)sub[i] * pld;
            const double *bri = Brow + (size_t)i * w;
            for (int32_t r = r0; r < r1; ++r) {
              const double l = lbi[r];
              double *tr = T + (size_t)(r - r0) * w;
              for (int32_t c = 0; c < w; ++c) tr[c] += l * bri[c];
            }
          }
        }
        for (int32_t r = r0; r < r1; ++r) {
          const int32_t dst = dsts[r];
          const double *tr = T + (size_t)(r - r0) * w;
          if (dst >= 0) {
            double *Wd = W + dst;
            for (int32_t c = 0; c < w; ++c) {
              Wd[(size_t)c * H] -= tr[c];
            }
          }
        }
      }
    }
    kls_snb_spin_barrier_wait(&ctx->coop, crew);
  }
  /* fused getrf: i-steps row-sliced with a barrier per step */
  {
    double *diag = W + uw;
    for (int32_t i = 0; i < w; ++i) {
      double *ci = diag + (size_t)i * H;
      const double pivot = ci[i];
      if (tid == 0) {
        udiag[k1 + sn->start + i] = pivot;
        if (pivot == 0.0) {
          kls_snb_record_singular(st, solver, k1 + sn->start + i);
        }
      }
      const int32_t hrem = H - uw - i - 1;
      const int32_t chunk = (hrem + crew - 1) / crew;
      const int32_t r0 = tid * chunk;
      int32_t r1 = r0 + chunk;
      if (r1 > hrem) r1 = hrem;
      double *lsub = ci + i + 1;
      for (int32_t r = r0; r < r1; ++r) lsub[r] /= pivot;

      kls_snb_spin_barrier_wait(&ctx->coop, crew);
      if (tid == 0) {
        atomic_store_explicit(&blk->published[s], i + 1,
                              memory_order_release);
      }
      for (int32_t j = i + 1; j < w; ++j) {
        double *cj = diag + (size_t)j * H;
        const double u = cj[i];
        if (u == 0.0) continue;
        double *dstv = cj + i + 1;
        for (int32_t r = r0; r < r1; ++r) {
          dstv[r] -= lsub[r] * u;
        }
      }
      kls_snb_spin_barrier_wait(&ctx->coop, crew);
    }
  }
  /* writeback: column-sliced */
  {
    trilinos_klu_l_numeric *numeric = solver->numeric;
    const UF_long *lip = numeric->Lip + k1;
    const UF_long *llen = numeric->Llen + k1;
    const UF_long *uip = numeric->Uip + k1;
    const UF_long *ulen = numeric->Ulen + k1;
    double *lu = (double *)numeric->LUbx[blk - st->blocks];
    const double *base = panels + sn->panel;
    for (int64_t k = sn->start + tid; k < sn->start + sn->width;
         k += crew) {
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
  kls_snb_spin_barrier_wait(&ctx->coop, crew);
}

typedef struct kls_snb_par_arg {
  kls_snb_par_ctx *ctx;
  int tid;
} kls_snb_par_arg;

static void *kls_snb_par_worker(void *argp) {
  kls_snb_par_arg *arg = (kls_snb_par_arg *)argp;
  kls_snb_par_ctx *ctx = arg->ctx;
  kls_snb_state *st = ctx->st;
  int go;
  while ((go = atomic_load_explicit(&ctx->go, memory_order_acquire)) == 0) {
    kls_cpu_relax();
  }
  if (go != 1) {
    return NULL;
  }
  for (int64_t seg = 0; seg < st->seg_count; ++seg) {
    const int64_t begin = st->seg_ptr[2 * seg];
    const int64_t end = st->seg_ptr[2 * seg + 1];
    if (st->seg_coop[seg]) {
      /* single big snode: whole crew cooperates (rendezvous barrier) */
      (void)pthread_barrier_wait(&ctx->barrier);
      kls_snb_block *blk = &st->blocks[st->slevel_block[begin]];
      const int64_t s = (int64_t)st->slevel_snode[begin];
      kls_snb_process_target_coop(ctx, blk, s, arg->tid);
      if (arg->tid == 0) {
        atomic_store_explicit(&blk->published[s],
                              blk->snodes[s].width,
                              memory_order_release);
      }
      (void)pthread_barrier_wait(&ctx->barrier);
      continue;
    }
    /* pipelined segment: claims in topological order, dependencies
       enforced by per-edge spins on published getrf columns */
    if (arg->tid == 0) {
      atomic_store_explicit(&ctx->cursor, begin, memory_order_relaxed);
    }
    (void)pthread_barrier_wait(&ctx->barrier);
    for (;;) {
      const int64_t idx = atomic_fetch_add_explicit(&ctx->cursor, 1,
                                                    memory_order_relaxed);
      if (idx >= end) {
        break;
      }
      kls_snb_block *blk = &st->blocks[st->slevel_block[idx]];
      kls_snb_process_target(ctx->solver, st, blk,
                             (int64_t)st->slevel_snode[idx], ctx->values,
                             ctx->rs, ctx->scale, st->ub_w[arg->tid],
                             st->gemm_w[arg->tid], 1);
    }
    (void)pthread_barrier_wait(&ctx->barrier);
  }
  return NULL;
}

/* Level-scheduled parallel supernodal refactor: workers claim whole
 * snodes within a level of the snode DAG (blocks merged; they are
 * independent), with immediate per-snode writeback. */
static UF_long kls_snb_refactor_parallel(kls_solver *solver,
                                         const double *numeric_values,
                                         int nthreads) {
  kls_snb_state *st = solver->snb;
  if (st->slevel_ptr == NULL || st->slevel_count <= 0 ||
      st->seg_ptr == NULL || st->seg_count <= 0) {
    return kls_snb_refactor(solver, numeric_values);
  }
  for (UF_long b = 0; b < st->nblocks; ++b) {
    kls_snb_block *blk = &st->blocks[b];
    if (blk->published == NULL) continue;
    for (int64_t s = 0; s < blk->snode_count; ++s) {
      atomic_store_explicit(&blk->published[s], 0, memory_order_relaxed);
    }
  }
  if (nthreads > st->nworkers) {
    nthreads = st->nworkers;
  }
  if (nthreads < 2) {
    return kls_snb_refactor(solver, numeric_values);
  }
  if (!kls_snb_refactor_prelude(solver, numeric_values)) {
    return 0;
  }
  kls_snb_par_ctx ctx;
  ctx.solver = solver;
  ctx.st = st;
  ctx.values = numeric_values;
  ctx.rs = solver->numeric->Rs;
  ctx.scale = (int)solver->common.scale;
  ctx.nthreads = nthreads;
  atomic_init(&ctx.cursor, 0);
  atomic_init(&ctx.go, 0);
  atomic_init(&ctx.coop.arrived, 0);
  atomic_init(&ctx.coop.generation, 0u);
  if (pthread_barrier_init(&ctx.barrier, NULL, (unsigned)nthreads) != 0) {
    return kls_snb_refactor(solver, numeric_values);
  }
  pthread_t tids[64];
  kls_snb_par_arg args[64];
  int spawned = 0;
  for (int t = 1; t < nthreads; ++t) {
    args[t].ctx = &ctx;
    args[t].tid = t;
    if (pthread_create(&tids[t], NULL, kls_snb_par_worker, &args[t]) != 0) {
      break;
    }
    spawned = t;
  }
  if (spawned != nthreads - 1) {
    /* could not build the crew: release the spawned workers without
       letting anyone touch the barrier, then run serially */
    atomic_store_explicit(&ctx.go, 2, memory_order_release);
    for (int t = 1; t <= spawned; ++t) {
      pthread_join(tids[t], NULL);
    }
    pthread_barrier_destroy(&ctx.barrier);
    return kls_snb_refactor(solver, numeric_values);
  }
  atomic_store_explicit(&ctx.go, 1, memory_order_release);
  args[0].ctx = &ctx;
  args[0].tid = 0;
  (void)kls_snb_par_worker(&args[0]);
  for (int t = 1; t <= spawned; ++t) {
    pthread_join(tids[t], NULL);
  }
  pthread_barrier_destroy(&ctx.barrier);
  if (atomic_load_explicit(&st->failed, memory_order_acquire)) {
    return 0;
  }
  return kls_snb_refactor_epilogue(solver);
}

static int kls_snb_env_disabled(void) {
  const char *value = getenv("KLS_DISABLE_SNB_REFACTOR");
  return value != NULL && value[0] == '1';
}

/* Run the supernodal acceptance at factor time so its one-time A/B cost
 * lands in the once-charged factor term instead of inflating the first
 * (99x-charged) refactor.  Also lets snb sort the numeric below the
 * snode-panel flops floor: panels have no batch-run floor of their own,
 * and TSOPF-class fused-block rows live down there. */
static void kls_snb_maybe_accept(kls_solver *solver,
                                 const double *numeric_values,
                                 double *elapsed) {
  if (solver == NULL || solver->snb_decision != 0 || solver->snb_declined ||
      solver->snb_trial_verdict < 0 ||
      solver->snb != NULL || kls_snb_env_disabled() ||
      solver->symbolic == NULL ||
      solver->numeric == NULL || solver->numeric->LUbx == NULL ||
      (solver->fp32_decision > 0) ||
      solver->common.status == TRILINOS_KLU_SINGULAR) {
    return;
  }
  if (solver->symbolic->maxblock > 1000000u &&
      getenv("KLS_SNB_FORCE_TRIAL") == NULL) {
    /* SNB's retained schedule is deliberately rich: it performs several
       max-block-sized passes and materializes per-supernode unions and edge
       maps.  On million-column blocks that setup cannot amortize over H100;
       worse, sparse giant circuits can spend the entire process timeout in
       preparation before reaching a single refactor (circuit5M and rajat31).
       The union audit's largest SNB adopter has a 321K-column block. */
    solver->snb_declined = 1;
    solver->snb_decision = -1;
    return;
  }
  if (solver->options.backend == KLS_BACKEND_SERIAL &&
      getenv("KLS_SNB_FORCE_TRIAL") == NULL) {
    /* Keep the opt-in serial backend lean on cache-resident factors.  Its
       mapped refactor needs only a cheap position map, whereas an SNB A/B
       trial materializes panels and executes an extra refactor.  Once the
       numeric value arrays alone exceed 32 MiB, the mapped walk is firmly
       memory-bound.  Require a pivot-stable numeric and a dominant BTF block
       as well: otherwise panel materialization is scattered work and its
       trial is not safely amortizable.  These are factor-shape invariants,
       independent of matrix names or benchmark membership. */
    const double factor_value_bytes =
      ((double)solver->numeric->lnz + (double)solver->numeric->unz) *
      (double)sizeof(double);
    if (factor_value_bytes < 32.0 * 1024.0 * 1024.0 ||
        solver->common.noffdiag != 0u ||
        (double)solver->symbolic->maxblock < 0.75 * (double)solver->n) {
      return;
    }
  }
  const double accept_start = kls_now_seconds();
  if (getenv("KLS_TRACE_SNB") != NULL) {
    fprintf(stderr,
            "KLS snb gate: prep=%d decl=%d verd=%d dec=%d fp32=%d wall=%.1f\n",
            solver->snode_prepared, solver->snb_declined,
            solver->snb_trial_verdict, solver->snb_decision,
            solver->fp32_decision,
            solver->snb_factor_start > 0.0
              ? kls_now_seconds() - solver->snb_factor_start : -1.0);
  }
  if (!solver->snode_prepared) {
    /* relaxed sort: same safety guards as the snode-panel prep, lower
       size floor (no retained position-based structures may exist) */
    if (solver->n < 2000u || solver->common.flops < 5.0e5 ||
        solver->refactor_col_ptr != NULL ||
        solver->refactor_level_ptr != NULL ||
        solver->refactor_l_indices32 != NULL ||
        solver->row_refactor_group_count != 0u ||
        solver->row_refactor_values_ready ||
        solver->stats.last_factor_path == KLS_FACTOR_PATH_KLS_FIRST) {
      *elapsed += kls_now_seconds() - accept_start;
      return;
    }
    solver->snode_prepared = 1;
    if (!solver->snode_numeric_pre_sorted &&
        !kls_parallel_lu_sort(solver)) {
      solver->snb_declined = 1;
      *elapsed += kls_now_seconds() - accept_start;
      return;
    }
    solver->snode_numeric_pre_sorted = 1;
  }
  if (solver->options.threads > 1 &&
      getenv("KLS_SNB_FORCE_TRIAL") == NULL &&
      (solver->snb_factor_start <= 0.0 ||
       kls_now_seconds() - solver->snb_factor_start <
         (30.0))) {
    /* t4: prep+trial cost only amortizes on matrices whose factor is
       already expensive (the adopters measure 18s+); cheap-factor rows
       are the pipelined incumbent's territory and the trial would tax
       the once-charged term for nothing.  Not a permanent decline: a
       later accept site in the same (or a future) factor may sit past
       the threshold. */
    *elapsed += kls_now_seconds() - accept_start;
    return;
  }
  solver->snb_trial_budget = solver->options.threads > 1 ? 1.0 : 0.0;
  if (!kls_snb_prepare(solver)) {
    if (solver->snb_declined) {
      solver->snb_decision = -1;
    }
    *elapsed += kls_now_seconds() - accept_start;
    return;
  }
  solver->snb_incumbent_seconds = 0.0;
  solver->snb_trial_seconds = 0.0;
  if (numeric_values == NULL) {
    numeric_values = solver->values;
  }
  if (numeric_values != NULL &&
      solver->snb_trial_stamp == solver->snb_factor_start) {
    /* a numeric rebuild inside this same kls_factor already paid for a
       trial; do not pay again (the refactor-time hooks simply see no
       engine and fall through to the incumbent) */
    numeric_values = NULL;
  }
  if (numeric_values != NULL) {
    /* run the engine once now: the one-time trial cost lands in the
       once-charged factor term instead of the 99x-charged refactor
       average; the incumbent is captured on the first refactor call and
       the decision falls on the second. */
    const int threads = solver->options.threads;
    const double t0 = kls_now_seconds();
    const UF_long snb_ok =
      threads > 1
        ? kls_snb_refactor_parallel(solver, numeric_values, threads)
        : kls_snb_refactor(solver, numeric_values);
    const double trial = kls_now_seconds() - t0;
    solver->snb_trial_stamp = solver->snb_factor_start;
    if (snb_ok) {
      solver->snb_trial_seconds = trial;
    } else {
      solver->snb_decision = -1;
      kls_snb_free(solver);
      solver->snb_declined = 1;
    }
    if (kls_snb_trace_enabled()) {
      fprintf(stderr, "KLS snb: factor-time trial %.3fms t%d ok=%ld\n",
              trial * 1e3, threads, (long)snb_ok);
    }
  }
  *elapsed += kls_now_seconds() - accept_start;
}

/* Try the supernodal refactor.  Returns 1 when it handled the call (with
 * *ok_out as the refactor status), 0 when the caller should proceed with
 * the incumbent path.  On the first eligible call the incumbent serial
 * path runs first and both are timed; the engine is adopted per matrix
 * only on a clear win. */
static int kls_snb_try_refactor(kls_solver *solver,
                                const double *numeric_values,
                                int check_pivots,
                                UF_long *ok_out) {
  if (!check_pivots && solver->snb_retrial == 1 && solver->snb != NULL &&
      solver->snb_decision < 0 && !kls_snb_env_disabled() &&
      solver->numeric != NULL && solver->symbolic != NULL &&
      ++solver->snb_retrial_wait >= 8 &&
      solver->adaptive_refactor_seconds > 0.0) {
    /* steady re-trial of a near-miss rejection: the incumbent's rate
       is now its true steady (updated every refactor) */
    solver->snb_retrial = 2;
    const double t0 = kls_now_seconds();
    const UF_long rok = solver->options.threads > 1
      ? kls_snb_refactor_parallel(solver, numeric_values,
                                  solver->options.threads)
      : kls_snb_refactor(solver, numeric_values);
    const double dt = kls_now_seconds() - t0;
    if (kls_snb_trace_enabled()) {
      fprintf(stderr,
              "KLS snb steady re-trial: snb %.3fms vs incumbent %.3fms"
              " -> %s\n", dt * 1e3,
              solver->adaptive_refactor_seconds * 1e3,
              rok && dt < 0.85 * solver->adaptive_refactor_seconds
                ? "adopted" : "rejected");
    }
    if (rok && dt < 0.85 * solver->adaptive_refactor_seconds) {
      solver->snb_decision = 1;
      *ok_out = rok;
      kls_set_last_refactor_path(solver, KLS_REFACTOR_PATH_SNB);
      return 1;
    }
    solver->snb_trial_verdict = -1;
    kls_snb_free(solver);
    solver->snb_declined = 1;
    if (rok) {
      /* the re-trial computed valid values; use them for this call */
      *ok_out = rok;
      kls_set_last_refactor_path(solver, KLS_REFACTOR_PATH_SNB);
      return 1;
    }
    return 0;
  }
  if (check_pivots || kls_snb_env_disabled() ||
      solver->numeric == NULL ||
      solver->symbolic == NULL || solver->snb == NULL ||
      solver->snb_decision < 0 || solver->snb_trial_verdict < 0) {
    /* No fp32 conjunct here: a fp64 snb build is strictly more accurate
       and the 0.60 timed acceptance already arbitrates speed (pre2: the
       fp32-decided pipeline sat at 8.5s while the fp64 snb trial
       measured 3.1s and was vetoed by the flag alone). */
    if (kls_snb_trace_enabled()) {
      fprintf(stderr,
              "KLS snb try gate: cp=%d dis=%d num=%d sym=%d snb=%d fp32=%d"
              " dec=%d verd=%d\n",
              check_pivots, kls_snb_env_disabled(),
              solver->numeric == NULL, solver->symbolic == NULL,
              solver->snb == NULL,
              solver->fp32_decision > 0,
              solver->snb_decision, solver->snb_trial_verdict);
    }
    return 0;
  }
  const int threads = solver->options.threads;
  if (solver->snb_decision == 0) {
    /* Event-driven acceptance: the incumbent path runs (and is timed by
       kls_refactor's wrapper) on one call, the engine runs timed on the
       next, and the better one keeps the matrix.  No call ever runs two
       engines, so incumbent statistics stay untouched and the one-time
       cost is bounded by a single refactor. */
    if (solver->snb_incumbent_seconds == 0.0) {
      solver->snb_incumbent_seconds = -1.0;
      return 0; /* let the incumbent run this call */
    }
    if (solver->snb_incumbent_seconds < 0.0) {
      if (solver->adaptive_refactor_seconds <= 0.0) {
        return 0; /* wrapper timing unavailable; keep waiting */
      }
      solver->snb_incumbent_seconds = solver->adaptive_refactor_seconds;
    }
    UF_long snb_ok = 1;
    double snb_seconds = solver->snb_trial_seconds;
    if (snb_seconds <= 0.0) {
      /* no factor-time trial (values not retained there): trial now */
      const double t0 = kls_now_seconds();
      snb_ok =
        threads > 1 ? kls_snb_refactor_parallel(solver, numeric_values,
                                                threads)
                    : kls_snb_refactor(solver, numeric_values);
      snb_seconds = kls_now_seconds() - t0;
    }
    /* Single-sample timings on a busy box swing ~20-30%.  The parallel
       policy therefore retains its conservative 0.60 margin.  The serial
       backend has already passed a narrow structural gate before paying for
       this trial, so it retains any clear measured win (0.90). */
    const double adoption_ratio =
      solver->options.backend == KLS_BACKEND_SERIAL ? 0.90 : 0.60;
    solver->snb_decision =
      snb_ok && snb_seconds < adoption_ratio * solver->snb_incumbent_seconds
        ? 1 : -1;
    if (kls_snb_trace_enabled()) {
      fprintf(stderr,
              "KLS snb: acceptance incumbent %.3fms snb %.3fms t%d -> %s\n",
              solver->snb_incumbent_seconds * 1e3, snb_seconds * 1e3,
              threads, solver->snb_decision > 0 ? "adopted" : "rejected");
    }
    if (solver->snb_decision < 0) {
      if (snb_ok && snb_seconds < solver->snb_incumbent_seconds &&
          solver->snb_retrial == 0) {
        /* Near miss: the engine measured absolutely faster but not by
           the conservative margin - and both samples are warm-up era
           (bcircuit: trial incumbent 4.9ms vs its 1.57ms steady). Keep
           the engine and re-run the comparison once at steady state,
           when adaptive_refactor_seconds carries the incumbent's true
           rate. */
        solver->snb_retrial = 1;
        solver->snb_retrial_wait = 0;
      } else {
        solver->snb_trial_verdict = -1;
        kls_snb_free(solver);
        solver->snb_declined = 1;
      }
      if (solver->snb_trial_seconds > 0.0 || !snb_ok) {
        return 0; /* nothing computed this call; incumbent runs it */
      }
    }
    if (solver->snb_decision > 0 && solver->snb_trial_seconds > 0.0) {
      /* decided from the stored trial: compute this call with the engine */
      snb_ok =
        threads > 1 ? kls_snb_refactor_parallel(solver, numeric_values,
                                                threads)
                    : kls_snb_refactor(solver, numeric_values);
      if (!snb_ok) {
        solver->snb_decision = -1;
        kls_snb_free(solver);
        solver->snb_declined = 1;
        return 0;
      }
    }
    *ok_out = snb_ok;
    kls_set_last_refactor_path(solver, KLS_REFACTOR_PATH_SNB);
    return 1;
  }
  *ok_out = threads > 1
    ? kls_snb_refactor_parallel(solver, numeric_values, threads)
    : kls_snb_refactor(solver, numeric_values);
  kls_set_last_refactor_path(solver, KLS_REFACTOR_PATH_SNB);
  return 1;
}

static int kls_snode_panel_env_disabled(void) {
  const char *value = getenv("KLS_DISABLE_SNODE_PANEL_REFACTOR");
  return value != NULL && value[0] == '1';
}

/* ------------------------------------------------------------------ */
/* Parallel packed-LU column sort (hard core D, brick 1).  The         */
/* vendored TRILINOS_KLU_sort is a serial double transpose; columns    */
/* sort independently, so a crew with an atomic column cursor does the */
/* same job at memory speed.  Semantics identical: every L and U       */
/* column of every block ends ascending by row index with values       */
/* permuted along.                                                     */
/* ------------------------------------------------------------------ */

typedef struct kls_par_sort_pair {
  UF_long row;
  double val;
} kls_par_sort_pair;

static int kls_par_sort_pair_cmp(const void *a, const void *b) {
  const UF_long ra = ((const kls_par_sort_pair *)a)->row;
  const UF_long rb = ((const kls_par_sort_pair *)b)->row;
  return ra < rb ? -1 : (ra > rb ? 1 : 0);
}

static int kls_par_sort_column(UF_long *xi, double *xx, UF_long len,
                               kls_par_sort_pair **buf, UF_long *cap) {
  if (len < 2) {
    return 1;
  }
  int sorted = 1;
  for (UF_long p = 1; p < len; ++p) {
    if (xi[p - 1] > xi[p]) {
      sorted = 0;
      break;
    }
  }
  if (sorted) {
    return 1;
  }
  if (len <= 40u) {
    for (UF_long p = 1; p < len; ++p) {
      const UF_long ri = xi[p];
      const double vi = xx[p];
      UF_long q = p;
      while (q > 0 && xi[q - 1] > ri) {
        xi[q] = xi[q - 1];
        xx[q] = xx[q - 1];
        q--;
      }
      xi[q] = ri;
      xx[q] = vi;
    }
    return 1;
  }
  if (*cap < len) {
    UF_long next = *cap > 0 ? *cap : 4096u;
    while (next < len) next *= 2u;
    kls_par_sort_pair *grown =
      (kls_par_sort_pair *)realloc(*buf, (size_t)next * sizeof(**buf));
    if (grown == NULL) {
      return 0;
    }
    *buf = grown;
    *cap = next;
  }
  kls_par_sort_pair *pairs = *buf;
  for (UF_long p = 0; p < len; ++p) {
    pairs[p].row = xi[p];
    pairs[p].val = xx[p];
  }
  qsort(pairs, (size_t)len, sizeof(*pairs), kls_par_sort_pair_cmp);
  for (UF_long p = 0; p < len; ++p) {
    xi[p] = pairs[p].row;
    xx[p] = pairs[p].val;
  }
  return 1;
}

typedef struct kls_par_sort_ctx {
  trilinos_klu_l_symbolic *symbolic;
  trilinos_klu_l_numeric *numeric;
  _Atomic int64_t cursor;   /* global column index across all blocks */
  _Atomic int failed;
} kls_par_sort_ctx;

static void *kls_par_sort_worker(void *argp) {
  kls_par_sort_ctx *ctx = (kls_par_sort_ctx *)argp;
  trilinos_klu_l_symbolic *symbolic = ctx->symbolic;
  trilinos_klu_l_numeric *numeric = ctx->numeric;
  const UF_long n = symbolic->n;
  kls_par_sort_pair *buf = NULL;
  UF_long cap = 0;
  for (;;) {
    const int64_t start =
      atomic_fetch_add_explicit(&ctx->cursor, 256, memory_order_relaxed);
    if (start >= (int64_t)n) {
      break;
    }
    int64_t end = start + 256;
    if (end > (int64_t)n) end = (int64_t)n;
    for (int64_t gk = start; gk < end; ++gk) {
      /* locate the block (blocks are few; linear from R would be slow —
         binary search) */
      UF_long lo = 0, hi = (UF_long)symbolic->nblocks;
      while (lo + 1u < hi) {
        const UF_long mid = (lo + hi) >> 1;
        if (symbolic->R[mid] <= (UF_long)gk) lo = mid;
        else hi = mid;
      }
      const UF_long block = lo;
      const UF_long k1 = symbolic->R[block];
      const UF_long nk = symbolic->R[block + 1u] - k1;
      if (nk < 2u || numeric->LUbx[block] == NULL) {
        continue;
      }
      double *lu = (double *)numeric->LUbx[block];
      const UF_long k = (UF_long)gk - k1;
      UF_long *xi; double *xx; UF_long len;
      kls_klu_get_pointer(lu, numeric->Lip + k1, numeric->Llen + k1, k,
                          &xi, &xx, &len);
      if (!kls_par_sort_column(xi, xx, len, &buf, &cap)) {
        atomic_store_explicit(&ctx->failed, 1, memory_order_release);
        break;
      }
      kls_klu_get_pointer(lu, numeric->Uip + k1, numeric->Ulen + k1, k,
                          &xi, &xx, &len);
      if (!kls_par_sort_column(xi, xx, len, &buf, &cap)) {
        atomic_store_explicit(&ctx->failed, 1, memory_order_release);
        break;
      }
    }
    if (atomic_load_explicit(&ctx->failed, memory_order_acquire)) {
      break;
    }
  }
  free(buf);
  return NULL;
}

/* Drop-in for trilinos_klu_l_sort: parallel when it can be, vendored
 * serial otherwise.  Returns nonzero on success like the original. */
static UF_long kls_parallel_lu_sort(kls_solver *solver) {
  /* The compact solve streams copy the row-index order of each packed L/U
     column.  A deferred prep may reach this sort after an earlier solve has
     already built those streams, so invalidate them at the mutation boundary
     rather than relying on every sort caller to predict that lifecycle. */
  kls_invalidate_i32_solve(solver);
  int nt = solver->options.threads;
  if (nt > 16) nt = 16;
  if (nt < 2 ||
      (solver->n < 100000u) ||
      getenv("KLS_DISABLE_PARALLEL_SORT") != NULL) {
    return trilinos_klu_l_sort(solver->symbolic, solver->numeric,
                               &solver->common);
  }
  kls_par_sort_ctx ctx;
  ctx.symbolic = solver->symbolic;
  ctx.numeric = solver->numeric;
  atomic_init(&ctx.cursor, 0);
  atomic_init(&ctx.failed, 0);
  pthread_t tids[16];
  int spawned = 0;
  for (int t = 1; t < nt; ++t) {
    if (pthread_create(&tids[t], NULL, kls_par_sort_worker, &ctx) != 0) {
      break;
    }
    spawned = t;
  }
  (void)kls_par_sort_worker(&ctx);
  for (int t = 1; t <= spawned; ++t) {
    pthread_join(tids[t], NULL);
  }
  if (atomic_load_explicit(&ctx.failed, memory_order_acquire)) {
    /* allocation failure mid-sort leaves some columns sorted, some not;
       the vendored sort is order-insensitive, run it to finish */
    return trilinos_klu_l_sort(solver->symbolic, solver->numeric,
                               &solver->common);
  }
  return 1;
}

/* B'-CORE stage 2: dense consume over a padded run. The consumer's U
   positions must hold the run contiguously from producer j to the run
   end. Multipliers are computed sequentially using the panel's
   within-run slots (the sorted union's first (len-1) rows are exactly
   s+1..e by the chain property, so row r's slot is r-s-1), then the
   tail accumulates chunked over contiguous panel rows - the strict
   batch kernel's shape at relaxed coverage, pad slots contributing
   exact zeros. Returns producers consumed or 0. */
#if defined(__GNUC__) && defined(__x86_64__) && !defined(__clang__)
__attribute__((target_clones("default", "arch=x86-64-v3", "arch=x86-64-v4")))
#endif
static UF_long kls_padded_run_consume(const kls_solver *ps,
                                      UF_long k1,
                                      UF_long j,
                                      const UF_long *ui,
                                      double *ux,
                                      UF_long ucol_len,
                                      UF_long up,
                                      double *restrict x) {
  const UF_long run1 = ps->padded_run_of[k1 + j];
  if (run1 == 0u) {
    return 0;
  }
  const UF_long run = run1 - 1u;
  const UF_long s = ps->padded_run_start[run] - k1;
  const UF_long len = ps->padded_run_len[run];
  const UF_long e = s + len - 1u;
  UF_long tmax = e - j + 1u;
  if (tmax > ucol_len - up) {
    tmax = ucol_len - up;
  }
  if (tmax > 64u) {
    tmax = 64u;
  }
  UF_long t = 1;
  while (t < tmax && ui[up + t] == j + t) {
    t++;
  }
  if (t < 2u) {
    return 0;
  }
  const UF_long ubase = ps->padded_union_ptr[run];
  const UF_long ulen_run = ps->padded_union_ptr[run + 1u] - ubase;
  if (t * ulen_run < 192u) {
    return 0;
  }
  const UF_long *urows = ps->padded_union_rows + ubase;
  const double *panel =
    ps->padded_panel_values + ps->padded_panel_ptr[run];
  const UF_long within = len - 1u; /* union rows s+1..e occupy [0,within) */
  double xs[64];
  for (UF_long i = 0; i < t; ++i) {
    xs[i] = x[j + i];
    x[j + i] = 0.0;
  }
  for (UF_long i = 0; i < t; ++i) {
    const double u = xs[i];
    ux[up + i] = u;
    if (u == 0.0) {
      continue;
    }
    const double *row = panel + (j + i - s) * ulen_run;
    /* within-run updates: slots of rows j+i+1..e are contiguous;
       rows beyond the consumed prefix update the live workspace */
    for (UF_long r = j + i + 1u; r <= e; ++r) {
      const double v = u * row[r - s - 1u];
      if (r < j + t) {
        xs[r - j] -= v;
      } else {
        x[r] -= v;
      }
    }
  }
  /* shared tail: union slots [within, ulen_run) */
  for (UF_long q0 = within; q0 < ulen_run; q0 += 32u) {
    const UF_long qc = ulen_run - q0 < 32u ? ulen_run - q0 : 32u;
    double acc[32];
    for (UF_long q = 0; q < qc; ++q) {
      acc[q] = 0.0;
    }
    for (UF_long i = 0; i < t; ++i) {
      const double u = xs[i];
      if (u == 0.0) {
        continue;
      }
      const double *row = panel + (j + i - s) * ulen_run + q0;
      for (UF_long q = 0; q < qc; ++q) {
        acc[q] += row[q] * u;
      }
    }
    for (UF_long q = 0; q < qc; ++q) {
      x[urows[q0 + q]] -= acc[q];
    }
  }
  return t;
}

/* B'-CORE stage 1 build: relaxed runs (column k+1 within z symmetric-
   difference entries of column k's tail, chain entry required), their
   sorted union patterns, per-column slot maps, and zeroed panel value
   storage. Write-side only: the mapped kernel refreshes panel rows at
   column finalize; the stage-2 consume reads them. Census (09bd408):
   bcircuit z=2 covers 79% of columns at +9.9%% pad. */
static void kls_build_padded_panels(kls_solver *solver) {
  if (solver == NULL || solver->padded_run_of != NULL ||
      (getenv("KLS_ENABLE_PADDED_PANELS") == NULL &&
       !solver->padded_probe_build) ||
      solver->symbolic == NULL || solver->numeric == NULL) {
    return;
  }
  const char *zenv = getenv("KLS_PADDED_Z");
  const UF_long z = zenv != NULL ? (UF_long)atol(zenv) : 2u;
  const UF_long n = solver->n;
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  const UF_long union_cap = 4096u;
  UF_long *run_of = (UF_long *)calloc((size_t)n, sizeof(*run_of));
  UF_long *slot_ptr = (UF_long *)calloc((size_t)n + 1u, sizeof(*slot_ptr));
  UF_long *mark = (UF_long *)calloc((size_t)n, sizeof(*mark));
  UF_long *ulist = (UF_long *)malloc((size_t)n * sizeof(*ulist));
  if (run_of == NULL || slot_ptr == NULL || mark == NULL || ulist == NULL) {
    goto fail;
  }
  {
    /* pass A: detect runs, count unions/slots/panel size */
    UF_long nruns = 0, union_total = 0, slot_total = 0;
    double panel_total = 0.0;
    const double lnz_total = (double)numeric->lnz + 1.0;
    for (UF_long block = 0; block < symbolic->nblocks; ++block) {
      const UF_long k1 = symbolic->R[block];
      const UF_long nk = symbolic->R[block + 1u] - k1;
      if (nk < 2u) {
        continue;
      }
      UF_long *lip = numeric->Lip + k1;
      UF_long *llen = numeric->Llen + k1;
      double *lu = (double *)numeric->LUbx[block];
      if (lu == NULL) {
        continue;
      }
      UF_long s = 0;
      for (UF_long k = 0; k < nk; ++k) {
        int ext = 0;
        if (k + 1u < nk && llen[k] >= 1) {
          UF_long *li, *li2;
          double *lx, *lx2;
          UF_long l1, l2;
          kls_klu_get_pointer(lu, lip, llen, k, &li, &lx, &l1);
          kls_klu_get_pointer(lu, lip, llen, k + 1u, &li2, &lx2, &l2);
          UF_long a = li[0] == k + 1u ? 1u : 0u, b = 0, mism = 0;
          if (li[0] == k + 1u) {
            while (a < l1 && b < l2) {
              if (li[a] == li2[b]) {
                a++, b++;
              } else if (li[a] < li2[b]) {
                a++, mism++;
              } else {
                b++, mism++;
              }
            }
            mism += (l1 - a) + (l2 - b);
            ext = mism <= z;
          }
        }
        if (!ext) {
          if (k > s) {
            /* run [s..k]: union via marks */
            UF_long ulen_run = 0;
            int over = 0;
            for (UF_long c = s; c <= k && !over; ++c) {
              UF_long *li;
              double *lx;
              UF_long l1;
              kls_klu_get_pointer(lu, lip, llen, c, &li, &lx, &l1);
              for (UF_long p = 0; p < l1; ++p) {
                if (mark[li[p]] == 0u) {
                  mark[li[p]] = 1u;
                  if (ulen_run >= union_cap) {
                    over = 1;
                    break;
                  }
                  ulist[ulen_run++] = li[p];
                }
              }
            }
            for (UF_long u = 0; u < ulen_run; ++u) {
              mark[ulist[u]] = 0u;
            }
            if (!over) {
              const UF_long len = k - s + 1u;
              static UF_long min_work = 0;
              if (min_work == 0u) {
                const char *we = getenv("KLS_PADDED_MIN_WORK");
                min_work = we != NULL ? (UF_long)atol(we) : 1024u;
              }
              /* materialize only runs whose consume work amortizes
                 the per-refactor panel refresh (the refresh tax
                 cancels scrap-coverage wins: bcircuit all-runs read
                 +20%) */
              if (len * ulen_run >= min_work) {
                panel_total += (double)len * (double)ulen_run;
                if (panel_total > 2.5 * lnz_total) {
                  goto fail; /* padding economics out of budget */
                }
                nruns++;
                union_total += ulen_run;
                for (UF_long c = s; c <= k; ++c) {
                  run_of[k1 + c] = nruns; /* 1-based */
                  slot_total += llen[c];
                }
              }
            }
          }
          s = k + 1u;
        }
      }
    }
    if (nruns == 0u) {
      goto fail;
    }
    solver->padded_run_count = nruns;
    solver->padded_run_start =
      (UF_long *)malloc((size_t)nruns * sizeof(UF_long));
    solver->padded_run_len =
      (UF_long *)malloc((size_t)nruns * sizeof(UF_long));
    solver->padded_union_ptr =
      (UF_long *)calloc((size_t)nruns + 1u, sizeof(UF_long));
    solver->padded_union_rows =
      (UF_long *)malloc((size_t)union_total * sizeof(UF_long));
    solver->padded_slots =
      (UF_long *)malloc((size_t)slot_total * sizeof(UF_long));
    solver->padded_panel_ptr =
      (UF_long *)calloc((size_t)nruns + 1u, sizeof(UF_long));
    solver->padded_panel_values =
      (double *)calloc((size_t)panel_total, sizeof(double));
    if (solver->padded_run_start == NULL || solver->padded_run_len == NULL ||
        solver->padded_union_ptr == NULL ||
        solver->padded_union_rows == NULL || solver->padded_slots == NULL ||
        solver->padded_panel_ptr == NULL ||
        solver->padded_panel_values == NULL) {
      goto fail_arrays;
    }
    /* pass B: fill (same detection; runs enumerate identically) */
    UF_long run = 0, upos = 0, spos = 0;
    double ppos = 0.0;
    for (UF_long block = 0; block < symbolic->nblocks; ++block) {
      const UF_long k1 = symbolic->R[block];
      const UF_long nk = symbolic->R[block + 1u] - k1;
      if (nk < 2u) {
        continue;
      }
      UF_long *lip = numeric->Lip + k1;
      UF_long *llen = numeric->Llen + k1;
      double *lu = (double *)numeric->LUbx[block];
      if (lu == NULL) {
        continue;
      }
      for (UF_long k = 0; k < nk; ++k) {
        if (run_of[k1 + k] != run + 1u || (k > 0 && run_of[k1 + k - 1u] == run + 1u)) {
          continue; /* only handle each run at its first column */
        }
        UF_long e = k;
        while (e + 1u < nk && run_of[k1 + e + 1u] == run + 1u) {
          e++;
        }
        UF_long ulen_run = 0;
        for (UF_long c = k; c <= e; ++c) {
          UF_long *li;
          double *lx;
          UF_long l1;
          kls_klu_get_pointer(lu, lip, llen, c, &li, &lx, &l1);
          for (UF_long p = 0; p < l1; ++p) {
            if (mark[li[p]] == 0u) {
              mark[li[p]] = 1u;
              ulist[ulen_run++] = li[p];
            }
          }
        }
        /* sort the union (insertion into sorted output via simple qsort) */
        for (UF_long i = 1; i < ulen_run; ++i) {
          const UF_long v = ulist[i];
          UF_long j2 = i;
          while (j2 > 0 && ulist[j2 - 1u] > v) {
            ulist[j2] = ulist[j2 - 1u];
            j2--;
          }
          ulist[j2] = v;
        }
        for (UF_long u = 0; u < ulen_run; ++u) {
          mark[ulist[u]] = (UF_long)(u + 1u); /* slot+1 while filling */
          solver->padded_union_rows[upos + u] = ulist[u];
        }
        solver->padded_run_start[run] = k1 + k;
        solver->padded_run_len[run] = e - k + 1u;
        solver->padded_union_ptr[run + 1u] = upos + ulen_run;
        solver->padded_panel_ptr[run + 1u] =
          (UF_long)(ppos + (double)(e - k + 1u) * (double)ulen_run);
        for (UF_long c = k; c <= e; ++c) {
          UF_long *li;
          double *lx;
          UF_long l1;
          kls_klu_get_pointer(lu, lip, llen, c, &li, &lx, &l1);
          slot_ptr[k1 + c] = spos;
          for (UF_long p = 0; p < l1; ++p) {
            solver->padded_slots[spos++] = mark[li[p]] - 1u;
          }
        }
        for (UF_long u = 0; u < ulen_run; ++u) {
          mark[ulist[u]] = 0u;
        }
        upos += ulen_run;
        ppos += (double)(e - k + 1u) * (double)ulen_run;
        run++;
      }
    }
    solver->padded_slot_ptr = slot_ptr;
    solver->padded_run_of = run_of;
    free(mark);
    free(ulist);
    if (getenv("KLS_TRACE_SNODE") != NULL) {
      fprintf(stderr,
              "KLS padded panels: z=%ld runs=%ld unions=%ld slots=%ld"
              " panel_doubles=%.3g\n",
              (long)z, (long)nruns, (long)union_total, (long)slot_total,
              ppos);
    }
    return;
  }
fail_arrays:
  free(solver->padded_run_start);
  free(solver->padded_run_len);
  free(solver->padded_union_ptr);
  free(solver->padded_union_rows);
  free(solver->padded_slots);
  free(solver->padded_panel_ptr);
  free(solver->padded_panel_values);
  solver->padded_run_start = NULL;
  solver->padded_run_len = NULL;
  solver->padded_union_ptr = NULL;
  solver->padded_union_rows = NULL;
  solver->padded_slots = NULL;
  solver->padded_panel_ptr = NULL;
  solver->padded_panel_values = NULL;
  solver->padded_run_count = 0;
fail:
  free(run_of);
  free(slot_ptr);
  free(mark);
  free(ulist);
}

static void kls_maybe_prepare_snode_panels(kls_solver *solver,
                                           double *elapsed) {
  if (solver == NULL || elapsed == NULL || solver->snode_prepared ||
      solver->symbolic == NULL || solver->numeric == NULL ||
      solver->numeric->LUbx == NULL || solver->symbolic->R == NULL ||
      solver->n < 4000u || solver->common.flops < 5.0e6 ||
      kls_snode_panel_env_disabled()) {
    if (solver != NULL && getenv("KLS_TRACE_SNODE") != NULL) {
      fprintf(stderr,
              "KLS snode: prep declined A prepared=%d flops=%.3g\n",
              solver->snode_prepared, solver->common.flops);
    }
    return;
  }
  /* Sorting reorders every L/U column's packed storage, so no
     position-based retained structure may exist yet.  KLS-first assembled
     numerics keep their own seeded row metadata, so leave them unsorted. */
  if (solver->refactor_col_ptr != NULL ||
      solver->refactor_level_ptr != NULL ||
      solver->refactor_l_indices32 != NULL ||
      solver->row_refactor_group_count != 0u ||
      solver->row_refactor_values_ready ||
      solver->stats.last_factor_path == KLS_FACTOR_PATH_KLS_FIRST) {
    if (getenv("KLS_TRACE_SNODE") != NULL) {
      fprintf(stderr,
              "KLS snode: prep declined B map=%d lvl=%d l32=%d rowg=%ld "
              "rowv=%d path=%d\n",
              solver->refactor_col_ptr != NULL,
              solver->refactor_level_ptr != NULL,
              solver->refactor_l_indices32 != NULL,
              (long)solver->row_refactor_group_count,
              solver->row_refactor_values_ready,
              (int)solver->stats.last_factor_path);
    }
    return;
  }
  const double start = kls_now_seconds();
  solver->snode_prepared = 1;
  if (!solver->snode_numeric_pre_sorted &&
      !kls_parallel_lu_sort(solver)) {
    *elapsed += kls_now_seconds() - start;
    return;
  }
  UF_long *run_end = (UF_long *)calloc((size_t)solver->n, sizeof(*run_end));
  if (run_end == NULL) {
    *elapsed += kls_now_seconds() - start;
    return;
  }
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  trilinos_klu_l_numeric *numeric = solver->numeric;
  UF_long covered = 0;
  for (UF_long block = 0; block < symbolic->nblocks; ++block) {
    const UF_long k1 = symbolic->R[block];
    const UF_long k2 = symbolic->R[block + 1u];
    const UF_long nk = k2 - k1;
    if (nk < 2u) {
      continue;
    }
    UF_long *lip = numeric->Lip + k1;
    UF_long *llen = numeric->Llen + k1;
    double *lu = (double *)numeric->LUbx[block];
    if (lu == NULL) {
      continue;
    }
    /* Relaxed-run census (KLS_TRACE_RELAXED_SNODE=z): how much of the
       scattered residue becomes batchable if column k+1 may differ
       from column k's tail by up to z entries (padded to the union)?
       Sizes the padded-supernode build's payoff. Measurement only. */
    static long relaxed_z = -2;
    if (relaxed_z == -2) {
      const char *env = getenv("KLS_TRACE_RELAXED_SNODE");
      relaxed_z = env != NULL ? atol(env) : -1;
    }
    if (relaxed_z >= 0) {
      UF_long rel_runs = 0, rel_cols = 0, rel_run_len = 0, pad = 0,
              entries = 0;
      for (UF_long k = 0; k < nk; ++k) {
        entries += llen[k];
        int ext = 0;
        UF_long mism = 0;
        if (k + 1u < nk && llen[k] >= 1) {
          UF_long *li, *li2;
          double *lx, *lx2;
          UF_long l1, l2;
          kls_klu_get_pointer(lu, lip, llen, k, &li, &lx, &l1);
          kls_klu_get_pointer(lu, lip, llen, k + 1u, &li2, &lx2, &l2);
          /* merge walk over li tail (past the k+1 chain entry if
             present) and li2: count symmetric difference */
          UF_long a = li[0] == k + 1u ? 1u : 0u, b = 0;
          while (a < l1 && b < l2) {
            if (li[a] == li2[b]) {
              a++;
              b++;
            } else if (li[a] < li2[b]) {
              a++;
              mism++;
            } else {
              b++;
              mism++;
            }
          }
          mism += (l1 - a) + (l2 - b);
          ext = li[0] == k + 1u && mism <= (UF_long)relaxed_z;
        }
        if (ext) {
          rel_run_len++;
          pad += mism;
        } else {
          if (rel_run_len > 0) {
            rel_runs++;
            rel_cols += rel_run_len + 1u;
          }
          rel_run_len = 0;
        }
      }
      fprintf(stderr,
              "KLS relaxed-snode z=%ld block=%ld nk=%ld: runs=%ld"
              " cols_in=%ld (%.1f%%) pad=%ld entries=%ld\n",
              relaxed_z, (long)block, (long)nk, (long)rel_runs,
              (long)rel_cols,
              nk > 0 ? 100.0 * (double)rel_cols / (double)nk : 0.0,
              (long)pad, (long)entries);
    }
    UF_long start_col = 0;
    for (UF_long k = 0; k < nk; ++k) {
      int extends = 0;
      if (k + 1u < nk && llen[k] >= 1 && llen[k + 1u] == llen[k] - 1) {
        UF_long *li = NULL;
        UF_long *li2 = NULL;
        double *lx = NULL;
        double *lx2 = NULL;
        UF_long l1 = 0;
        UF_long l2 = 0;
        kls_klu_get_pointer(lu, lip, llen, k, &li, &lx, &l1);
        kls_klu_get_pointer(lu, lip, llen, k + 1u, &li2, &lx2, &l2);
        /* Sorted columns make strict supernode nesting a prefix test: the
           next column's whole pattern must be this column's tail. */
        if (li[0] == k + 1u &&
            memcmp(li + 1, li2, (size_t)l2 * sizeof(*li)) == 0) {
          extends = 1;
        }
      }
      if (!extends) {
        if (k > start_col) {
          for (UF_long c = start_col; c <= k; ++c) {
            run_end[k1 + c] = k1 + k + 1u;
          }
          covered += k - start_col + 1u;
        }
        start_col = k + 1u;
      }
    }
  }
  if (covered < 64u) {
    free(run_end);
  } else {
    solver->snode_run_end = run_end;
  }
  kls_build_padded_panels(solver);
  if (getenv("KLS_TRACE_SNODE") != NULL) {
    long hist[7] = {0, 0, 0, 0, 0, 0, 0};
    long nrun = 0;
    if (solver->snode_run_end != NULL) {
      for (UF_long c = 0; c < solver->n;) {
        const UF_long re = run_end[c];
        if (re <= c + 1u) {
          c++;
          continue;
        }
        const UF_long w = re - c;
        int b = w <= 2 ? 0 : w <= 4 ? 1 : w <= 8 ? 2 : w <= 16 ? 3
                : w <= 32 ? 4 : w <= 64 ? 5 : 6;
        hist[b]++;
        nrun++;
        c = re;
      }
    }
    fprintf(stderr,
            "KLS snode: sorted=1 covered_cols=%ld accepted=%d runs=%ld "
            "w2=%ld w3_4=%ld w5_8=%ld w9_16=%ld w17_32=%ld w33_64=%ld "
            "w65+=%ld\n",
            (long)covered, solver->snode_run_end != NULL, nrun, hist[0],
            hist[1], hist[2], hist[3], hist[4], hist[5], hist[6]);
    /* relaxed-run simulation: continue while the next column's pattern
       is a SUBSET of the current union tail (dropped rows = padding) */
    {
      UF_long rcovered = 0;
      int64_t rentries = 0;
      int64_t rpad = 0;
      int64_t rruns = 0;
      int64_t total_entries = 0;
      for (UF_long block = 0; block < symbolic->nblocks; ++block) {
        const UF_long k1 = symbolic->R[block];
        const UF_long k2 = symbolic->R[block + 1u];
        const UF_long nk = k2 - k1;
        if (nk < 2u || numeric->LUbx[block] == NULL) {
          continue;
        }
        UF_long *lip = numeric->Lip + k1;
        UF_long *llen = numeric->Llen + k1;
        double *lu = (double *)numeric->LUbx[block];
        UF_long rstart = 0;
        UF_long *uni = NULL;
        UF_long uni_len = 0;
        for (UF_long k = 0; k < nk; ++k) {
          UF_long *li = NULL;
          double *lx = NULL;
          UF_long ll = 0;
          kls_klu_get_pointer(lu, lip, llen, k, &li, &lx, &ll);
          total_entries += ll;
          int extends = 0;
          if (k + 1u < nk && k > rstart - 1u) {
            /* union tail = uni minus leading diag hits; next col must
               start with k+1 and be a subset */
            UF_long *li2 = NULL;
            double *lx2 = NULL;
            UF_long l2 = 0;
            kls_klu_get_pointer(lu, lip, llen, k + 1u, &li2, &lx2, &l2);
            const UF_long *base = uni != NULL ? uni : li;
            const UF_long blen = uni != NULL ? uni_len : ll;
            /* advance past rows <= k+1 in base; li2 must all be found */
            UF_long bi = 0;
            while (bi < blen && base[bi] <= k + 1u) bi++;
            /* subset test: also require li2[0..] to include nothing new
               and next col to be reached (base contains k+1) */
            int has_diag = 0;
            for (UF_long q = 0; q < blen && base[q] <= k + 1u; ++q) {
              if (base[q] == k + 1u) has_diag = 1;
            }
            if (has_diag && l2 <= blen - bi) {
              UF_long p2 = 0;
              UF_long q = bi;
              while (p2 < l2 && q < blen) {
                if (base[q] == li2[p2]) { p2++; q++; }
                else if (base[q] < li2[p2]) q++;
                else break;
              }
              if (p2 == l2) extends = 1;
            }
            if (extends && uni == NULL) {
              uni = li;
              uni_len = ll;
            }
          }
          if (!extends) {
            if (uni != NULL && k > rstart) {
              const UF_long w = k - rstart + 1u;
              rruns++;
              rcovered += w;
              /* panel entries = w * (uni tail length per col approx) */
              for (UF_long c = rstart; c <= k; ++c) {
                UF_long *lic; double *lxc; UF_long llc;
                kls_klu_get_pointer(lu, lip, llen, c, &lic, &lxc, &llc);
                rentries += llc;
              }
              /* pad = union positions minus real entries */
              UF_long below = 0;
              for (UF_long q = 0; q < uni_len; ++q) {
                if (uni[q] > k) below++;
              }
              rpad += (int64_t)(below + w) * w - ((int64_t)0);
            }
            uni = NULL;
            uni_len = 0;
            rstart = k + 1u;
          }
        }
      }
      fprintf(stderr,
              "KLS snode relaxed: runs=%lld cols=%ld entries=%lld/%lld"
              " (%.1f%%) panel_cells~%lld\n",
              (long long)rruns, (long)rcovered, (long long)rentries,
              (long long)total_entries,
              total_entries > 0
                ? 100.0 * (double)rentries / (double)total_entries : 0.0,
              (long long)rpad);
    }
  }
  *elapsed += kls_now_seconds() - start;
}

/* The deferred-preparation path can overlap run detection with the compact
   solve and column-schedule builders, but all of them must observe the final
   packed-LU order.  Perform only that mutating sort on the caller first;
   kls_maybe_prepare_snode_panels then sees the pre-sorted flag on its worker
   and completes the read-only run census concurrently. */
static int kls_prepare_snode_sort_for_overlap(kls_solver *solver,
                                              double *elapsed) {
  if (solver == NULL || elapsed == NULL || solver->snode_prepared ||
      solver->symbolic == NULL || solver->numeric == NULL ||
      solver->numeric->LUbx == NULL || solver->symbolic->R == NULL ||
      solver->n < 4000u || solver->common.flops < 5.0e6 ||
      kls_snode_panel_env_disabled() || solver->refactor_col_ptr != NULL ||
      solver->refactor_level_ptr != NULL ||
      solver->refactor_l_indices32 != NULL ||
      solver->row_refactor_group_count != 0u ||
      solver->row_refactor_values_ready ||
      solver->stats.last_factor_path == KLS_FACTOR_PATH_KLS_FIRST) {
    return 0;
  }
  if (solver->snode_numeric_pre_sorted) {
    return 1;
  }
  const double start = kls_now_seconds();
  const int sorted = kls_parallel_lu_sort(solver) != 0;
  *elapsed += kls_now_seconds() - start;
  if (sorted) {
    solver->snode_numeric_pre_sorted = 1;
  }
  return sorted;
}

static void maybe_prepare_refactor_map(kls_solver *solver,
                                       double *elapsed) {
  const int pool_map = kls_refactor_pool_map_is_worthwhile(solver);
  if (solver == NULL || elapsed == NULL ||
      kls_low_work_single_block_direct_public_capable(solver) ||
      solver->refactor_col_ptr != NULL ||
      solver->numeric == NULL ||
      (!pool_map && (solver->common.scale > 0 ||
                     solver->numeric->Rs != NULL ||
                     kls_parallel_refactor_is_eligible(solver))) ||
      !kls_refactor_map_is_eligible(solver)) {
    return;
  }
  const double start = kls_now_seconds();
  (void)kls_build_refactor_map(solver);
  *elapsed += kls_now_seconds() - start;
}

static int kls_refactor_schedule_is_eligible(const kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      solver->options.threads <= 1 || solver->n == 0) {
    return 0;
  }
  if (solver->numeric->Uip == NULL || solver->numeric->Ulen == NULL ||
      solver->numeric->LUbx == NULL || solver->symbolic->R == NULL) {
    return 0;
  }
  if (solver->dense_tail_cols > 0) {
    /* pipe-emitted dense-tail numerics: the fragmented-BTF gate below
       misses ss1 by hairs (nblocks 95977 vs the 100000 floor, coverage
       0.53 vs the 0.50 cap, and a scale exclusion the egraph body does
       not need).  The tail columns no-op in the dispatch and the BLAS3
       refresh rebuilds them after the parallel pass. */
    return solver->common.flops >= 5.0e8;
  }
  {
    /* Building the exact dependency graph is worthwhile when the retained
       numeric is large enough to dispatch and contains at least one pool's
       worth of fixed-pivot work.  Its measured schedule makes the final
       admission decision; no benchmark-family signature is needed here. */
    return solver->n >= kls_egraph_refactor_size_floor() &&
           solver->common.flops >=
             KLS_EGRAPH_REFACTOR_MIN_FLOPS_PER_THREAD *
               (double)solver->options.threads;
  }
}

static int kls_build_refactor_pipeline_graph(
  kls_solver *solver,
  const UF_long *levels,
  UF_long cluster_levels,
  UF_long level_count,
  UF_long pipeline_columns,
  UF_long **successor_ptr_out,
  UF_long **successors_out,
  UF_long **pred_count_out,
  UF_long *edge_count_out,
  UF_long **ready_cols_out,
  atomic_uint **ready_slots_out,
  atomic_ulong **remaining_preds_out) {
  KLS_SET_OPTIONAL_OUTPUT(successor_ptr_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(successors_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(pred_count_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(edge_count_out, 0u);
  KLS_SET_OPTIONAL_OUTPUT(ready_cols_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(ready_slots_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(remaining_preds_out, NULL);
  if (solver == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL || levels == NULL ||
      successor_ptr_out == NULL || successors_out == NULL ||
      pred_count_out == NULL || edge_count_out == NULL ||
      ready_cols_out == NULL || ready_slots_out == NULL ||
      remaining_preds_out == NULL) {
    return 0;
  }
  if (pipeline_columns == 0u || cluster_levels >= level_count) {
    return 1;
  }

  UF_long *successor_counts =
    (UF_long *)calloc((size_t)solver->n, sizeof(*successor_counts));
  UF_long *pred_count =
    (UF_long *)calloc((size_t)solver->n, sizeof(*pred_count));
  UF_long *successor_ptr =
    (UF_long *)calloc((size_t)solver->n + 1u, sizeof(*successor_ptr));
  if (successor_counts == NULL || pred_count == NULL ||
      successor_ptr == NULL) {
    free(successor_counts);
    free(pred_count);
    free(successor_ptr);
    return 0;
  }

  UF_long edge_count = 0u;
  for (UF_long block = 0; block < solver->symbolic->nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    const UF_long nk = k2 - k1;
    if (nk <= 1u) {
      continue;
    }
    double *lu = (double *)solver->numeric->LUbx[block];
    if (lu == NULL) {
      free(successor_counts);
      free(pred_count);
      free(successor_ptr);
      return 0;
    }
    const UF_long *uip = solver->numeric->Uip + k1;
    const UF_long *ulen = solver->numeric->Ulen + k1;
    for (UF_long k = 0; k < nk; ++k) {
      const UF_long global_col = k1 + k;
      if (levels[global_col] < cluster_levels) {
        continue;
      }
      UF_long *ui = NULL;
      double *ux = NULL;
      UF_long ucol_len = 0;
      kls_klu_get_pointer(lu, uip, ulen, k, &ui, &ux, &ucol_len);
      (void)ux;
      if (ui == NULL && ucol_len > 0u) {
        free(successor_counts);
        free(pred_count);
        free(successor_ptr);
        return 0;
      }
      for (UF_long p = 0; p < ucol_len; ++p) {
        const UF_long dep = ui[p];
        if (dep >= k) {
          free(successor_counts);
          free(pred_count);
          free(successor_ptr);
          return 0;
        }
        const UF_long global_dep = k1 + dep;
        if (levels[global_dep] < cluster_levels) {
          continue;
        }
        pred_count[global_col]++;
        successor_counts[global_dep]++;
        edge_count++;
      }
    }
  }

  for (UF_long col = 0; col < solver->n; ++col) {
    successor_ptr[col + 1u] = successor_ptr[col] + successor_counts[col];
  }
  if (successor_ptr[solver->n] != edge_count) {
    free(successor_counts);
    free(pred_count);
    free(successor_ptr);
    return 0;
  }

  UF_long *successors =
    (UF_long *)malloc((size_t)(edge_count > 0u ? edge_count : 1u) *
                      sizeof(*successors));
  UF_long *cursor =
    (UF_long *)malloc((size_t)solver->n * sizeof(*cursor));
  UF_long *ready_cols =
    (UF_long *)malloc((size_t)pipeline_columns * sizeof(*ready_cols));
  atomic_uint *ready_slots =
    (atomic_uint *)malloc((size_t)pipeline_columns * sizeof(*ready_slots));
  atomic_ulong *remaining_preds =
    (atomic_ulong *)malloc((size_t)solver->n * sizeof(*remaining_preds));
  if (successors == NULL || cursor == NULL || ready_cols == NULL ||
      ready_slots == NULL || remaining_preds == NULL) {
    free(successor_counts);
    free(pred_count);
    free(successor_ptr);
    free(successors);
    free(cursor);
    free(ready_cols);
    free(ready_slots);
    free(remaining_preds);
    return 0;
  }
  for (UF_long col = 0; col < solver->n; ++col) {
    cursor[col] = successor_ptr[col];
    atomic_init(&remaining_preds[col], 0ul);
  }
  for (UF_long pos = 0; pos < pipeline_columns; ++pos) {
    ready_cols[pos] = 0u;
    atomic_init(&ready_slots[pos], 0u);
  }

  for (UF_long block = 0; block < solver->symbolic->nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    const UF_long nk = k2 - k1;
    if (nk <= 1u) {
      continue;
    }
    double *lu = (double *)solver->numeric->LUbx[block];
    const UF_long *uip = solver->numeric->Uip + k1;
    const UF_long *ulen = solver->numeric->Ulen + k1;
    for (UF_long k = 0; k < nk; ++k) {
      const UF_long global_col = k1 + k;
      if (levels[global_col] < cluster_levels) {
        continue;
      }
      UF_long *ui = NULL;
      double *ux = NULL;
      UF_long ucol_len = 0;
      kls_klu_get_pointer(lu, uip, ulen, k, &ui, &ux, &ucol_len);
      (void)ux;
      for (UF_long p = 0; p < ucol_len; ++p) {
        const UF_long global_dep = k1 + ui[p];
        if (levels[global_dep] < cluster_levels) {
          continue;
        }
        const UF_long dst = cursor[global_dep]++;
        if (dst >= edge_count) {
          free(successor_counts);
          free(pred_count);
          free(successor_ptr);
          free(successors);
          free(cursor);
          free(ready_cols);
          free(ready_slots);
          free(remaining_preds);
          return 0;
        }
        successors[dst] = global_col;
      }
    }
  }
  for (UF_long col = 0; col < solver->n; ++col) {
    if (cursor[col] != successor_ptr[col + 1u]) {
      free(successor_counts);
      free(pred_count);
      free(successor_ptr);
      free(successors);
      free(cursor);
      free(ready_cols);
      free(ready_slots);
      free(remaining_preds);
      return 0;
    }
  }

  free(successor_counts);
  free(cursor);
  *successor_ptr_out = successor_ptr;
  *successors_out = successors;
  *pred_count_out = pred_count;
  *edge_count_out = edge_count;
  *ready_cols_out = ready_cols;
  *ready_slots_out = ready_slots;
  *remaining_preds_out = remaining_preds;
  return 1;
}

typedef struct kls_refactor_separator_component_entry {
  unsigned int component;
  UF_long columns;
  double work;
} kls_refactor_separator_component_entry;

static int kls_compare_refactor_separator_component_work_desc(
  const void *lhs_ptr, const void *rhs_ptr) {
  const kls_refactor_separator_component_entry *lhs =
    (const kls_refactor_separator_component_entry *)lhs_ptr;
  const kls_refactor_separator_component_entry *rhs =
    (const kls_refactor_separator_component_entry *)rhs_ptr;
  if (lhs->work > rhs->work) {
    return -1;
  }
  if (lhs->work < rhs->work) {
    return 1;
  }
  if (lhs->columns > rhs->columns) {
    return -1;
  }
  if (lhs->columns < rhs->columns) {
    return 1;
  }
  return lhs->component < rhs->component
    ? -1 : lhs->component > rhs->component ? 1 : 0;
}

/* Build a barrier-free cluster prefix from the separator forest retained by
   METIS.  A leaf component is admitted only after checking the actual U
   dependency pattern: every predecessor of every admitted cluster column
   must belong to the same component.  Whole components can then be assigned
   to thread-private queues, while the (normally small) rejected/separator
   remainder stays in natural topological order after one join. */
static int kls_build_refactor_separator_private_plan(
  const kls_solver *solver,
  const UF_long *levels,
  UF_long cluster_levels,
  const double *column_work,
  int thread_count,
  UF_long **private_cols_out,
  UF_long **thread_ptr_out,
  UF_long **tail_cols_out,
  UF_long **tail_level_ptr_out,
  UF_long **tail_level_thread_ptr_out,
  UF_long *private_column_count_out,
  UF_long *tail_column_count_out,
  UF_long *private_component_count_out,
  UF_long *unsafe_component_count_out) {
  KLS_SET_OPTIONAL_OUTPUT(private_cols_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(thread_ptr_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(tail_cols_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(tail_level_ptr_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(tail_level_thread_ptr_out, NULL);
  KLS_SET_OPTIONAL_OUTPUT(private_column_count_out, 0u);
  KLS_SET_OPTIONAL_OUTPUT(tail_column_count_out, 0u);
  KLS_SET_OPTIONAL_OUTPUT(private_component_count_out, 0u);
  KLS_SET_OPTIONAL_OUTPUT(unsafe_component_count_out, 0u);
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      levels == NULL || column_work == NULL || thread_count <= 1 ||
      cluster_levels == 0u || private_cols_out == NULL ||
      thread_ptr_out == NULL || tail_cols_out == NULL ||
      tail_level_ptr_out == NULL || tail_level_thread_ptr_out == NULL ||
      private_column_count_out == NULL || tail_column_count_out == NULL ||
      private_component_count_out == NULL ||
      unsafe_component_count_out == NULL ||
      !kls_separator_analysis_has_global_range(&solver->separator) ||
      solver->separator.global_begin != 0u ||
      solver->separator.global_end != solver->n ||
      solver->separator.n != solver->n ||
      solver->separator.component_count == 0u ||
      solver->separator.component_count > (UF_long)UINT_MAX ||
      solver->separator.component_kind == NULL ||
      solver->separator.order_component == NULL) {
    return 0;
  }

  const UF_long component_count = solver->separator.component_count;
  if (component_count > (UF_long)(SIZE_MAX / sizeof(UF_long)) ||
      component_count > (UF_long)(SIZE_MAX / sizeof(double)) ||
      component_count > (UF_long)(SIZE_MAX / sizeof(int)) ||
      component_count >
        (UF_long)(SIZE_MAX / sizeof(kls_refactor_separator_component_entry)) ||
      solver->n > (UF_long)(SIZE_MAX / sizeof(UF_long)) ||
      (size_t)thread_count > SIZE_MAX / sizeof(UF_long) - 1u ||
      (size_t)thread_count > SIZE_MAX / sizeof(double)) {
    return 0;
  }

  unsigned char *unsafe =
    (unsigned char *)calloc((size_t)component_count, sizeof(*unsafe));
  UF_long *component_columns =
    (UF_long *)calloc((size_t)component_count, sizeof(*component_columns));
  double *component_work =
    (double *)calloc((size_t)component_count, sizeof(*component_work));
  int *component_thread =
    (int *)malloc((size_t)component_count * sizeof(*component_thread));
  kls_refactor_separator_component_entry *entries =
    (kls_refactor_separator_component_entry *)malloc(
      (size_t)component_count * sizeof(*entries));
  UF_long *thread_counts =
    (UF_long *)calloc((size_t)thread_count, sizeof(*thread_counts));
  double *thread_work =
    (double *)calloc((size_t)thread_count, sizeof(*thread_work));
  UF_long *thread_ptr = NULL;
  UF_long *private_cols = NULL;
  UF_long *tail_cols = NULL;
  UF_long *tail_level_ptr = NULL;
  UF_long *tail_level_thread_ptr = NULL;
  UF_long *tail_next = NULL;
  int ok = 0;
  if (unsafe == NULL || component_columns == NULL ||
      component_work == NULL || component_thread == NULL ||
      entries == NULL || thread_counts == NULL || thread_work == NULL) {
    goto cleanup;
  }
  for (UF_long component = 0u; component < component_count; ++component) {
    component_thread[component] = -1;
  }

  /* Check independence against the post-pivot numeric pattern rather than
     assuming that the ordering-time graph remained unchanged. */
  for (UF_long block = 0u; block < solver->symbolic->nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    const UF_long nk = k2 - k1;
    if (nk <= 1u) {
      continue;
    }
    double *lu = (double *)solver->numeric->LUbx[block];
    if (lu == NULL) {
      goto cleanup;
    }
    const UF_long *uip = solver->numeric->Uip + k1;
    const UF_long *ulen = solver->numeric->Ulen + k1;
    for (UF_long k = 0u; k < nk; ++k) {
      const UF_long col = k1 + k;
      if (levels[col] >= cluster_levels) {
        continue;
      }
      const unsigned int component =
        solver->separator.order_component[col];
      if ((UF_long)component >= component_count) {
        goto cleanup;
      }
      if (solver->separator.component_kind[component] != 0u) {
        continue;
      }
      UF_long *ui = NULL;
      double *ux = NULL;
      UF_long ucol_len = 0u;
      kls_klu_get_pointer(lu, uip, ulen, k, &ui, &ux, &ucol_len);
      (void)ux;
      if (ui == NULL && ucol_len > 0u) {
        goto cleanup;
      }
      for (UF_long p = 0u; p < ucol_len; ++p) {
        if (ui[p] >= k) {
          goto cleanup;
        }
        const UF_long dep = k1 + ui[p];
        const unsigned int dep_component =
          solver->separator.order_component[dep];
        if ((UF_long)dep_component >= component_count) {
          goto cleanup;
        }
        if (dep_component != component) {
          unsafe[component] = 1u;
          break;
        }
      }
    }
  }

  UF_long cluster_columns = 0u;
  for (UF_long col = 0u; col < solver->n; ++col) {
    if (levels[col] >= cluster_levels) {
      continue;
    }
    cluster_columns++;
    const unsigned int component = solver->separator.order_component[col];
    if ((UF_long)component >= component_count) {
      goto cleanup;
    }
    if (solver->separator.component_kind[component] == 0u &&
        !unsafe[component]) {
      component_columns[component]++;
      component_work[component] += column_work[col];
    }
  }

  UF_long private_component_count = 0u;
  UF_long unsafe_component_count = 0u;
  UF_long private_column_count = 0u;
  for (UF_long component = 0u; component < component_count; ++component) {
    if (unsafe[component]) {
      unsafe_component_count++;
    }
    if (component_columns[component] == 0u) {
      continue;
    }
    if (component > (UF_long)UINT_MAX) {
      goto cleanup;
    }
    entries[private_component_count].component = (unsigned int)component;
    entries[private_component_count].columns = component_columns[component];
    entries[private_component_count].work = component_work[component];
    private_component_count++;
    private_column_count += component_columns[component];
  }
  if (private_component_count == 0u ||
      private_column_count < (UF_long)(2 * thread_count) ||
      private_column_count > cluster_columns) {
    goto cleanup;
  }
  const UF_long tail_column_count = cluster_columns - private_column_count;
  qsort(entries, (size_t)private_component_count, sizeof(*entries),
        kls_compare_refactor_separator_component_work_desc);

  for (UF_long pos = 0u; pos < private_component_count; ++pos) {
    int target = 0;
    for (int tid = 1; tid < thread_count; ++tid) {
      if (thread_work[tid] < thread_work[target]) {
        target = tid;
      }
    }
    const unsigned int component = entries[pos].component;
    component_thread[component] = target;
    thread_counts[target] += entries[pos].columns;
    thread_work[target] += entries[pos].work;
  }

  thread_ptr = (UF_long *)calloc((size_t)thread_count + 1u,
                                 sizeof(*thread_ptr));
  private_cols =
    (UF_long *)malloc((size_t)private_column_count * sizeof(*private_cols));
  tail_cols = tail_column_count > 0u
    ? (UF_long *)malloc((size_t)tail_column_count * sizeof(*tail_cols))
    : NULL;
  if (thread_ptr == NULL || private_cols == NULL ||
      (tail_column_count > 0u && tail_cols == NULL)) {
    goto cleanup;
  }
  for (int tid = 0; tid < thread_count; ++tid) {
    thread_ptr[tid + 1] = thread_ptr[tid] + thread_counts[tid];
    thread_counts[tid] = thread_ptr[tid];
  }
  if (thread_ptr[thread_count] != private_column_count) {
    goto cleanup;
  }

  if (tail_column_count > 0u) {
    tail_level_ptr = (UF_long *)calloc((size_t)cluster_levels + 1u,
                                       sizeof(*tail_level_ptr));
    if (tail_level_ptr == NULL) {
      goto cleanup;
    }
    for (UF_long col = 0u; col < solver->n; ++col) {
      if (levels[col] < cluster_levels) {
        const unsigned int component = solver->separator.order_component[col];
        if ((UF_long)component >= component_count) {
          goto cleanup;
        }
        if (component_thread[component] < 0) {
          tail_level_ptr[levels[col] + 1u]++;
        }
      }
    }
    for (UF_long level = 0u; level < cluster_levels; ++level) {
      tail_level_ptr[level + 1u] += tail_level_ptr[level];
    }
    if (tail_level_ptr[cluster_levels] != tail_column_count) {
      goto cleanup;
    }
    tail_next = (UF_long *)malloc((size_t)cluster_levels * sizeof(*tail_next));
    if (tail_next == NULL) {
      goto cleanup;
    }
    memcpy(tail_next, tail_level_ptr,
           (size_t)cluster_levels * sizeof(*tail_next));
  }
  for (UF_long col = 0u; col < solver->n; ++col) {
    if (levels[col] >= cluster_levels) {
      continue;
    }
    const unsigned int component = solver->separator.order_component[col];
    const int target = (UF_long)component < component_count
      ? component_thread[component] : -1;
    if (target >= 0) {
      const UF_long dst = thread_counts[target]++;
      if (dst >= private_column_count) {
        goto cleanup;
      }
      private_cols[dst] = col;
    } else {
      const UF_long dst = tail_next != NULL ? tail_next[levels[col]] : 0u;
      if (tail_next == NULL || dst >= tail_column_count) {
        free(tail_next);
        tail_next = NULL;
        goto cleanup;
      }
      tail_cols[dst] = col;
      tail_next[levels[col]]++;
    }
  }
  for (UF_long level = 0u; level < cluster_levels; ++level) {
    if (tail_next != NULL && tail_next[level] != tail_level_ptr[level + 1u]) {
      free(tail_next);
      tail_next = NULL;
      goto cleanup;
    }
  }
  free(tail_next);
  tail_next = NULL;
  if (tail_column_count > 0u) {
    if (cluster_levels >
          (UF_long)(SIZE_MAX / (size_t)(thread_count + 1)) ||
        (size_t)cluster_levels * (size_t)(thread_count + 1) >
          SIZE_MAX / sizeof(*tail_level_thread_ptr)) {
      goto cleanup;
    }
    const size_t tail_part_count =
      (size_t)cluster_levels * (size_t)(thread_count + 1);
    tail_level_thread_ptr =
      (UF_long *)malloc(tail_part_count * sizeof(*tail_level_thread_ptr));
    if (tail_level_thread_ptr == NULL) {
      goto cleanup;
    }
    for (UF_long level = 0u; level < cluster_levels; ++level) {
      const UF_long begin = tail_level_ptr[level];
      const UF_long end = tail_level_ptr[level + 1u];
      UF_long *parts = tail_level_thread_ptr +
        level * (UF_long)(thread_count + 1);
      parts[0] = begin;
      parts[thread_count] = end;
      double total_level_work = 0.0;
      for (UF_long pos = begin; pos < end; ++pos) {
        total_level_work += column_work[tail_cols[pos]];
      }
      UF_long pos = begin;
      double prefix_work = 0.0;
      for (int tid = 1; tid < thread_count; ++tid) {
        const double target = total_level_work * (double)tid /
                              (double)thread_count;
        while (pos < end && prefix_work < target) {
          prefix_work += column_work[tail_cols[pos]];
          pos++;
        }
        parts[tid] = pos;
      }
    }
  }
  for (int tid = 0; tid < thread_count; ++tid) {
    if (thread_counts[tid] != thread_ptr[tid + 1]) {
      goto cleanup;
    }
  }

  if (getenv("KLS_TRACE_EGRAPH_SEPARATOR_PRIVATE") != NULL) {
    double min_work = thread_work[0];
    double max_work = thread_work[0];
    for (int tid = 1; tid < thread_count; ++tid) {
      if (thread_work[tid] < min_work) {
        min_work = thread_work[tid];
      }
      if (thread_work[tid] > max_work) {
        max_work = thread_work[tid];
      }
    }
    fprintf(stderr,
            "KLS EGraph separator-private plan: cluster=%" PRIu64
            " private=%" PRIu64 " tail=%" PRIu64
            " components=%" PRIu64 " unsafe=%" PRIu64
            " work_imbalance=%.3f\n",
            (uint64_t)cluster_columns, (uint64_t)private_column_count,
            (uint64_t)tail_column_count, (uint64_t)private_component_count,
            (uint64_t)unsafe_component_count,
            min_work > 0.0 ? max_work / min_work : 0.0);
  }

  *private_cols_out = private_cols;
  *thread_ptr_out = thread_ptr;
  *tail_cols_out = tail_cols;
  *tail_level_ptr_out = tail_level_ptr;
  *tail_level_thread_ptr_out = tail_level_thread_ptr;
  *private_column_count_out = private_column_count;
  *tail_column_count_out = tail_column_count;
  *private_component_count_out = private_component_count;
  *unsafe_component_count_out = unsafe_component_count;
  private_cols = NULL;
  thread_ptr = NULL;
  tail_cols = NULL;
  tail_level_ptr = NULL;
  tail_level_thread_ptr = NULL;
  ok = 1;

cleanup:
  free(tail_next);
  free(unsafe);
  free(component_columns);
  free(component_work);
  free(component_thread);
  free(entries);
  free(thread_counts);
  free(thread_work);
  free(thread_ptr);
  free(private_cols);
  free(tail_cols);
  free(tail_level_ptr);
  free(tail_level_thread_ptr);
  return ok;
}

static int kls_build_refactor_schedule(kls_solver *solver) {
  if (!kls_refactor_schedule_is_eligible(solver)) {
    return 0;
  }
  const int build_separator_private =
    kls_egraph_separator_private_requested(solver);
  const int need_pipeline_ready_graph =
    kls_egraph_ready_queue_env_enabled();

  if (((solver->refactor_level_ptr != NULL &&
      solver->refactor_supernode_pipeline_end != NULL &&
      (!need_pipeline_ready_graph ||
       solver->refactor_pipeline_column_count == 0u ||
       solver->refactor_pipeline_successor_ptr != NULL))) &&
      (!build_separator_private ||
       (solver->refactor_separator_private_plan_attempted &&
        (solver->refactor_separator_private_cols == NULL ||
         (solver->refactor_separator_private_thread_ptr != NULL &&
          solver->refactor_separator_private_thread_count ==
            solver->options.threads &&
          solver->refactor_separator_private_cluster_level_count ==
            solver->refactor_cluster_level_count)))) &&
      solver->refactor_level_cols != NULL) {
    return 1;
  }

  free_refactor_schedule(solver);
  UF_long *levels = NULL;
  UF_long *successor_counts = NULL;
  UF_long *first_successor = NULL;
  uint64_t *row_hash1 = NULL;
  uint64_t *row_hash2 = NULL;

  double *column_work = NULL;
  UF_long *counts = NULL;
  UF_long *level_ptr = NULL;
  UF_long *level_cols = NULL;

  UF_long *level_thread_ptr = NULL;
  UF_long *separator_private_cols = NULL;
  UF_long *separator_private_thread_ptr = NULL;
  UF_long *separator_cluster_tail_cols = NULL;
  UF_long *separator_cluster_tail_level_ptr = NULL;
  UF_long *separator_cluster_tail_level_thread_ptr = NULL;
  UF_long *separator_private_cols_alpha4 = NULL;
  UF_long *separator_private_thread_ptr_alpha4 = NULL;
  UF_long *separator_cluster_tail_cols_alpha4 = NULL;
  UF_long *separator_cluster_tail_level_ptr_alpha4 = NULL;
  UF_long *separator_cluster_tail_level_thread_ptr_alpha4 = NULL;
  UF_long separator_private_column_count = 0u;
  UF_long separator_cluster_tail_column_count = 0u;
  UF_long separator_private_component_count = 0u;
  UF_long separator_private_unsafe_component_count = 0u;
  UF_long separator_private_column_count_alpha4 = 0u;
  UF_long separator_cluster_tail_column_count_alpha4 = 0u;
  int separator_private_plan_attempted = 0;
  UF_long *pipeline_successor_ptr = NULL;
  UF_long *pipeline_successors = NULL;
  UF_long *pipeline_pred_count = NULL;
  UF_long pipeline_edge_count = 0u;
  UF_long *pipeline_ready_cols = NULL;
  atomic_uint *pipeline_ready_slots = NULL;
  atomic_ulong *pipeline_remaining_preds = NULL;
  UF_long *supernode_pipeline_end = NULL;

#define KLS_FREE_REFACTOR_SCHEDULE_TEMP() \
  do { \
    free(levels); \
    free(successor_counts); \
    free(first_successor); \
    free(row_hash1); \
    free(row_hash2); \
    free(column_work); \
    free(counts); \
    free(level_ptr); \
    free(level_cols); \
    free(level_thread_ptr); \
    free(separator_private_cols); \
    free(separator_private_thread_ptr); \
    free(separator_cluster_tail_cols); \
    free(separator_cluster_tail_level_ptr); \
    free(separator_cluster_tail_level_thread_ptr); \
    free(separator_private_cols_alpha4); \
    free(separator_private_thread_ptr_alpha4); \
    free(separator_cluster_tail_cols_alpha4); \
    free(separator_cluster_tail_level_ptr_alpha4); \
    free(separator_cluster_tail_level_thread_ptr_alpha4); \
    free(pipeline_successor_ptr); \
    free(pipeline_successors); \
    free(pipeline_pred_count); \
    free(pipeline_ready_cols); \
    free(pipeline_ready_slots); \
    free(pipeline_remaining_preds); \
    free(supernode_pipeline_end); \
  } while (0)

  levels = (UF_long *)calloc((size_t)solver->n, sizeof(*levels));
  successor_counts =
    (UF_long *)calloc((size_t)solver->n, sizeof(*successor_counts));
  first_successor =
    (UF_long *)malloc((size_t)solver->n * sizeof(*first_successor));
  row_hash1 = (uint64_t *)calloc((size_t)solver->n, sizeof(*row_hash1));
  row_hash2 = (uint64_t *)calloc((size_t)solver->n, sizeof(*row_hash2));

  column_work =
    (double *)calloc((size_t)solver->n, sizeof(*column_work));
  supernode_pipeline_end =
    (UF_long *)calloc((size_t)solver->n, sizeof(*supernode_pipeline_end));
  if (levels == NULL || successor_counts == NULL ||
      first_successor == NULL || row_hash1 == NULL ||
      row_hash2 == NULL || column_work == NULL ||
      supernode_pipeline_end == NULL) {
    KLS_FREE_REFACTOR_SCHEDULE_TEMP();
    return 0;
  }

  for (UF_long k = 0; k < solver->n; ++k) {
    first_successor[k] = KLS_KLU_EMPTY;

  }

  UF_long max_level = 0;
  UF_long edges = 0;
  UF_long root_columns = 0;
  double max_column_work = 0.0;
  double total_work = 0.0;

  for (UF_long block = 0; block < solver->symbolic->nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    const UF_long nk = k2 - k1;
    if (nk <= 1u) {
      if (nk == 1u) {
        /* A singleton still loads its input value, writes its diagonal, and
           participates in the level-zero dispatch.  Giving it zero work
           makes a fragmented level look free and can place hundreds of
           thousands of real updates in one worker's nominally empty slice.
           Account for that irreducible work for every retained BTF, rather
           than enabling the correction through an input-shape policy. */
        column_work[k1] = 1.0;
        total_work += 1.0;
        if (max_column_work < 1.0) {
          max_column_work = 1.0;
        }
        root_columns++;
      }
      continue;
    }
    double *lu = (double *)solver->numeric->LUbx[block];
    if (lu == NULL) {
      KLS_FREE_REFACTOR_SCHEDULE_TEMP();
      return 0;
    }
    const UF_long *uip = solver->numeric->Uip + k1;
    const UF_long *ulen = solver->numeric->Ulen + k1;

    for (UF_long k = 0; k < nk; ++k) {
      UF_long *ui = NULL;
      double *ux = NULL;
      UF_long ucol_len = 0;
      kls_klu_get_pointer(lu, uip, ulen, k, &ui, &ux, &ucol_len);

      UF_long level = 0;
      double work = 1.0;
      for (UF_long p = 0; p < ucol_len; ++p) {
        const UF_long dep = ui[p];
        if (dep >= k) {
          KLS_FREE_REFACTOR_SCHEDULE_TEMP();
          return 0;
        }
        const UF_long global_dep = k1 + dep;
        const UF_long global_col = k1 + k;
        successor_counts[global_dep]++;
        if (first_successor[global_dep] == KLS_KLU_EMPTY) {
          first_successor[global_dep] = global_col;
        }
        row_hash1[global_dep] += kls_supernode_hash1(global_col);
        row_hash2[global_dep] += kls_supernode_hash2(global_col);

        const UF_long dep_level = levels[global_dep] + 1u;
        if (dep_level > level) {
          level = dep_level;
        }
        work += 1.0 + (double)solver->numeric->Llen[k1 + dep];

        edges++;
      }
      if (ucol_len == 0u) {
        root_columns++;
      }
      levels[k1 + k] = level;
      column_work[k1 + k] = work;
      if (work > max_column_work) {
        max_column_work = work;
      }
      total_work += work;
      if (level > max_level) {
        max_level = level;
      }
    }
  }

  UF_long leaf_columns = 0;
  UF_long max_fanout = 0;
  for (UF_long k = 0; k < solver->n; ++k) {
    const UF_long fanout = successor_counts[k];
    if (fanout == 0u) {
      leaf_columns++;
    }
    if (fanout > max_fanout) {
      max_fanout = fanout;
    }
  }

  UF_long level_count = max_level + 1u;
  counts = (UF_long *)calloc((size_t)level_count + 1u, sizeof(*counts));
  level_ptr =
    (UF_long *)malloc(((size_t)level_count + 1u) * sizeof(*level_ptr));
  const int thread_count = solver->options.threads;
  level_cols =
    (UF_long *)malloc((size_t)solver->n * sizeof(*level_cols));

  if (thread_count > 1) {
    level_thread_ptr =
      (UF_long *)malloc((size_t)level_count *
                        ((size_t)thread_count + 1u) *
                        sizeof(*level_thread_ptr));
  }
  if (counts == NULL || level_ptr == NULL ||
      level_cols == NULL) {
    KLS_FREE_REFACTOR_SCHEDULE_TEMP();
    return 0;
  }
  if (thread_count > 1 && level_thread_ptr == NULL) {
    KLS_FREE_REFACTOR_SCHEDULE_TEMP();
    return 0;
  }

  for (UF_long k = 0; k < solver->n; ++k) {
    if (levels[k] >= level_count) {
      KLS_FREE_REFACTOR_SCHEDULE_TEMP();
      return 0;
    }
    counts[levels[k] + 1u]++;
  }
  for (UF_long level = 0; level < level_count; ++level) {
    counts[level + 1u] += counts[level];
    level_ptr[level] = counts[level];
  }
  level_ptr[level_count] = counts[level_count];
  UF_long max_width = 0;
  UF_long cluster_levels = level_count;
  /* Stop globally synchronizing before a level narrows to the point where
     barrier latency rivals its useful parallel work.  A 2.5-worker width
     leaves enough independent columns for the fixed slices while handing the
     narrowing wave to the exact dependency pipeline a few levels earlier.
     This boundary depends only on available parallelism; the retained alpha-3
     and alpha-4 plans below still cover schedules with a different balance. */
  double cluster_width_alpha = 2.5;
  int explicit_cluster_width_alpha = 0;

  {
    const char *env = getenv("KLS_CLUSTER_WIDTH_ALPHA");
    if (env != NULL && env[0] != '\0') {
      const double parsed = atof(env);
      if (parsed > 0.0) {
        cluster_width_alpha = parsed;
        explicit_cluster_width_alpha = 1;
      }
    }
  }
  const double cluster_width_limit =
    cluster_width_alpha * (double)solver->options.threads;
  for (UF_long level = 0; level < level_count; ++level) {
    const UF_long width = level_ptr[level + 1u] - level_ptr[level];
    if (width > max_width) {
      max_width = width;
    }
    if (cluster_levels == level_count &&
        (double)width < cluster_width_limit) {
      cluster_levels = level;
    }
  }
  UF_long cluster_levels_alpha3 = cluster_levels;
  UF_long cluster_levels_alpha4 = cluster_levels;
  if (!explicit_cluster_width_alpha) {
    cluster_levels_alpha3 = level_count;
    const double alpha3_limit = 3.0 * (double)solver->options.threads;
    for (UF_long level = 0u; level < level_count; ++level) {
      const UF_long width = level_ptr[level + 1u] - level_ptr[level];
      if (cluster_levels_alpha3 == level_count &&
          (double)width < alpha3_limit) {
        cluster_levels_alpha3 = level;
      }
    }
    cluster_levels_alpha4 = level_count;
    const double alpha4_limit = 4.0 * (double)solver->options.threads;
    for (UF_long level = 0u; level < level_count; ++level) {
      const UF_long width = level_ptr[level + 1u] - level_ptr[level];
      if (cluster_levels_alpha4 == level_count &&
          (double)width < alpha4_limit) {
        cluster_levels_alpha4 = level;
      }
    }
  }

  if (level_cols != NULL) {
    UF_long *next = counts;

    for (UF_long k = 0; k < solver->n; ++k) {
      const UF_long level = levels[k];
      const UF_long dst = next[level]++;
      if (dst >= solver->n) {
        KLS_FREE_REFACTOR_SCHEDULE_TEMP();
        return 0;
      }
      level_cols[dst] = k;

    }
  }

  if (level_thread_ptr != NULL) {
    for (UF_long level = 0; level < level_count; ++level) {
      const UF_long begin = level_ptr[level];
      const UF_long end = level_ptr[level + 1u];
      UF_long *parts =
        level_thread_ptr + level * (UF_long)(thread_count + 1);
      parts[0] = begin;
      parts[thread_count] = end;
      /* Cluster levels keep barriers, but split contiguous slices by the
       * existing no-pivot column-work estimate instead of column count. */
      double level_work = 0.0;
      for (UF_long pos = begin; pos < end; ++pos) {
        const UF_long col = level_cols[pos];
        if (col >= solver->n) {
          KLS_FREE_REFACTOR_SCHEDULE_TEMP();
          return 0;
        }
        level_work += column_work[col];
      }
      UF_long pos = begin;
      double prefix_work = 0.0;
      for (int t = 1; t < thread_count; ++t) {
        const double target =
          level_work * (double)t / (double)thread_count;
        while (pos < end && prefix_work < target) {
          prefix_work += column_work[level_cols[pos]];
          pos++;
        }
        parts[t] = pos;
      }
    }
  }

  UF_long pipeline_columns = 0;
  double pipeline_work = 0.0;
  double pipeline_max_column_work = 0.0;
  if (cluster_levels < level_count) {
    for (UF_long pos = level_ptr[cluster_levels]; pos < solver->n; ++pos) {
      const UF_long col = level_cols[pos];
      if (col >= solver->n) {
        KLS_FREE_REFACTOR_SCHEDULE_TEMP();
        return 0;
      }
      pipeline_columns++;
      const double work = column_work[col];
      pipeline_work += work;
      if (work > pipeline_max_column_work) {
        pipeline_max_column_work = work;
      }
    }
  }

  if (build_separator_private && cluster_levels > 0u) {
    separator_private_plan_attempted = 1;
    (void)kls_build_refactor_separator_private_plan(
      solver, levels, cluster_levels, column_work, thread_count,
      &separator_private_cols, &separator_private_thread_ptr,
      &separator_cluster_tail_cols, &separator_cluster_tail_level_ptr,
      &separator_cluster_tail_level_thread_ptr,
      &separator_private_column_count,
      &separator_cluster_tail_column_count,
      &separator_private_component_count,
      &separator_private_unsafe_component_count);
    if (kls_repeated_update_workload(&solver->options) &&
        separator_private_cols != NULL &&
        separator_private_thread_ptr != NULL &&
        cluster_levels_alpha4 > 0u &&
        cluster_levels_alpha4 != cluster_levels) {
      UF_long alpha4_component_count = 0u;
      UF_long alpha4_unsafe_component_count = 0u;
      (void)kls_build_refactor_separator_private_plan(
        solver, levels, cluster_levels_alpha4, column_work, thread_count,
        &separator_private_cols_alpha4,
        &separator_private_thread_ptr_alpha4,
        &separator_cluster_tail_cols_alpha4,
        &separator_cluster_tail_level_ptr_alpha4,
        &separator_cluster_tail_level_thread_ptr_alpha4,
        &separator_private_column_count_alpha4,
        &separator_cluster_tail_column_count_alpha4,
        &alpha4_component_count,
        &alpha4_unsafe_component_count);
    }
  }

  if (need_pipeline_ready_graph) {
    if (!kls_build_refactor_pipeline_graph(
          solver, levels, cluster_levels, level_count, pipeline_columns,
          &pipeline_successor_ptr, &pipeline_successors,
          &pipeline_pred_count, &pipeline_edge_count,
          &pipeline_ready_cols, &pipeline_ready_slots,
          &pipeline_remaining_preds)) {
      KLS_FREE_REFACTOR_SCHEDULE_TEMP();
      return 0;
    }
  }

  for (UF_long block = 0; block < solver->symbolic->nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    if (k2 - k1 <= 1u) {
      continue;
    }
    UF_long start = k1;
    UF_long width = 1u;
    for (UF_long k = k1; k + 1u < k2; ++k) {
      const UF_long next = k + 1u;
      const uint64_t next_hash1 = kls_supernode_hash1(next);
      const uint64_t next_hash2 = kls_supernode_hash2(next);
      const int extends =
        first_successor[k] == next &&
        successor_counts[k] > 0u &&
        successor_counts[k] - 1u == successor_counts[next] &&
        row_hash1[k] - next_hash1 == row_hash1[next] &&
        row_hash2[k] - next_hash2 == row_hash2[next];
      if (extends) {
        width++;
      } else {
        kls_record_supernode_candidate(solver, start, width,
                                       successor_counts,
                                       supernode_pipeline_end);

        start = next;
        width = 1u;
      }
    }
    kls_record_supernode_candidate(solver, start, width, successor_counts,
                                   supernode_pipeline_end);

  }

  free(levels);
  free(successor_counts);
  free(first_successor);
  free(row_hash1);
  free(row_hash2);

  free(column_work);
  free(counts);

#undef KLS_FREE_REFACTOR_SCHEDULE_TEMP
  solver->refactor_level_ptr = level_ptr;
  solver->refactor_level_cols = level_cols;

  solver->refactor_level_thread_ptr = level_thread_ptr;
  solver->refactor_separator_private_cols = separator_private_cols;
  solver->refactor_separator_private_thread_ptr =
    separator_private_thread_ptr;
  solver->refactor_separator_cluster_tail_cols =
    separator_cluster_tail_cols;
  solver->refactor_separator_cluster_tail_level_ptr =
    separator_cluster_tail_level_ptr;
  solver->refactor_separator_cluster_tail_level_thread_ptr =
    separator_cluster_tail_level_thread_ptr;
  solver->refactor_separator_private_cols_alpha4 =
    separator_private_cols_alpha4;
  solver->refactor_separator_private_thread_ptr_alpha4 =
    separator_private_thread_ptr_alpha4;
  solver->refactor_separator_cluster_tail_cols_alpha4 =
    separator_cluster_tail_cols_alpha4;
  solver->refactor_separator_cluster_tail_level_ptr_alpha4 =
    separator_cluster_tail_level_ptr_alpha4;
  solver->refactor_separator_cluster_tail_level_thread_ptr_alpha4 =
    separator_cluster_tail_level_thread_ptr_alpha4;
  solver->refactor_pipeline_successor_ptr = pipeline_successor_ptr;
  solver->refactor_pipeline_successors = pipeline_successors;
  solver->refactor_pipeline_pred_count = pipeline_pred_count;
  solver->refactor_pipeline_edge_count = pipeline_edge_count;
  solver->refactor_pipeline_ready_cols = pipeline_ready_cols;
  solver->refactor_pipeline_ready_slots = pipeline_ready_slots;
  solver->refactor_pipeline_remaining_preds = pipeline_remaining_preds;
  solver->refactor_supernode_pipeline_end = supernode_pipeline_end;
  solver->refactor_level_thread_count = thread_count;
  solver->refactor_separator_private_thread_count =
    separator_private_cols != NULL ? thread_count : 0;
  solver->refactor_separator_private_plan_attempted =
    separator_private_plan_attempted;
  solver->refactor_level_count = level_count;
  solver->refactor_cluster_level_count_alpha3 = cluster_levels_alpha3;
  solver->refactor_cluster_level_count_alpha4 = cluster_levels_alpha4;
  solver->refactor_separator_private_cluster_level_count =
    separator_private_cols != NULL ? cluster_levels : 0u;
  solver->refactor_separator_private_column_count =
    separator_private_column_count;
  solver->refactor_separator_cluster_tail_column_count =
    separator_cluster_tail_column_count;
  solver->refactor_separator_private_column_count_alpha4 =
    separator_private_column_count_alpha4;
  solver->refactor_separator_cluster_tail_column_count_alpha4 =
    separator_cluster_tail_column_count_alpha4;
  solver->refactor_separator_private_component_count =
    separator_private_component_count;
  solver->refactor_separator_private_unsafe_component_count =
    separator_private_unsafe_component_count;
  solver->refactor_level_max_width = max_width;
  solver->refactor_dependency_edges = edges;
  solver->refactor_dependency_root_columns = root_columns;
  solver->refactor_dependency_leaf_columns = leaf_columns;
  solver->refactor_dependency_max_fanout = max_fanout;
  solver->refactor_dependency_max_column_work = max_column_work;
  solver->refactor_dependency_pipeline_max_column_work =
    pipeline_max_column_work;
  solver->refactor_cluster_level_count = cluster_levels;
  solver->refactor_pipeline_column_count = pipeline_columns;
  solver->refactor_dependency_work = total_work;
  solver->refactor_pipeline_work = pipeline_work;

  level_ptr = NULL;
  level_cols = NULL;

  level_thread_ptr = NULL;
  separator_private_cols = NULL;
  separator_private_thread_ptr = NULL;
  separator_cluster_tail_cols = NULL;
  separator_cluster_tail_level_ptr = NULL;
  separator_cluster_tail_level_thread_ptr = NULL;
  separator_private_cols_alpha4 = NULL;
  separator_private_thread_ptr_alpha4 = NULL;
  separator_cluster_tail_cols_alpha4 = NULL;
  separator_cluster_tail_level_ptr_alpha4 = NULL;
  separator_cluster_tail_level_thread_ptr_alpha4 = NULL;
  pipeline_successor_ptr = NULL;
  pipeline_successors = NULL;
  pipeline_pred_count = NULL;
  pipeline_ready_cols = NULL;
  pipeline_ready_slots = NULL;
  pipeline_remaining_preds = NULL;

  supernode_pipeline_end = NULL;
  return 1;
}

static void maybe_prepare_refactor_schedule(kls_solver *solver,
                                            double *elapsed) {
  if (solver == NULL || elapsed == NULL) {
    return;
  }
  if (!kls_refactor_schedule_is_eligible(solver)) {
    return;
  }
  const double start = kls_now_seconds();
  if (kls_build_refactor_schedule(solver)) {

  }
  *elapsed += kls_now_seconds() - start;
}

static double kls_egraph_refactor_floor(void) {
  double floor = KLS_EGRAPH_REFACTOR_MIN_FLOPS_PER_THREAD;
  const char *floor_env = getenv("KLS_EGRAPH_REFACTOR_FLOOR");
  if (floor_env != NULL && floor_env[0] != '\0') {
    const double parsed = atof(floor_env);
    if (parsed >= 0.0) {
      floor = parsed;
    }
  }
  return floor;
}

static UF_long kls_parallel_refactor(kls_solver *solver,
                                     double *numeric_values,
                                     int check_pivots) {
  /* Any numeric engine below may replace Udiag values.  An EGraph pass can
     republish the reciprocal mirror from its worker pool; all other paths
     leave the public wrapper to refresh it serially. */
  solver->i32solve_udiag_recip_fresh = 0;
  solver->tiny_singleton_rs_recip_fresh = 0;
  solver->parallel_refine_values_copied = 0;
  solver->compact_amf_two_block_exact_recip_fresh = 0;
  /* A tolerance selected below the caller's requested tolerance carries a
     solve-time residual contract.  Allocate its plain-frame snapshot before
     entering a parallel numeric engine so that the active crew can copy the
     current values while workers finish their refactor slices.  Engines that
     decline or run serially leave the flag clear and the public wrapper falls
     back to its existing serial snapshot. */
  if (solver->solve_refine_values == NULL &&
      solver->options.threads > 1 &&
      solver->stats.selected_pivot_tolerance > 0.0 &&
      solver->stats.selected_pivot_tolerance <
        solver->options.pivot_tolerance &&
      (!kls_uses_structural_initial_pivot_tolerance(solver)) &&
      solver->row_perm == NULL && solver->row_scale == NULL &&
      solver->col_scale == NULL && numeric_values != NULL &&
      solver->nnz >= (UF_long)(16384 * solver->options.threads) &&
      getenv("KLS_DISABLE_PARALLEL_REFINE_COPY") == NULL) {
    solver->solve_refine_values = (double *)malloc(
      (size_t)solver->nnz * sizeof(*solver->solve_refine_values));
  }
  /* Refinement computes residuals against the true input values, which for
     pass-through CSC inputs are not otherwise retained past this call. */
  if ((solver->numeric_needs_refinement || kls_fp32_refactor_wanted(solver) ||
       solver->pivot_nudge_count > 0) &&
      numeric_values != solver->pivot_nudge_values &&
      numeric_values != solver->solve_refine_values) {
    if (solver->solve_refine_values == NULL) {
      solver->solve_refine_values = (double *)malloc(
        (size_t)solver->nnz * sizeof(*solver->solve_refine_values));
    }
    if (solver->solve_refine_values != NULL) {
      memcpy(solver->solve_refine_values, numeric_values,
             (size_t)solver->nnz * sizeof(*numeric_values));
      solver->parallel_refine_values_copied = 1;
    }
  }
  if (solver->pivot_nudge_count > 0 && solver->pivot_nudge_values != NULL &&
      numeric_values != solver->pivot_nudge_values) {
    memcpy(solver->pivot_nudge_values, numeric_values,
           (size_t)solver->nnz * sizeof(*numeric_values));
    for (UF_long i = 0; i < solver->pivot_nudge_count; ++i) {
      solver->pivot_nudge_values[solver->pivot_nudge_pos[i]] +=
        solver->pivot_nudge_sigma[i];
    }
    numeric_values = solver->pivot_nudge_values;
  }
  solver->fast_reject_refresh_state = KLS_FAST_REJECT_REFRESH_UNKNOWN;
  kls_clear_egraph_refactor_last_stats(solver);
  kls_set_last_refactor_path(solver, KLS_REFACTOR_PATH_NONE);
  /* Once EGraph has published the authoritative packed factor, keep values
     in that representation across recurring updates.  Re-entering the
     generic engine tournament would require a full synchronization into
     KLU's interleaved storage before even timing another arm.  Structure and
     pivot policy are unchanged for a values-only refactor; an EGraph failure
     synchronizes the partial state before returning to the repair path. */
  if (!check_pivots && solver->refactor_l_packed_valid &&
      solver->refactor_u_packed_valid) {
    const int packed =
      kls_egraph_mapped_refactor(solver, numeric_values, check_pivots);
    if (packed >= 0) {
      kls_set_last_refactor_path(solver, KLS_REFACTOR_PATH_EGRAPH);
      return (UF_long)packed;
    }
    kls_sync_authoritative_packed_l_values(solver);
    kls_sync_authoritative_packed_u_values(solver);
    solver->refactor_l_packed_valid = 0;
    solver->refactor_u_packed_valid = 0;
  }
  /* Enter the lean walk early only after the first changed numeric has
     established its accuracy contract.  A direct user-position map or a
     live partial-BTF plan is durable evidence that a more specialized
     input-consuming route owns this factor; preserve it without naming its
     matrix. */
  if (!check_pivots && kls_low_work_btf_map32_policy_enabled(solver) &&
      solver->solve_contract_probe != 0 &&
      solver->refactor_input_user_pos32 == NULL &&
      solver->snb == NULL && solver->lean_probe_arm <= 0 &&
      !solver->row_refactor_values_ready &&
      !solver->row_refactor_values_dirty &&
      !kls_row_refactor_env_enabled() &&
      kls_build_refactor_map(solver)) {
    const int lean =
      kls_lean_btf_map32_refactor(solver, numeric_values);
    if (lean >= 0) {
      solver->row_refactor_solve_direct_ready = 0;
      solver->row_refactor_solve_validated = 0;
      kls_set_last_refactor_path(solver, KLS_REFACTOR_PATH_MAPPED);
      return (UF_long)lean;
    }
  }
  const int auto_row_refactor =
    kls_auto_row_refactor_should_run(solver);
  if (!check_pivots && solver->lean_probe_arm > 0) {
    solver->lean_pair_active = solver->lean_probe_arm == 2;
    const int lean = kls_lean_row_refactor_numeric(solver, numeric_values);
    if (lean >= 0) {
      kls_set_last_refactor_path(solver, KLS_REFACTOR_PATH_ROW);
      return (UF_long)lean;
    }
    return 0;
  }
  if (check_pivots && kls_checked_row_refactor_should_run()) {
    if (solver->options.threads > 1) {
      const int parallel_row_status =
        kls_threaded_row_refactor_numeric(solver, numeric_values, 1);
      if (parallel_row_status >= 0) {
        kls_set_last_refactor_path(solver, KLS_REFACTOR_PATH_ROW);
        return (UF_long)parallel_row_status;
      }
    }
    const int row_status =
      kls_serial_row_refactor_numeric(solver, numeric_values, 1);
    if (row_status >= 0) {
      kls_set_last_refactor_path(solver, KLS_REFACTOR_PATH_ROW);
      return (UF_long)row_status;
    }
  } else if (!check_pivots) {
    {
      UF_long snb_ok = 0;
      if (kls_snb_try_refactor(solver, numeric_values, 0, &snb_ok)) {
        return snb_ok;
      }
    }
    const int row_acceptance_ready =
      kls_row_refactor_acceptance_structurally_ready(solver);
    const int row_acceptance_wants_row =
      kls_row_refactor_acceptance_wants_row(solver);
    if (row_acceptance_wants_row ||
        (auto_row_refactor && !row_acceptance_ready)) {
    if (solver->options.threads > 1) {
      const int parallel_row_status =
        kls_threaded_row_refactor_numeric(solver, numeric_values, 0);
      if (parallel_row_status >= 0) {
        kls_set_last_refactor_path(solver, KLS_REFACTOR_PATH_ROW);
        return (UF_long)parallel_row_status;
      }
      /* structural decline: never re-prep (~1.4s on ASIC_680ks's
         583K-group schedule) just to decline again on this numeric */
      solver->row_accept_decision = -1;
    }
    if (kls_row_refactor_env_enabled() || auto_row_refactor) {
      const int row_status =
        kls_serial_row_refactor_numeric(solver, numeric_values, 0);
      if (row_status >= 0) {
        kls_set_last_refactor_path(solver, KLS_REFACTOR_PATH_ROW);
        return (UF_long)row_status;
      }
    } else if (solver->row_accept_decision == 0) {
      /* threaded trial infeasible: settle on the column engine */
      solver->row_accept_decision = -1;
    }
    }
  }
  if (!kls_publish_row_refactor_values(solver)) {
    return 0;
  }
  solver->row_refactor_values_ready = 0;
  solver->row_refactor_solve_direct_ready = 0;
  solver->row_refactor_solve_validated = 0;
  if (!check_pivots &&
      (solver->lean_probe_arm > 0 ||
       getenv("KLS_ENABLE_LEAN_ROW_CONSUME") != NULL)) {
    /* A positive arm is an already-settled engine decision.  Dispatch it
       independently of the generic EGraph work floor: fused supernodes can
       make the lean walk the faster high-flop engine as well (rajat03), and
       forcing such a choice through the environment should mean what it
       says. */
    solver->lean_pair_active = solver->lean_probe_arm == 2 ||
      getenv("KLS_LEAN_ROW_PAIR") != NULL;
    const int lean = kls_lean_row_refactor_numeric(solver, numeric_values);
    if (lean >= 0) {
      kls_set_last_refactor_path(solver, KLS_REFACTOR_PATH_ROW);
      return (UF_long)lean;
    }
  }
  const double egraph_floor = kls_egraph_refactor_floor();
  if (!check_pivots && solver->common.flops > 0.0 &&
      solver->common.flops <
        egraph_floor * (double)solver->options.threads) {
    /* Below the per-thread flop floor the pool dispatch costs more than
       it saves; the busy-wait pool moved the measured four-thread
       crossover down to ~1e6 total flops (re-measured 2026-07).  Prefer
       the serial mapped refactor under it and keep the EGraph for
       shapes it cannot cover. */
    if (kls_low_work_single_block_direct_public_capable(solver)) {
      UF_long serial_ok = 0u;
      int compact_used = 0;
      if (!kls_try_compact_map32_klu_tournament(
            solver, numeric_values, &serial_ok, &compact_used)) {
        serial_ok =
          kls_serial_klu_refactor_dense_tail(solver, numeric_values, 0);
      }
      kls_set_last_refactor_path(
        solver, compact_used
          ? KLS_REFACTOR_PATH_MAPPED : KLS_REFACTOR_PATH_KLU);
      return serial_ok;
    }
    {
      UF_long snb_ok = 0;
      if (kls_snb_try_refactor(solver, numeric_values, check_pivots,
                               &snb_ok)) {
        return snb_ok;
      }
    }
    if (!check_pivots && solver->options.threads > 1) {
      const int ptsref = kls_pts_mapped_refactor(solver, numeric_values);
      if (ptsref >= 0) {
        kls_set_last_refactor_path(solver, KLS_REFACTOR_PATH_MAPPED);
        return (UF_long)ptsref;
      }
    }
    const int mapped = kls_mapped_refactor(solver, numeric_values,
                                           check_pivots);
    if (mapped >= 0) {
      kls_set_last_refactor_path(solver, KLS_REFACTOR_PATH_MAPPED);
      return (UF_long)mapped;
    }
    UF_long serial_ok = 0u;
    int compact_used = 0;
    if (!kls_try_compact_map32_klu_tournament(
          solver, numeric_values, &serial_ok, &compact_used)) {
      serial_ok = kls_serial_klu_refactor_dense_tail(
        solver, numeric_values, check_pivots);
    }
    kls_set_last_refactor_path(
      solver, compact_used
        ? KLS_REFACTOR_PATH_MAPPED : KLS_REFACTOR_PATH_KLU);
    if (!serial_ok) {
      kls_record_fast_factor_failure(
        solver, KLS_FAST_FACTOR_FAIL_KLU_REFACTOR_FAILED,
        solver->common.status);
    }
    return serial_ok;
  }
  if (!check_pivots) {
    UF_long snb_ok = 0;
    if (kls_snb_try_refactor(solver, numeric_values, check_pivots,
                             &snb_ok)) {
      return snb_ok;
    }
    {
      UF_long pts_ok = 0;
      if (kls_pts_try_refactor_timed(solver, numeric_values, &pts_ok)) {
        kls_set_last_refactor_path(solver, KLS_REFACTOR_PATH_MAPPED);
        return pts_ok;
      }
    }
  }
  pthread_t deferred_pts_thread;
  const int deferred_pts_active =
    kls_pts_start_deferred_build(solver, &deferred_pts_thread);
  int egraph =
    kls_egraph_mapped_refactor(solver, numeric_values, check_pivots);
  if (deferred_pts_active) {
    pthread_join(deferred_pts_thread, NULL);
  }
  if (egraph > 0 && solver->fp32_last_used && !solver->fp32_validated) {
    /* First reduced-precision pass on this numeric: ill conditioning can
       amplify the storage rounding beyond what refinement recovers, so
       prove the factorization on A e before trusting the mode; on failure
       redo this pass in double precision and lock the mode off. */
    int probe_ok = 0;
    /* The mode's contract is float factors plus one refinement shot per
       solve, so probe what it delivers: solve A x = A e, apply one
       correction, and judge the corrected residual.  Judging the raw
       solve rejects float-cancellation columns the correction absorbs
       (ASIC_320k: six of 6.9M L entries off at 1e-2 raw, corrected
       solves at 2e-15). */
    double *probe =
      (double *)malloc(3u * (size_t)solver->n * sizeof(*probe));
    if (probe != NULL) {
      double *px = probe + solver->n;
      double *pc = px + solver->n;
      for (UF_long i = 0; i < solver->n; ++i) {
        probe[i] = 0.0;
      }
      for (UF_long j = 0; j < solver->n; ++j) {
        for (UF_long p = solver->col_ptr[j]; p < solver->col_ptr[j + 1u];
             ++p) {
          probe[solver->row_idx[p]] += numeric_values[p];
        }
      }
      double bmax = 0.0;
      for (UF_long i = 0; i < solver->n; ++i) {
        const double v = fabs(probe[i]);
        bmax = bmax < v ? v : bmax;
      }
      memcpy(px, probe, (size_t)solver->n * sizeof(*px));
      if (trilinos_klu_l_solve(solver->symbolic, solver->numeric, solver->n,
                               1, px, &solver->common) &&
          solver->common.status >= TRILINOS_KLU_OK) {
        int step_ok = 1;
        double rmax = 0.0;
        for (int step = 0; step < 2 && step_ok; ++step) {
          memcpy(pc, probe, (size_t)solver->n * sizeof(*pc));
          for (UF_long j = 0; j < solver->n; ++j) {
            const double xv = px[j];
            if (xv == 0.0) {
              continue;
            }
            for (UF_long p = solver->col_ptr[j];
                 p < solver->col_ptr[j + 1u]; ++p) {
              pc[solver->row_idx[p]] -= numeric_values[p] * xv;
            }
          }
          rmax = 0.0;
          for (UF_long i = 0; i < solver->n; ++i) {
            const double v = fabs(pc[i]);
            rmax = rmax < v ? v : rmax;
          }
          if (getenv("KLS_TRACE_REFINE") != NULL) {
            fprintf(stderr,
                    "KLS fp32 probe: step %d rmax %.3e bmax %.3e\n",
                    step, rmax, bmax);
          }
          if (step == 0) {
            step_ok =
              trilinos_klu_l_solve(solver->symbolic, solver->numeric,
                                   solver->n, 1, pc, &solver->common) &&
              solver->common.status >= TRILINOS_KLU_OK;
            if (step_ok) {
              for (UF_long i = 0; i < solver->n; ++i) {
                px[i] += pc[i];
              }
            }
          }
        }
        /* Bar chosen against the fleet: healthy fp32 adopters land
           their corrected residual at 1e-14..1e-9 relative (memchip
           7.7e-12/5e2, Freescale1 1.4e-8/8.1); the ASIC_320k failure
           mode sits at 3e-3.  1e-8 keeps the adopters and rejects the
           unrecoverable state by five orders. */
        probe_ok =
          step_ok && rmax <= (bmax > 0.0 ? bmax : 1.0) * 1.0e-8;
      }
      free(probe);
    }
    if (getenv("KLS_TRACE_REFINE") != NULL) {
      fprintf(stderr, "KLS fp32 probe: ok=%d\n", probe_ok);
    }
    if (probe_ok) {
      solver->fp32_validated = 1;
      solver->fp32_decision = 1;
    } else {
      /* The first float pass over a factor-built numeric fails this
         probe deterministically (ASIC_320k: rmax 3e-2) while the next
         pass validates at float accuracy; the failed pass is replaced
         by the double redo below, so give the mode one retry per
         numeric before locking it off. */
      if (solver->fp32_probe_retry_used) {
        solver->fp32_decision = -1;
      }
      solver->fp32_probe_retry_used = 1;
      solver->numeric_needs_refinement = 0;
      solver->fp32_redo_pass = 1;
      egraph = kls_egraph_mapped_refactor(solver, numeric_values,
                                          check_pivots);
      solver->fp32_redo_pass = 0;
    }
  }
  if (egraph >= 0) {
#ifdef KLS_HAVE_CBLAS
    if (egraph > 0 && solver->dense_tail_cols > 0 &&
        solver->common.status >= 0 &&
        getenv("KLS_DISABLE_DENSE_TAIL_REFACTOR") == NULL) {
      /* the egraph pass no-opped the dense-tail columns; rebuild them
         from the refreshed prefix (fallback: one full scalar walk) */
      if (!kls_dense_tail_refresh(solver, numeric_values)) {
        egraph = (int)trilinos_klu_l_refactor(solver->col_ptr,
                                              solver->row_idx,
                                              numeric_values,
                                              solver->symbolic,
                                              solver->numeric,
                                              &solver->common);
      }
    }
#endif
    kls_set_last_refactor_path(solver, KLS_REFACTOR_PATH_EGRAPH);
    if (egraph && check_pivots &&
        !kls_checked_refactor_accepts_rowwise_u(solver)) {
      return 0;
    }
    return (UF_long)egraph;
  }

  if (!kls_parallel_refactor_is_eligible(solver)) {
    {
      UF_long snb_ok = 0;
      if (kls_snb_try_refactor(solver, numeric_values, check_pivots,
                               &snb_ok)) {
        return snb_ok;
      }
    }
    const int mapped = kls_mapped_refactor(solver, numeric_values, check_pivots);
    if (mapped >= 0) {
      kls_set_last_refactor_path(solver, KLS_REFACTOR_PATH_MAPPED);
      if (mapped && check_pivots &&
          !solver->row_refactor_last_run &&
          !kls_checked_refactor_accepts_rowwise_u(solver)) {
        return 0;
      }
      return (UF_long)mapped;
    }
    if (check_pivots && solver->common.scale > 0) {
      const int checked = kls_serial_checked_scaled_refactor(solver,
                                                            numeric_values);
      if (checked >= 0) {
        if (checked &&
            !kls_checked_refactor_accepts_rowwise_u(solver)) {
          return 0;
        }
        return (UF_long)checked;
      }
    }
    UF_long ok = 0u;
    int compact_used = 0;
    if (!kls_try_compact_map32_klu_tournament(
          solver, numeric_values, &ok, &compact_used)) {
      ok = kls_serial_klu_refactor_dense_tail(
        solver, numeric_values, check_pivots);
    }
    kls_set_last_refactor_path(
      solver, compact_used
        ? KLS_REFACTOR_PATH_MAPPED : KLS_REFACTOR_PATH_KLU);
    if (!ok) {
      kls_record_fast_factor_failure(
        solver, KLS_FAST_FACTOR_FAIL_KLU_REFACTOR_FAILED,
        solver->common.status);
    }
    if (ok && check_pivots && solver->common.status >= 0 &&
        solver->common.status != TRILINOS_KLU_SINGULAR) {
      UF_long rejected_pivot = KLS_KLU_EMPTY;
      UF_long rejected_pivot_col = KLS_KLU_EMPTY;
      if (!kls_numeric_pivots_pass_threshold(solver, &rejected_pivot,
                                             &rejected_pivot_col)) {
        solver->fast_reject_refresh_state = KLS_FAST_REJECT_REFRESH_ALL;
        kls_record_fast_reject(solver, rejected_pivot, rejected_pivot_col);
        return 0;
      }
      if (!kls_checked_refactor_accepts_rowwise_u(solver)) {
        return 0;
      }
    }
    return ok;
  }

  trilinos_klu_l_common *common = &solver->common;
  trilinos_klu_l_symbolic *symbolic = solver->symbolic;

  common->status = TRILINOS_KLU_OK;
  common->numerical_rank = KLS_KLU_EMPTY;
  common->singular_col = KLS_KLU_EMPTY;
  common->nrealloc = 0;
  if (common->scale > 0 &&
      !trilinos_klu_l_scale((UF_long)common->scale, solver->n, solver->col_ptr,
                            solver->row_idx, numeric_values, solver->numeric->Rs,
                            NULL, common)) {
    return 0;
  }

  int thread_count = solver->options.threads;
  if ((UF_long)thread_count > symbolic->nblocks) {
    thread_count = (int)symbolic->nblocks;
  }
  if (thread_count < 2) {
    {
      UF_long snb_ok = 0;
      if (kls_snb_try_refactor(solver, numeric_values, check_pivots,
                               &snb_ok)) {
        return snb_ok;
      }
    }
    UF_long ok = 0u;
    int compact_used = 0;
    if (!kls_try_compact_map32_klu_tournament(
          solver, numeric_values, &ok, &compact_used)) {
      ok = kls_serial_klu_refactor_dense_tail(
        solver, numeric_values, check_pivots);
    }
    kls_set_last_refactor_path(
      solver, compact_used
        ? KLS_REFACTOR_PATH_MAPPED : KLS_REFACTOR_PATH_KLU);
    if (!ok) {
      kls_record_fast_factor_failure(
        solver, KLS_FAST_FACTOR_FAIL_KLU_REFACTOR_FAILED,
        solver->common.status);
    }
    if (ok && check_pivots && solver->common.status >= 0 &&
        solver->common.status != TRILINOS_KLU_SINGULAR) {
      UF_long rejected_pivot = KLS_KLU_EMPTY;
      UF_long rejected_pivot_col = KLS_KLU_EMPTY;
      if (!kls_numeric_pivots_pass_threshold(solver, &rejected_pivot,
                                             &rejected_pivot_col)) {
        solver->fast_reject_refresh_state = KLS_FAST_REJECT_REFRESH_ALL;
        kls_record_fast_reject(solver, rejected_pivot, rejected_pivot_col);
        return 0;
      }
      if (!kls_checked_refactor_accepts_rowwise_u(solver)) {
        return 0;
      }
    }
    return ok;
  }

  if (kls_refactor_pool_map_is_worthwhile(solver) &&
      solver->refactor_col_ptr == NULL) {
    (void)kls_build_refactor_map(solver);
  }
  if ((kls_refactor_l_index32_env_enabled() ||
       kls_refactor_u_index32_env_enabled()) &&
      solver->n <= (UF_long)INT32_MAX) {
    (void)kls_build_refactor_lu_pointer_cache(solver);
  }

  int invalid = 0;
  int pivot_rejected = 0;
  UF_long rejected_pivot = KLS_KLU_EMPTY;
  UF_long rejected_pivot_col = KLS_KLU_EMPTY;
  UF_long rejected_row = KLS_KLU_EMPTY;
  double rejected_multiplier_abs = -1.0;
  double rejected_pivot_abs = -1.0;
  double rejected_candidate_abs = -1.0;
  int singular = 0;
  UF_long numerical_rank = UF_long_max;
  UF_long singular_col = KLS_KLU_EMPTY;
  int pool_prefix_current = 0;
  if (!run_refactor_pool(solver, numeric_values, thread_count, 0u,
                         check_pivots,
                         &invalid, &pivot_rejected,
                         &rejected_pivot, &rejected_pivot_col,
                         &rejected_row, &rejected_multiplier_abs,
                         &rejected_pivot_abs, &rejected_candidate_abs,
                         &singular,
                         &numerical_rank, &singular_col,
                         &pool_prefix_current)) {
    kls_set_last_refactor_path(solver, KLS_REFACTOR_PATH_POOL);
    kls_record_fast_factor_failure(
      solver, KLS_FAST_FACTOR_FAIL_POOL_INVALID,
      TRILINOS_KLU_OUT_OF_MEMORY);
    common->status = TRILINOS_KLU_OUT_OF_MEMORY;
    return 0;
  }
  kls_set_last_refactor_path(solver, KLS_REFACTOR_PATH_POOL);

  if (invalid) {
    kls_record_fast_factor_failure(
      solver, KLS_FAST_FACTOR_FAIL_POOL_INVALID, TRILINOS_KLU_INVALID);
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }
  if (pivot_rejected) {
    solver->fast_reject_refresh_state =
      pool_prefix_current ? KLS_FAST_REJECT_REFRESH_PREFIX
                          : KLS_FAST_REJECT_REFRESH_UNKNOWN;
    if (pool_prefix_current) {
      kls_record_fast_reject_block_suffix_tail_seed(solver, rejected_pivot);
    } else {
      kls_clear_fast_reject_tail_seed(solver);
    }
    kls_record_fast_reject_detail(solver, rejected_pivot,
                                  rejected_pivot_col, rejected_row,
                                  rejected_multiplier_abs,
                                  rejected_pivot_abs,
                                  rejected_candidate_abs);
    common->status = TRILINOS_KLU_OK;
    return 0;
  }

  if (singular) {
    common->status = TRILINOS_KLU_SINGULAR;
    common->numerical_rank = numerical_rank;
    common->singular_col = singular_col;
    if (common->halt_if_singular) {
      return 0;
    }
  }

  if (common->scale > 0 && !kls_parallel_refactor_permute_scale(solver)) {
    common->status = TRILINOS_KLU_INVALID;
    return 0;
  }

  if (!singular) {
    common->status = TRILINOS_KLU_OK;
  }
  if (check_pivots && !singular &&
      !kls_checked_refactor_accepts_rowwise_u(solver)) {
    return 0;
  }
  return 1;
}

/* A rejected checked refactor has already measured the fixed-pattern engine
   on this exact retained representation.  Before launching a dynamic-pivot
   block repair, compare its remaining stored-factor work with the measured
   cost of rebuilding the same representation through KLU.  The repair's
   lower-bound model has four factor-work passes: dependency discovery, row
   materialization, dynamic pivot selection, and repacking.  This replaces
   dimension/block-count exception lists with an engine/work comparison. */
static int kls_measured_full_factor_beats_repair(
  const kls_solver *solver,
  UF_long rejected_pivot,
  double checked_seconds) {
  if (solver == NULL ||
      solver->symbolic == NULL || solver->numeric == NULL ||
      !(solver->numeric_full_factor_seconds > 0.0) ||
      !(checked_seconds > 0.0) || rejected_pivot == KLS_KLU_EMPTY) {
    return 0;
  }
  const UF_long block = kls_block_for_pivot(solver, rejected_pivot);
  if (block == KLS_KLU_EMPTY || block >= solver->symbolic->nblocks) {
    return 0;
  }
  const UF_long k1 = solver->symbolic->R[block];
  const UF_long nk = solver->symbolic->R[block + 1u] - k1;
  if (rejected_pivot < k1 || rejected_pivot >= k1 + nk) {
    return 0;
  }

  double repair_work = solver->stats.fast_rejected_pivoting_tail_work;
  if (!(repair_work > 0.0)) {
    UF_long repair_columns = 0u;
    if (!kls_fast_reject_suffix_work(
          solver, block, k1, nk, rejected_pivot - k1,
          &repair_columns, &repair_work)) {
      repair_work = 0.0;
    }
  }
  double total_work = solver->refactor_dependency_work;
  if (!(total_work > 0.0)) {
    double block_work = 0.0;
    if (solver->symbolic->nblocks == 1u &&
        kls_fast_reject_block_work(solver, block, k1, nk, &block_work)) {
      total_work = block_work;
    }
  }
  if (!(repair_work > 0.0) || !(total_work > 0.0)) {
    return 0;
  }
  double repair_fraction = repair_work / total_work;
  if (repair_fraction > 1.0) {
    repair_fraction = 1.0;
  }
  const double projected_repair_seconds =
    4.0 * checked_seconds * repair_fraction;
  const int rebuild =
    projected_repair_seconds >= 0.90 * solver->numeric_full_factor_seconds;
  if (getenv("KLS_TRACE_FAST_FACTOR_BUDGET") != NULL) {
    fprintf(stderr,
            "KLS fast-factor budget: checked=%.6f repair_work=%.3e"
            " total_work=%.3e projected=%.6f rebuild=%.6f choose=%s\n",
            checked_seconds, repair_work, total_work,
            projected_repair_seconds, solver->numeric_full_factor_seconds,
            rebuild ? "rebuild" : "repair");
  }
  return rebuild;
}

static UF_long kls_fast_factor_with_block_restarts(kls_solver *solver,
                                                   double *numeric_values) {
  const int max_restarts = 4;
  for (int attempt = 0; attempt <= max_restarts; ++attempt) {
    const double checked_start = kls_now_seconds();
    const UF_long ok = kls_parallel_refactor(solver, numeric_values, 1);
    const double checked_seconds = kls_now_seconds() - checked_start;
    if (ok || solver->common.status < 0 ||
        solver->common.status == TRILINOS_KLU_SINGULAR) {
      if (!ok && solver->common.status < 0) {
        kls_record_fast_factor_failure(
          solver, KLS_FAST_FACTOR_FAIL_INVALID_STATUS,
          solver->common.status);
      } else if (!ok && solver->common.status == TRILINOS_KLU_SINGULAR) {
        kls_record_fast_factor_failure(
          solver, KLS_FAST_FACTOR_FAIL_SINGULAR_STATUS,
          solver->common.status);
      }
      return ok;
    }
    if (solver->stats.fast_rejected_pivot < 0 ||
        solver->numeric == NULL) {
      kls_record_fast_factor_failure(
        solver, KLS_FAST_FACTOR_FAIL_NO_REJECT, solver->common.status);
      return 0;
    }
    if (attempt == max_restarts) {
      kls_record_fast_factor_failure(
        solver, KLS_FAST_FACTOR_FAIL_NO_REJECT, solver->common.status);
      return 0;
    }
    const UF_long rejected_pivot = (UF_long)solver->stats.fast_rejected_pivot;
    const kls_fast_reject_repair_mode repair_mode =
      kls_fast_reject_repair_mode_from_env();
    if (repair_mode == KLS_FAST_REJECT_REPAIR_FULL ||
        (repair_mode == KLS_FAST_REJECT_REPAIR_DEFAULT &&
         kls_measured_full_factor_beats_repair(
           solver, rejected_pivot, checked_seconds))) {
      solver->full_factor_preferred = 1;
      return 0;
    }
    const int repair_covers_remaining =
      kls_fast_repair_covers_remaining_columns(solver, rejected_pivot);
    if (solver->common.scale > 0) {
      if (!repair_covers_remaining &&
           solver->fast_reject_refresh_state !=
             KLS_FAST_REJECT_REFRESH_PREFIX) {
        return 0;
      }
    }
    const int kls_rebuild_restarts_before =
      solver->fast_kls_rebuild_restarts;
    if (!kls_pivot_restart_rejected_block(
          solver, numeric_values, rejected_pivot, repair_mode)) {
      return 0;
    }
    /* The restart re-pivoted and reallocated block storage (the
       row-first rebuild path frees the old LU and packs a fresh one
       with its own dynamic column pivots): the COLUMN-engine caches
       built for the old layout — LU pointer cache, refactor map,
       dependency schedule, snode panels, dense-tail bookkeeping —
       dangle into freed memory, and the egraph consuming them
       measured raj1 at e-4 residuals with values bit-exact to the
       WRONG factorization (2026-07-15).  Free exactly those; the
       retry below runs schedule-less (serial/row-checked) and the
       deferred consult re-preps on the next refactorization.  The
       ROW structures deliberately survive: they are the repair
       vehicle (the tail refresh below runs the checked row engine on
       them — the dense checked-row smoke encodes that contract) and
       re-derive their pointers from the live numeric per run.  The
       full replacement hammer is NOT safe here for the same reason.
       Contract state is reset: the rebuilt pivot sequence's accuracy
       is unclassified, and any captured refine values describe the
       dead numeric. */
    free_refactor_lu_pointer_cache(solver);
    free_refactor_map(solver);
    free_refactor_schedule(solver);
    free_snode_panels(solver);
    solver->dense_tail_cols = 0;
    solver->dense_tail_block = 0;
    solver->solve_contract_probe = 0;
    solver->solve_contract_verified = 0;
    free(solver->solve_refine_values);
    solver->solve_refine_values = NULL;
    if (solver->n >= 512 && getenv("KLS_SYNC_FACTOR_PREPS") == NULL) {
      solver->factor_preps_deferred = 1;
    }
    if (solver->fast_kls_rebuild_restarts >
        kls_rebuild_restarts_before) {
      UF_long next_rejected_pivot = KLS_KLU_EMPTY;
      UF_long rejected_pivot_col = KLS_KLU_EMPTY;
      if (!kls_numeric_pivots_pass_threshold(
            solver, &next_rejected_pivot, &rejected_pivot_col)) {
        kls_record_fast_reject(solver, next_rejected_pivot,
                               rejected_pivot_col);
        return 0;
      }
      if (!kls_checked_refactor_accepts_rowwise_u(solver)) {
        return 0;
      }
      return 1;
    }
    const UF_long repaired_block = kls_block_for_pivot(solver, rejected_pivot);
    if (repaired_block == KLS_KLU_EMPTY) {
      return 0;
    }
    if (repair_covers_remaining) {
      UF_long next_rejected_pivot = KLS_KLU_EMPTY;
      UF_long rejected_pivot_col = KLS_KLU_EMPTY;
      if (!kls_numeric_pivots_pass_threshold_from_block(
            solver, repaired_block, &next_rejected_pivot,
            &rejected_pivot_col)) {
        kls_record_fast_reject(solver, next_rejected_pivot,
                               rejected_pivot_col);
        return 0;
      }
      if (!kls_checked_refactor_accepts_rowwise_u_from_block(
            solver, repaired_block)) {
        return 0;
      }
      return 1;
    }
    if (solver->fast_reject_refresh_state == KLS_FAST_REJECT_REFRESH_PREFIX) {
      const UF_long tail_start = repaired_block + 1u;
      int tail_ok = -1;
      tail_ok =
        kls_parallel_refactor_tail_from_block(solver, numeric_values,
                                              tail_start, 1);
      if (tail_ok < 0) {
        tail_ok =
          solver->common.scale > 0
            ? kls_serial_checked_scaled_refactor_from_block(
                solver, numeric_values, tail_start)
            : kls_serial_refactor_tail_from_block(solver, numeric_values,
                                                  tail_start, 1);
      }
      if (tail_ok < 0 || solver->common.status < 0 ||
          solver->common.status == TRILINOS_KLU_SINGULAR) {
        return 0;
      }
      if (!tail_ok) {
        if (solver->stats.fast_rejected_pivot >= 0) {
          continue;
        }
        return 0;
      }
      UF_long next_rejected_pivot = KLS_KLU_EMPTY;
      UF_long rejected_pivot_col = KLS_KLU_EMPTY;
      if (!kls_numeric_pivots_pass_threshold_from_block(
            solver, repaired_block, &next_rejected_pivot,
            &rejected_pivot_col)) {
        kls_record_fast_reject(solver, next_rejected_pivot,
                               rejected_pivot_col);
        continue;
      }
      if (!kls_checked_refactor_accepts_rowwise_u_from_block(
            solver, repaired_block)) {
        if (solver->common.status < 0) {
          return 0;
        }
        continue;
      }
      return 1;
    }
  }
  return 0;
}
