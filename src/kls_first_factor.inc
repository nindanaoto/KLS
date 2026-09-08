/* Row-first and predicted first-factor construction.
 * Included by kls_numeric_engines.inc in the original definition order.
 */

static size_t kls_initial_block_lusize(const kls_solver *solver,
                                       UF_long k1,
                                       UF_long nk,
                                       double lsize_estimate) {
  if (solver == NULL || solver->col_ptr == NULL || nk == 0u ||
      k1 > solver->n || nk > solver->n - k1) {
    return 0u;
  }

  double lsize = lsize_estimate;
  const double n = (double)nk;
  if (lsize <= 0.0) {
    const UF_long p1 = solver->col_ptr[k1];
    const UF_long p2 = solver->col_ptr[k1 + nk];
    if (p2 < p1) {
      return 0u;
    }
    lsize = -lsize;
    lsize = lsize > 1.0 ? lsize : 1.0;
    lsize = lsize * (double)(p2 - p1) + n;
  }

  lsize = lsize > n + 1.0 ? lsize : n + 1.0;
  double max_lnz = (n * n + n) / 2.0;
  max_lnz = max_lnz < (double)INT_MAX ? max_lnz : (double)INT_MAX;
  lsize = lsize < max_lnz ? lsize : max_lnz;
  if (INT_OVERFLOW(lsize)) {
    return 0u;
  }

  const Int entries = (Int)lsize;
  const double units = DUNITS(Int, entries) + DUNITS(Entry, entries) +
                       DUNITS(Int, entries) + DUNITS(Entry, entries);
  if (INT_OVERFLOW(units) || units <= 0.0) {
    return 0u;
  }
  return (size_t)units;
}

static trilinos_klu_l_numeric *kls_allocate_numeric_skeleton(
  kls_solver *solver) {
  if (solver == NULL || solver->symbolic == NULL) {
    return NULL;
  }
  TRILINOS_KLU_common *common = &solver->common;
  const UF_long n = solver->symbolic->n;
  const UF_long nblocks = solver->symbolic->nblocks;
  const UF_long nzoff = solver->symbolic->nzoff;
  const UF_long maxblock = solver->symbolic->maxblock;

  Int ok = TRUE;
  const size_t n1 = TRILINOS_KLU_add_size_t((size_t)n, 1u, &ok);
  const size_t nzoff1 = TRILINOS_KLU_add_size_t((size_t)nzoff, 1u, &ok);
  trilinos_klu_l_numeric *numeric =
    (trilinos_klu_l_numeric *)TRILINOS_KLU_malloc(
      1u, sizeof(*numeric), common);
  if (numeric == NULL || common->status < TRILINOS_KLU_OK || !ok) {
    common->status =
      ok ? TRILINOS_KLU_OUT_OF_MEMORY : TRILINOS_KLU_TOO_LARGE;
    return NULL;
  }
  memset(numeric, 0, sizeof(*numeric));
  numeric->n = n;
  numeric->nblocks = nblocks;
  numeric->nzoff = nzoff;
  numeric->Pnum = (UF_long *)TRILINOS_KLU_malloc(n, sizeof(UF_long), common);
  numeric->Offp = (UF_long *)TRILINOS_KLU_malloc(n1, sizeof(UF_long), common);
  numeric->Offi =
    (UF_long *)TRILINOS_KLU_malloc(nzoff1, sizeof(UF_long), common);
  numeric->Offx = TRILINOS_KLU_malloc(nzoff1, sizeof(Entry), common);
  numeric->Lip = (UF_long *)TRILINOS_KLU_malloc(n, sizeof(UF_long), common);
  numeric->Uip = (UF_long *)TRILINOS_KLU_malloc(n, sizeof(UF_long), common);
  numeric->Llen = (UF_long *)TRILINOS_KLU_malloc(n, sizeof(UF_long), common);
  numeric->Ulen = (UF_long *)TRILINOS_KLU_malloc(n, sizeof(UF_long), common);
  numeric->LUsize =
    (size_t *)TRILINOS_KLU_malloc(nblocks, sizeof(size_t), common);
  numeric->LUbx = (void **)TRILINOS_KLU_malloc(nblocks, sizeof(Unit *), common);
  if (numeric->LUbx != NULL) {
    for (UF_long block = 0; block < nblocks; ++block) {
      numeric->LUbx[block] = NULL;
    }
  }
  if (numeric->LUsize != NULL) {
    for (UF_long block = 0; block < nblocks; ++block) {
      numeric->LUsize[block] = 0u;
    }
  }
  numeric->Udiag = TRILINOS_KLU_malloc(n, sizeof(Entry), common);
  numeric->Rs = NULL;
  if (common->scale > 0) {
    numeric->Rs =
      (double *)TRILINOS_KLU_malloc(n, sizeof(double), common);
  }
  numeric->Pinv = (UF_long *)TRILINOS_KLU_malloc(n, sizeof(UF_long), common);

  const size_t xwork_size = TRILINOS_KLU_mult_size_t(n, sizeof(Entry), &ok);
  const size_t solve_work_size =
    TRILINOS_KLU_mult_size_t(n, 3u * sizeof(Entry), &ok);
  const size_t factor_iwork_size =
    TRILINOS_KLU_mult_size_t(maxblock, 6u * sizeof(UF_long), &ok);
  numeric->worksize =
    TRILINOS_KLU_add_size_t(xwork_size,
                            solve_work_size > factor_iwork_size
                              ? solve_work_size : factor_iwork_size,
                            &ok);
  numeric->Work = TRILINOS_KLU_malloc(numeric->worksize, 1u, common);
  numeric->Xwork = numeric->Work;
  numeric->Iwork = (UF_long *)(((Entry *)numeric->Xwork) + n);

  if (!ok || common->status < TRILINOS_KLU_OK ||
      numeric->Pnum == NULL || numeric->Offp == NULL ||
      numeric->Offi == NULL || numeric->Offx == NULL ||
      numeric->Lip == NULL || numeric->Uip == NULL ||
      numeric->Llen == NULL || numeric->Ulen == NULL ||
      numeric->LUsize == NULL || numeric->LUbx == NULL ||
      numeric->Udiag == NULL || numeric->Pinv == NULL ||
      numeric->Work == NULL || (common->scale > 0 && numeric->Rs == NULL)) {
    common->status =
      ok ? TRILINOS_KLU_OUT_OF_MEMORY : TRILINOS_KLU_TOO_LARGE;
    trilinos_klu_l_free_numeric(&numeric, common);
    return NULL;
  }
  return numeric;
}

typedef struct kls_row_first_entries {
  UF_long *row;
  UF_long *col;
  UF_long *next;
  UF_long *col_head;
  UF_long *col_count;
  double *value;
  UF_long count;
  UF_long capacity;
  UF_long col_head_len;
  UF_long col_count_len;
} kls_row_first_entries;

static void kls_row_first_entries_free(kls_row_first_entries *entries) {
  if (entries == NULL) {
    return;
  }
  free(entries->row);
  free(entries->col);
  free(entries->next);
  free(entries->col_head);
  free(entries->col_count);
  free(entries->value);
  entries->row = NULL;
  entries->col = NULL;
  entries->next = NULL;
  entries->col_head = NULL;
  entries->col_count = NULL;
  entries->value = NULL;
  entries->count = 0;
  entries->capacity = 0;
  entries->col_head_len = 0;
  entries->col_count_len = 0;
}

static int kls_row_first_entries_enable_col_counts(
  kls_row_first_entries *entries,
  UF_long local_n) {
  if (entries == NULL || entries->count != 0u ||
      entries->capacity != 0u || entries->col_count != NULL) {
    return 0;
  }
  entries->col_count =
    (UF_long *)calloc((size_t)local_n, sizeof(*entries->col_count));
  if (entries->col_count == NULL && local_n > 0u) {
    return 0;
  }
  entries->col_count_len = local_n;
  return 1;
}

static int kls_row_first_entries_enable_col_links(
  kls_row_first_entries *entries,
  UF_long local_n) {
  if (entries == NULL || entries->count != 0u ||
      entries->capacity != 0u || entries->col_head != NULL ||
      entries->next != NULL) {
    return 0;
  }
  entries->col_head =
    (UF_long *)malloc((size_t)local_n * sizeof(*entries->col_head));
  if (entries->col_head == NULL && local_n > 0u) {
    return 0;
  }
  for (UF_long k = 0; k < local_n; ++k) {
    entries->col_head[k] = KLS_KLU_EMPTY;
  }
  entries->col_head_len = local_n;
  return 1;
}

static int kls_row_first_entries_reserve(kls_row_first_entries *entries,
                                         UF_long capacity) {
  if (entries == NULL || capacity <= entries->capacity) {
    return 1;
  }
  if (capacity > (UF_long)(SIZE_MAX / sizeof(*entries->row)) ||
      capacity > (UF_long)(SIZE_MAX / sizeof(*entries->value))) {
    return 0;
  }

  UF_long *new_row =
    (UF_long *)malloc((size_t)capacity * sizeof(*new_row));
  UF_long *new_col =
    (UF_long *)malloc((size_t)capacity * sizeof(*new_col));
  double *new_value =
    (double *)malloc((size_t)capacity * sizeof(*new_value));
  UF_long *new_next = entries->col_head != NULL
    ? (UF_long *)malloc((size_t)capacity * sizeof(*new_next))
    : NULL;
  if (new_row == NULL || new_col == NULL || new_value == NULL ||
      (entries->col_head != NULL && new_next == NULL)) {
    free(new_row);
    free(new_col);
    free(new_next);
    free(new_value);
    return 0;
  }

  if (entries->count > 0u) {
    memcpy(new_row, entries->row,
           (size_t)entries->count * sizeof(*new_row));
    memcpy(new_col, entries->col,
           (size_t)entries->count * sizeof(*new_col));
    memcpy(new_value, entries->value,
           (size_t)entries->count * sizeof(*new_value));
    if (entries->col_head != NULL) {
      memcpy(new_next, entries->next,
             (size_t)entries->count * sizeof(*new_next));
    }
  }

  free(entries->row);
  free(entries->col);
  free(entries->next);
  free(entries->value);
  entries->row = new_row;
  entries->col = new_col;
  entries->next = new_next;
  entries->value = new_value;
  entries->capacity = capacity;
  return 1;
}

static int kls_row_first_entries_reserve_append(
  kls_row_first_entries *entries,
  UF_long append_count) {
  if (entries == NULL) {
    return 0;
  }
  if (append_count == 0u) {
    return 1;
  }
  if (append_count > UF_long_max - entries->count) {
    return 0;
  }
  const UF_long needed = entries->count + append_count;
  if (needed <= entries->capacity) {
    return 1;
  }
  UF_long grown = entries->capacity == 0u ? 64u : entries->capacity;
  while (grown < needed) {
    if (grown > UF_long_max / 2u) {
      grown = needed;
      break;
    }
    grown *= 2u;
  }
  return kls_row_first_entries_reserve(entries, grown);
}

static int kls_row_first_entries_append_reserved(
  kls_row_first_entries *entries,
  UF_long row,
  UF_long col,
  double value) {
  if (entries == NULL || entries->count >= entries->capacity ||
      (entries->col_head != NULL && col >= entries->col_head_len) ||
      (entries->col_count != NULL && col >= entries->col_count_len)) {
    return 0;
  }
  const UF_long pos = entries->count;
  entries->row[pos] = row;
  entries->col[pos] = col;
  entries->value[pos] = value;
  if (entries->col_count != NULL) {
    entries->col_count[col]++;
  }
  if (entries->col_head != NULL) {
    entries->next[pos] = entries->col_head[col];
    entries->col_head[col] = pos;
  }
  entries->count++;
  return 1;
}

static UF_long kls_row_first_symbolic_reserve(
  const kls_solver *solver,
  UF_long block,
  UF_long nk) {
  if (solver == NULL || solver->symbolic == NULL ||
      solver->symbolic->Lnz == NULL || block >= solver->symbolic->nblocks) {
    return 0;
  }
  const double estimate = solver->symbolic->Lnz[block];
  if (!isfinite(estimate) || estimate <= (double)nk ||
      estimate > (double)UF_long_max) {
    return 0;
  }
  const UF_long entries = (UF_long)(estimate - (double)nk);
  return entries >= 64u ? entries : 0u;
}

static int kls_row_first_entries_append(kls_row_first_entries *entries,
                                        UF_long row,
                                        UF_long col,
                                        double value) {
  if (!kls_row_first_entries_reserve_append(entries, 1u)) {
    return 0;
  }
  return kls_row_first_entries_append_reserved(entries, row, col, value);
}

static int kls_row_first_entries_append_consecutive(
  kls_row_first_entries *entries,
  UF_long row,
  UF_long col_begin,
  const double *values,
  UF_long count) {
  if (count == 0u) {
    return 1;
  }
  if (entries == NULL || values == NULL ||
      count > UF_long_max - col_begin ||
      count > UF_long_max - entries->count ||
      !kls_row_first_entries_reserve_append(entries, count)) {
    return 0;
  }
  const UF_long col_end = col_begin + count;
  if ((entries->col_head != NULL && col_end > entries->col_head_len) ||
      (entries->col_count != NULL && col_end > entries->col_count_len)) {
    return 0;
  }
  const UF_long pos_begin = entries->count;
  for (UF_long local = 0; local < count; ++local) {
    const UF_long pos = pos_begin + local;
    const UF_long col = col_begin + local;
    entries->row[pos] = row;
    entries->col[pos] = col;
    entries->value[pos] = values[local];
    if (entries->col_count != NULL) {
      entries->col_count[col]++;
    }
    if (entries->col_head != NULL) {
      entries->next[pos] = entries->col_head[col];
      entries->col_head[col] = pos;
    }
  }
  entries->count = pos_begin + count;
  return 1;
}

static int kls_row_first_entries_truncate(kls_row_first_entries *entries,
                                          UF_long count) {
  if (entries == NULL || count > entries->count) {
    return 0;
  }
  if (count == entries->count) {
    return 1;
  }
  if (entries->col_count != NULL) {
    for (UF_long p = count; p < entries->count; ++p) {
      const UF_long col = entries->col[p];
      if (col >= entries->col_count_len || entries->col_count[col] == 0u) {
        return 0;
      }
      entries->col_count[col]--;
    }
  }
  entries->count = count;
  if (entries->col_head != NULL) {
    if (entries->next == NULL) {
      return 0;
    }
    for (UF_long col = 0; col < entries->col_head_len; ++col) {
      entries->col_head[col] = KLS_KLU_EMPTY;
    }
    for (UF_long p = 0; p < entries->count; ++p) {
      const UF_long col = entries->col[p];
      if (col >= entries->col_head_len) {
        return 0;
      }
      entries->next[p] = entries->col_head[col];
      entries->col_head[col] = p;
    }
  }
  return 1;
}

static int kls_row_first_entries_swap_columns(
  kls_row_first_entries *u_entries,
  UF_long a,
  UF_long b,
  UF_long local_n) {
  if (a == b) {
    return 1;
  }
  if (u_entries == NULL || a >= local_n || b >= local_n) {
    return 0;
  }
  if (u_entries->col_head != NULL) {
    if ((u_entries->count > 0u && u_entries->next == NULL) ||
        a >= u_entries->col_head_len || b >= u_entries->col_head_len) {
      return 0;
    }
    const UF_long head_a = u_entries->col_head[a];
    const UF_long head_b = u_entries->col_head[b];
    UF_long guard = 0;
    for (UF_long p = head_a; p != KLS_KLU_EMPTY; p = u_entries->next[p]) {
      if (p >= u_entries->count || ++guard > u_entries->count) {
        return 0;
      }
      u_entries->col[p] = b;
      if (u_entries->col[p] <= u_entries->row[p] ||
          u_entries->row[p] >= local_n || u_entries->col[p] >= local_n) {
        return 0;
      }
    }
    guard = 0;
    for (UF_long p = head_b; p != KLS_KLU_EMPTY; p = u_entries->next[p]) {
      if (p >= u_entries->count || ++guard > u_entries->count) {
        return 0;
      }
      u_entries->col[p] = a;
      if (u_entries->col[p] <= u_entries->row[p] ||
          u_entries->row[p] >= local_n || u_entries->col[p] >= local_n) {
        return 0;
      }
    }
    u_entries->col_head[a] = head_b;
    u_entries->col_head[b] = head_a;
    if (u_entries->col_count != NULL) {
      if (a >= u_entries->col_count_len || b >= u_entries->col_count_len) {
        return 0;
      }
      const UF_long count_a = u_entries->col_count[a];
      u_entries->col_count[a] = u_entries->col_count[b];
      u_entries->col_count[b] = count_a;
    }
    return 1;
  }
  for (UF_long p = 0; p < u_entries->count; ++p) {
    if (u_entries->col[p] == a) {
      u_entries->col[p] = b;
    } else if (u_entries->col[p] == b) {
      u_entries->col[p] = a;
    }
    if (u_entries->col[p] <= u_entries->row[p] ||
        u_entries->row[p] >= local_n || u_entries->col[p] >= local_n) {
      return 0;
    }
  }
  if (u_entries->col_count != NULL) {
    if (a >= u_entries->col_count_len || b >= u_entries->col_count_len) {
      return 0;
    }
    const UF_long count_a = u_entries->col_count[a];
    u_entries->col_count[a] = u_entries->col_count[b];
    u_entries->col_count[b] = count_a;
  }
  return 1;
}

static void kls_row_first_heap_push(UF_long *heap,
                                    UF_long *size,
                                    UF_long value) {
  UF_long pos = (*size)++;
  while (pos > 0u) {
    const UF_long parent = (pos - 1u) / 2u;
    if (heap[parent] <= value) {
      break;
    }
    heap[pos] = heap[parent];
    pos = parent;
  }
  heap[pos] = value;
}

static UF_long kls_row_first_heap_pop(UF_long *heap, UF_long *size) {
  const UF_long result = heap[0];
  const UF_long value = heap[--(*size)];
  UF_long pos = 0;
  while (1) {
    const UF_long left = 2u * pos + 1u;
    const UF_long right = left + 1u;
    if (left >= *size) {
      break;
    }
    UF_long child = left;
    if (right < *size && heap[right] < heap[left]) {
      child = right;
    }
    if (heap[child] >= value) {
      break;
    }
    heap[pos] = heap[child];
    pos = child;
  }
  if (*size > 0u) {
    heap[pos] = value;
  }
  return result;
}

static int kls_row_first_heap_consume_ready_dep(UF_long *heap,
                                                UF_long *size,
                                                UF_long dep) {
  if (heap == NULL || size == NULL) {
    return 0;
  }
  if (*size == 0u) {
    return 1;
  }
  const UF_long root = heap[0];
  if (root == dep) {
    const UF_long popped = kls_row_first_heap_pop(heap, size);
    return popped == dep;
  }
  return root > dep;
}

static int kls_row_first_exchange_columns(kls_row_first_entries *u_entries,
                                          UF_long *q_order,
                                          UF_long *col_pos,
                                          UF_long a,
                                          UF_long b,
                                          UF_long local_n,
                                          UF_long global_n) {
  if (a == b) {
    return 1;
  }
  if (u_entries == NULL || q_order == NULL || col_pos == NULL ||
      a >= local_n || b >= local_n ||
      q_order[a] >= global_n || q_order[b] >= global_n) {
    return 0;
  }
  const UF_long old_a = q_order[a];
  const UF_long old_b = q_order[b];
  q_order[a] = old_b;
  q_order[b] = old_a;
  col_pos[old_b] = a;
  col_pos[old_a] = b;
  return kls_row_first_entries_swap_columns(u_entries, a, b, local_n);
}

static int kls_pack_row_first_block_numeric(
  kls_solver *solver,
  trilinos_klu_l_numeric *numeric,
  UF_long block,
  UF_long k1,
  UF_long nk,
  const UF_long *local_row_order,
  const kls_row_first_entries *l_entries,
  const kls_row_first_entries *u_entries,
  const double *udiag_values,
  UF_long *total_lnz,
  UF_long *total_unz,
  UF_long *max_lnz_block,
  UF_long *max_unz_block) {
  if (solver == NULL || numeric == NULL || l_entries == NULL ||
      u_entries == NULL || udiag_values == NULL ||
      numeric->Lip == NULL || numeric->Llen == NULL ||
      numeric->Uip == NULL || numeric->Ulen == NULL ||
      numeric->Udiag == NULL || numeric->LUbx == NULL ||
      numeric->LUsize == NULL || numeric->Pnum == NULL ||
      solver->symbolic == NULL || solver->symbolic->P == NULL ||
      block >= solver->symbolic->nblocks || k1 > solver->n ||
      nk == 0u || nk > solver->n - k1 ||
      total_lnz == NULL || total_unz == NULL ||
      max_lnz_block == NULL || max_unz_block == NULL) {
    return 0;
  }

  if (l_entries->col_count == NULL || u_entries->col_count == NULL ||
      l_entries->col_count_len < nk || u_entries->col_count_len < nk) {
    return 0;
  }
  const UF_long *l_count = l_entries->col_count;
  const UF_long *u_count = u_entries->col_count;

  size_t lusize = 0u;
  for (UF_long k = 0; k < nk; ++k) {
    const size_t l_units = kls_klu_units_for_indices(l_count[k]);
    const size_t u_units = kls_klu_units_for_indices(u_count[k]);
    if (l_units > SIZE_MAX - lusize ||
        (size_t)l_count[k] > SIZE_MAX - lusize - l_units) {
      return 0;
    }
    lusize += l_units + (size_t)l_count[k];
    if (u_units > SIZE_MAX - lusize ||
        (size_t)u_count[k] > SIZE_MAX - lusize - u_units) {
      return 0;
    }
    lusize += u_units + (size_t)u_count[k];
  }
  if (lusize == 0u) {
    lusize = 1u;
  }

  Unit *lu =
    (Unit *)TRILINOS_KLU_malloc(lusize, sizeof(Unit), &solver->common);
  if (lu == NULL) {
    return 0;
  }

  size_t lup = 0u;
  for (UF_long k = 0; k < nk; ++k) {
    numeric->Lip[k1 + k] = (UF_long)lup;
    numeric->Llen[k1 + k] = l_count[k];
    lup += kls_klu_units_for_indices(l_count[k]) + (size_t)l_count[k];
    numeric->Uip[k1 + k] = (UF_long)lup;
    numeric->Ulen[k1 + k] = u_count[k];
    lup += kls_klu_units_for_indices(u_count[k]) + (size_t)u_count[k];
    ((Entry *)numeric->Udiag)[k1 + k] = udiag_values[k];
    const UF_long symbolic_row = local_row_order != NULL
      ? local_row_order[k] : k;
    if (symbolic_row >= nk) {
      (void)TRILINOS_KLU_free(lu, lusize, sizeof(Unit), &solver->common);
      return 0;
    }
    numeric->Pnum[k1 + k] = solver->symbolic->P[k1 + symbolic_row];
  }

  UF_long *l_write =
    (UF_long *)calloc((size_t)nk, sizeof(*l_write));
  UF_long *u_write =
    (UF_long *)calloc((size_t)nk, sizeof(*u_write));
  if (l_write == NULL || u_write == NULL) {
    (void)TRILINOS_KLU_free(lu, lusize, sizeof(Unit), &solver->common);
    free(l_write);
    free(u_write);
    return 0;
  }

  for (UF_long p = 0; p < l_entries->count; ++p) {
    const UF_long col = l_entries->col[p];
    if (col >= nk || l_entries->row[p] >= nk) {
      (void)TRILINOS_KLU_free(lu, lusize, sizeof(Unit), &solver->common);
      free(l_write);
      free(u_write);
      return 0;
    }
    const UF_long offset = l_write[col]++;
    if (offset >= l_count[col]) {
      (void)TRILINOS_KLU_free(lu, lusize, sizeof(Unit), &solver->common);
      free(l_write);
      free(u_write);
      return 0;
    }
    UF_long *li = (UF_long *)(lu + numeric->Lip[k1 + col]);
    Entry *lx =
      (Entry *)(lu + numeric->Lip[k1 + col] +
                kls_klu_units_for_indices(l_count[col]));
    li[offset] = l_entries->row[p];
    lx[offset] = l_entries->value[p];
  }
  for (UF_long p = 0; p < u_entries->count; ++p) {
    const UF_long col = u_entries->col[p];
    if (col >= nk || u_entries->row[p] >= nk) {
      (void)TRILINOS_KLU_free(lu, lusize, sizeof(Unit), &solver->common);
      free(l_write);
      free(u_write);
      return 0;
    }
    const UF_long offset = u_write[col]++;
    if (offset >= u_count[col]) {
      (void)TRILINOS_KLU_free(lu, lusize, sizeof(Unit), &solver->common);
      free(l_write);
      free(u_write);
      return 0;
    }
    UF_long *ui = (UF_long *)(lu + numeric->Uip[k1 + col]);
    Entry *ux =
      (Entry *)(lu + numeric->Uip[k1 + col] +
                kls_klu_units_for_indices(u_count[col]));
    ui[offset] = u_entries->row[p];
    ux[offset] = u_entries->value[p];
  }

  const UF_long lnz_block = l_entries->count + nk;
  const UF_long unz_block = u_entries->count + nk;
  numeric->LUbx[block] = lu;
  numeric->LUsize[block] = lusize;
  *total_lnz += lnz_block;
  *total_unz += unz_block;
  if (*max_lnz_block < lnz_block) {
    *max_lnz_block = lnz_block;
  }
  if (*max_unz_block < unz_block) {
    *max_unz_block = unz_block;
  }
  free(l_write);
  free(u_write);
  return 1;
}

typedef struct kls_row_first_lu_seed {
  UF_long *row;
  UF_long *col;
  double **numeric_value;
  double *value;
  UF_long count;
  UF_long capacity;
} kls_row_first_lu_seed;

typedef struct kls_row_first_input_seed {
  UF_long *row;
  UF_long *col;
  UF_long *pos;
  UF_long count;
  UF_long capacity;
} kls_row_first_input_seed;

typedef struct kls_row_first_refactor_seed {
  UF_long n;
  int enabled;
  kls_row_first_lu_seed l;
  kls_row_first_lu_seed u;
  kls_row_first_input_seed input;
} kls_row_first_refactor_seed;

typedef struct kls_first_separator_component_entry {
  unsigned int component;
  UF_long rows;
  double work;
} kls_first_separator_component_entry;

static int kls_compare_first_separator_component_work_desc(
  const void *a,
  const void *b) {
  const kls_first_separator_component_entry *left =
    (const kls_first_separator_component_entry *)a;
  const kls_first_separator_component_entry *right =
    (const kls_first_separator_component_entry *)b;
  if (left->work < right->work) {
    return 1;
  }
  if (left->work > right->work) {
    return -1;
  }
  if (left->rows < right->rows) {
    return 1;
  }
  if (left->rows > right->rows) {
    return -1;
  }
  return (left->component > right->component) -
         (left->component < right->component);
}

static void kls_first_separator_queue_plan_free(
  kls_first_separator_queue_plan *plan) {
  if (plan == NULL) {
    return;
  }
  free(plan->private_rows);
  free(plan->pipeline_rows);
  free(plan->thread_ptr);
  memset(plan, 0, sizeof(*plan));
}

static double kls_first_separator_row_input_work(UF_long input_count) {
  const double count = (double)input_count;
  const char *model = getenv("KLS_FIRST_SEPARATOR_WORK_MODEL");
  if (model != NULL) {
    if (strcmp(model, "rows") == 0) {
      return 1.0;
    }
    if (strcmp(model, "linear") == 0) {
      return count;
    }
  }
  return count * count;
}

static int kls_build_first_separator_queue_plan_partitioned(
  const kls_solver *solver,
  UF_long k1,
  UF_long k2,
  int thread_count,
  const UF_long *row_input_count,
  kls_first_separator_queue_plan *plan) {
  if (solver == NULL || plan == NULL || row_input_count == NULL ||
      thread_count <= 1 || k2 <= k1 || k2 > solver->n ||
      !kls_separator_analysis_has_global_range(&solver->separator) ||
      solver->separator.component_count <= 1u ||
      solver->separator.order_component == NULL ||
      solver->separator.component_left_child == NULL ||
      solver->separator.component_right_child == NULL ||
      k1 < solver->separator.global_begin ||
      k2 > solver->separator.global_end) {
    return 0;
  }

  const kls_separator_analysis *separator = &solver->separator;
  const UF_long component_count = separator->component_count;
  const UF_long local_begin = k1 - separator->global_begin;
  const UF_long local_end = k2 - separator->global_begin;
  if (local_end > separator->n || component_count > (UF_long)UINT_MAX) {
    return 0;
  }
  const UF_long block_rows = k2 - k1;
  const size_t component_count_size = (size_t)component_count;
  const size_t thread_count_size = (size_t)thread_count;
  if ((UF_long)component_count_size != component_count ||
      (UF_long)thread_count_size != (UF_long)thread_count ||
      component_count_size > SIZE_MAX / sizeof(UF_long) ||
      component_count_size > SIZE_MAX / sizeof(double) ||
      component_count_size > SIZE_MAX / sizeof(unsigned char) ||
      component_count_size > SIZE_MAX / sizeof(int) ||
      component_count_size >
        SIZE_MAX / sizeof(kls_first_separator_component_entry) ||
      block_rows > (UF_long)(SIZE_MAX / sizeof(UF_long)) ||
      thread_count_size > (SIZE_MAX / sizeof(UF_long)) - 1u ||
      thread_count_size > SIZE_MAX / sizeof(double)) {
    return 0;
  }

  int ok = 0;
  UF_long *component_rows =
    (UF_long *)calloc(component_count_size, sizeof(*component_rows));
  double *component_work =
    (double *)calloc(component_count_size, sizeof(*component_work));
  UF_long *subtree_begin =
    (UF_long *)malloc(component_count_size * sizeof(*subtree_begin));
  UF_long *subtree_end =
    (UF_long *)malloc(component_count_size * sizeof(*subtree_end));
  double *subtree_work =
    (double *)calloc(component_count_size, sizeof(*subtree_work));
  UF_long *candidates =
    (UF_long *)malloc(component_count_size * sizeof(*candidates));
  unsigned char *component_pipeline =
    (unsigned char *)calloc(component_count_size, sizeof(*component_pipeline));
  int *component_thread =
    (int *)malloc(component_count_size * sizeof(*component_thread));
  kls_first_separator_component_entry *candidate_entries =
    (kls_first_separator_component_entry *)
      malloc(component_count_size * sizeof(*candidate_entries));
  UF_long *thread_counts =
    (UF_long *)calloc(thread_count_size, sizeof(*thread_counts));
  double *thread_work =
    (double *)calloc(thread_count_size, sizeof(*thread_work));
  UF_long *thread_ptr = NULL;
  UF_long *private_queue = NULL;
  UF_long *pipeline_queue = NULL;
  if (component_rows == NULL || component_work == NULL ||
      subtree_begin == NULL || subtree_end == NULL ||
      subtree_work == NULL || candidates == NULL ||
      component_pipeline == NULL || component_thread == NULL ||
      candidate_entries == NULL || thread_counts == NULL ||
      thread_work == NULL) {
    goto cleanup;
  }
  for (UF_long component = 0; component < component_count; ++component) {
    component_thread[component] = -1;
  }

  for (UF_long local = local_begin; local < local_end; ++local) {
    const unsigned int component = separator->order_component[local];
    const UF_long global_row = separator->global_begin + local;
    if ((UF_long)component >= component_count ||
        global_row < k1 || global_row >= k2) {
      goto cleanup;
    }
    component_rows[component]++;
    component_work[component] +=
      kls_first_separator_row_input_work(row_input_count[global_row - k1]);
  }

  for (UF_long component = 0; component < component_count; ++component) {
    UF_long begin = component;
    UF_long end = component + 1u;
    double work = component_work[component];
    const UF_long left = separator->component_left_child[component];
    const UF_long right = separator->component_right_child[component];
    if (left != KLS_KLU_EMPTY) {
      if (left >= component || left >= component_count) {
        goto cleanup;
      }
      begin = subtree_begin[left] < begin ? subtree_begin[left] : begin;
      end = subtree_end[left] > end ? subtree_end[left] : end;
      work += subtree_work[left];
    }
    if (right != KLS_KLU_EMPTY) {
      if (right >= component || right >= component_count) {
        goto cleanup;
      }
      begin = subtree_begin[right] < begin ? subtree_begin[right] : begin;
      end = subtree_end[right] > end ? subtree_end[right] : end;
      work += subtree_work[right];
    }
    subtree_begin[component] = begin;
    subtree_end[component] = end;
    subtree_work[component] = work;
  }

  const UF_long root_component = component_count - 1u;
  if (subtree_work[root_component] <= 0.0) {
    goto cleanup;
  }
  UF_long candidate_count = 1u;
  UF_long pipeline_component_count = 0;
  candidates[0] = root_component;
  double candidate_work_sum = subtree_work[root_component];
  while (candidate_count > 0u && candidate_work_sum > 0.0) {
    UF_long max_pos = 0;
    double max_work = subtree_work[candidates[0]];
    for (UF_long pos = 1u; pos < candidate_count; ++pos) {
      const double work = subtree_work[candidates[pos]];
      if (work > max_work) {
        max_pos = pos;
        max_work = work;
      }
    }
    if ((double)thread_count * max_work <
        KLS_ROW_REFACTOR_SEPARATOR_BALANCE_BETA * candidate_work_sum) {
      break;
    }
    const UF_long component = candidates[max_pos];
    const UF_long left = separator->component_left_child[component];
    const UF_long right = separator->component_right_child[component];
    const int split_left =
      left != KLS_KLU_EMPTY && left < component_count &&
      subtree_work[left] > 0.0;
    const int split_right =
      right != KLS_KLU_EMPTY && right < component_count &&
      subtree_work[right] > 0.0;
    if (!split_left && !split_right) {
      break;
    }
    candidates[max_pos] = candidates[candidate_count - 1u];
    candidate_count--;
    candidate_work_sum -= subtree_work[component];
    component_pipeline[component] = 1u;
    pipeline_component_count++;
    if (split_left) {
      candidates[candidate_count++] = left;
      candidate_work_sum += subtree_work[left];
    }
    if (split_right) {
      candidates[candidate_count++] = right;
      candidate_work_sum += subtree_work[right];
    }
  }
  if (candidate_count == 0u || pipeline_component_count == 0u) {
    goto cleanup;
  }

  for (UF_long pos = 0; pos < candidate_count; ++pos) {
    if (candidates[pos] > (UF_long)UINT_MAX) {
      goto cleanup;
    }
    candidate_entries[pos].component = (unsigned int)candidates[pos];
    candidate_entries[pos].rows = 0;
    candidate_entries[pos].work = subtree_work[candidates[pos]];
  }
  qsort(candidate_entries, (size_t)candidate_count,
        sizeof(*candidate_entries),
        kls_compare_first_separator_component_work_desc);

  for (UF_long pos = 0; pos < candidate_count; ++pos) {
    const UF_long component = candidate_entries[pos].component;
    int target = 0;
    double target_work = thread_work[0];
    for (int tid = 1; tid < thread_count; ++tid) {
      if (thread_work[tid] < target_work) {
        target = tid;
        target_work = thread_work[tid];
      }
    }
    for (UF_long c = subtree_begin[component]; c < subtree_end[component];
         ++c) {
      if (c >= component_count || component_pipeline[c]) {
        continue;
      }
      component_thread[c] = target;
    }
    thread_work[target] += subtree_work[component];
  }

  UF_long private_rows = 0;
  UF_long pipeline_rows = 0;
  for (UF_long local = local_begin; local < local_end; ++local) {
    const unsigned int component = separator->order_component[local];
    if ((UF_long)component >= component_count) {
      goto cleanup;
    }
    if (component_pipeline[component]) {
      pipeline_rows++;
      continue;
    }
    const int target = component_thread[component];
    if (target < 0 || target >= thread_count) {
      goto cleanup;
    }
    thread_counts[target]++;
    private_rows++;
  }
  if (private_rows == 0u || pipeline_rows == 0u ||
      private_rows + pipeline_rows != block_rows) {
    goto cleanup;
  }

  UF_long nonempty_threads = 0;
  UF_long max_thread_rows = 0;
  double min_thread_work = 0.0;
  double max_thread_work = 0.0;
  for (int tid = 0; tid < thread_count; ++tid) {
    if (thread_counts[tid] == 0u) {
      continue;
    }
    nonempty_threads++;
    if (thread_counts[tid] > max_thread_rows) {
      max_thread_rows = thread_counts[tid];
    }
    if (min_thread_work == 0.0 || thread_work[tid] < min_thread_work) {
      min_thread_work = thread_work[tid];
    }
    if (thread_work[tid] > max_thread_work) {
      max_thread_work = thread_work[tid];
    }
  }
  if (nonempty_threads == 0u) {
    goto cleanup;
  }

  thread_ptr =
    (UF_long *)calloc(thread_count_size + 1u, sizeof(*thread_ptr));
  private_queue =
    (UF_long *)malloc((size_t)private_rows * sizeof(*private_queue));
  pipeline_queue =
    (UF_long *)malloc((size_t)pipeline_rows * sizeof(*pipeline_queue));
  if (thread_ptr == NULL || private_queue == NULL ||
      pipeline_queue == NULL) {
    goto cleanup;
  }
  for (int tid = 0; tid < thread_count; ++tid) {
    thread_ptr[tid + 1] = thread_ptr[tid] + thread_counts[tid];
    thread_counts[tid] = thread_ptr[tid];
  }
  if (thread_ptr[thread_count] != private_rows) {
    goto cleanup;
  }

  UF_long pipeline_next = 0;
  for (UF_long local = local_begin; local < local_end; ++local) {
    const unsigned int component = separator->order_component[local];
    const UF_long global_row = separator->global_begin + local;
    if (component_pipeline[component]) {
      if (pipeline_next >= pipeline_rows) {
        goto cleanup;
      }
      pipeline_queue[pipeline_next++] = global_row;
    } else {
      const int target = component_thread[component];
      if (target < 0 || target >= thread_count) {
        goto cleanup;
      }
      const UF_long dst = thread_counts[target]++;
      if (dst >= private_rows) {
        goto cleanup;
      }
      private_queue[dst] = global_row;
    }
  }
  if (pipeline_next != pipeline_rows) {
    goto cleanup;
  }
  for (int tid = 0; tid < thread_count; ++tid) {
    if (thread_counts[tid] != thread_ptr[tid + 1]) {
      goto cleanup;
    }
  }

  plan->private_rows = private_queue;
  plan->pipeline_rows = pipeline_queue;
  plan->thread_ptr = thread_ptr;
  plan->private_rows_count = private_rows;
  plan->pipeline_rows_count = pipeline_rows;
  plan->private_component_count = candidate_count;
  plan->pipeline_component_count = pipeline_component_count;
  plan->nonempty_threads = nonempty_threads;
  plan->max_thread_rows = max_thread_rows;
  plan->min_thread_work = min_thread_work;
  plan->max_thread_work = max_thread_work;
  plan->partitioned = 1;
  plan->split_component_count = pipeline_component_count;
  private_queue = NULL;
  pipeline_queue = NULL;
  thread_ptr = NULL;
  ok = 1;

cleanup:
  free(component_rows);
  free(component_work);
  free(subtree_begin);
  free(subtree_end);
  free(subtree_work);
  free(candidates);
  free(component_pipeline);
  free(component_thread);
  free(candidate_entries);
  free(thread_counts);
  free(thread_work);
  free(thread_ptr);
  free(private_queue);
  free(pipeline_queue);
  return ok;
}

static int kls_build_first_separator_queue_plan_legacy_from_counts(
  const kls_solver *solver,
  UF_long k1,
  UF_long k2,
  int thread_count,
  const UF_long *row_input_count,
  kls_first_separator_queue_plan *plan) {
  if (plan != NULL) {
    memset(plan, 0, sizeof(*plan));
  }
  if (solver == NULL || plan == NULL || row_input_count == NULL ||
      thread_count <= 1 || k2 <= k1 || k2 > solver->n ||
      !kls_separator_analysis_has_global_range(&solver->separator) ||
      solver->separator.component_count <= 1u ||
      solver->separator.component_kind == NULL ||
      solver->separator.order_component == NULL ||
      k1 < solver->separator.global_begin ||
      k2 > solver->separator.global_end) {
    return 0;
  }

  const kls_separator_analysis *separator = &solver->separator;
  const UF_long component_count = separator->component_count;
  const UF_long local_begin = k1 - separator->global_begin;
  const UF_long local_end = k2 - separator->global_begin;
  if (local_end > separator->n || component_count > (UF_long)UINT_MAX) {
    return 0;
  }
  const size_t component_count_size = (size_t)component_count;
  const size_t thread_count_size = (size_t)thread_count;
  const UF_long block_rows = k2 - k1;
  if ((UF_long)component_count_size != component_count ||
      (UF_long)thread_count_size != (UF_long)thread_count ||
      component_count_size > SIZE_MAX / sizeof(UF_long) ||
      component_count_size > SIZE_MAX / sizeof(int) ||
      component_count_size > SIZE_MAX / sizeof(double) ||
      component_count_size >
        SIZE_MAX / sizeof(kls_first_separator_component_entry) ||
      block_rows > (UF_long)(SIZE_MAX / sizeof(UF_long)) ||
      thread_count_size > (SIZE_MAX / sizeof(UF_long)) - 1u ||
      thread_count_size > SIZE_MAX / sizeof(double)) {
    return 0;
  }

  int ok = 0;
  UF_long *component_private_rows =
    (UF_long *)calloc(component_count_size, sizeof(*component_private_rows));
  UF_long *component_pipeline_rows =
    (UF_long *)calloc(component_count_size, sizeof(*component_pipeline_rows));
  double *component_work =
    (double *)calloc(component_count_size, sizeof(*component_work));
  int *component_thread =
    (int *)malloc(component_count_size * sizeof(*component_thread));
  kls_first_separator_component_entry *private_components =
    (kls_first_separator_component_entry *)
      malloc(component_count_size * sizeof(*private_components));
  UF_long *thread_counts =
    (UF_long *)calloc(thread_count_size, sizeof(*thread_counts));
  double *thread_work =
    (double *)calloc(thread_count_size, sizeof(*thread_work));
  UF_long *thread_ptr = NULL;
  UF_long *private_queue = NULL;
  UF_long *pipeline_queue = NULL;
  if (component_private_rows == NULL || component_pipeline_rows == NULL ||
      component_work == NULL || component_thread == NULL ||
      private_components == NULL || thread_counts == NULL ||
      thread_work == NULL) {
    goto cleanup;
  }
  for (UF_long component = 0; component < component_count; ++component) {
    component_thread[component] = -1;
  }

  UF_long private_rows = 0;
  UF_long pipeline_rows = 0;
  UF_long private_component_count = 0;
  UF_long pipeline_component_count = 0;
  for (UF_long local = local_begin; local < local_end; ++local) {
    const unsigned int component = separator->order_component[local];
    const UF_long global_row = separator->global_begin + local;
    if ((UF_long)component >= component_count ||
        global_row < k1 || global_row >= k2) {
      goto cleanup;
    }
    const double row_work =
      kls_first_separator_row_input_work(row_input_count[global_row - k1]);
    if (separator->component_kind[component] == 0u) {
      if (component_private_rows[component] == 0u) {
        private_components[private_component_count].component = component;
        private_component_count++;
      }
      component_private_rows[component]++;
      component_work[component] += row_work;
      private_rows++;
    } else {
      if (component_pipeline_rows[component] == 0u) {
        pipeline_component_count++;
      }
      component_pipeline_rows[component]++;
      pipeline_rows++;
    }
  }
  if (private_rows == 0u || private_component_count == 0u) {
    goto cleanup;
  }

  for (UF_long pos = 0; pos < private_component_count; ++pos) {
    const unsigned int component = private_components[pos].component;
    private_components[pos].rows = component_private_rows[component];
    private_components[pos].work = component_work[component];
  }
  qsort(private_components, (size_t)private_component_count,
        sizeof(*private_components),
        kls_compare_first_separator_component_work_desc);

  for (UF_long pos = 0; pos < private_component_count; ++pos) {
    const unsigned int component = private_components[pos].component;
    int target = 0;
    double target_work = thread_work[0];
    for (int tid = 1; tid < thread_count; ++tid) {
      if (thread_work[tid] < target_work) {
        target = tid;
        target_work = thread_work[tid];
      }
    }
    component_thread[component] = target;
    thread_counts[target] += component_private_rows[component];
    thread_work[target] += component_work[component];
  }

  UF_long nonempty_threads = 0;
  UF_long max_thread_rows = 0;
  double min_thread_work = 0.0;
  double max_thread_work = 0.0;
  for (int tid = 0; tid < thread_count; ++tid) {
    if (thread_counts[tid] > 0u) {
      if (nonempty_threads == 0u || thread_work[tid] < min_thread_work) {
        min_thread_work = thread_work[tid];
      }
      if (thread_work[tid] > max_thread_work) {
        max_thread_work = thread_work[tid];
      }
      nonempty_threads++;
      if (thread_counts[tid] > max_thread_rows) {
        max_thread_rows = thread_counts[tid];
      }
    }
  }
  if (nonempty_threads == 0u) {
    goto cleanup;
  }

  thread_ptr =
    (UF_long *)calloc(thread_count_size + 1u, sizeof(*thread_ptr));
  private_queue =
    (UF_long *)malloc((size_t)private_rows * sizeof(*private_queue));
  pipeline_queue = pipeline_rows > 0u
    ? (UF_long *)malloc((size_t)pipeline_rows * sizeof(*pipeline_queue))
    : NULL;
  if (thread_ptr == NULL || private_queue == NULL ||
      (pipeline_rows > 0u && pipeline_queue == NULL)) {
    goto cleanup;
  }
  for (int tid = 0; tid < thread_count; ++tid) {
    thread_ptr[tid + 1] = thread_ptr[tid] + thread_counts[tid];
    thread_counts[tid] = thread_ptr[tid];
  }
  if (thread_ptr[thread_count] != private_rows) {
    goto cleanup;
  }

  UF_long pipeline_next = 0;
  for (UF_long local = local_begin; local < local_end; ++local) {
    const unsigned int component = separator->order_component[local];
    const UF_long global_row = separator->global_begin + local;
    if (separator->component_kind[component] == 0u) {
      const int target = component_thread[component];
      if (target < 0 || target >= thread_count) {
        goto cleanup;
      }
      const UF_long dst = thread_counts[target]++;
      if (dst >= private_rows) {
        goto cleanup;
      }
      private_queue[dst] = global_row;
    } else {
      if (pipeline_next >= pipeline_rows) {
        goto cleanup;
      }
      pipeline_queue[pipeline_next++] = global_row;
    }
  }
  for (int tid = 0; tid < thread_count; ++tid) {
    if (thread_counts[tid] != thread_ptr[tid + 1]) {
      goto cleanup;
    }
  }
  if (pipeline_next != pipeline_rows) {
    goto cleanup;
  }

  plan->private_rows = private_queue;
  plan->pipeline_rows = pipeline_queue;
  plan->thread_ptr = thread_ptr;
  plan->private_rows_count = private_rows;
  plan->pipeline_rows_count = pipeline_rows;
  plan->private_component_count = private_component_count;
  plan->pipeline_component_count = pipeline_component_count;
  plan->nonempty_threads = nonempty_threads;
  plan->max_thread_rows = max_thread_rows;
  plan->min_thread_work = min_thread_work;
  plan->max_thread_work = max_thread_work;
  private_queue = NULL;
  pipeline_queue = NULL;
  thread_ptr = NULL;
  ok = 1;

cleanup:
  free(component_private_rows);
  free(component_pipeline_rows);
  free(component_work);
  free(component_thread);
  free(private_components);
  free(thread_counts);
  free(thread_work);
  free(thread_ptr);
  free(private_queue);
  free(pipeline_queue);
  return ok;
}

static int kls_build_first_separator_queue_plan_legacy(
  const kls_solver *solver,
  UF_long k1,
  UF_long k2,
  int thread_count,
  kls_first_separator_queue_plan *plan) {
  if (plan != NULL) {
    memset(plan, 0, sizeof(*plan));
  }
  if (solver == NULL || solver->symbolic == NULL || plan == NULL ||
      thread_count <= 1 || k2 <= k1 || k2 > solver->n ||
      solver->col_ptr == NULL || solver->row_idx == NULL ||
      solver->symbolic->P == NULL || solver->symbolic->Q == NULL ||
      !kls_separator_analysis_has_global_range(&solver->separator)) {
    return 0;
  }

  const UF_long n = solver->n;
  const UF_long block_rows = k2 - k1;
  const size_t n_size = (size_t)n;
  if ((UF_long)n_size != n ||
      n_size > SIZE_MAX / sizeof(UF_long) ||
      block_rows > (UF_long)(SIZE_MAX / sizeof(UF_long))) {
    return 0;
  }

  int ok = 0;
  UF_long *symbolic_pos_by_row =
    (UF_long *)malloc(n_size * sizeof(*symbolic_pos_by_row));
  UF_long *row_input_count =
    (UF_long *)calloc((size_t)block_rows, sizeof(*row_input_count));
  if (symbolic_pos_by_row == NULL || row_input_count == NULL) {
    goto cleanup;
  }
  for (UF_long row = 0; row < n; ++row) {
    symbolic_pos_by_row[row] = KLS_KLU_EMPTY;
  }
  for (UF_long pos = 0; pos < n; ++pos) {
    const UF_long row = solver->symbolic->P[pos];
    if (row >= n || symbolic_pos_by_row[row] != KLS_KLU_EMPTY) {
      goto cleanup;
    }
    symbolic_pos_by_row[row] = pos;
  }

  for (UF_long local_col = 0; local_col < block_rows; ++local_col) {
    const UF_long oldcol = solver->symbolic->Q[k1 + local_col];
    if (oldcol >= n ||
        solver->col_ptr[oldcol] > solver->col_ptr[oldcol + 1u] ||
        solver->col_ptr[oldcol + 1u] > solver->nnz) {
      goto cleanup;
    }
    for (UF_long p = solver->col_ptr[oldcol];
         p < solver->col_ptr[oldcol + 1u]; ++p) {
      const UF_long oldrow = solver->row_idx[p];
      if (oldrow >= n) {
        goto cleanup;
      }
      const UF_long row_pos = symbolic_pos_by_row[oldrow];
      if (row_pos >= k1 && row_pos < k2) {
        row_input_count[row_pos - k1]++;
      }
    }
  }

  ok = kls_build_first_separator_queue_plan_legacy_from_counts(
    solver, k1, k2, thread_count, row_input_count, plan);

cleanup:
  free(symbolic_pos_by_row);
  free(row_input_count);
  return ok;
}

static int kls_build_first_separator_queue_plan(
  const kls_solver *solver,
  UF_long k1,
  UF_long k2,
  int thread_count,
  kls_first_separator_queue_plan *plan) {
  if (plan != NULL) {
    memset(plan, 0, sizeof(*plan));
  }
  if (solver == NULL || solver->symbolic == NULL || plan == NULL ||
      thread_count <= 1 || k2 <= k1 || k2 > solver->n ||
      solver->col_ptr == NULL || solver->row_idx == NULL ||
      solver->symbolic->P == NULL || solver->symbolic->Q == NULL ||
      !kls_separator_analysis_has_global_range(&solver->separator) ||
      solver->separator.component_count <= 1u ||
      solver->separator.order_component == NULL ||
      k1 < solver->separator.global_begin ||
      k2 > solver->separator.global_end) {
    return 0;
  }

  if (solver->separator.component_left_child == NULL ||
      solver->separator.component_right_child == NULL) {
    return kls_build_first_separator_queue_plan_legacy(
      solver, k1, k2, thread_count, plan);
  }

  const UF_long n = solver->n;
  const UF_long block_rows = k2 - k1;
  const size_t n_size = (size_t)n;
  if ((UF_long)n_size != n ||
      n_size > SIZE_MAX / sizeof(UF_long) ||
      block_rows > (UF_long)(SIZE_MAX / sizeof(UF_long))) {
    return kls_build_first_separator_queue_plan_legacy(
      solver, k1, k2, thread_count, plan);
  }

  int ok = 0;
  int counted = 0;
  UF_long *symbolic_pos_by_row =
    (UF_long *)malloc(n_size * sizeof(*symbolic_pos_by_row));
  UF_long *row_input_count =
    (UF_long *)calloc((size_t)block_rows, sizeof(*row_input_count));
  if (symbolic_pos_by_row == NULL || row_input_count == NULL) {
    free(symbolic_pos_by_row);
    free(row_input_count);
    return kls_build_first_separator_queue_plan_legacy(
      solver, k1, k2, thread_count, plan);
  }

  for (UF_long row = 0; row < n; ++row) {
    symbolic_pos_by_row[row] = KLS_KLU_EMPTY;
  }
  for (UF_long pos = 0; pos < n; ++pos) {
    const UF_long row = solver->symbolic->P[pos];
    if (row >= n || symbolic_pos_by_row[row] != KLS_KLU_EMPTY) {
      goto cleanup;
    }
    symbolic_pos_by_row[row] = pos;
  }

  for (UF_long local_col = 0; local_col < block_rows; ++local_col) {
    const UF_long oldcol = solver->symbolic->Q[k1 + local_col];
    if (oldcol >= n ||
        solver->col_ptr[oldcol] > solver->col_ptr[oldcol + 1u] ||
        solver->col_ptr[oldcol + 1u] > solver->nnz) {
      goto cleanup;
    }
    for (UF_long p = solver->col_ptr[oldcol];
         p < solver->col_ptr[oldcol + 1u]; ++p) {
      const UF_long oldrow = solver->row_idx[p];
      if (oldrow >= n) {
        goto cleanup;
      }
      const UF_long row_pos = symbolic_pos_by_row[oldrow];
      if (row_pos >= k1 && row_pos < k2) {
        row_input_count[row_pos - k1]++;
      }
    }
  }
  counted = 1;

  ok = kls_build_first_separator_queue_plan_partitioned(
    solver, k1, k2, thread_count, row_input_count, plan);

cleanup:
  if (ok) {
    free(symbolic_pos_by_row);
    free(row_input_count);
    return 1;
  }
  kls_first_separator_queue_plan_free(plan);
  if (counted) {
    ok = kls_build_first_separator_queue_plan_legacy_from_counts(
      solver, k1, k2, thread_count, row_input_count, plan);
    free(symbolic_pos_by_row);
    free(row_input_count);
    return ok;
  }
  free(symbolic_pos_by_row);
  free(row_input_count);
  return kls_build_first_separator_queue_plan_legacy(
    solver, k1, k2, thread_count, plan);
}

static void kls_record_first_separator_queue_plans(
  kls_solver *solver,
  int thread_count) {
  if (solver == NULL || solver->symbolic == NULL ||
      solver->symbolic->R == NULL || thread_count <= 1 ||
      solver->symbolic->nblocks == 0u) {
    return;
  }
  const UF_long nblocks = solver->symbolic->nblocks;
  const int record_only_large_blocks =
    nblocks > 1024u &&
    solver->symbolic->maxblock >=
      KLS_FIRST_FACTOR_DOMINANT_BTF_PIPELINE_MIN_BLOCK;
  /*
   * Queue-plan recording is diagnostic; the executor rebuilds the accepted
   * block plan before running it.  Huge BTF forests such as pre2 can have
   * tens of thousands of fringe blocks, and each plan needs an n-row map.
   * Keep diagnostics focused on the dominant/large executable blocks there.
   */

  UF_long private_components = 0;
  UF_long pipeline_components = 0;
  UF_long private_rows = 0;
  UF_long pipeline_rows = 0;
  UF_long nonempty_threads = 0;
  UF_long max_thread_rows = 0;
  UF_long partitioned_count = 0;
  UF_long split_components = 0;
  double min_thread_work = 0.0;
  double max_thread_work = 0.0;
  int recorded = 0;
  for (UF_long block = 0; block < nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    if (k2 < k1) {
      continue;
    }
    const UF_long block_rows = k2 - k1;
    if (record_only_large_blocks &&
        block_rows != solver->symbolic->maxblock &&
        block_rows < KLS_FIRST_FACTOR_DOMINANT_BTF_PIPELINE_MIN_BLOCK) {
      continue;
    }
    kls_first_separator_queue_plan plan;
    if (!kls_build_first_separator_queue_plan(solver, k1, k2, thread_count,
                                              &plan)) {
      continue;
    }
    recorded = 1;
    private_components += plan.private_component_count;
    pipeline_components += plan.pipeline_component_count;
    private_rows += plan.private_rows_count;
    pipeline_rows += plan.pipeline_rows_count;
    if (plan.nonempty_threads > nonempty_threads) {
      nonempty_threads = plan.nonempty_threads;
    }
    if (plan.max_thread_rows > max_thread_rows) {
      max_thread_rows = plan.max_thread_rows;
    }
    if (plan.partitioned) {
      partitioned_count++;
      if (split_components > UF_long_max - plan.split_component_count) {
        split_components = UF_long_max;
      } else {
        split_components += plan.split_component_count;
      }
    }
    if (plan.nonempty_threads > 0u) {
      if (min_thread_work == 0.0 ||
          (plan.min_thread_work > 0.0 &&
           plan.min_thread_work < min_thread_work)) {
        min_thread_work = plan.min_thread_work;
      }
      if (plan.max_thread_work > max_thread_work) {
        max_thread_work = plan.max_thread_work;
      }
    }
    kls_first_separator_queue_plan_free(&plan);
  }

  if (!recorded || private_rows == 0u) {
    return;
  }
  solver->kls_first_last_separator_queue = 1;
  solver->kls_first_separator_queue_run_count++;
  solver->kls_first_last_separator_queue_private_components =
    private_components;
  solver->kls_first_last_separator_queue_pipeline_components =
    pipeline_components;
  solver->kls_first_last_separator_queue_private_rows = private_rows;
  solver->kls_first_last_separator_queue_pipeline_rows = pipeline_rows;
  solver->kls_first_last_separator_queue_nonempty_threads = nonempty_threads;
  solver->kls_first_last_separator_queue_max_thread_rows = max_thread_rows;
  solver->kls_first_last_separator_queue_min_thread_work = min_thread_work;
  solver->kls_first_last_separator_queue_max_thread_work = max_thread_work;
  solver->kls_first_last_separator_queue_partitioned =
    partitioned_count > 0u ? 1 : 0;
  if (solver->kls_first_separator_queue_partitioned_count >
      UF_long_max - partitioned_count) {
    solver->kls_first_separator_queue_partitioned_count = UF_long_max;
  } else {
    solver->kls_first_separator_queue_partitioned_count +=
      partitioned_count;
  }
  solver->kls_first_last_separator_queue_split_components =
    split_components;
}

typedef struct kls_row_first_block_context {
  kls_solver *solver;
  trilinos_klu_l_numeric *numeric;
  double *numeric_values;
  const UF_long *row_ptr;
  const UF_long *row_cols;
  const UF_long *row_input_pos;
  UF_long *q_order;
  UF_long *col_pos;
  const UF_long *separator_component_last;
  const UF_long *separator_local_order;
  UF_long k1;
  UF_long nk;
  UF_long n;
  UF_long separator_block_base;
  UF_long pivot_col_limit;
  const unsigned char *pivot_col_mask;
  int scaled;
  int scale;
  int use_separator_for_block;
  int use_pivot_col_limit;
  int strict_separator_pivot_scope;
  double tol;
} kls_row_first_block_context;

typedef struct kls_row_first_supernode_panel_cache {
  UF_long *panel_id_by_row;
  UF_long *start;
  UF_long *width;
  UF_long *tail_begin;
  UF_long *tail_len;
  UF_long *dense_begin;
  UF_long *tail_value_begin;
  unsigned char *active;
  UF_long *tail_cols;
  double *dense_values;
  double *tail_values;
  UF_long row_count;
  UF_long panel_count;
  UF_long panel_capacity;
  UF_long tail_col_count;
  UF_long tail_col_capacity;
  UF_long dense_value_count;
  UF_long dense_value_capacity;
  UF_long tail_value_count;
  UF_long tail_value_capacity;
  int enabled;
} kls_row_first_supernode_panel_cache;

typedef struct kls_row_first_workspace {
  double *x;
  double *supernode_workspace;
  unsigned int *mark;
  UF_long *pattern;
  UF_long *dep_heap;
  UF_long *u_row_ptr;
  UF_long *u_row_end;
  UF_long supernode_workspace_capacity;
  kls_row_first_supernode_panel_cache *supernode_panel_cache;
  unsigned char *supernode_fit_status;
} kls_row_first_workspace;

typedef struct kls_row_first_partial_row {
  UF_long row;
  unsigned int generation;
  UF_long pattern_count;
  UF_long dep_heap_size;
  double pivot;
} kls_row_first_partial_row;

typedef struct kls_row_first_row_stats {
  UF_long dynamic_column_pivots;
  UF_long separator_dynamic_column_pivots;
  UF_long separator_extent_dynamic_column_pivots;
  UF_long separator_dynamic_column_fallbacks;
  UF_long separator_dynamic_column_rejects;
  UF_long supernode_update_groups;
  UF_long supernode_update_rows;
  UF_long supernode_panel_update_groups;
  UF_long supernode_panel_update_rows;
  UF_long pipeline_prefix_panel_rebuilds;
  UF_long pipeline_prefix_panel_rebuild_rows;
  UF_long active_rank_pivot_resets;
  UF_long active_rank_pivot_reset_rows;
  UF_long active_rank_pivot_panel_rebuilds;
  UF_long active_rank_pivot_panel_rebuild_rows;
  UF_long panel_cache_builds;
  UF_long panel_cache_build_panels;
  UF_long panel_cache_build_entries;
  UF_long panel_cache_appends;
  UF_long panel_cache_append_panels;
  UF_long panel_cache_append_entries;
} kls_row_first_row_stats;

typedef enum kls_row_first_factor_row_failure {
  KLS_ROW_FIRST_FACTOR_ROW_FAIL_NONE = 0,
  KLS_ROW_FIRST_FACTOR_ROW_FAIL_SETUP,
  KLS_ROW_FIRST_FACTOR_ROW_FAIL_INPUT,
  KLS_ROW_FIRST_FACTOR_ROW_FAIL_VALUE,
  KLS_ROW_FIRST_FACTOR_ROW_FAIL_DEPENDENCY,
  KLS_ROW_FIRST_FACTOR_ROW_FAIL_UPDATE,
  KLS_ROW_FIRST_FACTOR_ROW_FAIL_PIVOT_SCOPE,
  KLS_ROW_FIRST_FACTOR_ROW_FAIL_PIVOT_REJECT,
  KLS_ROW_FIRST_FACTOR_ROW_FAIL_PIVOT_OWNER,
  KLS_ROW_FIRST_FACTOR_ROW_FAIL_PIVOT_EXCHANGE,
  KLS_ROW_FIRST_FACTOR_ROW_FAIL_PIVOT_VERIFY,
  KLS_ROW_FIRST_FACTOR_ROW_FAIL_APPEND
} kls_row_first_factor_row_failure;

typedef struct kls_row_first_factor_row_failure_info {
  kls_row_first_factor_row_failure reason;
  UF_long row;
  UF_long related_col;
  int owner;
  int actual_owner;
} kls_row_first_factor_row_failure_info;

static int kls_row_first_pivot_failure_col_state(
  const kls_row_first_block_context *ctx,
  const int *row_owner,
  UF_long col) {
  if (col == KLS_KLU_EMPTY || ctx == NULL || col >= ctx->nk) {
    return -1;
  }
  if (row_owner != NULL) {
    return row_owner[col];
  }
  if (ctx->pivot_col_mask != NULL) {
    return ctx->pivot_col_mask[col] ? 1 : 0;
  }
  return -1;
}

static const char *kls_row_first_factor_row_failure_name(
  kls_row_first_factor_row_failure reason) {
  switch (reason) {
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_NONE:
      return "none";
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_SETUP:
      return "setup";
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_INPUT:
      return "input";
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_VALUE:
      return "value";
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_DEPENDENCY:
      return "dependency";
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_UPDATE:
      return "update";
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_PIVOT_SCOPE:
      return "pivot-scope";
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_PIVOT_REJECT:
      return "pivot-reject";
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_PIVOT_OWNER:
      return "pivot-owner";
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_PIVOT_EXCHANGE:
      return "pivot-exchange";
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_PIVOT_VERIFY:
      return "pivot-verify";
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_APPEND:
      return "append";
  }
  return "unknown";
}

static void kls_row_first_factor_row_failure_set(
  kls_row_first_factor_row_failure_info *failure,
  kls_row_first_factor_row_failure reason,
  UF_long row,
  UF_long related_col,
  int owner,
  int actual_owner) {
  if (failure == NULL) {
    return;
  }
  failure->reason = reason;
  failure->row = row;
  failure->related_col = related_col;
  failure->owner = owner;
  failure->actual_owner = actual_owner;
}

static void kls_row_first_factor_row_failure_clear(
  kls_row_first_factor_row_failure_info *failure,
  UF_long row,
  int owner) {
  kls_row_first_factor_row_failure_set(
    failure, KLS_ROW_FIRST_FACTOR_ROW_FAIL_NONE, row, KLS_KLU_EMPTY,
    owner, -1);
}

static void kls_row_first_stats_add(UF_long *target, UF_long value) {
  if (target == NULL || value == 0u) {
    return;
  }
  *target = *target > UF_long_max - value ? UF_long_max : *target + value;
}

static UF_long kls_row_first_panel_cache_stored_entries(UF_long width,
                                                        UF_long tail_len) {
  UF_long entries = 0;
  if (width > 0u && width <= UF_long_max / width) {
    entries = width * width;
  } else if (width > 0u) {
    return UF_long_max;
  }
  const UF_long tail_values =
    tail_len > 0u && width > UF_long_max / tail_len
      ? UF_long_max
      : width * tail_len;
  if (entries > UF_long_max - tail_values) {
    return UF_long_max;
  }
  entries += tail_values;
  return entries > UF_long_max - tail_len ? UF_long_max : entries + tail_len;
}

static UF_long kls_row_first_panel_cache_current_entries(
  const kls_row_first_supernode_panel_cache *cache) {
  if (cache == NULL) {
    return 0u;
  }
  UF_long entries = cache->dense_value_count;
  if (entries > UF_long_max - cache->tail_value_count) {
    return UF_long_max;
  }
  entries += cache->tail_value_count;
  if (entries > UF_long_max - cache->tail_col_count) {
    return UF_long_max;
  }
  return entries + cache->tail_col_count;
}

static int kls_row_first_panel_cache_allows_missed_append(
  const kls_row_first_supernode_panel_cache *cache,
  UF_long width,
  UF_long tail_len,
  UF_long max_stored_entries) {
  if (max_stored_entries == UF_long_max) {
    return 1;
  }
  const UF_long panel_entries =
    kls_row_first_panel_cache_stored_entries(width, tail_len);
  if (panel_entries == UF_long_max || panel_entries > max_stored_entries) {
    return 0;
  }
  const UF_long current_entries =
    kls_row_first_panel_cache_current_entries(cache);
  return current_entries <= max_stored_entries &&
         panel_entries <= max_stored_entries - current_entries;
}

static void kls_row_first_record_panel_cache_build(
  kls_row_first_row_stats *stats,
  UF_long panels,
  UF_long dense_values,
  UF_long tail_values,
  UF_long tail_cols) {
  if (stats == NULL) {
    return;
  }
  kls_row_first_stats_add(&stats->panel_cache_builds, 1u);
  kls_row_first_stats_add(&stats->panel_cache_build_panels, panels);
  kls_row_first_stats_add(&stats->panel_cache_build_entries, dense_values);
  kls_row_first_stats_add(&stats->panel_cache_build_entries, tail_values);
  kls_row_first_stats_add(&stats->panel_cache_build_entries, tail_cols);
}

static void kls_row_first_record_panel_cache_append(
  kls_row_first_row_stats *stats,
  UF_long width,
  UF_long tail_len) {
  if (stats == NULL) {
    return;
  }
  kls_row_first_stats_add(&stats->panel_cache_appends, 1u);
  kls_row_first_stats_add(&stats->panel_cache_append_panels, 1u);
  kls_row_first_stats_add(
    &stats->panel_cache_append_entries,
    kls_row_first_panel_cache_stored_entries(width, tail_len));
}

static void kls_row_first_accumulate_panel_cache_stats(
  kls_row_first_row_stats *dst,
  const kls_row_first_row_stats *src) {
  if (dst == NULL || src == NULL) {
    return;
  }
  kls_row_first_stats_add(&dst->panel_cache_builds,
                          src->panel_cache_builds);
  kls_row_first_stats_add(&dst->panel_cache_build_panels,
                          src->panel_cache_build_panels);
  kls_row_first_stats_add(&dst->panel_cache_build_entries,
                          src->panel_cache_build_entries);
  kls_row_first_stats_add(&dst->panel_cache_appends,
                          src->panel_cache_appends);
  kls_row_first_stats_add(&dst->panel_cache_append_panels,
                          src->panel_cache_append_panels);
  kls_row_first_stats_add(&dst->panel_cache_append_entries,
                          src->panel_cache_append_entries);
  kls_row_first_stats_add(&dst->active_rank_pivot_resets,
                          src->active_rank_pivot_resets);
  kls_row_first_stats_add(&dst->active_rank_pivot_reset_rows,
                          src->active_rank_pivot_reset_rows);
  kls_row_first_stats_add(&dst->active_rank_pivot_panel_rebuilds,
                          src->active_rank_pivot_panel_rebuilds);
  kls_row_first_stats_add(&dst->active_rank_pivot_panel_rebuild_rows,
                          src->active_rank_pivot_panel_rebuild_rows);
}

typedef struct kls_row_first_pivot_choice {
  UF_long global_col;
  double global_abs;
  UF_long limited_col;
  double limited_abs;
  UF_long scoped_col;
  double scoped_abs;
  int scoped_col_exact;
  int scope_enabled;
  int limit_enabled;
} kls_row_first_pivot_choice;

static int kls_row_first_partial_apply_one_dep(
  const kls_row_first_block_context *ctx,
  kls_row_first_workspace *workspace,
  kls_row_first_entries *local_l_entries,
  const kls_row_first_entries *published_u_entries,
  const double *udiag_values,
  const UF_long *u_row_ptr,
  const UF_long *u_row_end,
  kls_row_first_partial_row *state,
  UF_long dep);

static int kls_row_first_partial_apply_supernode_run(
  const kls_row_first_block_context *ctx,
  kls_row_first_workspace *workspace,
  kls_row_first_entries *local_l_entries,
  const kls_row_first_entries *published_u_entries,
  const double *udiag_values,
  const UF_long *u_row_ptr,
  const UF_long *u_row_end,
  kls_row_first_partial_row *state,
  UF_long dep_begin,
  UF_long dep_end,
  UF_long *run_rows_out,
  int *used_panel_out,
  kls_row_first_row_stats *stats);

static void kls_row_first_supernode_panel_cache_free(
  kls_row_first_supernode_panel_cache *cache);

static int kls_row_first_supernode_panel_cache_init_rows(
  kls_row_first_supernode_panel_cache *cache,
  UF_long nk);

static void kls_row_first_supernode_panel_cache_reset_active(
  kls_row_first_supernode_panel_cache *cache,
  kls_row_first_workspace *workspace,
  UF_long nk,
  UF_long *open_start_io,
  UF_long *open_end_io);

static void kls_row_first_supernode_panel_cache_publish_completed(
  kls_row_first_supernode_panel_cache *cache,
  const kls_row_first_entries *u_entries,
  const UF_long *u_row_ptr,
  const UF_long *u_row_end,
  const double *udiag_values,
  UF_long nk,
  kls_row_first_row_stats *stats,
  UF_long *open_start_io,
  UF_long *open_end_io,
  UF_long row);

static UF_long kls_row_first_owned_supernode_end(
  const kls_row_first_entries *u_entries,
  const UF_long *u_row_ptr,
  const UF_long *u_row_end,
  const unsigned char *row_done,
  const int *row_owner,
  int owner,
  UF_long nk,
  UF_long current_row,
  unsigned char *fit_status,
  UF_long dep);

static int kls_row_first_separator_position(
  const kls_row_first_block_context *ctx,
  UF_long local_col,
  UF_long *separator_pos_out) {
  if (separator_pos_out != NULL) {
    *separator_pos_out = KLS_KLU_EMPTY;
  }
  if (ctx == NULL || separator_pos_out == NULL ||
      !ctx->use_separator_for_block || ctx->solver == NULL ||
      ctx->solver->separator.order_component == NULL ||
      local_col >= ctx->nk) {
    return 0;
  }
  const UF_long separator_local =
    ctx->separator_local_order != NULL
      ? ctx->separator_local_order[local_col] : local_col;
  if (separator_local >= ctx->nk ||
      separator_local >
        UF_long_max - ctx->separator_block_base) {
    return 0;
  }
  const UF_long separator_pos =
    ctx->separator_block_base + separator_local;
  if (separator_pos >= ctx->solver->separator.n) {
    return 0;
  }
  *separator_pos_out = separator_pos;
  return 1;
}

static int kls_row_first_pivot_scope(
  const kls_row_first_block_context *ctx,
  UF_long row,
  unsigned int *component_out,
  UF_long *component_last_out) {
  if (component_out != NULL) {
    *component_out = UINT_MAX;
  }
  if (component_last_out != NULL) {
    *component_last_out = KLS_KLU_EMPTY;
  }
  if (ctx == NULL || component_out == NULL || component_last_out == NULL ||
      !ctx->use_separator_for_block || ctx->solver == NULL ||
      ctx->separator_component_last == NULL ||
      ctx->solver->separator.order_component == NULL) {
    return 0;
  }
  UF_long separator_pos = KLS_KLU_EMPTY;
  if (!kls_row_first_separator_position(ctx, row, &separator_pos)) {
    return 0;
  }
  const unsigned int component =
    ctx->solver->separator.order_component[separator_pos];
  if ((UF_long)component >= ctx->solver->separator.component_count) {
    return 0;
  }
  const UF_long component_last = ctx->separator_component_last[component];
  if (component_last == KLS_KLU_EMPTY ||
      component_last >= ctx->solver->separator.n ||
      component_last < separator_pos) {
    return 0;
  }
  *component_out = component;
  *component_last_out = component_last;
  return 1;
}

static int kls_row_first_validate_tail_separator_scope(
  const kls_row_first_block_context *ctx,
  const UF_long *row_order,
  UF_long row_count,
  UF_long nk,
  const unsigned char *active_mask,
  UF_long *scoped_rows_out) {
  KLS_SET_OPTIONAL_OUTPUT(scoped_rows_out, 0);
  if (ctx == NULL || row_order == NULL || active_mask == NULL ||
      scoped_rows_out == NULL || nk > ctx->nk ||
      !ctx->use_separator_for_block) {
    return 0;
  }
  UF_long scoped_rows = 0;
  for (UF_long pos = 0; pos < row_count; ++pos) {
    const UF_long row = row_order[pos];
    unsigned int component = UINT_MAX;
    UF_long component_last = KLS_KLU_EMPTY;
    if (row >= nk || !active_mask[row] ||
        !kls_row_first_pivot_scope(ctx, row, &component,
                                   &component_last)) {
      return 0;
    }
    (void)component;
    (void)component_last;
    scoped_rows++;
  }
  *scoped_rows_out = scoped_rows;
  return 1;
}

static int kls_row_first_collect_pivot_choice(
  const kls_row_first_block_context *ctx,
  const UF_long *pattern,
  UF_long pattern_count,
  const double *x,
  UF_long row,
  const int *row_owner,
  int owner,
  kls_row_first_pivot_choice *choice) {
  if (choice != NULL) {
    choice->global_col = KLS_KLU_EMPTY;
    choice->global_abs = 0.0;
    choice->limited_col = KLS_KLU_EMPTY;
    choice->limited_abs = 0.0;
    choice->scoped_col = KLS_KLU_EMPTY;
    choice->scoped_abs = 0.0;
    choice->scoped_col_exact = 0;
    choice->scope_enabled = 0;
    choice->limit_enabled = 0;
  }
  if (ctx == NULL || pattern == NULL || x == NULL || choice == NULL ||
      row >= ctx->nk) {
    return 0;
  }

  unsigned int pivot_component = UINT_MAX;
  UF_long pivot_component_last = KLS_KLU_EMPTY;
  if (ctx->use_separator_for_block) {
    if (!kls_row_first_pivot_scope(
          ctx, row, &pivot_component, &pivot_component_last)) {
      return 0;
    }
    choice->scope_enabled = 1;
  }
  const int limit_enabled =
    (ctx->use_pivot_col_limit && ctx->pivot_col_limit <= ctx->nk) ||
    ctx->pivot_col_mask != NULL;
  if (limit_enabled) {
    choice->limit_enabled = 1;
  }

  double scoped_best_abs = -1.0;
  double limited_best_abs = -1.0;
  const int owner_restricted = row_owner != NULL && owner >= 0;
  for (UF_long p = 0; p < pattern_count; ++p) {
    const UF_long col = pattern[p];
    if (col <= row) {
      continue;
    }
    if (owner_restricted && row_owner[col] != owner) {
      continue;
    }
    const double abs_value = fabs(x[col]);
    if (abs_value > choice->global_abs) {
      choice->global_abs = abs_value;
      choice->global_col = col;
    }
    int within_limit = 1;
    if (ctx->use_pivot_col_limit && ctx->pivot_col_limit <= ctx->nk &&
        col >= ctx->pivot_col_limit) {
      within_limit = 0;
    }
    if (ctx->pivot_col_mask != NULL && !ctx->pivot_col_mask[col]) {
      within_limit = 0;
    }
    if (within_limit && abs_value > limited_best_abs) {
      limited_best_abs = abs_value;
      choice->limited_abs = abs_value;
      choice->limited_col = col;
    }
    if (choice->scope_enabled) {
      UF_long separator_pos = KLS_KLU_EMPTY;
      if (!kls_row_first_separator_position(ctx, col, &separator_pos)) {
        return 0;
      }
      const unsigned int col_component =
        ctx->solver->separator.order_component[separator_pos];
      const int exact_component = col_component == pivot_component;
      if (within_limit &&
          separator_pos <= pivot_component_last &&
          (!ctx->strict_separator_pivot_scope || exact_component) &&
          abs_value > scoped_best_abs) {
        scoped_best_abs = abs_value;
        choice->scoped_abs = abs_value;
        choice->scoped_col = col;
        choice->scoped_col_exact = exact_component;
      }
    }
  }
  return 1;
}

static int kls_row_first_select_pivot_exchange(
  const kls_row_first_pivot_choice *choice,
  double pivot,
  double tol,
  int *pivot_needed_out,
  UF_long *selected_col_out,
  int *selected_separator_exact_out,
  int *selected_separator_extent_out) {
  KLS_SET_OPTIONAL_OUTPUT(pivot_needed_out, 0);
  if (selected_col_out != NULL) {
    *selected_col_out = KLS_KLU_EMPTY;
  }
  KLS_SET_OPTIONAL_OUTPUT(selected_separator_exact_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(selected_separator_extent_out, 0);
  if (choice == NULL || pivot_needed_out == NULL ||
      selected_col_out == NULL) {
    return 0;
  }

  const double check_abs =
    choice->limit_enabled ? choice->global_abs
                          : (choice->scope_enabled ? choice->scoped_abs
                                                   : choice->global_abs);
  const double candidate_abs =
    choice->scope_enabled ? choice->scoped_abs
                          : (choice->limit_enabled ? choice->limited_abs
                                                   : choice->global_abs);
  const UF_long best_col =
    choice->scope_enabled ? choice->scoped_col
                          : (choice->limit_enabled ? choice->limited_col
                                                   : choice->global_col);
  if (pivot != 0.0 &&
      (check_abs <= 0.0 || fabs(pivot) >= tol * check_abs)) {
    return 1;
  }
  if (best_col == KLS_KLU_EMPTY || candidate_abs <= 0.0 ||
      (choice->limit_enabled && candidate_abs < tol * check_abs)) {
    return 0;
  }

  *pivot_needed_out = 1;
  *selected_col_out = best_col;
  if (choice->scope_enabled) {
    if (selected_separator_exact_out != NULL) {
      *selected_separator_exact_out = choice->scoped_col_exact;
    }
    if (selected_separator_extent_out != NULL) {
      *selected_separator_extent_out = !choice->scoped_col_exact;
    }
  }
  return 1;
}

static int kls_row_first_factor_one_row(
  const kls_row_first_block_context *ctx,
  kls_row_first_workspace *workspace,
  kls_row_first_entries *l_entries,
  kls_row_first_entries *u_entries,
  double *udiag_values,
  UF_long i,
  const int *row_owner,
  int owner,
  unsigned char *row_done,
  kls_row_first_row_stats *stats,
  kls_row_first_factor_row_failure_info *failure) {
  kls_row_first_factor_row_failure_clear(failure, i, owner);
  if (ctx == NULL || workspace == NULL || l_entries == NULL ||
      u_entries == NULL || udiag_values == NULL ||
      workspace->x == NULL || workspace->mark == NULL ||
      workspace->pattern == NULL || workspace->dep_heap == NULL ||
      workspace->u_row_ptr == NULL || workspace->u_row_end == NULL ||
      ctx->solver == NULL || ctx->numeric == NULL ||
      ctx->numeric_values == NULL || ctx->row_ptr == NULL ||
      ctx->row_cols == NULL || ctx->row_input_pos == NULL ||
      ctx->q_order == NULL || ctx->col_pos == NULL || i >= ctx->nk) {
    kls_row_first_factor_row_failure_set(
      failure, KLS_ROW_FIRST_FACTOR_ROW_FAIL_SETUP, i, KLS_KLU_EMPTY,
      owner, -1);
    return 0;
  }
  const unsigned int generation = (unsigned int)(i + 1u);
  if (generation == 0u) {
    kls_row_first_factor_row_failure_set(
      failure, KLS_ROW_FIRST_FACTOR_ROW_FAIL_SETUP, i, KLS_KLU_EMPTY,
      owner, -1);
    return 0;
  }

  kls_solver *solver = ctx->solver;
  double *x = workspace->x;
  unsigned int *mark = workspace->mark;
  UF_long *pattern = workspace->pattern;
  UF_long *dep_heap = workspace->dep_heap;
  UF_long pattern_count = 0;
  UF_long dep_heap_size = 0;
  for (UF_long p = ctx->row_ptr[i]; p < ctx->row_ptr[i + 1u]; ++p) {
    const UF_long oldcol = ctx->row_cols[p];
    if (oldcol >= ctx->n || ctx->col_pos[oldcol] >= ctx->nk ||
        ctx->row_input_pos[p] >= solver->nnz) {
      kls_row_first_factor_row_failure_set(
        failure, KLS_ROW_FIRST_FACTOR_ROW_FAIL_INPUT, i, oldcol,
        owner, -1);
      return 0;
    }
    const UF_long col = ctx->col_pos[oldcol];
    if (mark[col] != generation) {
      mark[col] = generation;
      pattern[pattern_count++] = col;
      if (col < i) {
        kls_row_first_heap_push(dep_heap, &dep_heap_size, col);
      }
      if (!kls_refactor_input_value(
            solver, ctx->numeric_values,
            ctx->scaled ? ctx->numeric->Rs : NULL,
            ctx->scale, ctx->row_input_pos[p], x + col)) {
        kls_row_first_factor_row_failure_set(
          failure, KLS_ROW_FIRST_FACTOR_ROW_FAIL_VALUE, i, col,
          owner, -1);
        return 0;
      }
    } else {
      double value = 0.0;
      if (!kls_refactor_input_value(
            solver, ctx->numeric_values,
            ctx->scaled ? ctx->numeric->Rs : NULL,
            ctx->scale, ctx->row_input_pos[p], &value)) {
        kls_row_first_factor_row_failure_set(
          failure, KLS_ROW_FIRST_FACTOR_ROW_FAIL_VALUE, i, col,
          owner, -1);
        return 0;
      }
      x[col] += value;
    }
  }

  workspace->u_row_ptr[i] = u_entries->count;
  kls_row_first_partial_row state;
  memset(&state, 0, sizeof(state));
  state.row = i;
  state.generation = generation;
  state.pattern_count = pattern_count;
  state.dep_heap_size = dep_heap_size;
  while (state.dep_heap_size > 0u) {
    const UF_long dep =
      kls_row_first_heap_pop(dep_heap, &state.dep_heap_size);
    if (dep >= i || mark[dep] != generation ||
        (row_owner != NULL && row_owner[dep] != owner) ||
        (row_done != NULL && !row_done[dep])) {
      kls_row_first_factor_row_failure_set(
        failure, KLS_ROW_FIRST_FACTOR_ROW_FAIL_DEPENDENCY, i, dep,
        owner, row_owner != NULL ? row_owner[dep] : -1);
      return 0;
    }
    const UF_long run_end =
      kls_row_first_owned_supernode_end(
        u_entries, workspace->u_row_ptr, workspace->u_row_end, row_done,
        row_owner, owner, ctx->nk, i, workspace->supernode_fit_status, dep);
    UF_long run_rows = 1;
    int used_panel = 0;
    if (run_end > dep) {
      if (!kls_row_first_partial_apply_supernode_run(
            ctx, workspace, l_entries, u_entries, udiag_values,
            workspace->u_row_ptr, workspace->u_row_end, &state, dep,
            run_end, &run_rows, &used_panel, stats)) {
        kls_row_first_factor_row_failure_set(
          failure, KLS_ROW_FIRST_FACTOR_ROW_FAIL_UPDATE, i, dep,
          owner, row_owner != NULL ? row_owner[dep] : -1);
        return 0;
      }
    } else if (!kls_row_first_partial_apply_one_dep(
                 ctx, workspace, l_entries, u_entries, udiag_values,
                 workspace->u_row_ptr, workspace->u_row_end, &state,
                 dep)) {
      kls_row_first_factor_row_failure_set(
        failure, KLS_ROW_FIRST_FACTOR_ROW_FAIL_UPDATE, i, dep,
        owner, row_owner != NULL ? row_owner[dep] : -1);
      return 0;
    }
    if (run_rows > 1u && stats != NULL) {
      stats->supernode_update_groups++;
      stats->supernode_update_rows += run_rows;
      if (used_panel) {
        stats->supernode_panel_update_groups++;
        stats->supernode_panel_update_rows += run_rows;
      }
    }
  }
  pattern_count = state.pattern_count;

  double pivot = mark[i] == generation ? x[i] : 0.0;
  UF_long u_count = 0;
  for (UF_long p = 0; p < pattern_count; ++p) {
    const UF_long col = pattern[p];
    if (col > i) {
      dep_heap[u_count++] = col;
    }
  }
  kls_row_first_pivot_choice pivot_choice;
  if (!kls_row_first_collect_pivot_choice(
        ctx, pattern, pattern_count, x, i, row_owner, owner,
        &pivot_choice)) {
    kls_row_first_factor_row_failure_set(
      failure, KLS_ROW_FIRST_FACTOR_ROW_FAIL_PIVOT_SCOPE, i,
      KLS_KLU_EMPTY, owner, -1);
    return 0;
  }
  int pivot_needed = 0;
  UF_long selected_col = KLS_KLU_EMPTY;
  int selected_separator_exact = 0;
  int selected_separator_extent = 0;
  if (!kls_row_first_select_pivot_exchange(
        &pivot_choice, pivot, ctx->tol, &pivot_needed, &selected_col,
        &selected_separator_exact, &selected_separator_extent)) {
    if (ctx->use_separator_for_block && stats != NULL) {
      stats->separator_dynamic_column_fallbacks++;
      stats->separator_dynamic_column_rejects++;
    }
    kls_row_first_factor_row_failure_set(
      failure, KLS_ROW_FIRST_FACTOR_ROW_FAIL_PIVOT_REJECT, i,
      pivot_choice.global_col, owner,
      kls_row_first_pivot_failure_col_state(
        ctx, row_owner, pivot_choice.global_col));
    return 0;
  }
  if (pivot_needed) {
    if (selected_col == KLS_KLU_EMPTY) {
      kls_row_first_factor_row_failure_set(
        failure, KLS_ROW_FIRST_FACTOR_ROW_FAIL_PIVOT_REJECT, i,
        pivot_choice.global_col, owner,
        kls_row_first_pivot_failure_col_state(
          ctx, row_owner, pivot_choice.global_col));
      return 0;
    }
    if (row_owner != NULL && row_owner[selected_col] != owner) {
      kls_row_first_factor_row_failure_set(
        failure, KLS_ROW_FIRST_FACTOR_ROW_FAIL_PIVOT_OWNER, i,
        selected_col, owner, row_owner[selected_col]);
      return 0;
    }
    if (ctx->use_separator_for_block) {
      if (selected_separator_exact) {
        if (stats != NULL) {
          stats->separator_dynamic_column_pivots++;
        }
      } else if (selected_separator_extent) {
        if (stats != NULL) {
          stats->separator_extent_dynamic_column_pivots++;
        }
      } else {
        if (stats != NULL) {
          stats->separator_dynamic_column_fallbacks++;
          stats->separator_dynamic_column_rejects++;
        }
        kls_row_first_factor_row_failure_set(
          failure, KLS_ROW_FIRST_FACTOR_ROW_FAIL_PIVOT_REJECT, i,
          selected_col, owner, row_owner != NULL ? row_owner[selected_col] : -1);
        return 0;
      }
    }
    if (!kls_row_first_exchange_columns(
          u_entries, ctx->q_order, ctx->col_pos, i, selected_col,
          ctx->nk, ctx->n)) {
      kls_row_first_factor_row_failure_set(
        failure, KLS_ROW_FIRST_FACTOR_ROW_FAIL_PIVOT_EXCHANGE, i,
        selected_col, owner, row_owner != NULL ? row_owner[selected_col] : -1);
      return 0;
    }
    const double swapped_pivot = x[selected_col];
    x[selected_col] = pivot;
    x[i] = swapped_pivot;
    mark[i] = generation;
    if (stats != NULL) {
      stats->dynamic_column_pivots++;
    }
    pivot = x[i];
    u_count = 0;
    for (UF_long p = 0; p < pattern_count; ++p) {
      const UF_long col = pattern[p];
      if (col > i) {
        dep_heap[u_count++] = col;
      }
    }
    if (!kls_row_first_collect_pivot_choice(
          ctx, pattern, pattern_count, x, i, row_owner, owner,
          &pivot_choice) ||
        !kls_row_first_select_pivot_exchange(
          &pivot_choice, pivot, ctx->tol, &pivot_needed, &selected_col,
          NULL, NULL) ||
        pivot_needed) {
      kls_row_first_factor_row_failure_set(
        failure, KLS_ROW_FIRST_FACTOR_ROW_FAIL_PIVOT_VERIFY, i,
        selected_col, owner,
        row_owner != NULL && selected_col < ctx->nk ? row_owner[selected_col] : -1);
      return 0;
    }
  }
  if (u_count > 1u) {
    qsort(dep_heap, (size_t)u_count, sizeof(*dep_heap),
          kls_compare_uf_long);
  }
  udiag_values[i] = pivot;
  for (UF_long p = 0; p < u_count; ++p) {
    const UF_long col = dep_heap[p];
    if (!kls_row_first_entries_append(u_entries, i, col, x[col])) {
      kls_row_first_factor_row_failure_set(
        failure, KLS_ROW_FIRST_FACTOR_ROW_FAIL_APPEND, i, col,
        owner, -1);
      return 0;
    }
  }
  workspace->u_row_end[i] = u_entries->count;
  for (UF_long p = 0; p < pattern_count; ++p) {
    x[pattern[p]] = 0.0;
  }
  if (mark[i] == generation) {
    x[i] = 0.0;
  }
  if (row_done != NULL) {
    row_done[i] = 1u;
  }
  return 1;
}

typedef struct kls_row_first_private_worker {
  const kls_row_first_block_context *ctx;
  const kls_first_separator_queue_plan *plan;
  const int *row_owner;
  kls_row_first_workspace workspace;
  kls_row_first_entries l_entries;
  kls_row_first_entries u_entries;
  kls_row_first_row_stats stats;
  kls_row_first_factor_row_failure_info failure;
  double *udiag_values;
  UF_long begin;
  UF_long end;
  UF_long completed_rows;
  UF_long failed_pos;
  UF_long failed_row;
  UF_long failed_global_row;
  int tid;
  int ok;
} kls_row_first_private_worker;

static void kls_row_first_workspace_free(
  kls_row_first_workspace *workspace) {
  if (workspace == NULL) {
    return;
  }
  free(workspace->x);
  free(workspace->supernode_workspace);
  free(workspace->mark);
  free(workspace->pattern);
  free(workspace->dep_heap);
  free(workspace->u_row_ptr);
  free(workspace->u_row_end);
  free(workspace->supernode_fit_status);
  memset(workspace, 0, sizeof(*workspace));
}

static int kls_row_first_workspace_enable_supernode_fit_cache(
  kls_row_first_workspace *workspace,
  UF_long nk) {
  if (workspace == NULL ||
      nk > (UF_long)(SIZE_MAX / sizeof(*workspace->supernode_fit_status))) {
    return 0;
  }
  free(workspace->supernode_fit_status);
  workspace->supernode_fit_status =
    (unsigned char *)calloc((size_t)nk,
                            sizeof(*workspace->supernode_fit_status));
  return workspace->supernode_fit_status != NULL || nk == 0u;
}

static int kls_row_first_workspace_reserve_supernode(
  kls_row_first_workspace *workspace,
  UF_long entry_count) {
  if (entry_count == 0u) {
    return 1;
  }
  if (workspace == NULL ||
      entry_count > (UF_long)(SIZE_MAX / sizeof(*workspace->supernode_workspace))) {
    return 0;
  }
  if (workspace->supernode_workspace_capacity >= entry_count) {
    return 1;
  }
  double *grown =
    (double *)realloc(workspace->supernode_workspace,
                      (size_t)entry_count *
                        sizeof(*workspace->supernode_workspace));
  if (grown == NULL) {
    return 0;
  }
  /* consumers accumulate into this buffer assuming zeroed lanes; the
     realloc'd tail is recycled heap whose stale bytes (another
     solver's freed values) entered sums silently - the task #21
     wrong-solve vehicle.  Zero the grown region. */
  memset(grown + workspace->supernode_workspace_capacity, 0,
         (size_t)(entry_count - workspace->supernode_workspace_capacity) *
           sizeof(*grown));
  workspace->supernode_workspace = grown;
  workspace->supernode_workspace_capacity = entry_count;
  return 1;
}

static int kls_row_first_workspace_init(
  kls_row_first_workspace *workspace,
  UF_long nk) {
  if (workspace == NULL ||
      nk > (UF_long)(SIZE_MAX / sizeof(*workspace->x)) ||
      nk > (UF_long)(SIZE_MAX / sizeof(*workspace->mark)) ||
      nk > (UF_long)(SIZE_MAX / sizeof(*workspace->pattern))) {
    return 0;
  }
  /* re-init used to memset the struct BEFORE freeing: every second
     init leaked the old buffers as live memory, and any retained view
     of them read stale-but-live data (the cross-solver wrong-solve
     family, task #21).  Free first - retained views then dangle to
     freed memory, which ASAN reports instead of silently mis-solving. */
  kls_row_first_workspace_free(workspace);
  workspace->x = (double *)calloc((size_t)nk, sizeof(*workspace->x));
  workspace->mark =
    (unsigned int *)calloc((size_t)nk, sizeof(*workspace->mark));
  workspace->pattern =
    (UF_long *)malloc((size_t)nk * sizeof(*workspace->pattern));
  workspace->dep_heap =
    (UF_long *)malloc((size_t)nk * sizeof(*workspace->dep_heap));
  workspace->u_row_ptr =
    (UF_long *)calloc((size_t)nk + 1u, sizeof(*workspace->u_row_ptr));
  workspace->u_row_end =
    (UF_long *)calloc((size_t)nk, sizeof(*workspace->u_row_end));
  if (workspace->x == NULL || workspace->mark == NULL ||
      workspace->pattern == NULL || workspace->dep_heap == NULL ||
      workspace->u_row_ptr == NULL || workspace->u_row_end == NULL) {
    kls_row_first_workspace_free(workspace);
    return 0;
  }
  return 1;
}

static void kls_row_first_private_worker_free(
  kls_row_first_private_worker *worker) {
  if (worker == NULL) {
    return;
  }
  kls_row_first_workspace_free(&worker->workspace);
  kls_row_first_entries_free(&worker->l_entries);
  kls_row_first_entries_free(&worker->u_entries);
}

static void *kls_row_first_private_worker_main(void *arg) {
  kls_row_first_private_worker *worker =
    (kls_row_first_private_worker *)arg;
  if (worker == NULL || worker->ctx == NULL || worker->plan == NULL ||
      worker->row_owner == NULL || worker->udiag_values == NULL ||
      worker->tid < 0) {
    return NULL;
  }
  const UF_long nk = worker->ctx->nk;
  worker->ok = 0;
  worker->completed_rows = 0;
  worker->failed_pos = KLS_KLU_EMPTY;
  worker->failed_row = KLS_KLU_EMPTY;
  worker->failed_global_row = KLS_KLU_EMPTY;
  kls_row_first_factor_row_failure_clear(&worker->failure, KLS_KLU_EMPTY,
                                         worker->tid);
  kls_row_first_supernode_panel_cache panel_cache;
  memset(&panel_cache, 0, sizeof(panel_cache));
  UF_long open_panel_start = KLS_KLU_EMPTY;
  UF_long open_panel_end = KLS_KLU_EMPTY;
  if (!kls_row_first_workspace_init(&worker->workspace, nk) ||
      !kls_row_first_entries_enable_col_counts(&worker->l_entries, nk) ||
      !kls_row_first_entries_enable_col_counts(&worker->u_entries, nk) ||
      !kls_row_first_entries_enable_col_links(&worker->u_entries, nk)) {
    kls_row_first_factor_row_failure_set(
      &worker->failure, KLS_ROW_FIRST_FACTOR_ROW_FAIL_SETUP,
      KLS_KLU_EMPTY, KLS_KLU_EMPTY, worker->tid, -1);
    goto done;
  }
  (void)kls_row_first_workspace_enable_supernode_fit_cache(
    &worker->workspace, nk);
  if (kls_row_first_supernode_panel_cache_init_rows(&panel_cache, nk)) {
    worker->workspace.supernode_panel_cache = &panel_cache;
  }

  for (UF_long pos = worker->begin; pos < worker->end; ++pos) {
    const UF_long global_row = worker->plan->private_rows[pos];
    const UF_long l_count_before = worker->l_entries.count;
    const UF_long u_count_before = worker->u_entries.count;
    if (global_row < worker->ctx->k1 ||
        global_row >= worker->ctx->k1 + nk) {
      worker->failed_pos = pos;
      worker->failed_global_row = global_row;
      (void)kls_row_first_entries_truncate(&worker->l_entries,
                                           l_count_before);
      (void)kls_row_first_entries_truncate(&worker->u_entries,
                                           u_count_before);
      kls_row_first_factor_row_failure_set(
        &worker->failure, KLS_ROW_FIRST_FACTOR_ROW_FAIL_INPUT,
        KLS_KLU_EMPTY, global_row, worker->tid, -1);
      goto done;
    }
    const UF_long row = global_row - worker->ctx->k1;
    const UF_long pivots_before = worker->stats.dynamic_column_pivots;
    if (worker->row_owner[row] != worker->tid) {
      worker->failed_pos = pos;
      worker->failed_row = row;
      worker->failed_global_row = global_row;
      (void)kls_row_first_entries_truncate(&worker->l_entries,
                                           l_count_before);
      (void)kls_row_first_entries_truncate(&worker->u_entries,
                                           u_count_before);
      worker->workspace.u_row_ptr[row] = 0u;
      worker->workspace.u_row_end[row] = 0u;
      worker->udiag_values[row] = 0.0;
      kls_row_first_factor_row_failure_set(
        &worker->failure, KLS_ROW_FIRST_FACTOR_ROW_FAIL_DEPENDENCY,
        row, row, worker->tid, worker->row_owner[row]);
      goto done;
    }
    if (!kls_row_first_factor_one_row(
          worker->ctx, &worker->workspace, &worker->l_entries,
          &worker->u_entries, worker->udiag_values, row,
          worker->row_owner, worker->tid, NULL, &worker->stats,
          &worker->failure)) {
      worker->failed_pos = pos;
      worker->failed_row = row;
      worker->failed_global_row = global_row;
      (void)kls_row_first_entries_truncate(&worker->l_entries,
                                           l_count_before);
      (void)kls_row_first_entries_truncate(&worker->u_entries,
                                           u_count_before);
      worker->workspace.u_row_ptr[row] = 0u;
      worker->workspace.u_row_end[row] = 0u;
      worker->udiag_values[row] = 0.0;
      goto done;
    }
    if (worker->stats.dynamic_column_pivots != pivots_before) {
      kls_row_first_supernode_panel_cache_reset_active(
        &panel_cache, &worker->workspace, nk, &open_panel_start,
        &open_panel_end);
      free(worker->workspace.supernode_fit_status);
      worker->workspace.supernode_fit_status = NULL;
    }
    kls_row_first_supernode_panel_cache_publish_completed(
      worker->workspace.supernode_panel_cache != NULL ? &panel_cache : NULL,
      &worker->u_entries, worker->workspace.u_row_ptr,
      worker->workspace.u_row_end, worker->udiag_values, nk,
      &worker->stats, &open_panel_start, &open_panel_end, row);
    worker->completed_rows++;
  }
  worker->ok = 1;

done:
  worker->workspace.supernode_panel_cache = NULL;
  kls_row_first_supernode_panel_cache_free(&panel_cache);
  return NULL;
}

static int kls_first_separator_private_ownership_valid(
  UF_long k1,
  UF_long nk,
  UF_long n,
  const UF_long *row_ptr,
  const UF_long *row_cols,
  const UF_long *col_pos,
  const kls_first_separator_queue_plan *plan,
  int thread_count) {
  if (nk == 0u || k1 > UF_long_max - nk || row_ptr == NULL ||
      row_cols == NULL || col_pos == NULL || plan == NULL ||
      plan->thread_ptr == NULL || plan->private_rows == NULL ||
      thread_count <= 1 || nk > (UF_long)(SIZE_MAX / sizeof(int))) {
    return 0;
  }

  int *row_owner = (int *)malloc((size_t)nk * sizeof(*row_owner));
  UF_long *post_rank =
    (UF_long *)malloc((size_t)nk * sizeof(*post_rank));
  if (row_owner == NULL || post_rank == NULL) {
    free(row_owner);
    free(post_rank);
    return 0;
  }
  for (UF_long row = 0; row < nk; ++row) {
    row_owner[row] = -1;
  }

  int ok = 1;
  for (int tid = 0; tid < thread_count && ok; ++tid) {
    const UF_long begin = plan->thread_ptr[tid];
    const UF_long end = plan->thread_ptr[tid + 1];
    if (begin > end || end > plan->private_rows_count) {
      ok = 0;
      break;
    }
    for (UF_long pos = begin; pos < end; ++pos) {
      const UF_long global_row = plan->private_rows[pos];
      if (global_row < k1 || global_row >= k1 + nk) {
        ok = 0;
        break;
      }
      const UF_long row = global_row - k1;
      if (row_owner[row] != -1) {
        ok = 0;
        break;
      }
      row_owner[row] = tid;
    }
  }
  for (UF_long pos = 0; pos < plan->pipeline_rows_count && ok; ++pos) {
    const UF_long global_row = plan->pipeline_rows[pos];
    if (global_row < k1 || global_row >= k1 + nk) {
      ok = 0;
      break;
    }
    const UF_long row = global_row - k1;
    if (row_owner[row] != -1) {
      ok = 0;
      break;
    }
    row_owner[row] = -2;
  }
  for (UF_long row = 0; row < nk && ok; ++row) {
    if (row_owner[row] == -1) {
      ok = 0;
    }
  }

  /* Validate the separator schedule in the coordinates in which it is
     executed.  The caller stably moves every private row before every
     pipeline/separator row and applies the same permutation to columns.  A
     separator entry that is below a private row in the incoming ordering is
     therefore an upper entry after that permutation, not a private-phase
     dependency.  Testing the incoming coordinates rejected otherwise valid
     collapsed-subtree plans and fell back to a much larger synchronized
     separator suffix. */
  if (ok) {
    UF_long rank = 0u;
    for (UF_long row = 0; row < nk; ++row) {
      if (row_owner[row] >= 0) {
        post_rank[row] = rank++;
      }
    }
    for (UF_long row = 0; row < nk; ++row) {
      if (row_owner[row] == -2) {
        post_rank[row] = rank++;
      }
    }
    if (rank != nk) {
      ok = 0;
    }
  }

  for (int tid = 0; tid < thread_count && ok; ++tid) {
    for (UF_long pos = plan->thread_ptr[tid];
         pos < plan->thread_ptr[tid + 1] && ok; ++pos) {
      const UF_long row = plan->private_rows[pos] - k1;
      if (row >= nk || row_owner[row] != tid ||
          row_ptr[row] > row_ptr[row + 1u]) {
        ok = 0;
        break;
      }
      for (UF_long p = row_ptr[row]; p < row_ptr[row + 1u]; ++p) {
        const UF_long oldcol = row_cols[p];
        if (oldcol >= n || col_pos[oldcol] >= nk) {
          ok = 0;
          break;
        }
        const UF_long col = col_pos[oldcol];
        const int owner = row_owner[col];
        if (owner == -1 ||
            (post_rank[col] <= post_rank[row] && owner != tid)) {
          if (getenv("KLS_TRACE_FIRST_SEPARATOR_OWNERSHIP") != NULL) {
            fprintf(stderr,
                    "KLS first separator ownership conflict:"
                    " tid=%d row=%" PRIu64 " col=%" PRIu64
                    " post_row=%" PRIu64 " post_col=%" PRIu64
                    " owner=%d oldcol=%" PRIu64 "\n",
                    tid, (uint64_t)row, (uint64_t)col,
                    (uint64_t)post_rank[row], (uint64_t)post_rank[col], owner,
                    (uint64_t)oldcol);
          }
          ok = 0;
          break;
        }
      }
    }
  }

  free(row_owner);
  free(post_rank);
  return ok;
}

static int kls_row_first_build_private_row_owners(
  const kls_row_first_block_context *ctx,
  const kls_first_separator_queue_plan *plan,
  int thread_count,
  int *row_owner) {
  if (ctx == NULL || plan == NULL || row_owner == NULL ||
      plan->thread_ptr == NULL || plan->private_rows == NULL ||
      thread_count <= 1 || ctx->nk == 0u) {
    return 0;
  }
  for (UF_long row = 0; row < ctx->nk; ++row) {
    row_owner[row] = -1;
  }
  for (int tid = 0; tid < thread_count; ++tid) {
    const UF_long begin = plan->thread_ptr[tid];
    const UF_long end = plan->thread_ptr[tid + 1];
    if (begin > end || end > plan->private_rows_count) {
      return 0;
    }
    for (UF_long pos = begin; pos < end; ++pos) {
      const UF_long global_row = plan->private_rows[pos];
      if (global_row < ctx->k1 || global_row >= ctx->k1 + ctx->nk) {
        return 0;
      }
      const UF_long row = global_row - ctx->k1;
      if (row_owner[row] != -1) {
        return 0;
      }
      row_owner[row] = tid;
    }
  }
  for (UF_long pos = 0; pos < plan->pipeline_rows_count; ++pos) {
    const UF_long global_row = plan->pipeline_rows[pos];
    if (global_row < ctx->k1 || global_row >= ctx->k1 + ctx->nk) {
      return 0;
    }
    const UF_long row = global_row - ctx->k1;
    if (row_owner[row] != -1) {
      return 0;
    }
    row_owner[row] = -2;
  }
  for (UF_long row = 0; row < ctx->nk; ++row) {
    if (row_owner[row] == -1) {
      return 0;
    }
  }

  for (int tid = 0; tid < thread_count; ++tid) {
    for (UF_long pos = plan->thread_ptr[tid];
         pos < plan->thread_ptr[tid + 1]; ++pos) {
      const UF_long row = plan->private_rows[pos] - ctx->k1;
      if (row_owner[row] != tid) {
        return 0;
      }
      for (UF_long p = ctx->row_ptr[row]; p < ctx->row_ptr[row + 1u]; ++p) {
        const UF_long oldcol = ctx->row_cols[p];
        if (oldcol >= ctx->n || ctx->col_pos[oldcol] >= ctx->nk) {
          return 0;
        }
        const UF_long col = ctx->col_pos[oldcol];
        const int owner = row_owner[col];
        if (owner == -1 || (col <= row && owner != tid)) {
          return 0;
        }
      }
    }
  }
  return 1;
}

static int kls_row_first_append_entries(
  kls_row_first_entries *dst,
  const kls_row_first_entries *src) {
  if (dst == NULL || src == NULL) {
    return 0;
  }
  if (src->count == 0u) {
    return 1;
  }
  if (dst == src ||
      !kls_row_first_entries_reserve_append(dst, src->count)) {
    return 0;
  }
  if (dst->col_count != NULL) {
    for (UF_long p = 0; p < src->count; ++p) {
      if (src->col[p] >= dst->col_count_len) {
        return 0;
      }
    }
  }
  if (dst->col_head == NULL) {
    const UF_long pos = dst->count;
    memcpy(dst->row + pos, src->row,
           (size_t)src->count * sizeof(*dst->row));
    memcpy(dst->col + pos, src->col,
           (size_t)src->count * sizeof(*dst->col));
    memcpy(dst->value + pos, src->value,
           (size_t)src->count * sizeof(*dst->value));
    if (dst->col_count != NULL) {
      for (UF_long p = 0; p < src->count; ++p) {
        dst->col_count[src->col[p]]++;
      }
    }
    dst->count += src->count;
    return 1;
  }
  for (UF_long p = 0; p < src->count; ++p) {
    if (!kls_row_first_entries_append_reserved(
          dst, src->row[p], src->col[p], src->value[p])) {
      return 0;
    }
  }
  return 1;
}

static int kls_row_first_private_failure_allows_partial_commit(
  kls_row_first_factor_row_failure reason) {
  switch (reason) {
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_NONE:
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_SETUP:
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_INPUT:
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_VALUE:
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_DEPENDENCY:
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_UPDATE:
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_PIVOT_SCOPE:
      return 1;
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_PIVOT_REJECT:
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_PIVOT_OWNER:
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_PIVOT_EXCHANGE:
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_PIVOT_VERIFY:
    case KLS_ROW_FIRST_FACTOR_ROW_FAIL_APPEND:
      return 0;
  }
  return 0;
}

static int kls_row_first_run_parallel_private_phase(
  const kls_row_first_block_context *ctx,
  const kls_first_separator_queue_plan *plan,
  int thread_count,
  const int *row_owner,
  kls_row_first_entries *l_entries,
  kls_row_first_entries *u_entries,
  double *udiag_values,
  UF_long *u_row_ptr,
  UF_long *u_row_end,
  unsigned char *row_done,
  kls_row_first_row_stats *stats,
  UF_long *private_rows_out,
  UF_long *thread_rows_out,
  int allow_partial,
  int *all_private_completed_out) {
  KLS_SET_OPTIONAL_OUTPUT(private_rows_out, 0u);
  KLS_SET_OPTIONAL_OUTPUT(thread_rows_out, 0u);
  KLS_SET_OPTIONAL_OUTPUT(all_private_completed_out, 0);
  if (ctx == NULL || plan == NULL || thread_count <= 1 ||
      row_owner == NULL || l_entries == NULL || u_entries == NULL ||
      udiag_values == NULL || u_row_ptr == NULL || u_row_end == NULL ||
      row_done == NULL || plan->thread_ptr == NULL ||
      plan->private_rows == NULL || plan->private_rows_count == 0u ||
      plan->nonempty_threads <= 1u) {
    return 0;
  }
  pthread_t *threads =
    (pthread_t *)calloc((size_t)thread_count, sizeof(*threads));
  kls_row_first_private_worker *workers =
    (kls_row_first_private_worker *)
      calloc((size_t)thread_count, sizeof(*workers));
  if (threads == NULL || workers == NULL) {
    free(threads);
    free(workers);
    return 0;
  }

  int created = 0;
  for (int tid = 0; tid < thread_count; ++tid) {
    workers[tid].ctx = ctx;
    workers[tid].plan = plan;
    workers[tid].row_owner = row_owner;
    workers[tid].udiag_values = udiag_values;
    workers[tid].begin = plan->thread_ptr[tid];
    workers[tid].end = plan->thread_ptr[tid + 1];
    workers[tid].tid = tid;
    if (pthread_create(&threads[tid], NULL,
                       kls_row_first_private_worker_main,
                       workers + tid) != 0) {
      break;
    }
    created++;
  }
  for (int tid = 0; tid < created; ++tid) {
    (void)pthread_join(threads[tid], NULL);
  }

  int ok = created == thread_count;
  UF_long thread_rows = 0;
  UF_long commit_private_rows = 0;
  UF_long commit_thread_rows = 0;
  int partial_safe = created == thread_count;
  for (int tid = 0; tid < thread_count; ++tid) {
    const UF_long rows = plan->thread_ptr[tid + 1] - plan->thread_ptr[tid];
    if (rows > 0u) {
      thread_rows++;
    }
    if (workers[tid].completed_rows > rows) {
      ok = 0;
      partial_safe = 0;
    }
    if (rows > 0u && !workers[tid].ok) {
      ok = 0;
      if (!kls_row_first_private_failure_allows_partial_commit(
            workers[tid].failure.reason)) {
        partial_safe = 0;
      }
      if (workers[tid].failed_row == KLS_KLU_EMPTY ||
          workers[tid].failed_row >= ctx->nk) {
        partial_safe = 0;
      }
    }
  }
  for (int tid = 0; tid < thread_count; ++tid) {
    const UF_long completed_end =
      workers[tid].begin + workers[tid].completed_rows;
    UF_long thread_commit_rows = 0;
    for (UF_long pos = workers[tid].begin; pos < completed_end; ++pos) {
      if (plan->private_rows[pos] < ctx->k1 ||
          plan->private_rows[pos] >= ctx->k1 + ctx->nk) {
        partial_safe = 0;
        continue;
      }
      commit_private_rows++;
      thread_commit_rows++;
    }
    if (thread_commit_rows > 0u) {
      commit_thread_rows++;
    }
  }
  if (!ok && kls_first_factor_trace_env_enabled()) {
    for (int tid = 0; tid < thread_count; ++tid) {
      const UF_long rows =
        plan->thread_ptr[tid + 1] - plan->thread_ptr[tid];
      if (rows == 0u || workers[tid].ok) {
        continue;
      }
      const kls_row_first_factor_row_failure_info *failure =
        &workers[tid].failure;
      fprintf(stderr,
              "KLS first-factor separator-private failure:"
              " tid=%d thread_begin=%" PRIu64
              " thread_end=%" PRIu64
              " failed_pos=%" PRIu64
              " failed_row=%" PRIu64
              " failed_global_row=%" PRIu64
              " reason=%s related_col=%" PRIu64
              " owner=%d actual_owner=%d"
              " dynamic_pivots=%" PRIu64
              " separator_pivots=%" PRIu64
              " separator_extent_pivots=%" PRIu64
              " separator_fallbacks=%" PRIu64
              " separator_rejects=%" PRIu64
              " completed_rows=%" PRIu64
              " l_entries=%" PRIu64
              " u_entries=%" PRIu64 "\n",
              tid, (uint64_t)plan->thread_ptr[tid],
              (uint64_t)plan->thread_ptr[tid + 1],
              (uint64_t)workers[tid].failed_pos,
              (uint64_t)workers[tid].failed_row,
              (uint64_t)workers[tid].failed_global_row,
              kls_row_first_factor_row_failure_name(failure->reason),
              (uint64_t)failure->related_col, failure->owner,
              failure->actual_owner,
              (uint64_t)workers[tid].stats.dynamic_column_pivots,
              (uint64_t)workers[tid].stats.separator_dynamic_column_pivots,
              (uint64_t)
                workers[tid].stats.separator_extent_dynamic_column_pivots,
              (uint64_t)
                workers[tid].stats.separator_dynamic_column_fallbacks,
              (uint64_t)workers[tid].stats.separator_dynamic_column_rejects,
              (uint64_t)workers[tid].completed_rows,
              (uint64_t)workers[tid].l_entries.count,
              (uint64_t)workers[tid].u_entries.count);
      break;
    }
  }

  const int commit_private =
    ok || (allow_partial && partial_safe && commit_private_rows > 0u);
  if (!commit_private) {
    ok = 0;
    goto cleanup_workers;
  }

  const UF_long l_count_before = l_entries->count;
  const UF_long u_count_before = u_entries->count;
  int commit_ok = 1;
  for (int tid = 0; tid < thread_count && commit_ok; ++tid) {
    commit_ok = kls_row_first_append_entries(l_entries,
                                             &workers[tid].l_entries);
  }
  if (commit_ok) {
    for (int tid = 0; tid < thread_count && commit_ok; ++tid) {
      const UF_long completed_end =
        workers[tid].begin + workers[tid].completed_rows;
      for (UF_long pos = workers[tid].begin;
           pos < completed_end && commit_ok; ++pos) {
        const UF_long row = plan->private_rows[pos] - ctx->k1;
        if (row >= ctx->nk || row_owner[row] != tid) {
          commit_ok = 0;
          break;
        }
        u_row_ptr[row] = u_entries->count;
        for (UF_long up = workers[tid].workspace.u_row_ptr[row];
             up < workers[tid].workspace.u_row_end[row]; ++up) {
          if (up >= workers[tid].u_entries.count ||
              !kls_row_first_entries_append(
                u_entries, workers[tid].u_entries.row[up],
                workers[tid].u_entries.col[up],
                workers[tid].u_entries.value[up])) {
            commit_ok = 0;
            break;
          }
        }
        u_row_end[row] = u_entries->count;
        row_done[row] = 1u;
      }
    }
  }
  if (!commit_ok) {
    (void)kls_row_first_entries_truncate(l_entries, l_count_before);
    (void)kls_row_first_entries_truncate(u_entries, u_count_before);
    for (int tid = 0; tid < thread_count; ++tid) {
      const UF_long completed_end =
        workers[tid].begin + workers[tid].completed_rows;
      for (UF_long pos = workers[tid].begin; pos < completed_end; ++pos) {
        if (plan->private_rows[pos] < ctx->k1 ||
            plan->private_rows[pos] >= ctx->k1 + ctx->nk) {
          continue;
        }
        const UF_long row = plan->private_rows[pos] - ctx->k1;
        row_done[row] = 0u;
        u_row_ptr[row] = 0u;
        u_row_end[row] = 0u;
        udiag_values[row] = 0.0;
      }
    }
    ok = 0;
    goto cleanup_workers;
  }

  if (stats != NULL) {
    for (int tid = 0; tid < thread_count; ++tid) {
      stats->dynamic_column_pivots +=
        workers[tid].stats.dynamic_column_pivots;
      stats->separator_dynamic_column_pivots +=
        workers[tid].stats.separator_dynamic_column_pivots;
      stats->separator_extent_dynamic_column_pivots +=
        workers[tid].stats.separator_extent_dynamic_column_pivots;
      stats->separator_dynamic_column_fallbacks +=
        workers[tid].stats.separator_dynamic_column_fallbacks;
      stats->separator_dynamic_column_rejects +=
        workers[tid].stats.separator_dynamic_column_rejects;
      stats->supernode_update_groups +=
        workers[tid].stats.supernode_update_groups;
      stats->supernode_update_rows +=
        workers[tid].stats.supernode_update_rows;
      stats->supernode_panel_update_groups +=
        workers[tid].stats.supernode_panel_update_groups;
      stats->supernode_panel_update_rows +=
        workers[tid].stats.supernode_panel_update_rows;
      kls_row_first_accumulate_panel_cache_stats(stats,
                                                 &workers[tid].stats);
    }
  }
  if (private_rows_out != NULL) {
    *private_rows_out = commit_private_rows;
  }
  if (thread_rows_out != NULL) {
    *thread_rows_out = ok ? thread_rows : commit_thread_rows;
  }
  if (all_private_completed_out != NULL) {
    *all_private_completed_out =
      ok && commit_private_rows == plan->private_rows_count;
  }
  ok = 1;

cleanup_workers:
  for (int tid = 0; tid < thread_count; ++tid) {
    kls_row_first_private_worker_free(workers + tid);
  }
  free(threads);
  free(workers);
  return ok;
}

typedef enum kls_row_first_pipeline_failure {
  KLS_ROW_FIRST_PIPELINE_FAIL_NONE = 0,
  KLS_ROW_FIRST_PIPELINE_FAIL_SETUP,
  KLS_ROW_FIRST_PIPELINE_FAIL_DEPENDENCY,
  KLS_ROW_FIRST_PIPELINE_FAIL_PIVOT,
  KLS_ROW_FIRST_PIPELINE_FAIL_APPEND,
  KLS_ROW_FIRST_PIPELINE_FAIL_THREAD
} kls_row_first_pipeline_failure;

struct kls_row_first_pipeline_worker;
typedef struct kls_row_first_pipeline_worker kls_row_first_pipeline_worker;

typedef struct kls_row_first_pipeline_batch_target {
  kls_row_first_pipeline_worker *worker;
  kls_row_first_partial_row *state;
  double lij;
} kls_row_first_pipeline_batch_target;

static void kls_row_first_partial_row_clear(
  const kls_row_first_partial_row *state,
  kls_row_first_workspace *workspace) {
  if (state == NULL || workspace == NULL || workspace->x == NULL ||
      workspace->pattern == NULL) {
    return;
  }
  for (UF_long p = 0; p < state->pattern_count; ++p) {
    const UF_long col = workspace->pattern[p];
    workspace->x[col] = 0.0;
    if (workspace->mark != NULL &&
        workspace->mark[col] == state->generation) {
      workspace->mark[col] = 0u;
    }
  }
  if (state->row < UF_long_max) {
    workspace->x[state->row] = 0.0;
    if (workspace->mark != NULL &&
        workspace->mark[state->row] == state->generation) {
      workspace->mark[state->row] = 0u;
    }
  }
}

static int kls_row_first_partial_row_begin(
  const kls_row_first_block_context *ctx,
  kls_row_first_workspace *workspace,
  kls_row_first_partial_row *state,
  UF_long i) {
  if (ctx == NULL || workspace == NULL || state == NULL ||
      workspace->x == NULL || workspace->mark == NULL ||
      workspace->pattern == NULL || workspace->dep_heap == NULL ||
      ctx->solver == NULL || ctx->numeric == NULL ||
      ctx->numeric_values == NULL || ctx->row_ptr == NULL ||
      ctx->row_cols == NULL || ctx->row_input_pos == NULL ||
      ctx->col_pos == NULL || i >= ctx->nk) {
    return 0;
  }
  memset(state, 0, sizeof(*state));
  state->row = i;
  state->generation = (unsigned int)(i + 1u);
  if (state->generation == 0u) {
    return 0;
  }

  kls_solver *solver = ctx->solver;
  double *x = workspace->x;
  unsigned int *mark = workspace->mark;
  UF_long *pattern = workspace->pattern;
  UF_long *dep_heap = workspace->dep_heap;
  for (UF_long p = ctx->row_ptr[i]; p < ctx->row_ptr[i + 1u]; ++p) {
    const UF_long oldcol = ctx->row_cols[p];
    if (oldcol >= ctx->n || ctx->col_pos[oldcol] >= ctx->nk ||
        ctx->row_input_pos[p] >= solver->nnz) {
      return 0;
    }
    const UF_long col = ctx->col_pos[oldcol];
    if (mark[col] != state->generation) {
      mark[col] = state->generation;
      pattern[state->pattern_count++] = col;
      if (col < i) {
        kls_row_first_heap_push(dep_heap, &state->dep_heap_size, col);
      }
      if (!kls_refactor_input_value(
            solver, ctx->numeric_values,
            ctx->scaled ? ctx->numeric->Rs : NULL,
            ctx->scale, ctx->row_input_pos[p], x + col)) {
        return 0;
      }
    } else {
      double value = 0.0;
      if (!kls_refactor_input_value(
            solver, ctx->numeric_values,
            ctx->scaled ? ctx->numeric->Rs : NULL,
            ctx->scale, ctx->row_input_pos[p], &value)) {
        return 0;
      }
      x[col] += value;
    }
  }
  return 1;
}

static int kls_row_first_u_rows_fit_supernode(
  const kls_row_first_entries *u_entries,
  const UF_long *u_row_ptr,
  const UF_long *u_row_end,
  UF_long prev,
  UF_long row,
  UF_long nk) {
  if (u_entries == NULL || u_row_ptr == NULL || u_row_end == NULL ||
      row == 0u || prev + 1u != row || row >= nk ||
      u_row_ptr[prev] > u_row_end[prev] ||
      u_row_ptr[row] > u_row_end[row] ||
      u_row_end[prev] > u_entries->count ||
      u_row_end[row] > u_entries->count) {
    return 0;
  }
  UF_long p = u_row_ptr[prev];
  const UF_long pend = u_row_end[prev];
  UF_long q = u_row_ptr[row];
  const UF_long qend = u_row_end[row];
  if (p >= pend || u_entries->col[p] != row) {
    return 0;
  }
  p++;
  const UF_long prev_tail_len = pend - p;
  if (prev_tail_len != qend - q) {
    return 0;
  }
  while (p < pend) {
    const UF_long pcol = u_entries->col[p];
    const UF_long qcol = u_entries->col[q];
    if (pcol <= row || qcol <= row || pcol != qcol) {
      return 0;
    }
    p++;
    q++;
  }
  return p == pend && q == qend;
}

static void kls_row_first_supernode_panel_cache_free(
  kls_row_first_supernode_panel_cache *cache) {
  if (cache == NULL) {
    return;
  }
  free(cache->panel_id_by_row);
  free(cache->start);
  free(cache->width);
  free(cache->tail_begin);
  free(cache->tail_len);
  free(cache->dense_begin);
  free(cache->tail_value_begin);
  free(cache->active);
  free(cache->tail_cols);
  free(cache->dense_values);
  free(cache->tail_values);
  memset(cache, 0, sizeof(*cache));
}

static int kls_row_first_supernode_panel_shape(
  const kls_row_first_entries *u_entries,
  const UF_long *u_row_ptr,
  const UF_long *u_row_end,
  UF_long nk,
  UF_long start,
  UF_long end,
  UF_long *tail_len_out) {
  KLS_SET_OPTIONAL_OUTPUT(tail_len_out, 0);
  if (u_entries == NULL || u_row_ptr == NULL || u_row_end == NULL ||
      start >= nk || end <= start || end >= nk ||
      u_row_ptr[end] > u_row_end[end] ||
      u_row_end[end] > u_entries->count) {
    return 0;
  }
  const UF_long tail_begin = u_row_ptr[end];
  const UF_long tail_len = u_row_end[end] - tail_begin;
  for (UF_long dep = start; dep <= end; ++dep) {
    if (u_row_ptr[dep] > u_row_end[dep] ||
        u_row_end[dep] > u_entries->count) {
      return 0;
    }
    const UF_long row_begin = u_row_ptr[dep];
    const UF_long row_end = u_row_end[dep];
    const UF_long internal_len = end - dep;
    if (row_end - row_begin < internal_len ||
        row_end - row_begin - internal_len != tail_len) {
      return 0;
    }
    for (UF_long offset = 0; offset < internal_len; ++offset) {
      const UF_long col = u_entries->col[row_begin + offset];
      if (col != dep + offset + 1u || col >= nk) {
        return 0;
      }
    }
    for (UF_long offset = 0; offset < tail_len; ++offset) {
      const UF_long col =
        u_entries->col[row_begin + internal_len + offset];
      if (tail_begin + offset >= u_entries->count ||
          col != u_entries->col[tail_begin + offset] ||
          col <= end || col >= nk) {
        return 0;
      }
    }
  }
  if (tail_len_out != NULL) {
    *tail_len_out = tail_len;
  }
  return 1;
}

static int kls_row_first_supernode_panel_cache_add_sizes(
  UF_long width,
  UF_long tail_len,
  UF_long *dense_values_out,
  UF_long *tail_values_out,
  UF_long *tail_cols_out) {
  if (dense_values_out == NULL || tail_values_out == NULL ||
      tail_cols_out == NULL ||
      width == 0u ||
      width > UF_long_max / width ||
      (tail_len > 0u && width > UF_long_max / tail_len)) {
    return 0;
  }
  const UF_long dense_values = width * width;
  const UF_long tail_values = width * tail_len;
  if (*dense_values_out > UF_long_max - dense_values ||
      *tail_values_out > UF_long_max - tail_values ||
      *tail_cols_out > UF_long_max - tail_len) {
    return 0;
  }
  *dense_values_out += dense_values;
  *tail_values_out += tail_values;
  *tail_cols_out += tail_len;
  return 1;
}

static int kls_row_first_supernode_panel_cache_size_fits(
  UF_long count,
  size_t entry_size) {
  return count == 0u || entry_size == 0u ||
    count <= (UF_long)(SIZE_MAX / entry_size);
}

static int kls_row_first_supernode_panel_cache_init_rows(
  kls_row_first_supernode_panel_cache *cache,
  UF_long nk) {
  if (cache == NULL ||
      !kls_row_first_supernode_panel_cache_size_fits(nk, sizeof(UF_long))) {
    return 0;
  }
  kls_row_first_supernode_panel_cache_free(cache);
  cache->panel_id_by_row =
    (UF_long *)malloc((size_t)nk * sizeof(*cache->panel_id_by_row));
  if (cache->panel_id_by_row == NULL && nk > 0u) {
    kls_row_first_supernode_panel_cache_free(cache);
    return 0;
  }
  for (UF_long row = 0; row < nk; ++row) {
    cache->panel_id_by_row[row] = KLS_KLU_EMPTY;
  }
  cache->row_count = nk;
  return 1;
}

static void kls_row_first_supernode_panel_cache_clear(
  kls_row_first_supernode_panel_cache *cache) {
  if (cache == NULL) {
    return;
  }
  if (cache->panel_id_by_row != NULL) {
    for (UF_long panel = 0; panel < cache->panel_count; ++panel) {
      if (cache->width == NULL || cache->start == NULL ||
          panel >= cache->panel_capacity) {
        break;
      }
      const UF_long start = cache->start[panel];
      const UF_long width = cache->width[panel];
      if (start >= cache->row_count ||
          width > cache->row_count - start) {
        continue;
      }
      for (UF_long local = 0; local < width; ++local) {
        const UF_long row = start + local;
        if (cache->panel_id_by_row[row] == panel) {
          cache->panel_id_by_row[row] = KLS_KLU_EMPTY;
        }
      }
    }
  }
  cache->panel_count = 0;
  cache->dense_value_count = 0;
  cache->tail_value_count = 0;
  cache->tail_col_count = 0;
  cache->enabled = 0;
}

static void kls_row_first_supernode_panel_cache_reset_active(
  kls_row_first_supernode_panel_cache *cache,
  kls_row_first_workspace *workspace,
  UF_long nk,
  UF_long *open_start_io,
  UF_long *open_end_io) {
  if (open_start_io != NULL) {
    *open_start_io = KLS_KLU_EMPTY;
  }
  if (open_end_io != NULL) {
    *open_end_io = KLS_KLU_EMPTY;
  }
  if (workspace != NULL) {
    workspace->supernode_panel_cache = NULL;
  }
  if (cache == NULL) {
    return;
  }
  if (cache->panel_id_by_row == NULL || cache->row_count != nk) {
    kls_row_first_supernode_panel_cache_free(cache);
    if (workspace != NULL &&
        kls_row_first_supernode_panel_cache_init_rows(cache, nk)) {
      workspace->supernode_panel_cache = cache;
    }
    return;
  }
  kls_row_first_supernode_panel_cache_clear(cache);
  if (workspace != NULL) {
    workspace->supernode_panel_cache = cache;
  }
}

static UF_long kls_row_first_supernode_panel_cache_growth(
  UF_long capacity,
  UF_long needed) {
  UF_long grown = capacity == 0u ? 4u : capacity;
  while (grown < needed) {
    const UF_long next = 2u * grown;
    if (next <= grown) {
      return needed;
    }
    grown = next;
  }
  return grown;
}

static int kls_row_first_supernode_panel_cache_reserve_uf(
  UF_long **array_io,
  UF_long used,
  UF_long *capacity_io,
  UF_long needed) {
  if (array_io == NULL || capacity_io == NULL) {
    return 0;
  }
  if (needed <= *capacity_io) {
    return 1;
  }
  const UF_long capacity =
    kls_row_first_supernode_panel_cache_growth(*capacity_io, needed);
  if (!kls_row_first_supernode_panel_cache_size_fits(
        capacity, sizeof(UF_long))) {
    return 0;
  }
  UF_long *grown =
    (UF_long *)realloc(*array_io, (size_t)capacity * sizeof(*grown));
  if (grown == NULL) {
    return 0;
  }
  (void)used;
  *array_io = grown;
  *capacity_io = capacity;
  return 1;
}

static int kls_row_first_supernode_panel_cache_reserve_double(
  double **array_io,
  UF_long used,
  UF_long *capacity_io,
  UF_long needed) {
  if (array_io == NULL || capacity_io == NULL) {
    return 0;
  }
  if (needed <= *capacity_io) {
    return 1;
  }
  const UF_long capacity =
    kls_row_first_supernode_panel_cache_growth(*capacity_io, needed);
  if (!kls_row_first_supernode_panel_cache_size_fits(
        capacity, sizeof(double))) {
    return 0;
  }
  double *grown =
    (double *)realloc(*array_io, (size_t)capacity * sizeof(*grown));
  if (grown == NULL) {
    return 0;
  }
  (void)used;
  *array_io = grown;
  *capacity_io = capacity;
  return 1;
}

static int kls_row_first_supernode_panel_cache_reserve_panels(
  kls_row_first_supernode_panel_cache *cache,
  UF_long needed) {
  if (cache == NULL) {
    return 0;
  }
  if (needed <= cache->panel_capacity) {
    return 1;
  }
  const UF_long capacity =
    kls_row_first_supernode_panel_cache_growth(cache->panel_capacity,
                                               needed);
  if (!kls_row_first_supernode_panel_cache_size_fits(
        capacity, sizeof(UF_long))) {
    return 0;
  }
  UF_long *start =
    (UF_long *)realloc(cache->start, (size_t)capacity * sizeof(*start));
  if (start != NULL) {
    cache->start = start;
  }
  UF_long *width =
    (UF_long *)realloc(cache->width, (size_t)capacity * sizeof(*width));
  if (width != NULL) {
    cache->width = width;
  }
  UF_long *tail_begin =
    (UF_long *)realloc(cache->tail_begin,
                       (size_t)capacity * sizeof(*tail_begin));
  if (tail_begin != NULL) {
    cache->tail_begin = tail_begin;
  }
  UF_long *tail_len =
    (UF_long *)realloc(cache->tail_len,
                       (size_t)capacity * sizeof(*tail_len));
  if (tail_len != NULL) {
    cache->tail_len = tail_len;
  }
  UF_long *dense_begin =
    (UF_long *)realloc(cache->dense_begin,
                       (size_t)capacity * sizeof(*dense_begin));
  if (dense_begin != NULL) {
    cache->dense_begin = dense_begin;
  }
  UF_long *tail_value_begin =
    (UF_long *)realloc(cache->tail_value_begin,
                       (size_t)capacity * sizeof(*tail_value_begin));
  if (tail_value_begin != NULL) {
    cache->tail_value_begin = tail_value_begin;
  }
  unsigned char *active =
    (unsigned char *)realloc(cache->active,
                             (size_t)capacity * sizeof(*active));
  if (active != NULL) {
    cache->active = active;
  }
  if (start == NULL || width == NULL || tail_begin == NULL ||
      tail_len == NULL || dense_begin == NULL ||
      tail_value_begin == NULL || active == NULL) {
    return 0;
  }
  cache->panel_capacity = capacity;
  return 1;
}

static int kls_row_first_supernode_panel_cache_append(
  kls_row_first_supernode_panel_cache *cache,
  const kls_row_first_entries *u_entries,
  const UF_long *u_row_ptr,
  const UF_long *u_row_end,
  const double *udiag_values,
  UF_long nk,
  UF_long start,
  UF_long end,
  UF_long max_stored_entries,
  kls_row_first_row_stats *stats) {
  if (cache == NULL || cache->panel_id_by_row == NULL ||
      u_entries == NULL || u_row_ptr == NULL || u_row_end == NULL ||
      udiag_values == NULL || cache->row_count != nk ||
      start >= nk || end <= start || end >= nk) {
    return 0;
  }
  const UF_long panel_width = end - start + 1u;
  if (panel_width < KLS_EGRAPH_CACHED_SUPERNODE_MIN_ROWS) {
    return 0;
  }
  UF_long tail_len = 0;
  if (!kls_row_first_supernode_panel_shape(
        u_entries, u_row_ptr, u_row_end, nk, start, end, &tail_len)) {
    return 0;
  }
  if (!kls_egraph_cached_supernode_update_allows(panel_width, 0u,
                                                 tail_len)) {
    return 0;
  }
  if (!kls_row_first_panel_cache_allows_missed_append(
        cache, panel_width, tail_len, max_stored_entries)) {
    return 0;
  }
  const UF_long panel_count_needed = cache->panel_count + 1u;
  if (panel_count_needed <= cache->panel_count ||
      !kls_row_first_supernode_panel_cache_reserve_panels(
        cache, panel_count_needed)) {
    return 0;
  }
  UF_long dense_needed = cache->dense_value_count;
  UF_long tail_value_needed = cache->tail_value_count;
  UF_long tail_col_needed = cache->tail_col_count;
  if (!kls_row_first_supernode_panel_cache_add_sizes(
        panel_width, tail_len, &dense_needed, &tail_value_needed,
        &tail_col_needed)) {
    return 0;
  }
  if (!kls_row_first_supernode_panel_cache_reserve_double(
        &cache->dense_values, cache->dense_value_count,
        &cache->dense_value_capacity, dense_needed) ||
      !kls_row_first_supernode_panel_cache_reserve_double(
        &cache->tail_values, cache->tail_value_count,
        &cache->tail_value_capacity, tail_value_needed) ||
      !kls_row_first_supernode_panel_cache_reserve_uf(
        &cache->tail_cols, cache->tail_col_count,
        &cache->tail_col_capacity, tail_col_needed)) {
    return 0;
  }

  const UF_long panel = cache->panel_count++;
  const UF_long dense_begin = cache->dense_value_count;
  const UF_long tail_begin = cache->tail_col_count;
  const UF_long tail_value_begin = cache->tail_value_count;
  cache->start[panel] = start;
  cache->width[panel] = panel_width;
  cache->tail_begin[panel] = tail_begin;
  cache->tail_len[panel] = tail_len;
  cache->dense_begin[panel] = dense_begin;
  cache->tail_value_begin[panel] = tail_value_begin;
  cache->active[panel] = 1u;
  const UF_long source_tail_begin = u_row_ptr[end];
  for (UF_long offset = 0; offset < tail_len; ++offset) {
    cache->tail_cols[tail_begin + offset] =
      u_entries->col[source_tail_begin + offset];
  }
  for (UF_long local = 0; local < panel_width; ++local) {
    const UF_long dep = start + local;
    cache->panel_id_by_row[dep] = panel;
    const UF_long row_begin = u_row_ptr[dep];
    const UF_long internal_len = end - dep;
    double *dense_row =
      cache->dense_values + dense_begin + local * panel_width;
    dense_row[local] = udiag_values[dep];
    for (UF_long offset = 0; offset < internal_len; ++offset) {
      dense_row[local + offset + 1u] =
        u_entries->value[row_begin + offset];
    }
    double *tail_values =
      tail_len > 0u
        ? cache->tail_values + tail_value_begin + local * tail_len
        : NULL;
    for (UF_long offset = 0; offset < tail_len; ++offset) {
      tail_values[offset] =
        u_entries->value[row_begin + internal_len + offset];
    }
  }
  cache->dense_value_count = dense_needed;
  cache->tail_value_count = tail_value_needed;
  cache->tail_col_count = tail_col_needed;
  cache->enabled = 1;
  kls_row_first_record_panel_cache_append(stats, panel_width, tail_len);
  return 1;
}

static void kls_row_first_supernode_panel_cache_publish_completed(
  kls_row_first_supernode_panel_cache *cache,
  const kls_row_first_entries *u_entries,
  const UF_long *u_row_ptr,
  const UF_long *u_row_end,
  const double *udiag_values,
  UF_long nk,
  kls_row_first_row_stats *stats,
  UF_long *open_start_io,
  UF_long *open_end_io,
  UF_long row) {
  if (open_start_io == NULL || open_end_io == NULL) {
    return;
  }
  if (cache == NULL || cache->panel_id_by_row == NULL ||
      u_entries == NULL || u_row_ptr == NULL || u_row_end == NULL ||
      udiag_values == NULL || row >= nk) {
    *open_start_io = row;
    *open_end_io = row;
    return;
  }
  if (*open_start_io == KLS_KLU_EMPTY ||
      *open_end_io == KLS_KLU_EMPTY ||
      *open_start_io >= nk || *open_end_io >= nk || row == 0u) {
    *open_start_io = row;
    *open_end_io = row;
    return;
  }
  if (*open_end_io + 1u == row &&
      kls_row_first_u_rows_fit_supernode(
        u_entries, u_row_ptr, u_row_end, *open_end_io, row, nk)) {
    *open_end_io = row;
    return;
  }
  if (*open_end_io > *open_start_io) {
    (void)kls_row_first_supernode_panel_cache_append(
      cache, u_entries, u_row_ptr, u_row_end, udiag_values, nk,
      *open_start_io, *open_end_io, UF_long_max, stats);
  }
  *open_start_io = row;
  *open_end_io = row;
}

static int kls_row_first_supernode_panel_cache_build(
  kls_row_first_supernode_panel_cache *cache,
  const kls_row_first_entries *u_entries,
  const UF_long *u_row_ptr,
  const UF_long *u_row_end,
  const double *udiag_values,
  const unsigned char *row_done,
  const UF_long *supernode_start,
  const UF_long *supernode_end,
  UF_long nk,
  UF_long ready_count,
  kls_row_first_row_stats *stats) {
  if (cache == NULL) {
    return 0;
  }
  kls_row_first_supernode_panel_cache_free(cache);
  if (u_entries == NULL || u_row_ptr == NULL || u_row_end == NULL ||
      udiag_values == NULL ||
      supernode_start == NULL || supernode_end == NULL ||
      ready_count > nk ||
      !kls_row_first_supernode_panel_cache_size_fits(nk, sizeof(UF_long))) {
    return 0;
  }

  UF_long panel_count = 0;
  UF_long dense_value_count = 0;
  UF_long tail_value_count = 0;
  UF_long tail_col_count = 0;
  for (UF_long row = 0; row < ready_count; ++row) {
    if (supernode_start[row] != row) {
      continue;
    }
    const UF_long end = supernode_end[row];
    if (end <= row || end >= ready_count) {
      continue;
    }
    int all_done = 1;
    if (row_done != NULL) {
      for (UF_long dep = row; dep <= end; ++dep) {
        if (row_done[dep] == 0u) {
          all_done = 0;
          break;
        }
      }
    }
    if (!all_done) {
      continue;
    }
    UF_long tail_len = 0;
    if (!kls_row_first_supernode_panel_shape(
          u_entries, u_row_ptr, u_row_end, nk, row, end, &tail_len)) {
      continue;
    }
    const UF_long width = end - row + 1u;
    if (!kls_row_first_supernode_panel_cache_add_sizes(
          width, tail_len, &dense_value_count, &tail_value_count,
          &tail_col_count)) {
      return 0;
    }
    panel_count++;
    row = end;
  }
  cache->panel_id_by_row =
    (UF_long *)malloc((size_t)nk * sizeof(*cache->panel_id_by_row));
  if (cache->panel_id_by_row == NULL && nk > 0u) {
    kls_row_first_supernode_panel_cache_free(cache);
    return 0;
  }
  for (UF_long row = 0; row < nk; ++row) {
    cache->panel_id_by_row[row] = KLS_KLU_EMPTY;
  }
  cache->row_count = nk;
  if (panel_count == 0u) {
    kls_row_first_record_panel_cache_build(stats, 0u, 0u, 0u, 0u);
    return 1;
  }
  if (!kls_row_first_supernode_panel_cache_size_fits(
        panel_count, sizeof(UF_long)) ||
      !kls_row_first_supernode_panel_cache_size_fits(
        dense_value_count, sizeof(double)) ||
      !kls_row_first_supernode_panel_cache_size_fits(
        tail_value_count, sizeof(double)) ||
      !kls_row_first_supernode_panel_cache_size_fits(
        tail_col_count, sizeof(UF_long))) {
    kls_row_first_supernode_panel_cache_free(cache);
    return 0;
  }

  cache->start =
    (UF_long *)malloc((size_t)panel_count * sizeof(*cache->start));
  cache->width =
    (UF_long *)malloc((size_t)panel_count * sizeof(*cache->width));
  cache->tail_begin =
    (UF_long *)malloc((size_t)panel_count * sizeof(*cache->tail_begin));
  cache->tail_len =
    (UF_long *)malloc((size_t)panel_count * sizeof(*cache->tail_len));
  cache->dense_begin =
    (UF_long *)malloc((size_t)panel_count * sizeof(*cache->dense_begin));
  cache->tail_value_begin =
    (UF_long *)malloc((size_t)panel_count *
                      sizeof(*cache->tail_value_begin));
  cache->active =
    (unsigned char *)malloc((size_t)panel_count * sizeof(*cache->active));
  cache->tail_cols = tail_col_count > 0u
    ? (UF_long *)malloc((size_t)tail_col_count * sizeof(*cache->tail_cols))
    : NULL;
  cache->dense_values =
    (double *)calloc((size_t)dense_value_count,
                     sizeof(*cache->dense_values));
  cache->tail_values = tail_value_count > 0u
    ? (double *)malloc((size_t)tail_value_count *
                       sizeof(*cache->tail_values))
    : NULL;
  if (cache->start == NULL || cache->width == NULL ||
      cache->tail_begin == NULL || cache->tail_len == NULL ||
      cache->dense_begin == NULL || cache->tail_value_begin == NULL ||
      cache->active == NULL ||
      (tail_col_count > 0u && cache->tail_cols == NULL) ||
      cache->dense_values == NULL ||
      (tail_value_count > 0u && cache->tail_values == NULL)) {
    kls_row_first_supernode_panel_cache_free(cache);
    return 0;
  }

  UF_long next_panel = 0;
  UF_long next_dense = 0;
  UF_long next_tail_col = 0;
  UF_long next_tail_value = 0;
  for (UF_long row = 0; row < ready_count; ++row) {
    if (supernode_start[row] != row) {
      continue;
    }
    const UF_long end = supernode_end[row];
    if (end <= row || end >= ready_count) {
      continue;
    }
    int all_done = 1;
    if (row_done != NULL) {
      for (UF_long dep = row; dep <= end; ++dep) {
        if (row_done[dep] == 0u) {
          all_done = 0;
          break;
        }
      }
    }
    if (!all_done) {
      continue;
    }
    UF_long tail_len = 0;
    if (!kls_row_first_supernode_panel_shape(
          u_entries, u_row_ptr, u_row_end, nk, row, end, &tail_len)) {
      continue;
    }
    const UF_long width = end - row + 1u;
    const UF_long panel = next_panel++;
    cache->start[panel] = row;
    cache->width[panel] = width;
    cache->tail_begin[panel] = next_tail_col;
    cache->tail_len[panel] = tail_len;
    cache->dense_begin[panel] = next_dense;
    cache->tail_value_begin[panel] = next_tail_value;
    cache->active[panel] = 1u;
    const UF_long source_tail_begin = u_row_ptr[end];
    for (UF_long offset = 0; offset < tail_len; ++offset) {
      cache->tail_cols[next_tail_col + offset] =
        u_entries->col[source_tail_begin + offset];
    }
    for (UF_long local = 0; local < width; ++local) {
      const UF_long dep = row + local;
      cache->panel_id_by_row[dep] = panel;
      const UF_long row_begin = u_row_ptr[dep];
      const UF_long internal_len = end - dep;
      double *dense_row =
        cache->dense_values + next_dense + local * width;
      dense_row[local] = udiag_values[dep];
      for (UF_long offset = 0; offset < internal_len; ++offset) {
        dense_row[local + offset + 1u] =
          u_entries->value[row_begin + offset];
      }
      double *tail_values =
        tail_len > 0u
          ? cache->tail_values + next_tail_value + local * tail_len
          : NULL;
      for (UF_long offset = 0; offset < tail_len; ++offset) {
        tail_values[offset] =
          u_entries->value[row_begin + internal_len + offset];
      }
    }
    next_dense += width * width;
    next_tail_col += tail_len;
    next_tail_value += width * tail_len;
    row = end;
  }
  cache->panel_count = next_panel;
  cache->panel_capacity = panel_count;
  cache->dense_value_count = next_dense;
  cache->dense_value_capacity = dense_value_count;
  cache->tail_col_count = next_tail_col;
  cache->tail_col_capacity = tail_col_count;
  cache->tail_value_count = next_tail_value;
  cache->tail_value_capacity = tail_value_count;
  cache->enabled = next_panel > 0u ? 1 : 0;
  kls_row_first_record_panel_cache_build(
    stats, next_panel, next_dense, next_tail_value, next_tail_col);
  return 1;
}

static void kls_row_first_supernode_panel_cache_deactivate_panel(
  kls_row_first_supernode_panel_cache *cache,
  UF_long panel) {
  if (cache == NULL || cache->active == NULL ||
      cache->start == NULL || cache->width == NULL ||
      cache->panel_id_by_row == NULL ||
      panel >= cache->panel_count ||
      cache->active[panel] == 0u) {
    return;
  }
  cache->active[panel] = 0u;
  const UF_long start = cache->start[panel];
  const UF_long width = cache->width[panel];
  if (start >= cache->row_count) {
    return;
  }
  UF_long local_limit = cache->row_count - start;
  if (local_limit > width) {
    local_limit = width;
  }
  for (UF_long local = 0; local < local_limit; ++local) {
    const UF_long row = start + local;
    if (cache->panel_id_by_row[row] == panel) {
      cache->panel_id_by_row[row] = KLS_KLU_EMPTY;
    }
  }
}

static UF_long kls_row_first_supernode_panel_cache_exchange_columns(
  kls_row_first_supernode_panel_cache *cache,
  UF_long col_a,
  UF_long col_b) {
  if (col_a == col_b) {
    return 0;
  }
  if (cache == NULL || cache->active == NULL ||
      cache->panel_id_by_row == NULL ||
      cache->start == NULL || cache->width == NULL ||
      cache->tail_begin == NULL || cache->tail_len == NULL) {
    return 0;
  }
  UF_long invalidated = 0;
  UF_long still_active = 0;
  for (UF_long panel = 0; panel < cache->panel_count; ++panel) {
    if (cache->active[panel] == 0u) {
      continue;
    }
    const UF_long start = cache->start[panel];
    const UF_long width = cache->width[panel];
    const UF_long panel_end =
      width > 0u && width - 1u <= UF_long_max - start
        ? start + width - 1u : UF_long_max;
    const int a_in_dense =
      width > 0u && col_a >= start && col_a - start < width;
    const int b_in_dense =
      width > 0u && col_b >= start && col_b - start < width;
    if (a_in_dense || b_in_dense) {
      kls_row_first_supernode_panel_cache_deactivate_panel(cache, panel);
      invalidated++;
      continue;
    }
    /* Tail-only swaps preserve the cached triangular panel. */
    const UF_long tail_begin = cache->tail_begin[panel];
    const UF_long tail_len = cache->tail_len[panel];
    if (tail_len > 0u && cache->tail_cols == NULL) {
      kls_row_first_supernode_panel_cache_deactivate_panel(cache, panel);
      invalidated++;
      continue;
    }
    if (tail_begin > cache->tail_col_count ||
        tail_len > cache->tail_col_count - tail_begin) {
      kls_row_first_supernode_panel_cache_deactivate_panel(cache, panel);
      invalidated++;
      continue;
    }
    for (UF_long offset = 0; offset < tail_len; ++offset) {
      UF_long *tail_col = cache->tail_cols + tail_begin + offset;
      if (*tail_col == col_a) {
        if (col_b <= panel_end) {
          kls_row_first_supernode_panel_cache_deactivate_panel(cache,
                                                              panel);
          invalidated++;
          goto next_panel;
        }
        *tail_col = col_b;
      } else if (*tail_col == col_b) {
        if (col_a <= panel_end) {
          kls_row_first_supernode_panel_cache_deactivate_panel(cache,
                                                              panel);
          invalidated++;
          goto next_panel;
        }
        *tail_col = col_a;
      }
    }
    still_active++;
next_panel:
    ;
  }
  cache->enabled = still_active > 0u ? 1 : 0;
  return invalidated;
}

static UF_long kls_row_first_owned_supernode_end(
  const kls_row_first_entries *u_entries,
  const UF_long *u_row_ptr,
  const UF_long *u_row_end,
  const unsigned char *row_done,
  const int *row_owner,
  int owner,
  UF_long nk,
  UF_long current_row,
  unsigned char *fit_status,
  UF_long dep) {
  if (u_entries == NULL || u_row_ptr == NULL || u_row_end == NULL ||
      dep >= nk || dep >= current_row) {
    return dep;
  }

  UF_long end = dep;
  while (end + 1u < current_row) {
    const UF_long next = end + 1u;
    if (row_owner != NULL && row_owner[next] != owner) {
      break;
    }
    if (row_done != NULL && row_done[next] == 0u) {
      break;
    }
    int rows_fit = 0;
    if (fit_status != NULL && next < nk && fit_status[next] != 0u) {
      rows_fit = fit_status[next] == 2u;
    } else {
      rows_fit =
        kls_row_first_u_rows_fit_supernode(
          u_entries, u_row_ptr, u_row_end, end, next, nk);
      if (fit_status != NULL && next < nk) {
        fit_status[next] = rows_fit ? 2u : 1u;
      }
    }
    if (!rows_fit) {
      break;
    }
    end = next;
  }
  return end;
}

static void kls_row_first_supernodes_reset(
  const kls_row_first_entries *u_entries,
  const UF_long *u_row_ptr,
  const UF_long *u_row_end,
  const unsigned char *row_done,
  UF_long nk,
  UF_long ready_count,
  UF_long *supernode_start,
  UF_long *supernode_end) {
  if (supernode_start == NULL || supernode_end == NULL) {
    return;
  }
  for (UF_long row = 0; row < nk; ++row) {
    supernode_start[row] = row;
    supernode_end[row] = row;
  }
  if (u_entries == NULL || u_row_ptr == NULL || u_row_end == NULL ||
      ready_count > nk) {
    return;
  }
  for (UF_long row = 1; row < ready_count; ++row) {
    if (row_done != NULL &&
        (row_done[row] == 0u || row_done[row - 1u] == 0u)) {
      continue;
    }
    if (!kls_row_first_u_rows_fit_supernode(
          u_entries, u_row_ptr, u_row_end, row - 1u, row, nk)) {
      continue;
    }
    const UF_long start = supernode_start[row - 1u];
    if (start >= nk) {
      continue;
    }
    supernode_start[row] = start;
    supernode_end[start] = row;
  }
}

static void kls_row_first_supernodes_reset_prefix(
  const kls_row_first_entries *u_entries,
  const UF_long *u_row_ptr,
  const UF_long *u_row_end,
  const unsigned char *row_done,
  UF_long nk,
  UF_long ready_count,
  UF_long *supernode_start,
  UF_long *supernode_end) {
  if (supernode_start == NULL || supernode_end == NULL ||
      ready_count > nk) {
    return;
  }
  for (UF_long row = 0; row < ready_count; ++row) {
    supernode_start[row] = row;
    supernode_end[row] = row;
  }
  if (u_entries == NULL || u_row_ptr == NULL || u_row_end == NULL) {
    return;
  }
  for (UF_long row = 1; row < ready_count; ++row) {
    if (row_done != NULL &&
        (row_done[row] == 0u || row_done[row - 1u] == 0u)) {
      continue;
    }
    if (!kls_row_first_u_rows_fit_supernode(
          u_entries, u_row_ptr, u_row_end, row - 1u, row, nk)) {
      continue;
    }
    const UF_long start = supernode_start[row - 1u];
    if (start >= ready_count) {
      continue;
    }
    supernode_start[row] = start;
    supernode_end[start] = row;
  }
}

static void kls_row_first_supernodes_publish_row(
  const kls_row_first_entries *u_entries,
  const UF_long *u_row_ptr,
  const UF_long *u_row_end,
  const unsigned char *row_done,
  UF_long nk,
  UF_long merge_floor,
  UF_long row,
  UF_long *supernode_start,
  UF_long *supernode_end) {
  if (supernode_start == NULL || supernode_end == NULL || row >= nk) {
    return;
  }
  supernode_start[row] = row;
  supernode_end[row] = row;
  if (row == 0u || row <= merge_floor || row_done == NULL ||
      row_done[row - 1u] == 0u) {
    return;
  }
  if (!kls_row_first_u_rows_fit_supernode(
        u_entries, u_row_ptr, u_row_end, row - 1u, row, nk)) {
    return;
  }
  const UF_long start = supernode_start[row - 1u];
  if (start >= nk) {
    return;
  }
  supernode_start[row] = start;
  supernode_end[start] = row;
}

static int kls_row_first_dependency_ready(
  const unsigned char *row_done,
  const UF_long *active_rank,
  UF_long ready_limit,
  UF_long dep);

static UF_long kls_row_first_ready_supernode_end(
  const UF_long *supernode_start,
  const UF_long *supernode_end,
  const unsigned char *row_done,
  const UF_long *active_rank,
  UF_long nk,
  UF_long supernode_valid_begin,
  UF_long ready_limit,
  UF_long current_row,
  UF_long dep) {
  if (supernode_start == NULL || supernode_end == NULL ||
      row_done == NULL || dep >= nk || dep >= current_row) {
    return dep;
  }
  if (dep < supernode_valid_begin) {
    return dep;
  }
  const UF_long start = supernode_start[dep];
  if (start >= nk || start < supernode_valid_begin || dep < start) {
    return dep;
  }
  const UF_long end = supernode_end[start];
  if (end <= dep || end >= nk) {
    return dep;
  }
  UF_long ready_end = end;
  if (active_rank == NULL && ready_end >= ready_limit) {
    if (ready_limit == 0u) {
      return dep;
    }
    ready_end = ready_limit - 1u;
  }
  if (ready_end >= current_row) {
    if (current_row == 0u) {
      return dep;
    }
    ready_end = current_row - 1u;
  }
  if (ready_end <= dep) {
    return dep;
  }
  for (UF_long row = dep + 1u; row <= ready_end; ++row) {
    if (!kls_row_first_dependency_ready(row_done, active_rank,
                                        ready_limit, row)) {
      return row - 1u;
    }
  }
  return ready_end;
}

static UF_long kls_row_first_ready_cached_panel_end(
  const kls_row_first_supernode_panel_cache *cache,
  const unsigned char *row_done,
  const UF_long *active_rank,
  UF_long nk,
  UF_long supernode_valid_begin,
  UF_long ready_limit,
  UF_long current_row,
  UF_long dep) {
  if (cache == NULL || !cache->enabled ||
      cache->panel_id_by_row == NULL ||
      cache->start == NULL || cache->width == NULL ||
      cache->active == NULL ||
      row_done == NULL || cache->row_count != nk ||
      dep >= nk || dep >= current_row) {
    return dep;
  }
  if (dep < supernode_valid_begin) {
    return dep;
  }
  const UF_long panel = cache->panel_id_by_row[dep];
  if (panel == KLS_KLU_EMPTY || panel >= cache->panel_count ||
      cache->active[panel] == 0u) {
    return dep;
  }
  const UF_long start = cache->start[panel];
  const UF_long width = cache->width[panel];
  if (start < supernode_valid_begin) {
    return dep;
  }
  if (width <= 1u || dep < start ||
      width - 1u > UF_long_max - start) {
    return dep;
  }
  UF_long ready_end = start + width - 1u;
  if (dep >= ready_end || ready_end >= nk ||
      ready_end >= cache->row_count) {
    return dep;
  }
  if (active_rank == NULL && ready_end >= ready_limit) {
    if (ready_limit == 0u) {
      return dep;
    }
    ready_end = ready_limit - 1u;
  }
  if (ready_end >= current_row) {
    if (current_row == 0u) {
      return dep;
    }
    ready_end = current_row - 1u;
  }
  if (ready_end <= dep) {
    return dep;
  }
  for (UF_long row = dep + 1u; row <= ready_end; ++row) {
    if (!kls_row_first_dependency_ready(row_done, active_rank,
                                        ready_limit, row)) {
      return row - 1u;
    }
  }
  return ready_end;
}

static int kls_row_first_dependency_ready(
  const unsigned char *row_done,
  const UF_long *active_rank,
  UF_long ready_limit,
  UF_long dep) {
  if (row_done == NULL || row_done[dep] == 0u) {
    return 0;
  }
  if (active_rank == NULL) {
    return dep < ready_limit;
  }
  const UF_long rank = active_rank[dep];
  return rank == KLS_KLU_EMPTY || rank < ready_limit;
}

static int kls_row_first_partial_apply_one_dep(
  const kls_row_first_block_context *ctx,
  kls_row_first_workspace *workspace,
  kls_row_first_entries *local_l_entries,
  const kls_row_first_entries *published_u_entries,
  const double *udiag_values,
  const UF_long *u_row_ptr,
  const UF_long *u_row_end,
  kls_row_first_partial_row *state,
  UF_long dep) {
  if (ctx == NULL || workspace == NULL || local_l_entries == NULL ||
      published_u_entries == NULL || udiag_values == NULL ||
      u_row_ptr == NULL || u_row_end == NULL || state == NULL ||
      workspace->x == NULL || workspace->mark == NULL ||
      workspace->pattern == NULL || workspace->dep_heap == NULL ||
      dep >= state->row || state->row >= ctx->nk ||
      workspace->mark[dep] != state->generation) {
    return 0;
  }

  double *x = workspace->x;
  unsigned int *mark = workspace->mark;
  UF_long *pattern = workspace->pattern;
  UF_long *dep_heap = workspace->dep_heap;
  const double dep_pivot = udiag_values[dep];
  if (dep_pivot == 0.0) {
    return 0;
  }
  const double lij = x[dep] / dep_pivot;
  if (!kls_row_first_entries_append(local_l_entries, state->row, dep,
                                    lij)) {
    return 0;
  }
  x[dep] = 0.0;
  if (u_row_ptr[dep] > u_row_end[dep] ||
      u_row_end[dep] > published_u_entries->count) {
    return 0;
  }
  for (UF_long up = u_row_ptr[dep]; up < u_row_end[dep]; ++up) {
    const UF_long col = published_u_entries->col[up];
    if (col <= dep || col >= ctx->nk) {
      return 0;
    }
    if (mark[col] != state->generation) {
      mark[col] = state->generation;
      pattern[state->pattern_count++] = col;
      if (col < state->row) {
        kls_row_first_heap_push(dep_heap, &state->dep_heap_size, col);
      }
      x[col] = 0.0;
    }
    x[col] -= lij * published_u_entries->value[up];
  }
  return 1;
}

static int kls_row_first_partial_apply_supernode_run_scalar(
  const kls_row_first_block_context *ctx,
  kls_row_first_workspace *workspace,
  kls_row_first_entries *local_l_entries,
  const kls_row_first_entries *published_u_entries,
  const double *udiag_values,
  const UF_long *u_row_ptr,
  const UF_long *u_row_end,
  kls_row_first_partial_row *state,
  UF_long dep_begin,
  UF_long dep_end,
  UF_long *run_rows_out) {
  KLS_SET_OPTIONAL_OUTPUT(run_rows_out, 0);
  if (ctx == NULL || workspace == NULL || local_l_entries == NULL ||
      published_u_entries == NULL || udiag_values == NULL ||
      u_row_ptr == NULL || u_row_end == NULL || state == NULL ||
      workspace->x == NULL || workspace->mark == NULL ||
      workspace->pattern == NULL || workspace->dep_heap == NULL ||
      dep_begin > dep_end || dep_end >= state->row ||
      state->row >= ctx->nk ||
      workspace->mark[dep_begin] != state->generation) {
    return 0;
  }

  double *x = workspace->x;
  unsigned int *mark = workspace->mark;
  UF_long *pattern = workspace->pattern;
  UF_long *dep_heap = workspace->dep_heap;
  const unsigned int generation = state->generation;
  UF_long run_rows = 0;
  if (!kls_row_first_entries_reserve_append(
        local_l_entries, dep_end - dep_begin + 1u)) {
    return 0;
  }

  /* Later rows in the run may be discovered by the triangular prefix update. */
  for (UF_long dep = dep_begin; dep <= dep_end; ++dep) {
    if (dep != dep_begin) {
      if (mark[dep] != generation) {
        break;
      }
      if (state->dep_heap_size == 0u || dep_heap[0] != dep) {
        return 0;
      }
      const UF_long popped =
        kls_row_first_heap_pop(dep_heap, &state->dep_heap_size);
      if (popped != dep) {
        return 0;
      }
    }
    const double dep_pivot = udiag_values[dep];
    if (dep_pivot == 0.0 ||
        u_row_ptr[dep] > u_row_end[dep] ||
        u_row_end[dep] > published_u_entries->count) {
      return 0;
    }
    const double lij = x[dep] / dep_pivot;
    if (!kls_row_first_entries_append_reserved(
          local_l_entries, state->row, dep, lij)) {
      return 0;
    }
    x[dep] = 0.0;
    for (UF_long up = u_row_ptr[dep]; up < u_row_end[dep]; ++up) {
      const UF_long col = published_u_entries->col[up];
      if (col <= dep || col >= ctx->nk) {
        return 0;
      }
      if (mark[col] != generation) {
        mark[col] = generation;
        pattern[state->pattern_count++] = col;
        if (col < state->row) {
          kls_row_first_heap_push(dep_heap, &state->dep_heap_size, col);
        }
        x[col] = 0.0;
      }
      x[col] -= lij * published_u_entries->value[up];
    }
    run_rows++;
  }
  if (run_rows_out != NULL) {
    *run_rows_out = run_rows;
  }
  return 1;
}

#ifdef KLS_HAVE_CBLAS
static int kls_row_first_partial_apply_supernode_run_cached_cblas(
  const kls_row_first_block_context *ctx,
  kls_row_first_workspace *workspace,
  kls_row_first_entries *local_l_entries,
  kls_row_first_partial_row *state,
  UF_long dep_begin,
  UF_long start,
  UF_long width,
  UF_long panel_offset,
  UF_long available_end,
  UF_long run_limit,
  UF_long dense_suffix_len,
  UF_long tail_len,
  const double *dense_panel,
  const double *tail_values,
  const UF_long *tail_cols,
  UF_long *run_rows_out) {
  KLS_SET_OPTIONAL_OUTPUT(run_rows_out, 0);
  if (ctx == NULL || workspace == NULL || local_l_entries == NULL ||
      state == NULL || dense_panel == NULL ||
      workspace->x == NULL || workspace->mark == NULL ||
      workspace->pattern == NULL || workspace->dep_heap == NULL ||
      run_limit <= 1u || panel_offset >= width ||
      available_end > width ||
      run_limit != available_end - panel_offset ||
      run_limit > (UF_long)INT_MAX ||
      width > (UF_long)INT_MAX ||
      dense_suffix_len > (UF_long)INT_MAX ||
      tail_len > (UF_long)INT_MAX ||
      (tail_len > 0u && (tail_values == NULL || tail_cols == NULL)) ||
      !kls_cblas_supernode_env_enabled()) {
    return 0;
  }

  if (!kls_cblas_row_first_supernode_vector_update_allows(
        run_limit, dense_suffix_len, tail_len)) {
    return 0;
  }
  if (tail_len > UF_long_max - run_limit) {
    return 0;
  }
  const UF_long workspace_len = run_limit + tail_len;
  if (!kls_row_first_workspace_reserve_supernode(workspace,
                                                 workspace_len)) {
    return 0;
  }

  double *x = workspace->x;
  unsigned int *mark = workspace->mark;
  UF_long *pattern = workspace->pattern;
  UF_long *dep_heap = workspace->dep_heap;
  double *rhs = workspace->supernode_workspace;
  double *tail_workspace = tail_len > 0u ? rhs + run_limit : NULL;
  const unsigned int generation = state->generation;

  for (UF_long local = 0; local < run_limit; ++local) {
    const UF_long dep = dep_begin + local;
    const UF_long panel_local = panel_offset + local;
    const double pivot = dense_panel[panel_local * width + panel_local];
    if (dep >= state->row || pivot == 0.0) {
      return -1;
    }
    rhs[local] = mark[dep] == generation ? x[dep] : 0.0;
    if (mark[dep] != generation) {
      if (local == 0u) {
        return -1;
      }
      mark[dep] = generation;
      pattern[state->pattern_count++] = dep;
      x[dep] = 0.0;
    }
  }

  cblas_dtrsv(CblasRowMajor, CblasUpper, CblasTrans, CblasNonUnit,
              (int)run_limit,
              dense_panel + panel_offset * width + panel_offset,
              (int)width, rhs, 1);

  for (UF_long local = 1u; local < run_limit; ++local) {
    const UF_long dep = dep_begin + local;
    if (!kls_row_first_heap_consume_ready_dep(
          dep_heap, &state->dep_heap_size, dep)) {
      return -1;
    }
  }
  if (!kls_row_first_entries_append_consecutive(
        local_l_entries, state->row, dep_begin, rhs, run_limit)) {
    return -1;
  }
  for (UF_long local = 0; local < run_limit; ++local) {
    const UF_long dep = dep_begin + local;
    x[dep] = 0.0;
  }

  if (dense_suffix_len > 0u) {
    for (UF_long offset = 0; offset < dense_suffix_len; ++offset) {
      const UF_long col = start + available_end + offset;
      if (col >= ctx->nk) {
        return -1;
      }
      if (mark[col] != generation) {
        mark[col] = generation;
        pattern[state->pattern_count++] = col;
        if (col < state->row) {
          kls_row_first_heap_push(dep_heap, &state->dep_heap_size, col);
        }
        x[col] = 0.0;
      }
    }
    cblas_dgemv(CblasRowMajor, CblasTrans, (int)run_limit,
                (int)dense_suffix_len, -1.0,
                dense_panel + panel_offset * width + available_end,
                (int)width, rhs, 1, 1.0,
                x + start + available_end, 1);
  }

  if (tail_len > 0u) {
    cblas_dgemv(CblasRowMajor, CblasTrans, (int)run_limit,
                (int)tail_len, 1.0,
                tail_values + panel_offset * tail_len, (int)tail_len,
                rhs, 1, 0.0, tail_workspace, 1);
    for (UF_long offset = 0; offset < tail_len; ++offset) {
      const UF_long col = tail_cols[offset];
      if (col >= ctx->nk) {
        return -1;
      }
      if (mark[col] != generation) {
        mark[col] = generation;
        pattern[state->pattern_count++] = col;
        if (col < state->row) {
          kls_row_first_heap_push(dep_heap, &state->dep_heap_size, col);
        }
        x[col] = 0.0;
      }
      x[col] -= tail_workspace[offset];
    }
  }

  if (run_rows_out != NULL) {
    *run_rows_out = run_limit;
  }
  return 1;
}
#endif

static int kls_row_first_partial_apply_supernode_run_cached_portable(
  const kls_row_first_block_context *ctx,
  kls_row_first_workspace *workspace,
  kls_row_first_entries *local_l_entries,
  kls_row_first_partial_row *state,
  UF_long dep_begin,
  UF_long start,
  UF_long width,
  UF_long panel_offset,
  UF_long available_end,
  UF_long run_limit,
  UF_long dense_suffix_len,
  UF_long tail_len,
  const double *dense_panel,
  const double *tail_values,
  const UF_long *tail_cols,
  UF_long *run_rows_out) {
  KLS_SET_OPTIONAL_OUTPUT(run_rows_out, 0);
  if (ctx == NULL || workspace == NULL || local_l_entries == NULL ||
      state == NULL || dense_panel == NULL ||
      workspace->x == NULL || workspace->mark == NULL ||
      workspace->pattern == NULL || workspace->dep_heap == NULL ||
      run_limit <= 1u || panel_offset >= width ||
      available_end > width ||
      run_limit != available_end - panel_offset ||
      (tail_len > 0u && (tail_values == NULL || tail_cols == NULL))) {
    return 0;
  }
  if (!kls_row_first_workspace_reserve_supernode(workspace, run_limit)) {
    return 0;
  }

  double *x = workspace->x;
  unsigned int *mark = workspace->mark;
  UF_long *pattern = workspace->pattern;
  UF_long *dep_heap = workspace->dep_heap;
  double *rhs = workspace->supernode_workspace;
  const unsigned int generation = state->generation;

  for (UF_long local = 0; local < run_limit; ++local) {
    const UF_long dep = dep_begin + local;
    const UF_long panel_local = panel_offset + local;
    if (dep >= state->row) {
      return -1;
    }
    const double pivot = dense_panel[panel_local * width + panel_local];
    if (pivot == 0.0) {
      return -1;
    }
    rhs[local] = mark[dep] == generation ? x[dep] : 0.0;
    if (mark[dep] != generation) {
      if (local == 0u) {
        return -1;
      }
      mark[dep] = generation;
      pattern[state->pattern_count++] = dep;
      x[dep] = 0.0;
    }
  }

  for (UF_long local = 0; local < run_limit; ++local) {
    const UF_long panel_local = panel_offset + local;
    double lij = rhs[local];
    for (UF_long prev = 0; prev < local; ++prev) {
      const UF_long prev_panel_local = panel_offset + prev;
      lij -= rhs[prev] *
             dense_panel[prev_panel_local * width + panel_local];
    }
    lij /= dense_panel[panel_local * width + panel_local];
    rhs[local] = lij;
  }

  for (UF_long local = 1u; local < run_limit; ++local) {
    const UF_long dep = dep_begin + local;
    if (!kls_row_first_heap_consume_ready_dep(
          dep_heap, &state->dep_heap_size, dep)) {
      return -1;
    }
  }
  if (!kls_row_first_entries_append_consecutive(
        local_l_entries, state->row, dep_begin, rhs, run_limit)) {
    return -1;
  }
  for (UF_long local = 0; local < run_limit; ++local) {
    const UF_long dep = dep_begin + local;
    x[dep] = 0.0;
  }

  for (UF_long offset = 0; offset < dense_suffix_len; ++offset) {
    const UF_long panel_col = available_end + offset;
    const UF_long col = start + panel_col;
    if (col >= ctx->nk) {
      return -1;
    }
    if (mark[col] != generation) {
      mark[col] = generation;
      pattern[state->pattern_count++] = col;
      if (col < state->row) {
        kls_row_first_heap_push(dep_heap, &state->dep_heap_size, col);
      }
      x[col] = 0.0;
    }
    double update = 0.0;
    for (UF_long local = 0; local < run_limit; ++local) {
      const UF_long panel_local = panel_offset + local;
      update += rhs[local] *
                dense_panel[panel_local * width + panel_col];
    }
    x[col] -= update;
  }

  for (UF_long offset = 0; offset < tail_len; ++offset) {
    const UF_long col = tail_cols[offset];
    if (col >= ctx->nk) {
      return -1;
    }
    if (mark[col] != generation) {
      mark[col] = generation;
      pattern[state->pattern_count++] = col;
      if (col < state->row) {
        kls_row_first_heap_push(dep_heap, &state->dep_heap_size, col);
      }
      x[col] = 0.0;
    }
    double update = 0.0;
    for (UF_long local = 0; local < run_limit; ++local) {
      const UF_long panel_local = panel_offset + local;
      update += rhs[local] *
                tail_values[panel_local * tail_len + offset];
    }
    x[col] -= update;
  }

  if (run_rows_out != NULL) {
    *run_rows_out = run_limit;
  }
  return 1;
}

static int kls_row_first_partial_apply_supernode_run_cached(
  const kls_row_first_block_context *ctx,
  kls_row_first_workspace *workspace,
  kls_row_first_entries *local_l_entries,
  const double *udiag_values,
  kls_row_first_partial_row *state,
  UF_long dep_begin,
  UF_long dep_end,
  UF_long *run_rows_out,
  int *used_panel_out) {
  KLS_SET_OPTIONAL_OUTPUT(run_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(used_panel_out, 0);
  if (ctx == NULL || workspace == NULL || local_l_entries == NULL ||
      udiag_values == NULL || state == NULL ||
      workspace->x == NULL || workspace->mark == NULL ||
      workspace->pattern == NULL || workspace->dep_heap == NULL ||
      dep_begin > dep_end || dep_end >= state->row ||
      state->row >= ctx->nk ||
      workspace->mark[dep_begin] != state->generation) {
    return -1;
  }
  const kls_row_first_supernode_panel_cache *cache =
    workspace->supernode_panel_cache;
  if (cache == NULL || !cache->enabled || cache->panel_id_by_row == NULL ||
      dep_begin >= cache->row_count) {
    return 0;
  }
  const UF_long panel = cache->panel_id_by_row[dep_begin];
  if (panel == KLS_KLU_EMPTY || panel >= cache->panel_count) {
    return 0;
  }
  if (cache->active == NULL || cache->active[panel] == 0u) {
    return 0;
  }
  const UF_long start = cache->start[panel];
  const UF_long width = cache->width[panel];
  const UF_long tail_len = cache->tail_len[panel];
  if (width <= 1u || dep_begin < start || dep_end < dep_begin ||
      width - 1u > UF_long_max - start) {
    return 0;
  }
  const UF_long panel_end = start + width - 1u;
  const UF_long effective_dep_end =
    dep_end < panel_end ? dep_end : panel_end;
  if (effective_dep_end < dep_begin ||
      effective_dep_end >= cache->row_count) {
    return 0;
  }
  const UF_long panel_offset = dep_begin - start;
  const UF_long available_end = effective_dep_end - start + 1u;
  const UF_long run_limit = effective_dep_end - dep_begin + 1u;
  if (available_end <= panel_offset ||
      available_end > width ||
      run_limit <= 1u ||
      run_limit != available_end - panel_offset) {
    return 0;
  }
  const UF_long dense_suffix_len = width - available_end;
  if (tail_len > UF_long_max - dense_suffix_len) {
    return 0;
  }

  double *x = workspace->x;
  unsigned int *mark = workspace->mark;
  UF_long *pattern = workspace->pattern;
  UF_long *dep_heap = workspace->dep_heap;
  const unsigned int generation = state->generation;
  const double *dense_panel =
    cache->dense_values + cache->dense_begin[panel];
  const double *tail_values =
    tail_len > 0u
      ? cache->tail_values + cache->tail_value_begin[panel]
      : NULL;
  const UF_long *tail_cols =
    tail_len > 0u ? cache->tail_cols + cache->tail_begin[panel] : NULL;
  if (dense_panel == NULL ||
      (tail_len > 0u && (tail_values == NULL || tail_cols == NULL))) {
    return 0;
  }

#ifdef KLS_HAVE_CBLAS
  const int cblas_status =
    kls_row_first_partial_apply_supernode_run_cached_cblas(
      ctx, workspace, local_l_entries, state, dep_begin, start, width,
      panel_offset, available_end, run_limit, dense_suffix_len, tail_len,
      dense_panel, tail_values, tail_cols, run_rows_out);
  if (cblas_status > 0) {
    if (used_panel_out != NULL) {
      *used_panel_out = 1;
    }
    return 1;
  }
  if (cblas_status < 0) {
    return -1;
  }
#endif
  const int portable_status =
    kls_row_first_partial_apply_supernode_run_cached_portable(
      ctx, workspace, local_l_entries, state, dep_begin, start, width,
      panel_offset, available_end, run_limit, dense_suffix_len, tail_len,
      dense_panel, tail_values, tail_cols, run_rows_out);
  if (portable_status > 0) {
    if (used_panel_out != NULL) {
      *used_panel_out = 1;
    }
    return 1;
  }
  if (portable_status < 0) {
    return -1;
  }

  const UF_long workspace_len = dense_suffix_len + tail_len;
  if (!kls_row_first_workspace_reserve_supernode(workspace,
                                                 workspace_len)) {
    return 0;
  }
  double *trailing_workspace = workspace->supernode_workspace;
  if (workspace_len > 0u) {
    memset(trailing_workspace, 0,
           (size_t)workspace_len * sizeof(*trailing_workspace));
  }

  UF_long run_rows = 0;
  if (!kls_row_first_entries_reserve_append(local_l_entries, run_limit)) {
    return -1;
  }
  for (UF_long local = 0; local < run_limit; ++local) {
    const UF_long dep = dep_begin + local;
    const UF_long panel_local = panel_offset + local;
    if (local > 0u) {
      if (mark[dep] != generation ||
          !kls_row_first_heap_consume_ready_dep(
            dep_heap, &state->dep_heap_size, dep)) {
        return -1;
      }
    }
    const double dep_pivot = udiag_values[dep];
    if (dep_pivot == 0.0) {
      return -1;
    }
    const double lij = x[dep] / dep_pivot;
    if (!kls_row_first_entries_append_reserved(
          local_l_entries, state->row, dep, lij)) {
      return -1;
    }
    x[dep] = 0.0;

    const double *dense_row = dense_panel + panel_local * width;
    for (UF_long target = panel_local + 1u; target < width; ++target) {
      const double value = dense_row[target];
      if (target < available_end) {
        const UF_long col = start + target;
        if (mark[col] != generation) {
          mark[col] = generation;
          pattern[state->pattern_count++] = col;
          x[col] = 0.0;
        }
        x[col] -= lij * value;
      } else {
        trailing_workspace[target - available_end] += lij * value;
      }
    }
    for (UF_long offset = 0; offset < tail_len; ++offset) {
      trailing_workspace[dense_suffix_len + offset] +=
        lij * tail_values[panel_local * tail_len + offset];
    }
    run_rows++;
  }

  for (UF_long offset = 0; offset < dense_suffix_len; ++offset) {
    const UF_long col = start + available_end + offset;
    if (mark[col] != generation) {
      mark[col] = generation;
      pattern[state->pattern_count++] = col;
      if (col < state->row) {
        kls_row_first_heap_push(dep_heap, &state->dep_heap_size, col);
      }
      x[col] = 0.0;
    }
    x[col] -= trailing_workspace[offset];
  }
  for (UF_long offset = 0; offset < tail_len; ++offset) {
    const UF_long col = tail_cols[offset];
    if (mark[col] != generation) {
      mark[col] = generation;
      pattern[state->pattern_count++] = col;
      if (col < state->row) {
        kls_row_first_heap_push(dep_heap, &state->dep_heap_size, col);
      }
      x[col] = 0.0;
    }
    x[col] -= trailing_workspace[dense_suffix_len + offset];
  }

  if (run_rows_out != NULL) {
    *run_rows_out = run_rows;
  }
  if (used_panel_out != NULL) {
    *used_panel_out = 1;
  }
  return 1;
}

/* Callers supply runs already identified by the supernode machinery. */
static int kls_row_first_partial_apply_supernode_run_compact(
  const kls_row_first_block_context *ctx,
  kls_row_first_workspace *workspace,
  kls_row_first_entries *local_l_entries,
  const kls_row_first_entries *published_u_entries,
  const double *udiag_values,
  const UF_long *u_row_ptr,
  const UF_long *u_row_end,
  kls_row_first_partial_row *state,
  UF_long dep_begin,
  UF_long dep_end,
  UF_long *run_rows_out) {
  KLS_SET_OPTIONAL_OUTPUT(run_rows_out, 0);
  if (ctx == NULL || workspace == NULL || local_l_entries == NULL ||
      published_u_entries == NULL || udiag_values == NULL ||
      u_row_ptr == NULL || u_row_end == NULL || state == NULL ||
      workspace->x == NULL || workspace->mark == NULL ||
      workspace->pattern == NULL || workspace->dep_heap == NULL ||
      dep_begin > dep_end || dep_end >= state->row ||
      state->row >= ctx->nk ||
      workspace->mark[dep_begin] != state->generation) {
    return -1;
  }
  const UF_long run_limit = dep_end - dep_begin + 1u;
  if (run_limit <= 1u) {
    return 0;
  }
  if (u_row_ptr[dep_end] > u_row_end[dep_end] ||
      u_row_end[dep_end] > published_u_entries->count) {
    return -1;
  }
  const UF_long tail_begin = u_row_ptr[dep_end];
  const UF_long tail_len = u_row_end[dep_end] - tail_begin;

  for (UF_long local = 0; local < run_limit; ++local) {
    const UF_long dep = dep_begin + local;
    if (udiag_values[dep] == 0.0 ||
        u_row_ptr[dep] > u_row_end[dep] ||
        u_row_end[dep] > published_u_entries->count) {
      return -1;
    }
    const UF_long row_begin = u_row_ptr[dep];
    const UF_long row_end = u_row_end[dep];
    const UF_long internal_len = dep_end - dep;
    if (row_end - row_begin < internal_len ||
        row_end - row_begin - internal_len != tail_len) {
      return 0;
    }
  }

  if (!kls_row_first_workspace_reserve_supernode(workspace, tail_len)) {
    return 0;
  }

  double *x = workspace->x;
  unsigned int *mark = workspace->mark;
  UF_long *pattern = workspace->pattern;
  UF_long *dep_heap = workspace->dep_heap;
  double *trailing_workspace = workspace->supernode_workspace;
  const unsigned int generation = state->generation;
  if (tail_len > 0u) {
    memset(trailing_workspace, 0,
           (size_t)tail_len * sizeof(*trailing_workspace));
  }

  UF_long run_rows = 0;
  if (!kls_row_first_entries_reserve_append(local_l_entries, run_limit)) {
    return -1;
  }
  for (UF_long local = 0; local < run_limit; ++local) {
    const UF_long dep = dep_begin + local;
    if (local > 0u) {
      if (mark[dep] != generation ||
          !kls_row_first_heap_consume_ready_dep(
            dep_heap, &state->dep_heap_size, dep)) {
        return -1;
      }
    }
    const double dep_pivot = udiag_values[dep];
    const double lij = x[dep] / dep_pivot;
    if (!kls_row_first_entries_append_reserved(
          local_l_entries, state->row, dep, lij)) {
      return -1;
    }
    x[dep] = 0.0;

    const UF_long row_begin = u_row_ptr[dep];
    const UF_long internal_len = dep_end - dep;
    for (UF_long offset = 0; offset < internal_len; ++offset) {
      const UF_long col = dep + offset + 1u;
      if (mark[col] != generation) {
        mark[col] = generation;
        pattern[state->pattern_count++] = col;
        x[col] = 0.0;
      }
      x[col] -= lij * published_u_entries->value[row_begin + offset];
    }
    for (UF_long offset = 0; offset < tail_len; ++offset) {
      trailing_workspace[offset] +=
        lij * published_u_entries->value[row_begin + internal_len + offset];
    }
    run_rows++;
  }

  for (UF_long offset = 0; offset < tail_len; ++offset) {
    const UF_long col = published_u_entries->col[tail_begin + offset];
    if (mark[col] != generation) {
      mark[col] = generation;
      pattern[state->pattern_count++] = col;
      if (col < state->row) {
        kls_row_first_heap_push(dep_heap, &state->dep_heap_size, col);
      }
      x[col] = 0.0;
    }
    x[col] -= trailing_workspace[offset];
  }

  if (run_rows_out != NULL) {
    *run_rows_out = run_rows;
  }
  return 1;
}

static int kls_row_first_partial_apply_supernode_run(
  const kls_row_first_block_context *ctx,
  kls_row_first_workspace *workspace,
  kls_row_first_entries *local_l_entries,
  const kls_row_first_entries *published_u_entries,
  const double *udiag_values,
  const UF_long *u_row_ptr,
  const UF_long *u_row_end,
  kls_row_first_partial_row *state,
  UF_long dep_begin,
  UF_long dep_end,
  UF_long *run_rows_out,
  int *used_panel_out,
  kls_row_first_row_stats *stats) {
  const int cached_status =
    kls_row_first_partial_apply_supernode_run_cached(
      ctx, workspace, local_l_entries, udiag_values, state, dep_begin,
      dep_end, run_rows_out, used_panel_out);
  if (cached_status > 0) {
    return 1;
  }
  if (cached_status < 0) {
    return 0;
  }
  if (workspace != NULL && workspace->supernode_panel_cache != NULL &&
      kls_row_first_supernode_panel_cache_append(
        workspace->supernode_panel_cache, published_u_entries,
        u_row_ptr, u_row_end, udiag_values, ctx->nk, dep_begin, dep_end,
        KLS_ROW_FIRST_MISSED_PANEL_CACHE_MAX_STORED_ENTRIES,
        stats)) {
    const int published_cached_status =
      kls_row_first_partial_apply_supernode_run_cached(
        ctx, workspace, local_l_entries, udiag_values, state, dep_begin,
        dep_end, run_rows_out, used_panel_out);
    if (published_cached_status > 0) {
      return 1;
    }
    if (published_cached_status < 0) {
      return 0;
    }
  }
  KLS_SET_OPTIONAL_OUTPUT(used_panel_out, 0);
  const int compact_status =
    kls_row_first_partial_apply_supernode_run_compact(
      ctx, workspace, local_l_entries, published_u_entries, udiag_values,
      u_row_ptr, u_row_end, state, dep_begin, dep_end, run_rows_out);
  if (compact_status > 0) {
    return 1;
  }
  if (compact_status < 0) {
    return 0;
  }
  return kls_row_first_partial_apply_supernode_run_scalar(
    ctx, workspace, local_l_entries, published_u_entries, udiag_values,
    u_row_ptr, u_row_end, state, dep_begin, dep_end, run_rows_out);
}

static int kls_row_first_partial_apply_ready_until(
  const kls_row_first_block_context *ctx,
  kls_row_first_workspace *workspace,
  kls_row_first_entries *local_l_entries,
  const kls_row_first_entries *published_u_entries,
  const double *udiag_values,
  const UF_long *u_row_ptr,
  const UF_long *u_row_end,
  const UF_long *supernode_start,
  const UF_long *supernode_end,
  UF_long supernode_valid_begin,
  const unsigned char *row_done,
  const UF_long *active_rank,
  UF_long ready_limit,
  kls_row_first_partial_row *state,
  UF_long stop_dep,
  UF_long *supernode_groups_out,
  UF_long *supernode_rows_out,
  UF_long *supernode_panel_groups_out,
  UF_long *supernode_panel_rows_out,
  kls_row_first_row_stats *stats,
  int *blocked_out);

static int kls_row_first_partial_apply_ready(
  const kls_row_first_block_context *ctx,
  kls_row_first_workspace *workspace,
  kls_row_first_entries *local_l_entries,
  const kls_row_first_entries *published_u_entries,
  const double *udiag_values,
  const UF_long *u_row_ptr,
  const UF_long *u_row_end,
  const UF_long *supernode_start,
  const UF_long *supernode_end,
  UF_long supernode_valid_begin,
  const unsigned char *row_done,
  const UF_long *active_rank,
  UF_long ready_limit,
  kls_row_first_partial_row *state,
  UF_long *supernode_groups_out,
  UF_long *supernode_rows_out,
  UF_long *supernode_panel_groups_out,
  UF_long *supernode_panel_rows_out,
  kls_row_first_row_stats *stats,
  int *blocked_out) {
  return kls_row_first_partial_apply_ready_until(
    ctx, workspace, local_l_entries, published_u_entries, udiag_values,
    u_row_ptr, u_row_end, supernode_start, supernode_end,
    supernode_valid_begin, row_done, active_rank, ready_limit, state,
    KLS_KLU_EMPTY, supernode_groups_out, supernode_rows_out,
    supernode_panel_groups_out, supernode_panel_rows_out, stats,
    blocked_out);
}

static int kls_row_first_partial_apply_ready_until(
  const kls_row_first_block_context *ctx,
  kls_row_first_workspace *workspace,
  kls_row_first_entries *local_l_entries,
  const kls_row_first_entries *published_u_entries,
  const double *udiag_values,
  const UF_long *u_row_ptr,
  const UF_long *u_row_end,
  const UF_long *supernode_start,
  const UF_long *supernode_end,
  UF_long supernode_valid_begin,
  const unsigned char *row_done,
  const UF_long *active_rank,
  UF_long ready_limit,
  kls_row_first_partial_row *state,
  UF_long stop_dep,
  UF_long *supernode_groups_out,
  UF_long *supernode_rows_out,
  UF_long *supernode_panel_groups_out,
  UF_long *supernode_panel_rows_out,
  kls_row_first_row_stats *stats,
  int *blocked_out) {
  KLS_SET_OPTIONAL_OUTPUT(blocked_out, 0);
  if (ctx == NULL || workspace == NULL || local_l_entries == NULL ||
      published_u_entries == NULL || udiag_values == NULL ||
      u_row_ptr == NULL || u_row_end == NULL || row_done == NULL ||
      state == NULL || workspace->x == NULL || workspace->mark == NULL ||
      workspace->pattern == NULL || workspace->dep_heap == NULL ||
      state->row >= ctx->nk) {
    return 0;
  }

  UF_long *dep_heap = workspace->dep_heap;
  unsigned int *mark = workspace->mark;
  const UF_long i = state->row;
  while (state->dep_heap_size > 0u) {
    const UF_long peek = dep_heap[0];
    if (stop_dep != KLS_KLU_EMPTY && peek >= stop_dep) {
      return 1;
    }
    if (!kls_row_first_dependency_ready(row_done, active_rank,
                                        ready_limit, peek)) {
      if (blocked_out != NULL) {
        *blocked_out = 1;
      }
      return 1;
    }
    const UF_long dep =
      kls_row_first_heap_pop(dep_heap, &state->dep_heap_size);
    if (dep >= i || mark[dep] != state->generation) {
      return 0;
    }
    UF_long run_end =
      kls_row_first_ready_supernode_end(
        supernode_start, supernode_end, row_done, active_rank, ctx->nk,
        supernode_valid_begin, ready_limit, i, dep);
    if (stop_dep != KLS_KLU_EMPTY && run_end >= stop_dep) {
      run_end = stop_dep > 0u ? stop_dep - 1u : dep;
    }
    UF_long run_rows = 1;
    int used_panel = 0;
    if (run_end > dep) {
      if (!kls_row_first_partial_apply_supernode_run(
            ctx, workspace, local_l_entries, published_u_entries,
            udiag_values, u_row_ptr, u_row_end, state, dep, run_end,
            &run_rows, &used_panel, stats)) {
        return 0;
      }
    } else {
      UF_long cached_run_end =
        kls_row_first_ready_cached_panel_end(
          workspace->supernode_panel_cache, row_done, active_rank,
          ctx->nk, supernode_valid_begin, ready_limit, i, dep);
      if (stop_dep != KLS_KLU_EMPTY && cached_run_end >= stop_dep) {
        cached_run_end = stop_dep > 0u ? stop_dep - 1u : dep;
      }
      if (cached_run_end > dep) {
        const int cached_status =
          kls_row_first_partial_apply_supernode_run_cached(
            ctx, workspace, local_l_entries, udiag_values, state, dep,
            cached_run_end, &run_rows, &used_panel);
        if (cached_status < 0) {
          return 0;
        }
      }
      if (run_rows <= 1u &&
          !kls_row_first_partial_apply_one_dep(
            ctx, workspace, local_l_entries, published_u_entries,
            udiag_values, u_row_ptr, u_row_end, state, dep)) {
        return 0;
      }
    }
    if (run_rows > 1u) {
      if (supernode_groups_out != NULL) {
        (*supernode_groups_out)++;
      }
      if (supernode_rows_out != NULL) {
        *supernode_rows_out += run_rows;
      }
      if (used_panel) {
        if (supernode_panel_groups_out != NULL) {
          (*supernode_panel_groups_out)++;
        }
        if (supernode_panel_rows_out != NULL) {
          *supernode_panel_rows_out += run_rows;
        }
      }
    }
  }
  return 1;
}

/* test-only fault injection: treat every k-th pipeline commit row's
   diagonal as unusable so the separator pipeline's pivot tail/restart
   machinery is exercised deterministically on any ordering.  The
   injected pivot is the row's legitimate strongest candidate, so the
   factorization stays numerically valid (residuals police it); row-id
   modulo keeps the injection independent of worker interleaving.
   Zero cost when the env is unset. */
static UF_long kls_row_first_debug_force_pivot_every(void) {
  static _Thread_local UF_long cached = -1;
  if (cached < 0) {
    const char *value = getenv("KLS_DEBUG_FORCE_ROW_PIPELINE_PIVOT");
    cached = 0;
    if (value != NULL && value[0] != '\0') {
      char *end = NULL;
      const long parsed = strtol(value, &end, 10);
      if (end != value && (end == NULL || *end == '\0') && parsed > 0) {
        cached = (UF_long)parsed;
      }
    }
  }
  return cached;
}

static int kls_row_first_partial_finish(
  const kls_row_first_block_context *ctx,
  kls_row_first_workspace *workspace,
  kls_row_first_entries *local_u_entries,
  kls_row_first_partial_row *state,
  int allow_pivot,
  int *pivot_needed_out,
  int *pivoted_out,
  UF_long *selected_col_out,
  int *selected_separator_exact_out,
  int *selected_separator_extent_out) {
  KLS_SET_OPTIONAL_OUTPUT(pivot_needed_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(pivoted_out, 0);
  if (selected_col_out != NULL) {
    *selected_col_out = KLS_KLU_EMPTY;
  }
  KLS_SET_OPTIONAL_OUTPUT(selected_separator_exact_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(selected_separator_extent_out, 0);
  if (ctx == NULL || workspace == NULL || local_u_entries == NULL ||
      state == NULL || workspace->x == NULL || workspace->mark == NULL ||
      workspace->pattern == NULL || workspace->dep_heap == NULL ||
      state->row >= ctx->nk) {
    return 0;
  }

  double *x = workspace->x;
  UF_long *pattern = workspace->pattern;
  UF_long *dep_heap = workspace->dep_heap;
  const UF_long i = state->row;
  double pivot =
    workspace->mark[i] == state->generation ? x[i] : 0.0;
  kls_row_first_pivot_choice pivot_choice;
  int pivot_needed = 0;
  UF_long selected_col = KLS_KLU_EMPTY;
  int selected_separator_exact = 0;
  int selected_separator_extent = 0;
  if (!kls_row_first_collect_pivot_choice(
        ctx, pattern, state->pattern_count, x, i, NULL, -1,
        &pivot_choice) ||
      !kls_row_first_select_pivot_exchange(
        &pivot_choice, pivot, ctx->tol, &pivot_needed, &selected_col,
        &selected_separator_exact, &selected_separator_extent)) {
    if (pivot_needed_out != NULL) {
      *pivot_needed_out = 1;
    }
    return 0;
  }
  if (allow_pivot && !pivot_needed && ctx->use_separator_for_block) {
    const UF_long force_every = kls_row_first_debug_force_pivot_every();
    if (force_every > 0 && i % force_every == 0) {
      const UF_long alt_col =
        pivot_choice.scope_enabled
          ? pivot_choice.scoped_col
          : (pivot_choice.limit_enabled ? pivot_choice.limited_col
                                        : pivot_choice.global_col);
      const double alt_abs =
        pivot_choice.scope_enabled
          ? pivot_choice.scoped_abs
          : (pivot_choice.limit_enabled ? pivot_choice.limited_abs
                                        : pivot_choice.global_abs);
      if (alt_col != KLS_KLU_EMPTY && alt_col != i && alt_abs > 0.0) {
        pivot_needed = 1;
        selected_col = alt_col;
        if (pivot_choice.scope_enabled) {
          selected_separator_exact = pivot_choice.scoped_col_exact;
          selected_separator_extent = !pivot_choice.scoped_col_exact;
        }
      }
    }
  }
  if (pivot_needed) {
    if (pivot_needed_out != NULL) {
      *pivot_needed_out = 1;
    }
    if (!allow_pivot || selected_col == KLS_KLU_EMPTY) {
      return 0;
    }
    const UF_long pivot_col = selected_col;
    const int pivot_separator_exact = selected_separator_exact;
    const int pivot_separator_extent = selected_separator_extent;
    const double swapped_pivot = x[selected_col];
    x[selected_col] = pivot;
    x[i] = swapped_pivot;
    workspace->mark[i] = state->generation;
    pivot = x[i];
    UF_long verify_selected_col = KLS_KLU_EMPTY;
    if (!kls_row_first_collect_pivot_choice(
          ctx, pattern, state->pattern_count, x, i, NULL, -1,
          &pivot_choice) ||
        !kls_row_first_select_pivot_exchange(
          &pivot_choice, pivot, ctx->tol, &pivot_needed,
          &verify_selected_col, NULL, NULL) ||
        pivot_needed) {
      return 0;
    }
    if (pivoted_out != NULL) {
      *pivoted_out = 1;
    }
    if (selected_col_out != NULL) {
      *selected_col_out = pivot_col;
    }
    if (selected_separator_exact_out != NULL) {
      *selected_separator_exact_out = pivot_separator_exact;
    }
    if (selected_separator_extent_out != NULL) {
      *selected_separator_extent_out = pivot_separator_extent;
    }
  }
  UF_long u_count = 0;
  for (UF_long p = 0; p < state->pattern_count; ++p) {
    const UF_long col = pattern[p];
    if (col > i) {
      dep_heap[u_count++] = col;
    }
  }
  if (u_count > 1u) {
    qsort(dep_heap, (size_t)u_count, sizeof(*dep_heap),
          kls_compare_uf_long);
  }
  state->pivot = pivot;
  for (UF_long p = 0; p < u_count; ++p) {
    const UF_long col = dep_heap[p];
    if (!kls_row_first_entries_append(local_u_entries, i, col, x[col])) {
      return 0;
    }
  }
  return 1;
}

typedef struct kls_row_first_pipeline_shared {
  const kls_row_first_block_context *ctx;
  const UF_long *row_order;
  const UF_long *active_rank;
  UF_long begin;
  UF_long end;
  kls_row_first_workspace *workspace;
  kls_row_first_entries *l_entries;
  kls_row_first_entries *u_entries;
  kls_row_first_entries *private_u_entries;
  const UF_long *private_u_row_ptr;
  const UF_long *private_u_row_end;
  UF_long *supernode_start;
  UF_long *supernode_end;
  UF_long supernode_valid_begin;
  kls_row_first_supernode_panel_cache *private_supernode_panel_cache;
  UF_long open_panel_start;
  UF_long open_panel_end;
  double *udiag_values;
  unsigned char *row_done;
  kls_row_first_row_stats *stats;
  kls_row_first_pipeline_worker *workers;

  kls_row_first_pipeline_batch_target *batch_targets;

  int worker_count;

  int batch_target_capacity;

  atomic_ulong next_pos;
  pthread_mutex_t lock;
  pthread_cond_t cond;
  UF_long completed_pos;

  UF_long order_epoch;
  UF_long pivot_tail_rows;
  UF_long pivot_restarts;

  int trace_enabled;
  int producer_batch_enabled;
  int active_catchup_batch_enabled;

  int pivot_supernode_rebase_enabled;
  int supernode_producer_batch_enabled;
  int supernode_producer_batch_active;
  int failed;
  kls_row_first_pipeline_failure failure_reason;
} kls_row_first_pipeline_shared;

struct kls_row_first_pipeline_worker {
  kls_row_first_pipeline_shared *shared;
  kls_row_first_workspace workspace;
  kls_row_first_entries l_entries;
  kls_row_first_entries u_entries;
  kls_row_first_partial_row *active_state;
  UF_long active_order_epoch;
  kls_row_first_partial_row state_storage;

  int external_update_active;
  UF_long rows;
  UF_long partial_rows;
  UF_long pipeline_wait_partial_rows;
  UF_long pipeline_wait_partial_deps;
  UF_long pipeline_supernode_update_groups;
  UF_long pipeline_supernode_update_rows;
  UF_long pipeline_supernode_panel_update_groups;
  UF_long pipeline_supernode_panel_update_rows;
};

static void kls_row_first_pipeline_worker_free(
  kls_row_first_pipeline_worker *worker) {
  if (worker == NULL) {
    return;
  }
  kls_row_first_workspace_free(&worker->workspace);
  kls_row_first_entries_free(&worker->l_entries);
  kls_row_first_entries_free(&worker->u_entries);

}

static int kls_row_first_pipeline_worker_scalar_batch_candidate(
  const kls_row_first_pipeline_shared *shared,
  const kls_row_first_pipeline_worker *worker,
  UF_long dep,
  int *not_root_out) {
  KLS_SET_OPTIONAL_OUTPUT(not_root_out, 0);
  if (shared == NULL || worker == NULL || worker->active_state == NULL ||
      worker->workspace.x == NULL || worker->workspace.mark == NULL ||
      worker->workspace.dep_heap == NULL ||
      shared->row_done == NULL || shared->ctx == NULL) {
    return 0;
  }
  const kls_row_first_partial_row *state = worker->active_state;
  if (worker->active_order_epoch != shared->order_epoch) {
    return 0;
  }

  if (worker->external_update_active) {
    return 0;
  }
  if (state->row >= shared->ctx->nk || dep >= state->row ||
      state->dep_heap_size == 0u) {
    return 0;
  }
  if (worker->workspace.mark[dep] != state->generation) {
    return 0;
  }
  if (worker->workspace.dep_heap[0] != dep) {
    KLS_SET_OPTIONAL_OUTPUT(not_root_out, 1);
    return 0;
  }
  if (!kls_row_first_dependency_ready(shared->row_done,
                                      shared->active_rank,
                                      shared->completed_pos, dep)) {
    return 0;
  }
  const UF_long run_end =
    kls_row_first_ready_supernode_end(
      shared->supernode_start, shared->supernode_end, shared->row_done,
      shared->active_rank, shared->ctx->nk, shared->supernode_valid_begin,
      shared->completed_pos, state->row, dep);
  if (run_end > dep) {
    return 0;
  }
  const UF_long cached_run_end =
    kls_row_first_ready_cached_panel_end(
      worker->workspace.supernode_panel_cache, shared->row_done,
      shared->active_rank, shared->ctx->nk, shared->supernode_valid_begin,
      shared->completed_pos, state->row, dep);
  if (cached_run_end > dep) {
    return 0;
  }
  return 1;
}

static int kls_row_first_pipeline_worker_supernode_batch_candidate(
  const kls_row_first_pipeline_shared *shared,
  const kls_row_first_pipeline_worker *worker,
  UF_long dep_begin,
  UF_long dep_end) {
  if (shared == NULL || worker == NULL || worker->active_state == NULL ||
      worker->workspace.x == NULL || worker->workspace.mark == NULL ||
      worker->workspace.dep_heap == NULL ||
      shared->row_done == NULL || shared->ctx == NULL ||
      dep_begin > dep_end || dep_end >= shared->ctx->nk) {
    return 0;
  }
  const kls_row_first_partial_row *state = worker->active_state;
  if (worker->active_order_epoch != shared->order_epoch) {
    return 0;
  }

  if (worker->external_update_active) {
    return 0;
  }
  if (state->row >= shared->ctx->nk || dep_end >= state->row ||
      state->dep_heap_size == 0u) {
    return 0;
  }
  if (worker->workspace.mark[dep_begin] != state->generation) {
    return 0;
  }
  if (worker->workspace.dep_heap[0] != dep_begin) {
    return 0;
  }
  if (!kls_row_first_dependency_ready(shared->row_done,
                                      shared->active_rank,
                                      shared->completed_pos,
                                      dep_begin)) {
    return 0;
  }
  const UF_long ready_end =
    kls_row_first_ready_supernode_end(
      shared->supernode_start, shared->supernode_end, shared->row_done,
      shared->active_rank, shared->ctx->nk, shared->supernode_valid_begin,
      shared->completed_pos, state->row, dep_begin);
  if (ready_end < dep_end) {
    return 0;
  }
  return 1;
}

static int kls_row_first_pipeline_apply_supernode_producer_batch(
  kls_row_first_pipeline_shared *shared,
  kls_row_first_pipeline_worker *current_worker,
  UF_long dep_begin,
  UF_long dep_end) {
  if (shared == NULL || current_worker == NULL ||
      !shared->producer_batch_enabled ||
      !shared->supernode_producer_batch_enabled ||
      shared->workers == NULL || shared->batch_targets == NULL ||
      shared->batch_target_capacity <= 1 || shared->ctx == NULL ||
      shared->u_entries == NULL || shared->workspace == NULL ||
      shared->workspace->u_row_ptr == NULL ||
      shared->workspace->u_row_end == NULL ||
      shared->udiag_values == NULL ||
      dep_begin >= dep_end || dep_end >= shared->ctx->nk) {
    return 1;
  }
  if (shared->supernode_start == NULL || shared->supernode_end == NULL ||
      dep_begin < shared->supernode_valid_begin ||
      shared->supernode_start[dep_begin] != dep_begin ||
      shared->supernode_end[dep_begin] < dep_end) {
    return 1;
  }

  const UF_long run_rows = dep_end - dep_begin + 1u;
  UF_long run_entries = 0u;
  for (UF_long dep = dep_begin; dep <= dep_end; ++dep) {
    const UF_long row_begin = shared->workspace->u_row_ptr[dep];
    const UF_long row_end = shared->workspace->u_row_end[dep];
    if (row_begin > row_end || row_end > shared->u_entries->count ||
        shared->udiag_values[dep] == 0.0) {
      return 1;
    }
    kls_row_first_stats_add(&run_entries, row_end - row_begin);
  }
  if (run_entries == 0u) {
    return 1;
  }

  UF_long target_count = 0u;
  for (int tid = 0; tid < shared->worker_count; ++tid) {
    kls_row_first_pipeline_worker *worker = shared->workers + tid;
    if (worker == current_worker) {
      continue;
    }
    if (!kls_row_first_pipeline_worker_supernode_batch_candidate(
          shared, worker, dep_begin, dep_end)) {
      continue;
    }
    if (target_count >= (UF_long)shared->batch_target_capacity) {
      return 0;
    }
    shared->batch_targets[target_count].worker = worker;
    shared->batch_targets[target_count].state = worker->active_state;
    shared->batch_targets[target_count].lij = 0.0;
    target_count++;
  }

  if (target_count < 2u) {
    return 1;
  }
  const UF_long saved_targets = target_count - 1u;
  const UF_long saved_stream =
    saved_targets != 0u && run_entries > UF_long_max / saved_targets
      ? UF_long_max
      : run_entries * saved_targets;
  if (saved_stream < KLS_ROW_FIRST_PRODUCER_BATCH_MIN_SAVED_STREAM) {
    return 1;
  }

  /*
   * The opt-in supernode producer can be expensive; claim target states and
   * run it outside the pipeline mutex while owner threads wait on the gate.
   */
  shared->supernode_producer_batch_active = 1;
  for (UF_long t = 0; t < target_count; ++t) {
    kls_row_first_pipeline_worker *worker =
      shared->batch_targets[t].worker;
    worker->external_update_active = 1;
  }

  pthread_mutex_unlock(&shared->lock);
  int ok = 1;
  for (UF_long t = 0; t < target_count; ++t) {
    kls_row_first_pipeline_worker *worker =
      shared->batch_targets[t].worker;
    kls_row_first_partial_row *state = shared->batch_targets[t].state;
    UF_long applied_rows = 0u;
    if (!kls_row_first_entries_reserve_append(&worker->l_entries,
                                              run_rows)) {
      ok = 0;
      break;
    }
    const int compact_status =
      kls_row_first_partial_apply_supernode_run_compact(
        shared->ctx, &worker->workspace, &worker->l_entries,
        shared->u_entries, shared->udiag_values,
        shared->workspace->u_row_ptr, shared->workspace->u_row_end,
        state, dep_begin, dep_end, &applied_rows);
    if (compact_status < 0) {
      ok = 0;
      break;
    }
    if (compact_status == 0 &&
        !kls_row_first_partial_apply_supernode_run_scalar(
          shared->ctx, &worker->workspace, &worker->l_entries,
          shared->u_entries, shared->udiag_values,
          shared->workspace->u_row_ptr, shared->workspace->u_row_end,
          state, dep_begin, dep_end, &applied_rows)) {
      ok = 0;
      break;
    }
    if (applied_rows != run_rows) {
      ok = 0;
      break;
    }
    worker->pipeline_supernode_update_groups++;
    worker->pipeline_supernode_update_rows += run_rows;
  }
  pthread_mutex_lock(&shared->lock);
  for (UF_long t = 0; t < target_count; ++t) {
    kls_row_first_pipeline_worker *worker =
      shared->batch_targets[t].worker;
    worker->external_update_active = 0;
  }
  shared->supernode_producer_batch_active = 0;
  pthread_cond_broadcast(&shared->cond);
  if (!ok) {
    return 0;
  }

  return 1;
}

static int kls_row_first_pipeline_apply_scalar_producer_targets(
  const kls_row_first_pipeline_shared *shared,
  kls_row_first_pipeline_batch_target *targets,
  UF_long target_count,
  UF_long dep,
  const UF_long *row_cols,
  const double *row_values,
  UF_long row_entries) {
  if (shared == NULL || targets == NULL || target_count == 0u ||
      shared->ctx == NULL || row_cols == NULL || row_values == NULL ||
      dep >= shared->ctx->nk) {
    return 0;
  }

  for (UF_long t = 0; t < target_count; ++t) {
    kls_row_first_pipeline_worker *worker = targets[t].worker;
    if (worker == NULL ||
        !kls_row_first_entries_reserve_append(&worker->l_entries, 1u)) {
      return 0;
    }
  }

  for (UF_long t = 0; t < target_count; ++t) {
    kls_row_first_pipeline_worker *worker = targets[t].worker;
    kls_row_first_partial_row *state = targets[t].state;
    if (worker == NULL || state == NULL ||
        worker->workspace.dep_heap == NULL ||
        worker->workspace.x == NULL) {
      return 0;
    }
    const UF_long popped =
      kls_row_first_heap_pop(worker->workspace.dep_heap,
                             &state->dep_heap_size);
    if (popped != dep ||
        !kls_row_first_entries_append_reserved(
          &worker->l_entries, state->row, dep, targets[t].lij)) {
      return 0;
    }
    worker->workspace.x[dep] = 0.0;
  }

  for (UF_long up = 0u; up < row_entries; ++up) {
    const UF_long col = row_cols[up];
    if (col <= dep || col >= shared->ctx->nk) {
      return 0;
    }
    const double uvalue = row_values[up];
    for (UF_long t = 0; t < target_count; ++t) {
      kls_row_first_pipeline_worker *worker = targets[t].worker;
      kls_row_first_partial_row *state = targets[t].state;
      if (worker == NULL || state == NULL ||
          worker->workspace.x == NULL ||
          worker->workspace.mark == NULL ||
          worker->workspace.pattern == NULL ||
          worker->workspace.dep_heap == NULL) {
        return 0;
      }
      double *x = worker->workspace.x;
      unsigned int *mark = worker->workspace.mark;
      UF_long *pattern = worker->workspace.pattern;
      UF_long *dep_heap = worker->workspace.dep_heap;
      if (mark[col] != state->generation) {
        mark[col] = state->generation;
        pattern[state->pattern_count++] = col;
        if (col < state->row) {
          kls_row_first_heap_push(dep_heap, &state->dep_heap_size, col);
        }
        x[col] = 0.0;
      }
      x[col] -= targets[t].lij * uvalue;
    }
  }

  return 1;
}

static int kls_row_first_pipeline_apply_producer_batch(
  kls_row_first_pipeline_shared *shared,
  kls_row_first_pipeline_worker *current_worker,
  UF_long dep) {
  if (shared == NULL || current_worker == NULL ||
      !shared->producer_batch_enabled ||
      shared->workers == NULL || shared->batch_targets == NULL ||
      shared->batch_target_capacity <= 1 || shared->ctx == NULL ||
      shared->u_entries == NULL || shared->workspace == NULL ||
      shared->workspace->u_row_ptr == NULL ||
      shared->workspace->u_row_end == NULL ||
      shared->udiag_values == NULL || dep >= shared->ctx->nk) {
    return 1;
  }
  const UF_long row_begin = shared->workspace->u_row_ptr[dep];
  const UF_long row_end = shared->workspace->u_row_end[dep];
  if (row_begin > row_end || row_end > shared->u_entries->count ||
      row_begin == row_end || shared->udiag_values[dep] == 0.0) {
    return 1;
  }
  const UF_long row_entries = row_end - row_begin;

  UF_long target_count = 0;
  for (int tid = 0; tid < shared->worker_count; ++tid) {
    kls_row_first_pipeline_worker *worker = shared->workers + tid;
    if (worker == current_worker) {
      continue;
    }
    int not_root = 0;
    if (!kls_row_first_pipeline_worker_scalar_batch_candidate(
          shared, worker, dep, &not_root)) {
      int retry_after_catchup = 0;
      if (shared->active_catchup_batch_enabled &&
          not_root &&
          worker->active_state != NULL &&
          worker->active_order_epoch == shared->order_epoch &&
          worker->workspace.mark != NULL &&
          worker->workspace.dep_heap != NULL &&
          dep < worker->active_state->row &&
          worker->workspace.mark[dep] ==
            worker->active_state->generation &&
          worker->active_state->dep_heap_size > 0u &&
          worker->workspace.dep_heap[0] < dep) {
        if (!kls_row_first_partial_apply_ready_until(
              shared->ctx, &worker->workspace, &worker->l_entries,
              shared->u_entries, shared->udiag_values,
              shared->workspace->u_row_ptr, shared->workspace->u_row_end,
              shared->supernode_start, shared->supernode_end,
              shared->supernode_valid_begin, shared->row_done,
              shared->active_rank, shared->completed_pos,
              worker->active_state, dep,
              &worker->pipeline_supernode_update_groups,
              &worker->pipeline_supernode_update_rows,
              &worker->pipeline_supernode_panel_update_groups,
              &worker->pipeline_supernode_panel_update_rows,
              shared->stats,
              NULL)) {
          return 0;
        }
        retry_after_catchup = 1;
      }
      if (!retry_after_catchup ||
          !kls_row_first_pipeline_worker_scalar_batch_candidate(
            shared, worker, dep, NULL)) {
        continue;
      }
    }
    if (target_count >= (UF_long)shared->batch_target_capacity) {
      return 0;
    }
    kls_row_first_partial_row *state = worker->active_state;
    const double lij = worker->workspace.x[dep] / shared->udiag_values[dep];
    shared->batch_targets[target_count].worker = worker;
    shared->batch_targets[target_count].state = state;
    shared->batch_targets[target_count].lij = lij;
    target_count++;
  }
  if (target_count < 2u) {
    return 1;
  }
  const UF_long saved_targets = target_count - 1u;
  const UF_long saved_stream =
    saved_targets != 0u && row_entries > UF_long_max / saved_targets
      ? UF_long_max
      : row_entries * saved_targets;
  if (saved_stream < KLS_ROW_FIRST_PRODUCER_BATCH_MIN_SAVED_STREAM) {
    return 1;
  }

  {
    if (target_count >
          (UF_long)(SIZE_MAX / sizeof(*shared->batch_targets)) ||
        row_entries > (UF_long)(SIZE_MAX / sizeof(*shared->u_entries->col)) ||
        row_entries > (UF_long)(SIZE_MAX /
                                sizeof(*shared->u_entries->value))) {
      return 0;
    }
    kls_row_first_pipeline_batch_target *detached_targets =
      (kls_row_first_pipeline_batch_target *)malloc(
        (size_t)target_count * sizeof(*detached_targets));
    if (detached_targets == NULL) {
      return 0;
    }
    UF_long *detached_cols =
      (UF_long *)malloc((size_t)row_entries *
                        sizeof(*detached_cols));
    double *detached_values =
      (double *)malloc((size_t)row_entries *
                       sizeof(*detached_values));
    if (detached_cols == NULL || detached_values == NULL) {
      free(detached_values);
      free(detached_cols);
      free(detached_targets);
      return 0;
    }
    memcpy(detached_cols, shared->u_entries->col + row_begin,
           (size_t)row_entries * sizeof(*detached_cols));
    memcpy(detached_values, shared->u_entries->value + row_begin,
           (size_t)row_entries * sizeof(*detached_values));
    for (UF_long t = 0; t < target_count; ++t) {
      kls_row_first_pipeline_worker *worker =
        shared->batch_targets[t].worker;
      detached_targets[t] = shared->batch_targets[t];
      if (worker != NULL) {
        worker->external_update_active = 1;
      }
    }
    pthread_mutex_unlock(&shared->lock);
    const int ok = kls_row_first_pipeline_apply_scalar_producer_targets(
      shared, detached_targets, target_count, dep,
      detached_cols, detached_values, row_entries);
    pthread_mutex_lock(&shared->lock);
    for (UF_long t = 0; t < target_count; ++t) {
      kls_row_first_pipeline_worker *worker = detached_targets[t].worker;
      if (worker != NULL) {
        worker->external_update_active = 0;
      }
    }
    pthread_cond_broadcast(&shared->cond);
    if (!ok) {
      free(detached_values);
      free(detached_cols);
      free(detached_targets);
      return 0;
    }
    memcpy(shared->batch_targets, detached_targets,
           (size_t)target_count * sizeof(*shared->batch_targets));
    free(detached_values);
    free(detached_cols);
    free(detached_targets);
  }

  return 1;
}

static void kls_row_first_pipeline_trace_print(
  const kls_row_first_pipeline_shared *shared,
  const char *event) {
  if (shared == NULL || !shared->trace_enabled || shared->ctx == NULL ||
      event == NULL) {
    return;
  }
  const UF_long completed =
    shared->completed_pos >= shared->begin
      ? shared->completed_pos - shared->begin
      : 0u;
  const UF_long total =
    shared->end >= shared->begin ? shared->end - shared->begin : 0u;
  fprintf(stderr,
          "KLS row-pipeline trace: event=%s completed=%" PRIu64
          "/%" PRIu64 " begin=%" PRIu64 " end=%" PRIu64 "\n",
          event, (uint64_t)completed, (uint64_t)total,
          (uint64_t)shared->begin, (uint64_t)shared->end);

}

static void kls_row_first_pipeline_mark_failed(
  kls_row_first_pipeline_shared *shared,
  kls_row_first_pipeline_failure reason) {
  if (shared == NULL) {
    return;
  }
  shared->failed = 1;
  if (shared->failure_reason == KLS_ROW_FIRST_PIPELINE_FAIL_NONE) {
    shared->failure_reason = reason;
  }
}

static void kls_row_first_pipeline_disable_prefix_supernodes(
  kls_row_first_pipeline_shared *shared) {
  if (shared == NULL) {
    return;
  }
  shared->supernode_start = NULL;
  shared->supernode_end = NULL;
  shared->open_panel_start = KLS_KLU_EMPTY;
  shared->open_panel_end = KLS_KLU_EMPTY;
}

static void kls_row_first_pipeline_disable_prefix_supernodes_after_pivot(
  kls_row_first_pipeline_shared *shared,
  UF_long pivot_row,
  UF_long pivot_col) {
  if (shared == NULL) {
    return;
  }
  shared->open_panel_start = KLS_KLU_EMPTY;
  shared->open_panel_end = KLS_KLU_EMPTY;
  if (shared->pivot_supernode_rebase_enabled &&
      shared->ctx != NULL &&
      shared->supernode_start != NULL &&
      shared->supernode_end != NULL) {
    const UF_long floor =
      pivot_row >= shared->ctx->nk ? shared->ctx->nk : pivot_row + 1u;
    if (floor > shared->supernode_valid_begin) {
      shared->supernode_valid_begin = floor;
    }
  } else {
    kls_row_first_pipeline_disable_prefix_supernodes(shared);
  }
  if (shared->private_supernode_panel_cache == NULL || shared->ctx == NULL) {
    return;
  }
  kls_row_first_supernode_panel_cache *cache =
    shared->private_supernode_panel_cache;
  if (cache->row_count == shared->ctx->nk &&
      cache->panel_id_by_row != NULL) {
    (void)kls_row_first_supernode_panel_cache_exchange_columns(
      cache, pivot_row, pivot_col);
  } else {
    kls_row_first_supernode_panel_cache_free(cache);
    if (!kls_row_first_supernode_panel_cache_init_rows(
          cache, shared->ctx->nk)) {
      shared->private_supernode_panel_cache = NULL;
    }
  }
}

static void kls_row_first_pipeline_rebuild_prefix_panel_cache(
  kls_row_first_pipeline_shared *shared,
  UF_long pivot_col_a,
  UF_long pivot_col_b) {
  if (shared == NULL || shared->ctx == NULL) {
    return;
  }
  if (shared->active_rank != NULL) {
    return;
  }
  kls_row_first_pipeline_disable_prefix_supernodes_after_pivot(
    shared, pivot_col_a, pivot_col_b);
}

static void kls_row_first_pipeline_publish_completed_panel(
  kls_row_first_pipeline_shared *shared,
  UF_long row) {
  if (shared == NULL ||
      shared->private_supernode_panel_cache == NULL ||
      shared->workspace == NULL ||
      shared->u_entries == NULL ||
      shared->workspace->u_row_ptr == NULL ||
      shared->workspace->u_row_end == NULL ||
      shared->udiag_values == NULL ||
      shared->ctx == NULL) {
    return;
  }
  kls_row_first_supernode_panel_cache *cache =
    shared->private_supernode_panel_cache;
  kls_row_first_supernode_panel_cache_publish_completed(
    cache, shared->u_entries,
    shared->workspace->u_row_ptr, shared->workspace->u_row_end,
    shared->udiag_values, shared->ctx->nk, shared->stats,
    &shared->open_panel_start, &shared->open_panel_end, row);
}

static void *kls_row_first_pipeline_worker_main(void *arg) {
  kls_row_first_pipeline_worker *worker =
    (kls_row_first_pipeline_worker *)arg;
  if (worker == NULL || worker->shared == NULL) {
    return NULL;
  }
  kls_row_first_pipeline_shared *shared = worker->shared;
  if (shared->ctx == NULL || shared->row_order == NULL ||
      shared->workspace == NULL || shared->l_entries == NULL ||
      shared->u_entries == NULL || shared->private_u_entries == NULL ||
      shared->private_u_row_ptr == NULL ||
      shared->private_u_row_end == NULL ||
      shared->supernode_start == NULL ||
      shared->supernode_end == NULL ||
      shared->udiag_values == NULL ||
      shared->row_done == NULL || shared->begin > shared->end) {
    return NULL;
  }
  if (!kls_row_first_workspace_init(&worker->workspace,
                                    shared->ctx->nk)) {
    pthread_mutex_lock(&shared->lock);
    kls_row_first_pipeline_mark_failed(
      shared, KLS_ROW_FIRST_PIPELINE_FAIL_SETUP);
    pthread_cond_broadcast(&shared->cond);
    pthread_mutex_unlock(&shared->lock);
    return NULL;
  }
  worker->workspace.supernode_panel_cache =
    shared->private_supernode_panel_cache;

  for (;;) {
    UF_long pos = KLS_KLU_EMPTY;

    {
      pos = (UF_long)atomic_fetch_add_explicit(&shared->next_pos, 1ul,
                                               memory_order_acq_rel);
    }
    UF_long row = KLS_KLU_EMPTY;
    UF_long adopted_order_epoch = 0u;
    kls_row_first_partial_row adopted_state;
    memset(&adopted_state, 0, sizeof(adopted_state));
    adopted_state.row = UF_long_max;
    int adopted_lookahead = 0;
    if (pos >= shared->end) {
      break;
    }
    if (pos < shared->begin) {
      pthread_mutex_lock(&shared->lock);
      kls_row_first_pipeline_mark_failed(
        shared, KLS_ROW_FIRST_PIPELINE_FAIL_DEPENDENCY);
      pthread_cond_broadcast(&shared->cond);
      pthread_mutex_unlock(&shared->lock);
      break;
    }

    pthread_mutex_lock(&shared->lock);
    while (!shared->failed &&
           shared->supernode_producer_batch_active) {
      if (pthread_cond_wait(&shared->cond, &shared->lock) != 0) {
        kls_row_first_pipeline_mark_failed(
          shared, KLS_ROW_FIRST_PIPELINE_FAIL_THREAD);
        pthread_cond_broadcast(&shared->cond);
        break;
      }
    }
    if (!shared->failed) {
      {
        adopted_lookahead =
          0;
        if (!adopted_lookahead) {
          adopted_lookahead =
            0;
        }
      }
    }
    pthread_mutex_unlock(&shared->lock);
    if (!adopted_lookahead) {
      row = shared->row_order[pos];
    }
    int stop_worker = 0;
    int retry_row = 1;
    while (retry_row) {
      retry_row = 0;
      kls_row_first_partial_row state;
      memset(&state, 0, sizeof(state));
      state.row = UF_long_max;
      if (adopted_lookahead) {
        state = adopted_state;
      } else {
        worker->l_entries.count = 0;
        worker->u_entries.count = 0;
      }
      int blocked = 0;
      int begin_ok = 1;
      UF_long order_epoch = 0;

      pthread_mutex_lock(&shared->lock);
      if (shared->failed) {
        stop_worker = 1;
        begin_ok = 0;
      } else if (adopted_lookahead) {
        order_epoch = adopted_order_epoch;
        if (order_epoch != shared->order_epoch) {
          retry_row = 1;
          begin_ok = 0;
        }
      } else {
        order_epoch = shared->order_epoch;
        if (row >= shared->ctx->nk ||
            !kls_row_first_partial_row_begin(shared->ctx,
                                             &worker->workspace,
                                             &state, row) ||
            !kls_row_first_partial_apply_ready(
              shared->ctx, &worker->workspace, &worker->l_entries,
              shared->private_u_entries, shared->udiag_values,
              shared->private_u_row_ptr, shared->private_u_row_end,
              shared->supernode_start, shared->supernode_end,
              shared->supernode_valid_begin, shared->row_done,
              shared->active_rank, shared->begin, &state,
              &worker->pipeline_supernode_update_groups,
              &worker->pipeline_supernode_update_rows,
              &worker->pipeline_supernode_panel_update_groups,
              &worker->pipeline_supernode_panel_update_rows,
              shared->stats,
              &blocked)) {
          kls_row_first_pipeline_mark_failed(
            shared, KLS_ROW_FIRST_PIPELINE_FAIL_DEPENDENCY);
          pthread_cond_broadcast(&shared->cond);
          stop_worker = 1;
          begin_ok = 0;
        }
      }
      if (begin_ok) {
        worker->active_state = &state;
        worker->active_order_epoch = order_epoch;
      }
      pthread_mutex_unlock(&shared->lock);
      if (!begin_ok) {
        kls_row_first_partial_row_clear(&state, &worker->workspace);
        if (retry_row) {
          adopted_lookahead = 0;
          worker->l_entries.count = 0;
          worker->u_entries.count = 0;
          continue;
        }
        break;
      }
      adopted_lookahead = 0;
      worker->partial_rows++;

      pthread_mutex_lock(&shared->lock);
      int wait_partial = 0;
      while (!shared->failed && shared->completed_pos != pos) {
        if (shared->supernode_producer_batch_active ||
            worker->external_update_active) {
          if (pthread_cond_wait(&shared->cond, &shared->lock) != 0) {
            kls_row_first_pipeline_mark_failed(
              shared, KLS_ROW_FIRST_PIPELINE_FAIL_THREAD);
            pthread_cond_broadcast(&shared->cond);
            break;
          }
          continue;
        }
        if (order_epoch != shared->order_epoch) {
          retry_row = 1;
          break;
        }
        const UF_long before_l_count = worker->l_entries.count;
        if (!kls_row_first_partial_apply_ready(
              shared->ctx, &worker->workspace, &worker->l_entries,
              shared->u_entries, shared->udiag_values,
              shared->workspace->u_row_ptr,
              shared->workspace->u_row_end,
              shared->supernode_start, shared->supernode_end,
              shared->supernode_valid_begin, shared->row_done,
              shared->active_rank, shared->completed_pos,
              &state,
              &worker->pipeline_supernode_update_groups,
              &worker->pipeline_supernode_update_rows,
              &worker->pipeline_supernode_panel_update_groups,
              &worker->pipeline_supernode_panel_update_rows,
              shared->stats,
              NULL)) {
          kls_row_first_pipeline_mark_failed(
            shared, KLS_ROW_FIRST_PIPELINE_FAIL_DEPENDENCY);
          pthread_cond_broadcast(&shared->cond);
          break;
        }
        if (worker->l_entries.count > before_l_count) {
          wait_partial = 1;
          worker->pipeline_wait_partial_deps +=
            worker->l_entries.count - before_l_count;
          continue;
        }
        if (pthread_cond_wait(&shared->cond, &shared->lock) != 0) {
          kls_row_first_pipeline_mark_failed(
            shared, KLS_ROW_FIRST_PIPELINE_FAIL_THREAD);
          pthread_cond_broadcast(&shared->cond);
          break;
        }
      }
      if (!retry_row && !shared->failed &&
          order_epoch != shared->order_epoch) {
        retry_row = 1;
      }
      while (!retry_row && !shared->failed &&
             (shared->supernode_producer_batch_active ||
              worker->external_update_active)) {
        if (pthread_cond_wait(&shared->cond, &shared->lock) != 0) {
          kls_row_first_pipeline_mark_failed(
            shared, KLS_ROW_FIRST_PIPELINE_FAIL_THREAD);
          pthread_cond_broadcast(&shared->cond);
          break;
        }
      }
      if (wait_partial) {
        worker->pipeline_wait_partial_rows++;
      }
      if (retry_row) {
        worker->active_state = NULL;
        worker->active_order_epoch = 0u;
        pthread_mutex_unlock(&shared->lock);
        kls_row_first_partial_row_clear(&state, &worker->workspace);
        continue;
      }
      if (!shared->failed) {
        blocked = 0;
        int pivot_needed = 0;
        int pivoted = 0;
        UF_long selected_col = KLS_KLU_EMPTY;
        int selected_separator_exact = 0;
        int selected_separator_extent = 0;
        int dependencies_ok = 1;
        const int unlock_commit_drain =
          ((shared->supernode_start == NULL ||
            shared->supernode_end == NULL) &&
           shared->active_rank == NULL) ||
          (shared->active_rank != NULL &&
           shared->ctx->nk >=
             KLS_ROW_FIRST_PIPELINE_PREFIX_CACHE_REBUILD_MAX_ROWS);
        if (unlock_commit_drain) {
          /*
           * This row is now the commit cursor. Later rows cannot append to the
           * published U storage until this row finishes, so the dependency drain
           * can run without holding the pipeline mutex. Large active-rank phases
           * still let waiting workers touch the shared panel cache under the
           * mutex, so the unlocked cursor consumes only read-only supernode
           * metadata and scalar/compact published-U rows.
           */
          kls_row_first_workspace detached_workspace;
          kls_row_first_workspace *drain_workspace = &worker->workspace;
          kls_row_first_row_stats *drain_stats = shared->stats;
          if (shared->active_rank != NULL) {
            detached_workspace = worker->workspace;
            detached_workspace.supernode_panel_cache = NULL;
            drain_workspace = &detached_workspace;
            drain_stats = NULL;
          }
          pthread_mutex_unlock(&shared->lock);
          dependencies_ok =
            kls_row_first_partial_apply_ready(
              shared->ctx, drain_workspace, &worker->l_entries,
              shared->u_entries, shared->udiag_values,
              shared->workspace->u_row_ptr, shared->workspace->u_row_end,
              shared->supernode_start, shared->supernode_end,
              shared->supernode_valid_begin, shared->row_done,
              shared->active_rank, shared->ctx->nk, &state,
              &worker->pipeline_supernode_update_groups,
              &worker->pipeline_supernode_update_rows,
              &worker->pipeline_supernode_panel_update_groups,
              &worker->pipeline_supernode_panel_update_rows,
              drain_stats,
              &blocked);
          if (drain_workspace != &worker->workspace) {
            worker->workspace.supernode_workspace =
              drain_workspace->supernode_workspace;
            worker->workspace.supernode_workspace_capacity =
              drain_workspace->supernode_workspace_capacity;
          }
          pthread_mutex_lock(&shared->lock);
        } else {
          dependencies_ok =
            kls_row_first_partial_apply_ready(
              shared->ctx, &worker->workspace, &worker->l_entries,
              shared->u_entries, shared->udiag_values,
              shared->workspace->u_row_ptr, shared->workspace->u_row_end,
              shared->supernode_start, shared->supernode_end,
              shared->supernode_valid_begin, shared->row_done,
              shared->active_rank, shared->ctx->nk, &state,
              &worker->pipeline_supernode_update_groups,
              &worker->pipeline_supernode_update_rows,
              &worker->pipeline_supernode_panel_update_groups,
              &worker->pipeline_supernode_panel_update_rows,
              shared->stats,
              &blocked);
        }
        if (!dependencies_ok || shared->failed || blocked) {
          kls_row_first_pipeline_mark_failed(
            shared, KLS_ROW_FIRST_PIPELINE_FAIL_DEPENDENCY);
        } else if (!kls_row_first_partial_finish(
                     shared->ctx, &worker->workspace, &worker->u_entries,
                     &state, 1, &pivot_needed, &pivoted, &selected_col,
                     &selected_separator_exact,
                     &selected_separator_extent)) {
          kls_row_first_pipeline_mark_failed(
            shared, pivot_needed ? KLS_ROW_FIRST_PIPELINE_FAIL_PIVOT
                                 : KLS_ROW_FIRST_PIPELINE_FAIL_APPEND);
        } else {
          const UF_long l_need = worker->l_entries.count;
          const UF_long u_need = worker->u_entries.count;
          if (l_need > UF_long_max - shared->l_entries->count ||
              u_need > UF_long_max - shared->u_entries->count ||
              !kls_row_first_entries_reserve(
                shared->l_entries, shared->l_entries->count + l_need) ||
              !kls_row_first_entries_reserve(
                shared->u_entries, shared->u_entries->count + u_need)) {
            kls_row_first_pipeline_mark_failed(
              shared, KLS_ROW_FIRST_PIPELINE_FAIL_APPEND);
          } else if (pivoted &&
                     (selected_col == KLS_KLU_EMPTY ||
                      !kls_row_first_exchange_columns(
                        shared->u_entries, shared->ctx->q_order,
                        shared->ctx->col_pos, row, selected_col,
                        shared->ctx->nk, shared->ctx->n) ||
                      !kls_row_first_entries_swap_columns(
                        shared->private_u_entries, row, selected_col,
                        shared->ctx->nk))) {
            kls_row_first_pipeline_mark_failed(
              shared, KLS_ROW_FIRST_PIPELINE_FAIL_PIVOT);
          } else if (!kls_row_first_append_entries(shared->l_entries,
                                                   &worker->l_entries)) {
            kls_row_first_pipeline_mark_failed(
              shared, KLS_ROW_FIRST_PIPELINE_FAIL_APPEND);
          } else {
            shared->workspace->u_row_ptr[row] = shared->u_entries->count;
            if (!kls_row_first_append_entries(shared->u_entries,
                                              &worker->u_entries)) {
              kls_row_first_pipeline_mark_failed(
                shared, KLS_ROW_FIRST_PIPELINE_FAIL_APPEND);
            } else {
              shared->workspace->u_row_end[row] = shared->u_entries->count;
              shared->udiag_values[row] = state.pivot;
              if (pivoted) {
                shared->open_panel_start = KLS_KLU_EMPTY;
                shared->open_panel_end = KLS_KLU_EMPTY;
              }
              if (shared->active_rank != NULL) {
                shared->row_done[row] = 1u;
                if (pivoted) {
                  kls_row_first_pipeline_disable_prefix_supernodes_after_pivot(
                    shared, row, selected_col);
                } else {
                  kls_row_first_supernodes_publish_row(
                    shared->u_entries, shared->workspace->u_row_ptr,
                    shared->workspace->u_row_end, shared->row_done,
                    shared->ctx->nk, shared->supernode_valid_begin, row,
                    shared->supernode_start, shared->supernode_end);
                }
              } else {
                const UF_long merge_floor =
                  shared->begin > shared->supernode_valid_begin
                    ? shared->begin
                    : shared->supernode_valid_begin;
                kls_row_first_supernodes_publish_row(
                  shared->u_entries, shared->workspace->u_row_ptr,
                  shared->workspace->u_row_end, shared->row_done,
                  shared->ctx->nk, merge_floor, row,
                  shared->supernode_start, shared->supernode_end);
                shared->row_done[row] = 1u;
              }

              kls_row_first_pipeline_publish_completed_panel(shared, row);
              if (!pivoted) {

              }

              shared->completed_pos = pos + 1u;
              if (!pivoted && shared->trace_enabled) {

              }
              if (!pivoted &&
                  shared->supernode_start != NULL &&
                  row < shared->ctx->nk) {
                const UF_long supernode_dep =
                  shared->supernode_start[row];
                if (supernode_dep < row &&
                    supernode_dep >= shared->supernode_valid_begin &&
                    shared->supernode_end[supernode_dep] >= row &&
                    !kls_row_first_pipeline_apply_supernode_producer_batch(
                      shared, worker, supernode_dep, row)) {
                  kls_row_first_pipeline_mark_failed(
                    shared, KLS_ROW_FIRST_PIPELINE_FAIL_DEPENDENCY);
                }
              }
              if (!pivoted &&
                  !kls_row_first_pipeline_apply_producer_batch(
                    shared, worker, row)) {
                kls_row_first_pipeline_mark_failed(
                  shared, KLS_ROW_FIRST_PIPELINE_FAIL_DEPENDENCY);
              }
              if (pivoted) {
                if (shared->stats != NULL) {
                  shared->stats->dynamic_column_pivots++;
                  if (selected_separator_exact) {
                    shared->stats->separator_dynamic_column_pivots++;
                  } else if (selected_separator_extent) {
                    shared->stats->separator_extent_dynamic_column_pivots++;
                  }
                }
                shared->pivot_tail_rows += shared->end - pos;
                shared->pivot_restarts++;

                shared->order_epoch++;
                if (shared->active_rank == NULL) {
                  kls_row_first_pipeline_rebuild_prefix_panel_cache(
                    shared, row, selected_col);
                }

              }
              worker->rows++;
              if (shared->trace_enabled) {
                if (((shared->completed_pos - shared->begin) %
                       KLS_ROW_FIRST_PIPELINE_TRACE_INTERVAL) == 0u ||
                    shared->completed_pos == shared->end) {
                  kls_row_first_pipeline_trace_print(shared, "progress");
                }
              }
            }
          }
        }
        pthread_cond_broadcast(&shared->cond);
      }
      const int failed = shared->failed;
      worker->active_state = NULL;
      worker->active_order_epoch = 0u;
      pthread_mutex_unlock(&shared->lock);
      kls_row_first_partial_row_clear(&state, &worker->workspace);
      if (failed) {
        stop_worker = 1;
        break;
      }
    }
    if (stop_worker) {
      break;
    }
  }
  return NULL;
}

static int kls_row_first_run_parallel_pipeline_phase(
  const kls_row_first_block_context *ctx,
  const UF_long *row_order,
  const UF_long *active_rank,
  UF_long begin,
  UF_long end,
  UF_long supernode_valid_begin,
  int thread_count,
  kls_row_first_workspace *workspace,
  kls_row_first_entries *l_entries,
  kls_row_first_entries *u_entries,
  double *udiag_values,
  unsigned char *row_done,
  kls_row_first_row_stats *stats,
  UF_long *pipeline_rows_out,
  UF_long *thread_rows_out,
  UF_long *partial_rows_out,
  UF_long *wait_partial_rows_out,
  UF_long *wait_partial_deps_out,
  UF_long *supernode_update_groups_out,
  UF_long *supernode_update_rows_out,
  UF_long *supernode_panel_update_groups_out,
  UF_long *supernode_panel_update_rows_out,
  UF_long *completed_pos_out,
  int *pivot_tail_ready_out,
  UF_long *pivot_tail_rows_out,
  UF_long *pivot_restarts_out) {
  if (completed_pos_out != NULL) {
    *completed_pos_out = begin;
  }
  KLS_SET_OPTIONAL_OUTPUT(pivot_tail_ready_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(pivot_tail_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(pivot_restarts_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(supernode_panel_update_groups_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(supernode_panel_update_rows_out, 0);
  if (ctx == NULL || row_order == NULL || begin >= end ||
      thread_count <= 1 || workspace == NULL || l_entries == NULL ||
      u_entries == NULL || udiag_values == NULL || row_done == NULL) {
    return 0;
  }
  const UF_long pipeline_rows = end - begin;
  int worker_count = thread_count;
  if ((UF_long)worker_count > pipeline_rows) {
    worker_count = (int)pipeline_rows;
  }
  if (worker_count <= 0) {
    return 0;
  }

  pthread_t *threads =
    (pthread_t *)calloc((size_t)worker_count, sizeof(*threads));
  kls_row_first_pipeline_worker *workers =
    (kls_row_first_pipeline_worker *)
      calloc((size_t)worker_count, sizeof(*workers));

  const int batch_target_capacity = worker_count;
  kls_row_first_pipeline_batch_target *batch_targets =
    (kls_row_first_pipeline_batch_target *)
      calloc((size_t)batch_target_capacity, sizeof(*batch_targets));
  kls_row_first_entries private_u_entries;
  memset(&private_u_entries, 0, sizeof(private_u_entries));
  UF_long *private_u_row_ptr =
    (UF_long *)malloc((size_t)ctx->nk * sizeof(*private_u_row_ptr));
  UF_long *private_u_row_end =
    (UF_long *)malloc((size_t)ctx->nk * sizeof(*private_u_row_end));
  UF_long *supernode_start =
    (UF_long *)malloc((size_t)ctx->nk * sizeof(*supernode_start));
  UF_long *supernode_end =
    (UF_long *)malloc((size_t)ctx->nk * sizeof(*supernode_end));
  if ((((((((((threads == NULL || workers == NULL || batch_targets == NULL)))))))))) {
    free(threads);
    free(workers);

    free(batch_targets);
    free(private_u_row_ptr);
    free(private_u_row_end);
    free(supernode_start);
    free(supernode_end);
    return 0;
  }
  if (private_u_row_ptr == NULL || private_u_row_end == NULL ||
      supernode_start == NULL || supernode_end == NULL ||
      !kls_row_first_append_entries(&private_u_entries, u_entries)) {
    free(threads);
    free(workers);

    free(batch_targets);
    kls_row_first_entries_free(&private_u_entries);
    free(private_u_row_ptr);
    free(private_u_row_end);
    free(supernode_start);
    free(supernode_end);
    return 0;
  }
  memcpy(private_u_row_ptr, workspace->u_row_ptr,
         (size_t)ctx->nk * sizeof(*private_u_row_ptr));
  memcpy(private_u_row_end, workspace->u_row_end,
         (size_t)ctx->nk * sizeof(*private_u_row_end));
  kls_row_first_supernode_panel_cache private_supernode_panel_cache;
  memset(&private_supernode_panel_cache, 0,
         sizeof(private_supernode_panel_cache));
  if (active_rank != NULL) {
    kls_row_first_supernodes_reset(
      &private_u_entries, private_u_row_ptr, private_u_row_end, row_done,
      ctx->nk, ctx->nk, supernode_start, supernode_end);
    (void)kls_row_first_supernode_panel_cache_build(
      &private_supernode_panel_cache, &private_u_entries, private_u_row_ptr,
      private_u_row_end, udiag_values, row_done, supernode_start,
      supernode_end, ctx->nk, ctx->nk, stats);
  } else {
    kls_row_first_supernodes_reset_prefix(
      &private_u_entries, private_u_row_ptr, private_u_row_end, row_done,
      ctx->nk, begin, supernode_start, supernode_end);
    (void)kls_row_first_supernode_panel_cache_build(
      &private_supernode_panel_cache, &private_u_entries, private_u_row_ptr,
      private_u_row_end, udiag_values, row_done, supernode_start,
      supernode_end, ctx->nk, begin, stats);
  }

  kls_row_first_pipeline_shared shared;
  memset(&shared, 0, sizeof(shared));
  shared.ctx = ctx;
  shared.row_order = row_order;
  shared.active_rank = active_rank;
  shared.begin = begin;
  shared.end = end;
  shared.workspace = workspace;
  shared.l_entries = l_entries;
  shared.u_entries = u_entries;
  shared.private_u_entries = &private_u_entries;
  shared.private_u_row_ptr = private_u_row_ptr;
  shared.private_u_row_end = private_u_row_end;
  shared.supernode_start = supernode_start;
  shared.supernode_end = supernode_end;
  shared.supernode_valid_begin =
    supernode_valid_begin <= ctx->nk ? supernode_valid_begin : ctx->nk;
  shared.private_supernode_panel_cache =
    private_supernode_panel_cache.panel_id_by_row != NULL
      ? &private_supernode_panel_cache : NULL;
  shared.open_panel_start = KLS_KLU_EMPTY;
  shared.open_panel_end = KLS_KLU_EMPTY;
  shared.udiag_values = udiag_values;
  shared.row_done = row_done;
  shared.stats = stats;
  shared.completed_pos = begin;

  shared.order_epoch = 0;

  shared.trace_enabled =
    (kls_row_pipeline_trace_env_enabled());
  shared.workers = workers;

  shared.batch_targets = batch_targets;

  shared.worker_count = worker_count;

  shared.batch_target_capacity = batch_target_capacity;

  shared.producer_batch_enabled =
    !kls_row_pipeline_producer_batch_env_disabled();
  shared.active_catchup_batch_enabled =
    shared.producer_batch_enabled &&
    kls_row_pipeline_active_catchup_env_enabled();

  shared.pivot_supernode_rebase_enabled =
    !kls_row_pipeline_pivot_supernode_rebase_env_disabled();
  shared.supernode_producer_batch_enabled =
    shared.producer_batch_enabled &&
    kls_row_pipeline_supernode_producer_batch_env_enabled();

  {

  }

  atomic_init(&shared.next_pos, (unsigned long)begin);
  if (pthread_mutex_init(&shared.lock, NULL) != 0) {
    kls_row_first_supernode_panel_cache_free(
      &private_supernode_panel_cache);
    free(threads);
    free(workers);

    free(batch_targets);
    kls_row_first_entries_free(&private_u_entries);
    free(private_u_row_ptr);
    free(private_u_row_end);
    free(supernode_start);
    free(supernode_end);

    return 0;
  }
  if (pthread_cond_init(&shared.cond, NULL) != 0) {
    pthread_mutex_destroy(&shared.lock);
    kls_row_first_supernode_panel_cache_free(
      &private_supernode_panel_cache);
    free(threads);
    free(workers);

    free(batch_targets);
    kls_row_first_entries_free(&private_u_entries);
    free(private_u_row_ptr);
    free(private_u_row_end);
    free(supernode_start);
    free(supernode_end);

    return 0;
  }

  kls_row_first_pipeline_trace_print(&shared, "start");

  int created = 0;
  for (int tid = 0; tid < worker_count; ++tid) {
    workers[tid].shared = &shared;
    if (pthread_create(&threads[tid], NULL,
                       kls_row_first_pipeline_worker_main,
                       workers + tid) != 0) {
      pthread_mutex_lock(&shared.lock);
      shared.failed = 1;
      pthread_cond_broadcast(&shared.cond);
      pthread_mutex_unlock(&shared.lock);
      break;
    }
    created++;
  }
  for (int tid = 0; tid < created; ++tid) {
    (void)pthread_join(threads[tid], NULL);
  }

  UF_long processed_rows = 0;
  UF_long partial_rows = 0;
  UF_long wait_partial_rows = 0;
  UF_long wait_partial_deps = 0;
  UF_long supernode_update_groups = 0;
  UF_long supernode_update_rows = 0;
  UF_long supernode_panel_update_groups = 0;
  UF_long supernode_panel_update_rows = 0;
  UF_long active_threads = 0;
  for (int tid = 0; tid < created; ++tid) {
    processed_rows += workers[tid].rows;
    partial_rows += workers[tid].partial_rows;
    wait_partial_rows += workers[tid].pipeline_wait_partial_rows;
    wait_partial_deps += workers[tid].pipeline_wait_partial_deps;
    supernode_update_groups +=
      workers[tid].pipeline_supernode_update_groups;
    supernode_update_rows += workers[tid].pipeline_supernode_update_rows;
    supernode_panel_update_groups +=
      workers[tid].pipeline_supernode_panel_update_groups;
    supernode_panel_update_rows +=
      workers[tid].pipeline_supernode_panel_update_rows;
    if (workers[tid].rows > 0u) {
      active_threads++;
    }
  }
  const int ok = created == worker_count && !shared.failed &&
    shared.completed_pos == end && processed_rows == pipeline_rows;
  const int pivot_tail_ready = created == worker_count && shared.failed &&
    shared.failure_reason == KLS_ROW_FIRST_PIPELINE_FAIL_PIVOT &&
    shared.completed_pos >= begin &&
    shared.completed_pos < end &&
    processed_rows == shared.completed_pos - begin;
  if (shared.trace_enabled) {
    kls_row_first_pipeline_trace_print(
      &shared, ok ? "done" : (pivot_tail_ready ? "pivot-tail" : "failed"));
  }
  if (completed_pos_out != NULL) {
    *completed_pos_out = shared.completed_pos;
  }
  if (pivot_tail_ready_out != NULL) {
    *pivot_tail_ready_out = pivot_tail_ready ? 1 : 0;
  }
  if ((ok || pivot_tail_ready) && pipeline_rows_out != NULL) {
    *pipeline_rows_out = processed_rows;
  }
  if ((ok || pivot_tail_ready) && thread_rows_out != NULL) {
    *thread_rows_out = active_threads;
  }
  if ((ok || pivot_tail_ready) && partial_rows_out != NULL) {
    *partial_rows_out = partial_rows;
  }
  if ((ok || pivot_tail_ready) && wait_partial_rows_out != NULL) {
    *wait_partial_rows_out = wait_partial_rows;
  }
  if ((ok || pivot_tail_ready) && wait_partial_deps_out != NULL) {
    *wait_partial_deps_out = wait_partial_deps;
  }
  if ((ok || pivot_tail_ready) && supernode_update_groups_out != NULL) {
    *supernode_update_groups_out = supernode_update_groups;
  }
  if ((ok || pivot_tail_ready) && supernode_update_rows_out != NULL) {
    *supernode_update_rows_out = supernode_update_rows;
  }
  if ((ok || pivot_tail_ready) &&
      supernode_panel_update_groups_out != NULL) {
    *supernode_panel_update_groups_out = supernode_panel_update_groups;
  }
  if ((ok || pivot_tail_ready) &&
      supernode_panel_update_rows_out != NULL) {
    *supernode_panel_update_rows_out = supernode_panel_update_rows;
  }
  if ((ok || pivot_tail_ready) && pivot_tail_rows_out != NULL) {
    *pivot_tail_rows_out = shared.pivot_tail_rows;
  }
  if ((ok || pivot_tail_ready) && pivot_restarts_out != NULL) {
    *pivot_restarts_out = shared.pivot_restarts;
  }

  pthread_cond_destroy(&shared.cond);
  pthread_mutex_destroy(&shared.lock);
  for (int tid = 0; tid < worker_count; ++tid) {
    kls_row_first_pipeline_worker_free(workers + tid);
  }

  kls_row_first_supernode_panel_cache_free(&private_supernode_panel_cache);
  kls_row_first_entries_free(&private_u_entries);
  free(private_u_row_ptr);
  free(private_u_row_end);
  free(supernode_start);
  free(supernode_end);
  free(threads);
  free(workers);

  free(batch_targets);
  return ok;
}

static int kls_row_first_run_restartable_pipeline_suffix(
  const kls_row_first_block_context *ctx,
  const UF_long *row_order,
  const UF_long *active_rank,
  UF_long begin,
  UF_long end,
  int thread_count,
  kls_row_first_workspace *workspace,
  kls_row_first_entries *l_entries,
  kls_row_first_entries *u_entries,
  double *udiag_values,
  unsigned char *row_done,
  kls_row_first_row_stats *stats,
  UF_long *pipeline_rows_out,
  UF_long *thread_rows_out,
  UF_long *partial_rows_out,
  UF_long *wait_partial_rows_out,
  UF_long *wait_partial_deps_out,
  UF_long *supernode_update_groups_out,
  UF_long *supernode_update_rows_out,
  UF_long *supernode_panel_update_groups_out,
  UF_long *supernode_panel_update_rows_out,
  UF_long *pivot_tail_rows_out,
  UF_long *pivot_restarts_out,
  UF_long *pivot_serial_rows_out,
  kls_row_first_supernode_panel_cache *pivot_serial_panel_cache,
  UF_long *pivot_serial_open_panel_start_io,
  UF_long *pivot_serial_open_panel_end_io) {
  KLS_SET_OPTIONAL_OUTPUT(pipeline_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(thread_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(partial_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(wait_partial_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(wait_partial_deps_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(supernode_update_groups_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(supernode_update_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(supernode_panel_update_groups_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(supernode_panel_update_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(pivot_tail_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(pivot_restarts_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(pivot_serial_rows_out, 0);
  if (ctx == NULL || row_order == NULL || begin > end ||
      thread_count <= 1 || workspace == NULL || l_entries == NULL ||
      u_entries == NULL || udiag_values == NULL || row_done == NULL) {
    return 0;
  }
  if (begin == end) {
    return 1;
  }

  UF_long cursor = begin;
  UF_long pipeline_rows = 0;
  UF_long max_threads = 0;
  UF_long partial_rows = 0;
  UF_long wait_partial_rows = 0;
  UF_long wait_partial_deps = 0;
  UF_long supernode_update_groups = 0;
  UF_long supernode_update_rows = 0;
  UF_long supernode_panel_update_groups = 0;
  UF_long supernode_panel_update_rows = 0;
  UF_long pivot_tail_rows = 0;
  UF_long pivot_restarts = 0;
  UF_long pivot_serial_rows = 0;
  UF_long supernode_valid_begin = 0u;
  while (cursor < end) {
    UF_long phase_rows = 0;
    UF_long phase_threads = 0;
    UF_long phase_partial_rows = 0;
    UF_long phase_wait_partial_rows = 0;
    UF_long phase_wait_partial_deps = 0;
    UF_long phase_supernode_groups = 0;
    UF_long phase_supernode_rows = 0;
    UF_long phase_supernode_panel_groups = 0;
    UF_long phase_supernode_panel_rows = 0;
    UF_long phase_pivot_tail_rows = 0;
    UF_long phase_pivot_restarts = 0;
    UF_long completed_pos = cursor;
    int pivot_tail_ready = 0;
    if (kls_row_first_run_parallel_pipeline_phase(
          ctx, row_order, active_rank, cursor, end, supernode_valid_begin,
          thread_count, workspace, l_entries, u_entries, udiag_values,
          row_done, stats, &phase_rows, &phase_threads,
          &phase_partial_rows, &phase_wait_partial_rows,
          &phase_wait_partial_deps, &phase_supernode_groups,
          &phase_supernode_rows, &phase_supernode_panel_groups,
          &phase_supernode_panel_rows, &completed_pos,
          &pivot_tail_ready, &phase_pivot_tail_rows,
          &phase_pivot_restarts)) {
      pipeline_rows += phase_rows;
      if (phase_threads > max_threads) {
        max_threads = phase_threads;
      }
      partial_rows += phase_partial_rows;
      wait_partial_rows += phase_wait_partial_rows;
      wait_partial_deps += phase_wait_partial_deps;
      supernode_update_groups += phase_supernode_groups;
      supernode_update_rows += phase_supernode_rows;
      supernode_panel_update_groups += phase_supernode_panel_groups;
      supernode_panel_update_rows += phase_supernode_panel_rows;
      pivot_tail_rows += phase_pivot_tail_rows;
      pivot_restarts += phase_pivot_restarts;
      cursor = end;
      break;
    }
    if (!pivot_tail_ready || completed_pos < cursor ||
        completed_pos >= end) {
      return 0;
    }
    pipeline_rows += phase_rows;
    if (phase_threads > max_threads) {
      max_threads = phase_threads;
    }
    partial_rows += phase_partial_rows;
    wait_partial_rows += phase_wait_partial_rows;
    wait_partial_deps += phase_wait_partial_deps;
    supernode_update_groups += phase_supernode_groups;
    supernode_update_rows += phase_supernode_rows;
    supernode_panel_update_groups += phase_supernode_panel_groups;
    supernode_panel_update_rows += phase_supernode_panel_rows;
    pivot_tail_rows += phase_pivot_tail_rows;
    pivot_restarts += phase_pivot_restarts;
    pivot_tail_rows += end - completed_pos;
    pivot_restarts++;
    const UF_long pivot_row = row_order[completed_pos];
    const UF_long pivots_before =
      stats != NULL ? stats->dynamic_column_pivots : 0u;
    kls_row_first_factor_row_failure_info pivot_failure;
    kls_row_first_factor_row_failure_clear(&pivot_failure, pivot_row, -1);
    if (pivot_row >= ctx->nk ||
        !kls_row_first_factor_one_row(
          ctx, workspace, l_entries, u_entries, udiag_values, pivot_row,
          NULL, -1, row_done, stats, &pivot_failure)) {
      if (kls_first_factor_trace_env_enabled()) {
        fprintf(stderr,
                "KLS first-factor row-pipeline pivot-tail failure:"
                " pos=%" PRIu64 " row=%" PRIu64
                " reason=%s related_col=%" PRIu64
                " actual_owner=%d active_rank=%d\n",
                (uint64_t)completed_pos, (uint64_t)pivot_row,
                kls_row_first_factor_row_failure_name(pivot_failure.reason),
                (uint64_t)pivot_failure.related_col,
                pivot_failure.actual_owner, active_rank != NULL ? 1 : 0);
      }
      return 0;
    }
    if (stats != NULL && stats->dynamic_column_pivots != pivots_before) {
      const UF_long floor =
        pivot_row >= ctx->nk ? ctx->nk : pivot_row + 1u;
      if (floor > supernode_valid_begin) {
        supernode_valid_begin = floor;
      }
      kls_row_first_supernode_panel_cache_reset_active(
        pivot_serial_panel_cache, workspace, ctx->nk,
        pivot_serial_open_panel_start_io,
        pivot_serial_open_panel_end_io);
    }
    kls_row_first_supernode_panel_cache_publish_completed(
      workspace->supernode_panel_cache != NULL
        ? pivot_serial_panel_cache : NULL,
      u_entries, workspace->u_row_ptr, workspace->u_row_end,
      udiag_values, ctx->nk, stats, pivot_serial_open_panel_start_io,
      pivot_serial_open_panel_end_io, pivot_row);
    pivot_serial_rows++;
    cursor = completed_pos + 1u;
  }

  if (pipeline_rows_out != NULL) {
    *pipeline_rows_out = pipeline_rows;
  }
  if (thread_rows_out != NULL) {
    *thread_rows_out = max_threads;
  }
  if (partial_rows_out != NULL) {
    *partial_rows_out = partial_rows;
  }
  if (wait_partial_rows_out != NULL) {
    *wait_partial_rows_out = wait_partial_rows;
  }
  if (wait_partial_deps_out != NULL) {
    *wait_partial_deps_out = wait_partial_deps;
  }
  if (supernode_update_groups_out != NULL) {
    *supernode_update_groups_out = supernode_update_groups;
  }
  if (supernode_update_rows_out != NULL) {
    *supernode_update_rows_out = supernode_update_rows;
  }
  if (supernode_panel_update_groups_out != NULL) {
    *supernode_panel_update_groups_out = supernode_panel_update_groups;
  }
  if (supernode_panel_update_rows_out != NULL) {
    *supernode_panel_update_rows_out = supernode_panel_update_rows;
  }
  if (pivot_tail_rows_out != NULL) {
    *pivot_tail_rows_out = pivot_tail_rows;
  }
  if (pivot_restarts_out != NULL) {
    *pivot_restarts_out = pivot_restarts;
  }
  if (pivot_serial_rows_out != NULL) {
    *pivot_serial_rows_out = pivot_serial_rows;
  }
  return 1;
}

static int kls_row_first_run_etree_prefactor_tail_phase(
  const kls_row_first_block_context *ctx,
  const UF_long *row_order,
  const unsigned char *active_mask,
  const UF_long *tail_parent,
  UF_long active_rows,
  UF_long supernode_valid_begin,
  int thread_count,
  kls_row_first_workspace *workspace,
  kls_row_first_entries *l_entries,
  kls_row_first_entries *u_entries,
  double *udiag_values,
  unsigned char *row_done,
  kls_row_first_row_stats *stats,
  UF_long *thread_rows_out,
  UF_long *prefactor_rows_out,
  UF_long *wait_partial_rows_out,
  UF_long *wait_partial_deps_out,
  UF_long *supernode_update_groups_out,
  UF_long *supernode_update_rows_out,
  UF_long *supernode_panel_update_groups_out,
  UF_long *supernode_panel_update_rows_out,
  UF_long *pivot_tail_rows_out,
  UF_long *pivot_restarts_out) {
  KLS_SET_OPTIONAL_OUTPUT(thread_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(prefactor_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(wait_partial_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(wait_partial_deps_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(supernode_update_groups_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(supernode_update_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(supernode_panel_update_groups_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(supernode_panel_update_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(pivot_tail_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(pivot_restarts_out, 0);
  if (ctx == NULL || row_order == NULL || active_mask == NULL ||
      tail_parent == NULL ||
      active_rows == 0u || thread_count <= 1 || workspace == NULL ||
      l_entries == NULL || u_entries == NULL || udiag_values == NULL ||
      row_done == NULL || workspace->u_row_ptr == NULL ||
      workspace->u_row_end == NULL) {
    return 0;
  }

  UF_long *active_rank =
    (UF_long *)malloc((size_t)ctx->nk * sizeof(*active_rank));
  if (active_rank == NULL) {
    return 0;
  }
  for (UF_long row = 0; row < ctx->nk; ++row) {
    active_rank[row] = KLS_KLU_EMPTY;
  }
  for (UF_long pos = 0; pos < active_rows; ++pos) {
    const UF_long row = row_order[pos];
    if (row >= ctx->nk || active_mask[row] == 0u ||
        active_rank[row] != KLS_KLU_EMPTY) {
      free(active_rank);
      return 0;
    }
    active_rank[row] = pos;
  }
  int active_tail_is_chain = 1;
  for (UF_long pos = 0; pos < active_rows; ++pos) {
    const UF_long row = row_order[pos];
    const UF_long parent = tail_parent[row];
    if (parent == row) {
      free(active_rank);
      return 0;
    }
    if (parent >= 0 && parent < ctx->nk && active_mask[parent] != 0u) {
      const UF_long parent_rank = active_rank[parent];
      if (parent_rank == KLS_KLU_EMPTY || parent_rank <= pos) {
        free(active_rank);
        return 0;
      }
      if (parent_rank != pos + 1u) {
        active_tail_is_chain = 0;
      }
    } else if (pos + 1u < active_rows) {
      active_tail_is_chain = 0;
    }
  }

  /* A width-one affected subtree has no exploitable row concurrency.  The
     direct topological walk avoids ready-queue, worker, and speculative-state
     machinery; the force switch is useful for new-architecture experiments,
     while the default gate is the exact structural certificate. */
  const int serial_etree_tail =
    getenv("KLS_ENABLE_FAST_REJECT_SERIAL_ETREE_TAIL") != NULL ||
    (getenv("KLS_DISABLE_FAST_REJECT_SERIAL_ETREE_TAIL") == NULL &&
     active_tail_is_chain);
  if (serial_etree_tail) {
    UF_long open_panel_start = KLS_KLU_EMPTY;
    UF_long open_panel_end = KLS_KLU_EMPTY;
    for (UF_long pos = 0; pos < active_rows; ++pos) {
      const UF_long row = row_order[pos];
      const UF_long pivots_before =
        stats != NULL ? stats->dynamic_column_pivots : 0u;
      if (row >= ctx->nk ||
          !kls_row_first_factor_one_row(
            ctx, workspace, l_entries, u_entries, udiag_values, row,
            NULL, -1, row_done, stats, NULL)) {
        free(active_rank);
        return 0;
      }
      if (stats != NULL && stats->dynamic_column_pivots != pivots_before) {
        const UF_long floor = row >= ctx->nk ? ctx->nk : row + 1u;
        if (floor > supernode_valid_begin) {
          supernode_valid_begin = floor;
        }
        if (workspace->supernode_panel_cache != NULL) {
          kls_row_first_supernode_panel_cache_reset_active(
            workspace->supernode_panel_cache, workspace, ctx->nk,
            &open_panel_start, &open_panel_end);
        }
      }
      kls_row_first_supernode_panel_cache_publish_completed(
        workspace->supernode_panel_cache, u_entries,
        workspace->u_row_ptr, workspace->u_row_end, udiag_values,
        ctx->nk, stats, &open_panel_start, &open_panel_end, row);
    }
    if (thread_rows_out != NULL) {
      *thread_rows_out = 1u;
    }
    if (prefactor_rows_out != NULL) {
      *prefactor_rows_out = active_rows;
    }
    free(active_rank);
    return 1;
  }

  UF_long cursor = 0;
  UF_long max_threads = 0;
  UF_long prefactor_rows = 0;
  UF_long wait_partial_rows = 0;
  UF_long wait_partial_deps = 0;
  UF_long supernode_update_groups = 0;
  UF_long supernode_update_rows = 0;
  UF_long supernode_panel_update_groups = 0;
  UF_long supernode_panel_update_rows = 0;
  UF_long pivot_tail_rows = 0;
  UF_long pivot_restarts = 0;
  UF_long open_panel_start = KLS_KLU_EMPTY;
  UF_long open_panel_end = KLS_KLU_EMPTY;
  while (cursor < active_rows) {
    UF_long phase_pipeline_rows = 0;
    UF_long phase_threads = 0;
    UF_long phase_prefactor_rows = 0;
    UF_long phase_wait_partial_rows = 0;
    UF_long phase_wait_partial_deps = 0;
    UF_long phase_supernode_update_groups = 0;
    UF_long phase_supernode_update_rows = 0;
    UF_long phase_supernode_panel_update_groups = 0;
    UF_long phase_supernode_panel_update_rows = 0;
    UF_long phase_pivot_tail_rows = 0;
    UF_long phase_pivot_restarts = 0;
    UF_long completed_pos = cursor;
    int pivot_tail_ready = 0;
    if (kls_row_first_run_parallel_pipeline_phase(
          ctx, row_order, active_rank, cursor, active_rows,
          supernode_valid_begin, thread_count, workspace, l_entries,
          u_entries, udiag_values, row_done, stats, &phase_pipeline_rows,
          &phase_threads, &phase_prefactor_rows, &phase_wait_partial_rows,
          &phase_wait_partial_deps, &phase_supernode_update_groups,
          &phase_supernode_update_rows, &phase_supernode_panel_update_groups,
          &phase_supernode_panel_update_rows, &completed_pos,
          &pivot_tail_ready, &phase_pivot_tail_rows,
          &phase_pivot_restarts)) {
      (void)phase_pipeline_rows;
      if (phase_threads > max_threads) {
        max_threads = phase_threads;
      }
      prefactor_rows += phase_prefactor_rows;
      wait_partial_rows += phase_wait_partial_rows;
      wait_partial_deps += phase_wait_partial_deps;
      supernode_update_groups += phase_supernode_update_groups;
      supernode_update_rows += phase_supernode_update_rows;
      supernode_panel_update_groups +=
        phase_supernode_panel_update_groups;
      supernode_panel_update_rows += phase_supernode_panel_update_rows;
      pivot_tail_rows += phase_pivot_tail_rows;
      pivot_restarts += phase_pivot_restarts;
      cursor = active_rows;
      break;
    }
    if (!pivot_tail_ready || completed_pos < cursor ||
        completed_pos >= active_rows) {
      free(active_rank);
      return 0;
    }
    (void)phase_pipeline_rows;
    if (phase_threads > max_threads) {
      max_threads = phase_threads;
    }
    prefactor_rows += phase_prefactor_rows;
    wait_partial_rows += phase_wait_partial_rows;
    wait_partial_deps += phase_wait_partial_deps;
    supernode_update_groups += phase_supernode_update_groups;
    supernode_update_rows += phase_supernode_update_rows;
    supernode_panel_update_groups +=
      phase_supernode_panel_update_groups;
    supernode_panel_update_rows += phase_supernode_panel_update_rows;
    pivot_tail_rows += phase_pivot_tail_rows;
    pivot_restarts += phase_pivot_restarts;
    pivot_tail_rows += active_rows - completed_pos;
    pivot_restarts++;

    const UF_long pivot_row = row_order[completed_pos];
    const UF_long pivots_before =
      stats != NULL ? stats->dynamic_column_pivots : 0u;
    if (pivot_row >= ctx->nk ||
        !kls_row_first_factor_one_row(
          ctx, workspace, l_entries, u_entries, udiag_values, pivot_row,
          NULL, -1, row_done, stats, NULL)) {
      free(active_rank);
      return 0;
    }
    if (stats != NULL && stats->dynamic_column_pivots != pivots_before) {
      const UF_long floor =
        pivot_row >= ctx->nk ? ctx->nk : pivot_row + 1u;
      if (floor > supernode_valid_begin) {
        supernode_valid_begin = floor;
      }
      if (workspace->supernode_panel_cache != NULL) {
        kls_row_first_supernode_panel_cache_reset_active(
          workspace->supernode_panel_cache, workspace, ctx->nk,
          &open_panel_start, &open_panel_end);
      }
    }
    kls_row_first_supernode_panel_cache_publish_completed(
      workspace->supernode_panel_cache, u_entries, workspace->u_row_ptr,
      workspace->u_row_end, udiag_values, ctx->nk, stats, &open_panel_start,
      &open_panel_end, pivot_row);
    cursor = completed_pos + 1u;
  }
  if (thread_rows_out != NULL) {
    *thread_rows_out = max_threads;
  }
  if (prefactor_rows_out != NULL) {
    *prefactor_rows_out = prefactor_rows;
  }
  if (wait_partial_rows_out != NULL) {
    *wait_partial_rows_out = wait_partial_rows;
  }
  if (wait_partial_deps_out != NULL) {
    *wait_partial_deps_out = wait_partial_deps;
  }
  if (supernode_update_groups_out != NULL) {
    *supernode_update_groups_out = supernode_update_groups;
  }
  if (supernode_update_rows_out != NULL) {
    *supernode_update_rows_out = supernode_update_rows;
  }
  if (supernode_panel_update_groups_out != NULL) {
    *supernode_panel_update_groups_out = supernode_panel_update_groups;
  }
  if (supernode_panel_update_rows_out != NULL) {
    *supernode_panel_update_rows_out = supernode_panel_update_rows;
  }
  if (pivot_tail_rows_out != NULL) {
    *pivot_tail_rows_out = pivot_tail_rows;
  }
  if (pivot_restarts_out != NULL) {
    *pivot_restarts_out = pivot_restarts;
  }
  free(active_rank);
  return 1;
}

static void kls_row_first_lu_seed_free(kls_row_first_lu_seed *seed) {
  if (seed == NULL) {
    return;
  }
  free(seed->row);
  free(seed->col);
  free(seed->numeric_value);
  free(seed->value);
  memset(seed, 0, sizeof(*seed));
}

static void kls_row_first_input_seed_free(kls_row_first_input_seed *seed) {
  if (seed == NULL) {
    return;
  }
  free(seed->row);
  free(seed->col);
  free(seed->pos);
  memset(seed, 0, sizeof(*seed));
}

static void kls_row_first_refactor_seed_clear(
  kls_row_first_refactor_seed *seed) {
  if (seed == NULL) {
    return;
  }
  kls_row_first_lu_seed_free(&seed->l);
  kls_row_first_lu_seed_free(&seed->u);
  kls_row_first_input_seed_free(&seed->input);
  seed->n = 0;
  seed->enabled = 0;
}

static int kls_row_first_refactor_seed_init(
  kls_row_first_refactor_seed *seed,
  UF_long n) {
  if (seed == NULL) {
    return 0;
  }
  memset(seed, 0, sizeof(*seed));
  seed->n = n;
  seed->enabled = 1;
  return 1;
}

static int kls_row_first_lu_seed_reserve(kls_row_first_lu_seed *seed,
                                         UF_long extra) {
  if (seed == NULL || extra == 0u) {
    return seed != NULL;
  }
  if (extra > UF_long_max - seed->count) {
    return 0;
  }
  const UF_long needed = seed->count + extra;
  if (needed <= seed->capacity) {
    return 1;
  }
  UF_long capacity = seed->capacity == 0u ? 64u : seed->capacity;
  while (capacity < needed) {
    const UF_long grown = 2u * capacity;
    if (grown <= capacity) {
      capacity = needed;
      break;
    }
    capacity = grown;
  }
  if (capacity > (UF_long)(SIZE_MAX / sizeof(*seed->row)) ||
      capacity > (UF_long)(SIZE_MAX / sizeof(*seed->numeric_value)) ||
      capacity > (UF_long)(SIZE_MAX / sizeof(*seed->value))) {
    return 0;
  }
  UF_long *row =
    (UF_long *)realloc(seed->row, (size_t)capacity * sizeof(*row));
  if (row == NULL) {
    return 0;
  }
  seed->row = row;
  UF_long *col =
    (UF_long *)realloc(seed->col, (size_t)capacity * sizeof(*col));
  if (col == NULL) {
    return 0;
  }
  seed->col = col;
  double **numeric_value =
    (double **)realloc(seed->numeric_value,
                       (size_t)capacity * sizeof(*numeric_value));
  if (numeric_value == NULL) {
    return 0;
  }
  seed->numeric_value = numeric_value;
  double *value =
    (double *)realloc(seed->value, (size_t)capacity * sizeof(*value));
  if (value == NULL) {
    return 0;
  }
  seed->value = value;
  seed->capacity = capacity;
  return 1;
}

static int kls_row_first_input_seed_reserve(kls_row_first_input_seed *seed,
                                            UF_long extra) {
  if (seed == NULL || extra == 0u) {
    return seed != NULL;
  }
  if (extra > UF_long_max - seed->count) {
    return 0;
  }
  const UF_long needed = seed->count + extra;
  if (needed <= seed->capacity) {
    return 1;
  }
  UF_long capacity = seed->capacity == 0u ? 64u : seed->capacity;
  while (capacity < needed) {
    const UF_long grown = 2u * capacity;
    if (grown <= capacity) {
      capacity = needed;
      break;
    }
    capacity = grown;
  }
  if (capacity > (UF_long)(SIZE_MAX / sizeof(*seed->row))) {
    return 0;
  }
  UF_long *row =
    (UF_long *)realloc(seed->row, (size_t)capacity * sizeof(*row));
  if (row == NULL) {
    return 0;
  }
  seed->row = row;
  UF_long *col =
    (UF_long *)realloc(seed->col, (size_t)capacity * sizeof(*col));
  if (col == NULL) {
    return 0;
  }
  seed->col = col;
  UF_long *pos =
    (UF_long *)realloc(seed->pos, (size_t)capacity * sizeof(*pos));
  if (pos == NULL) {
    return 0;
  }
  seed->pos = pos;
  seed->capacity = capacity;
  return 1;
}

static int kls_row_first_lu_seed_append(kls_row_first_lu_seed *seed,
                                        UF_long row,
                                        UF_long col,
                                        double *numeric_value,
                                        double value) {
  if (seed == NULL || numeric_value == NULL ||
      !kls_row_first_lu_seed_reserve(seed, 1u)) {
    return 0;
  }
  const UF_long p = seed->count++;
  seed->row[p] = row;
  seed->col[p] = col;
  seed->numeric_value[p] = numeric_value;
  seed->value[p] = value;
  return 1;
}

static int kls_row_first_input_seed_append(kls_row_first_input_seed *seed,
                                           UF_long row,
                                           UF_long col,
                                           UF_long pos) {
  if (seed == NULL || !kls_row_first_input_seed_reserve(seed, 1u)) {
    return 0;
  }
  const UF_long p = seed->count++;
  seed->row[p] = row;
  seed->col[p] = col;
  seed->pos[p] = pos;
  return 1;
}

static int kls_row_first_refactor_seed_capture_block(
  kls_row_first_refactor_seed *seed,
  const kls_solver *solver,
  const trilinos_klu_l_numeric *numeric,
  UF_long block,
  UF_long k1,
  UF_long nk,
  const UF_long *row_ptr,
  const UF_long *row_cols,
  const UF_long *row_input_pos,
  const UF_long *col_pos,
  const kls_row_first_entries *l_entries,
  const kls_row_first_entries *u_entries) {
  if (seed == NULL || !seed->enabled) {
    return 1;
  }
  if (solver == NULL || numeric == NULL || row_ptr == NULL ||
      (row_ptr[nk] > 0u && (row_cols == NULL || row_input_pos == NULL)) ||
      col_pos == NULL || l_entries == NULL || u_entries == NULL ||
      numeric->LUbx == NULL || numeric->Lip == NULL ||
      numeric->Llen == NULL || numeric->Uip == NULL ||
      numeric->Ulen == NULL || block >= solver->symbolic->nblocks ||
      k1 > seed->n || nk > seed->n - k1) {
    return 0;
  }

  Unit *lu = numeric->LUbx[block];
  if (lu == NULL) {
    return 0;
  }

  for (UF_long row = 0; row < nk; ++row) {
    for (UF_long p = row_ptr[row]; p < row_ptr[row + 1u]; ++p) {
      const UF_long oldcol = row_cols[p];
      if (oldcol >= seed->n || col_pos[oldcol] >= nk ||
          row_input_pos[p] >= solver->nnz) {
        return 0;
      }
      if (!kls_row_first_input_seed_append(&seed->input, k1 + row,
                                           k1 + col_pos[oldcol],
                                           row_input_pos[p])) {
        return 0;
      }
    }
  }

  UF_long *l_write =
    nk > 0u ? (UF_long *)calloc((size_t)nk, sizeof(*l_write)) : NULL;
  UF_long *u_write =
    nk > 0u ? (UF_long *)calloc((size_t)nk, sizeof(*u_write)) : NULL;
  if ((nk > 0u && l_write == NULL) || (nk > 0u && u_write == NULL)) {
    free(l_write);
    free(u_write);
    return 0;
  }

  for (UF_long p = 0; p < l_entries->count; ++p) {
    const UF_long row = l_entries->row[p];
    const UF_long col = l_entries->col[p];
    if (row >= nk || col >= nk || row <= col) {
      free(l_write);
      free(u_write);
      return 0;
    }
    const UF_long len = numeric->Llen[k1 + col];
    const UF_long offset = l_write[col]++;
    if (offset >= len) {
      free(l_write);
      free(u_write);
      return 0;
    }
    double *lx =
      (double *)(lu + numeric->Lip[k1 + col] +
                 kls_klu_units_for_indices(len));
    if (!kls_row_first_lu_seed_append(&seed->l, k1 + row, k1 + col,
                                      lx + offset, l_entries->value[p])) {
      free(l_write);
      free(u_write);
      return 0;
    }
  }

  for (UF_long p = 0; p < u_entries->count; ++p) {
    const UF_long row = u_entries->row[p];
    const UF_long col = u_entries->col[p];
    if (row >= nk || col >= nk || col <= row) {
      free(l_write);
      free(u_write);
      return 0;
    }
    const UF_long len = numeric->Ulen[k1 + col];
    const UF_long offset = u_write[col]++;
    if (offset >= len) {
      free(l_write);
      free(u_write);
      return 0;
    }
    double *ux =
      (double *)(lu + numeric->Uip[k1 + col] +
                 kls_klu_units_for_indices(len));
    if (!kls_row_first_lu_seed_append(&seed->u, k1 + row, k1 + col,
                                      ux + offset, u_entries->value[p])) {
      free(l_write);
      free(u_write);
      return 0;
    }
  }

  free(l_write);
  free(u_write);
  return 1;
}

static int kls_row_first_refactor_seed_install(
  kls_solver *solver,
  kls_row_first_refactor_seed *seed) {
  if (solver == NULL || seed == NULL || !seed->enabled ||
      seed->n != solver->n) {
    kls_row_first_refactor_seed_clear(seed);
    return 0;
  }
  const UF_long n = seed->n;
  UF_long *l_ptr = (UF_long *)calloc((size_t)n + 1u, sizeof(*l_ptr));
  UF_long *u_ptr = (UF_long *)calloc((size_t)n + 1u, sizeof(*u_ptr));
  UF_long *input_ptr =
    (UF_long *)calloc((size_t)n + 1u, sizeof(*input_ptr));
  if (l_ptr == NULL || u_ptr == NULL || input_ptr == NULL) {
    free(l_ptr);
    free(u_ptr);
    free(input_ptr);
    kls_row_first_refactor_seed_clear(seed);
    return 0;
  }

  for (UF_long p = 0; p < seed->l.count; ++p) {
    if (seed->l.row[p] >= n || seed->l.col[p] >= seed->l.row[p]) {
      free(l_ptr);
      free(u_ptr);
      free(input_ptr);
      kls_row_first_refactor_seed_clear(seed);
      return 0;
    }
    l_ptr[seed->l.row[p] + 1u]++;
  }
  for (UF_long p = 0; p < seed->u.count; ++p) {
    if (seed->u.row[p] >= n || seed->u.col[p] <= seed->u.row[p] ||
        seed->u.col[p] >= n) {
      free(l_ptr);
      free(u_ptr);
      free(input_ptr);
      kls_row_first_refactor_seed_clear(seed);
      return 0;
    }
    u_ptr[seed->u.row[p] + 1u]++;
  }
  for (UF_long p = 0; p < seed->input.count; ++p) {
    if (seed->input.row[p] >= n || seed->input.col[p] >= n ||
        seed->input.pos[p] >= solver->nnz) {
      free(l_ptr);
      free(u_ptr);
      free(input_ptr);
      kls_row_first_refactor_seed_clear(seed);
      return 0;
    }
    input_ptr[seed->input.row[p] + 1u]++;
  }
  for (UF_long row = 0; row < n; ++row) {
    l_ptr[row + 1u] += l_ptr[row];
    u_ptr[row + 1u] += u_ptr[row];
    input_ptr[row + 1u] += input_ptr[row];
  }

  const UF_long l_nnz = l_ptr[n];
  const UF_long u_nnz = u_ptr[n];
  const UF_long input_nnz = input_ptr[n];
  UF_long *l_cols = l_nnz > 0u
    ? (UF_long *)malloc((size_t)l_nnz * sizeof(*l_cols)) : NULL;
  double **l_values = l_nnz > 0u
    ? (double **)malloc((size_t)l_nnz * sizeof(*l_values)) : NULL;
  double *l_row_values = l_nnz > 0u
    ? (double *)malloc((size_t)l_nnz * sizeof(*l_row_values)) : NULL;
  UF_long *u_cols = u_nnz > 0u
    ? (UF_long *)malloc((size_t)u_nnz * sizeof(*u_cols)) : NULL;
  double **u_values = u_nnz > 0u
    ? (double **)malloc((size_t)u_nnz * sizeof(*u_values)) : NULL;
  double *u_row_values = kls_aligned_double_values(u_nnz);
  UF_long *input_cols = input_nnz > 0u
    ? (UF_long *)malloc((size_t)input_nnz * sizeof(*input_cols)) : NULL;
  UF_long *input_pos = input_nnz > 0u
    ? (UF_long *)malloc((size_t)input_nnz * sizeof(*input_pos)) : NULL;
  UF_long *next = n > 0u ? (UF_long *)malloc((size_t)n * sizeof(*next))
                         : NULL;
  if ((l_nnz > 0u &&
       (l_cols == NULL || l_values == NULL || l_row_values == NULL)) ||
      (u_nnz > 0u &&
       (u_cols == NULL || u_values == NULL || u_row_values == NULL)) ||
      (input_nnz > 0u && (input_cols == NULL || input_pos == NULL)) ||
      (n > 0u && next == NULL)) {
    free(l_ptr);
    free(u_ptr);
    free(input_ptr);
    free(l_cols);
    free(l_values);
    free(l_row_values);
    free(u_cols);
    free(u_values);
    free(u_row_values);
    free(input_cols);
    free(input_pos);
    free(next);
    kls_row_first_refactor_seed_clear(seed);
    return 0;
  }

  if (n > 0u) {
    memcpy(next, l_ptr, (size_t)n * sizeof(*next));
  }
  for (UF_long p = 0; p < seed->l.count; ++p) {
    const UF_long row = seed->l.row[p];
    const UF_long dst = next[row]++;
    l_cols[dst] = seed->l.col[p];
    l_values[dst] = seed->l.numeric_value[p];
    l_row_values[dst] = seed->l.value[p];
  }
  if (n > 0u) {
    memcpy(next, u_ptr, (size_t)n * sizeof(*next));
  }
  for (UF_long p = 0; p < seed->u.count; ++p) {
    const UF_long row = seed->u.row[p];
    const UF_long dst = next[row]++;
    u_cols[dst] = seed->u.col[p];
    u_values[dst] = seed->u.numeric_value[p];
    u_row_values[dst] = seed->u.value[p];
  }
  if (n > 0u) {
    memcpy(next, input_ptr, (size_t)n * sizeof(*next));
  }
  for (UF_long p = 0; p < seed->input.count; ++p) {
    const UF_long row = seed->input.row[p];
    const UF_long dst = next[row]++;
    input_cols[dst] = seed->input.col[p];
    input_pos[dst] = seed->input.pos[p];
  }

  kls_row_first_refactor_seed_clear(seed);
  free_row_refactor_pattern(solver);
  if (!kls_finish_row_refactor_pattern_from_arrays(
        solver, l_ptr, l_cols, l_values, l_row_values,
        u_ptr, u_cols, u_values, u_row_values,
        input_ptr, input_cols, input_pos, next)) {
    return 0;
  }
  solver->row_refactor_values_ready = 1;
  solver->row_refactor_values_dirty = 0;
  solver->row_refactor_solve_direct_ready = 0;
  solver->row_refactor_solve_validated = 0;
  solver->kls_first_last_row_refactor_seeded_rows = n;
  solver->kls_first_row_refactor_seeded_row_count += n;
  return 1;
}

typedef struct kls_row_first_parallel_shared {
  kls_solver *solver;
  trilinos_klu_l_numeric *numeric;
  double *numeric_values;
  const UF_long *psinv;
  UF_long *q_order;
  kls_row_first_refactor_seed *row_refactor_seed;
  const UF_long *separator_component_last;
  atomic_ulong next_block;
  atomic_int failed;
  pthread_mutex_t commit_lock;
  int scaled;
  int scale;
  int use_separator_pivot_domains;
  int force_block_pipeline;
  int block_pipeline_thread_count;
  UF_long block_pipeline_target_block;
  UF_long block_pipeline_begin;
  UF_long block_pipeline_end;
  const unsigned char *block_pipeline_active_mask;
  UF_long block_pipeline_active_rows;
  const UF_long *block_pipeline_row_order;
  UF_long block_pipeline_order_rows;
  int block_pipeline_etree_tail_order;
  int block_pipeline_separator_tail_scope;
  UF_long block_pipeline_separator_tail_scope_rows;
  double tol;
  UF_long total_lnz;
  UF_long total_unz;
  UF_long max_lnz_block;
  UF_long max_unz_block;
  UF_long dynamic_column_pivots;
  UF_long separator_dynamic_column_pivots;
  UF_long separator_extent_dynamic_column_pivots;
  UF_long separator_dynamic_column_fallbacks;
  UF_long row_supernode_update_groups;
  UF_long row_supernode_update_rows;
  UF_long row_supernode_panel_update_groups;
  UF_long row_supernode_panel_update_rows;
  UF_long row_panel_cache_builds;
  UF_long row_panel_cache_build_panels;
  UF_long row_panel_cache_build_entries;
  UF_long row_panel_cache_appends;
  UF_long row_panel_cache_append_panels;
  UF_long row_panel_cache_append_entries;
  UF_long active_rank_pivot_resets;
  UF_long active_rank_pivot_reset_rows;
  UF_long active_rank_pivot_panel_rebuilds;
  UF_long active_rank_pivot_panel_rebuild_rows;
  UF_long row_pipeline_runs;
  UF_long row_pipeline_rows;
  UF_long row_pipeline_threads;
  UF_long row_pipeline_prefix_rows;
  UF_long row_pipeline_suffix_rows;
  UF_long row_pipeline_gap_rows;
  UF_long row_pipeline_pivot_tail_rows;
  UF_long row_pipeline_pivot_restarts;
  UF_long row_pipeline_etree_ready_runs;
  UF_long row_pipeline_etree_ready_rows;
  UF_long row_pipeline_etree_ready_threads;
  UF_long row_pipeline_etree_prefactor_runs;
  UF_long row_pipeline_etree_prefactor_rows;
  UF_long row_pipeline_etree_prefactor_threads;
  UF_long row_pipeline_etree_prefactor_wait_rows;
  UF_long row_pipeline_etree_prefactor_wait_deps;
  UF_long row_pipeline_separator_queue_runs;
  UF_long row_pipeline_separator_private_rows;
  UF_long row_pipeline_separator_pipeline_rows;
  UF_long row_pipeline_separator_private_threads;
  int row_pipeline_separator_partitioned;
  UF_long row_pipeline_separator_split_components;
  atomic_ulong separator_dynamic_column_rejects;
  UF_long committed_blocks;
} kls_row_first_parallel_shared;

typedef struct kls_row_first_parallel_worker {
  kls_row_first_parallel_shared *shared;
  UF_long *row_counts;
  UF_long *row_ptr;
  double *x;
  double *supernode_workspace;
  unsigned int *mark;
  UF_long *pattern;
  UF_long *dep_heap;
  UF_long *u_row_ptr;
  UF_long *u_row_end;
  double *udiag_values;
  UF_long *col_pos;
  UF_long *block_q_order;
  UF_long supernode_workspace_capacity;
} kls_row_first_parallel_worker;

static void kls_row_first_parallel_worker_free(
  kls_row_first_parallel_worker *worker) {
  if (worker == NULL) {
    return;
  }
  free(worker->row_counts);
  free(worker->row_ptr);
  free(worker->x);
  free(worker->supernode_workspace);
  free(worker->mark);
  free(worker->pattern);
  free(worker->dep_heap);
  free(worker->u_row_ptr);
  free(worker->u_row_end);
  free(worker->udiag_values);
  free(worker->col_pos);
  free(worker->block_q_order);
  memset(worker, 0, sizeof(*worker));
}

static int kls_row_first_parallel_worker_init(
  kls_row_first_parallel_worker *worker,
  kls_row_first_parallel_shared *shared) {
  if (worker == NULL || shared == NULL || shared->solver == NULL ||
      shared->solver->symbolic == NULL) {
    return 0;
  }
  memset(worker, 0, sizeof(*worker));
  worker->shared = shared;
  const UF_long n = shared->solver->n;
  const UF_long maxblock = shared->solver->symbolic->maxblock;
  worker->row_counts =
    (UF_long *)calloc((size_t)maxblock, sizeof(*worker->row_counts));
  worker->row_ptr =
    (UF_long *)calloc((size_t)maxblock + 1u, sizeof(*worker->row_ptr));
  worker->x = (double *)calloc((size_t)maxblock, sizeof(*worker->x));
  worker->mark =
    (unsigned int *)calloc((size_t)maxblock, sizeof(*worker->mark));
  worker->pattern =
    (UF_long *)malloc((size_t)maxblock * sizeof(*worker->pattern));
  worker->dep_heap =
    (UF_long *)malloc((size_t)maxblock * sizeof(*worker->dep_heap));
  worker->u_row_ptr =
    (UF_long *)calloc((size_t)maxblock + 1u,
                      sizeof(*worker->u_row_ptr));
  worker->u_row_end =
    (UF_long *)calloc((size_t)maxblock, sizeof(*worker->u_row_end));
  worker->udiag_values =
    (double *)malloc((size_t)maxblock * sizeof(*worker->udiag_values));
  worker->col_pos = (UF_long *)malloc((size_t)n * sizeof(*worker->col_pos));
  worker->block_q_order =
    (UF_long *)malloc((size_t)maxblock * sizeof(*worker->block_q_order));
  if (worker->row_counts == NULL || worker->row_ptr == NULL ||
      worker->x == NULL || worker->mark == NULL ||
      worker->pattern == NULL || worker->dep_heap == NULL ||
      worker->u_row_ptr == NULL || worker->u_row_end == NULL ||
      worker->udiag_values == NULL || worker->col_pos == NULL ||
      worker->block_q_order == NULL) {
    kls_row_first_parallel_worker_free(worker);
    return 0;
  }
  for (UF_long k = 0; k < n; ++k) {
    worker->col_pos[k] = KLS_KLU_EMPTY;
  }
  return 1;
}

static int kls_row_first_parallel_commit_block(
  kls_row_first_parallel_worker *worker,
  UF_long block,
  UF_long k1,
  UF_long nk,
  const UF_long *row_ptr,
  const UF_long *row_cols,
  const UF_long *row_input_pos,
  const UF_long *local_row_order,
  const kls_row_first_entries *l_entries,
  const kls_row_first_entries *u_entries,
  UF_long dynamic_column_pivots,
  UF_long separator_dynamic_column_pivots,
  UF_long separator_extent_dynamic_column_pivots,
  UF_long separator_dynamic_column_fallbacks,
  UF_long row_supernode_update_groups,
  UF_long row_supernode_update_rows,
  UF_long row_supernode_panel_update_groups,
  UF_long row_supernode_panel_update_rows,
  UF_long row_panel_cache_builds,
  UF_long row_panel_cache_build_panels,
  UF_long row_panel_cache_build_entries,
  UF_long row_panel_cache_appends,
  UF_long row_panel_cache_append_panels,
  UF_long row_panel_cache_append_entries,
  UF_long active_rank_pivot_resets,
  UF_long active_rank_pivot_reset_rows,
  UF_long active_rank_pivot_panel_rebuilds,
  UF_long active_rank_pivot_panel_rebuild_rows,
  UF_long row_pipeline_rows,
  UF_long row_pipeline_threads,
  UF_long row_pipeline_prefix_rows,
  UF_long row_pipeline_suffix_rows,
  UF_long row_pipeline_gap_rows,
  UF_long row_pipeline_pivot_tail_rows,
  UF_long row_pipeline_pivot_restarts,
  UF_long row_pipeline_etree_ready_rows,
  UF_long row_pipeline_etree_ready_threads,
  UF_long row_pipeline_etree_prefactor_rows,
  UF_long row_pipeline_etree_prefactor_threads,
  UF_long row_pipeline_etree_prefactor_wait_rows,
  UF_long row_pipeline_etree_prefactor_wait_deps,
  int row_pipeline_separator_queue,
  UF_long row_pipeline_separator_private_rows,
  UF_long row_pipeline_separator_pipeline_rows,
  UF_long row_pipeline_separator_private_threads,
  int row_pipeline_separator_partitioned,
  UF_long row_pipeline_separator_split_components) {
  if (worker == NULL || worker->shared == NULL ||
      worker->shared->solver == NULL || row_ptr == NULL ||
      l_entries == NULL || u_entries == NULL ||
      worker->udiag_values == NULL || worker->block_q_order == NULL ||
      worker->col_pos == NULL) {
    return 0;
  }
  kls_row_first_parallel_shared *shared = worker->shared;
  int ok = 1;
  pthread_mutex_lock(&shared->commit_lock);
  if (atomic_load_explicit(&shared->failed, memory_order_acquire)) {
    ok = 0;
  } else if (!kls_pack_row_first_block_numeric(
               shared->solver, shared->numeric, block, k1, nk,
               local_row_order, l_entries, u_entries, worker->udiag_values,
               &shared->total_lnz, &shared->total_unz,
               &shared->max_lnz_block, &shared->max_unz_block)) {
    ok = 0;
  } else {
    memcpy(shared->q_order + k1, worker->block_q_order,
           (size_t)nk * sizeof(*worker->block_q_order));
    if (shared->row_refactor_seed != NULL &&
        shared->row_refactor_seed->enabled &&
        !kls_row_first_refactor_seed_capture_block(
          shared->row_refactor_seed, shared->solver, shared->numeric,
          block, k1, nk, row_ptr, row_cols, row_input_pos,
          worker->col_pos, l_entries, u_entries)) {
      kls_row_first_refactor_seed_clear(shared->row_refactor_seed);
    }
    shared->dynamic_column_pivots += dynamic_column_pivots;
    shared->separator_dynamic_column_pivots +=
      separator_dynamic_column_pivots;
    shared->separator_extent_dynamic_column_pivots +=
      separator_extent_dynamic_column_pivots;
    shared->separator_dynamic_column_fallbacks +=
      separator_dynamic_column_fallbacks;
    shared->row_supernode_update_groups +=
      row_supernode_update_groups;
    shared->row_supernode_update_rows += row_supernode_update_rows;
    shared->row_supernode_panel_update_groups +=
      row_supernode_panel_update_groups;
    shared->row_supernode_panel_update_rows +=
      row_supernode_panel_update_rows;
    kls_row_first_stats_add(&shared->row_panel_cache_builds,
                            row_panel_cache_builds);
    kls_row_first_stats_add(&shared->row_panel_cache_build_panels,
                            row_panel_cache_build_panels);
    kls_row_first_stats_add(&shared->row_panel_cache_build_entries,
                            row_panel_cache_build_entries);
    kls_row_first_stats_add(&shared->row_panel_cache_appends,
                            row_panel_cache_appends);
    kls_row_first_stats_add(&shared->row_panel_cache_append_panels,
                            row_panel_cache_append_panels);
    kls_row_first_stats_add(&shared->row_panel_cache_append_entries,
                            row_panel_cache_append_entries);
    kls_row_first_stats_add(&shared->active_rank_pivot_resets,
                            active_rank_pivot_resets);
    kls_row_first_stats_add(&shared->active_rank_pivot_reset_rows,
                            active_rank_pivot_reset_rows);
    kls_row_first_stats_add(&shared->active_rank_pivot_panel_rebuilds,
                            active_rank_pivot_panel_rebuilds);
    kls_row_first_stats_add(&shared->active_rank_pivot_panel_rebuild_rows,
                            active_rank_pivot_panel_rebuild_rows);
    if (row_pipeline_rows > 0u) {
      shared->row_pipeline_runs++;
      shared->row_pipeline_rows += row_pipeline_rows;
      if (row_pipeline_threads > shared->row_pipeline_threads) {
        shared->row_pipeline_threads = row_pipeline_threads;
      }
      shared->row_pipeline_prefix_rows += row_pipeline_prefix_rows;
      shared->row_pipeline_suffix_rows += row_pipeline_suffix_rows;
      shared->row_pipeline_gap_rows += row_pipeline_gap_rows;
      shared->row_pipeline_pivot_tail_rows += row_pipeline_pivot_tail_rows;
      shared->row_pipeline_pivot_restarts += row_pipeline_pivot_restarts;
      if (row_pipeline_etree_ready_rows > 0u) {
        shared->row_pipeline_etree_ready_runs++;
        shared->row_pipeline_etree_ready_rows +=
          row_pipeline_etree_ready_rows;
        if (row_pipeline_etree_ready_threads >
            shared->row_pipeline_etree_ready_threads) {
          shared->row_pipeline_etree_ready_threads =
            row_pipeline_etree_ready_threads;
        }
      }
      if (row_pipeline_etree_prefactor_rows > 0u) {
        shared->row_pipeline_etree_prefactor_runs++;
        shared->row_pipeline_etree_prefactor_rows +=
          row_pipeline_etree_prefactor_rows;
        if (row_pipeline_etree_prefactor_threads >
            shared->row_pipeline_etree_prefactor_threads) {
          shared->row_pipeline_etree_prefactor_threads =
            row_pipeline_etree_prefactor_threads;
        }
        shared->row_pipeline_etree_prefactor_wait_rows +=
          row_pipeline_etree_prefactor_wait_rows;
        shared->row_pipeline_etree_prefactor_wait_deps +=
          row_pipeline_etree_prefactor_wait_deps;
      }
    }
    if (row_pipeline_separator_queue) {
      shared->row_pipeline_separator_queue_runs++;
      shared->row_pipeline_separator_private_rows +=
        row_pipeline_separator_private_rows;
      shared->row_pipeline_separator_pipeline_rows +=
        row_pipeline_separator_pipeline_rows;
      if (row_pipeline_separator_private_threads >
          shared->row_pipeline_separator_private_threads) {
        shared->row_pipeline_separator_private_threads =
          row_pipeline_separator_private_threads;
      }
      if (row_pipeline_separator_partitioned) {
        shared->row_pipeline_separator_partitioned = 1;
      }
      shared->row_pipeline_separator_split_components +=
        row_pipeline_separator_split_components;
    }
    shared->committed_blocks++;
  }
  if (!ok) {
    atomic_store_explicit(&shared->failed, 1, memory_order_release);
  }
  pthread_mutex_unlock(&shared->commit_lock);
  return ok;
}

static int kls_row_first_parallel_seed_preserved_rows(
  kls_row_first_parallel_worker *worker,
  UF_long block,
  UF_long k1,
  UF_long nk,
  UF_long active_begin,
  UF_long active_end,
  const unsigned char *active_mask,
  kls_row_first_entries *l_entries,
  kls_row_first_entries *u_entries,
  kls_row_first_workspace *workspace,
  unsigned char *row_done,
  kls_row_first_supernode_panel_cache *panel_cache,
  UF_long *open_panel_start_io,
  UF_long *open_panel_end_io) {
  if (worker == NULL || worker->shared == NULL ||
      worker->shared->solver == NULL || worker->shared->numeric == NULL ||
      l_entries == NULL || u_entries == NULL || workspace == NULL ||
      row_done == NULL || worker->udiag_values == NULL ||
      worker->shared->psinv == NULL || active_begin > active_end ||
      active_end > nk ||
      block >= worker->shared->numeric->nblocks ||
      nk > (UF_long)(KLS_MAX_ALLOCATION / sizeof(UF_long))) {
    return 0;
  }
  if (active_begin == 0u && active_end == nk) {
    if (active_mask == NULL) {
      return 1;
    }
  }

  kls_row_first_parallel_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  trilinos_klu_l_numeric *numeric = shared->numeric;
  if (numeric->LUbx == NULL || numeric->LUbx[block] == NULL ||
      numeric->Lip == NULL || numeric->Llen == NULL ||
      numeric->Uip == NULL || numeric->Ulen == NULL ||
      numeric->Udiag == NULL || numeric->Pnum == NULL ||
      workspace->u_row_ptr == NULL || workspace->u_row_end == NULL ||
      k1 > solver->n || nk > solver->n - k1) {
    return 0;
  }

  for (UF_long k = 0; k < nk; ++k) {
    const int active = active_mask != NULL
      ? active_mask[k] != 0u
      : (k >= active_begin && k < active_end);
    if (active) {
      continue;
    }
    const UF_long old_row = numeric->Pnum[k1 + k];
    if (old_row >= solver->n || shared->psinv[old_row] != k1 + k) {
      return 0;
    }
  }

  Unit *lu = (Unit *)numeric->LUbx[block];
  const UF_long *lip = numeric->Lip + k1;
  const UF_long *llen = numeric->Llen + k1;
  const UF_long *uip = numeric->Uip + k1;
  const UF_long *ulen = numeric->Ulen + k1;
  const Entry *udiag = ((Entry *)numeric->Udiag) + k1;

  for (UF_long col = 0; col < nk; ++col) {
    const int col_active = active_mask != NULL
      ? active_mask[col] != 0u
      : (col >= active_begin && col < active_end);
    UF_long *li = NULL;
    Entry *lx = NULL;
    UF_long l_len = 0;
    kls_klu_get_pointer((double *)lu, lip, llen, col, &li, &lx, &l_len);
    for (UF_long p = 0; p < l_len; ++p) {
      const UF_long row = li[p];
      if (row >= nk) {
        return 0;
      }
      const int row_active = active_mask != NULL
        ? active_mask[row] != 0u
        : (row >= active_begin && row < active_end);
      if (row_active) {
        continue;
      }
      if (col_active) {
        return 0;
      }
      if (row <= col ||
          !kls_row_first_entries_append(l_entries, row, col, lx[p])) {
        return 0;
      }
    }
  }

  /* Transpose preserved strict-U entries in two linear passes.  The older
     row-at-a-time seed rescanned every later U column for every preserved
     row, making a small late ETree tail more expensive than a whole-block
     factor.  Count row buckets once, reserve once, then fill the same
     contiguous row representation while maintaining the column links used
     by dynamic pivot updates. */
  const size_t nk_size = (size_t)nk;
  UF_long *preserved_u_count =
    (UF_long *)malloc(nk_size * sizeof(*preserved_u_count));
  UF_long *preserved_u_offset =
    (UF_long *)malloc(nk_size * sizeof(*preserved_u_offset));
  if ((nk > 0u && preserved_u_count == NULL) ||
      (nk > 0u && preserved_u_offset == NULL)) {
    free(preserved_u_count);
    free(preserved_u_offset);
    return 0;
  }
  if (nk_size > 0u) {
    memset(preserved_u_count, 0,
           nk_size * sizeof(*preserved_u_count));
  }
  int preserved_u_ok = 1;
  for (UF_long col = 0; col < nk && preserved_u_ok; ++col) {
    UF_long *ui = NULL;
    Entry *ux = NULL;
    UF_long u_len = 0;
    kls_klu_get_pointer((double *)lu, uip, ulen, col, &ui, &ux, &u_len);
    for (UF_long p = 0; p < u_len; ++p) {
      const UF_long row = ui[p];
      if (row >= nk || row >= col) {
        preserved_u_ok = 0;
        break;
      }
      const int row_active = active_mask != NULL
        ? active_mask[row] != 0u
        : (row >= active_begin && row < active_end);
      if (!row_active) {
        if (preserved_u_count[row] == UF_long_max) {
          preserved_u_ok = 0;
          break;
        }
        preserved_u_count[row]++;
      }
    }
  }
  UF_long preserved_u_total = 0u;
  for (UF_long row = 0; row < nk && preserved_u_ok; ++row) {
    preserved_u_offset[row] = preserved_u_total;
    if (preserved_u_count[row] > UF_long_max - preserved_u_total) {
      preserved_u_ok = 0;
      break;
    }
    preserved_u_total += preserved_u_count[row];
  }
  const UF_long preserved_u_base = u_entries->count;
  if (!preserved_u_ok ||
      !kls_row_first_entries_reserve_append(u_entries,
                                            preserved_u_total)) {
    free(preserved_u_count);
    free(preserved_u_offset);
    return 0;
  }
  for (UF_long col = 0; col < nk; ++col) {
    UF_long *ui = NULL;
    Entry *ux = NULL;
    UF_long u_len = 0;
    kls_klu_get_pointer((double *)lu, uip, ulen, col, &ui, &ux, &u_len);
    for (UF_long p = 0; p < u_len; ++p) {
      const UF_long row = ui[p];
      const int row_active = active_mask != NULL
        ? active_mask[row] != 0u
        : (row >= active_begin && row < active_end);
      if (row_active) {
        continue;
      }
      const UF_long pos = preserved_u_base + preserved_u_offset[row]++;
      u_entries->row[pos] = row;
      u_entries->col[pos] = col;
      u_entries->value[pos] = ux[p];
      if (u_entries->col_count != NULL) {
        u_entries->col_count[col]++;
      }
      if (u_entries->col_head != NULL) {
        u_entries->next[pos] = u_entries->col_head[col];
        u_entries->col_head[col] = pos;
      }
    }
  }
  u_entries->count = preserved_u_base + preserved_u_total;
  for (UF_long row = 0; row < nk; ++row) {
    const int row_active = active_mask != NULL
      ? active_mask[row] != 0u
      : (row >= active_begin && row < active_end);
    if (row_active) {
      continue;
    }
    const UF_long row_end = preserved_u_base + preserved_u_offset[row];
    workspace->u_row_ptr[row] = row_end - preserved_u_count[row];
    workspace->u_row_end[row] = row_end;
    worker->udiag_values[row] = udiag[row];
    row_done[row] = 1u;
    kls_row_first_supernode_panel_cache_publish_completed(
      workspace->supernode_panel_cache != NULL ? panel_cache : NULL,
      u_entries, workspace->u_row_ptr, workspace->u_row_end,
      worker->udiag_values, nk, NULL, open_panel_start_io,
      open_panel_end_io, row);
  }
  free(preserved_u_count);
  free(preserved_u_offset);
  return 1;
}

static int kls_try_row_first_etree_ready_tail_repair(
  kls_row_first_parallel_worker *worker,
  UF_long block,
  UF_long k1,
  UF_long nk,
  const kls_row_first_block_context *ctx,
  const UF_long *row_order,
  const unsigned char *active_mask,
  UF_long active_rows,
  UF_long active_begin,
  UF_long active_end,
  UF_long reserve_entries,
  kls_row_first_entries *l_entries,
  kls_row_first_entries *u_entries,
  kls_row_first_row_stats *row_stats,
  UF_long *row_pipeline_ready_rows_out,
  UF_long *row_pipeline_threads_out,
  UF_long *row_pipeline_prefactor_rows_out,
  UF_long *row_pipeline_prefactor_wait_rows_out,
  UF_long *row_pipeline_prefactor_wait_deps_out,
  UF_long *supernode_update_groups_out,
  UF_long *supernode_update_rows_out,
  UF_long *supernode_panel_update_groups_out,
  UF_long *supernode_panel_update_rows_out,
  UF_long *pivot_tail_rows_out,
  UF_long *pivot_restarts_out) {
  KLS_SET_OPTIONAL_OUTPUT(row_pipeline_ready_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(row_pipeline_threads_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(row_pipeline_prefactor_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(row_pipeline_prefactor_wait_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(row_pipeline_prefactor_wait_deps_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(supernode_update_groups_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(supernode_update_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(supernode_panel_update_groups_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(supernode_panel_update_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(pivot_tail_rows_out, 0);
  KLS_SET_OPTIONAL_OUTPUT(pivot_restarts_out, 0);
  if (worker == NULL || worker->shared == NULL ||
      worker->udiag_values == NULL || worker->u_row_ptr == NULL ||
      worker->u_row_end == NULL || ctx == NULL || ctx->solver == NULL ||
      row_order == NULL || active_mask == NULL || active_rows == 0u ||
      active_begin >= active_end || active_end > nk ||
      l_entries == NULL || u_entries == NULL || row_stats == NULL ||
      ctx->solver->fast_reject_tail_parent == NULL ||
      ctx->solver->fast_reject_tail_plan_block != block ||
      ctx->solver->fast_reject_tail_plan_k1 != k1 ||
      ctx->solver->fast_reject_tail_plan_nk != nk ||
      worker->shared->block_pipeline_thread_count <= 1) {
    return 0;
  }

  kls_row_first_entries trial_l_entries;
  kls_row_first_entries trial_u_entries;
  memset(&trial_l_entries, 0, sizeof(trial_l_entries));
  memset(&trial_u_entries, 0, sizeof(trial_u_entries));
  kls_row_first_workspace trial_workspace;
  memset(&trial_workspace, 0, sizeof(trial_workspace));
  kls_row_first_supernode_panel_cache trial_panel_cache;
  memset(&trial_panel_cache, 0, sizeof(trial_panel_cache));
  double *trial_udiag =
    (double *)malloc((size_t)nk * sizeof(*trial_udiag));
  unsigned char *trial_row_done =
    (unsigned char *)calloc((size_t)nk, sizeof(*trial_row_done));
  trial_workspace.x = (double *)calloc((size_t)nk, sizeof(*trial_workspace.x));
  trial_workspace.mark =
    (unsigned int *)calloc((size_t)nk, sizeof(*trial_workspace.mark));
  trial_workspace.pattern =
    (UF_long *)malloc((size_t)nk * sizeof(*trial_workspace.pattern));
  trial_workspace.dep_heap =
    (UF_long *)malloc((size_t)nk * sizeof(*trial_workspace.dep_heap));
  trial_workspace.u_row_ptr =
    (UF_long *)calloc((size_t)nk + 1u, sizeof(*trial_workspace.u_row_ptr));
  trial_workspace.u_row_end =
    (UF_long *)calloc((size_t)nk, sizeof(*trial_workspace.u_row_end));
  int ok = 0;
  if (trial_udiag == NULL || trial_row_done == NULL ||
      trial_workspace.x == NULL || trial_workspace.mark == NULL ||
      trial_workspace.pattern == NULL || trial_workspace.dep_heap == NULL ||
      trial_workspace.u_row_ptr == NULL ||
      trial_workspace.u_row_end == NULL ||
      !kls_row_first_entries_enable_col_counts(&trial_l_entries, nk) ||
      !kls_row_first_entries_enable_col_counts(&trial_u_entries, nk) ||
      !kls_row_first_entries_enable_col_links(&trial_u_entries, nk)) {
    goto cleanup;
  }
  if (reserve_entries > 0u) {
    (void)kls_row_first_entries_reserve(&trial_l_entries, reserve_entries);
    (void)kls_row_first_entries_reserve(&trial_u_entries, reserve_entries);
  }
  if (kls_row_first_supernode_panel_cache_init_rows(
        &trial_panel_cache, nk)) {
    trial_workspace.supernode_panel_cache = &trial_panel_cache;
  }

  double *saved_udiag = worker->udiag_values;
  worker->udiag_values = trial_udiag;
  const int seed_ok = kls_row_first_parallel_seed_preserved_rows(
    worker, block, k1, nk, active_begin, active_end, active_mask,
    &trial_l_entries, &trial_u_entries, &trial_workspace, trial_row_done,
    &trial_panel_cache, NULL, NULL);
  worker->udiag_values = saved_udiag;
  if (!seed_ok) {
    goto cleanup;
  }

  unsigned char *ready_mask =
    (unsigned char *)malloc((size_t)nk * sizeof(*ready_mask));
  if (ready_mask == NULL) {
    goto cleanup;
  }
  memcpy(ready_mask, active_mask, (size_t)nk * sizeof(*ready_mask));
  const UF_long boundary_row = row_order[0];
  UF_long supernode_valid_begin = 0u;
  const UF_long boundary_pivots_before =
    row_stats != NULL ? row_stats->dynamic_column_pivots : 0u;
  if (boundary_row >= nk || ready_mask[boundary_row] == 0u ||
      !kls_row_first_factor_one_row(
        ctx, &trial_workspace, &trial_l_entries, &trial_u_entries,
        trial_udiag, boundary_row, NULL, -1, trial_row_done, row_stats,
        NULL)) {
    free(ready_mask);
    goto cleanup;
  }
  if (row_stats != NULL &&
      row_stats->dynamic_column_pivots != boundary_pivots_before) {
    supernode_valid_begin =
      boundary_row >= nk ? nk : boundary_row + 1u;
  }
  ready_mask[boundary_row] = 0u;
  if (trial_workspace.supernode_panel_cache != NULL) {
    kls_row_first_supernode_panel_cache_free(&trial_panel_cache);
    if (kls_row_first_supernode_panel_cache_init_rows(
          &trial_panel_cache, nk)) {
      trial_workspace.supernode_panel_cache = &trial_panel_cache;
    } else {
      trial_workspace.supernode_panel_cache = NULL;
    }
  }

  UF_long ready_threads = 0;
  UF_long ready_rows = active_rows - 1u;
  UF_long prefactor_rows = 0;
  UF_long prefactor_wait_rows = 0;
  UF_long prefactor_wait_deps = 0;
  UF_long supernode_update_groups = 0;
  UF_long supernode_update_rows = 0;
  UF_long supernode_panel_update_groups = 0;
  UF_long supernode_panel_update_rows = 0;
  UF_long pivot_tail_rows = 0;
  UF_long pivot_restarts = 0;
  if (ready_rows > 0u &&
      !kls_row_first_run_etree_prefactor_tail_phase(
        ctx, row_order + 1u, ready_mask,
        ctx->solver->fast_reject_tail_parent, ready_rows,
        supernode_valid_begin,
        worker->shared->block_pipeline_thread_count, &trial_workspace,
        &trial_l_entries, &trial_u_entries, trial_udiag, trial_row_done,
        row_stats, &ready_threads, &prefactor_rows, &prefactor_wait_rows,
        &prefactor_wait_deps, &supernode_update_groups,
        &supernode_update_rows, &supernode_panel_update_groups,
        &supernode_panel_update_rows, &pivot_tail_rows,
        &pivot_restarts)) {
    free(ready_mask);
    goto cleanup;
  }
  free(ready_mask);

  memcpy(worker->u_row_ptr, trial_workspace.u_row_ptr,
         ((size_t)nk + 1u) * sizeof(*worker->u_row_ptr));
  memcpy(worker->u_row_end, trial_workspace.u_row_end,
         (size_t)nk * sizeof(*worker->u_row_end));
  memcpy(worker->udiag_values, trial_udiag,
         (size_t)nk * sizeof(*worker->udiag_values));
  kls_row_first_entries_free(l_entries);
  kls_row_first_entries_free(u_entries);
  *l_entries = trial_l_entries;
  *u_entries = trial_u_entries;
  memset(&trial_l_entries, 0, sizeof(trial_l_entries));
  memset(&trial_u_entries, 0, sizeof(trial_u_entries));
  if (row_pipeline_ready_rows_out != NULL) {
    *row_pipeline_ready_rows_out = ready_rows;
  }
  if (row_pipeline_threads_out != NULL) {
    *row_pipeline_threads_out = ready_threads;
  }
  if (row_pipeline_prefactor_rows_out != NULL) {
    *row_pipeline_prefactor_rows_out = prefactor_rows;
  }
  if (row_pipeline_prefactor_wait_rows_out != NULL) {
    *row_pipeline_prefactor_wait_rows_out = prefactor_wait_rows;
  }
  if (row_pipeline_prefactor_wait_deps_out != NULL) {
    *row_pipeline_prefactor_wait_deps_out = prefactor_wait_deps;
  }
  if (supernode_update_groups_out != NULL) {
    *supernode_update_groups_out = supernode_update_groups;
  }
  if (supernode_update_rows_out != NULL) {
    *supernode_update_rows_out = supernode_update_rows;
  }
  if (supernode_panel_update_groups_out != NULL) {
    *supernode_panel_update_groups_out = supernode_panel_update_groups;
  }
  if (supernode_panel_update_rows_out != NULL) {
    *supernode_panel_update_rows_out = supernode_panel_update_rows;
  }
  if (pivot_tail_rows_out != NULL) {
    *pivot_tail_rows_out = pivot_tail_rows;
  }
  if (pivot_restarts_out != NULL) {
    *pivot_restarts_out = pivot_restarts;
  }
  ok = 1;

cleanup:
  kls_row_first_entries_free(&trial_l_entries);
  kls_row_first_entries_free(&trial_u_entries);
  kls_row_first_supernode_panel_cache_free(&trial_panel_cache);
  /* the manual list missed supernode_fit_status (leaked) and left the
     struct un-zeroed; use the canonical free */
  kls_row_first_workspace_free(&trial_workspace);
  free(trial_udiag);
  free(trial_row_done);
  return ok;
}

static int kls_row_first_parallel_factor_block(
  kls_row_first_parallel_worker *worker,
  UF_long block) {
  if (worker == NULL || worker->shared == NULL ||
      worker->shared->solver == NULL ||
      worker->shared->solver->symbolic == NULL ||
      worker->shared->numeric == NULL) {
    return 0;
  }
  kls_row_first_parallel_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  const UF_long n = solver->n;
  const UF_long k1 = solver->symbolic->R[block];
  const UF_long k2 = solver->symbolic->R[block + 1u];
  if (k2 <= k1 || k2 > n ||
      k2 - k1 > solver->symbolic->maxblock) {
    return 0;
  }
  const UF_long nk = k2 - k1;
  const int use_separator_for_block =
    shared->use_separator_pivot_domains &&
    shared->separator_component_last != NULL &&
    k1 >= solver->separator.global_begin &&
    k2 <= solver->separator.global_end;
  const UF_long separator_block_base = use_separator_for_block
    ? k1 - solver->separator.global_begin
    : 0u;

  memset(worker->row_counts, 0,
         (size_t)nk * sizeof(*worker->row_counts));
  memset(worker->u_row_ptr, 0,
         ((size_t)nk + 1u) * sizeof(*worker->u_row_ptr));
  memset(worker->u_row_end, 0, (size_t)nk * sizeof(*worker->u_row_end));

  for (UF_long local_col = 0; local_col < nk; ++local_col) {
    const UF_long oldcol = solver->symbolic->Q[k1 + local_col];
    if (oldcol >= n || worker->col_pos[oldcol] != KLS_KLU_EMPTY) {
      return 0;
    }
    worker->block_q_order[local_col] = oldcol;
    worker->col_pos[oldcol] = local_col;
  }

  int ok = 0;
  UF_long *row_cols = NULL;
  UF_long *row_input_pos = NULL;
  UF_long *row_next = NULL;
  UF_long *row_order = NULL;
  UF_long *local_row_order = NULL;
  UF_long *old_to_new_row = NULL;
  unsigned char *row_seen = NULL;
  unsigned char *private_active_mask = NULL;
  UF_long *active_rank = NULL;
  int *row_owner = NULL;
  unsigned char *row_done = NULL;
  kls_row_first_entries l_entries;
  kls_row_first_entries u_entries;
  memset(&l_entries, 0, sizeof(l_entries));
  memset(&u_entries, 0, sizeof(u_entries));
  kls_first_separator_queue_plan separator_plan;
  memset(&separator_plan, 0, sizeof(separator_plan));
  int use_separator_row_order = 0;
  int separator_private_executed = 0;
  int separator_private_owner_valid = 0;
  int separator_private_run_ok = 0;
  int separator_private_complete = 0;
  int separator_private_partial_active = 0;
  int row_pipeline_separator_queue = 0;
  int row_pipeline_separator_partitioned = 0;
  UF_long row_pipeline_separator_private_rows = 0;
  UF_long row_pipeline_separator_pipeline_rows = 0;
  UF_long row_pipeline_separator_private_threads = 0;
  UF_long row_pipeline_separator_split_components = 0;
  UF_long block_separator_private_rows = 0;
  UF_long block_separator_pipeline_rows = 0;
  UF_long dynamic_column_pivots = 0;
  UF_long separator_dynamic_column_pivots = 0;
  UF_long separator_extent_dynamic_column_pivots = 0;
  UF_long separator_dynamic_column_fallbacks = 0;
  UF_long row_supernode_update_groups = 0;
  UF_long row_supernode_update_rows = 0;
  UF_long row_supernode_panel_update_groups = 0;
  UF_long row_supernode_panel_update_rows = 0;
  UF_long row_panel_cache_builds = 0;
  UF_long row_panel_cache_build_panels = 0;
  UF_long row_panel_cache_build_entries = 0;
  UF_long row_panel_cache_appends = 0;
  UF_long row_panel_cache_append_panels = 0;
  UF_long row_panel_cache_append_entries = 0;
  UF_long active_rank_pivot_resets = 0;
  UF_long active_rank_pivot_reset_rows = 0;
  UF_long active_rank_pivot_panel_rebuilds = 0;
  UF_long active_rank_pivot_panel_rebuild_rows = 0;
  UF_long row_pipeline_rows = 0;
  UF_long row_pipeline_threads = 0;
  UF_long row_pipeline_prefix_rows = 0;
  UF_long row_pipeline_suffix_rows = 0;
  UF_long row_pipeline_gap_rows = 0;
  UF_long row_pipeline_pivot_tail_rows = 0;
  UF_long row_pipeline_pivot_restarts = 0;
  UF_long row_pipeline_etree_ready_rows = 0;
  UF_long row_pipeline_etree_ready_threads = 0;
  UF_long row_pipeline_etree_prefactor_rows = 0;
  UF_long row_pipeline_etree_prefactor_threads = 0;
  UF_long row_pipeline_etree_prefactor_wait_rows = 0;
  UF_long row_pipeline_etree_prefactor_wait_deps = 0;
  UF_long open_panel_start = KLS_KLU_EMPTY;
  UF_long open_panel_end = KLS_KLU_EMPTY;
  UF_long trace_active_rows = 0;
  int trace_etree_attempted = 0;
  int trace_etree_succeeded = 0;
  int trace_seed_succeeded = 0;
  int trace_pipeline_succeeded = 0;
  int trace_commit_attempted = 0;
  kls_row_first_supernode_panel_cache row_supernode_panel_cache;
  memset(&row_supernode_panel_cache, 0, sizeof(row_supernode_panel_cache));
  kls_row_first_workspace row_workspace;
  memset(&row_workspace, 0, sizeof(row_workspace));
  row_workspace.x = worker->x;
  row_workspace.supernode_workspace = worker->supernode_workspace;
  row_workspace.mark = worker->mark;
  row_workspace.pattern = worker->pattern;
  row_workspace.dep_heap = worker->dep_heap;
  row_workspace.u_row_ptr = worker->u_row_ptr;
  row_workspace.u_row_end = worker->u_row_end;
  row_workspace.supernode_workspace_capacity =
    worker->supernode_workspace_capacity;
  if (kls_row_first_supernode_panel_cache_init_rows(
        &row_supernode_panel_cache, nk)) {
    row_workspace.supernode_panel_cache = &row_supernode_panel_cache;
  }

  for (UF_long local_col = 0; local_col < nk; ++local_col) {
    const UF_long oldcol = worker->block_q_order[local_col];
    for (UF_long p = solver->col_ptr[oldcol];
         p < solver->col_ptr[oldcol + 1u]; ++p) {
      const UF_long oldrow = solver->row_idx[p];
      if (oldrow >= n || shared->psinv[oldrow] >= n) {
        goto cleanup;
      }
      const UF_long row_pos = shared->psinv[oldrow];
      if (row_pos < k1) {
        continue;
      }
      if (row_pos >= k2) {
        goto cleanup;
      }
      worker->row_counts[row_pos - k1]++;
    }
  }
  worker->row_ptr[0] = 0;
  for (UF_long i = 0; i < nk; ++i) {
    worker->row_ptr[i + 1u] =
      worker->row_ptr[i] + worker->row_counts[i];
  }
  const UF_long row_nnz = worker->row_ptr[nk];
  row_cols = row_nnz > 0u
    ? (UF_long *)malloc((size_t)row_nnz * sizeof(*row_cols)) : NULL;
  row_input_pos = row_nnz > 0u
    ? (UF_long *)malloc((size_t)row_nnz * sizeof(*row_input_pos)) : NULL;
  row_next = (UF_long *)malloc((size_t)nk * sizeof(*row_next));
  if ((row_nnz > 0u && (row_cols == NULL || row_input_pos == NULL)) ||
      row_next == NULL) {
    goto cleanup;
  }
  memcpy(row_next, worker->row_ptr, (size_t)nk * sizeof(*row_next));
  for (UF_long local_col = 0; local_col < nk; ++local_col) {
    const UF_long oldcol = worker->block_q_order[local_col];
    for (UF_long p = solver->col_ptr[oldcol];
         p < solver->col_ptr[oldcol + 1u]; ++p) {
      const UF_long row_pos = shared->psinv[solver->row_idx[p]];
      if (row_pos < k1) {
        continue;
      }
      const UF_long row = row_pos - k1;
      const UF_long dst = row_next[row]++;
      row_cols[dst] = oldcol;
      row_input_pos[dst] = p;
    }
  }
  free(row_next);
  row_next = NULL;

  if (!kls_row_first_entries_enable_col_counts(&l_entries, nk) ||
      !kls_row_first_entries_enable_col_counts(&u_entries, nk) ||
      !kls_row_first_entries_enable_col_links(&u_entries, nk)) {
    goto cleanup;
  }
  const UF_long reserve_entries =
    kls_row_first_symbolic_reserve(solver, block, nk);
  if (reserve_entries > 0u) {
    (void)kls_row_first_entries_reserve(&l_entries, reserve_entries);
    (void)kls_row_first_entries_reserve(&u_entries, reserve_entries);
  }

  const int block_pipeline_allowed =
    shared->force_block_pipeline &&
    shared->block_pipeline_thread_count > 1 &&
    nk > 1u &&
    (shared->block_pipeline_target_block == KLS_KLU_EMPTY ||
     shared->block_pipeline_target_block == block);

  if (block_pipeline_allowed &&
      shared->block_pipeline_active_mask == NULL &&
      use_separator_for_block) {
    int have_separator_plan = kls_build_first_separator_queue_plan(
      solver, k1, k2, shared->block_pipeline_thread_count,
      &separator_plan);
    const int initial_separator_plan = have_separator_plan;
    const int initial_separator_partitioned =
      have_separator_plan && separator_plan.partitioned;
    int initial_separator_ownership_valid =
      initial_separator_partitioned ? 0 : 1;
    if (have_separator_plan && separator_plan.partitioned) {
      initial_separator_ownership_valid =
        kls_first_separator_private_ownership_valid(
          k1, nk, n, worker->row_ptr, row_cols, worker->col_pos,
          &separator_plan, shared->block_pipeline_thread_count);
    }
    if (have_separator_plan && separator_plan.partitioned &&
        !initial_separator_ownership_valid) {
      kls_first_separator_queue_plan_free(&separator_plan);
      have_separator_plan =
        kls_build_first_separator_queue_plan_legacy_from_counts(
          solver, k1, k2, shared->block_pipeline_thread_count,
          worker->row_counts, &separator_plan);
    }
    if (kls_first_factor_trace_env_enabled()) {
      fprintf(stderr,
              "KLS first-factor separator-queue decision:"
              " block=%" PRIu64 " rows=%" PRIu64
              " initial_plan=%d initial_partitioned=%d"
              " initial_ownership_valid=%d final_plan=%d"
              " final_partitioned=%d final_private_rows=%" PRIu64
              " final_pipeline_rows=%" PRIu64
              " final_nonempty_threads=%" PRIu64 "\n",
              (uint64_t)block, (uint64_t)nk,
              initial_separator_plan, initial_separator_partitioned,
              initial_separator_ownership_valid, have_separator_plan,
              have_separator_plan ? separator_plan.partitioned : 0,
              (uint64_t)(have_separator_plan
                           ? separator_plan.private_rows_count : 0u),
              (uint64_t)(have_separator_plan
                           ? separator_plan.pipeline_rows_count : 0u),
              (uint64_t)(have_separator_plan
                           ? separator_plan.nonempty_threads : 0u));
    }
    if (have_separator_plan) {
      row_order = (UF_long *)malloc((size_t)nk * sizeof(*row_order));
      row_seen = (unsigned char *)calloc((size_t)nk, sizeof(*row_seen));
      if (row_order == NULL || row_seen == NULL) {
        goto cleanup;
      }
      UF_long row_order_count = 0;
      for (UF_long pos = 0; pos < separator_plan.private_rows_count; ++pos) {
        const UF_long global_row = separator_plan.private_rows[pos];
        if (global_row < k1 || global_row >= k2) {
          goto cleanup;
        }
        const UF_long local_row = global_row - k1;
        if (row_seen[local_row]) {
          goto cleanup;
        }
        row_seen[local_row] = 1u;
        row_order[row_order_count++] = local_row;
      }
      block_separator_private_rows = row_order_count;
      qsort(row_order, (size_t)block_separator_private_rows,
            sizeof(*row_order), kls_compare_uf_long);
      for (UF_long pos = 0; pos < separator_plan.pipeline_rows_count; ++pos) {
        const UF_long global_row = separator_plan.pipeline_rows[pos];
        if (global_row < k1 || global_row >= k2) {
          goto cleanup;
        }
        const UF_long local_row = global_row - k1;
        if (row_seen[local_row]) {
          goto cleanup;
        }
        row_seen[local_row] = 1u;
        row_order[row_order_count++] = local_row;
      }
      block_separator_pipeline_rows =
        row_order_count - block_separator_private_rows;
      qsort(row_order + block_separator_private_rows,
            (size_t)block_separator_pipeline_rows, sizeof(*row_order),
            kls_compare_uf_long);
      if (row_order_count == nk &&
          block_separator_private_rows > 0u &&
          block_separator_pipeline_rows > 0u) {
        use_separator_row_order = 1;
      }
    }
  }

  if (use_separator_row_order) {
    local_row_order = row_order;
    row_order = (UF_long *)malloc((size_t)nk * sizeof(*row_order));
    old_to_new_row =
      (UF_long *)malloc((size_t)nk * sizeof(*old_to_new_row));
    if (row_order == NULL || old_to_new_row == NULL) {
      goto cleanup;
    }
    for (UF_long new_row = 0; new_row < nk; ++new_row) {
      const UF_long old_row = local_row_order[new_row];
      if (old_row >= nk) {
        goto cleanup;
      }
      row_order[new_row] = new_row;
      old_to_new_row[old_row] = new_row;
    }
    for (UF_long pos = 0; pos < separator_plan.private_rows_count; ++pos) {
      const UF_long old_global = separator_plan.private_rows[pos];
      if (old_global < k1 || old_global >= k2) {
        goto cleanup;
      }
      separator_plan.private_rows[pos] =
        k1 + old_to_new_row[old_global - k1];
    }
    for (UF_long pos = 0; pos < separator_plan.pipeline_rows_count; ++pos) {
      const UF_long old_global = separator_plan.pipeline_rows[pos];
      if (old_global < k1 || old_global >= k2) {
        goto cleanup;
      }
      separator_plan.pipeline_rows[pos] =
        k1 + old_to_new_row[old_global - k1];
    }

    UF_long *old_block_q =
      (UF_long *)malloc((size_t)nk * sizeof(*old_block_q));
    if (old_block_q == NULL) {
      goto cleanup;
    }
    for (UF_long local_col = 0; local_col < nk; ++local_col) {
      old_block_q[local_col] = worker->block_q_order[local_col];
      if (old_block_q[local_col] < n) {
        worker->col_pos[old_block_q[local_col]] = KLS_KLU_EMPTY;
      }
    }
    for (UF_long local_col = 0; local_col < nk; ++local_col) {
      const UF_long old_local_col = local_row_order[local_col];
      if (old_local_col >= nk) {
        free(old_block_q);
        goto cleanup;
      }
      worker->block_q_order[local_col] = old_block_q[old_local_col];
    }
    free(old_block_q);

    memset(worker->row_counts, 0,
           (size_t)nk * sizeof(*worker->row_counts));
    memset(worker->row_ptr, 0,
           ((size_t)nk + 1u) * sizeof(*worker->row_ptr));
    for (UF_long local_col = 0; local_col < nk; ++local_col) {
      const UF_long oldcol = worker->block_q_order[local_col];
      if (oldcol >= n || worker->col_pos[oldcol] != KLS_KLU_EMPTY) {
        goto cleanup;
      }
      worker->col_pos[oldcol] = local_col;
    }

    for (UF_long local_col = 0; local_col < nk; ++local_col) {
      const UF_long oldcol = worker->block_q_order[local_col];
      for (UF_long p = solver->col_ptr[oldcol];
           p < solver->col_ptr[oldcol + 1u]; ++p) {
        const UF_long oldrow = solver->row_idx[p];
        if (oldrow >= n || shared->psinv[oldrow] >= n) {
          goto cleanup;
        }
        const UF_long row_pos = shared->psinv[oldrow];
        if (row_pos < k1) {
          continue;
        }
        if (row_pos >= k2) {
          goto cleanup;
        }
        const UF_long old_local_row = row_pos - k1;
        if (old_local_row >= nk) {
          goto cleanup;
        }
        worker->row_counts[old_to_new_row[old_local_row]]++;
      }
    }
    worker->row_ptr[0] = 0;
    for (UF_long i = 0; i < nk; ++i) {
      worker->row_ptr[i + 1u] =
        worker->row_ptr[i] + worker->row_counts[i];
    }
    const UF_long rebuilt_row_nnz = worker->row_ptr[nk];
    UF_long *rebuilt_row_cols = rebuilt_row_nnz > 0u
      ? (UF_long *)malloc((size_t)rebuilt_row_nnz *
                          sizeof(*rebuilt_row_cols)) : NULL;
    UF_long *rebuilt_row_input_pos = rebuilt_row_nnz > 0u
      ? (UF_long *)malloc((size_t)rebuilt_row_nnz *
                          sizeof(*rebuilt_row_input_pos)) : NULL;
    row_next = (UF_long *)malloc((size_t)nk * sizeof(*row_next));
    if ((rebuilt_row_nnz > 0u &&
         (rebuilt_row_cols == NULL || rebuilt_row_input_pos == NULL)) ||
        row_next == NULL) {
      free(rebuilt_row_cols);
      free(rebuilt_row_input_pos);
      goto cleanup;
    }
    memcpy(row_next, worker->row_ptr, (size_t)nk * sizeof(*row_next));
    for (UF_long local_col = 0; local_col < nk; ++local_col) {
      const UF_long oldcol = worker->block_q_order[local_col];
      for (UF_long p = solver->col_ptr[oldcol];
           p < solver->col_ptr[oldcol + 1u]; ++p) {
        const UF_long oldrow = solver->row_idx[p];
        if (oldrow >= n || shared->psinv[oldrow] >= n) {
          free(rebuilt_row_cols);
          free(rebuilt_row_input_pos);
          goto cleanup;
        }
        const UF_long row_pos = shared->psinv[oldrow];
        if (row_pos < k1) {
          continue;
        }
        if (row_pos >= k2) {
          free(rebuilt_row_cols);
          free(rebuilt_row_input_pos);
          goto cleanup;
        }
        const UF_long old_local_row = row_pos - k1;
        if (old_local_row >= nk) {
          free(rebuilt_row_cols);
          free(rebuilt_row_input_pos);
          goto cleanup;
        }
        const UF_long row = old_to_new_row[old_local_row];
        const UF_long dst = row_next[row]++;
        rebuilt_row_cols[dst] = oldcol;
        rebuilt_row_input_pos[dst] = p;
      }
    }
    free(row_next);
    row_next = NULL;
    free(row_cols);
    free(row_input_pos);
    row_cols = rebuilt_row_cols;
    row_input_pos = rebuilt_row_input_pos;
    row_pipeline_separator_partitioned = separator_plan.partitioned ? 1 : 0;
    row_pipeline_separator_split_components =
      separator_plan.split_component_count;
  }

  kls_row_first_block_context row_ctx;
  memset(&row_ctx, 0, sizeof(row_ctx));
  row_ctx.solver = solver;
  row_ctx.numeric = shared->numeric;
  row_ctx.numeric_values = shared->numeric_values;
  row_ctx.row_ptr = worker->row_ptr;
  row_ctx.row_cols = row_cols;
  row_ctx.row_input_pos = row_input_pos;
  row_ctx.q_order = worker->block_q_order;
  row_ctx.col_pos = worker->col_pos;
  row_ctx.separator_component_last = shared->separator_component_last;
  row_ctx.separator_local_order = local_row_order;
  row_ctx.k1 = k1;
  row_ctx.nk = nk;
  row_ctx.n = n;
  row_ctx.separator_block_base = separator_block_base;
  row_ctx.pivot_col_limit = nk;
  row_ctx.scaled = shared->scaled;
  row_ctx.scale = shared->scale;
  row_ctx.use_separator_for_block = use_separator_for_block;
  row_ctx.strict_separator_pivot_scope =
    kls_row_first_strict_separator_pivot_scope_env_enabled();
  row_ctx.tol = shared->tol;

  kls_row_first_row_stats row_stats;
  memset(&row_stats, 0, sizeof(row_stats));
  if (block_pipeline_allowed) {
    UF_long pipeline_begin =
      shared->block_pipeline_begin <= nk ? shared->block_pipeline_begin : 0u;
    UF_long pipeline_end =
      shared->block_pipeline_end > pipeline_begin &&
          shared->block_pipeline_end <= nk
        ? shared->block_pipeline_end : nk;
    const unsigned char *active_mask = shared->block_pipeline_active_mask;
    const int use_active_mask =
      active_mask != NULL && shared->block_pipeline_active_rows > 0u;
    UF_long active_rows = 0;
    if (use_active_mask) {
      row_ctx.pivot_col_mask = active_mask;
      row_ctx.pivot_col_limit = nk;
      row_ctx.use_pivot_col_limit = 0;
    } else {
      row_ctx.pivot_col_limit = pipeline_end;
      row_ctx.use_pivot_col_limit = pipeline_end < nk;
    }
    if (row_order == NULL) {
      row_order = (UF_long *)malloc((size_t)nk * sizeof(*row_order));
    }
    row_done = (unsigned char *)calloc((size_t)nk, sizeof(*row_done));
    if (use_active_mask) {
      active_rank = (UF_long *)malloc((size_t)nk * sizeof(*active_rank));
    }
    if (row_order == NULL || row_done == NULL) {
      goto cleanup;
    }
    if (use_active_mask) {
      if (active_rank == NULL) {
        goto cleanup;
      }
      for (UF_long i = 0; i < nk; ++i) {
        active_rank[i] = KLS_KLU_EMPTY;
      }
      if (shared->block_pipeline_row_order != NULL &&
          shared->block_pipeline_order_rows > 0u) {
        if (shared->block_pipeline_order_rows !=
            shared->block_pipeline_active_rows) {
          goto cleanup;
        }
        for (UF_long pos = 0; pos < shared->block_pipeline_order_rows;
             ++pos) {
          const UF_long row = shared->block_pipeline_row_order[pos];
          if (row >= nk || row < pipeline_begin || row >= pipeline_end ||
              !active_mask[row] || active_rank[row] != KLS_KLU_EMPTY) {
            goto cleanup;
          }
          row_order[active_rows] = row;
          active_rank[row] = active_rows;
          active_rows++;
        }
      } else {
        for (UF_long i = 0; i < nk; ++i) {
          if (!active_mask[i]) {
            continue;
          }
          if (i < pipeline_begin || i >= pipeline_end) {
            goto cleanup;
          }
          row_order[active_rows] = i;
          active_rank[i] = active_rows;
          active_rows++;
        }
      }
      if (active_rows == 0u ||
          active_rows != shared->block_pipeline_active_rows ||
          pipeline_end <= pipeline_begin ||
          active_rows > pipeline_end - pipeline_begin) {
        goto cleanup;
      }
      trace_active_rows = active_rows;
      if (use_separator_for_block) {
        UF_long scoped_rows = 0;
        if (!kls_row_first_validate_tail_separator_scope(
              &row_ctx, row_order, active_rows, nk, active_mask,
              &scoped_rows) ||
            scoped_rows != active_rows) {
          goto cleanup;
        }
        shared->block_pipeline_separator_tail_scope = 1;
        shared->block_pipeline_separator_tail_scope_rows = scoped_rows;
      }
      row_pipeline_gap_rows = pipeline_end - pipeline_begin - active_rows;
      pipeline_begin = 0;
      pipeline_end = active_rows;
    } else {
      for (UF_long i = 0; i < nk; ++i) {
        row_order[i] = i;
      }
      if (use_separator_row_order) {
        row_owner = (int *)malloc((size_t)nk * sizeof(*row_owner));
        if (row_owner != NULL) {
          separator_private_owner_valid =
            kls_row_first_build_private_row_owners(
              &row_ctx, &separator_plan, shared->block_pipeline_thread_count,
              row_owner);
        }
        if (row_owner == NULL || !separator_private_owner_valid) {
          free(row_owner);
          row_owner = NULL;
          use_separator_row_order = 0;
          row_pipeline_separator_partitioned = 0;
          row_pipeline_separator_split_components = 0;
        } else if ((separator_private_run_ok =
                      kls_row_first_run_parallel_private_phase(
                        &row_ctx, &separator_plan,
                        shared->block_pipeline_thread_count,
                        row_owner, &l_entries, &u_entries,
                        worker->udiag_values, worker->u_row_ptr,
                        worker->u_row_end, row_done, &row_stats,
                        &row_pipeline_separator_private_rows,
                        &row_pipeline_separator_private_threads, 1,
                        &separator_private_complete))) {
          separator_private_executed = 1;
          row_pipeline_separator_queue = 1;
          if (separator_private_complete) {
            pipeline_begin = block_separator_private_rows;
            pipeline_end = nk;
          } else {
            private_active_mask =
              (unsigned char *)calloc((size_t)nk,
                                      sizeof(*private_active_mask));
            active_rank = (UF_long *)malloc((size_t)nk *
                                            sizeof(*active_rank));
            if (private_active_mask == NULL || active_rank == NULL) {
              goto cleanup;
            }
            active_rows = 0u;
            for (UF_long i = 0; i < nk; ++i) {
              active_rank[i] = KLS_KLU_EMPTY;
              if (row_done[i]) {
                continue;
              }
              private_active_mask[i] = 1u;
              row_order[active_rows] = i;
              active_rank[i] = active_rows;
              active_rows++;
            }
            if (active_rows == 0u ||
                active_rows > nk ||
                row_pipeline_separator_private_rows != nk - active_rows) {
              goto cleanup;
            }
            row_ctx.pivot_col_mask = private_active_mask;
            row_ctx.pivot_col_limit = nk;
            row_ctx.use_pivot_col_limit = 0;
            row_ctx.use_separator_for_block = 0;
            separator_private_partial_active = 1;
            pipeline_begin = 0u;
            pipeline_end = active_rows;
          }
        } else {
          free(row_owner);
          row_owner = NULL;
          use_separator_row_order = 0;
          separator_private_complete = 0;
          separator_private_partial_active = 0;
          row_pipeline_separator_private_rows = 0;
          row_pipeline_separator_private_threads = 0;
          row_pipeline_separator_partitioned = 0;
          row_pipeline_separator_split_components = 0;
        }
        if (kls_first_factor_trace_env_enabled()) {
          fprintf(stderr,
                  "KLS first-factor separator-private decision:"
                  " block=%" PRIu64 " rows=%" PRIu64
                  " use_order=%d owner_valid=%d run_ok=%d"
                  " executed=%d complete=%d partial_active=%d"
                  " private_rows=%" PRIu64
                  " pipeline_begin=%" PRIu64
                  " pipeline_end=%" PRIu64 "\n",
                  (uint64_t)block, (uint64_t)nk,
                  use_separator_row_order, separator_private_owner_valid,
                  separator_private_run_ok, separator_private_executed,
                  separator_private_complete,
                  separator_private_partial_active,
                  (uint64_t)row_pipeline_separator_private_rows,
                  (uint64_t)pipeline_begin, (uint64_t)pipeline_end);
        }
      }
    }
    row_pipeline_prefix_rows = use_active_mask
      ? shared->block_pipeline_begin : 0u;
    row_pipeline_suffix_rows = use_active_mask
      ? nk - shared->block_pipeline_end
      : (separator_private_partial_active ? 0u : nk - pipeline_end);
    UF_long pipeline_partial_rows = 0;
    UF_long pipeline_wait_partial_rows = 0;
    UF_long pipeline_wait_partial_deps = 0;
    UF_long pipeline_supernode_update_groups = 0;
    UF_long pipeline_supernode_update_rows = 0;
    UF_long pipeline_supernode_panel_update_groups = 0;
    UF_long pipeline_supernode_panel_update_rows = 0;
    UF_long pipeline_pivot_serial_rows = 0;
    UF_long pipeline_etree_ready_rows = 0;
    UF_long pipeline_etree_prefactor_rows = 0;
    UF_long pipeline_etree_prefactor_wait_rows = 0;
    UF_long pipeline_etree_prefactor_wait_deps = 0;
    trace_etree_attempted =
      use_active_mask && active_rows > 1u &&
      shared->block_pipeline_etree_tail_order;
    if (trace_etree_attempted) {
      trace_etree_succeeded = kls_try_row_first_etree_ready_tail_repair(
        worker, block, k1, nk, &row_ctx, row_order, active_mask,
        active_rows, shared->block_pipeline_begin,
        shared->block_pipeline_end, reserve_entries, &l_entries,
        &u_entries, &row_stats, &pipeline_etree_ready_rows,
        &row_pipeline_etree_ready_threads,
        &pipeline_etree_prefactor_rows,
        &pipeline_etree_prefactor_wait_rows,
        &pipeline_etree_prefactor_wait_deps,
        &pipeline_supernode_update_groups,
        &pipeline_supernode_update_rows,
        &pipeline_supernode_panel_update_groups,
        &pipeline_supernode_panel_update_rows,
        &row_pipeline_pivot_tail_rows,
        &row_pipeline_pivot_restarts);
    }
    if (trace_etree_succeeded) {
      row_pipeline_rows = active_rows;
      row_pipeline_threads = row_pipeline_etree_ready_threads;
      row_pipeline_etree_ready_rows = pipeline_etree_ready_rows;
      row_pipeline_etree_prefactor_rows = pipeline_etree_prefactor_rows;
      row_pipeline_etree_prefactor_threads =
        row_pipeline_etree_ready_threads;
      row_pipeline_etree_prefactor_wait_rows =
        pipeline_etree_prefactor_wait_rows;
      row_pipeline_etree_prefactor_wait_deps =
        pipeline_etree_prefactor_wait_deps;
    } else {
      kls_row_first_supernode_panel_cache_reset_active(
        &row_supernode_panel_cache, &row_workspace, nk, &open_panel_start,
        &open_panel_end);
      if (!separator_private_executed) {
        trace_seed_succeeded = kls_row_first_parallel_seed_preserved_rows(
              worker, block, k1, nk,
              use_active_mask ? shared->block_pipeline_begin : pipeline_begin,
              use_active_mask ? shared->block_pipeline_end : pipeline_end,
              active_mask,
              &l_entries, &u_entries,
              &row_workspace, row_done, &row_supernode_panel_cache,
              &open_panel_start, &open_panel_end);
        if (!trace_seed_succeeded) {
          goto cleanup;
        }
      }
      trace_pipeline_succeeded = kls_row_first_run_restartable_pipeline_suffix(
            &row_ctx, row_order, active_rank, pipeline_begin, pipeline_end,
            shared->block_pipeline_thread_count,
            &row_workspace, &l_entries, &u_entries, worker->udiag_values,
            row_done, &row_stats, &row_pipeline_rows, &row_pipeline_threads,
            &pipeline_partial_rows, &pipeline_wait_partial_rows,
            &pipeline_wait_partial_deps, &pipeline_supernode_update_groups,
            &pipeline_supernode_update_rows,
            &pipeline_supernode_panel_update_groups,
            &pipeline_supernode_panel_update_rows,
            &row_pipeline_pivot_tail_rows,
            &row_pipeline_pivot_restarts, &pipeline_pivot_serial_rows,
            row_workspace.supernode_panel_cache != NULL
              ? &row_supernode_panel_cache : NULL,
            &open_panel_start, &open_panel_end);
      if (!trace_pipeline_succeeded) {
        goto cleanup;
      }
      if (separator_private_executed) {
        row_pipeline_separator_pipeline_rows = row_pipeline_rows;
        row_pipeline_rows += row_pipeline_separator_private_rows;
        if (row_pipeline_separator_private_threads > row_pipeline_threads) {
          row_pipeline_threads = row_pipeline_separator_private_threads;
        }
      }
    }
    (void)pipeline_partial_rows;
    (void)pipeline_wait_partial_rows;
    (void)pipeline_wait_partial_deps;
    (void)pipeline_pivot_serial_rows;
    row_supernode_update_groups += pipeline_supernode_update_groups;
    row_supernode_update_rows += pipeline_supernode_update_rows;
    row_supernode_panel_update_groups +=
      pipeline_supernode_panel_update_groups;
    row_supernode_panel_update_rows += pipeline_supernode_panel_update_rows;
  } else {
    for (UF_long i = 0; i < nk; ++i) {
      const UF_long dynamic_pivots_before =
        row_stats.dynamic_column_pivots;
      const UF_long separator_rejects_before =
        row_stats.separator_dynamic_column_rejects;
      if (!kls_row_first_factor_one_row(
            &row_ctx, &row_workspace, &l_entries, &u_entries,
            worker->udiag_values, i, NULL, -1, NULL, &row_stats, NULL)) {
        if (row_stats.separator_dynamic_column_rejects >
            separator_rejects_before) {
          const UF_long reject_delta =
            row_stats.separator_dynamic_column_rejects -
            separator_rejects_before;
          atomic_fetch_add_explicit(
            &shared->separator_dynamic_column_rejects,
            (unsigned long)reject_delta, memory_order_relaxed);
        }
        goto cleanup;
      }
      worker->u_row_ptr[i + 1u] = u_entries.count;
      if (row_stats.dynamic_column_pivots != dynamic_pivots_before) {
        kls_row_first_supernode_panel_cache_reset_active(
          &row_supernode_panel_cache, &row_workspace, nk,
          &open_panel_start, &open_panel_end);
      }
      kls_row_first_supernode_panel_cache_publish_completed(
        row_workspace.supernode_panel_cache != NULL ? &row_supernode_panel_cache
                                                    : NULL,
        &u_entries, worker->u_row_ptr,
        worker->u_row_end, worker->udiag_values, nk, &row_stats,
        &open_panel_start,
        &open_panel_end, i);
    }
  }
  dynamic_column_pivots = row_stats.dynamic_column_pivots;
  separator_dynamic_column_pivots =
    row_stats.separator_dynamic_column_pivots;
  separator_extent_dynamic_column_pivots =
    row_stats.separator_extent_dynamic_column_pivots;
  separator_dynamic_column_fallbacks =
    row_stats.separator_dynamic_column_fallbacks;
  row_supernode_update_groups += row_stats.supernode_update_groups;
  row_supernode_update_rows += row_stats.supernode_update_rows;
  row_supernode_panel_update_groups +=
    row_stats.supernode_panel_update_groups;
  row_supernode_panel_update_rows +=
    row_stats.supernode_panel_update_rows;
  kls_row_first_stats_add(&row_panel_cache_builds,
                          row_stats.panel_cache_builds);
  kls_row_first_stats_add(&row_panel_cache_build_panels,
                          row_stats.panel_cache_build_panels);
  kls_row_first_stats_add(&row_panel_cache_build_entries,
                          row_stats.panel_cache_build_entries);
  kls_row_first_stats_add(&row_panel_cache_appends,
                          row_stats.panel_cache_appends);
  kls_row_first_stats_add(&row_panel_cache_append_panels,
                          row_stats.panel_cache_append_panels);
  kls_row_first_stats_add(&row_panel_cache_append_entries,
                          row_stats.panel_cache_append_entries);
  kls_row_first_stats_add(&active_rank_pivot_resets,
                          row_stats.active_rank_pivot_resets);
  kls_row_first_stats_add(&active_rank_pivot_reset_rows,
                          row_stats.active_rank_pivot_reset_rows);
  kls_row_first_stats_add(&active_rank_pivot_panel_rebuilds,
                          row_stats.active_rank_pivot_panel_rebuilds);
  kls_row_first_stats_add(&active_rank_pivot_panel_rebuild_rows,
                          row_stats.active_rank_pivot_panel_rebuild_rows);

  trace_commit_attempted = 1;
  ok = kls_row_first_parallel_commit_block(
    worker, block, k1, nk, worker->row_ptr, row_cols, row_input_pos,
    local_row_order, &l_entries, &u_entries, dynamic_column_pivots,
    separator_dynamic_column_pivots, separator_extent_dynamic_column_pivots,
    separator_dynamic_column_fallbacks, row_supernode_update_groups,
    row_supernode_update_rows, row_supernode_panel_update_groups,
    row_supernode_panel_update_rows, row_panel_cache_builds,
    row_panel_cache_build_panels, row_panel_cache_build_entries,
    row_panel_cache_appends, row_panel_cache_append_panels,
    row_panel_cache_append_entries, active_rank_pivot_resets,
    active_rank_pivot_reset_rows, active_rank_pivot_panel_rebuilds,
    active_rank_pivot_panel_rebuild_rows, row_pipeline_rows,
    row_pipeline_threads,
    row_pipeline_prefix_rows, row_pipeline_suffix_rows,
    row_pipeline_gap_rows,
    row_pipeline_pivot_tail_rows, row_pipeline_pivot_restarts,
    row_pipeline_etree_ready_rows, row_pipeline_etree_ready_threads,
    row_pipeline_etree_prefactor_rows,
    row_pipeline_etree_prefactor_threads,
    row_pipeline_etree_prefactor_wait_rows,
    row_pipeline_etree_prefactor_wait_deps,
    row_pipeline_separator_queue,
    row_pipeline_separator_private_rows,
    row_pipeline_separator_pipeline_rows,
    row_pipeline_separator_private_threads,
    row_pipeline_separator_partitioned,
    row_pipeline_separator_split_components);

cleanup:
  if (!ok && getenv("KLS_TRACE_FAST_REJECT_REPAIR") != NULL) {
    fprintf(stderr,
            "KLS row repair block stage: block=%ld active=%ld"
            " etree=%d/%d seed=%d pipeline=%d commit=%d"
            " L=%ld U=%ld pivots=%ld\n",
            (long)block, (long)trace_active_rows, trace_etree_attempted,
            trace_etree_succeeded, trace_seed_succeeded,
            trace_pipeline_succeeded, trace_commit_attempted,
            (long)l_entries.count, (long)u_entries.count,
            (long)row_stats.dynamic_column_pivots);
  }
  kls_row_first_supernode_panel_cache_free(&row_supernode_panel_cache);
  worker->supernode_workspace = row_workspace.supernode_workspace;
  worker->supernode_workspace_capacity =
    row_workspace.supernode_workspace_capacity;
  kls_row_first_entries_free(&l_entries);
  kls_row_first_entries_free(&u_entries);
  kls_first_separator_queue_plan_free(&separator_plan);
  free(row_cols);
  free(row_input_pos);
  free(row_next);
  free(row_order);
  free(local_row_order);
  free(old_to_new_row);
  free(row_seen);
  free(private_active_mask);
  free(active_rank);
  free(row_owner);
  free(row_done);
  for (UF_long local_col = 0; local_col < nk; ++local_col) {
    const UF_long oldcol = worker->block_q_order[local_col];
    if (oldcol < n) {
      worker->col_pos[oldcol] = KLS_KLU_EMPTY;
    }
  }
  return ok;
}

static int kls_row_first_repaired_prefix_unchanged(
  const kls_solver *solver,
  UF_long k1,
  UF_long prefix_cols,
  const UF_long *saved_pnum,
  const UF_long *saved_q_block,
  const UF_long *q_order) {
  if (solver == NULL || solver->numeric == NULL ||
      solver->numeric->Pnum == NULL || saved_pnum == NULL ||
      saved_q_block == NULL || q_order == NULL ||
      k1 > solver->n || prefix_cols > solver->n - k1) {
    return 0;
  }
  for (UF_long k = 0; k < prefix_cols; ++k) {
    if (solver->numeric->Pnum[k1 + k] != saved_pnum[k] ||
        q_order[k1 + k] != saved_q_block[k]) {
      return 0;
    }
  }
  return 1;
}

static int kls_try_row_first_rebuild_rejected_block(
  kls_solver *solver,
  double *numeric_values,
  UF_long block,
  UF_long k1,
  UF_long nk,
  UF_long rejected_pivot,
  const UF_long *old_pblock,
  UF_long old_lnz_block,
  UF_long old_unz_block,
  UF_long *pblock_out) {
  if (solver == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL || numeric_values == NULL ||
      old_pblock == NULL || pblock_out == NULL ||
      solver->symbolic->P == NULL || solver->symbolic->Q == NULL ||
      solver->symbolic->R == NULL || solver->symbolic->Lnz == NULL ||
      solver->numeric->Pinv == NULL || solver->numeric->Pnum == NULL ||
      solver->numeric->Lip == NULL || solver->numeric->Llen == NULL ||
      solver->numeric->Uip == NULL || solver->numeric->Ulen == NULL ||
      solver->numeric->Udiag == NULL || solver->numeric->LUbx == NULL ||
      solver->numeric->LUsize == NULL || block >= solver->symbolic->nblocks ||
      k1 > solver->n || nk == 0u || nk > solver->n - k1 ||
      rejected_pivot < k1 || rejected_pivot >= k1 + nk ||
      nk > solver->symbolic->maxblock ||
      (solver->common.scale > 0 && solver->numeric->Rs == NULL) ||
      (solver->common.scale <= 0 && solver->numeric->Rs != NULL) ||
      solver->numeric->lnz < old_lnz_block ||
      solver->numeric->unz < old_unz_block) {
    if (getenv("KLS_TRACE_FAST_REJECT_REPAIR") != NULL) {
      fprintf(stderr,
              "KLS rejected-block row repair ineligible: block=%ld k1=%ld"
              " nk=%ld reject=%ld n=%ld nblocks=%ld maxblock=%ld"
              " scale=%ld Rs=%d lnz=%ld/%ld unz=%ld/%ld\n",
              (long)block, (long)k1, (long)nk, (long)rejected_pivot,
              (long)(solver != NULL ? solver->n : 0u),
              (long)(solver != NULL && solver->symbolic != NULL
                       ? solver->symbolic->nblocks : 0u),
              (long)(solver != NULL && solver->symbolic != NULL
                       ? solver->symbolic->maxblock : 0u),
              (long)(solver != NULL ? solver->common.scale : 0L),
              solver != NULL && solver->numeric != NULL &&
                solver->numeric->Rs != NULL,
              (long)(solver != NULL && solver->numeric != NULL
                       ? solver->numeric->lnz : 0u),
              (long)old_lnz_block,
              (long)(solver != NULL && solver->numeric != NULL
                       ? solver->numeric->unz : 0u),
              (long)old_unz_block);
    }
    return 0;
  }

  const UF_long n = solver->n;
  const int scaled = solver->common.scale > 0;
  const size_t nk_size = (size_t)nk;
  if ((UF_long)nk_size != nk) {
    return 0;
  }

  UF_long *saved_ptrs = NULL;
  UF_long *saved_pnum = NULL;
  UF_long *saved_q_block = NULL;
  Entry *saved_udiag = NULL;
  UF_long *q_order = NULL;
  UF_long *separator_component_last = NULL;
  unsigned char *pipeline_active_mask = NULL;
  UF_long *pipeline_row_order = NULL;
  int mutex_initialized = 0;
  int snapshot_ready = 0;
  Unit *old_lu = NULL;
  size_t old_lusize = 0u;
  UF_long saved_lnz = 0u;
  UF_long saved_unz = 0u;
  UF_long saved_max_lnz = 0u;
  UF_long saved_max_unz = 0u;
  trilinos_klu_l_common saved_common;
  memset(&saved_common, 0, sizeof(saved_common));
  kls_row_first_parallel_worker worker;
  memset(&worker, 0, sizeof(worker));
  kls_row_first_parallel_shared shared;
  memset(&shared, 0, sizeof(shared));
  shared.block_pipeline_target_block = KLS_KLU_EMPTY;
  atomic_init(&shared.next_block, 0ul);
  atomic_init(&shared.failed, 0);
  atomic_init(&shared.separator_dynamic_column_rejects, 0ul);

  if (nk_size > SIZE_MAX / (4u * sizeof(*saved_ptrs))) {
    return 0;
  }
  saved_ptrs = (UF_long *)malloc(4u * nk_size * sizeof(*saved_ptrs));
  saved_pnum = (UF_long *)malloc(nk_size * sizeof(*saved_pnum));
  saved_q_block = (UF_long *)malloc(nk_size * sizeof(*saved_q_block));
  saved_udiag = (Entry *)malloc(nk_size * sizeof(*saved_udiag));
  q_order = n > 0u
    ? (UF_long *)malloc((size_t)n * sizeof(*q_order)) : NULL;
  if (saved_ptrs == NULL || saved_pnum == NULL ||
      saved_q_block == NULL || saved_udiag == NULL ||
      (n > 0u && q_order == NULL)) {
    goto fail;
  }

  UF_long *saved_lip = saved_ptrs;
  UF_long *saved_llen = saved_ptrs + nk;
  UF_long *saved_uip = saved_ptrs + 2u * nk;
  UF_long *saved_ulen = saved_ptrs + 3u * nk;
  memcpy(saved_lip, solver->numeric->Lip + k1,
         nk_size * sizeof(*saved_lip));
  memcpy(saved_llen, solver->numeric->Llen + k1,
         nk_size * sizeof(*saved_llen));
  memcpy(saved_uip, solver->numeric->Uip + k1,
         nk_size * sizeof(*saved_uip));
  memcpy(saved_ulen, solver->numeric->Ulen + k1,
         nk_size * sizeof(*saved_ulen));
  memcpy(saved_pnum, solver->numeric->Pnum + k1,
         nk_size * sizeof(*saved_pnum));
  memcpy(saved_q_block, solver->symbolic->Q + k1,
         nk_size * sizeof(*saved_q_block));
  memcpy(saved_udiag, ((Entry *)solver->numeric->Udiag) + k1,
         nk_size * sizeof(*saved_udiag));
  for (UF_long k = 0; k < n; ++k) {
    const UF_long oldcol = solver->symbolic->Q[k];
    if (oldcol >= n) {
      goto fail;
    }
    q_order[k] = oldcol;
  }

  old_lu = solver->numeric->LUbx[block];
  old_lusize = solver->numeric->LUsize[block];
  saved_lnz = solver->numeric->lnz;
  saved_unz = solver->numeric->unz;
  saved_max_lnz = solver->numeric->max_lnz_block;
  saved_max_unz = solver->numeric->max_unz_block;
  saved_common = solver->common;
  snapshot_ready = 1;

  UF_long *psinv = solver->numeric->Pinv;
  for (UF_long k = 0; k < n; ++k) {
    psinv[k] = KLS_KLU_EMPTY;
  }
  for (UF_long k = 0; k < n; ++k) {
    const UF_long row = solver->symbolic->P[k];
    if (row >= n || psinv[row] != KLS_KLU_EMPTY) {
      goto fail;
    }
    psinv[row] = k;
  }
  for (UF_long k = 0; k < n; ++k) {
    if (psinv[k] == KLS_KLU_EMPTY) {
      goto fail;
    }
  }

  const int use_separator_pivot_domains =
    kls_separator_analysis_has_global_range(&solver->separator) &&
    solver->separator.component_count > 0u &&
    solver->separator.order_component != NULL;
  const int separator_covers_rejected_block =
    use_separator_pivot_domains &&
    k1 >= solver->separator.global_begin &&
    k1 + nk <= solver->separator.global_end;
  if (use_separator_pivot_domains) {
    if (solver->separator.component_count >
        (UF_long)(SIZE_MAX / sizeof(*separator_component_last))) {
      goto fail;
    }
    separator_component_last =
      (UF_long *)malloc((size_t)solver->separator.component_count *
                        sizeof(*separator_component_last));
    if (separator_component_last == NULL) {
      goto fail;
    }
    for (UF_long component = 0; component < solver->separator.component_count;
         ++component) {
      separator_component_last[component] = KLS_KLU_EMPTY;
    }
    for (UF_long pos = 0; pos < solver->separator.n; ++pos) {
      const unsigned int component = solver->separator.order_component[pos];
      if ((UF_long)component >= solver->separator.component_count) {
        goto fail;
      }
      separator_component_last[component] = pos;
    }
  }

  if (pthread_mutex_init(&shared.commit_lock, NULL) != 0) {
    goto fail;
  }
  mutex_initialized = 1;
  shared.solver = solver;
  shared.numeric = solver->numeric;
  shared.numeric_values = numeric_values;
  shared.psinv = psinv;
  shared.q_order = q_order;
  shared.separator_component_last = separator_component_last;
  shared.scaled = scaled;
  shared.scale = scaled ? (int)solver->common.scale : 0;
  shared.use_separator_pivot_domains = use_separator_pivot_domains;
  shared.force_block_pipeline = solver->options.threads > 1 && nk > 1u;
  shared.block_pipeline_thread_count = solver->options.threads;
  const UF_long local_reject = rejected_pivot - k1;
  UF_long tail_begin = 0;
  UF_long tail_end = 0;
  UF_long tail_columns = 0;
  int tail_suffix_exact = 0;
  UF_long repaired_tail_begin = local_reject;
  int repaired_tail_etree_mask = 0;
  int repaired_tail_etree_order = 0;
  if (shared.force_block_pipeline &&
      kls_pivot_tail_plan_local_range(
        solver, k1, nk, local_reject, &tail_begin, &tail_end,
        &tail_columns, NULL, &tail_suffix_exact) &&
      (!separator_covers_rejected_block || tail_suffix_exact)) {
    const size_t nk_size = (size_t)nk;
    if ((UF_long)nk_size == nk) {
      pipeline_active_mask =
        (unsigned char *)malloc(nk_size * sizeof(*pipeline_active_mask));
      pipeline_row_order =
        (UF_long *)malloc(nk_size * sizeof(*pipeline_row_order));
    }
    const int mask_ok =
      pipeline_active_mask != NULL && pipeline_row_order != NULL &&
      kls_build_pivot_tail_restart_mask_range(
        solver, k1, nk, local_reject, tail_begin, tail_end,
        pipeline_active_mask);
    const int mask_etree_ok =
      mask_ok &&
      kls_pivot_tail_restart_mask_matches_etree_forest(
        solver, block, k1, nk, local_reject, tail_begin, tail_end,
        pipeline_active_mask);
    UF_long active_rows = 0;
    const int order_ok =
      mask_etree_ok &&
      kls_build_pivot_tail_restart_order_range(
        solver, k1, nk, local_reject, tail_begin, tail_end,
        pipeline_active_mask, pipeline_row_order, &active_rows);
    const int refresh_ok =
      order_ok && active_rows == (UF_long)tail_columns &&
      active_rows > 0u &&
      (solver->fast_reject_refresh_state !=
         KLS_FAST_REJECT_REFRESH_PREFIX ||
       tail_suffix_exact ||
       kls_refresh_pivot_tail_preserved_block_columns(
         solver, numeric_values, block, k1, nk, pipeline_active_mask));
    if (refresh_ok) {
      shared.block_pipeline_begin = tail_begin;
      shared.block_pipeline_end = tail_end;
      repaired_tail_begin = tail_begin;
      repaired_tail_etree_mask = 1;
      repaired_tail_etree_order = 1;
      shared.block_pipeline_active_mask = pipeline_active_mask;
      shared.block_pipeline_active_rows = active_rows;
      shared.block_pipeline_row_order = pipeline_row_order;
      shared.block_pipeline_order_rows = active_rows;
      shared.block_pipeline_etree_tail_order = 1;
    } else {
      free(pipeline_active_mask);
      free(pipeline_row_order);
      pipeline_active_mask = NULL;
      pipeline_row_order = NULL;
    }
  }
  shared.tol = solver->common.tol;

  if (!kls_row_first_parallel_worker_init(&worker, &shared)) {
    goto fail;
  }

  solver->common.status = TRILINOS_KLU_OK;
  solver->common.numerical_rank = KLS_KLU_EMPTY;
  solver->common.singular_col = KLS_KLU_EMPTY;
  solver->common.noffdiag = saved_common.noffdiag;
  if (!kls_row_first_parallel_factor_block(&worker, block) ||
      atomic_load_explicit(&shared.failed, memory_order_acquire) ||
      shared.committed_blocks != 1u ||
      solver->common.status < TRILINOS_KLU_OK ||
      solver->common.status == TRILINOS_KLU_SINGULAR ||
      solver->numeric->LUbx[block] == NULL ||
      shared.total_lnz == 0u || shared.total_unz == 0u) {
    goto fail;
  }

  for (UF_long k = 0; k < nk; ++k) {
    const UF_long row = solver->numeric->Pnum[k1 + k];
    if (row >= n || psinv[row] < k1 || psinv[row] >= k1 + nk) {
      goto fail;
    }
    pblock_out[k] = psinv[row] - k1;
  }
  memcpy(solver->symbolic->Q + k1, q_order + k1,
         nk_size * sizeof(*q_order));

  if (!kls_rebuild_numeric_pinv(solver)) {
    goto fail;
  }
  int offdiag_refreshed = 0;
  solver->fast_repaired_last_offdiag_suffix_refresh = 0;
  if (kls_row_first_repaired_prefix_unchanged(
        solver, k1, repaired_tail_begin, saved_pnum, saved_q_block,
        q_order)) {
    offdiag_refreshed =
      kls_recompute_offdiag_suffix_from_pinv(solver, numeric_values,
                                             k1 + repaired_tail_begin);
    if (offdiag_refreshed) {
      solver->fast_repaired_last_offdiag_suffix_refresh = 1;
      solver->fast_repaired_offdiag_suffix_refresh_count++;
    }
  }
  if (!offdiag_refreshed) {
    offdiag_refreshed =
      kls_recompute_offdiag_from_pinv(solver, numeric_values);
    if (offdiag_refreshed) {
      solver->fast_repaired_offdiag_full_refresh_count++;
    }
  }
  if (!offdiag_refreshed) {
    goto fail;
  }
  if (scaled && !kls_parallel_refactor_permute_scale(solver)) {
    goto fail;
  }

  solver->numeric->lnz = saved_lnz - old_lnz_block + shared.total_lnz;
  solver->numeric->unz = saved_unz - old_unz_block + shared.total_unz;
  solver->numeric->max_lnz_block =
    saved_max_lnz < shared.max_lnz_block
      ? shared.max_lnz_block : saved_max_lnz;
  solver->numeric->max_unz_block =
    saved_max_unz < shared.max_unz_block
      ? shared.max_unz_block : saved_max_unz;
  if (solver->symbolic->Lnz[block] < 0.0) {
    solver->symbolic->Lnz[block] =
      (double)(shared.total_lnz > shared.total_unz
                 ? shared.total_lnz : shared.total_unz);
  }
  UF_long offdiag = 0;
  for (UF_long k = 0; k < n; ++k) {
    if (solver->numeric->Pnum[k] != solver->symbolic->P[k]) {
      offdiag++;
    }
  }
  solver->common.noffdiag = offdiag;
  if (shared.block_pipeline_active_mask != NULL) {
    UF_long next_rejected_pivot = KLS_KLU_EMPTY;
    UF_long rejected_pivot_col = KLS_KLU_EMPTY;
    UF_long rejected_row = KLS_KLU_EMPTY;
    double rejected_ratio_abs = -1.0;
    double rejected_pivot_abs = -1.0;
    double rejected_row_max_abs = -1.0;
    if (!kls_numeric_pivots_pass_threshold_from_block(
          solver, block, &next_rejected_pivot, &rejected_pivot_col) ||
        kls_numeric_rowwise_u_pivots_pass_threshold_from_block(
          solver, block, &next_rejected_pivot, &rejected_pivot_col,
          &rejected_row, &rejected_ratio_abs, &rejected_pivot_abs,
          &rejected_row_max_abs) != 1) {
      goto fail;
    }
  }

  UF_long repaired_tail_end = nk;
  if (shared.block_pipeline_end > 0u && shared.block_pipeline_end <= nk) {
    repaired_tail_end = shared.block_pipeline_end;
  }
  UF_long repaired_tail_skipped_columns = 0;
  double repaired_tail_skipped_work = 0.0;
  int repaired_tail_exact_mask = 0;
  if (shared.block_pipeline_active_mask != NULL) {
    repaired_tail_end = nk;
    for (UF_long k = repaired_tail_begin; k < nk; ++k) {
      if (shared.block_pipeline_active_mask[k]) {
        continue;
      }
      double column_work = 0.0;
      if (!kls_fast_reject_column_work(solver, block, k1, nk, k,
                                       &column_work)) {
        goto fail;
      }
      repaired_tail_skipped_columns++;
      repaired_tail_skipped_work += column_work;
    }
    repaired_tail_exact_mask = 1;
  } else if (repaired_tail_end < nk) {
    for (UF_long k = repaired_tail_end; k < nk; ++k) {
      double column_work = 0.0;
      if (!kls_fast_reject_column_work(solver, block, k1, nk, k,
                                       &column_work)) {
        goto fail;
      }
      repaired_tail_skipped_columns++;
      repaired_tail_skipped_work += column_work;
    }
    repaired_tail_exact_mask = 1;
  }
  kls_record_fast_repaired_block_stats(
    solver, block, k1, nk, rejected_pivot, repaired_tail_begin,
    repaired_tail_end,
    repaired_tail_skipped_columns, repaired_tail_skipped_work,
    repaired_tail_exact_mask, repaired_tail_etree_mask, old_pblock,
    pblock_out);
  if (old_lu != NULL) {
    (void)TRILINOS_KLU_free(old_lu, old_lusize, sizeof(Unit),
                            &solver->common);
  }
  kls_invalidate_i32_solve(solver);
  solver->fast_block_restarts++;
  solver->fast_kls_block_restarts++;
  if (shared.row_pipeline_runs > 0u) {
    solver->fast_kls_block_restart_last_row_pipeline = 1;
    solver->fast_kls_block_restart_row_pipeline_count +=
      shared.row_pipeline_runs;
    solver->fast_kls_block_restart_last_row_pipeline_rows =
      shared.row_pipeline_rows;
    solver->fast_kls_block_restart_last_row_pipeline_threads =
      shared.row_pipeline_threads;
    solver->fast_kls_block_restart_last_row_pipeline_prefix_rows =
      shared.row_pipeline_prefix_rows;
    solver->fast_kls_block_restart_last_row_pipeline_suffix_rows =
      shared.row_pipeline_suffix_rows;
    solver->fast_kls_block_restart_last_row_pipeline_gap_rows =
      shared.row_pipeline_gap_rows;
    if (repaired_tail_etree_order &&
        shared.block_pipeline_etree_tail_order) {
      solver->fast_kls_block_restart_last_row_pipeline_etree_tail = 1;
      solver->fast_kls_block_restart_row_pipeline_etree_tail_count +=
        shared.row_pipeline_runs;
      solver->fast_kls_block_restart_last_row_pipeline_etree_tail_rows =
        shared.row_pipeline_rows;
      solver->fast_kls_block_restart_last_row_pipeline_etree_tail_gap_rows =
        shared.row_pipeline_gap_rows;
      solver
        ->fast_kls_block_restart_last_row_pipeline_etree_tail_exact_mask =
          repaired_tail_exact_mask;
    }
    if (shared.row_pipeline_etree_ready_runs > 0u) {
      solver->fast_kls_block_restart_last_row_pipeline_etree_ready = 1;
      solver->fast_kls_block_restart_row_pipeline_etree_ready_count +=
        shared.row_pipeline_etree_ready_runs;
      solver->fast_kls_block_restart_last_row_pipeline_etree_ready_rows =
        shared.row_pipeline_etree_ready_rows;
      solver->fast_kls_block_restart_last_row_pipeline_etree_ready_threads =
        shared.row_pipeline_etree_ready_threads;
    }
    if (shared.row_pipeline_etree_prefactor_runs > 0u) {
      solver->fast_kls_block_restart_last_row_pipeline_etree_prefactor = 1;
      solver->fast_kls_block_restart_row_pipeline_etree_prefactor_count +=
        shared.row_pipeline_etree_prefactor_runs;
      solver->fast_kls_block_restart_last_row_pipeline_etree_prefactor_rows =
        shared.row_pipeline_etree_prefactor_rows;
      solver
        ->fast_kls_block_restart_last_row_pipeline_etree_prefactor_threads =
          shared.row_pipeline_etree_prefactor_threads;
      solver
        ->fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_rows =
          shared.row_pipeline_etree_prefactor_wait_rows;
      solver
        ->fast_kls_block_restart_last_row_pipeline_etree_prefactor_wait_deps =
          shared.row_pipeline_etree_prefactor_wait_deps;
    }
    if (shared.block_pipeline_separator_tail_scope) {
      solver
        ->fast_kls_block_restart_last_row_pipeline_separator_tail_scope = 1;
      solver
        ->fast_kls_block_restart_row_pipeline_separator_tail_scope_count +=
          shared.row_pipeline_runs;
      solver
        ->fast_kls_block_restart_last_row_pipeline_separator_tail_scope_rows =
          shared.block_pipeline_separator_tail_scope_rows;
    }
    if (shared.row_pipeline_separator_queue_runs > 0u) {
      solver->fast_kls_block_restart_last_row_pipeline_separator_queue = 1;
      solver->fast_kls_block_restart_row_pipeline_separator_queue_count +=
        shared.row_pipeline_separator_queue_runs;
      solver
        ->fast_kls_block_restart_last_row_pipeline_separator_private_rows =
          shared.row_pipeline_separator_private_rows;
      solver
        ->fast_kls_block_restart_last_row_pipeline_separator_pipeline_rows =
          shared.row_pipeline_separator_pipeline_rows;
      solver
        ->fast_kls_block_restart_last_row_pipeline_separator_private_threads =
          shared.row_pipeline_separator_private_threads;
      solver
        ->fast_kls_block_restart_last_row_pipeline_separator_partitioned =
          shared.row_pipeline_separator_partitioned;
      solver
        ->fast_kls_block_restart_last_row_pipeline_separator_split_components =
          shared.row_pipeline_separator_split_components;
    }
    solver->fast_kls_block_restart_last_row_pipeline_pivot_tail_rows =
      shared.row_pipeline_pivot_tail_rows;
    solver->fast_kls_block_restart_last_row_pipeline_pivot_restarts =
      shared.row_pipeline_pivot_restarts;
    solver
      ->fast_kls_block_restart_last_row_pipeline_supernode_update_groups =
        shared.row_supernode_update_groups;
    solver->fast_kls_block_restart_last_row_pipeline_supernode_update_rows =
      shared.row_supernode_update_rows;
    solver
      ->fast_kls_block_restart_last_row_pipeline_supernode_panel_update_groups =
        shared.row_supernode_panel_update_groups;
    solver
      ->fast_kls_block_restart_last_row_pipeline_supernode_panel_update_rows =
        shared.row_supernode_panel_update_rows;
  }
  solver->stats.fast_block_restarts = solver->fast_block_restarts;
  solver->stats.fast_kls_block_restarts = solver->fast_kls_block_restarts;
  solver->stats.fast_kls_rebuild_restarts =
    solver->fast_kls_rebuild_restarts;
  kls_fill_fast_kls_block_restart_pipeline_stats(solver);
  solver->stats.fast_tail_restarts = solver->fast_tail_restarts;
  solver->stats.fast_repaired_last_offdiag_suffix_refresh =
    solver->fast_repaired_last_offdiag_suffix_refresh;
  solver->stats.fast_repaired_offdiag_suffix_refresh_count =
    (int64_t)solver->fast_repaired_offdiag_suffix_refresh_count;
  solver->stats.fast_repaired_offdiag_full_refresh_count =
    (int64_t)solver->fast_repaired_offdiag_full_refresh_count;
  free_refactor_map(solver);
  free_refactor_schedule(solver);
  free_row_refactor_pattern_preserve_diagnostics(solver);
  free_refactor_lu_pointer_cache(solver);
  solver->common.status = TRILINOS_KLU_OK;
  solver->common.numerical_rank = KLS_KLU_EMPTY;
  solver->common.singular_col = KLS_KLU_EMPTY;

  kls_row_first_parallel_worker_free(&worker);
  if (mutex_initialized) {
    (void)pthread_mutex_destroy(&shared.commit_lock);
  }
  free(pipeline_active_mask);
  free(pipeline_row_order);
  free(saved_ptrs);
  free(saved_pnum);
  free(saved_q_block);
  free(saved_udiag);
  free(q_order);
  free(separator_component_last);
  return 1;

fail:
  if (getenv("KLS_TRACE_FAST_REJECT_REPAIR") != NULL) {
    fprintf(stderr,
            "KLS rejected-block row repair failed: block=%ld k1=%ld nk=%ld"
            " reject=%ld status=%ld shared_failed=%d committed=%ld"
            " lnz=%ld unz=%ld etree_mask=%d active_rows=%ld"
            " pipeline_runs=%ld pivot_restarts=%ld\n",
            (long)block, (long)k1, (long)nk, (long)rejected_pivot,
            (long)(solver != NULL ? solver->common.status : -999L),
            atomic_load_explicit(&shared.failed, memory_order_acquire),
            (long)shared.committed_blocks, (long)shared.total_lnz,
            (long)shared.total_unz,
            shared.block_pipeline_active_mask != NULL,
            (long)shared.block_pipeline_active_rows,
            (long)shared.row_pipeline_runs,
            (long)shared.row_pipeline_pivot_restarts);
  }
  if (snapshot_ready) {
    Unit *current_lu = solver != NULL && solver->numeric != NULL &&
        solver->numeric->LUbx != NULL && block < solver->symbolic->nblocks
      ? solver->numeric->LUbx[block] : NULL;
    size_t current_lusize = solver != NULL && solver->numeric != NULL &&
        solver->numeric->LUsize != NULL && block < solver->symbolic->nblocks
      ? solver->numeric->LUsize[block] : 0u;
    if (current_lu != NULL && current_lu != old_lu) {
      (void)TRILINOS_KLU_free(current_lu, current_lusize, sizeof(Unit),
                              &solver->common);
    }
    if (saved_ptrs != NULL) {
      memcpy(solver->numeric->Lip + k1, saved_ptrs,
             nk_size * sizeof(*solver->numeric->Lip));
      memcpy(solver->numeric->Llen + k1, saved_ptrs + nk,
             nk_size * sizeof(*solver->numeric->Llen));
      memcpy(solver->numeric->Uip + k1, saved_ptrs + 2u * nk,
             nk_size * sizeof(*solver->numeric->Uip));
      memcpy(solver->numeric->Ulen + k1, saved_ptrs + 3u * nk,
             nk_size * sizeof(*solver->numeric->Ulen));
    }
    if (saved_pnum != NULL) {
      memcpy(solver->numeric->Pnum + k1, saved_pnum,
             nk_size * sizeof(*solver->numeric->Pnum));
    }
    if (saved_q_block != NULL) {
      memcpy(solver->symbolic->Q + k1, saved_q_block,
             nk_size * sizeof(*solver->symbolic->Q));
    }
    if (saved_udiag != NULL) {
      memcpy(((Entry *)solver->numeric->Udiag) + k1, saved_udiag,
             nk_size * sizeof(*saved_udiag));
    }
    solver->numeric->LUbx[block] = old_lu;
    solver->numeric->LUsize[block] = old_lusize;
    solver->numeric->lnz = saved_lnz;
    solver->numeric->unz = saved_unz;
    solver->numeric->max_lnz_block = saved_max_lnz;
    solver->numeric->max_unz_block = saved_max_unz;
    solver->common = saved_common;
    (void)kls_rebuild_numeric_pinv(solver);
    (void)kls_recompute_offdiag_from_pinv(solver, numeric_values);
    solver->common = saved_common;
  }

  kls_row_first_parallel_worker_free(&worker);
  if (mutex_initialized) {
    (void)pthread_mutex_destroy(&shared.commit_lock);
  }
  free(pipeline_active_mask);
  free(pipeline_row_order);
  free(saved_ptrs);
  free(saved_pnum);
  free(saved_q_block);
  free(saved_udiag);
  free(q_order);
  free(separator_component_last);
  return 0;
}

static void *kls_row_first_parallel_worker_main(void *arg) {
  kls_row_first_parallel_worker worker;
  kls_row_first_parallel_shared *shared =
    (kls_row_first_parallel_shared *)arg;
  if (!kls_row_first_parallel_worker_init(&worker, shared)) {
    if (shared != NULL) {
      atomic_store_explicit(&shared->failed, 1, memory_order_release);
    }
    return NULL;
  }
  const UF_long nblocks =
    shared != NULL && shared->solver != NULL && shared->solver->symbolic != NULL
      ? shared->solver->symbolic->nblocks : 0u;
  for (;;) {
    if (atomic_load_explicit(&shared->failed, memory_order_acquire)) {
      break;
    }
    const UF_long block = (UF_long)atomic_fetch_add_explicit(
      &shared->next_block, 1ul, memory_order_acq_rel);
    if (block >= nblocks) {
      break;
    }
    if (!kls_row_first_parallel_factor_block(&worker, block)) {
      atomic_store_explicit(&shared->failed, 1, memory_order_release);
      break;
    }
  }
  kls_row_first_parallel_worker_free(&worker);
  return NULL;
}

static int kls_first_factor_dominant_btf_pipeline_target(
  const kls_solver *solver,
  UF_long *block_out) {
  if (block_out != NULL) {
    *block_out = KLS_KLU_EMPTY;
  }
  if (solver == NULL || solver->symbolic == NULL ||
      solver->symbolic->R == NULL || solver->n == 0u ||
      solver->options.threads <= 1 || solver->symbolic->nblocks < 2u ||
      solver->symbolic->maxblock <
        KLS_FIRST_FACTOR_DOMINANT_BTF_PIPELINE_MIN_BLOCK) {
    return 0;
  }

  const double coverage =
    (double)solver->symbolic->maxblock / (double)solver->n;
  if (coverage < KLS_FIRST_FACTOR_DOMINANT_BTF_PIPELINE_MIN_COVERAGE) {
    return 0;
  }

  for (UF_long block = 0; block < solver->symbolic->nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    if (k2 >= k1 && k2 - k1 == solver->symbolic->maxblock) {
      if (block_out != NULL) {
        *block_out = block;
      }
      return 1;
    }
  }
  return 0;
}

static int kls_try_first_factor_row_uplooking_blocks_parallel(
  kls_solver *solver,
  double *numeric_values) {
  if (solver == NULL || solver->symbolic == NULL || numeric_values == NULL ||
      solver->numeric != NULL || solver->symbolic->nblocks < 2u ||
      solver->options.threads <= 1 ||
      solver->symbolic->P == NULL || solver->symbolic->Q == NULL ||
      solver->symbolic->R == NULL || solver->symbolic->maxblock == 0u ||
      solver->n == 0u) {
    return -1;
  }
  int thread_count = solver->options.threads;
  if ((UF_long)thread_count > solver->symbolic->nblocks) {
    thread_count = (int)solver->symbolic->nblocks;
  }
  if (thread_count < 2) {
    return -1;
  }

  TRILINOS_KLU_common *common = &solver->common;
  trilinos_klu_l_common saved_common = *common;
  int q_committed = 0;
  common->status = TRILINOS_KLU_OK;
  common->numerical_rank = KLS_KLU_EMPTY;
  common->singular_col = KLS_KLU_EMPTY;
  common->noffdiag = 0;
  common->nrealloc = 0;

  trilinos_klu_l_numeric *numeric = kls_allocate_numeric_skeleton(solver);
  if (numeric == NULL) {
    *common = saved_common;
    return 0;
  }
  solver->numeric = numeric;

  UF_long *q_order = NULL;
  UF_long *saved_q = NULL;
  UF_long *separator_component_last = NULL;
  kls_row_first_refactor_seed row_refactor_seed;
  memset(&row_refactor_seed, 0, sizeof(row_refactor_seed));
  pthread_t *threads = NULL;
  int mutex_initialized = 0;
  kls_row_first_parallel_shared shared;
  memset(&shared, 0, sizeof(shared));
  shared.block_pipeline_target_block = KLS_KLU_EMPTY;
  atomic_init(&shared.next_block, 0ul);
  atomic_init(&shared.failed, 0);
  atomic_init(&shared.separator_dynamic_column_rejects, 0ul);

  const UF_long n = solver->n;
  const int scaled = common->scale > 0;
  if (scaled &&
      !trilinos_klu_l_scale((UF_long)common->scale, n,
                            solver->col_ptr, solver->row_idx,
                            numeric_values, numeric->Rs, NULL, common)) {
    goto fail;
  }

  UF_long *psinv = numeric->Pinv;
  for (UF_long k = 0; k < n; ++k) {
    psinv[k] = KLS_KLU_EMPTY;
  }
  for (UF_long k = 0; k < n; ++k) {
    const UF_long row = solver->symbolic->P[k];
    if (row >= n || psinv[row] != KLS_KLU_EMPTY) {
      goto fail;
    }
    psinv[row] = k;
  }
  for (UF_long k = 0; k < n; ++k) {
    if (psinv[k] == KLS_KLU_EMPTY) {
      goto fail;
    }
  }

  q_order = (UF_long *)malloc((size_t)n * sizeof(*q_order));
  saved_q = (UF_long *)malloc((size_t)n * sizeof(*saved_q));
  if (q_order == NULL || saved_q == NULL) {
    goto fail;
  }
  for (UF_long k = 0; k < n; ++k) {
    const UF_long oldcol = solver->symbolic->Q[k];
    if (oldcol >= n) {
      goto fail;
    }
    q_order[k] = oldcol;
    saved_q[k] = oldcol;
  }

  const int use_separator_pivot_domains =
    kls_separator_analysis_has_global_range(&solver->separator) &&
    solver->separator.component_count > 0u &&
    solver->separator.order_component != NULL;
  if (use_separator_pivot_domains) {
    if (solver->separator.component_count >
        (UF_long)(SIZE_MAX / sizeof(*separator_component_last))) {
      goto fail;
    }
    separator_component_last =
      (UF_long *)malloc((size_t)solver->separator.component_count *
                        sizeof(*separator_component_last));
    if (separator_component_last == NULL) {
      goto fail;
    }
    for (UF_long component = 0; component < solver->separator.component_count;
         ++component) {
      separator_component_last[component] = KLS_KLU_EMPTY;
    }
    for (UF_long pos = 0; pos < solver->separator.n; ++pos) {
      const unsigned int component = solver->separator.order_component[pos];
      if ((UF_long)component >= solver->separator.component_count) {
        goto fail;
      }
      separator_component_last[component] = pos;
    }
  }

  (void)kls_row_first_refactor_seed_init(&row_refactor_seed, n);
  if (pthread_mutex_init(&shared.commit_lock, NULL) != 0) {
    goto fail;
  }
  mutex_initialized = 1;
  shared.solver = solver;
  shared.numeric = numeric;
  shared.numeric_values = numeric_values;
  shared.psinv = psinv;
  shared.q_order = q_order;
  shared.row_refactor_seed = &row_refactor_seed;
  shared.separator_component_last = separator_component_last;
  shared.scaled = scaled;
  shared.scale = (int)common->scale;
  shared.use_separator_pivot_domains = use_separator_pivot_domains;
  UF_long dominant_pipeline_block = KLS_KLU_EMPTY;
  if (kls_first_factor_dominant_btf_pipeline_target(
        solver, &dominant_pipeline_block)) {
    const UF_long block_begin =
      solver->symbolic->R[dominant_pipeline_block];
    const UF_long block_end =
      solver->symbolic->R[dominant_pipeline_block + 1u];
    const UF_long block_rows =
      block_end >= block_begin ? block_end - block_begin : 0u;
    const int dominant_has_separator =
      use_separator_pivot_domains &&
      block_begin >= solver->separator.global_begin &&
      block_end <= solver->separator.global_end;
    solver->kls_first_last_dominant_btf_pipeline = 1;
    solver->kls_first_dominant_btf_pipeline_count++;
    solver->kls_first_last_dominant_btf_pipeline_block =
      dominant_pipeline_block;
    solver->kls_first_last_dominant_btf_pipeline_rows = block_rows;
    solver->kls_first_last_dominant_btf_pipeline_has_separator =
      dominant_has_separator;
    if (!dominant_has_separator) {
      solver->kls_first_dominant_btf_pipeline_without_separator_count++;
    }
    kls_trace_first_factor_dominant_btf(
      solver, dominant_pipeline_block, block_rows, dominant_has_separator);
    shared.force_block_pipeline = 1;
    shared.block_pipeline_thread_count = solver->options.threads;
    shared.block_pipeline_target_block = dominant_pipeline_block;
  }
  shared.tol = common->tol;
  shared.max_lnz_block = 1u;
  shared.max_unz_block = 1u;

  threads =
    (pthread_t *)calloc((size_t)thread_count, sizeof(*threads));
  if (threads == NULL) {
    goto fail;
  }
  int created = 0;
  for (; created < thread_count; ++created) {
    if (pthread_create(&threads[created], NULL,
                       kls_row_first_parallel_worker_main, &shared) != 0) {
      atomic_store_explicit(&shared.failed, 1, memory_order_release);
      break;
    }
  }
  for (int i = 0; i < created; ++i) {
    (void)pthread_join(threads[i], NULL);
  }
  if (created != thread_count ||
      atomic_load_explicit(&shared.failed, memory_order_acquire) ||
      shared.committed_blocks != solver->symbolic->nblocks) {
    goto fail;
  }

  numeric->lnz = shared.total_lnz;
  numeric->unz = shared.total_unz;
  numeric->max_lnz_block = shared.max_lnz_block;
  numeric->max_unz_block = shared.max_unz_block;
  if (!kls_rebuild_numeric_pinv(solver)) {
    goto fail;
  }
  memcpy(solver->symbolic->Q, q_order, (size_t)n * sizeof(*q_order));
  q_committed = 1;
  if (!kls_recompute_offdiag_from_pinv(solver, numeric_values)) {
    goto fail;
  }
  if (scaled && !kls_parallel_refactor_permute_scale(solver)) {
    goto fail;
  }
  if (row_refactor_seed.enabled) {
    (void)kls_row_first_refactor_seed_install(solver, &row_refactor_seed);
  }
  if (common->status == TRILINOS_KLU_OK) {
    common->numerical_rank = n;
    common->singular_col = n;
  }
  solver->kls_first_last_row_uplooking_columns = n;
  solver->kls_first_row_uplooking_column_count += n;
  if (shared.row_pipeline_runs > 0u) {
    solver->kls_first_last_row_pipeline = 1;
    solver->kls_first_row_pipeline_run_count += shared.row_pipeline_runs;
    solver->kls_first_last_row_pipeline_rows = shared.row_pipeline_rows;
    solver->kls_first_last_row_pipeline_threads =
      shared.row_pipeline_threads;
    if (shared.row_pipeline_pivot_tail_rows > 0u) {
      solver->kls_first_last_row_pipeline_pivot_tail = 1;
      solver->kls_first_row_pipeline_pivot_tail_run_count +=
        shared.row_pipeline_runs;
      solver->kls_first_last_row_pipeline_pivot_tail_rows =
        shared.row_pipeline_pivot_tail_rows;
      solver->kls_first_last_row_pipeline_pivot_restarts =
        shared.row_pipeline_pivot_restarts;
      solver->kls_first_row_pipeline_pivot_restart_count +=
        shared.row_pipeline_pivot_restarts;
    }
    if (shared.row_pipeline_separator_queue_runs > 0u) {
      solver->kls_first_last_separator_queue_executed = 1;
      solver->kls_first_separator_queue_executed_run_count +=
        shared.row_pipeline_separator_queue_runs;
      solver->kls_first_last_separator_queue_executed_private_rows =
        shared.row_pipeline_separator_private_rows;
      solver->kls_first_last_separator_queue_executed_pipeline_rows =
        shared.row_pipeline_separator_pipeline_rows;
      solver->kls_first_last_separator_queue_parallel_private = 1;
      solver->kls_first_separator_queue_parallel_private_run_count +=
        shared.row_pipeline_separator_queue_runs;
      solver->kls_first_last_separator_queue_parallel_private_rows =
        shared.row_pipeline_separator_private_rows;
      solver->kls_first_last_separator_queue_parallel_private_threads =
        shared.row_pipeline_separator_private_threads;
      solver->kls_first_last_separator_queue_parallel_pipeline = 1;
      solver->kls_first_separator_queue_parallel_pipeline_run_count +=
        shared.row_pipeline_separator_queue_runs;
      solver->kls_first_last_separator_queue_parallel_pipeline_rows =
        shared.row_pipeline_separator_pipeline_rows;
      solver->kls_first_last_separator_queue_parallel_pipeline_threads =
        shared.row_pipeline_threads;
    }
  }
  if (shared.row_supernode_update_groups > 0u) {
    solver->kls_first_last_row_supernode_update = 1;
    solver->kls_first_row_supernode_update_run_count++;
    kls_row_first_stats_add(&solver->kls_first_row_supernode_update_groups,
                            shared.row_supernode_update_groups);
    kls_row_first_stats_add(&solver->kls_first_row_supernode_update_rows,
                            shared.row_supernode_update_rows);
    solver->kls_first_last_row_supernode_update_groups =
      shared.row_supernode_update_groups;
    solver->kls_first_last_row_supernode_update_rows =
      shared.row_supernode_update_rows;
  }
  if (shared.row_supernode_panel_update_groups > 0u) {
    solver->kls_first_last_row_supernode_panel_update = 1;
    solver->kls_first_row_supernode_panel_update_run_count++;
    kls_row_first_stats_add(
      &solver->kls_first_row_supernode_panel_update_groups,
      shared.row_supernode_panel_update_groups);
    kls_row_first_stats_add(
      &solver->kls_first_row_supernode_panel_update_rows,
      shared.row_supernode_panel_update_rows);
    solver->kls_first_last_row_supernode_panel_update_groups =
      shared.row_supernode_panel_update_groups;
    solver->kls_first_last_row_supernode_panel_update_rows =
      shared.row_supernode_panel_update_rows;
  }
  kls_row_first_stats_add(&solver->kls_first_row_panel_cache_build_count,
                          shared.row_panel_cache_builds);
  kls_row_first_stats_add(&solver->kls_first_row_panel_cache_build_panels,
                          shared.row_panel_cache_build_panels);
  kls_row_first_stats_add(&solver->kls_first_row_panel_cache_build_entries,
                          shared.row_panel_cache_build_entries);
  kls_row_first_stats_add(&solver->kls_first_row_panel_cache_append_count,
                          shared.row_panel_cache_appends);
  kls_row_first_stats_add(&solver->kls_first_row_panel_cache_append_panels,
                          shared.row_panel_cache_append_panels);
  kls_row_first_stats_add(&solver->kls_first_row_panel_cache_append_entries,
                          shared.row_panel_cache_append_entries);
  kls_row_first_stats_add(&solver->kls_first_active_rank_pivot_reset_count,
                          shared.active_rank_pivot_resets);
  kls_row_first_stats_add(&solver->kls_first_active_rank_pivot_reset_rows,
                          shared.active_rank_pivot_reset_rows);
  kls_row_first_stats_add(
    &solver->kls_first_active_rank_pivot_panel_rebuild_count,
    shared.active_rank_pivot_panel_rebuilds);
  kls_row_first_stats_add(
    &solver->kls_first_active_rank_pivot_panel_rebuild_rows,
    shared.active_rank_pivot_panel_rebuild_rows);
  solver->kls_first_last_dynamic_column_pivots =
    shared.dynamic_column_pivots;
  solver->kls_first_dynamic_column_pivot_count +=
    shared.dynamic_column_pivots;
  solver->kls_first_last_separator_dynamic_column_pivots =
    shared.separator_dynamic_column_pivots;
  solver->kls_first_separator_dynamic_column_pivot_count +=
    shared.separator_dynamic_column_pivots;
  solver->kls_first_last_separator_extent_dynamic_column_pivots =
    shared.separator_extent_dynamic_column_pivots;
  solver->kls_first_separator_extent_dynamic_column_pivot_count +=
    shared.separator_extent_dynamic_column_pivots;
  solver->kls_first_last_separator_dynamic_column_fallbacks =
    shared.separator_dynamic_column_fallbacks;
  solver->kls_first_separator_dynamic_column_fallback_count +=
    shared.separator_dynamic_column_fallbacks;
  solver->kls_first_last_separator_dynamic_column_rejects =
    (UF_long)atomic_load_explicit(&shared.separator_dynamic_column_rejects,
                                  memory_order_relaxed);
  solver->kls_first_separator_dynamic_column_reject_count +=
    solver->kls_first_last_separator_dynamic_column_rejects;
  solver->kls_first_last_parallel_btf_blocks = shared.committed_blocks;
  solver->kls_first_parallel_btf_block_count += shared.committed_blocks;

  if (mutex_initialized) {
    (void)pthread_mutex_destroy(&shared.commit_lock);
  }
  free(threads);
  free(q_order);
  free(saved_q);
  free(separator_component_last);
  return 1;

fail:
  kls_row_first_refactor_seed_clear(&row_refactor_seed);
  if (q_committed && saved_q != NULL && solver->symbolic != NULL &&
      solver->symbolic->Q != NULL) {
    memcpy(solver->symbolic->Q, saved_q, (size_t)n * sizeof(*saved_q));
  }
  if (threads != NULL) {
    free(threads);
  }
  if (mutex_initialized) {
    (void)pthread_mutex_destroy(&shared.commit_lock);
  }
  free(q_order);
  free(saved_q);
  free(separator_component_last);
  trilinos_klu_l_free_numeric(&numeric, common);
  solver->numeric = NULL;
  *common = saved_common;
  return 0;
}

static int kls_try_first_factor_row_uplooking_blocks_impl(
  kls_solver *solver,
  double *numeric_values,
  int allow_separator_schedule) {
  if (solver == NULL || solver->symbolic == NULL || numeric_values == NULL ||
      solver->numeric != NULL || solver->symbolic->nblocks == 0u ||
      solver->symbolic->P == NULL || solver->symbolic->Q == NULL ||
      solver->symbolic->R == NULL) {
    return 0;
  }

  const UF_long n = solver->n;
  const UF_long maxblock = solver->symbolic->maxblock;
  if (n == 0u || maxblock == 0u) {
    return 0;
  }

  const int parallel_status =
    kls_try_first_factor_row_uplooking_blocks_parallel(solver,
                                                       numeric_values);
  if (parallel_status > 0) {
    return 1;
  }

  TRILINOS_KLU_common *common = &solver->common;
  trilinos_klu_l_common saved_common = *common;
  int q_committed = 0;
  UF_long separator_dynamic_column_rejects = 0;
  common->status = TRILINOS_KLU_OK;
  common->numerical_rank = KLS_KLU_EMPTY;
  common->singular_col = KLS_KLU_EMPTY;
  common->noffdiag = 0;
  common->nrealloc = 0;

  trilinos_klu_l_numeric *numeric = kls_allocate_numeric_skeleton(solver);
  if (numeric == NULL) {
    *common = saved_common;
    return 0;
  }
  solver->numeric = numeric;
  UF_long *psinv = numeric->Pinv;
  UF_long *row_counts = NULL;
  UF_long *row_ptr = NULL;
  double *x = NULL;
  unsigned int *mark = NULL;
  UF_long *pattern = NULL;
  UF_long *dep_heap = NULL;
  UF_long *u_row_ptr = NULL;
  UF_long *u_row_end = NULL;
  unsigned char *row_done = NULL;
  double *udiag_values = NULL;
  UF_long *q_order = NULL;
  UF_long *saved_q = NULL;
  UF_long *col_pos = NULL;
  UF_long *separator_component_last = NULL;
  kls_row_first_refactor_seed row_refactor_seed;
  memset(&row_refactor_seed, 0, sizeof(row_refactor_seed));
  const int scaled = common->scale > 0;
  if (scaled &&
      !trilinos_klu_l_scale((UF_long)common->scale, n,
                            solver->col_ptr, solver->row_idx,
                            numeric_values, numeric->Rs, NULL, common)) {
    goto fail;
  }

  row_counts = (UF_long *)calloc((size_t)maxblock, sizeof(*row_counts));
  row_ptr = (UF_long *)calloc((size_t)maxblock + 1u, sizeof(*row_ptr));
  x = (double *)calloc((size_t)maxblock, sizeof(*x));
  mark = (unsigned int *)calloc((size_t)maxblock, sizeof(*mark));
  pattern = (UF_long *)malloc((size_t)maxblock * sizeof(*pattern));
  dep_heap = (UF_long *)malloc((size_t)maxblock * sizeof(*dep_heap));
  u_row_ptr = (UF_long *)calloc((size_t)maxblock + 1u,
                                sizeof(*u_row_ptr));
  u_row_end = (UF_long *)calloc((size_t)maxblock, sizeof(*u_row_end));
  row_done = (unsigned char *)calloc((size_t)maxblock, sizeof(*row_done));
  udiag_values = (double *)malloc((size_t)maxblock *
                                  sizeof(*udiag_values));
  q_order = (UF_long *)malloc((size_t)n * sizeof(*q_order));
  saved_q = (UF_long *)malloc((size_t)n * sizeof(*saved_q));
  col_pos = (UF_long *)malloc((size_t)n * sizeof(*col_pos));
  if (row_counts == NULL || row_ptr == NULL || x == NULL ||
      mark == NULL || pattern == NULL || dep_heap == NULL ||
      u_row_ptr == NULL || u_row_end == NULL || row_done == NULL ||
      udiag_values == NULL || q_order == NULL || saved_q == NULL ||
      col_pos == NULL) {
    goto fail;
  }
  (void)kls_row_first_refactor_seed_init(&row_refactor_seed, n);

  for (UF_long k = 0; k < n; ++k) {
    psinv[k] = KLS_KLU_EMPTY;
  }
  for (UF_long k = 0; k < n; ++k) {
    const UF_long row = solver->symbolic->P[k];
    if (row >= n || psinv[row] != KLS_KLU_EMPTY) {
      goto fail;
    }
    psinv[row] = k;
  }
  for (UF_long k = 0; k < n; ++k) {
    if (psinv[k] == KLS_KLU_EMPTY) {
      goto fail;
    }
  }

  for (UF_long k = 0; k < n; ++k) {
    col_pos[k] = KLS_KLU_EMPTY;
  }
  for (UF_long k = 0; k < n; ++k) {
    const UF_long oldcol = solver->symbolic->Q[k];
    if (oldcol >= n || col_pos[oldcol] != KLS_KLU_EMPTY) {
      goto fail;
    }
    q_order[k] = oldcol;
    saved_q[k] = oldcol;
    col_pos[oldcol] = k;
  }
  for (UF_long k = 0; k < n; ++k) {
    col_pos[k] = KLS_KLU_EMPTY;
  }

  const double tol = common->tol;
  UF_long dynamic_column_pivots = 0;
  UF_long separator_dynamic_column_pivots = 0;
  UF_long separator_extent_dynamic_column_pivots = 0;
  UF_long separator_dynamic_column_fallbacks = 0;
  UF_long row_supernode_update_groups = 0;
  UF_long row_supernode_update_rows = 0;
  UF_long row_supernode_panel_update_groups = 0;
  UF_long row_supernode_panel_update_rows = 0;
  UF_long row_panel_cache_builds = 0;
  UF_long row_panel_cache_build_panels = 0;
  UF_long row_panel_cache_build_entries = 0;
  UF_long row_panel_cache_appends = 0;
  UF_long row_panel_cache_append_panels = 0;
  UF_long row_panel_cache_append_entries = 0;
  UF_long active_rank_pivot_resets = 0;
  UF_long active_rank_pivot_reset_rows = 0;
  UF_long active_rank_pivot_panel_rebuilds = 0;
  UF_long active_rank_pivot_panel_rebuild_rows = 0;
  const int use_separator_pivot_domains =
    kls_separator_analysis_has_global_range(&solver->separator) &&
    solver->separator.component_count > 0u &&
    solver->separator.order_component != NULL;
  if (use_separator_pivot_domains) {
    if (solver->separator.component_count >
        (UF_long)(SIZE_MAX / sizeof(*separator_component_last))) {
      goto fail;
    }
    separator_component_last =
      (UF_long *)malloc((size_t)solver->separator.component_count *
                        sizeof(*separator_component_last));
    if (separator_component_last == NULL) {
      goto fail;
    }
    for (UF_long component = 0; component < solver->separator.component_count;
         ++component) {
      separator_component_last[component] = KLS_KLU_EMPTY;
    }
    for (UF_long pos = 0; pos < solver->separator.n; ++pos) {
      const unsigned int component = solver->separator.order_component[pos];
      if ((UF_long)component >= solver->separator.component_count) {
        goto fail;
      }
      separator_component_last[component] = pos;
    }
  }
  UF_long total_lnz = 0;
  UF_long total_unz = 0;
  UF_long max_lnz_block = 1;
  UF_long max_unz_block = 1;
  int separator_queue_executed = 0;
  UF_long separator_queue_executed_private_rows = 0;
  UF_long separator_queue_executed_pipeline_rows = 0;
  int separator_queue_parallel_private = 0;
  UF_long separator_queue_parallel_private_rows = 0;
  UF_long separator_queue_parallel_private_threads = 0;
  int separator_queue_parallel_pipeline = 0;
  UF_long separator_queue_parallel_pipeline_rows = 0;
  UF_long separator_queue_parallel_pipeline_threads = 0;
  int row_pipeline_executed = 0;
  UF_long row_pipeline_rows = 0;
  UF_long row_pipeline_threads = 0;
  int row_pipeline_partial = 0;
  UF_long row_pipeline_partial_rows = 0;
  UF_long row_pipeline_partial_threads = 0;
  int row_pipeline_pivot_tail = 0;
  UF_long row_pipeline_pivot_tail_rows = 0;
  UF_long row_pipeline_pivot_restarts = 0;
  UF_long row_pipeline_pivot_serial_rows = 0;
  int row_pipeline_prefix_panel_rebuild = 0;
  UF_long row_pipeline_prefix_panel_rebuilds = 0;
  UF_long row_pipeline_prefix_panel_rebuild_rows = 0;
  int separator_queue_pipeline_partial = 0;
  UF_long separator_queue_pipeline_partial_rows = 0;
  UF_long separator_queue_pipeline_partial_threads = 0;
  int separator_queue_pipeline_wait_partial = 0;
  UF_long separator_queue_pipeline_wait_partial_rows = 0;
  UF_long separator_queue_pipeline_wait_partial_deps = 0;
  int separator_queue_pipeline_supernode_update = 0;
  UF_long separator_queue_pipeline_supernode_update_groups = 0;
  UF_long separator_queue_pipeline_supernode_update_rows = 0;
  int separator_queue_pipeline_supernode_panel_update = 0;
  UF_long separator_queue_pipeline_supernode_panel_update_groups = 0;
  UF_long separator_queue_pipeline_supernode_panel_update_rows = 0;
  int separator_queue_pipeline_pivot_tail = 0;
  UF_long separator_queue_pipeline_pivot_tail_rows = 0;
  UF_long separator_queue_pipeline_pivot_restarts = 0;
  UF_long separator_queue_pipeline_pivot_serial_rows = 0;
  int separator_queue_pipeline_prefix_panel_rebuild = 0;
  UF_long separator_queue_pipeline_prefix_panel_rebuilds = 0;
  UF_long separator_queue_pipeline_prefix_panel_rebuild_rows = 0;

  for (UF_long block = 0; block < solver->symbolic->nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    if (k2 <= k1 || k2 > n || k2 - k1 > maxblock) {
      goto fail;
    }
    const UF_long nk = k2 - k1;
    const int use_separator_for_block =
      use_separator_pivot_domains &&
      separator_component_last != NULL &&
      k1 >= solver->separator.global_begin &&
      k2 <= solver->separator.global_end;
    const UF_long separator_block_base = use_separator_for_block
      ? k1 - solver->separator.global_begin
      : 0u;
    memset(row_counts, 0, (size_t)nk * sizeof(*row_counts));
    memset(u_row_ptr, 0, ((size_t)nk + 1u) * sizeof(*u_row_ptr));
    memset(u_row_end, 0, (size_t)nk * sizeof(*u_row_end));
    memset(row_done, 0, (size_t)nk * sizeof(*row_done));

    for (UF_long local_col = 0; local_col < nk; ++local_col) {
      const UF_long oldcol = q_order[k1 + local_col];
      if (oldcol >= n || col_pos[oldcol] != KLS_KLU_EMPTY) {
        goto fail;
      }
      col_pos[oldcol] = local_col;
    }

    for (UF_long local_col = 0; local_col < nk; ++local_col) {
      const UF_long oldcol = q_order[k1 + local_col];
      for (UF_long p = solver->col_ptr[oldcol];
           p < solver->col_ptr[oldcol + 1u]; ++p) {
        const UF_long oldrow = solver->row_idx[p];
        if (oldrow >= n || psinv[oldrow] >= n) {
          goto fail;
        }
        const UF_long row_pos = psinv[oldrow];
        if (row_pos < k1) {
          continue;
        }
        if (row_pos >= k2) {
          goto fail;
        }
        row_counts[row_pos - k1]++;
      }
    }
    row_ptr[0] = 0;
    for (UF_long i = 0; i < nk; ++i) {
      row_ptr[i + 1u] = row_ptr[i] + row_counts[i];
    }
    const UF_long row_nnz = row_ptr[nk];
    UF_long *row_cols = row_nnz > 0u
      ? (UF_long *)malloc((size_t)row_nnz * sizeof(*row_cols)) : NULL;
    UF_long *row_input_pos = row_nnz > 0u
      ? (UF_long *)malloc((size_t)row_nnz * sizeof(*row_input_pos)) : NULL;
    UF_long *row_next =
      (UF_long *)malloc((size_t)nk * sizeof(*row_next));
    if ((row_nnz > 0u && (row_cols == NULL || row_input_pos == NULL)) ||
        row_next == NULL) {
      free(row_cols);
      free(row_input_pos);
      free(row_next);
      goto fail;
    }
    memcpy(row_next, row_ptr, (size_t)nk * sizeof(*row_next));
    for (UF_long local_col = 0; local_col < nk; ++local_col) {
      const UF_long oldcol = q_order[k1 + local_col];
      for (UF_long p = solver->col_ptr[oldcol];
           p < solver->col_ptr[oldcol + 1u]; ++p) {
        const UF_long row_pos = psinv[solver->row_idx[p]];
        if (row_pos < k1) {
          continue;
        }
        const UF_long row = row_pos - k1;
        const UF_long dst = row_next[row]++;
        row_cols[dst] = oldcol;
        row_input_pos[dst] = p;
      }
    }
    free(row_next);
    row_next = NULL;

    kls_row_first_entries l_entries;
    kls_row_first_entries u_entries;
    memset(&l_entries, 0, sizeof(l_entries));
    memset(&u_entries, 0, sizeof(u_entries));
    UF_long *row_order = NULL;
    UF_long *local_row_order = NULL;
    UF_long *old_to_new_row = NULL;
    unsigned char *row_seen = NULL;
    int *row_owner = NULL;
    int use_separator_row_order = 0;
    int parallel_private_executed = 0;
    int parallel_pipeline_executed = 0;
    UF_long block_parallel_private_rows = 0;
    UF_long block_parallel_private_threads = 0;
    UF_long block_parallel_pipeline_rows = 0;
    UF_long block_parallel_pipeline_threads = 0;
    UF_long block_pipeline_partial_rows = 0;
    UF_long block_pipeline_wait_partial_rows = 0;
    UF_long block_pipeline_wait_partial_deps = 0;
    UF_long block_pipeline_supernode_update_groups = 0;
    UF_long block_pipeline_supernode_update_rows = 0;
    UF_long block_pipeline_supernode_panel_update_groups = 0;
    UF_long block_pipeline_supernode_panel_update_rows = 0;
    UF_long block_pipeline_pivot_tail_rows = 0;
    UF_long block_pipeline_pivot_restarts = 0;
    UF_long block_pipeline_pivot_serial_rows = 0;
    UF_long block_separator_private_rows = 0;
    UF_long block_separator_pipeline_rows = 0;
    kls_first_separator_queue_plan separator_plan;
    memset(&separator_plan, 0, sizeof(separator_plan));
    kls_row_first_supernode_panel_cache row_supernode_panel_cache;
    memset(&row_supernode_panel_cache, 0, sizeof(row_supernode_panel_cache));
    UF_long row_open_panel_start = KLS_KLU_EMPTY;
    UF_long row_open_panel_end = KLS_KLU_EMPTY;
    if (!kls_row_first_entries_enable_col_counts(&l_entries, nk) ||
        !kls_row_first_entries_enable_col_counts(&u_entries, nk) ||
        !kls_row_first_entries_enable_col_links(&u_entries, nk)) {
      goto fail_block_entries;
    }
    const UF_long reserve_entries =
      kls_row_first_symbolic_reserve(solver, block, nk);
    if (reserve_entries > 0u) {
      (void)kls_row_first_entries_reserve(&l_entries, reserve_entries);
      (void)kls_row_first_entries_reserve(&u_entries, reserve_entries);
    }

    if (allow_separator_schedule && use_separator_for_block &&
        solver->options.threads > 1) {
      int have_separator_plan = kls_build_first_separator_queue_plan(
        solver, k1, k2, solver->options.threads, &separator_plan);
      if (have_separator_plan && separator_plan.partitioned &&
          !kls_first_separator_private_ownership_valid(
            k1, nk, n, row_ptr, row_cols, col_pos, &separator_plan,
            solver->options.threads)) {
        kls_first_separator_queue_plan_free(&separator_plan);
        have_separator_plan =
          kls_build_first_separator_queue_plan_legacy_from_counts(
            solver, k1, k2, solver->options.threads, row_counts,
            &separator_plan);
      }
      if (have_separator_plan) {
        row_order = (UF_long *)malloc((size_t)nk * sizeof(*row_order));
        row_seen = (unsigned char *)calloc((size_t)nk, sizeof(*row_seen));
        if (row_order == NULL || row_seen == NULL) {
          goto fail_block_entries;
        }
        UF_long row_order_count = 0;
        for (UF_long pos = 0; pos < separator_plan.private_rows_count; ++pos) {
          const UF_long global_row = separator_plan.private_rows[pos];
          if (global_row < k1 || global_row >= k2) {
            goto fail_block_entries;
          }
          const UF_long local_row = global_row - k1;
          if (row_seen[local_row]) {
            goto fail_block_entries;
          }
          row_seen[local_row] = 1u;
          row_order[row_order_count++] = local_row;
        }
        block_separator_private_rows = row_order_count;
        qsort(row_order, (size_t)block_separator_private_rows,
              sizeof(*row_order), kls_compare_uf_long);
        for (UF_long pos = 0; pos < separator_plan.pipeline_rows_count; ++pos) {
          const UF_long global_row = separator_plan.pipeline_rows[pos];
          if (global_row < k1 || global_row >= k2) {
            goto fail_block_entries;
          }
          const UF_long local_row = global_row - k1;
          if (row_seen[local_row]) {
            goto fail_block_entries;
          }
          row_seen[local_row] = 1u;
          row_order[row_order_count++] = local_row;
        }
        block_separator_pipeline_rows =
          row_order_count - block_separator_private_rows;
        qsort(row_order + block_separator_private_rows,
              (size_t)block_separator_pipeline_rows, sizeof(*row_order),
              kls_compare_uf_long);
        if (row_order_count == nk && block_separator_private_rows > 0u) {
          use_separator_row_order = 1;
        }
      }
    }
    if (use_separator_row_order) {
      local_row_order = row_order;
      row_order = (UF_long *)malloc((size_t)nk * sizeof(*row_order));
      old_to_new_row =
        (UF_long *)malloc((size_t)nk * sizeof(*old_to_new_row));
      if (row_order == NULL || old_to_new_row == NULL) {
        goto fail_block_entries;
      }
      for (UF_long new_row = 0; new_row < nk; ++new_row) {
        const UF_long old_row = local_row_order[new_row];
        if (old_row >= nk) {
          goto fail_block_entries;
        }
        row_order[new_row] = new_row;
        old_to_new_row[old_row] = new_row;
      }
      for (UF_long pos = 0; pos < separator_plan.private_rows_count; ++pos) {
        const UF_long old_global = separator_plan.private_rows[pos];
        if (old_global < k1 || old_global >= k2) {
          goto fail_block_entries;
        }
        separator_plan.private_rows[pos] =
          k1 + old_to_new_row[old_global - k1];
      }
      for (UF_long pos = 0; pos < separator_plan.pipeline_rows_count; ++pos) {
        const UF_long old_global = separator_plan.pipeline_rows[pos];
        if (old_global < k1 || old_global >= k2) {
          goto fail_block_entries;
        }
        separator_plan.pipeline_rows[pos] =
          k1 + old_to_new_row[old_global - k1];
      }

      UF_long *old_block_q =
        (UF_long *)malloc((size_t)nk * sizeof(*old_block_q));
      if (old_block_q == NULL) {
        goto fail_block_entries;
      }
      for (UF_long local_col = 0; local_col < nk; ++local_col) {
        old_block_q[local_col] = q_order[k1 + local_col];
        if (old_block_q[local_col] < n) {
          col_pos[old_block_q[local_col]] = KLS_KLU_EMPTY;
        }
      }
      for (UF_long local_col = 0; local_col < nk; ++local_col) {
        const UF_long old_local_col = local_row_order[local_col];
        if (old_local_col >= nk) {
          free(old_block_q);
          goto fail_block_entries;
        }
        q_order[k1 + local_col] = old_block_q[old_local_col];
      }
      free(old_block_q);

      memset(row_counts, 0, (size_t)nk * sizeof(*row_counts));
      memset(row_ptr, 0, ((size_t)nk + 1u) * sizeof(*row_ptr));
      for (UF_long local_col = 0; local_col < nk; ++local_col) {
        const UF_long oldcol = q_order[k1 + local_col];
        if (oldcol >= n || col_pos[oldcol] != KLS_KLU_EMPTY) {
          goto fail_block_entries;
        }
        col_pos[oldcol] = local_col;
      }

      for (UF_long local_col = 0; local_col < nk; ++local_col) {
        const UF_long oldcol = q_order[k1 + local_col];
        for (UF_long p = solver->col_ptr[oldcol];
             p < solver->col_ptr[oldcol + 1u]; ++p) {
          const UF_long oldrow = solver->row_idx[p];
          if (oldrow >= n || psinv[oldrow] >= n) {
            goto fail_block_entries;
          }
          const UF_long row_pos = psinv[oldrow];
          if (row_pos < k1) {
            continue;
          }
          if (row_pos >= k2) {
            goto fail_block_entries;
          }
          const UF_long old_local_row = row_pos - k1;
          if (old_local_row >= nk) {
            goto fail_block_entries;
          }
          row_counts[old_to_new_row[old_local_row]]++;
        }
      }
      row_ptr[0] = 0;
      for (UF_long i = 0; i < nk; ++i) {
        row_ptr[i + 1u] = row_ptr[i] + row_counts[i];
      }
      const UF_long rebuilt_row_nnz = row_ptr[nk];
      UF_long *rebuilt_row_cols = rebuilt_row_nnz > 0u
        ? (UF_long *)malloc((size_t)rebuilt_row_nnz *
                            sizeof(*rebuilt_row_cols)) : NULL;
      UF_long *rebuilt_row_input_pos = rebuilt_row_nnz > 0u
        ? (UF_long *)malloc((size_t)rebuilt_row_nnz *
                            sizeof(*rebuilt_row_input_pos)) : NULL;
      row_next = (UF_long *)malloc((size_t)nk * sizeof(*row_next));
      if ((rebuilt_row_nnz > 0u &&
           (rebuilt_row_cols == NULL || rebuilt_row_input_pos == NULL)) ||
          row_next == NULL) {
        free(rebuilt_row_cols);
        free(rebuilt_row_input_pos);
        goto fail_block_entries;
      }
      memcpy(row_next, row_ptr, (size_t)nk * sizeof(*row_next));
      for (UF_long local_col = 0; local_col < nk; ++local_col) {
        const UF_long oldcol = q_order[k1 + local_col];
        for (UF_long p = solver->col_ptr[oldcol];
             p < solver->col_ptr[oldcol + 1u]; ++p) {
          const UF_long oldrow = solver->row_idx[p];
          if (oldrow >= n || psinv[oldrow] >= n) {
            free(rebuilt_row_cols);
            free(rebuilt_row_input_pos);
            goto fail_block_entries;
          }
          const UF_long row_pos = psinv[oldrow];
          if (row_pos < k1) {
            continue;
          }
          if (row_pos >= k2) {
            free(rebuilt_row_cols);
            free(rebuilt_row_input_pos);
            goto fail_block_entries;
          }
          const UF_long old_local_row = row_pos - k1;
          if (old_local_row >= nk) {
            free(rebuilt_row_cols);
            free(rebuilt_row_input_pos);
            goto fail_block_entries;
          }
          const UF_long row = old_to_new_row[old_local_row];
          const UF_long dst = row_next[row]++;
          rebuilt_row_cols[dst] = oldcol;
          rebuilt_row_input_pos[dst] = p;
        }
      }
      free(row_next);
      row_next = NULL;
      free(row_cols);
      free(row_input_pos);
      row_cols = rebuilt_row_cols;
      row_input_pos = rebuilt_row_input_pos;
    }

    kls_row_first_block_context row_ctx;
    memset(&row_ctx, 0, sizeof(row_ctx));
    row_ctx.solver = solver;
    row_ctx.numeric = numeric;
    row_ctx.numeric_values = numeric_values;
    row_ctx.row_ptr = row_ptr;
    row_ctx.row_cols = row_cols;
    row_ctx.row_input_pos = row_input_pos;
    row_ctx.q_order = q_order + k1;
    row_ctx.col_pos = col_pos;
    row_ctx.separator_component_last = separator_component_last;
    row_ctx.separator_local_order = local_row_order;
    row_ctx.k1 = k1;
    row_ctx.nk = nk;
    row_ctx.n = n;
    row_ctx.separator_block_base = separator_block_base;
    row_ctx.scaled = scaled;
    row_ctx.scale = (int)common->scale;
    row_ctx.use_separator_for_block = use_separator_for_block;
    row_ctx.strict_separator_pivot_scope =
      kls_row_first_strict_separator_pivot_scope_env_enabled();
    row_ctx.tol = tol;

    kls_row_first_workspace row_workspace;
    memset(&row_workspace, 0, sizeof(row_workspace));
    row_workspace.x = x;
    row_workspace.mark = mark;
    row_workspace.pattern = pattern;
    row_workspace.dep_heap = dep_heap;
    row_workspace.u_row_ptr = u_row_ptr;
    row_workspace.u_row_end = u_row_end;
    if (kls_row_first_supernode_panel_cache_init_rows(
          &row_supernode_panel_cache, nk)) {
      row_workspace.supernode_panel_cache = &row_supernode_panel_cache;
    }

    kls_row_first_row_stats row_stats;
    memset(&row_stats, 0, sizeof(row_stats));
    if (use_separator_row_order && separator_plan.nonempty_threads > 1u) {
      row_owner = (int *)malloc((size_t)nk * sizeof(*row_owner));
      if (row_owner == NULL ||
          !kls_row_first_build_private_row_owners(
            &row_ctx, &separator_plan, solver->options.threads,
            row_owner)) {
        goto fail_block_entries;
      }
      if (kls_row_first_run_parallel_private_phase(
            &row_ctx, &separator_plan, solver->options.threads, row_owner,
            &l_entries, &u_entries, udiag_values, u_row_ptr, u_row_end,
            row_done, &row_stats, &block_parallel_private_rows,
            &block_parallel_private_threads, 0, NULL)) {
        parallel_private_executed = 1;
      } else {
        goto fail_block_entries;
      }
    }
    if (use_separator_row_order) {
      if (!parallel_private_executed) {
        for (UF_long iter = 0; iter < block_separator_private_rows; ++iter) {
          const UF_long i = row_order[iter];
          const UF_long pivots_before = row_stats.dynamic_column_pivots;
          if (!kls_row_first_factor_one_row(
                &row_ctx, &row_workspace, &l_entries, &u_entries,
                udiag_values, i, NULL, -1, row_done, &row_stats, NULL)) {
            goto fail_block_entries;
          }
          if (row_stats.dynamic_column_pivots != pivots_before) {
            kls_row_first_supernode_panel_cache_reset_active(
              &row_supernode_panel_cache, &row_workspace, nk,
              &row_open_panel_start, &row_open_panel_end);
          }
          kls_row_first_supernode_panel_cache_publish_completed(
            row_workspace.supernode_panel_cache != NULL
              ? &row_supernode_panel_cache : NULL,
            &u_entries, u_row_ptr, u_row_end, udiag_values, nk,
            &row_stats, &row_open_panel_start, &row_open_panel_end, i);
        }
      }
      if (block_separator_pipeline_rows > 0u) {
        kls_row_first_supernode_panel_cache_reset_active(
          &row_supernode_panel_cache, &row_workspace, nk,
          &row_open_panel_start, &row_open_panel_end);
        if (kls_row_first_run_restartable_pipeline_suffix(
              &row_ctx, row_order, NULL, block_separator_private_rows, nk,
              solver->options.threads, &row_workspace, &l_entries,
              &u_entries, udiag_values, row_done, &row_stats,
              &block_parallel_pipeline_rows,
              &block_parallel_pipeline_threads,
              &block_pipeline_partial_rows,
              &block_pipeline_wait_partial_rows,
              &block_pipeline_wait_partial_deps,
              &block_pipeline_supernode_update_groups,
              &block_pipeline_supernode_update_rows,
              &block_pipeline_supernode_panel_update_groups,
              &block_pipeline_supernode_panel_update_rows,
              &block_pipeline_pivot_tail_rows,
              &block_pipeline_pivot_restarts,
              &block_pipeline_pivot_serial_rows,
              row_workspace.supernode_panel_cache != NULL
                ? &row_supernode_panel_cache : NULL,
              &row_open_panel_start, &row_open_panel_end)) {
          parallel_pipeline_executed = 1;
        } else {
          goto fail_block_entries;
        }
      }
    } else if (solver->options.threads > 1 && nk > 1u) {
      free(row_order);
      row_order = NULL;
      free(row_seen);
      row_seen = NULL;
      row_order = (UF_long *)malloc((size_t)nk * sizeof(*row_order));
      if (row_order == NULL) {
        goto fail_block_entries;
      }
      for (UF_long iter = 0; iter < nk; ++iter) {
        row_order[iter] = iter;
      }
      kls_row_first_supernode_panel_cache_reset_active(
        &row_supernode_panel_cache, &row_workspace, nk,
        &row_open_panel_start, &row_open_panel_end);
      if (kls_row_first_run_restartable_pipeline_suffix(
            &row_ctx, row_order, NULL, 0, nk, solver->options.threads,
            &row_workspace, &l_entries, &u_entries, udiag_values, row_done,
            &row_stats, &block_parallel_pipeline_rows,
            &block_parallel_pipeline_threads,
            &block_pipeline_partial_rows,
            &block_pipeline_wait_partial_rows,
            &block_pipeline_wait_partial_deps,
            &block_pipeline_supernode_update_groups,
            &block_pipeline_supernode_update_rows,
            &block_pipeline_supernode_panel_update_groups,
            &block_pipeline_supernode_panel_update_rows,
            &block_pipeline_pivot_tail_rows,
            &block_pipeline_pivot_restarts,
            &block_pipeline_pivot_serial_rows,
            row_workspace.supernode_panel_cache != NULL
              ? &row_supernode_panel_cache : NULL,
            &row_open_panel_start, &row_open_panel_end)) {
        parallel_pipeline_executed = 1;
      } else {
        goto fail_block_entries;
      }
    } else {
      for (UF_long iter = 0; iter < nk; ++iter) {
        const UF_long pivots_before = row_stats.dynamic_column_pivots;
        if (!kls_row_first_factor_one_row(
              &row_ctx, &row_workspace, &l_entries, &u_entries,
              udiag_values, iter, NULL, -1, NULL, &row_stats, NULL)) {
          goto fail_block_entries;
        }
        if (row_stats.dynamic_column_pivots != pivots_before) {
          kls_row_first_supernode_panel_cache_reset_active(
            &row_supernode_panel_cache, &row_workspace, nk,
            &row_open_panel_start, &row_open_panel_end);
        }
        kls_row_first_supernode_panel_cache_publish_completed(
          row_workspace.supernode_panel_cache != NULL
            ? &row_supernode_panel_cache : NULL,
          &u_entries, u_row_ptr, u_row_end, udiag_values, nk,
          &row_stats, &row_open_panel_start, &row_open_panel_end, iter);
      }
    }
    dynamic_column_pivots += row_stats.dynamic_column_pivots;
    separator_dynamic_column_pivots +=
      row_stats.separator_dynamic_column_pivots;
    separator_extent_dynamic_column_pivots +=
      row_stats.separator_extent_dynamic_column_pivots;
    separator_dynamic_column_fallbacks +=
      row_stats.separator_dynamic_column_fallbacks;
    separator_dynamic_column_rejects +=
      row_stats.separator_dynamic_column_rejects;
    row_supernode_update_groups += row_stats.supernode_update_groups;
    row_supernode_update_rows += row_stats.supernode_update_rows;
    row_supernode_panel_update_groups +=
      row_stats.supernode_panel_update_groups;
    row_supernode_panel_update_rows +=
      row_stats.supernode_panel_update_rows;
    if (use_separator_row_order) {
      separator_queue_executed = 1;
      separator_queue_executed_private_rows += block_separator_private_rows;
      separator_queue_executed_pipeline_rows += block_separator_pipeline_rows;
    }
    if (parallel_private_executed) {
      separator_queue_parallel_private = 1;
      separator_queue_parallel_private_rows += block_parallel_private_rows;
      if (block_parallel_private_threads >
          separator_queue_parallel_private_threads) {
        separator_queue_parallel_private_threads =
          block_parallel_private_threads;
      }
    }
    if (parallel_pipeline_executed && use_separator_row_order) {
      separator_queue_parallel_pipeline = 1;
      separator_queue_parallel_pipeline_rows += block_parallel_pipeline_rows;
      if (block_parallel_pipeline_threads >
          separator_queue_parallel_pipeline_threads) {
        separator_queue_parallel_pipeline_threads =
          block_parallel_pipeline_threads;
      }
      if (block_pipeline_partial_rows > 0u) {
        separator_queue_pipeline_partial = 1;
        separator_queue_pipeline_partial_rows += block_pipeline_partial_rows;
        if (block_parallel_pipeline_threads >
            separator_queue_pipeline_partial_threads) {
          separator_queue_pipeline_partial_threads =
            block_parallel_pipeline_threads;
        }
      }
      if (block_pipeline_wait_partial_rows > 0u) {
        separator_queue_pipeline_wait_partial = 1;
        separator_queue_pipeline_wait_partial_rows +=
          block_pipeline_wait_partial_rows;
        separator_queue_pipeline_wait_partial_deps +=
          block_pipeline_wait_partial_deps;
      }
      if (block_pipeline_supernode_update_groups > 0u) {
        separator_queue_pipeline_supernode_update = 1;
        separator_queue_pipeline_supernode_update_groups +=
          block_pipeline_supernode_update_groups;
        separator_queue_pipeline_supernode_update_rows +=
          block_pipeline_supernode_update_rows;
      }
      if (block_pipeline_supernode_panel_update_groups > 0u) {
        separator_queue_pipeline_supernode_panel_update = 1;
        separator_queue_pipeline_supernode_panel_update_groups +=
          block_pipeline_supernode_panel_update_groups;
        separator_queue_pipeline_supernode_panel_update_rows +=
          block_pipeline_supernode_panel_update_rows;
      }
      if (block_pipeline_pivot_tail_rows > 0u) {
        separator_queue_pipeline_pivot_tail = 1;
        separator_queue_pipeline_pivot_tail_rows +=
          block_pipeline_pivot_tail_rows;
        separator_queue_pipeline_pivot_restarts +=
          block_pipeline_pivot_restarts;
        separator_queue_pipeline_pivot_serial_rows +=
          block_pipeline_pivot_serial_rows;
      }
      if (row_stats.pipeline_prefix_panel_rebuilds > 0u) {
        separator_queue_pipeline_prefix_panel_rebuild = 1;
        separator_queue_pipeline_prefix_panel_rebuilds +=
          row_stats.pipeline_prefix_panel_rebuilds;
        separator_queue_pipeline_prefix_panel_rebuild_rows +=
          row_stats.pipeline_prefix_panel_rebuild_rows;
      }
    } else if (parallel_pipeline_executed) {
      row_pipeline_executed = 1;
      row_pipeline_rows += block_parallel_pipeline_rows;
      if (block_parallel_pipeline_threads > row_pipeline_threads) {
        row_pipeline_threads = block_parallel_pipeline_threads;
      }
      if (block_pipeline_partial_rows > 0u) {
        row_pipeline_partial = 1;
        row_pipeline_partial_rows += block_pipeline_partial_rows;
        if (block_parallel_pipeline_threads > row_pipeline_partial_threads) {
          row_pipeline_partial_threads = block_parallel_pipeline_threads;
        }
      }
      if (block_pipeline_pivot_tail_rows > 0u) {
        row_pipeline_pivot_tail = 1;
        row_pipeline_pivot_tail_rows += block_pipeline_pivot_tail_rows;
        row_pipeline_pivot_restarts += block_pipeline_pivot_restarts;
        row_pipeline_pivot_serial_rows += block_pipeline_pivot_serial_rows;
      }
      if (row_stats.pipeline_prefix_panel_rebuilds > 0u) {
        row_pipeline_prefix_panel_rebuild = 1;
        row_pipeline_prefix_panel_rebuilds +=
          row_stats.pipeline_prefix_panel_rebuilds;
        row_pipeline_prefix_panel_rebuild_rows +=
          row_stats.pipeline_prefix_panel_rebuild_rows;
      }
    }
    kls_row_first_stats_add(&row_panel_cache_builds,
                            row_stats.panel_cache_builds);
    kls_row_first_stats_add(&row_panel_cache_build_panels,
                            row_stats.panel_cache_build_panels);
    kls_row_first_stats_add(&row_panel_cache_build_entries,
                            row_stats.panel_cache_build_entries);
    kls_row_first_stats_add(&row_panel_cache_appends,
                            row_stats.panel_cache_appends);
    kls_row_first_stats_add(&row_panel_cache_append_panels,
                            row_stats.panel_cache_append_panels);
    kls_row_first_stats_add(&row_panel_cache_append_entries,
                            row_stats.panel_cache_append_entries);
    kls_row_first_stats_add(&active_rank_pivot_resets,
                            row_stats.active_rank_pivot_resets);
    kls_row_first_stats_add(&active_rank_pivot_reset_rows,
                            row_stats.active_rank_pivot_reset_rows);
    kls_row_first_stats_add(&active_rank_pivot_panel_rebuilds,
                            row_stats.active_rank_pivot_panel_rebuilds);
    kls_row_first_stats_add(&active_rank_pivot_panel_rebuild_rows,
                            row_stats.active_rank_pivot_panel_rebuild_rows);
    if (!kls_pack_row_first_block_numeric(
          solver, numeric, block, k1, nk, local_row_order,
          &l_entries, &u_entries,
          udiag_values, &total_lnz, &total_unz, &max_lnz_block,
          &max_unz_block)) {
      goto fail_block_entries;
    }
    if (row_refactor_seed.enabled &&
        !kls_row_first_refactor_seed_capture_block(
          &row_refactor_seed, solver, numeric, block, k1, nk,
          row_ptr, row_cols, row_input_pos, col_pos,
          &l_entries, &u_entries)) {
      kls_row_first_refactor_seed_clear(&row_refactor_seed);
    }
    /* The block-local workspace aliases the outer row arrays, but its
       supernode scratch is grown locally by the row/pipeline consumers. */
    free(row_workspace.supernode_workspace);
    row_workspace.supernode_workspace = NULL;
    row_workspace.supernode_workspace_capacity = 0u;
    kls_row_first_supernode_panel_cache_free(&row_supernode_panel_cache);
    kls_row_first_entries_free(&l_entries);
    kls_row_first_entries_free(&u_entries);
    kls_first_separator_queue_plan_free(&separator_plan);
    free(row_order);
    free(local_row_order);
    free(old_to_new_row);
    free(row_seen);
    free(row_owner);
    free(row_cols);
    free(row_input_pos);
    for (UF_long local_col = 0; local_col < nk; ++local_col) {
      const UF_long oldcol = q_order[k1 + local_col];
      if (oldcol < n) {
        col_pos[oldcol] = KLS_KLU_EMPTY;
      }
    }
    continue;

fail_block_entries:
    free(row_workspace.supernode_workspace);
    row_workspace.supernode_workspace = NULL;
    row_workspace.supernode_workspace_capacity = 0u;
    kls_row_first_supernode_panel_cache_free(&row_supernode_panel_cache);
    kls_row_first_entries_free(&l_entries);
    kls_row_first_entries_free(&u_entries);
    kls_first_separator_queue_plan_free(&separator_plan);
    free(row_order);
    free(local_row_order);
    free(old_to_new_row);
    free(row_seen);
    free(row_owner);
    free(row_cols);
    free(row_input_pos);
    free(row_next);
    goto fail;
  }

  numeric->lnz = total_lnz;
  numeric->unz = total_unz;
  numeric->max_lnz_block = max_lnz_block;
  numeric->max_unz_block = max_unz_block;
  if (!kls_rebuild_numeric_pinv(solver)) {
    goto fail;
  }
  memcpy(solver->symbolic->Q, q_order, (size_t)n * sizeof(*q_order));
  q_committed = 1;
  if (!kls_recompute_offdiag_from_pinv(solver, numeric_values)) {
    goto fail;
  }
  if (scaled && !kls_parallel_refactor_permute_scale(solver)) {
    goto fail;
  }
  if (row_refactor_seed.enabled) {
    (void)kls_row_first_refactor_seed_install(solver, &row_refactor_seed);
  }
  if (common->status == TRILINOS_KLU_OK) {
    common->numerical_rank = n;
    common->singular_col = n;
  }
  solver->kls_first_last_row_uplooking_columns = n;
  solver->kls_first_row_uplooking_column_count += n;
  if (row_pipeline_executed) {
    solver->kls_first_last_row_pipeline = 1;
    solver->kls_first_row_pipeline_run_count++;
    solver->kls_first_last_row_pipeline_rows = row_pipeline_rows;
    solver->kls_first_last_row_pipeline_threads = row_pipeline_threads;
  }
  if (row_pipeline_partial) {
    solver->kls_first_last_row_pipeline_partial = 1;
    solver->kls_first_row_pipeline_partial_run_count++;
    solver->kls_first_last_row_pipeline_partial_rows =
      row_pipeline_partial_rows;
    solver->kls_first_last_row_pipeline_partial_threads =
      row_pipeline_partial_threads;
  }
  if (row_pipeline_pivot_tail) {
    solver->kls_first_last_row_pipeline_pivot_tail = 1;
    solver->kls_first_row_pipeline_pivot_tail_run_count++;
    solver->kls_first_last_row_pipeline_pivot_tail_rows =
      row_pipeline_pivot_tail_rows;
    solver->kls_first_last_row_pipeline_pivot_restarts =
      row_pipeline_pivot_restarts;
    solver->kls_first_row_pipeline_pivot_restart_count +=
      row_pipeline_pivot_restarts;
    solver->kls_first_last_row_pipeline_pivot_serial_rows =
      row_pipeline_pivot_serial_rows;
  }
  if (row_pipeline_prefix_panel_rebuild) {
    solver->kls_first_last_row_pipeline_prefix_panel_rebuild = 1;
    solver->kls_first_row_pipeline_prefix_panel_rebuild_count +=
      row_pipeline_prefix_panel_rebuilds;
    solver->kls_first_last_row_pipeline_prefix_panel_rebuild_rows =
      row_pipeline_prefix_panel_rebuild_rows;
  }
  kls_row_first_stats_add(&solver->kls_first_row_panel_cache_build_count,
                          row_panel_cache_builds);
  kls_row_first_stats_add(&solver->kls_first_row_panel_cache_build_panels,
                          row_panel_cache_build_panels);
  kls_row_first_stats_add(&solver->kls_first_row_panel_cache_build_entries,
                          row_panel_cache_build_entries);
  kls_row_first_stats_add(&solver->kls_first_row_panel_cache_append_count,
                          row_panel_cache_appends);
  kls_row_first_stats_add(&solver->kls_first_row_panel_cache_append_panels,
                          row_panel_cache_append_panels);
  kls_row_first_stats_add(&solver->kls_first_row_panel_cache_append_entries,
                          row_panel_cache_append_entries);
  kls_row_first_stats_add(&solver->kls_first_active_rank_pivot_reset_count,
                          active_rank_pivot_resets);
  kls_row_first_stats_add(&solver->kls_first_active_rank_pivot_reset_rows,
                          active_rank_pivot_reset_rows);
  kls_row_first_stats_add(
    &solver->kls_first_active_rank_pivot_panel_rebuild_count,
    active_rank_pivot_panel_rebuilds);
  kls_row_first_stats_add(
    &solver->kls_first_active_rank_pivot_panel_rebuild_rows,
    active_rank_pivot_panel_rebuild_rows);
  if (row_supernode_update_groups > 0u) {
    solver->kls_first_last_row_supernode_update = 1;
    solver->kls_first_row_supernode_update_run_count++;
    kls_row_first_stats_add(&solver->kls_first_row_supernode_update_groups,
                            row_supernode_update_groups);
    kls_row_first_stats_add(&solver->kls_first_row_supernode_update_rows,
                            row_supernode_update_rows);
    solver->kls_first_last_row_supernode_update_groups =
      row_supernode_update_groups;
    solver->kls_first_last_row_supernode_update_rows =
      row_supernode_update_rows;
  }
  if (row_supernode_panel_update_groups > 0u) {
    solver->kls_first_last_row_supernode_panel_update = 1;
    solver->kls_first_row_supernode_panel_update_run_count++;
    kls_row_first_stats_add(
      &solver->kls_first_row_supernode_panel_update_groups,
      row_supernode_panel_update_groups);
    kls_row_first_stats_add(&solver->kls_first_row_supernode_panel_update_rows,
                            row_supernode_panel_update_rows);
    solver->kls_first_last_row_supernode_panel_update_groups =
      row_supernode_panel_update_groups;
    solver->kls_first_last_row_supernode_panel_update_rows =
      row_supernode_panel_update_rows;
  }
  solver->kls_first_last_dynamic_column_pivots = dynamic_column_pivots;
  solver->kls_first_dynamic_column_pivot_count += dynamic_column_pivots;
  solver->kls_first_last_separator_dynamic_column_pivots =
    separator_dynamic_column_pivots;
  solver->kls_first_separator_dynamic_column_pivot_count +=
    separator_dynamic_column_pivots;
  solver->kls_first_last_separator_extent_dynamic_column_pivots =
    separator_extent_dynamic_column_pivots;
  solver->kls_first_separator_extent_dynamic_column_pivot_count +=
    separator_extent_dynamic_column_pivots;
  solver->kls_first_last_separator_dynamic_column_fallbacks =
    separator_dynamic_column_fallbacks;
  solver->kls_first_separator_dynamic_column_fallback_count +=
    separator_dynamic_column_fallbacks;
  solver->kls_first_last_separator_dynamic_column_rejects =
    separator_dynamic_column_rejects;
  solver->kls_first_separator_dynamic_column_reject_count +=
    separator_dynamic_column_rejects;
  if (separator_queue_executed) {
    solver->kls_first_last_separator_queue_executed = 1;
    solver->kls_first_separator_queue_executed_run_count++;
    solver->kls_first_last_separator_queue_executed_private_rows =
      separator_queue_executed_private_rows;
    solver->kls_first_last_separator_queue_executed_pipeline_rows =
      separator_queue_executed_pipeline_rows;
  }
  if (separator_queue_parallel_private) {
    solver->kls_first_last_separator_queue_parallel_private = 1;
    solver->kls_first_separator_queue_parallel_private_run_count++;
    solver->kls_first_last_separator_queue_parallel_private_rows =
      separator_queue_parallel_private_rows;
    solver->kls_first_last_separator_queue_parallel_private_threads =
      separator_queue_parallel_private_threads;
  }
  if (separator_queue_parallel_pipeline) {
    solver->kls_first_last_separator_queue_parallel_pipeline = 1;
    solver->kls_first_separator_queue_parallel_pipeline_run_count++;
    solver->kls_first_last_separator_queue_parallel_pipeline_rows =
      separator_queue_parallel_pipeline_rows;
    solver->kls_first_last_separator_queue_parallel_pipeline_threads =
      separator_queue_parallel_pipeline_threads;
  }
  if (separator_queue_pipeline_partial) {
    solver->kls_first_last_separator_queue_pipeline_partial = 1;
    solver->kls_first_separator_queue_pipeline_partial_run_count++;
    solver->kls_first_last_separator_queue_pipeline_partial_rows =
      separator_queue_pipeline_partial_rows;
    solver->kls_first_last_separator_queue_pipeline_partial_threads =
      separator_queue_pipeline_partial_threads;
  }
  if (separator_queue_pipeline_wait_partial) {
    solver->kls_first_last_separator_queue_pipeline_wait_partial = 1;
    solver->kls_first_separator_queue_pipeline_wait_partial_run_count++;
    solver->kls_first_last_separator_queue_pipeline_wait_partial_rows =
      separator_queue_pipeline_wait_partial_rows;
    solver->kls_first_last_separator_queue_pipeline_wait_partial_deps =
      separator_queue_pipeline_wait_partial_deps;
  }
  if (separator_queue_pipeline_supernode_update) {
    solver->kls_first_last_separator_queue_pipeline_supernode_update = 1;
    solver->kls_first_separator_queue_pipeline_supernode_update_run_count++;
    solver->kls_first_last_separator_queue_pipeline_supernode_update_groups =
      separator_queue_pipeline_supernode_update_groups;
    solver->kls_first_last_separator_queue_pipeline_supernode_update_rows =
      separator_queue_pipeline_supernode_update_rows;
  }
  if (separator_queue_pipeline_supernode_panel_update) {
    solver->kls_first_last_separator_queue_pipeline_supernode_panel_update =
      1;
    solver->kls_first_separator_queue_pipeline_supernode_panel_update_run_count++;
    solver->kls_first_last_separator_queue_pipeline_supernode_panel_update_groups =
      separator_queue_pipeline_supernode_panel_update_groups;
    solver->kls_first_last_separator_queue_pipeline_supernode_panel_update_rows =
      separator_queue_pipeline_supernode_panel_update_rows;
  }
  if (separator_queue_pipeline_pivot_tail) {
    solver->kls_first_last_separator_queue_pipeline_pivot_tail = 1;
    solver->kls_first_separator_queue_pipeline_pivot_tail_run_count++;
    solver->kls_first_last_separator_queue_pipeline_pivot_tail_rows =
      separator_queue_pipeline_pivot_tail_rows;
    solver->kls_first_last_separator_queue_pipeline_pivot_restarts =
      separator_queue_pipeline_pivot_restarts;
    solver->kls_first_separator_queue_pipeline_pivot_restart_count +=
      separator_queue_pipeline_pivot_restarts;
    solver->kls_first_last_separator_queue_pipeline_pivot_serial_rows =
      separator_queue_pipeline_pivot_serial_rows;
  }
  if (separator_queue_pipeline_prefix_panel_rebuild) {
    solver->kls_first_last_separator_queue_pipeline_prefix_panel_rebuild = 1;
    solver->kls_first_separator_queue_pipeline_prefix_panel_rebuild_count +=
      separator_queue_pipeline_prefix_panel_rebuilds;
    solver->kls_first_last_separator_queue_pipeline_prefix_panel_rebuild_rows =
      separator_queue_pipeline_prefix_panel_rebuild_rows;
  }
  free(row_counts);
  free(row_ptr);
  free(x);
  free(mark);
  free(pattern);
  free(dep_heap);
  free(u_row_ptr);
  free(u_row_end);
  free(row_done);
  free(udiag_values);
  free(q_order);
  free(saved_q);
  free(col_pos);
  free(separator_component_last);
  return 1;

fail:
  if (separator_dynamic_column_rejects > 0u) {
    solver->kls_first_last_separator_dynamic_column_fallbacks =
      separator_dynamic_column_fallbacks;
    solver->kls_first_separator_dynamic_column_fallback_count +=
      separator_dynamic_column_fallbacks;
    solver->kls_first_last_separator_dynamic_column_rejects =
      separator_dynamic_column_rejects;
    solver->kls_first_separator_dynamic_column_reject_count +=
      separator_dynamic_column_rejects;
  }
  kls_row_first_refactor_seed_clear(&row_refactor_seed);
  if (q_committed && saved_q != NULL && solver->symbolic != NULL &&
      solver->symbolic->Q != NULL) {
    memcpy(solver->symbolic->Q, saved_q, (size_t)n * sizeof(*saved_q));
  }
  free(row_counts);
  free(row_ptr);
  free(x);
  free(mark);
  free(pattern);
  free(dep_heap);
  free(u_row_ptr);
  free(u_row_end);
  free(row_done);
  free(udiag_values);
  free(q_order);
  free(saved_q);
  free(col_pos);
  free(separator_component_last);
  trilinos_klu_l_free_numeric(&numeric, common);
  solver->numeric = NULL;
  *common = saved_common;
  return 0;
}

static int kls_try_first_factor_row_uplooking_blocks(
  kls_solver *solver,
  double *numeric_values) {
  if (solver == NULL || solver->symbolic == NULL || numeric_values == NULL ||
      solver->numeric != NULL || solver->symbolic->nblocks == 0u ||
      solver->symbolic->P == NULL || solver->symbolic->Q == NULL ||
      solver->symbolic->R == NULL) {
    return 0;
  }

  const UF_long n = solver->n;
  const UF_long maxblock = solver->symbolic->maxblock;
  if (n == 0u || maxblock == 0u) {
    return 0;
  }

  kls_record_first_separator_queue_plans(solver, solver->options.threads);
  if (solver->options.threads > 1 &&
      kls_try_first_factor_row_uplooking_blocks_impl(solver,
                                                     numeric_values, 1)) {
    return 1;
  }
  return kls_try_first_factor_row_uplooking_blocks_impl(solver,
                                                        numeric_values, 0);
}

typedef struct kls_pivot_first_parallel_shared {
  kls_solver *solver;
  trilinos_klu_l_numeric *numeric;
  double *numeric_values;
  const UF_long *psinv;
  const UF_long *offp_seed;
  TRILINOS_KLU_common common_template;
  atomic_ulong next_block;
  atomic_int failed;
  pthread_mutex_t commit_lock;
  UF_long total_lnz;
  UF_long total_unz;
  UF_long max_lnz_block;
  UF_long max_unz_block;
  UF_long noffdiag;
  UF_long nrealloc;
  UF_long committed_blocks;
  UF_long singular_rank;
  UF_long singular_col;
  int singular;
  int scaled;
} kls_pivot_first_parallel_shared;

typedef struct kls_pivot_first_parallel_worker {
  kls_pivot_first_parallel_shared *shared;
  double *x;
  UF_long *iwork;
  UF_long *pblock;
  UF_long *offp;
} kls_pivot_first_parallel_worker;

static void kls_pivot_first_parallel_worker_free(
  kls_pivot_first_parallel_worker *worker) {
  if (worker == NULL) {
    return;
  }
  free(worker->x);
  free(worker->iwork);
  free(worker->offp);
  memset(worker, 0, sizeof(*worker));
}

static int kls_pivot_first_parallel_worker_init(
  kls_pivot_first_parallel_worker *worker,
  kls_pivot_first_parallel_shared *shared) {
  if (worker == NULL || shared == NULL || shared->solver == NULL ||
      shared->solver->symbolic == NULL || shared->offp_seed == NULL) {
    return 0;
  }
  memset(worker, 0, sizeof(*worker));
  worker->shared = shared;
  const UF_long n = shared->solver->n;
  const UF_long maxblock = shared->solver->symbolic->maxblock;
  if (maxblock == 0u || n == KLS_KLU_EMPTY ||
      maxblock > (UF_long)(SIZE_MAX / (6u * sizeof(*worker->iwork))) ||
      n + 1u == 0u || n + 1u > (UF_long)(SIZE_MAX / sizeof(*worker->offp))) {
    return 0;
  }
  worker->x = (double *)calloc((size_t)maxblock, sizeof(*worker->x));
  worker->iwork =
    (UF_long *)malloc(6u * (size_t)maxblock * sizeof(*worker->iwork));
  worker->offp =
    (UF_long *)malloc(((size_t)n + 1u) * sizeof(*worker->offp));
  if (worker->x == NULL || worker->iwork == NULL ||
      worker->offp == NULL) {
    kls_pivot_first_parallel_worker_free(worker);
    return 0;
  }
  worker->pblock = worker->iwork + 5u * (size_t)maxblock;
  memcpy(worker->offp, shared->offp_seed,
         ((size_t)n + 1u) * sizeof(*worker->offp));
  return 1;
}

static int kls_build_pivot_first_offp_seed(
  const kls_solver *solver,
  const trilinos_klu_l_numeric *numeric,
  const UF_long *psinv,
  UF_long *offp_seed) {
  if (solver == NULL || solver->symbolic == NULL || numeric == NULL ||
      psinv == NULL || offp_seed == NULL || solver->symbolic->Q == NULL ||
      solver->symbolic->R == NULL) {
    return 0;
  }
  UF_long poff = 0;
  UF_long block = 0;
  for (UF_long k = 0; k < solver->n; ++k) {
    while (block + 1u < solver->symbolic->nblocks &&
           solver->symbolic->R[block + 1u] <= k) {
      block++;
    }
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    const UF_long oldcol = solver->symbolic->Q[k];
    if (oldcol >= solver->n || k < k1 || k >= k2) {
      return 0;
    }
    offp_seed[k] = poff;
    for (UF_long p = solver->col_ptr[oldcol];
         p < solver->col_ptr[oldcol + 1u]; ++p) {
      const UF_long oldrow = solver->row_idx[p];
      if (oldrow >= solver->n) {
        return 0;
      }
      const UF_long row_pos = psinv[oldrow];
      if (row_pos >= solver->n) {
        return 0;
      }
      if (row_pos < k1) {
        if (poff >= numeric->nzoff) {
          return 0;
        }
        poff++;
      } else if (row_pos >= k2) {
        return 0;
      }
    }
  }
  offp_seed[solver->n] = poff;
  return poff == numeric->nzoff;
}

static void kls_pivot_first_parallel_record_singular(
  kls_pivot_first_parallel_shared *shared,
  UF_long rank,
  UF_long col) {
  if (shared == NULL) {
    return;
  }
  if (!shared->singular || rank < shared->singular_rank) {
    shared->singular = 1;
    shared->singular_rank = rank;
    shared->singular_col = col;
  }
  if (shared->solver != NULL &&
      shared->solver->zero_pivot_collect != NULL) {
    const long slot = atomic_fetch_add_explicit(
      &shared->solver->zero_pivot_collect_count, 1, memory_order_relaxed);
    if (slot < shared->solver->zero_pivot_collect_cap) {
      shared->solver->zero_pivot_collect[slot] = rank;
    }
  }
}

static int kls_pivot_first_parallel_commit_block(
  kls_pivot_first_parallel_worker *worker,
  UF_long block,
  UF_long k1,
  UF_long nk,
  Unit *new_lu,
  size_t new_size,
  UF_long lnz_block,
  UF_long unz_block,
  const UF_long *pblock,
  TRILINOS_KLU_common *local_common) {
  if (worker == NULL || worker->shared == NULL ||
      worker->shared->solver == NULL || worker->shared->numeric == NULL ||
      pblock == NULL || local_common == NULL ||
      block >= worker->shared->numeric->nblocks) {
    if (new_lu != NULL && local_common != NULL) {
      (void)TRILINOS_KLU_free(new_lu, new_size, sizeof(Unit), local_common);
    }
    return 0;
  }
  kls_pivot_first_parallel_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  trilinos_klu_l_numeric *numeric = shared->numeric;
  int ok = 1;
  pthread_mutex_lock(&shared->commit_lock);
  if (atomic_load_explicit(&shared->failed, memory_order_acquire)) {
    ok = 0;
  } else if (new_lu == NULL || new_size == 0u) {
    ok = 0;
  } else {
    for (UF_long k = 0; k < nk; ++k) {
      const UF_long local_row = pblock[k];
      if (local_row >= nk || k1 + local_row >= solver->n) {
        ok = 0;
        break;
      }
    }
  }
  if (ok && local_common->memusage >
              SIZE_MAX - solver->common.memusage) {
    ok = 0;
  }
  if (ok) {
    const size_t base_mem = solver->common.memusage;
    numeric->LUbx[block] = new_lu;
    numeric->LUsize[block] = new_size;
    for (UF_long k = 0; k < nk; ++k) {
      numeric->Pnum[k1 + k] = solver->symbolic->P[k1 + pblock[k]];
    }
    shared->total_lnz += lnz_block;
    shared->total_unz += unz_block;
    shared->max_lnz_block =
      shared->max_lnz_block < lnz_block ? lnz_block : shared->max_lnz_block;
    shared->max_unz_block =
      shared->max_unz_block < unz_block ? unz_block : shared->max_unz_block;
    shared->noffdiag += (UF_long)local_common->noffdiag;
    shared->nrealloc += (UF_long)local_common->nrealloc;
    solver->common.memusage += local_common->memusage;
    if (solver->common.mempeak < solver->common.memusage) {
      solver->common.mempeak = solver->common.memusage;
    }
    if (local_common->mempeak > 0u &&
        local_common->mempeak <= SIZE_MAX - base_mem &&
        solver->common.mempeak < base_mem + local_common->mempeak) {
      solver->common.mempeak = base_mem + local_common->mempeak;
    }
    if (local_common->status == TRILINOS_KLU_SINGULAR) {
      kls_pivot_first_parallel_record_singular(
        shared, (UF_long)local_common->numerical_rank,
        (UF_long)local_common->singular_col);
    }
    if (solver->symbolic->Lnz[block] < 0.0) {
      solver->symbolic->Lnz[block] =
        (double)(lnz_block > unz_block ? lnz_block : unz_block);
    }
    shared->committed_blocks++;
    new_lu = NULL;
  }
  if (!ok) {
    atomic_store_explicit(&shared->failed, 1, memory_order_release);
  }
  pthread_mutex_unlock(&shared->commit_lock);
  if (new_lu != NULL) {
    (void)TRILINOS_KLU_free(new_lu, new_size, sizeof(Unit), local_common);
  }
  return ok;
}

static int kls_pivot_first_parallel_commit_singleton(
  kls_pivot_first_parallel_worker *worker,
  UF_long block,
  UF_long k1,
  double pivot,
  int singular,
  UF_long singular_col) {
  if (worker == NULL || worker->shared == NULL ||
      worker->shared->solver == NULL || worker->shared->numeric == NULL ||
      block >= worker->shared->numeric->nblocks) {
    return 0;
  }
  kls_pivot_first_parallel_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  trilinos_klu_l_numeric *numeric = shared->numeric;
  int ok = 1;
  pthread_mutex_lock(&shared->commit_lock);
  if (atomic_load_explicit(&shared->failed, memory_order_acquire)) {
    ok = 0;
  } else {
    ((Entry *)numeric->Udiag)[k1] = pivot;
    numeric->Pnum[k1] = solver->symbolic->P[k1];
    numeric->Lip[k1] = 0;
    numeric->Uip[k1] = 0;
    numeric->Llen[k1] = 0;
    numeric->Ulen[k1] = 0;
    shared->total_lnz++;
    shared->total_unz++;
    if (solver->symbolic->Lnz[block] < 0.0) {
      solver->symbolic->Lnz[block] = 1.0;
    }
    if (singular) {
      kls_pivot_first_parallel_record_singular(shared, k1, singular_col);
    }
    shared->committed_blocks++;
  }
  if (!ok) {
    atomic_store_explicit(&shared->failed, 1, memory_order_release);
  }
  pthread_mutex_unlock(&shared->commit_lock);
  return ok;
}

static int kls_pivot_first_parallel_factor_singleton(
  kls_pivot_first_parallel_worker *worker,
  UF_long block,
  UF_long k1) {
  if (worker == NULL || worker->shared == NULL ||
      worker->shared->solver == NULL || worker->shared->numeric == NULL) {
    return 0;
  }
  kls_pivot_first_parallel_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  trilinos_klu_l_numeric *numeric = shared->numeric;
  const UF_long oldcol = solver->symbolic->Q[k1];
  const UF_long oldrow = solver->symbolic->P[k1];
  if (oldcol >= solver->n || oldrow >= solver->n) {
    return 0;
  }
  double pivot = 0.0;
  for (UF_long p = solver->col_ptr[oldcol];
       p < solver->col_ptr[oldcol + 1u]; ++p) {
    const UF_long row = solver->row_idx[p];
    if (row >= solver->n) {
      return 0;
    }
    const UF_long row_pos = shared->psinv[row];
    if (row_pos == k1) {
      pivot = shared->numeric_values[p];
      if (shared->scaled) {
        const double rs = numeric->Rs[row];
        if (rs == 0.0) {
          return 0;
        }
        pivot /= rs;
      }
    } else if (row_pos > k1) {
      return 0;
    }
  }
  const int singular = pivot == 0.0;
  if (singular && shared->common_template.halt_if_singular) {
    return 0;
  }
  return kls_pivot_first_parallel_commit_singleton(
    worker, block, k1, pivot, singular, oldcol);
}

static int kls_pivot_first_parallel_factor_block(
  kls_pivot_first_parallel_worker *worker,
  UF_long block) {
  if (worker == NULL || worker->shared == NULL ||
      worker->shared->solver == NULL ||
      worker->shared->solver->symbolic == NULL ||
      worker->shared->numeric == NULL) {
    return 0;
  }
  kls_pivot_first_parallel_shared *shared = worker->shared;
  kls_solver *solver = shared->solver;
  trilinos_klu_l_numeric *numeric = shared->numeric;
  if (block >= solver->symbolic->nblocks) {
    return 0;
  }
  const UF_long k1 = solver->symbolic->R[block];
  const UF_long k2 = solver->symbolic->R[block + 1u];
  if (k2 <= k1 || k2 > solver->n ||
      k2 - k1 > solver->symbolic->maxblock) {
    return 0;
  }
  const UF_long nk = k2 - k1;
  if (nk == 1u) {
    return kls_pivot_first_parallel_factor_singleton(worker, block, k1);
  }

  double lsize =
    solver->symbolic->Lnz[block] < 0.0
      ? -(shared->common_template.initmem)
      : shared->common_template.initmem_amd * solver->symbolic->Lnz[block] +
          (double)nk;
  TRILINOS_KLU_common local_common = shared->common_template;
  local_common.status = TRILINOS_KLU_OK;
  local_common.numerical_rank = KLS_KLU_EMPTY;
  local_common.singular_col = KLS_KLU_EMPTY;
  local_common.noffdiag = 0;
  local_common.nrealloc = 0;
  local_common.memusage = 0u;
  local_common.mempeak = 0u;

  Unit *new_lu = NULL;
  UF_long lnz_block = 0;
  UF_long unz_block = 0;
  const size_t new_size =
    TRILINOS_KLU_kernel_factor(nk, solver->col_ptr, solver->row_idx,
                               shared->numeric_values, solver->symbolic->Q,
                               lsize, &new_lu,
                               ((Entry *)numeric->Udiag) + k1,
                               numeric->Llen + k1, numeric->Ulen + k1,
                               numeric->Lip + k1, numeric->Uip + k1,
                               worker->pblock, &lnz_block, &unz_block,
                               (Entry *)worker->x, worker->iwork, k1,
                               (UF_long *)shared->psinv,
                               shared->scaled ? numeric->Rs : NULL,
                               worker->offp, numeric->Offi,
                               (Entry *)numeric->Offx, &local_common);
  if (new_size == 0u || new_lu == NULL ||
      local_common.status < TRILINOS_KLU_OK ||
      (local_common.status == TRILINOS_KLU_SINGULAR &&
       local_common.halt_if_singular)) {
    if (new_lu != NULL) {
      (void)TRILINOS_KLU_free(new_lu, new_size, sizeof(Unit),
                              &local_common);
    }
    return 0;
  }
  return kls_pivot_first_parallel_commit_block(
    worker, block, k1, nk, new_lu, new_size, lnz_block, unz_block,
    worker->pblock, &local_common);
}

static void *kls_pivot_first_parallel_worker_main(void *arg) {
  kls_pivot_first_parallel_worker worker;
  kls_pivot_first_parallel_shared *shared =
    (kls_pivot_first_parallel_shared *)arg;
  if (!kls_pivot_first_parallel_worker_init(&worker, shared)) {
    if (shared != NULL) {
      atomic_store_explicit(&shared->failed, 1, memory_order_release);
    }
    return NULL;
  }
  const UF_long nblocks =
    shared != NULL && shared->solver != NULL && shared->solver->symbolic != NULL
      ? shared->solver->symbolic->nblocks : 0u;
  for (;;) {
    if (atomic_load_explicit(&shared->failed, memory_order_acquire)) {
      break;
    }
    const UF_long block = (UF_long)atomic_fetch_add_explicit(
      &shared->next_block, 1ul, memory_order_acq_rel);
    if (block >= nblocks) {
      break;
    }
    if (!kls_pivot_first_parallel_factor_block(&worker, block)) {
      atomic_store_explicit(&shared->failed, 1, memory_order_release);
      break;
    }
  }
  kls_pivot_first_parallel_worker_free(&worker);
  return NULL;
}

static int kls_try_first_factor_with_pivoted_blocks_parallel(
  kls_solver *solver,
  double *numeric_values,
  const UF_long *psinv) {
  if (solver == NULL || solver->symbolic == NULL || solver->numeric == NULL ||
      numeric_values == NULL || psinv == NULL ||
      solver->symbolic->nblocks < 2u || solver->options.threads <= 1 ||
      solver->symbolic->P == NULL || solver->symbolic->Q == NULL ||
      solver->symbolic->R == NULL || solver->symbolic->Lnz == NULL ||
      solver->symbolic->maxblock == 0u || solver->n == 0u) {
    return -1;
  }
  int thread_count = solver->options.threads;
  if ((UF_long)thread_count > solver->symbolic->nblocks) {
    thread_count = (int)solver->symbolic->nblocks;
  }
  if (thread_count < 2) {
    return -1;
  }

  UF_long *offp_seed = NULL;
  pthread_t *threads = NULL;
  int mutex_initialized = 0;
  kls_pivot_first_parallel_shared shared;
  memset(&shared, 0, sizeof(shared));
  shared.singular_rank = KLS_KLU_EMPTY;
  shared.singular_col = KLS_KLU_EMPTY;
  shared.max_lnz_block = 1u;
  shared.max_unz_block = 1u;

  if (solver->n + 1u == 0u ||
      solver->n + 1u > (UF_long)(SIZE_MAX / sizeof(*offp_seed))) {
    return 0;
  }
  offp_seed =
    (UF_long *)malloc(((size_t)solver->n + 1u) * sizeof(*offp_seed));
  if (offp_seed == NULL ||
      !kls_build_pivot_first_offp_seed(solver, solver->numeric, psinv,
                                       offp_seed)) {
    free(offp_seed);
    return 0;
  }
  if (pthread_mutex_init(&shared.commit_lock, NULL) != 0) {
    free(offp_seed);
    return 0;
  }
  mutex_initialized = 1;
  threads = (pthread_t *)calloc((size_t)thread_count, sizeof(*threads));
  if (threads == NULL) {
    goto fail;
  }

  atomic_init(&shared.next_block, 0ul);
  atomic_init(&shared.failed, 0);
  shared.solver = solver;
  shared.numeric = solver->numeric;
  shared.numeric_values = numeric_values;
  shared.psinv = psinv;
  shared.offp_seed = offp_seed;
  shared.common_template = solver->common;
  shared.scaled = solver->common.scale > 0;

  int created = 0;
  for (; created < thread_count; ++created) {
    if (pthread_create(&threads[created], NULL,
                       kls_pivot_first_parallel_worker_main, &shared) != 0) {
      atomic_store_explicit(&shared.failed, 1, memory_order_release);
      break;
    }
  }
  for (int i = 0; i < created; ++i) {
    (void)pthread_join(threads[i], NULL);
  }
  if (created != thread_count ||
      atomic_load_explicit(&shared.failed, memory_order_acquire) ||
      shared.committed_blocks != solver->symbolic->nblocks) {
    goto fail;
  }

  solver->numeric->lnz = shared.total_lnz;
  solver->numeric->unz = shared.total_unz;
  solver->numeric->max_lnz_block = shared.max_lnz_block;
  solver->numeric->max_unz_block = shared.max_unz_block;
  solver->common.noffdiag = shared.noffdiag;
  solver->common.nrealloc = shared.nrealloc;
  if (shared.singular) {
    solver->common.status = TRILINOS_KLU_SINGULAR;
    solver->common.numerical_rank = shared.singular_rank;
    solver->common.singular_col = shared.singular_col;
  } else {
    solver->common.status = TRILINOS_KLU_OK;
    solver->common.numerical_rank = solver->n;
    solver->common.singular_col = solver->n;
  }
  if (!kls_rebuild_numeric_pinv(solver) ||
      !kls_recompute_offdiag_from_pinv(solver, numeric_values)) {
    goto fail;
  }
  if (shared.scaled && !kls_parallel_refactor_permute_scale(solver)) {
    goto fail;
  }
  solver->kls_first_last_parallel_btf_blocks = shared.committed_blocks;
  solver->kls_first_parallel_btf_block_count += shared.committed_blocks;

  if (mutex_initialized) {
    (void)pthread_mutex_destroy(&shared.commit_lock);
  }
  free(threads);
  free(offp_seed);
  return 1;

fail:
  if (threads != NULL) {
    free(threads);
  }
  if (mutex_initialized) {
    (void)pthread_mutex_destroy(&shared.commit_lock);
  }
  free(offp_seed);
  return 0;
}

static int kls_try_first_factor_with_pivoted_blocks(kls_solver *solver,
                                                    double *numeric_values) {
  if (solver == NULL || solver->symbolic == NULL || numeric_values == NULL ||
      solver->numeric != NULL || solver->symbolic->P == NULL ||
      solver->symbolic->Q == NULL ||
      solver->symbolic->R == NULL || solver->symbolic->Lnz == NULL) {
    return 0;
  }

  TRILINOS_KLU_common *common = &solver->common;
  trilinos_klu_l_common saved_common = *common;
  common->initmem_amd = common->initmem_amd > 1.0 ? common->initmem_amd : 1.0;
  /* tiny systems: growth reallocs cost more than the memory (rajat03:
     the 4.3ms serial factor pays mid-factor LU grows; 3x of a <1MB
     estimate is free) */
  if (solver->n > 0 && solver->n <= 32768 && common->initmem_amd < 1.5) {
    common->initmem_amd = 1.5;
  }
  common->initmem = common->initmem > 1.0 ? common->initmem : 1.0;
  common->tol = common->tol < 1.0 ? common->tol : 1.0;
  common->tol = common->tol > 0.0 ? common->tol : 0.0;
  common->memgrow = common->memgrow > 1.0 ? common->memgrow : 1.0;
  common->status = TRILINOS_KLU_OK;
  common->numerical_rank = KLS_KLU_EMPTY;
  common->singular_col = KLS_KLU_EMPTY;
  common->noffdiag = 0;
  common->nrealloc = 0;

  trilinos_klu_l_numeric *numeric = kls_allocate_numeric_skeleton(solver);
  if (numeric == NULL) {
    *common = saved_common;
    return 0;
  }
  solver->numeric = numeric;
  const int scaled = common->scale > 0;
  if (common->scale >= 0 &&
      !trilinos_klu_l_scale((UF_long)common->scale, solver->n,
                            solver->col_ptr, solver->row_idx,
                            numeric_values, scaled ? numeric->Rs : NULL,
                            numeric->Pnum, common)) {
    goto fail;
  }

  UF_long *psinv = numeric->Pinv;
  for (UF_long k = 0; k < solver->n; ++k) {
    psinv[k] = KLS_KLU_EMPTY;
  }
  for (UF_long k = 0; k < solver->n; ++k) {
    const UF_long row = solver->symbolic->P[k];
    if (row >= solver->n) {
      goto fail;
    }
    psinv[row] = k;
  }
  for (UF_long k = 0; k < solver->n; ++k) {
    if (psinv[k] == KLS_KLU_EMPTY) {
      goto fail;
    }
  }
  (void)kls_build_refactor_map(solver);

  const int parallel_status =
    kls_try_first_factor_with_pivoted_blocks_parallel(solver,
                                                      numeric_values,
                                                      psinv);
  if (parallel_status > 0) {
    free_refactor_map(solver);
    return 1;
  }
  if (parallel_status == 0) {
    goto fail;
  }

  UF_long total_lnz = 0;
  UF_long total_unz = 0;
  UF_long max_lnz_block = 1;
  UF_long max_unz_block = 1;
  UF_long *old_pblock = numeric->Iwork;
  UF_long *pblock =
    numeric->Iwork + 5u * (size_t)solver->symbolic->maxblock;
  for (UF_long block = 0; block < solver->symbolic->nblocks; ++block) {
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    if (k2 < k1 || k2 > solver->n) {
      goto fail;
    }
    const UF_long nk = k2 - k1;
    if (nk == 0u || nk > solver->symbolic->maxblock) {
      goto fail;
    }
    if (nk == 1u) {
      const UF_long oldcol = solver->symbolic->Q[k1];
      const UF_long oldrow = solver->symbolic->P[k1];
      if (oldcol >= solver->n || oldrow >= solver->n) {
        goto fail;
      }
      double pivot = 0.0;
      for (UF_long p = solver->col_ptr[oldcol];
           p < solver->col_ptr[oldcol + 1u]; ++p) {
        const UF_long row = solver->row_idx[p];
        if (row >= solver->n) {
          goto fail;
        }
        const UF_long row_pos = psinv[row];
        if (row_pos == k1) {
          pivot = numeric_values[p];
          if (scaled) {
            const double rs = numeric->Rs[row];
            if (rs == 0.0) {
              goto fail;
            }
            pivot /= rs;
          }
        } else if (row_pos > k1) {
          goto fail;
        }
      }
      ((Entry *)numeric->Udiag)[k1] = pivot;
      numeric->Pnum[k1] = oldrow;
      numeric->Lip[k1] = 0;
      numeric->Uip[k1] = 0;
      numeric->Llen[k1] = 0;
      numeric->Ulen[k1] = 0;
      total_lnz++;
      total_unz++;
      if (solver->symbolic->Lnz[block] < 0.0) {
        solver->symbolic->Lnz[block] = 1.0;
      }
      if (pivot == 0.0) {
        common->status = TRILINOS_KLU_SINGULAR;
        if (common->numerical_rank == KLS_KLU_EMPTY) {
          common->numerical_rank = k1;
          common->singular_col = oldcol;
        }
        if (common->halt_if_singular) {
          goto fail;
        }
      }
      continue;
    }
    for (UF_long k = 0; k < nk; ++k) {
      old_pblock[k] = k;
    }

    Unit *new_lu = NULL;
    size_t new_size = 0u;
    UF_long lnz_block = 0;
    UF_long unz_block = 0;
    UF_long skipped_columns = 0;
    double skipped_work = 0.0;
    int exact_mask_used = 0;
    if (!kls_try_pivot_tail_restart_rejected_block(
          solver, numeric_values, block, k1, nk, 0u, nk, NULL, old_pblock,
          psinv, &new_lu, &new_size, &lnz_block, &unz_block, pblock,
          &skipped_columns, &skipped_work, &exact_mask_used) ||
        new_lu == NULL || new_size == 0u || common->status < 0 ||
        (common->status == TRILINOS_KLU_SINGULAR &&
         common->halt_if_singular)) {
      goto fail;
    }
    numeric->LUbx[block] = new_lu;
    numeric->LUsize[block] = new_size;
    for (UF_long k = 0; k < nk; ++k) {
      const UF_long local_row = pblock[k];
      if (local_row >= nk || k1 + local_row >= solver->n) {
        goto fail;
      }
      numeric->Pnum[k1 + k] = solver->symbolic->P[k1 + local_row];
    }
    total_lnz += lnz_block;
    total_unz += unz_block;
    max_lnz_block =
      max_lnz_block < lnz_block ? lnz_block : max_lnz_block;
    max_unz_block =
      max_unz_block < unz_block ? unz_block : max_unz_block;
    if (solver->symbolic->Lnz[block] < 0.0) {
      solver->symbolic->Lnz[block] =
        (double)(lnz_block > unz_block ? lnz_block : unz_block);
    }
  }

  numeric->lnz = total_lnz;
  numeric->unz = total_unz;
  numeric->max_lnz_block = max_lnz_block;
  numeric->max_unz_block = max_unz_block;
  free_refactor_map(solver);
  if (!kls_rebuild_numeric_pinv(solver) ||
      !kls_recompute_offdiag_from_pinv(solver, numeric_values)) {
    goto fail;
  }
  if (scaled && !kls_parallel_refactor_permute_scale(solver)) {
    goto fail;
  }
  if (common->status == TRILINOS_KLU_OK) {
    common->numerical_rank = solver->n;
    common->singular_col = solver->n;
  }
  return 1;

fail:
  free_refactor_map(solver);
  trilinos_klu_l_free_numeric(&numeric, common);
  solver->numeric = NULL;
  *common = saved_common;
  return 0;
}

static int kls_first_rebuild_quality_accepts(
  kls_solver *solver,
  const trilinos_klu_l_numeric *saved_numeric,
  const trilinos_klu_l_common *saved_common) {
  if (solver == NULL || solver->symbolic == NULL ||
      solver->numeric == NULL || saved_numeric == NULL ||
      saved_common == NULL) {
    return 0;
  }

  trilinos_klu_l_common saved_eval = *saved_common;
  (void)trilinos_klu_l_flops(solver->symbolic,
                             (trilinos_klu_l_numeric *)saved_numeric,
                             &saved_eval);
  (void)trilinos_klu_l_flops(solver->symbolic, solver->numeric,
                             &solver->common);

  const double saved_fill =
    (double)saved_numeric->lnz + (double)saved_numeric->unz;
  const double rebuilt_fill =
    (double)solver->numeric->lnz + (double)solver->numeric->unz;
  if (saved_fill > 0.0 && rebuilt_fill > 1.01 * saved_fill) {
    return 0;
  }

  const double saved_flops = saved_eval.flops;
  const double rebuilt_flops = solver->common.flops;
  if (saved_flops > 0.0 && rebuilt_flops > 1.02 * saved_flops) {
    return 0;
  }

  if (solver->common.noffdiag > saved_eval.noffdiag &&
      solver->common.noffdiag > 16u) {
    return 0;
  }
  return 1;
}

static int kls_try_rebuild_current_numeric_with_kls_first_mode(
  kls_solver *solver,
  double *numeric_values,
  double *elapsed,
  int force) {
  if (solver != NULL && solver->numeric_is_predicted) {
    /* Predicted numerics reject to the bounded KLU fallback; the row-first
       rebuild pipelines assume KLU-built state and churn unboundedly on
       these shapes. */
    return 0;
  }
  if (solver == NULL || numeric_values == NULL || solver->numeric == NULL ||
      (force ? kls_first_factor_env_disabled()
             : !kls_should_try_first_factor(solver))) {
    return 0;
  }

  trilinos_klu_l_numeric *saved_numeric = solver->numeric;
  trilinos_klu_l_common saved_common = solver->common;
  UF_long *saved_q = NULL;
  if (solver->symbolic->Q != NULL && solver->n > 0u) {
    saved_q = (UF_long *)malloc((size_t)solver->n * sizeof(*saved_q));
    if (saved_q == NULL) {
      return 0;
    }
    memcpy(saved_q, solver->symbolic->Q,
           (size_t)solver->n * sizeof(*saved_q));
  }
  free_refactor_map(solver);
  free_refactor_schedule(solver);
  free_refactor_lu_pointer_cache(solver);
  free_fast_reject_tail_plan(solver);
  free_row_refactor_pattern(solver);
  solver->row_refactor_auto_enabled = 0;
  solver->numeric = NULL;

  const double start = kls_now_seconds();
  int rebuilt = kls_try_first_factor_row_uplooking_blocks(solver,
                                                          numeric_values);
  if (!rebuilt) {
    rebuilt = kls_try_first_factor_with_pivoted_blocks(solver,
                                                       numeric_values);
  }
  if (elapsed != NULL) {
    *elapsed += kls_now_seconds() - start;
  }

  if (rebuilt && solver->numeric != NULL && solver->common.status >= 0 &&
      kls_first_rebuild_quality_accepts(solver, saved_numeric,
                                        &saved_common)) {
    free_snode_panels(solver);
    trilinos_klu_l_free_numeric(&saved_numeric, &saved_common);
    free(saved_q);
    return 1;
  }

  if (solver->numeric != NULL) {
    trilinos_klu_l_free_numeric(&solver->numeric, &solver->common);
  }
  free_refactor_map(solver);
  free_refactor_schedule(solver);
  free_refactor_lu_pointer_cache(solver);
  free_fast_reject_tail_plan(solver);
  free_row_refactor_pattern(solver);
  if (saved_q != NULL && solver->symbolic->Q != NULL) {
    memcpy(solver->symbolic->Q, saved_q,
           (size_t)solver->n * sizeof(*saved_q));
  }
  free(saved_q);
  solver->numeric = saved_numeric;
  solver->common = saved_common;
  return 0;
}

static int kls_try_rebuild_current_numeric_with_kls_first(
  kls_solver *solver,
  double *numeric_values,
  double *elapsed) {
  return kls_try_rebuild_current_numeric_with_kls_first_mode(
    solver, numeric_values, elapsed, 0);
}

/* Predicted-pattern first factorization: KLU's Gilbert-Peierls first factor
   spends two orders of magnitude more time discovering the pattern per
   column than the refactorization spends computing the same values (the
   NICSLU paper documents KLU at 462s on memchip whose refactor takes about
   a second).  For large structurally symmetric patterns with a full
   diagonal, the exact no-pivot LU pattern equals the symbolic Cholesky
   pattern of the block, which an elimination-tree reach computes in near
   linear time.  Build that pattern directly into a KLU numeric object and
   compute the first values with the ordinary refactorization; a pivot
   threshold check keeps the quality guarantees and any rejection falls
   back to the standard factorization. */
static int kls_compare_uf_long_asc(const void *a, const void *b) {
  const UF_long left = *(const UF_long *)a;
  const UF_long right = *(const UF_long *)b;
  return (left > right) - (left < right);
}

static int kls_predicted_pivoting_fill(kls_solver *solver,
                                       const double *numeric_values);
static int kls_predicted_pivot_rescue(kls_solver *solver,
                                      const double *numeric_values,
                                      double *elapsed);

/* Parallel predicted-pattern build for one large block: the etree reach
   walks and especially the 66M-scale random L-scatter of pass 2 are
   DRAM-latency-bound and column-independent.  Threads own contiguous
   column ranges with private mark arrays; pass 1 retains per-thread
   partial column counts so pass 2 threads write L entries at exclusive
   pre-scanned cursors - the ascending-k append order (and every byte of
   the result) is identical to the serial build. */
typedef struct kls_pred_build_ctx {
  UF_long nk;
  UF_long k1;
  const UF_long *ucol_ptr;
  const UF_long *urows;
  const UF_long *parent;
  UF_long *rowlen;             /* out: per-column U reach length */
  int32_t **count_partial;     /* [t][nk] pass-1 partial col counts */
  UF_long *thread_lnz;         /* [t] */
  double *lu;                  /* pass 2: block LU area */
  trilinos_klu_l_numeric *numeric;
  UF_long **lcursor;           /* [t][nk] pass-2 L write cursors */
  int nthreads;
  int phase;                   /* 1 = counts, 2 = fill */
  int failed;
} kls_pred_build_ctx;

typedef struct kls_pred_build_job {
  kls_pred_build_ctx *ctx;
  UF_long lo, hi;
  int tid;
  pthread_t thread;
} kls_pred_build_job;

static void *kls_pred_build_worker(void *arg) {
  kls_pred_build_job *job = (kls_pred_build_job *)arg;
  kls_pred_build_ctx *ctx = job->ctx;
  const UF_long nk = ctx->nk;
  const UF_long *ucol_ptr = ctx->ucol_ptr;
  const UF_long *urows = ctx->urows;
  const UF_long *parent = ctx->parent;
  int32_t *mark = (int32_t *)malloc((size_t)nk * sizeof(*mark));
  UF_long *urow = NULL;
  if (mark == NULL) {
    ctx->failed = 1;
    return NULL;
  }
  memset(mark, 0xff, (size_t)nk * sizeof(*mark));
  if (ctx->phase == 1) {
    int32_t *count = ctx->count_partial[job->tid];
    UF_long lnz = 0;
    for (UF_long k = job->lo; k < job->hi; ++k) {
      mark[k] = (int32_t)k;
      UF_long len = 0;
      for (UF_long p = ucol_ptr[k]; p < ucol_ptr[k + 1u]; ++p) {
        UF_long i = urows[p];
        while (i < k && mark[i] != (int32_t)k) {
          mark[i] = (int32_t)k;
          count[i]++;
          len++;
          i = parent[i];
          if (i == KLS_KLU_EMPTY) {
            break;
          }
        }
      }
      ctx->rowlen[k] = len;
      lnz += len;
    }
    ctx->thread_lnz[job->tid] = lnz;
  } else {
    double *lu = ctx->lu;
    trilinos_klu_l_numeric *numeric = ctx->numeric;
    const UF_long k1 = ctx->k1;
    UF_long *lcur = ctx->lcursor[job->tid];
    urow = (UF_long *)malloc((size_t)nk * sizeof(*urow));
    if (urow == NULL) {
      ctx->failed = 1;
      free(mark);
      return NULL;
    }
    for (UF_long k = job->lo; k < job->hi; ++k) {
      mark[k] = (int32_t)k;
      UF_long len = 0;
      for (UF_long p = ucol_ptr[k]; p < ucol_ptr[k + 1u]; ++p) {
        UF_long i = urows[p];
        while (i < k && mark[i] != (int32_t)k) {
          mark[i] = (int32_t)k;
          urow[len++] = i;
          i = parent[i];
          if (i == KLS_KLU_EMPTY) {
            break;
          }
        }
      }
      for (UF_long a = 1u; a < len; ++a) {
        const UF_long v = urow[a];
        UF_long b = a;
        while (b > 0u && urow[b - 1u] > v) {
          urow[b] = urow[b - 1u];
          b--;
        }
        urow[b] = v;
      }
      UF_long *ui = (UF_long *)(lu + numeric->Uip[k1 + k]);
      for (UF_long p = 0; p < len; ++p) {
        ui[p] = urow[p];
        const UF_long j = urow[p];
        UF_long *li = (UF_long *)(lu + numeric->Lip[k1 + j]);
        li[lcur[j]++] = k;
      }
    }
  }
  free(urow);
  free(mark);
  return NULL;
}

/* Run one phase across nthreads-1 spawned workers plus the caller.
   Returns 0 on any spawn or allocation failure (caller falls back). */
static int kls_pred_build_run(kls_pred_build_ctx *ctx,
                              kls_pred_build_job *jobs) {
  const int nt = ctx->nthreads;
  for (int t = 1; t < nt; ++t) {
    if (pthread_create(&jobs[t].thread, NULL, kls_pred_build_worker,
                       &jobs[t]) != 0) {
      ctx->failed = 1;
      for (int u = 1; u < t; ++u) {
        pthread_join(jobs[u].thread, NULL);
      }
      return 0;
    }
  }
  kls_pred_build_worker(&jobs[0]);
  for (int t = 1; t < nt; ++t) {
    pthread_join(jobs[t].thread, NULL);
  }
  return ctx->failed == 0;
}

/* Prediction can be worthwhile even when it merely matches KLU's cold
   factor time: a fixed-pivot numeric may remove substantially more work from
   every advertised update.  Admit that experiment from the symbolic work
   accumulated over the caller's lifecycle.  Construction, fill bounds, and
   the mandatory true-residual probe remain the acceptance authorities. */
static int kls_generic_repeated_predicted_lifecycle_candidate(
  const kls_solver *solver) {
  if (solver == NULL ||
      solver->symbolic == NULL || solver->options.threads <= 1 ||
      !kls_repeated_update_workload(&solver->options) ||
      solver->common.scale > 0) {
    return 0;
  }
  /* User-ordering symbolics do not publish fill/flop estimates.  A generic
     ND candidate with a bounded A union A' stream has already passed the
     same lifecycle-work budget before it is installed, so the missing
     estimate must not suppress its value census and exact-residual probe. */
  return (solver->generic_nd_portfolio_selected &&
          solver->generic_nd_bounded_symmetric_union) ||
    kls_symbolic_repeated_prediction_economic(
      solver->symbolic, &solver->options);
}

/* A lifecycle-funded prediction trial should not spend a full pattern build
   when the selected static diagonal is already far below its input-column
   scale.  This value census is deliberately permissive: unpressured pivots
   pass immediately, and a pressured set passes when at least three quarters
   remain above 1e-8 of their block-local column maximum.  The later factor
   and true-residual probe still make the final decision. */
static int kls_generic_predicted_diagonal_values_capable(
  const kls_solver *solver,
  const double *numeric_values) {
  if (solver == NULL || solver->symbolic == NULL || numeric_values == NULL ||
      solver->col_ptr == NULL || solver->row_idx == NULL ||
      solver->symbolic->P == NULL || solver->symbolic->Q == NULL ||
      solver->symbolic->R == NULL ||
      (solver->common.scale > 0 &&
       (solver->options.scale != KLS_SCALE_AUTO ||
        getenv("KLS_PREDICTED_TRY_SCALED") != NULL))) {
    return 0;
  }
  UF_long *pinv = (UF_long *)malloc((size_t)solver->n * sizeof(*pinv));
  if (pinv == NULL && solver->n > 0u) {
    return 0;
  }
  int valid = 1;
  for (UF_long k = 0u; k < solver->n; ++k) {
    pinv[k] = KLS_KLU_EMPTY;
  }
  for (UF_long k = 0u; k < solver->n; ++k) {
    const UF_long row = solver->symbolic->P[k];
    if (row >= solver->n || pinv[row] != KLS_KLU_EMPTY) {
      valid = 0;
      break;
    }
    pinv[row] = k;
  }
  UF_long pressured = 0u;
  UF_long recovered = 0u;
  UF_long block = 0u;
  for (UF_long k = 0u; valid && k < solver->n; ++k) {
    while (block + 1u < solver->symbolic->nblocks &&
           solver->symbolic->R[block + 1u] <= k) {
      block++;
    }
    if (block >= solver->symbolic->nblocks) {
      valid = 0;
      break;
    }
    const UF_long k1 = solver->symbolic->R[block];
    const UF_long k2 = solver->symbolic->R[block + 1u];
    const UF_long oldcol = solver->symbolic->Q[k];
    const UF_long oldrow = solver->symbolic->P[k];
    if (oldcol >= solver->n || k < k1 || k >= k2) {
      valid = 0;
      break;
    }
    double diagonal = 0.0;
    double column_max = 0.0;
    for (UF_long p = solver->col_ptr[oldcol];
         p < solver->col_ptr[oldcol + 1u]; ++p) {
      const UF_long row = solver->row_idx[p];
      const double value = fabs(numeric_values[p]);
      if (row >= solver->n || !isfinite(value)) {
        valid = 0;
        break;
      }
      const UF_long newrow = pinv[row];
      if (newrow >= k1 && newrow < k2 && value > column_max) {
        column_max = value;
      }
      if (row == oldrow && value > diagonal) {
        diagonal = value;
      }
    }
    if (!valid || !(column_max > 0.0) ||
        diagonal >= solver->options.pivot_tolerance * column_max) {
      continue;
    }
    pressured++;
    recovered += diagonal >= 1.0e-8 * column_max;
  }
  free(pinv);
  const int capable = valid &&
    (pressured == 0u || recovered >= pressured - pressured / 4u);
  if (getenv("KLS_TRACE_PREDICTED") != NULL) {
    fprintf(stderr,
            "KLS predicted lifecycle value census: valid=%d pressure="
            "%ld recovered=%ld capable=%d\n",
            valid, (long)pressured, (long)recovered, capable);
  }
  return capable;
}

static int kls_predicted_probe_first_threshold_shape(
  const kls_solver *solver) {
  {
    /* Threshold violations describe a multiplier bound, not an inaccurate
       solution.  For a lifecycle-funded generic trial, probe the exact
       unmodified numeric first; accepting still requires the ordinary true
       residual below, while a failure discards the trial and returns to KLU.
       This avoids value-dependent diagonal nudges becoming an implicit input
       classifier. */
    return kls_generic_repeated_predicted_lifecycle_candidate(solver);
  }

  if (solver == NULL || solver->symbolic == NULL || solver->n < 100000u ||
      solver->n > 150000u || solver->nnz > 10u * solver->n ||
      solver->symbolic->nblocks < 2u ||
      solver->symbolic->nblocks > 64u ||
      solver->symbolic->maxblock < 100000u ||
      solver->common.scale > 0 || solver->common.flops < 2.0e7 ||
      solver->common.flops >= 1.0e8) {
    return 0;
  }
  /* A nearly whole-matrix BTF block with modest work can carry a handful of
     large-but-benign multipliers.  Diagonal nudges made that numeric less
     accurate on the dc/trans cohort, whereas the unmodified factor passed
     the mandatory true-residual probe and repeated changing-value checks.
     Let the probe arbitrate this narrow structural class first. */
  return 1000.0 * (double)solver->symbolic->maxblock >=
         995.0 * (double)solver->n;
}

static int kls_predicted_pattern_first_factor(kls_solver *solver,
                                               double *numeric_values,
                                               double *elapsed) {
  /* Prediction is a parallel pattern-construction optimization.  In serial
     it is slower than KLU on the ASIC class, and its suffix-closed storage
     is not a valid input to the later mapped-refactor portfolio. */
  if (solver == NULL || solver->symbolic == NULL || solver->numeric != NULL ||
      numeric_values == NULL ||
      solver->options.threads <= 1 ||
      (solver->n < 20000 && !solver->block_trial_active) ||
      solver->options.pivot_tolerance <= 0.0) {
    if (solver != NULL && solver->block_trial_active &&
        getenv("KLS_TRACE_PREDICTED") != NULL) {
      fprintf(stderr, "KLS block: predicted gate1\n");
    }
    return 0;
  }
  if (solver->common.scale > 0 &&
      solver->options.scale != KLS_SCALE_AUTO) {
    return 0;
  }
  trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  trilinos_klu_l_common *common = &solver->common;
  const UF_long n = solver->n;
  const UF_long nblocks = symbolic->nblocks;
  const UF_long maxblock = symbolic->maxblock;
  const UF_long *P = symbolic->P;
  const UF_long *Q = symbolic->Q;
  const UF_long *R = symbolic->R;
  if (P == NULL || Q == NULL || R == NULL || maxblock < 2) {
    if (solver->block_trial_active && getenv("KLS_TRACE_PREDICTED") != NULL) {
      fprintf(stderr, "KLS block: predicted gate2 maxblock=%ld nblocks=%ld\n",
              (long)maxblock, (long)nblocks);
    }
    return 0;
  }
  /* The KLU first factorization only becomes pattern-discovery bound at
     scale; below it the prediction attempt (build + fill + probe) costs
     more than it saves and the symmetrized pattern can carry extra fill
     into every refactorization. */
  {
    const char *dis = getenv("KLS_DISABLE_PREDICTED_FIRST");
    if (dis != NULL && dis[0] == '1') {
      return 0;
    }
  }
  const int generic_nd_unscored_bounded =
    solver->generic_nd_portfolio_selected &&
    solver->generic_nd_bounded_symmetric_union &&
    !(symbolic->lnz + symbolic->unz > 0.0);
  const int generic_nd_static_values_capable =
    generic_nd_unscored_bounded
      ? kls_generic_predicted_diagonal_values_capable(solver, numeric_values)
      : 1;
  const int generic_nd_fill_only =
    generic_nd_unscored_bounded &&
    solver->options.scale == KLS_SCALE_AUTO &&
    (solver->common.scale > 0 || !generic_nd_static_values_capable);
  const int generic_scaled_auto_prediction_census =
    solver->options.scale == KLS_SCALE_AUTO &&
    solver->common.scale > 0 &&
    getenv("KLS_PREDICTED_TRY_SCALED") == NULL &&
    solver->options.threads > 1 &&
    kls_repeated_update_workload(&solver->options) &&
    kls_symbolic_repeated_prediction_economic(
      solver->symbolic, &solver->options);
  if (solver->generic_nd_portfolio_selected &&
      !solver->generic_nd_bounded_symmetric_union) {
    /* On an ND frame with an unbounded A+A' stream, prediction's residual
       may be excellent while the extra union fill is paid by every later
       update.  Let ordinary pivoted KLU discover the actual LU pattern; the
       ND numeric-fill rollback remains active. */
    return 0;
  }
  if ((kls_generic_repeated_predicted_lifecycle_candidate(solver) ||
       generic_scaled_auto_prediction_census) &&
      !kls_generic_predicted_diagonal_values_capable(
        solver, numeric_values) &&
      !(solver->generic_nd_portfolio_selected &&
        solver->generic_nd_bounded_symmetric_union)) {
    /* For an ordinary minimum-degree proposal the census avoids a pattern
       build whose selected diagonal is predictably unusable.  AUTO may have
       selected a positive robust scale, but the prediction arm below still
       runs unscaled; apply the same current-value pressure certificate before
       paying for it.  A provisional bounded-union ND ordering is different:
       its otherwise-unavailable fill evidence and mandatory exact-residual
       probe are the numeric validation for the ordering itself.  Let that
       already-funded trial proceed; any failed probe still restores the
       recorded incumbent. */
    return 0;
  }

  if (n < 500000 && !(symbolic->lnz >= 5.0e6) &&
      !solver->generic_nd_portfolio_selected &&
      !solver->block_trial_active && !solver->prestatic_adopted_unfactored &&
      !kls_generic_repeated_predicted_lifecycle_candidate(solver) &&
      (solver->metis_race == NULL || n >= 1000000) &&
      getenv("KLS_FORCE_PIVOT_FILL") == NULL) {
    /* With a METIS race pending this numeric is a bootstrap the
       promotion will replace: the parallel predicted fill is far
       cheaper than a serial KLU first factorization, and a predicted
       numeric also keeps the unmeasured scale re-trials off. */
    if (getenv("KLS_TRACE_PREDICTED") != NULL) {
      fprintf(stderr, "KLS predicted: skipped n=%ld est_lnz=%.3e\n",
              (long)n, symbolic->lnz);
    }
    return 0;
  }

  {
    /* Static-diagonal prediction requires the pivot diagonal to be
       structurally PRESENT under the symbolic's P/Q: patterns whose
       diagonals are mostly missing (mac_econ: 89% under the raced
       METIS symbolic) factor to exact-cancellation zero pivots and the
       honest probe measures e+33 - the whole attempt (~1.5s of build,
       preps and passes) is a priori hopeless.  O(nnz) scan. */
    UF_long missing = 0;
    UF_long k, p;
    const UF_long *P = symbolic->P;
    UF_long *pinv_probe = (UF_long *)malloc((size_t)n * sizeof(UF_long));
    if (pinv_probe != NULL) {
      for (k = 0; k < n; ++k) {
        pinv_probe[P[k]] = k;
      }
      for (k = 0; k < n; ++k) {
        const UF_long oldcol = Q[k];
        int found = 0;
        for (p = solver->col_ptr[oldcol]; p < solver->col_ptr[oldcol + 1u];
             ++p) {
          if (pinv_probe[solver->row_idx[p]] == k) {
            found = 1;
            break;
          }
        }
        if (!found) {
          missing++;
        }
      }
      free(pinv_probe);
      if (missing * 4u > (UF_long)n &&
          getenv("KLS_PREDICTED_TRY_SCALED") == NULL) {
        if (getenv("KLS_TRACE_PREDICTED") != NULL) {
          fprintf(stderr,
                  "KLS predicted: skipped, structural diagonal missing "
                  "%ld of %ld\n", (long)missing, (long)n);
        }
        return 0;
      }
    }
  }
  const double start = kls_now_seconds();
  /* The predicted factorization runs unscaled: the solve probe proves the
     result, and every strong-diagonal candidate that reaches this path
     settles on the unscaled mode in the auto-scale trials anyway.  Restore
     the mode on any rejection. */
  const UF_long saved_scale = solver->common.scale;
  if (getenv("KLS_PREDICTED_TRY_SCALED") == NULL) {
    solver->common.scale = -1;
  } else if (solver->common.scale <= 0) {
    /* scaled-static experiment: row scaling changes the addends of the
       exact cancellations that produce the zero pivots (missing-
       diagonal columns' pivots are pure update sums) */
    solver->common.scale = 1;
  }

  UF_long *pinv = (UF_long *)malloc((size_t)n * sizeof(*pinv));
  UF_long *col_block = generic_nd_fill_only
    ? NULL : (UF_long *)malloc((size_t)n * sizeof(*col_block));
  if (pinv == NULL || (!generic_nd_fill_only && col_block == NULL)) {
    solver->common.scale = saved_scale;
    free(pinv);
    free(col_block);
    return 0;
  }
  if (!generic_nd_fill_only) {
    for (UF_long block = 0; block < nblocks; ++block) {
      for (UF_long k = R[block]; k < R[block + 1u]; ++k) {
        col_block[k] = block;
      }
    }
  }
  for (UF_long k = 0; k < n; ++k) {
    pinv[P[k]] = k;
  }

  /* Block-local pattern in CSC over final (pivot) indices, plus the
     off-diagonal entry count. */
  UF_long *bcol_ptr = (UF_long *)malloc(((size_t)maxblock + 1u) *
                                        sizeof(*bcol_ptr));
  UF_long *work_rows = NULL;
  UF_long work_rows_capacity = 0;
  UF_long *parent = (UF_long *)malloc((size_t)maxblock * sizeof(*parent));
  UF_long *ancestor = (UF_long *)malloc((size_t)maxblock * sizeof(*ancestor));
  UF_long *mark = (UF_long *)malloc((size_t)maxblock * sizeof(*mark));
  UF_long *rowlen = (UF_long *)malloc((size_t)maxblock * sizeof(*rowlen));
  UF_long *colcount = (UF_long *)malloc((size_t)maxblock * sizeof(*colcount));
  UF_long *rowpat_ptr = NULL;
  UF_long *rowpat = NULL;
  UF_long *bcol_ptr_union = NULL;
  UF_long *work_rows_union = NULL;
  int ok = bcol_ptr != NULL && parent != NULL && ancestor != NULL &&
           mark != NULL && rowlen != NULL && colcount != NULL;

  trilinos_klu_l_numeric *numeric = NULL;
  int numeric_adopted = 0;      /* solver->numeric aliases the local */
  UF_long total_lnz = 0;
  UF_long total_unz = 0;
  UF_long max_lnz_block = 1;
  UF_long max_unz_block = 1;

  if (ok && !generic_nd_fill_only) {
    numeric = (trilinos_klu_l_numeric *)trilinos_klu_l_malloc(
      sizeof(*numeric), 1, common);
    ok = numeric != NULL && common->status >= TRILINOS_KLU_OK;
  }
  if (ok && !generic_nd_fill_only) {
    memset(numeric, 0, sizeof(*numeric));
    numeric->n = n;
    numeric->nblocks = nblocks;
    numeric->nzoff = symbolic->nzoff;
    numeric->Pnum =
      (UF_long *)trilinos_klu_l_malloc((size_t)n, sizeof(UF_long), common);
    numeric->Pinv =
      (UF_long *)trilinos_klu_l_malloc((size_t)n, sizeof(UF_long), common);
    numeric->Offp = (UF_long *)trilinos_klu_l_malloc((size_t)n + 1u,
                                                     sizeof(UF_long), common);
    numeric->Offi = (UF_long *)trilinos_klu_l_malloc(
      (size_t)symbolic->nzoff + 1u, sizeof(UF_long), common);
    numeric->Offx = trilinos_klu_l_malloc((size_t)symbolic->nzoff + 1u,
                                          sizeof(double), common);
    numeric->Lip =
      (UF_long *)trilinos_klu_l_malloc((size_t)n, sizeof(UF_long), common);
    numeric->Uip =
      (UF_long *)trilinos_klu_l_malloc((size_t)n, sizeof(UF_long), common);
    numeric->Llen =
      (UF_long *)trilinos_klu_l_malloc((size_t)n, sizeof(UF_long), common);
    numeric->Ulen =
      (UF_long *)trilinos_klu_l_malloc((size_t)n, sizeof(UF_long), common);
    numeric->LUsize = (size_t *)trilinos_klu_l_malloc((size_t)nblocks,
                                                      sizeof(size_t), common);
    numeric->LUbx = (void **)trilinos_klu_l_malloc((size_t)nblocks,
                                                   sizeof(void *), common);
    numeric->Udiag =
      trilinos_klu_l_malloc((size_t)n, sizeof(double), common);
    numeric->Rs = NULL;
    {
      size_t s = (size_t)n * sizeof(double);
      size_t n3 = (size_t)n * 3u * sizeof(double);
      size_t b6 = (size_t)maxblock * 6u * sizeof(UF_long);
      numeric->worksize = s + (n3 > b6 ? n3 : b6);
      numeric->Work = trilinos_klu_l_malloc(numeric->worksize, 1, common);
      numeric->Xwork = numeric->Work;
      numeric->Iwork =
        numeric->Work != NULL
          ? (UF_long *)((double *)numeric->Work + n)
          : NULL;
    }
    ok = common->status >= TRILINOS_KLU_OK && numeric->Pnum != NULL &&
         numeric->Pinv != NULL && numeric->Offp != NULL &&
         numeric->Offi != NULL && numeric->Offx != NULL &&
         numeric->Lip != NULL && numeric->Uip != NULL &&
         numeric->Llen != NULL && numeric->Ulen != NULL &&
         numeric->LUsize != NULL && numeric->LUbx != NULL &&
         numeric->Udiag != NULL && numeric->Work != NULL;
    if (!ok && getenv("KLS_TRACE_PREDICTED") != NULL) {
      fprintf(stderr, "KLS predicted: skeleton alloc failed status=%d\n",
              (int)common->status);
    }
    if (ok) {
      for (UF_long block = 0; block < nblocks; ++block) {
        numeric->LUbx[block] = NULL;
        numeric->LUsize[block] = 0;
      }
      memcpy(numeric->Pnum, P, (size_t)n * sizeof(*P));
      memcpy(numeric->Pinv, pinv, (size_t)n * sizeof(*pinv));
    }
  }

  /* Off-diagonal pattern in input-column order, matching the sequential
     Offx writes the refactorization performs. */
  if (ok && !generic_nd_fill_only) {
    UF_long poff = 0;
    for (UF_long k = 0; k < n && ok; ++k) {
      numeric->Offp[k] = poff;
      const UF_long oldcol = Q[k];
      const UF_long k1 = R[col_block[k]];
      for (UF_long p = solver->col_ptr[oldcol];
           p < solver->col_ptr[oldcol + 1u]; ++p) {
        const UF_long newrow = pinv[solver->row_idx[p]];
        if (newrow < k1) {
          if (poff >= symbolic->nzoff) {
            if (getenv("KLS_TRACE_PREDICTED") != NULL) {
              fprintf(stderr, "KLS predicted: off overflow at col %ld\n",
                      (long)k);
            }
            ok = 0;
            break;
          }
          numeric->Offi[poff++] = newrow;
        }
      }
    }
    if (ok && getenv("KLS_TRACE_PREDICTED") != NULL &&
        solver->block_trial_active) {
      fprintf(stderr, "KLS predicted: off section done\n");
    }
    if (ok) {
      numeric->Offp[n] = poff;
      ok = poff == symbolic->nzoff;
      if (!ok && getenv("KLS_TRACE_PREDICTED") != NULL) {
        fprintf(stderr, "KLS predicted: off count %ld != nzoff %ld\n",
                (long)poff, (long)symbolic->nzoff);
      }
    }
  }

  double kls_bp_union = 0.0, kls_bp_p1 = 0.0, kls_bp_lay = 0.0,
         kls_bp_p2 = 0.0;
  for (UF_long block = 0; block < nblocks && ok; ++block) {
    const UF_long k1 = R[block];
    const UF_long k2 = R[block + 1u];
    const UF_long nk = k2 - k1;
    double kls_bp_t = kls_now_seconds();
    if (solver->block_trial_active && getenv("KLS_TRACE_PREDICTED") != NULL &&
        nk >= 2) {
      fprintf(stderr, "KLS predicted: constructing block %ld nk=%ld\n",
              (long)block, (long)nk);
    }
    if (nk < 2) {
      if (!generic_nd_fill_only) {
        numeric->Lip[k1] = 0;
        numeric->Uip[k1] = 0;
        numeric->Llen[k1] = 0;
        numeric->Ulen[k1] = 0;
      }
      total_lnz += 1;
      total_unz += 1;
      continue;
    }
    /* Local pattern: count entries, then fill row-index lists per local
       column; verify a full structural diagonal.  Reject any block whose
       pattern is not structurally symmetric: only there does the symmetric
       elimination-tree reach give the exact no-pivot LU pattern. */
    UF_long bnz = 0;
    for (UF_long k = 0; k < nk; ++k) {
      const UF_long oldcol = Q[k1 + k];
      for (UF_long p = solver->col_ptr[oldcol];
           p < solver->col_ptr[oldcol + 1u]; ++p) {
        const UF_long newrow = pinv[solver->row_idx[p]];
        if (newrow >= k1 && newrow < k2) {
          bnz++;
        }
      }
    }
    if (bnz > work_rows_capacity) {
      free(work_rows);
      work_rows = (UF_long *)malloc((size_t)bnz * sizeof(*work_rows));
      work_rows_capacity = work_rows != NULL ? bnz : 0;
      if (work_rows == NULL) {
        ok = 0;
        break;
      }
    }
    memset(bcol_ptr, 0, ((size_t)nk + 1u) * sizeof(*bcol_ptr));
    for (UF_long k = 0; k < nk; ++k) {
      const UF_long oldcol = Q[k1 + k];
      for (UF_long p = solver->col_ptr[oldcol];
           p < solver->col_ptr[oldcol + 1u]; ++p) {
        const UF_long newrow = pinv[solver->row_idx[p]];
        if (newrow >= k1 && newrow < k2) {
          bcol_ptr[k + 1u]++;
        }
      }
    }
    for (UF_long k = 0; k < nk; ++k) {
      bcol_ptr[k + 1u] += bcol_ptr[k];
    }
    for (UF_long k = 0; k < nk; ++k) {
      UF_long cursor = bcol_ptr[k];
      const UF_long oldcol = Q[k1 + k];
      for (UF_long p = solver->col_ptr[oldcol];
           p < solver->col_ptr[oldcol + 1u]; ++p) {
        const UF_long newrow = pinv[solver->row_idx[p]];
        if (newrow >= k1 && newrow < k2) {
          work_rows[cursor++] = newrow - k1;
        }
      }
    }
    /* The union diagonal exists exactly where the block pattern has one:
       check it before paying for the symmetrized union at all. */
    {
      UF_long diag_present = 0;
      for (UF_long k = 0; k < nk; ++k) {
        for (UF_long p = bcol_ptr[k]; p < bcol_ptr[k + 1u]; ++p) {
          if (work_rows[p] == k) {
            diag_present++;
            break;
          }
        }
      }
      if (diag_present != nk && !solver->block_trial_active) {
        if (getenv("KLS_TRACE_PREDICTED") != NULL) {
          fprintf(stderr, "KLS predicted: diag precheck %ld/%ld\n",
                  (long)diag_present, (long)nk);
        }
        ok = 0;
        break;
      }
    }
    /* Symmetrize the block pattern: the symbolic Cholesky pattern of
       B | B' is the exact no-pivot LU bound (superset positions compute
       zeros), so mild asymmetry is absorbed instead of rejected.  Bail on
       strongly unsymmetric blocks where the bound is loose, and require a
       full structural diagonal. */
    {
      UF_long *ucol_ptr =
        (UF_long *)calloc((size_t)nk + 1u, sizeof(*ucol_ptr));
      UF_long *urows = (UF_long *)malloc((size_t)(2u * bnz + nk) *
                                         sizeof(*urows));
      if (ucol_ptr == NULL || urows == NULL) {
        free(ucol_ptr);
        free(urows);
        ok = 0;
        break;
      }
      for (UF_long k = 0; k < nk; ++k) {
        for (UF_long p = bcol_ptr[k]; p < bcol_ptr[k + 1u]; ++p) {
          const UF_long i = work_rows[p];
          ucol_ptr[k + 1u]++;
          if (i != k) {
            ucol_ptr[i + 1u]++;
          }
        }
      }
      for (UF_long k = 0; k < nk; ++k) {
        ucol_ptr[k + 1u] += ucol_ptr[k];
      }
      {
        UF_long *cursor = rowlen; /* reuse as scratch */
        memcpy(cursor, ucol_ptr, (size_t)nk * sizeof(*cursor));
        for (UF_long k = 0; k < nk; ++k) {
          for (UF_long p = bcol_ptr[k]; p < bcol_ptr[k + 1u]; ++p) {
            const UF_long i = work_rows[p];
            urows[cursor[k]++] = i;
            if (i != k) {
              urows[cursor[i]++] = k;
            }
          }
        }
      }
      /* sort + dedupe each column of the union */
      UF_long out = 0;
      UF_long diag_seen = 0;
      const int force_diag = solver->block_trial_active;
      for (UF_long k = 0; k < nk; ++k) {
        const UF_long lo = ucol_ptr[k];
        const UF_long hi = rowlen[k];
        if (hi - lo > 48) {
          /* Insertion sort is quadratic; spike columns in the union reach
             hundreds of thousands of entries. */
          qsort(urows + lo, (size_t)(hi - lo), sizeof(*urows),
                kls_compare_uf_long_asc);
        } else {
        for (UF_long a = lo + 1u; a < hi; ++a) {
          const UF_long v = urows[a];
          UF_long b = a;
          while (b > lo && urows[b - 1u] > v) {
            urows[b] = urows[b - 1u];
            b--;
          }
          urows[b] = v;
        }
        }
        const UF_long dst = out;
        UF_long prev = KLS_KLU_EMPTY;
        int have_diag = 0;
        for (UF_long p = lo; p < hi; ++p) {
          if (urows[p] != prev) {
            if (force_diag && !have_diag && urows[p] > k) {
              urows[out++] = k;
              have_diag = 1;
            }
            prev = urows[p];
            urows[out++] = prev;
            if (prev == k) {
              diag_seen++;
              have_diag = 1;
            }
          }
        }
        if (force_diag && !have_diag) {
          urows[out++] = k;
        }
        ucol_ptr[k] = dst;
      }
      ucol_ptr[nk] = out;
      if (solver->block_trial_active) {
        /* The trial's symmetric relabel has no matched diagonal; the
           pattern carries it explicitly (a bound of A+A'+I) and the
           pivoting fill provides the numeric pivots. */
        diag_seen = nk;
      }
      if (diag_seen != nk || out > 2u * bnz + nk) {
        if (getenv("KLS_TRACE_PREDICTED") != NULL) {
          fprintf(stderr,
                  "KLS predicted: block reject nk=%ld diag=%ld/%ld union=%ld"
                  " bound=%ld\n",
                  (long)nk, (long)diag_seen, (long)nk, (long)out,
                  (long)(2u * bnz + nk));
        }
        free(ucol_ptr);
        free(urows);
        ok = 0;
        break;
      }
      free(bcol_ptr_union);
      free(work_rows_union);
      bcol_ptr_union = ucol_ptr;
      work_rows_union = urows;
    }
    /* Elimination tree over the sorted symmetric block pattern. */
    for (UF_long k = 0; k < nk; ++k) {
      parent[k] = KLS_KLU_EMPTY;
      ancestor[k] = KLS_KLU_EMPTY;
    }
    for (UF_long k = 0; k < nk; ++k) {
      for (UF_long p = bcol_ptr_union[k]; p < bcol_ptr_union[k + 1u]; ++p) {
        UF_long i = work_rows_union[p];
        while (i != KLS_KLU_EMPTY && i < k) {
          const UF_long next = ancestor[i];
          ancestor[i] = k;
          if (next == KLS_KLU_EMPTY) {
            parent[i] = k;
            i = KLS_KLU_EMPTY;
          } else {
            i = next;
          }
        }
      }
    }
    kls_bp_union += kls_now_seconds() - kls_bp_t;
    kls_bp_t = kls_now_seconds();
    /* Pass 1: row reach counts give Ulen (row pattern length) and column
       counts give Llen.  Large blocks run the byte-identical parallel
       build (reach walks and the pass-2 L-scatter are DRAM-latency-bound
       and column-independent). */
    UF_long block_lnz = 0;
    int pb_used = 0;
    kls_pred_build_ctx pb_ctx;
    kls_pred_build_job pb_jobs[16];
    int32_t *pb_partial[16] = {NULL};
    UF_long *pb_lcur[16] = {NULL};
    UF_long pb_tlnz[16] = {0};
    int pb_nt = solver->options.threads > 16 ? 16 : solver->options.threads;
    if (nk >= 65536 && pb_nt >= 2 &&
        getenv("KLS_DISABLE_PARALLEL_PREDICTED_BUILD") == NULL) {
      int alloc_ok = 1;
      for (int t = 0; t < pb_nt; ++t) {
        pb_partial[t] = (int32_t *)calloc((size_t)nk, sizeof(int32_t));
        if (pb_partial[t] == NULL) {
          alloc_ok = 0;
          break;
        }
      }
      if (alloc_ok) {
        memset(&pb_ctx, 0, sizeof(pb_ctx));
        pb_ctx.nk = nk;
        pb_ctx.k1 = k1;
        pb_ctx.ucol_ptr = bcol_ptr_union;
        pb_ctx.urows = work_rows_union;
        pb_ctx.parent = parent;
        pb_ctx.rowlen = rowlen;
        pb_ctx.count_partial = pb_partial;
        pb_ctx.thread_lnz = pb_tlnz;
        pb_ctx.nthreads = pb_nt;
        pb_ctx.phase = 1;
        pb_ctx.failed = 0;
        for (int t = 0; t < pb_nt; ++t) {
          pb_jobs[t].ctx = &pb_ctx;
          pb_jobs[t].lo = nk * (UF_long)t / pb_nt;
          pb_jobs[t].hi = nk * (UF_long)(t + 1) / pb_nt;
          pb_jobs[t].tid = t;
        }
        if (kls_pred_build_run(&pb_ctx, pb_jobs)) {
          memset(colcount, 0, (size_t)nk * sizeof(*colcount));
          for (int t = 0; t < pb_nt; ++t) {
            const int32_t *part = pb_partial[t];
            for (UF_long i = 0; i < nk; ++i) {
              colcount[i] += part[i];
            }
            block_lnz += pb_tlnz[t];
          }
          pb_used = 1;
        }
      }
      if (!pb_used) {
        for (int t = 0; t < pb_nt; ++t) {
          free(pb_partial[t]);
          pb_partial[t] = NULL;
        }
      }
    }
    if (!pb_used) {
    memset(colcount, 0, (size_t)nk * sizeof(*colcount));
    for (UF_long k = 0; k < nk; ++k) {
      mark[k] = -1;
    }
    for (UF_long k = 0; k < nk; ++k) {
      mark[k] = k;
      UF_long len = 0;
      for (UF_long p = bcol_ptr_union[k]; p < bcol_ptr_union[k + 1u]; ++p) {
        UF_long i = work_rows_union[p];
        while (i < k && mark[i] != k) {
          mark[i] = k;
          colcount[i]++;
          len++;
          i = parent[i];
          if (i == KLS_KLU_EMPTY) {
            break;
          }
        }
      }
      rowlen[k] = len;
      block_lnz += len;
    }
    }
    if (solver->block_trial_active && getenv("KLS_TRACE_PREDICTED") != NULL) {
      fprintf(stderr, "KLS predicted: block %ld exact lnz=%ld\n",
              (long)block, (long)block_lnz);
    }
    if (solver->block_trial_active &&
        block_lnz > 16 * (UF_long)solver->nnz) {
      /* The blocked ordering only pays when the bound stays near the
         input size; a blowup here means the structure was misjudged. */
      ok = 0;
      break;
    }
    kls_bp_p1 += kls_now_seconds() - kls_bp_t;
    kls_bp_t = kls_now_seconds();
    if (generic_nd_fill_only) {
      UF_long block_unz = 0;
      for (UF_long k = 0; k < nk; ++k) {
        block_unz += rowlen[k];
      }
      for (int t = 0; t < pb_nt; ++t) {
        free(pb_partial[t]);
        pb_partial[t] = NULL;
      }
      total_lnz += block_lnz + nk;
      total_unz += block_unz + nk;
      if (block_lnz + nk > max_lnz_block) {
        max_lnz_block = block_lnz + nk;
      }
      if (block_unz + nk > max_unz_block) {
        max_unz_block = block_unz + nk;
      }
      continue;
    }
    /* Layout and allocation for this block. */
    size_t lusize = 0;
    for (UF_long k = 0; k < nk; ++k) {
      numeric->Llen[k1 + k] = colcount[k];
      numeric->Ulen[k1 + k] = rowlen[k];
      numeric->Lip[k1 + k] = (UF_long)lusize;
      lusize += 2u * (size_t)colcount[k];
      numeric->Uip[k1 + k] = (UF_long)lusize;
      lusize += 2u * (size_t)rowlen[k];
    }
    if (lusize == 0) {
      lusize = 1;
    }
    numeric->LUbx[block] =
      trilinos_klu_l_malloc(lusize, sizeof(double), common);
    numeric->LUsize[block] = lusize;
    if (numeric->LUbx[block] == NULL || common->status < TRILINOS_KLU_OK) {
      ok = 0;
      break;
    }
    kls_bp_lay += kls_now_seconds() - kls_bp_t;
    kls_bp_t = kls_now_seconds();
    /* Pass 2: fill sorted index areas.  U column k receives the sorted row
       reach of k directly; L columns receive k appended in ascending order
       as rows are visited. */
    {
      double *lu = (double *)numeric->LUbx[block];
      if (rowpat_ptr == NULL) {
        rowpat_ptr = (UF_long *)malloc((size_t)maxblock * sizeof(*rowpat_ptr));
        if (rowpat_ptr == NULL) {
          ok = 0;
          break;
        }
      }
      /* per-column L write cursors reuse colcount as cursor after reset */
      if (!pb_used) {
      for (UF_long k = 0; k < nk; ++k) {
        rowpat_ptr[k] = 0;
        mark[k] = -1;
      }
      }
      if (rowpat == NULL || (UF_long)block_lnz > 0) {
        free(rowpat);
        rowpat = (UF_long *)malloc(
          (size_t)(block_lnz > 0 ? block_lnz : 1) * sizeof(*rowpat));
        if (rowpat == NULL) {
          ok = 0;
          break;
        }
      }
      (void)rowpat;
      if (pb_used) {
        /* exclusive per-thread cursor bases from the pass-1 partials:
           thread ranges ascend, so per-column appends keep the exact
           serial ascending-k order */
        int pb2_ok = 1;
        for (int t = 0; t < pb_nt && pb2_ok; ++t) {
          pb_lcur[t] = (UF_long *)malloc((size_t)nk * sizeof(UF_long));
          if (pb_lcur[t] == NULL) {
            pb2_ok = 0;
          }
        }
        if (pb2_ok) {
          for (UF_long i = 0; i < nk; ++i) {
            pb_lcur[0][i] = 0;
          }
          for (int t = 1; t < pb_nt; ++t) {
            const int32_t *part = pb_partial[t - 1];
            const UF_long *prevc = pb_lcur[t - 1];
            UF_long *cur = pb_lcur[t];
            for (UF_long i = 0; i < nk; ++i) {
              cur[i] = prevc[i] + (UF_long)part[i];
            }
          }
          pb_ctx.lu = lu;
          pb_ctx.numeric = numeric;
          pb_ctx.lcursor = pb_lcur;
          pb_ctx.phase = 2;
          pb_ctx.failed = 0;
          pb2_ok = kls_pred_build_run(&pb_ctx, pb_jobs);
        }
        for (int t = 0; t < pb_nt; ++t) {
          free(pb_partial[t]);
          free(pb_lcur[t]);
          pb_partial[t] = NULL;
          pb_lcur[t] = NULL;
        }
        if (!pb2_ok) {
          ok = 0;
          break;
        }
      } else
      for (UF_long k = 0; k < nk; ++k) {
        mark[k] = k;
        UF_long len = 0;
        UF_long *urow = rowpat;
        for (UF_long p = bcol_ptr_union[k]; p < bcol_ptr_union[k + 1u]; ++p) {
          UF_long i = work_rows_union[p];
          while (i < k && mark[i] != k) {
            mark[i] = k;
            urow[len++] = i;
            i = parent[i];
            if (i == KLS_KLU_EMPTY) {
              break;
            }
          }
        }
        /* sort the row pattern ascending */
        for (UF_long a = 1u; a < len; ++a) {
          const UF_long v = urow[a];
          UF_long b = a;
          while (b > 0u && urow[b - 1u] > v) {
            urow[b] = urow[b - 1u];
            b--;
          }
          urow[b] = v;
        }
        UF_long *ui = (UF_long *)(lu + numeric->Uip[k1 + k]);
        for (UF_long p = 0; p < len; ++p) {
          ui[p] = urow[p];
          const UF_long j = urow[p];
          UF_long *li = (UF_long *)(lu + numeric->Lip[k1 + j]);
          li[rowpat_ptr[j]++] = k;
        }
      }
      for (UF_long k = 0; !pb_used && k < nk && ok; ++k) {
        if (rowpat_ptr[k] != (UF_long)numeric->Llen[k1 + k]) {
          if (getenv("KLS_TRACE_PREDICTED") != NULL) {
            fprintf(stderr, "KLS predicted: rowpat mismatch col %ld\n",
                    (long)(k1 + k));
          }
          ok = 0;
        }
      }
    }
    kls_bp_p2 += kls_now_seconds() - kls_bp_t;
    total_lnz += block_lnz + nk;
    {
      UF_long block_unz = 0;
      for (UF_long k = 0; k < nk; ++k) {
        block_unz += numeric->Ulen[k1 + k];
      }
      total_unz += block_unz + nk;
      if (block_lnz + nk > max_lnz_block) {
        max_lnz_block = block_lnz + nk;
      }
      if (block_unz + nk > max_unz_block) {
        max_unz_block = block_unz + nk;
      }
    }
  }

  if (getenv("KLS_TRACE_PREDICTED") != NULL) {
    fprintf(stderr,
            "KLS predicted: build phases union=%.2f p1=%.2f layout=%.2f"
            " p2=%.2f\n", kls_bp_union, kls_bp_p1, kls_bp_lay, kls_bp_p2);
  }
  free(work_rows);
  free(bcol_ptr);
  free(parent);
  free(ancestor);
  free(mark);
  free(rowlen);
  free(colcount);
  free(rowpat_ptr);
  free(rowpat);
  free(bcol_ptr_union);
  free(work_rows_union);
  free(col_block);
  free(pinv);

  if (ok) {
    const double est_total = symbolic->lnz + symbolic->unz;
    const double exact_total = (double)total_lnz + (double)total_unz;
    const int provisional_nd = solver->generic_nd_portfolio_selected;
    const double fallback_factor =
      solver->generic_nd_lifecycle_near_tie ? 1.05 : 0.95;
    const int generic_nd_fill_reject =
      provisional_nd &&
      ((est_total > 0.0 &&
        exact_total > 1.25 * est_total + (double)n) ||
       (solver->generic_nd_fallback_fill > 0.0 &&
        exact_total > fallback_factor * solver->generic_nd_fallback_fill));
    if ((est_total > 0.0 &&
         exact_total > 2.0 * est_total + (double)n) ||
        generic_nd_fill_reject) {
      if (getenv("KLS_TRACE_PREDICTED") != NULL) {
        fprintf(stderr,
                "KLS predicted: fill reject exact=%.3e est=%.3e "
                "fallback=%.3e provisional_nd=%d\n",
                exact_total, est_total, solver->generic_nd_fallback_fill,
                provisional_nd);
      }
      ok = 0;
    }
  }
  if (!ok) {
    if (getenv("KLS_TRACE_PREDICTED") != NULL) {
      fprintf(stderr, "KLS predicted: construction rejected\n");
    }
    solver->common.scale = saved_scale;
    if (numeric_adopted) {
      /* the local aliases solver->numeric; failure paths inside the
         value machinery may have freed it already (free_numeric NULLs
         the solver field) - never re-free through the stale alias */
      if (solver->numeric != NULL) {
        free_numeric(solver);
      }
      numeric = NULL;
    }
    if (numeric != NULL) {
      if (getenv("KLS_TRACE_PREDICTED") != NULL) {
        fprintf(stderr,
                "KLS predicted reject dump: n=%ld nzoff=%ld nblocks=%ld"
                " Pnum=%p Offp=%p Offi=%p Offx=%p Lip=%p LUsize=%p"
                " LUbx=%p Udiag=%p Rs=%p Pinv=%p Work=%p\n",
                (long)numeric->n, (long)numeric->nzoff,
                (long)numeric->nblocks, (void *)numeric->Pnum,
                (void *)numeric->Offp, (void *)numeric->Offi,
                (void *)numeric->Offx, (void *)numeric->Lip,
                (void *)numeric->LUsize, (void *)numeric->LUbx,
                (void *)numeric->Udiag, (void *)numeric->Rs,
                (void *)numeric->Pinv, (void *)numeric->Work);
        for (UF_long b2 = 0; b2 < numeric->nblocks && b2 < 8; ++b2) {
          fprintf(stderr, "KLS predicted reject blk %ld LUbx=%p size=%zu\n",
                  (long)b2, numeric->LUbx ? numeric->LUbx[b2] : NULL,
                  numeric->LUsize ? numeric->LUsize[b2] : 0u);
        }
      }
      trilinos_klu_l_free_numeric(&numeric, common);
    }
    *elapsed += kls_now_seconds() - start;
    return 0;
  }

  if (solver->generic_nd_portfolio_selected &&
      solver->generic_nd_bounded_symmetric_union &&
      !(symbolic->lnz + symbolic->unz > 0.0)) {
    /* KLU's user-ordering analyze cannot publish a fill estimate.  The
       predictor has now constructed the exact static structural pattern,
       even if its value pass later rejects the fixed diagonal.  Retain that
       pattern count on the provisional symbolic so ordinary pivoted KLU may
       validate a decisively smaller ND proposal under the existing realized
       fill rollback instead of discarding the ordering merely because the
       static numeric was unsuitable. */
    symbolic->lnz = (double)total_lnz;
    symbolic->unz = (double)total_unz;
  }

  if (generic_nd_fill_only) {
    /* AUTO already selected a robust scaled numeric for this provisional
       user ordering.  The static predictor deliberately runs unscaled, so
       building its numeric, value streams, and residual probe cannot certify
       the representation we will retain.  The exact symmetric-union reach
       above supplies the missing fill bound; now let ordinary pivoted KLU
       validate the scaled numeric once under the realized-fill rollback. */
    solver->common.scale = saved_scale;
    if (getenv("KLS_TRACE_PREDICTED") != NULL) {
      fprintf(stderr,
              "KLS predicted: fill-only certificate %.0f+%.0f\n",
              symbolic->lnz, symbolic->unz);
    }
    *elapsed += kls_now_seconds() - start;
    return 0;
  }

  numeric->lnz = total_lnz;
  numeric->unz = total_unz;
  numeric->max_lnz_block = max_lnz_block;
  numeric->max_unz_block = max_unz_block;

  /* Fill the first values through the solver's own refactorization
     machinery: the pattern is already ascending so the supernode panel
     preparation and retained schedules apply directly, and the threaded
     panel refactor computes the identical values one to two orders of
     magnitude faster than the serial KLU refactorization. */
  common->status = TRILINOS_KLU_OK;
  solver->numeric = numeric;
  numeric_adopted = 1;
  solver->numeric_is_predicted = 1;
  /* The pattern was written in ascending row order, so the supernode
     preparation can skip its sorting pass. */
  solver->snode_numeric_pre_sorted =
    getenv("KLS_PREDICTED_FORCE_SORT") == NULL;
  /* The routing gates below need the factor work estimate, which KLU
     computes from the pattern alone. */
  kls_update_numeric_diagnostics(solver, 0);
  const double build_done = kls_now_seconds();
  kls_maybe_prepare_snode_panels(solver, elapsed);
  kls_snb_maybe_accept(solver, numeric_values, elapsed);
  (void)kls_i32_solve_ready(solver);
  const double snode_done = kls_now_seconds();
  maybe_prepare_refactor_map(solver, elapsed);
  const double map_done = kls_now_seconds();
  maybe_prepare_refactor_schedule(solver, elapsed);
  const double sched_done = kls_now_seconds();
  if (getenv("KLS_TRACE_PREDICTED") != NULL) {
    fprintf(stderr,
            "KLS predicted: preps snode=%.2fs map=%.2fs sched=%.2fs"
            " fill=%ld+%ld\n",
            snode_done - build_done, map_done - snode_done,
            sched_done - map_done, (long)total_lnz, (long)total_unz);
  }
  UF_long refactor_ok = 0;
  int saw_singular = 0;
  if (getenv("KLS_FORCE_PIVOT_FILL") != NULL || solver->block_trial_active) {
    /* The block trial's matched diagonal is structurally incomplete by
       design; the pivoting fill is its intended factorization. */
    refactor_ok = 0;
    common->status = TRILINOS_KLU_SINGULAR;
  }
  UF_long fill_round_cap = 3u;
  const UF_long nudge_slot_cap =
    getenv("KLS_PREDICTED_NUDGE_SLOTS_8K") != NULL ? 8192u : 1024u;
  const int batch_fast_rounds =
    getenv("KLS_PREDICTED_BATCH_FAST") != NULL;
  UF_long batch_window = 32u;
  {
    const char *we = getenv("KLS_PREDICTED_BATCH_WINDOW");
    if (we != NULL && we[0] != '\0') {
      const long parsed = atol(we);
      if (parsed > 0) {
        batch_window = (UF_long)parsed;
      }
    }
  }
  {
    const char *rounds_env = getenv("KLS_PREDICTED_NUDGE_ROUNDS");
    if (rounds_env != NULL && rounds_env[0] != '\0') {
      const long parsed = atol(rounds_env);
      if (parsed > 0 && parsed <= 2048) {
        fill_round_cap = (UF_long)parsed;
      }
    }
  }
  const int zp_collect_mode =
    getenv("KLS_PREDICTED_COLLECT_ZEROS") != NULL &&
    !solver->block_trial_active;
  const Int zp_saved_halt = common->halt_if_singular;
  if (zp_collect_mode) {
    /* Batch zero-pivot collection: run each fill pass with the halt
       off so every worker records its exact-zero pivots (downstream
       Inf/NaN never compares equal to 0.0, so one pass collects the
       frontier of independent zeros; mac_econ: 507 one-per-pass
       rounds collapse to the dependency depth). */
    solver->zero_pivot_collect_cap = 4096;
    solver->zero_pivot_collect = (UF_long *)malloc(
      (size_t)solver->zero_pivot_collect_cap *
      sizeof(*solver->zero_pivot_collect));
    if (solver->zero_pivot_collect != NULL) {
      common->halt_if_singular = 0;
    }
  }
  UF_long kls_prev_batch_nudged = 0;
  if (refactor_ok == 0 &&
      common->status == TRILINOS_KLU_SINGULAR &&
      (getenv("KLS_FORCE_PIVOT_FILL") != NULL ||
       solver->block_trial_active)) {
    /* forced-fill path set above: skip the value rounds */
  } else
  for (UF_long fill_round = 0; fill_round < fill_round_cap; ++fill_round) {
    const double first_values_t0 = kls_now_seconds();
    if (solver->zero_pivot_collect != NULL) {
      atomic_store_explicit(&solver->zero_pivot_collect_count, 0,
                            memory_order_relaxed);
    }
    refactor_ok = kls_parallel_refactor(solver, numeric_values, 0);
    if (solver->zero_pivot_collect != NULL &&
        atomic_load_explicit(&solver->zero_pivot_collect_count,
                             memory_order_relaxed) > 0) {
      /* halt-off passes complete and their join reports neither a
         failed ok nor a singular status; the recorded zeros are the
         verdict */
      refactor_ok = 0;
      common->status = TRILINOS_KLU_SINGULAR;
    }
    if (getenv("KLS_TRACE_PREDICTED") != NULL) {
      fprintf(stderr,
              "KLS predicted: value pass %.2fs path=%d ok=%ld\n",
              kls_now_seconds() - first_values_t0,
              (int)solver->stats.last_refactor_path, (long)refactor_ok);
    }
    UF_long k = KLS_KLU_EMPTY;
    if (refactor_ok && common->status >= TRILINOS_KLU_OK &&
        common->status != TRILINOS_KLU_SINGULAR) {
      UF_long threshold_pivot = KLS_KLU_EMPTY;
      UF_long threshold_col = KLS_KLU_EMPTY;
      if (kls_numeric_pivots_pass_threshold(solver, &threshold_pivot,
                                            &threshold_col)) {
        break;
      }
      /* Pivots too small for the 1/tol multiplier bound: nudge every
         violating column at once, each by the amount its own multipliers
         require, so the round count depends on escalation depth rather
         than on how many weak pivots the matrix has.  The solve probe
         stays the accuracy arbiter. */
      const double tolv = common->tol > 0.0 ? common->tol : 0.001;
      UF_long batch_nudged = 0;
      int batch_failed = 0;
      if (getenv("KLS_PREDICTED_THRESHOLD_WAIVE") != NULL ||
          kls_predicted_probe_first_threshold_shape(solver)) {
        if (getenv("KLS_TRACE_PREDICTED") != NULL) {
          fprintf(stderr,
                  "KLS predicted: threshold deferred, probe arbitrates\n");
        }
        break;
      }
      /* Front-first: violations cascade downstream, so fix only the
         leading window each round; later columns are re-judged on sane
         upstream values in the next round. */
      for (UF_long block = 0;
           block < nblocks && !batch_failed && batch_nudged < batch_window;
           ++block) {
        const UF_long bk1 = R[block];
        const UF_long bk2 = R[block + 1u];
        if (bk2 - bk1 <= 1u) {
          continue;
        }
        double *lu = (double *)numeric->LUbx[block];
        const UF_long *lip = numeric->Lip + bk1;
        const UF_long *llen = numeric->Llen + bk1;
        const double *udiag = (const double *)numeric->Udiag;
        for (UF_long kk = 0; kk < bk2 - bk1 && batch_nudged < batch_window;
             ++kk) {
          UF_long *li = NULL;
          double *lx = NULL;
          UF_long lcol_len = 0;
          kls_klu_get_pointer(lu, lip, llen, kk, &li, &lx, &lcol_len);
          (void)li;
          double maxmult = 0.0;
          int bad = 0;
          for (UF_long p = 0; p < lcol_len; ++p) {
            const double v = fabs(lx[p]);
            if (!isfinite(v)) {
              bad = 1;
              break;
            }
            maxmult = maxmult < v ? v : maxmult;
          }
          if (!bad && !(maxmult * tolv > 1.0 + 1.0e-12)) {
            continue;
          }
          const UF_long gk = bk1 + kk;
          const UF_long goldcol = Q[gk];
          UF_long gdiag_pos = KLS_KLU_EMPTY;
          double gcolmax = 0.0;
          for (UF_long p = solver->col_ptr[goldcol];
               p < solver->col_ptr[goldcol + 1u]; ++p) {
            const double av = fabs(numeric_values[p]);
            gdiag_pos = numeric->Pinv[solver->row_idx[p]] == gk
              ? p : gdiag_pos;
            gcolmax = gcolmax < av ? av : gcolmax;
          }
          const double gpivot = fabs(udiag[gk]);
          double sigma = bad || gpivot == 0.0
            ? tolv * (gcolmax > 0.0 ? gcolmax : 1.0)
            : 1.25 * tolv * maxmult * gpivot;
          if (gdiag_pos == KLS_KLU_EMPTY) {
            /* Diagonal present only through symmetrization: nothing in A to
               nudge.  Leave the column to the solve probe. */
            continue;
          }
          if (solver->pivot_nudge_count + batch_nudged >=
              (solver->pivot_nudge_capacity > 0
                 ? solver->pivot_nudge_capacity : nudge_slot_cap)) {
            continue;
          }
          if (solver->pivot_nudge_values == NULL) {
            solver->pivot_nudge_values = (double *)malloc(
              (size_t)solver->nnz * sizeof(*solver->pivot_nudge_values));
            solver->pivot_nudge_pos = (UF_long *)malloc(
              nudge_slot_cap * sizeof(*solver->pivot_nudge_pos));
            solver->pivot_nudge_sigma = (double *)malloc(
              nudge_slot_cap * sizeof(*solver->pivot_nudge_sigma));
            if (solver->pivot_nudge_values == NULL ||
                solver->pivot_nudge_pos == NULL ||
                solver->pivot_nudge_sigma == NULL) {
              free_pivot_nudges(solver);
              batch_failed = 1;
              break;
            }
            solver->pivot_nudge_capacity = nudge_slot_cap;
          }
          UF_long slot = KLS_KLU_EMPTY;
          for (UF_long i = 0; i < solver->pivot_nudge_count; ++i) {
            if (solver->pivot_nudge_pos[i] == gdiag_pos) {
              slot = i;
              break;
            }
          }
          if (slot != KLS_KLU_EMPTY) {
            const double cap = gcolmax > 0.0 ? gcolmax : 1.0;
            double escalated = solver->pivot_nudge_sigma[slot] * 32.0;
            escalated = sigma > escalated ? sigma : escalated;
            if (!(fabs(escalated) <= cap)) {
              /* Escalation exhausted: keep the capped correction and let
                 the refining solve probe judge the numeric as a whole. */
              continue;
            }
            solver->pivot_nudge_sigma[slot] = escalated;
          } else {
            solver->pivot_nudge_pos[solver->pivot_nudge_count] = gdiag_pos;
            solver->pivot_nudge_sigma[solver->pivot_nudge_count] = sigma;
            solver->pivot_nudge_count++;
          }
          batch_nudged++;
        }
      }
      if (batch_failed) {
        refactor_ok = 0;
        break;
      }
      if (batch_nudged == 0) {
        if (getenv("KLS_TRACE_PREDICTED") != NULL) {
          fprintf(stderr,
                  "KLS predicted: threshold waived, probe arbitrates\n");
        }
        break;
      }
      if (getenv("KLS_TRACE_PREDICTED") != NULL) {
        fprintf(stderr,
                "KLS predicted: batch nudged %ld columns (total %ld)"
                " round=%ld\n",
                (long)batch_nudged, (long)solver->pivot_nudge_count,
                (long)fill_round);
      }
      if (fill_round >= 1u &&
          (UF_long)batch_nudged >= kls_prev_batch_nudged) {
        /* non-converging nudges: each pass exposes as many new zero
           pivots as the last (rajat25: 32 -> 30 -> 32 then a fill
           reject anyway) - the class pivots numerically and the
           diagonal replay cannot hold; reject now instead of paying
           the remaining rounds and the fill check */
        if (getenv("KLS_TRACE_PREDICTED") != NULL) {
          fprintf(stderr, "KLS predicted: nudge divergence reject\n");
        }
        refactor_ok = 0;
        break;
      }
      kls_prev_batch_nudged = (UF_long)batch_nudged;
      if (!batch_fast_rounds) {
        /* nudges change values, not pattern; the rebuild costs ~14s per
           round at pre2 scale and is skipped under
           KLS_PREDICTED_BATCH_FAST */
        free_refactor_map(solver);
        free_refactor_schedule(solver);
        maybe_prepare_refactor_map(solver, elapsed);
        maybe_prepare_refactor_schedule(solver, elapsed);
      }
      common->status = TRILINOS_KLU_OK;
      refactor_ok = 0;
      continue;
    } else {
      /* Zero pivot: rescue it with a tiny explicit diagonal correction and
         refill.  The correction persists with the numeric so every later
         refactorization factors the same nudged matrix, and the solve probe
         below still measures the residual against the original values. */
      if (common->status != TRILINOS_KLU_SINGULAR) {
        if (getenv("KLS_TRACE_PREDICTED") != NULL) {
          fprintf(stderr,
                  "KLS predicted: pass failed status=%d path=%d "
                  "round=%ld (not singular)\n",
                  (int)common->status,
                  (int)solver->stats.last_refactor_path, (long)fill_round);
        }
        break;
      }
      saw_singular = 1;
      k = (UF_long)common->numerical_rank;
      if (k == KLS_KLU_EMPTY || k >= n) {
        /* not every refactor path records numerical_rank; the singular
           column is the same failing pivot */
        k = (UF_long)common->singular_col;
      }
    }
    long zp_count = 0;
    if (solver->zero_pivot_collect != NULL) {
      zp_count = atomic_load_explicit(&solver->zero_pivot_collect_count,
                                      memory_order_relaxed);
      if (zp_count > solver->zero_pivot_collect_cap) {
        zp_count = solver->zero_pivot_collect_cap;
      }
      if (zp_count > 1 && getenv("KLS_TRACE_PREDICTED") != NULL) {
        fprintf(stderr, "KLS predicted: zp batch %ld zeros round=%ld\n",
                zp_count, (long)fill_round);
      }
    }
    int zp_fatal = 0;
    for (long zp_i = 0; zp_i < (zp_count > 0 ? zp_count : 1) && !zp_fatal;
         ++zp_i) {
    if (zp_count > 0) {
      k = solver->zero_pivot_collect[zp_i];
    }
    if (k == KLS_KLU_EMPTY || k >= n ||
        solver->pivot_nudge_count >= nudge_slot_cap) {
      zp_fatal = 1;
      break;
    }
    const UF_long oldcol = Q[k];
    UF_long diag_pos = KLS_KLU_EMPTY;
    double colmax = 0.0;
    for (UF_long p = solver->col_ptr[oldcol]; p < solver->col_ptr[oldcol + 1u];
         ++p) {
      const double av = fabs(numeric_values[p]);
      colmax = colmax < av ? av : colmax;
      if (numeric->Pinv[solver->row_idx[p]] == k) {
        diag_pos = p;
      }
    }
    if (diag_pos == KLS_KLU_EMPTY) {
      if (getenv("KLS_TRACE_PREDICTED") != NULL) {
        fprintf(stderr,
                "KLS predicted: zp diag missing k=%ld oldcol=%ld "
                "symP[k]=%ld Pnum[k]=%ld collen=%ld\n",
                (long)k, (long)oldcol, (long)symbolic->P[k],
                (long)numeric->Pnum[k],
                (long)(solver->col_ptr[oldcol + 1u] -
                       solver->col_ptr[oldcol]));
        for (UF_long p = solver->col_ptr[oldcol];
             p < solver->col_ptr[oldcol + 1u] &&
             p < solver->col_ptr[oldcol] + 8; ++p) {
          fprintf(stderr,
                  "KLS predicted:   row %ld pinv %ld val %.3e\n",
                  (long)solver->row_idx[p],
                  (long)numeric->Pinv[solver->row_idx[p]],
                  numeric_values[p]);
        }
      }
      zp_fatal = 1;
      break;
    }
    if (solver->pivot_nudge_values == NULL) {
      solver->pivot_nudge_values = (double *)malloc(
        (size_t)solver->nnz * sizeof(*solver->pivot_nudge_values));
      solver->pivot_nudge_pos = (UF_long *)malloc(
        nudge_slot_cap * sizeof(*solver->pivot_nudge_pos));
      solver->pivot_nudge_sigma = (double *)malloc(
        nudge_slot_cap * sizeof(*solver->pivot_nudge_sigma));
      if (solver->pivot_nudge_values == NULL ||
          solver->pivot_nudge_pos == NULL ||
          solver->pivot_nudge_sigma == NULL) {
        free_pivot_nudges(solver);
        zp_fatal = 1;
        break;
      }
      solver->pivot_nudge_capacity = nudge_slot_cap;
    }
    const double base = colmax > 0.0 ? colmax : 1.0;
    UF_long existing = KLS_KLU_EMPTY;
    for (UF_long i = 0; i < solver->pivot_nudge_count; ++i) {
      if (solver->pivot_nudge_pos[i] == diag_pos) {
        existing = i;
        break;
      }
    }
    if (existing != KLS_KLU_EMPTY) {
      solver->pivot_nudge_sigma[existing] *= 32.0;
      if (!(fabs(solver->pivot_nudge_sigma[existing]) <= base)) {
        if (getenv("KLS_TRACE_PREDICTED") != NULL) {
          fprintf(stderr,
                  "KLS predicted: zp escalation exhausted k=%ld "
                  "sigma=%.3e base=%.3e udiag=%.6e a_kk=%.6e\n",
                  (long)k, solver->pivot_nudge_sigma[existing], base,
                  ((const double *)numeric->Udiag)[k],
                  numeric_values[diag_pos]);
        }
        zp_fatal = 1;
        break;
      }
    } else {
      double tol_sigma =
        (common->tol > 0.0 ? common->tol : 0.001) * base;
      {
        /* Smaller sigma lowers the perturbation floor the probe's
           refinement converges to (it scales with |E|); growth in the
           unpivoted factor scales inversely.  pre2 stalls at 3.8e-8
           with the 1e-3 default. */
        const char *se = getenv("KLS_PREDICTED_NUDGE_SIGMA_SCALE");
        if (se != NULL && se[0] != '\0') {
          const double parsed = atof(se);
          if (parsed > 0.0 && parsed <= 1.0) {
            tol_sigma *= parsed;
          }
        }
      }
      solver->pivot_nudge_pos[solver->pivot_nudge_count] = diag_pos;
      solver->pivot_nudge_sigma[solver->pivot_nudge_count] = tol_sigma;
      solver->pivot_nudge_count++;
    }
    if (getenv("KLS_TRACE_PREDICTED") != NULL) {
      fprintf(stderr,
              "KLS predicted: nudge col=%ld sigma=%.3e round=%ld\n",
              (long)k,
              solver->pivot_nudge_sigma[existing != KLS_KLU_EMPTY
                                          ? existing
                                          : solver->pivot_nudge_count - 1u],
              (long)fill_round);
    }
    }
    if (zp_fatal) {
      break;
    }
    common->status = TRILINOS_KLU_OK;
    common->numerical_rank = KLS_KLU_EMPTY;
    common->singular_col = KLS_KLU_EMPTY;
    refactor_ok = 0;
  }
  common->halt_if_singular = zp_saved_halt;
  if (solver->zero_pivot_collect != NULL) {
    free(solver->zero_pivot_collect);
    solver->zero_pivot_collect = NULL;
    solver->zero_pivot_collect_cap = 0;
  }
  if (getenv("KLS_TRACE_PREDICTED") != NULL) {
    fprintf(stderr, "KLS predicted: build=%.2fs first_values=%.2fs\n",
            build_done - start, kls_now_seconds() - build_done);
  }
  int pivot_fill_used = 0;
  if (!refactor_ok || common->status < TRILINOS_KLU_OK ||
      common->status == TRILINOS_KLU_SINGULAR) {
    if (getenv("KLS_TRACE_PREDICTED") != NULL) {
      fprintf(stderr,
              "KLS predicted: fill reject ok=%ld status=%d singular_col=%ld\n",
              (long)refactor_ok, (int)common->status,
              (long)common->singular_col);
    }
    const int rescue_worthwhile =
      common->status == TRILINOS_KLU_SINGULAR || saw_singular;
    common->status = TRILINOS_KLU_OK;
    common->numerical_rank = KLS_KLU_EMPTY;
    common->singular_col = KLS_KLU_EMPTY;
    free_pivot_nudges(solver);
    if (!rescue_worthwhile ||
        !kls_predicted_pivot_rescue(solver, numeric_values, elapsed)) {
      if (getenv("KLS_TRACE_PREDICTED") != NULL) {
        fprintf(stderr, "KLS predicted: pivoting fill reject status=%d\n",
                (int)common->status);
      }
      free_numeric(solver);
      common->status = TRILINOS_KLU_OK;
      common->scale = saved_scale;
      *elapsed += kls_now_seconds() - start;
      return 0;
    }
    pivot_fill_used = 1;
    if (getenv("KLS_TRACE_PREDICTED") != NULL) {
      fprintf(stderr, "KLS predicted: pivoting fill accepted structurally\n");
    }
    free_refactor_map(solver);
    free_refactor_schedule(solver);
    maybe_prepare_refactor_map(solver, elapsed);
    maybe_prepare_refactor_schedule(solver, elapsed);
  }
  /* Reject numerically unstable no-pivot factorizations directly: solve
     A x = A e and require a small true residual.  Per-column multiplier
     bounds cannot stop growth from compounding across elimination levels
     (Freescale1 passes them yet yields a 5e-3 residual). */
  for (int probe_attempt = 0; probe_attempt < 2; ++probe_attempt) {
    double *probe_b =
      (double *)malloc(2u * (size_t)n * sizeof(*probe_b));
    if (probe_b == NULL) {
      free_numeric(solver);
      common->status = TRILINOS_KLU_OK;
      common->scale = saved_scale;
      *elapsed += kls_now_seconds() - start;
      return 0;
    }
    double *probe_x = probe_b + n;
    for (UF_long i = 0; i < n; ++i) {
      probe_b[i] = 0.0;
    }
    /* Varied test solution: probing A x = A e is degenerate on matrices
       whose rows sum to zero exactly (mac_econ accounting identities:
       b = 0, x = 0, residual 0 for ANY nonsingular factor). */
    for (UF_long j = 0; j < n; ++j) {
      const double xj = 1.0 + (double)(j % 17) * 0.01;
      for (UF_long p = solver->col_ptr[j]; p < solver->col_ptr[j + 1u];
           ++p) {
        probe_b[solver->row_idx[p]] += numeric_values[p] * xj;
      }
    }
    double *probe_r = (double *)malloc((size_t)n * sizeof(*probe_r));
    memcpy(probe_x, probe_b, (size_t)n * sizeof(*probe_x));
    double residual = HUGE_VAL;
    double raw_residual = HUGE_VAL;
    int refine_iters = 0;
    if (probe_r != NULL &&
        trilinos_klu_l_solve(symbolic, numeric, n, 1, probe_x, common) &&
        common->status >= TRILINOS_KLU_OK) {
      double bmax = 0.0;
      for (UF_long i = 0; i < n; ++i) {
        const double v = fabs(probe_b[i]);
        bmax = bmax < v ? v : bmax;
      }
      const double bscale = bmax > 0.0 ? bmax : 1.0;
      /* Refine while it helps: a reduced-accuracy factorization (nudged or
         waived pivots) that converges under refinement is still usable,
         with production solves flagged to refine the same way. */
      int probe_iter_cap = solver->block_trial_active ? 10 : 3;
      double probe_stall_factor = 0.5;
      {
        const char *pe = getenv("KLS_PREDICTED_PROBE_ITERS");
        if (pe != NULL && pe[0] != '\0') {
          const long parsed = atol(pe);
          if (parsed > 0 && parsed <= 30) {
            probe_iter_cap = (int)parsed;
            probe_stall_factor = 0.75;
          }
        }
      }
      for (int iter = 0; iter <= probe_iter_cap; ++iter) {
        memcpy(probe_r, probe_b, (size_t)n * sizeof(*probe_r));
        for (UF_long j = 0; j < n; ++j) {
          const double xv = probe_x[j];
          if (xv == 0.0) {
            continue;
          }
          for (UF_long p = solver->col_ptr[j]; p < solver->col_ptr[j + 1u];
               ++p) {
            probe_r[solver->row_idx[p]] -= numeric_values[p] * xv;
          }
        }
        double rmax = 0.0;
        int r_finite = 1;
        for (UF_long i = 0; i < n; ++i) {
          const double v = fabs(probe_r[i]);
          if (!isfinite(v)) {
            /* max-abs with NaN comparisons would silently keep 0.0 and
               accept a poisoned solve as a perfect residual */
            r_finite = 0;
            break;
          }
          rmax = rmax < v ? v : rmax;
        }
        const double rel = r_finite ? rmax / bscale : HUGE_VAL;
        if (iter == 0) {
          raw_residual = rel;
        }
        if (!(rel < residual * probe_stall_factor) ||
            rel < (solver->block_trial_active ? 3.0e-10 : 1.0e-13) ||
            iter == probe_iter_cap) {
          residual = rel < residual ? rel : residual;
          break;
        }
        residual = rel;
        if (rel < 1.0e-9 && iter > 0) {
          break;
        }
        if (!trilinos_klu_l_solve(symbolic, numeric, n, 1, probe_r, common) ||
            common->status < TRILINOS_KLU_OK) {
          break;
        }
        for (UF_long i = 0; i < n; ++i) {
          probe_x[i] += probe_r[i];
        }
        refine_iters = iter + 1;
      }
    }
    free(probe_r);
    free(probe_b);
    if (getenv("KLS_TRACE_PREDICTED") != NULL) {
      fprintf(stderr, "KLS predicted: probe_residual=%.3e refine_iters=%d\n",
              residual, refine_iters);
    }
    if (!(residual < 1e-9)) {
      if ((getenv("KLS_ENABLE_PIVOT_FILL") != NULL ||
           solver->block_trial_active) &&
          probe_attempt == 0 && !pivot_fill_used &&
          kls_predicted_pivot_rescue(solver, numeric_values, elapsed)) {
        pivot_fill_used = 1;
        if (getenv("KLS_TRACE_PREDICTED") != NULL) {
          fprintf(stderr,
                  "KLS predicted: probe retry after pivoting fill\n");
        }
        free_refactor_map(solver);
        free_refactor_schedule(solver);
        maybe_prepare_refactor_map(solver, elapsed);
        maybe_prepare_refactor_schedule(solver, elapsed);
        continue;
      }
      free_numeric(solver);
      common->status = TRILINOS_KLU_OK;
      common->scale = saved_scale;
      *elapsed += kls_now_seconds() - start;
      return 0;
    }
    if (refine_iters > 0 &&
        (!(raw_residual < 1.0e-10) ||
         getenv("KLS_STRICT_PREDICTED_REFINEMENT") != NULL)) {
      /* Only keep per-solve refinement when the raw factorization is
         genuinely inaccurate: a raw probe residual at or below 1e-10 is
         SPICE-grade on its own, and the refinement tax lands on every
         one of the ~200 cycle-charged solves. */
      solver->numeric_needs_refinement = 1;
      /* the probe converging in one or two sweeps proves a single
         correction reaches 1e-13-class: skip the post-correction
         verification SpMV on production solves */
      solver->solve_refine_single_shot =
        refine_iters <= 2 && getenv("KLS_STRICT_PREDICTED_REFINEMENT") == NULL;
    }
    break;
  }
  solver->auto_scale_checked = 1;
  *elapsed += kls_now_seconds() - start;
  return 1;
}

/* Serial pivoting fill for a predicted-pattern numeric: one left-looking
   pass per block with threshold partial pivoting confined to the fixed
   symbolic pattern.  A row swap relabels two not-yet-pivotal pattern slots:
   the workspace entries, the pending input scatter (through cur_slot), and
   the already-stored L values at those slots in earlier columns (found
   through the pattern's row lists) all swap, exactly as in dense getrf.
   The solve probe remains the arbiter for the rare pivot that the bound
   cannot express. */

/* Suffix-close the predicted pattern over its supernode run spans: whenever
   a column's pattern hits a run row, extend it to the run's end.  On the U
   side this is a running merge of the run columns' U lists, which makes the
   construction linear in the closed output.  Closure removes the one-sided
   veto cases entirely, so the pivoting fill can swap any weak pivot with
   any in-run candidate.  Values are left uninitialized; the pivoting fill
   rewrites every pattern position. */
static int kls_predicted_suffix_close_numeric_mode(
  kls_solver *solver, const unsigned char *span_full);

static int kls_predicted_suffix_close_numeric(kls_solver *solver) {
  return kls_predicted_suffix_close_numeric_mode(solver, NULL);
}

static int kls_predicted_suffix_close_numeric_mode(
  kls_solver *solver, const unsigned char *span_full) {
  trilinos_klu_l_numeric *numeric = solver->numeric;
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  trilinos_klu_l_common *common = &solver->common;
  if (numeric == NULL || symbolic == NULL) {
    return 0;
  }
  if (solver->snode_run_end == NULL) {
    /* No supernode runs detected: nothing to close; the pivoting fill's
       per-swap legality vetoes carry the burden alone. */
    return 1;
  }
  const UF_long nblocks = symbolic->nblocks;
  const UF_long *R = symbolic->R;
  const UF_long maxblock = symbolic->maxblock;
  UF_long *merge = (UF_long *)malloc((size_t)maxblock * sizeof(*merge));
  UF_long *merge_next = (UF_long *)malloc((size_t)maxblock *
                                          sizeof(*merge_next));
  UF_long *closed_llen = (UF_long *)malloc((size_t)maxblock *
                                           sizeof(*closed_llen));
  UF_long *closed_ulen_buf = (UF_long *)malloc((size_t)maxblock *
                                               sizeof(*closed_ulen_buf));
  UF_long *cursor = (UF_long *)malloc((size_t)maxblock * sizeof(*cursor));
  if (merge == NULL || merge_next == NULL || closed_llen == NULL ||
      closed_ulen_buf == NULL || cursor == NULL) {
    free(merge);
    free(merge_next);
    free(closed_llen);
    free(closed_ulen_buf);
    free(cursor);
    return 0;
  }
  int ok = 1;
  UF_long total_lnz = 0;
  UF_long total_unz = 0;
  UF_long max_lnz_block = 1;
  UF_long max_unz_block = 1;

  for (UF_long block = 0; block < nblocks && ok; ++block) {
    const UF_long k1 = R[block];
    const UF_long k2 = R[block + 1u];
    const UF_long nk = k2 - k1;
    if (nk <= 1u) {
      total_lnz += 1;
      total_unz += 1;
      continue;
    }
    double *lu = (double *)numeric->LUbx[block];
    if (lu == NULL) {
      ok = 0;
      break;
    }
    UF_long *uip = (UF_long *)malloc((size_t)nk * sizeof(*uip));
    UF_long *ulen = (UF_long *)malloc((size_t)nk * sizeof(*ulen));
    if (uip == NULL || ulen == NULL) {
      free(uip);
      free(ulen);
      ok = 0;
      break;
    }
    memcpy(uip, numeric->Uip + k1, (size_t)nk * sizeof(*uip));
    memcpy(ulen, numeric->Ulen + k1, (size_t)nk * sizeof(*ulen));

    /* Pass 1: full-span union sizes.  For a span S=[s,e) the closed U list
       of column k is (union of the span's unpadded U lists, below k) plus
       the dense in-span range [s,k): every earlier column that touches the
       span then carries its whole remainder, so any in-span row swap keeps
       every stored value inside every pattern. */
    memset(closed_llen, 0, (size_t)nk * sizeof(*closed_llen));
    UF_long *closed_ulen = closed_ulen_buf;
    UF_long block_unz = 0;
    size_t lusize = 0;
    for (UF_long s = 0; s < nk;) {
      UF_long e = solver->snode_run_end[k1 + s] - k1;
      if (e <= s || e > nk) {
        e = s + 1u;
      }
      /* running (suffix) merge: column k sees the union of U lists of
         span columns up to and including itself, plus the dense range */
      UF_long union_len = 0;
      /* A full span pre-accumulates its whole union; the per-column
         merges below are then idempotent and every emit sees the full
         span, which is exactly the closure the vetoed swaps need. */
      if (span_full != NULL && span_full[k1 + s]) {
        UF_long full_len = 0;
        for (UF_long k = s; k < e; ++k) {
          UF_long *ui = NULL;
          double *ux = NULL;
          UF_long ucol_len = 0;
          kls_klu_get_pointer(lu, uip, ulen, k, &ui, &ux, &ucol_len);
          (void)ux;
          UF_long a = 0;
          UF_long b = 0;
          UF_long w = 0;
          while (a < full_len || b < ucol_len) {
            const UF_long va = a < full_len ? merge[a] : nk;
            const UF_long vb = b < ucol_len ? ui[b] : nk;
            UF_long v = va < vb ? va : vb;
            if (va == v) {
              a++;
            }
            if (vb == v) {
              b++;
            }
            merge_next[w++] = v;
          }
          memcpy(merge, merge_next, (size_t)w * sizeof(*merge));
          full_len = w;
        }
        union_len = full_len;
      }
      for (UF_long k = s; k < e; ++k) {
        UF_long *ui = NULL;
        double *ux = NULL;
        UF_long ucol_len = 0;
        kls_klu_get_pointer(lu, uip, ulen, k, &ui, &ux, &ucol_len);
        (void)ux;
        UF_long a = 0;
        UF_long b = 0;
        UF_long w = 0;
        while (a < union_len || b < ucol_len) {
          const UF_long va = a < union_len ? merge[a] : nk;
          const UF_long vb = b < ucol_len ? ui[b] : nk;
          UF_long v = va < vb ? va : vb;
          if (va == v) {
            a++;
          }
          if (vb == v) {
            b++;
          }
          merge_next[w++] = v;
        }
        memcpy(merge, merge_next, (size_t)w * sizeof(*merge));
        union_len = w;
        UF_long upos = 0;
        while (upos < union_len && merge[upos] < k) {
          upos++;
        }
        /* count = |union below k| + |[s,k)| - |union ∩ [s,k)| */
        UF_long in_span_dup = 0;
        for (UF_long p = 0; p < upos; ++p) {
          if (merge[p] >= s) {
            in_span_dup = upos - p;
            break;
          }
        }
        const UF_long cw = upos + (k - s) - in_span_dup;
        closed_ulen[k] = cw;
        block_unz += cw;
        lusize += 2u * (size_t)cw;
        /* L counts via the same membership */
        {
          UF_long p = 0;
          UF_long rr = s;
          while (p < upos || rr < k) {
            const UF_long vu = p < upos ? merge[p] : nk;
            const UF_long vr = rr < k ? rr : nk;
            UF_long v = vu < vr ? vu : vr;
            if (vu == v) {
              p++;
            }
            if (vr == v) {
              rr++;
            }
            closed_llen[v]++;
          }
        }
      }
      s = e;
    }
    UF_long block_lnz = 0;
    for (UF_long jj = 0; jj < nk; ++jj) {
      block_lnz += closed_llen[jj];
      lusize += 2u * (size_t)closed_llen[jj];
    }
    if (lusize == 0) {
      lusize = 1;
    }
    double *new_lu =
      (double *)trilinos_klu_l_malloc(lusize, sizeof(double), common);
    if (new_lu == NULL || common->status < TRILINOS_KLU_OK) {
      common->status = TRILINOS_KLU_OK;
      free(uip);
      free(ulen);
      ok = 0;
      break;
    }
    /* layout */
    UF_long *nlip = numeric->Lip + k1;
    UF_long *nllen = numeric->Llen + k1;
    UF_long *nuip = numeric->Uip + k1;
    UF_long *nulen = numeric->Ulen + k1;
    size_t off = 0;
    for (UF_long k = 0; k < nk; ++k) {
      nlip[k] = (UF_long)off;
      off += 2u * (size_t)closed_llen[k];
      nuip[k] = (UF_long)off;
      off += 2u * (size_t)closed_ulen[k];
    }
    memcpy(cursor, nlip, (size_t)nk * sizeof(*cursor));
    /* Pass 2: emit closed lists and scatter L rows. */
    for (UF_long s = 0; s < nk;) {
      UF_long e = solver->snode_run_end[k1 + s] - k1;
      if (e <= s || e > nk) {
        e = s + 1u;
      }
      UF_long union_len = 0;
      /* A full span pre-accumulates its whole union; the per-column
         merges below are then idempotent and every emit sees the full
         span, which is exactly the closure the vetoed swaps need. */
      if (span_full != NULL && span_full[k1 + s]) {
        UF_long full_len = 0;
        for (UF_long k = s; k < e; ++k) {
          UF_long *ui = NULL;
          double *ux = NULL;
          UF_long ucol_len = 0;
          kls_klu_get_pointer(lu, uip, ulen, k, &ui, &ux, &ucol_len);
          (void)ux;
          UF_long a = 0;
          UF_long b = 0;
          UF_long w = 0;
          while (a < full_len || b < ucol_len) {
            const UF_long va = a < full_len ? merge[a] : nk;
            const UF_long vb = b < ucol_len ? ui[b] : nk;
            UF_long v = va < vb ? va : vb;
            if (va == v) {
              a++;
            }
            if (vb == v) {
              b++;
            }
            merge_next[w++] = v;
          }
          memcpy(merge, merge_next, (size_t)w * sizeof(*merge));
          full_len = w;
        }
        union_len = full_len;
      }
      for (UF_long k = s; k < e; ++k) {
        UF_long *ui = NULL;
        double *ux = NULL;
        UF_long ucol_len = 0;
        kls_klu_get_pointer(lu, uip, ulen, k, &ui, &ux, &ucol_len);
        (void)ux;
        UF_long a = 0;
        UF_long b = 0;
        UF_long w = 0;
        while (a < union_len || b < ucol_len) {
          const UF_long va = a < union_len ? merge[a] : nk;
          const UF_long vb = b < ucol_len ? ui[b] : nk;
          UF_long v = va < vb ? va : vb;
          if (va == v) {
            a++;
          }
          if (vb == v) {
            b++;
          }
          merge_next[w++] = v;
        }
        memcpy(merge, merge_next, (size_t)w * sizeof(*merge));
        union_len = w;
        UF_long upos = 0;
        while (upos < union_len && merge[upos] < k) {
          upos++;
        }
        UF_long *nui = (UF_long *)(new_lu + nuip[k]);
        UF_long p = 0;
        UF_long rr = s;
        UF_long ew = 0;
        while (p < upos || rr < k) {
          const UF_long vu = p < upos ? merge[p] : nk;
          const UF_long vr = rr < k ? rr : nk;
          UF_long v = vu < vr ? vu : vr;
          if (vu == v) {
            p++;
          }
          if (vr == v) {
            rr++;
          }
          nui[ew++] = v;
          UF_long *nli = (UF_long *)(new_lu + nlip[v]);
          nli[cursor[v] - nlip[v]] = k;
          cursor[v] += 1u;
        }
        if (ew != closed_ulen[k]) {
          ok = 0;
          break;
        }
        nulen[k] = ew;
      }
      if (!ok) {
        break;
      }
      s = e;
    }
    if (ok) {
      for (UF_long k = 0; k < nk; ++k) {
        if (cursor[k] - nlip[k] != closed_llen[k]) {
          ok = 0;
          break;
        }
        nllen[k] = closed_llen[k];
      }
    }
    free(uip);
    free(ulen);
    if (!ok) {
      trilinos_klu_l_free(new_lu, lusize, sizeof(double), common);
      break;
    }
    trilinos_klu_l_free(numeric->LUbx[block], numeric->LUsize[block],
                        sizeof(double), common);
    numeric->LUbx[block] = new_lu;
    numeric->LUsize[block] = lusize;
    total_lnz += block_lnz + nk;
    total_unz += block_unz + nk;
    if (block_lnz + nk > max_lnz_block) {
      max_lnz_block = block_lnz + nk;
    }
    if (block_unz + nk > max_unz_block) {
      max_unz_block = block_unz + nk;
    }
  }
  free(merge);
  free(merge_next);
  free(closed_llen);
  free(closed_ulen_buf);
  free(cursor);
  if (!ok) {
    return 0;
  }
  numeric->lnz = total_lnz;
  numeric->unz = total_unz;
  numeric->max_lnz_block = max_lnz_block;
  numeric->max_unz_block = max_unz_block;
  return 1;
}

static int kls_predicted_pivoting_fill_hot(kls_solver *solver,
                                           const double *numeric_values,
                                           unsigned char *span_hot);

static int kls_predicted_pivoting_fill(kls_solver *solver,
                                       const double *numeric_values) {
  return kls_predicted_pivoting_fill_hot(solver, numeric_values, NULL);
}

static int kls_predicted_pivoting_fill_hot(kls_solver *solver,
                                           const double *numeric_values,
                                           unsigned char *span_hot) {
  trilinos_klu_l_numeric *numeric = solver->numeric;
  const trilinos_klu_l_symbolic *symbolic = solver->symbolic;
  trilinos_klu_l_common *common = &solver->common;
  if (numeric == NULL || symbolic == NULL ||
      (symbolic->maxblock > 200000 &&
       getenv("KLS_PREDICTED_PIVOT_FILL_LARGE") == NULL) ||
      common->scale > 0 ||
      numeric->Rs != NULL) {
    if (getenv("KLS_TRACE_PREDICTED") != NULL) {
      fprintf(stderr, "KLS predicted: pivot fill gate1 maxblock=%ld\n",
              symbolic != NULL ? (long)symbolic->maxblock : -1L);
    }
    return 0;
  }
  if (!kls_build_refactor_map(solver) ||
      solver->refactor_col_ptr == NULL ||
      (symbolic->nblocks > 1u && solver->refactor_block_start == NULL) ||
      solver->refactor_row_idx == NULL ||
      solver->refactor_input_pos == NULL) {
    if (getenv("KLS_TRACE_PREDICTED") != NULL) {
      fprintf(stderr, "KLS predicted: pivot fill gate2 map=%d %d %d %d\n",
              solver->refactor_col_ptr != NULL,
              solver->refactor_block_start != NULL,
              solver->refactor_row_idx != NULL,
              solver->refactor_input_pos != NULL);
    }
    return 0;
  }
  const UF_long *map_block_start =
    symbolic->nblocks == 1u ? solver->refactor_col_ptr
                            : solver->refactor_block_start;
  const UF_long nblocks = symbolic->nblocks;
  const UF_long *R = symbolic->R;
  const UF_long *Q = symbolic->Q;
  const UF_long maxblock = symbolic->maxblock;
  const double tol = common->tol > 0.0 ? common->tol : 0.001;

  double *x = (double *)calloc((size_t)maxblock, sizeof(*x));
  UF_long *cur_slot = (UF_long *)malloc((size_t)maxblock * sizeof(*cur_slot));
  UF_long *row_cols_ptr =
    (UF_long *)malloc(((size_t)maxblock + 1u) * sizeof(*row_cols_ptr));
  unsigned char *cand_tried =
    (unsigned char *)calloc((size_t)maxblock, sizeof(*cand_tried));
  UF_long *row_cols = NULL;
  UF_long row_cols_capacity = 0;
  if (x == NULL || cur_slot == NULL || row_cols_ptr == NULL ||
      cand_tried == NULL) {
    free(x);
    free(cur_slot);
    free(row_cols_ptr);
    free(cand_tried);
    return 0;
  }

  int ok = 1;
  UF_long swap_count = 0;
  UF_long veto_count = 0;
  for (UF_long block = 0; block < nblocks && ok; ++block) {
    const UF_long k1 = R[block];
    const UF_long k2 = R[block + 1u];
    const UF_long nk = k2 - k1;
    if (nk == 1u) {
      const UF_long gk = k1;
      UF_long poff = numeric->Offp[gk];
      const UF_long poff_end = numeric->Offp[gk + 1u];
      double s = 0.0;
      for (UF_long p = solver->refactor_col_ptr[gk];
           p < map_block_start[gk]; ++p) {
        if (poff >= poff_end) {
          ok = 0;
          break;
        }
        ((double *)numeric->Offx)[poff++] =
          numeric_values[solver->refactor_input_pos[p]];
      }
      for (UF_long p = map_block_start[gk];
           ok && p < solver->refactor_col_ptr[gk + 1u]; ++p) {
        s = numeric_values[solver->refactor_input_pos[p]];
      }
      ((double *)numeric->Udiag)[gk] = s;
      continue;
    }
    double *lu = (double *)numeric->LUbx[block];
    const UF_long *lip = numeric->Lip + k1;
    const UF_long *llen = numeric->Llen + k1;
    const UF_long *uip = numeric->Uip + k1;
    const UF_long *ulen = numeric->Ulen + k1;
    if (lu == NULL) {
      ok = 0;
      break;
    }
    /* Row lists of the L pattern: which columns contain each slot. */
    memset(row_cols_ptr, 0, ((size_t)nk + 1u) * sizeof(*row_cols_ptr));
    UF_long lnz_block = 0;
    for (UF_long k = 0; k < nk; ++k) {
      UF_long *li = NULL;
      double *lx = NULL;
      UF_long len = 0;
      kls_klu_get_pointer(lu, lip, llen, k, &li, &lx, &len);
      (void)lx;
      for (UF_long p = 0; p < len; ++p) {
        row_cols_ptr[li[p] + 1u]++;
      }
      lnz_block += len;
    }
    for (UF_long i = 0; i < nk; ++i) {
      row_cols_ptr[i + 1u] += row_cols_ptr[i];
    }
    if (lnz_block > row_cols_capacity) {
      free(row_cols);
      row_cols = (UF_long *)malloc((size_t)lnz_block * sizeof(*row_cols));
      row_cols_capacity = row_cols != NULL ? lnz_block : 0;
      if (row_cols == NULL) {
        ok = 0;
        break;
      }
    }
    {
      UF_long *cursor = cur_slot; /* scratch before its real use */
      memcpy(cursor, row_cols_ptr, (size_t)nk * sizeof(*cursor));
      for (UF_long k = 0; k < nk; ++k) {
        UF_long *li = NULL;
        double *lx = NULL;
        UF_long len = 0;
        kls_klu_get_pointer(lu, lip, llen, k, &li, &lx, &len);
        (void)lx;
        for (UF_long p = 0; p < len; ++p) {
          row_cols[cursor[li[p]]++] = k;
        }
      }
    }
    for (UF_long i = 0; i < nk; ++i) {
      cur_slot[i] = i;
    }

    for (UF_long k = 0; k < nk && ok; ++k) {
      const UF_long gk = k1 + k;
      /* Scatter the input column through the current slot relabeling. */
      UF_long poff = numeric->Offp[gk];
      const UF_long poff_end = numeric->Offp[gk + 1u];
      for (UF_long p = solver->refactor_col_ptr[gk];
           p < map_block_start[gk]; ++p) {
        if (poff >= poff_end) {
          ok = 0;
          break;
        }
        ((double *)numeric->Offx)[poff++] =
          numeric_values[solver->refactor_input_pos[p]];
      }
      if (!ok) {
        break;
      }
      for (UF_long p = map_block_start[gk];
           p < solver->refactor_col_ptr[gk + 1u]; ++p) {
        const UF_long slot0 = solver->refactor_row_idx[p] - k1;
        x[cur_slot[slot0]] =
          numeric_values[solver->refactor_input_pos[p]];
      }
      /* Left-looking updates over the fixed U pattern. */
      UF_long *ui = NULL;
      double *ux = NULL;
      UF_long ucol_len = 0;
      kls_klu_get_pointer(lu, uip, ulen, k, &ui, &ux, &ucol_len);
      for (UF_long up = 0; up < ucol_len; ++up) {
        const UF_long j = ui[up];
        const double ujk = x[j];
        x[j] = 0.0;
        ux[up] = ujk;
        if (ujk != 0.0) {
          UF_long *li = NULL;
          double *lx = NULL;
          UF_long lcol_len = 0;
          kls_klu_get_pointer(lu, lip, llen, j, &li, &lx, &lcol_len);
          for (UF_long p = 0; p < lcol_len; ++p) {
            x[li[p]] -= lx[p] * ujk;
          }
        }
      }
      /* Threshold pivot within the pattern column. */
      UF_long *li = NULL;
      double *lx = NULL;
      UF_long lcol_len = 0;
      kls_klu_get_pointer(lu, lip, llen, k, &li, &lx, &lcol_len);
      double ukk = x[k];
      double cmax = fabs(ukk);
      UF_long best = KLS_KLU_EMPTY;
      /* Candidates stay inside the supernode run: path closure of the
         symmetrized bound guarantees any earlier column containing one of
         the run's rows contains those above it, so an in-run swap never
         moves a stored value outside a column's pattern.  Rows beyond the
         run offer no such guarantee and measurably destroy the factor. */
      UF_long run_end_local = nk;
      if (solver->snode_run_end != NULL) {
        run_end_local = solver->snode_run_end[gk] - k1;
        if (run_end_local <= k || run_end_local > nk) {
          run_end_local = k + 1u;
        }
      } else {
        /* No runs: any pattern row may serve; the veto pass rejects the
           structurally illegal ones swap by swap. */
        run_end_local = nk;
      }
      UF_long in_run_len = 0;
      while (in_run_len < lcol_len && li[in_run_len] < run_end_local) {
        in_run_len++;
      }
      double colmax_all = fabs(ukk);
      for (UF_long p = 0; p < lcol_len; ++p) {
        const double av = fabs(x[li[p]]);
        colmax_all = colmax_all < av ? av : colmax_all;
      }
      cmax = colmax_all;
      /* Candidates by descending magnitude; a veto moves to the next one
         instead of giving up on the swap entirely. */
      if (fabs(ukk) < tol * cmax) {
        for (UF_long p = 0; p < in_run_len; ++p) {
          cand_tried[li[p]] = 0;
        }
        for (int attempt = 0; attempt < 8; ++attempt) {
          double best_av = fabs(ukk);
          UF_long pick = KLS_KLU_EMPTY;
          for (UF_long p = 0; p < in_run_len; ++p) {
            if (cand_tried[li[p]]) {
              continue;
            }
            const double av = fabs(x[li[p]]);
            if (av > best_av) {
              best_av = av;
              pick = li[p];
            }
          }
          if (pick == KLS_KLU_EMPTY) {
            break;
          }
          cand_tried[pick] = 1;
          best = pick;
          break;
        }
      }
      while (best != KLS_KLU_EMPTY && fabs(ukk) < tol * cmax) {
        /* A swap is only performed when every earlier column carries both
           slots or neither: a one-sided column would have to drop a stored
           value outside its fixed pattern, corrupting the factor. */
        int vetoed = 0;
        {
          UF_long va = row_cols_ptr[k];
          const UF_long va_end = row_cols_ptr[k + 1u];
          UF_long vb = row_cols_ptr[best];
          const UF_long vb_end = row_cols_ptr[best + 1u];
          while (va < va_end || vb < vb_end) {
            const UF_long ca = va < va_end ? row_cols[va] : nk;
            const UF_long cb = vb < vb_end ? row_cols[vb] : nk;
            const UF_long col = ca < cb ? ca : cb;
            if (col >= k) {
              break;
            }
            const int in_a = va < va_end && ca == col;
            const int in_b = vb < vb_end && cb == col;
            if (in_a != in_b) {
              vetoed = 1;
              break;
            }
            va += in_a;
            vb += in_b;
          }
        }
        if (vetoed) {
          veto_count++;
          if (span_hot != NULL && solver->snode_run_end != NULL) {
            /* mark the span so the targeted closure can open it fully */
            UF_long ss = k;
            while (ss > 0 &&
                   k1 + ss < solver->snode_run_end[k1 + ss - 1u]) {
              ss--;
            }
            span_hot[k1 + ss] = 1;
          }
          /* next candidate by magnitude */
          double best_av = fabs(ukk);
          best = KLS_KLU_EMPTY;
          for (UF_long p = 0; p < in_run_len; ++p) {
            if (cand_tried[li[p]]) {
              continue;
            }
            const double av = fabs(x[li[p]]);
            if (av > best_av) {
              best_av = av;
              best = li[p];
            }
          }
          if (best != KLS_KLU_EMPTY) {
            cand_tried[best] = 1;
          }
          continue;
        }
        break;
      }
      if (best != KLS_KLU_EMPTY && fabs(ukk) < tol * cmax) {
        /* Swap slots k and best: workspace, pending scatters, and the
           already-stored L values of earlier columns at those slots. */
        const double tmp = x[k];
        x[k] = x[best];
        x[best] = tmp;
        for (UF_long i = 0; i < nk; ++i) {
          if (cur_slot[i] == k) {
            cur_slot[i] = best;
          } else if (cur_slot[i] == best) {
            cur_slot[i] = k;
          }
        }
        const UF_long gbest = k1 + best;
        const UF_long pk = numeric->Pnum[gk];
        numeric->Pnum[gk] = numeric->Pnum[gbest];
        numeric->Pnum[gbest] = pk;
        numeric->Pinv[numeric->Pnum[gk]] = gk;
        numeric->Pinv[numeric->Pnum[gbest]] = gbest;
        /* Prior columns holding either slot swap their stored values. */
        UF_long pa = row_cols_ptr[k];
        const UF_long pa_end = row_cols_ptr[k + 1u];
        UF_long pb = row_cols_ptr[best];
        const UF_long pb_end = row_cols_ptr[best + 1u];
        while (pa < pa_end || pb < pb_end) {
          UF_long col;
          if (pa < pa_end && (pb >= pb_end || row_cols[pa] <= row_cols[pb])) {
            col = row_cols[pa];
          } else {
            col = row_cols[pb];
          }
          if (col >= k) {
            break;
          }
          UF_long *cli = NULL;
          double *clx = NULL;
          UF_long clen = 0;
          kls_klu_get_pointer(lu, lip, llen, col, &cli, &clx, &clen);
          UF_long ia = KLS_KLU_EMPTY;
          UF_long ib = KLS_KLU_EMPTY;
          if (pa < pa_end && row_cols[pa] == col) {
            UF_long lo = 0;
            UF_long hi = clen;
            while (lo < hi) {
              const UF_long mid = (lo + hi) / 2u;
              if (cli[mid] < k) {
                lo = mid + 1u;
              } else {
                hi = mid;
              }
            }
            if (lo < clen && cli[lo] == k) {
              ia = lo;
            }
            pa++;
          }
          if (pb < pb_end && row_cols[pb] == col) {
            UF_long lo = 0;
            UF_long hi = clen;
            while (lo < hi) {
              const UF_long mid = (lo + hi) / 2u;
              if (cli[mid] < best) {
                lo = mid + 1u;
              } else {
                hi = mid;
              }
            }
            if (lo < clen && cli[lo] == best) {
              ib = lo;
            }
            pb++;
          }
          if (ia != KLS_KLU_EMPTY && ib != KLS_KLU_EMPTY) {
            const double tv = clx[ia];
            clx[ia] = clx[ib];
            clx[ib] = tv;
          } else if (ia != KLS_KLU_EMPTY || ib != KLS_KLU_EMPTY) {
            /* The veto pass excluded one-sided columns; reaching here means
               the row lists and patterns disagree. */
            ok = 0;
            break;
          }
        }
        ukk = x[k];
        swap_count++;
      }
      x[k] = 0.0;
      if (ukk == 0.0) {
        /* No in-run candidate: fall back to a diagonal correction that the
           persistent nudge list replays on every later refactorization. */
        double sigma = tol * (colmax_all > 0.0 ? colmax_all : 1.0);
        UF_long diag_pos = KLS_KLU_EMPTY;
        for (UF_long p = map_block_start[gk];
             p < solver->refactor_col_ptr[gk + 1u]; ++p) {
          const UF_long slot0 = solver->refactor_row_idx[p] - k1;
          if (cur_slot[slot0] == k) {
            diag_pos = solver->refactor_input_pos[p];
            break;
          }
        }
        if (diag_pos == KLS_KLU_EMPTY ||
            solver->pivot_nudge_count >=
              (solver->pivot_nudge_capacity > 0
                 ? solver->pivot_nudge_capacity : 0)) {
          if (diag_pos == KLS_KLU_EMPTY && solver->pivot_nudge_values == NULL) {
            /* allocate below */
          }
        }
        if (diag_pos != KLS_KLU_EMPTY) {
          if (solver->pivot_nudge_values == NULL) {
            solver->pivot_nudge_values = (double *)malloc(
              (size_t)solver->nnz * sizeof(*solver->pivot_nudge_values));
            solver->pivot_nudge_pos = (UF_long *)malloc(
              1024u * sizeof(*solver->pivot_nudge_pos));
            solver->pivot_nudge_sigma = (double *)malloc(
              1024u * sizeof(*solver->pivot_nudge_sigma));
            if (solver->pivot_nudge_values == NULL ||
                solver->pivot_nudge_pos == NULL ||
                solver->pivot_nudge_sigma == NULL) {
              free_pivot_nudges(solver);
              diag_pos = KLS_KLU_EMPTY;
            } else {
              solver->pivot_nudge_capacity = 1024u;
            }
          }
        }
        if (diag_pos != KLS_KLU_EMPTY &&
            solver->pivot_nudge_count < solver->pivot_nudge_capacity) {
          solver->pivot_nudge_pos[solver->pivot_nudge_count] = diag_pos;
          solver->pivot_nudge_sigma[solver->pivot_nudge_count] = sigma;
          solver->pivot_nudge_count++;
          ukk = sigma;
          if (getenv("KLS_TRACE_PREDICTED") != NULL) {
            fprintf(stderr,
                    "KLS predicted: pivot fill nudge col=%ld sigma=%.3e\n",
                    (long)gk, sigma);
          }
        } else {
          common->status = TRILINOS_KLU_SINGULAR;
          common->numerical_rank = gk;
          common->singular_col = Q[gk];
          ok = 0;
          break;
        }
      }
      ((double *)numeric->Udiag)[gk] = ukk;
      for (UF_long p = 0; p < lcol_len; ++p) {
        lx[p] = x[li[p]] / ukk;
        x[li[p]] = 0.0;
      }
    }
  }
  free(x);
  free(cur_slot);
  free(row_cols_ptr);
  free(cand_tried);
  free(row_cols);
  if (getenv("KLS_TRACE_PREDICTED") != NULL) {
    fprintf(stderr, "KLS predicted: pivot fill swaps=%ld vetoes=%ld ok=%d\n",
            (long)swap_count, (long)veto_count, ok);
  }
  if (!ok) {
    return 0;
  }
  common->status = TRILINOS_KLU_OK;
  return 1;
}

typedef struct kls_block_trial_entry {
  UF_long row;
  UF_long base_pos;
  double value;
} kls_block_trial_entry;

static int kls_block_trial_entry_cmp(const void *a, const void *b) {
  const kls_block_trial_entry *ea = (const kls_block_trial_entry *)a;
  const kls_block_trial_entry *eb = (const kls_block_trial_entry *)b;
  return (ea->row > eb->row) - (ea->row < eb->row);
}

/* Block-structured ordering: circuit families such as TSOPF repeat a
   physics block (one power network per contingency) that generic orderings
   interleave, forcing the pivots that inflate fill.  Detect the blocks as
   connected components after setting aside the dense coupling vertices,
   order each component contiguously (reverse Cuthill-McKee inside), append
   the couplers last, and let the predicted factorization with restricted
   pivoting run on the aligned spans. */
static int kls_ufl_cmp(const void *a, const void *b) {
  const UF_long x = *(const UF_long *)a;
  const UF_long y = *(const UF_long *)b;
  return x < y ? -1 : x > y ? 1 : 0;
}

static int kls_build_block_structured_order(UF_long n,
                                            const UF_long *col_ptr,
                                            const UF_long *row_idx,
                                            UF_long **perm_out,
                                            UF_long **comp_out) {
  *perm_out = NULL;
  *comp_out = NULL;
  if (n < 10000 || n > 2000000) {
    return 0;
  }
  const UF_long nnz = col_ptr[n];
  UF_long *deg = (UF_long *)calloc((size_t)n, sizeof(*deg));
  if (deg == NULL) {
    return 0;
  }
  for (UF_long j = 0; j < n; ++j) {
    for (UF_long p = col_ptr[j]; p < col_ptr[j + 1u]; ++p) {
      const UF_long i = row_idx[p];
      if (i == j) {
        continue;
      }
      deg[i]++;
      deg[j]++;
    }
  }
  UF_long *xadj = (UF_long *)malloc(((size_t)n + 1u) * sizeof(*xadj));
  UF_long *adj = (UF_long *)malloc(2u * (size_t)nnz * sizeof(*adj));
  if (xadj == NULL || adj == NULL) {
    free(deg);
    free(xadj);
    free(adj);
    return 0;
  }
  xadj[0] = 0;
  for (UF_long v = 0; v < n; ++v) {
    xadj[v + 1u] = xadj[v] + deg[v];
    deg[v] = xadj[v];
  }
  for (UF_long j = 0; j < n; ++j) {
    for (UF_long p = col_ptr[j]; p < col_ptr[j + 1u]; ++p) {
      const UF_long i = row_idx[p];
      if (i == j) {
        continue;
      }
      adj[deg[i]++] = j;
      adj[deg[j]++] = i;
    }
  }
  /* degrees again (deg was reused as cursor); duplicates in adj are
     harmless for BFS/RCM */
  for (UF_long v = 0; v < n; ++v) {
    deg[v] = xadj[v + 1u] - xadj[v];
  }
  /* Coupler vertices sit orders of magnitude above the block-interior
     degrees (TSOPF b9: interiors p99=14, couplers 34..7210); a floor
     relative to the p99 degree separates them where absolute thresholds
     miss the low-degree couplers. */
  UF_long dense_floor = 24;
  {
    UF_long *deg_sorted = (UF_long *)malloc((size_t)n * sizeof(*deg_sorted));
    if (deg_sorted != NULL) {
      memcpy(deg_sorted, deg, (size_t)n * sizeof(*deg_sorted));
      qsort(deg_sorted, (size_t)n, sizeof(*deg_sorted), kls_ufl_cmp);
      const UF_long p99 = deg_sorted[(size_t)((double)(n - 1) * 0.99)];
      if (2u * p99 > dense_floor) {
        dense_floor = 2u * p99;
      }
      free(deg_sorted);
    } else {
      const double avg = (double)xadj[n] / (double)n;
      if ((UF_long)(avg * 8.0) > dense_floor) {
        dense_floor = (UF_long)(avg * 8.0);
      }
    }
  }
  {
    const char *floor_env = getenv("KLS_BLOCK_DENSE_FLOOR");
    if (floor_env != NULL && floor_env[0] != '\0') {
      const long parsed = atol(floor_env);
      if (parsed > 0) {
        dense_floor = (UF_long)parsed;
      }
    }
  }
  UF_long dense_count = 0;
  for (UF_long v = 0; v < n; ++v) {
    if (deg[v] > dense_floor) {
      dense_count++;
    }
  }
  if (dense_count == 0 || dense_count > n / 100u) {
    free(deg);
    free(xadj);
    free(adj);
    return 0;
  }
  UF_long *comp = (UF_long *)malloc((size_t)n * sizeof(*comp));
  UF_long *queue = (UF_long *)malloc((size_t)n * sizeof(*queue));
  UF_long *perm = (UF_long *)malloc((size_t)n * sizeof(*perm));
  if (comp == NULL || queue == NULL || perm == NULL) {
    free(deg);
    free(xadj);
    free(adj);
    free(comp);
    free(queue);
    free(perm);
    return 0;
  }
  for (UF_long v = 0; v < n; ++v) {
    comp[v] = KLS_KLU_EMPTY;
  }
  UF_long ncomp = 0;
  UF_long covered = 0;
  UF_long max_comp = 0;
  /* Pass 1: collect components contiguously into queue[] so they can be
     emitted largest-first; a big leading component anchors the band and
     the matched diagonal downstream. */
  UF_long *comp_start =
    (UF_long *)malloc(((size_t)n + 1u) * sizeof(*comp_start));
  UF_long *comp_order = (UF_long *)malloc((size_t)n * sizeof(*comp_order));
  unsigned char *placed = (unsigned char *)calloc((size_t)n, 1u);
  if (comp_start == NULL || comp_order == NULL || placed == NULL) {
    free(comp_start);
    free(comp_order);
    free(placed);
    free(deg);
    free(xadj);
    free(adj);
    free(comp);
    free(queue);
    free(perm);
    return 0;
  }
  UF_long fill_ptr = 0;
  for (UF_long v = 0; v < n; ++v) {
    if (comp[v] != KLS_KLU_EMPTY || deg[v] > dense_floor) {
      continue;
    }
    comp_start[ncomp] = fill_ptr;
    UF_long head = fill_ptr;
    queue[fill_ptr++] = v;
    comp[v] = ncomp;
    while (head < fill_ptr) {
      const UF_long u = queue[head++];
      for (UF_long p = xadj[u]; p < xadj[u + 1u]; ++p) {
        const UF_long w = adj[p];
        if (comp[w] != KLS_KLU_EMPTY || deg[w] > dense_floor) {
          continue;
        }
        comp[w] = ncomp;
        queue[fill_ptr++] = w;
      }
    }
    const UF_long sz = fill_ptr - comp_start[ncomp];
    covered += sz;
    if (sz > max_comp) {
      max_comp = sz;
    }
    ncomp++;
  }
  comp_start[ncomp] = fill_ptr;
  for (UF_long c = 0; c < ncomp; ++c) {
    comp_order[c] = c;
  }
  for (UF_long a = 1; a < ncomp; ++a) {
    const UF_long c = comp_order[a];
    const UF_long csz = comp_start[c + 1u] - comp_start[c];
    UF_long b = a;
    while (b > 0 &&
           comp_start[comp_order[b - 1u] + 1u] -
             comp_start[comp_order[b - 1u]] < csz) {
      comp_order[b] = comp_order[b - 1u];
      b--;
    }
    comp_order[b] = c;
  }
  /* Pass 2, largest component first: Reverse Cuthill-McKee from the
     component's minimum-degree member with each vertex's undiscovered
     neighbours enqueued in ascending degree, then the whole order
     reversed.  The degree sorting is what keeps the component banded;
     without it the no-pivot fill bound explodes. */
  UF_long write = 0;
  for (UF_long ci = 0; ci < ncomp; ++ci) {
    const UF_long c = comp_order[ci];
    const UF_long ms = comp_start[c];
    const UF_long me = comp_start[c + 1u];
    UF_long start = queue[ms];
    for (UF_long p = ms; p < me; ++p) {
      if (deg[queue[p]] < deg[start]) {
        start = queue[p];
      }
    }
    UF_long head = 0;
    UF_long tail = 0;
    UF_long *cq = perm + write; /* RCM order built in place, reversed after */
    cq[tail++] = start;
    placed[start] = 1;
    while (head < tail) {
      const UF_long u = cq[head++];
      const UF_long first = tail;
      for (UF_long p = xadj[u]; p < xadj[u + 1u]; ++p) {
        const UF_long w = adj[p];
        if (placed[w] || comp[w] != c) {
          continue;
        }
        placed[w] = 1;
        cq[tail++] = w;
      }
      for (UF_long a = first + 1u; a < tail; ++a) {
        const UF_long w = cq[a];
        UF_long b = a;
        while (b > first && deg[cq[b - 1u]] > deg[w]) {
          cq[b] = cq[b - 1u];
          b--;
        }
        cq[b] = w;
      }
    }
    for (UF_long p = ms; p < me; ++p) {
      if (!placed[queue[p]]) {
        placed[queue[p]] = 1;
        cq[tail++] = queue[p];
      }
    }
    for (UF_long a = 0; a < tail / 2u; ++a) {
      const UF_long tmp = cq[a];
      cq[a] = cq[tail - 1u - a];
      cq[tail - 1u - a] = tmp;
    }
    write += tail;
  }
  /* Repeated-physics-block signature: the near-equal-size components must
     carry most of the matrix.  Irregular component profiles (hvdc2's
     40632/34723/30926..., circuit_4's 28k fragments) factor better under
     the regular orderings and only waste a trial here. */
  UF_long equal_mass = 0;
  for (UF_long c = 0; c < ncomp; ++c) {
    const UF_long sz = comp_start[c + 1u] - comp_start[c];
    if ((double)sz >= 0.95 * (double)max_comp) {
      equal_mass += sz;
    }
  }
  free(comp_start);
  free(comp_order);
  free(placed);
  for (UF_long v = 0; v < n; ++v) {
    if (deg[v] > dense_floor) {
      perm[write++] = v;
    }
  }
  free(deg);
  free(xadj);
  free(adj);
  free(queue);
  const int usable =
    write == n && ncomp >= 3u && covered * 10u >= n * 6u &&
    max_comp * 2u <= n && equal_mass * 2u >= n;
  if (!usable) {
    free(perm);
    free(comp);
    return 0;
  }
  *perm_out = perm;
  *comp_out = comp;
  return 1;
}

#ifdef KLS_HAVE_SPRAL_SCALING
static void maybe_select_block_structured_ordering(kls_solver *solver,
                                                   double *elapsed,
                                                   const double *values) {
  if (solver == NULL || !solver->options.static_pivoting ||
      solver->numeric != NULL || solver->row_perm != NULL ||
      solver->user_col_perm != NULL ||
      solver->input_format != KLS_INPUT_CSC ||
      solver->options.ordering != KLS_ORDERING_AUTO ||
      solver->options.scale > 0) {
    return;
  }
  const char *block_ordering_env = getenv("KLS_ENABLE_BLOCK_ORDERING");
  if (block_ordering_env != NULL && block_ordering_env[0] == '0' &&
      block_ordering_env[1] == '\0') {
    /* Keep an explicit master disable for controlled A/B.  In generic AUTO
       this routine only proposes a candidate: a measured factor trial,
       residual probe, and fill bound below still decide adoption. */
    return;
  }
  /* Reordering thousands of already-tiny full-rank SCCs is more work than
     their bounded AMD symbolic.  This decision is structural; the stricter
     factor profile still controls every recurring numeric shortcut. */

  if (solver->n > 200000 ||
      (solver->col_ptr != NULL && solver->col_ptr[solver->n] > 4000000)) {
    /* class ceiling, same as the choose fast path and the race gate:
       every block-structured adopter sits under 200K/4M, and the
       detector walk has a bad worst case on large patterns (the
       rajat29 lesson: 1.9s hunting a class the matrix cannot be) */
    return;
  }
  const double start = kls_now_seconds();
  UF_long *block_perm = solver->block_order_perm;
  solver->block_order_perm = NULL;
  UF_long *block_comp = NULL;
  if (block_perm == NULL &&
      !kls_build_block_structured_order(solver->n, solver->col_ptr,
                                        solver->row_idx, &block_perm,
                                        &block_comp)) {
    *elapsed += kls_now_seconds() - start;
    return;
  }
  if (getenv("KLS_TRACE_PREDICTED") != NULL) {
    fprintf(stderr, "KLS block: structure detected %.2fs\n",
            kls_now_seconds() - start);
  }

  const UF_long n = solver->n;
  const UF_long nnz = solver->nnz;
  UF_long *row_pos = NULL;      /* user row -> internal position */
  UF_long *trial_col_ptr = NULL;
  UF_long *trial_row_idx = NULL;
  UF_long *trial_input = NULL;
  double *trial_values = NULL;
  UF_long *trial_row_perm = NULL;
  UF_long *trial_col_perm = NULL;
  trilinos_klu_l_symbolic *trial_symbolic = NULL;
  trilinos_klu_l_common trial_common;
  (void)trilinos_klu_l_defaults(&trial_common);
  int accepted = 0;

  /* Band-capped Hungarian: only entries within a fixed distance of the
     blocked diagonal participate (couplers exempt), with proximity
     weights.  The matched diagonal is then banded by construction - the
     degenerate value-optimal matchings that scramble the band are simply
     not in the graph - and the restricted pivoting fill supplies the
     numerics afterwards. */
  {
    int64_t *sp_col_ptr =
      (int64_t *)malloc(((size_t)n + 1u) * sizeof(*sp_col_ptr));
    int *sp_row_idx = (int *)malloc((size_t)nnz * sizeof(*sp_row_idx));
    double *sp_values = (double *)malloc((size_t)nnz * sizeof(*sp_values));
    int *sp_match = (int *)malloc((size_t)n * sizeof(*sp_match));
    double *sp_rs = (double *)malloc((size_t)n * sizeof(*sp_rs));
    double *sp_cs = (double *)malloc((size_t)n * sizeof(*sp_cs));
    UF_long *pos_of = (UF_long *)malloc((size_t)n * sizeof(*pos_of));
    if (sp_col_ptr == NULL || sp_row_idx == NULL || sp_values == NULL ||
        sp_match == NULL || sp_rs == NULL || sp_cs == NULL ||
        pos_of == NULL || n > (UF_long)INT_MAX) {
      free(sp_col_ptr);
      free(sp_row_idx);
      free(sp_values);
      free(sp_match);
      free(sp_rs);
      free(sp_cs);
      free(pos_of);
      goto done;
    }
    for (UF_long k = 0; k < n; ++k) {
      pos_of[block_perm[k]] = k;
    }
    UF_long count = 0;
    for (UF_long k = 0; k < n; ++k) {
      const UF_long oldcol = block_perm[k];
      sp_col_ptr[k] = (int64_t)count;
      for (UF_long p = solver->col_ptr[oldcol];
           p < solver->col_ptr[oldcol + 1u]; ++p) {
        const UF_long r = solver->row_idx[p];
        if (values[p] == 0.0) {
          continue;
        }
        /* Raw value-optimal matching: with the largest-first min-degree
           RCM column order the value-optimal diagonal follows the physics
           of the blocks, and the band emerges without a distance cap.
           (Band-capped proximity weights were compensating for a weaker
           column order and block the same matching now.) */
        sp_row_idx[count] = (int)r;
        sp_values[count] = fabs(values[p]);
        count++;
      }
    }
    sp_col_ptr[n] = (int64_t)count;
    free(pos_of);
    struct spral_scaling_hungarian_options sp_options;
    struct spral_scaling_hungarian_inform sp_inform;
    spral_scaling_hungarian_default_options(&sp_options);
    sp_options.array_base = 0;
    sp_options.scale_if_singular = false;
    spral_scaling_hungarian_unsym_long((int)n, (int)n, sp_col_ptr,
                                       sp_row_idx, sp_values, sp_rs, sp_cs,
                                       sp_match, &sp_options, &sp_inform);
    free(sp_col_ptr);
    free(sp_row_idx);
    free(sp_values);
    free(sp_rs);
    free(sp_cs);
    if (sp_inform.flag != 0 && sp_inform.flag != -2) {
      free(sp_match);
      goto done;
    }
    row_pos = (UF_long *)malloc((size_t)n * sizeof(*row_pos));
    if (row_pos == NULL) {
      free(sp_match);
      goto done;
    }
    UF_long matched = 0;
    for (UF_long r = 0; r < n; ++r) {
      const int pos = sp_match[r];
      if (pos >= 0 && pos < (int)n) {
        row_pos[r] = (UF_long)pos;
        matched++;
      } else {
        row_pos[r] = KLS_KLU_EMPTY;
      }
    }
    free(sp_match);
    if (matched < n - n / 50u) {
      if (getenv("KLS_TRACE_PREDICTED") != NULL) {
        fprintf(stderr, "KLS block: match failed (%ld/%ld)\n",
                (long)matched, (long)n);
      }
      goto done;
    }
    if (matched != n) {
      /* Phase two: the starved columns (unmatchable without a dense row)
         and the unmatched rows both move to the tail, keeping the band
         intact ahead of them. */
      unsigned char *col_matched =
        (unsigned char *)calloc((size_t)n, sizeof(*col_matched));
      UF_long *new_perm = (UF_long *)malloc((size_t)n * sizeof(*new_perm));
      if (col_matched == NULL || new_perm == NULL) {
        free(col_matched);
        free(new_perm);
        goto done;
      }
      for (UF_long r = 0; r < n; ++r) {
        if (row_pos[r] != KLS_KLU_EMPTY) {
          col_matched[row_pos[r]] = 1;
        }
      }
      UF_long w = 0;
      for (UF_long k = 0; k < n; ++k) {
        if (col_matched[k]) {
          new_perm[w++] = block_perm[k];
        }
      }
      const UF_long tail_start = w;
      for (UF_long k = 0; k < n; ++k) {
        if (!col_matched[k]) {
          new_perm[w++] = block_perm[k];
        }
      }
      /* remap matched rows to the compacted positions */
      UF_long *pos_map = (UF_long *)malloc((size_t)n * sizeof(*pos_map));
      if (pos_map == NULL) {
        free(col_matched);
        free(new_perm);
        goto done;
      }
      UF_long w2 = 0;
      for (UF_long k = 0; k < n; ++k) {
        if (col_matched[k]) {
          pos_map[k] = w2++;
        } else {
          pos_map[k] = KLS_KLU_EMPTY;
        }
      }
      for (UF_long r = 0; r < n; ++r) {
        if (row_pos[r] != KLS_KLU_EMPTY) {
          row_pos[r] = pos_map[row_pos[r]];
        }
      }
      /* unmatched rows fill the tail positions in order */
      UF_long next_tail = tail_start;
      for (UF_long r = 0; r < n; ++r) {
        if (row_pos[r] == KLS_KLU_EMPTY) {
          row_pos[r] = next_tail++;
        }
      }
      memcpy(block_perm, new_perm, (size_t)n * sizeof(*new_perm));
      free(col_matched);
      free(new_perm);
      free(pos_map);
      if (getenv("KLS_TRACE_PREDICTED") != NULL) {
        fprintf(stderr, "KLS block: %ld starved columns moved to tail\n",
                (long)(n - matched));
      }
    }
  }

  /* Physical position-space matrix: internal entry (row_pos[r], k) for the
     user entry (r, block_perm[k]) - byte for byte what the validated
     file-based prototype factored. */
  trial_col_ptr = (UF_long *)malloc(((size_t)n + 1u) * sizeof(*trial_col_ptr));
  trial_row_idx = (UF_long *)malloc((size_t)nnz * sizeof(*trial_row_idx));
  trial_input = (UF_long *)malloc((size_t)nnz * sizeof(*trial_input));
  trial_values = (double *)malloc((size_t)nnz * sizeof(*trial_values));
  if (trial_col_ptr == NULL || trial_row_idx == NULL || trial_input == NULL ||
      trial_values == NULL) {
    goto done;
  }
  {
    UF_long *base_to_trial =
      (UF_long *)malloc((size_t)nnz * sizeof(*base_to_trial));
    kls_block_trial_entry *col_buf = NULL;
    UF_long col_buf_capacity = 0;
    if (base_to_trial == NULL) {
      goto done;
    }
    UF_long count = 0;
    for (UF_long k = 0; k < n; ++k) {
      const UF_long oldcol = block_perm[k];
      trial_col_ptr[k] = count;
      const UF_long len =
        solver->col_ptr[oldcol + 1u] - solver->col_ptr[oldcol];
      if (len > col_buf_capacity) {
        free(col_buf);
        col_buf = (kls_block_trial_entry *)malloc(
          (size_t)len * sizeof(*col_buf));
        col_buf_capacity = col_buf != NULL ? len : 0;
        if (col_buf == NULL) {
          free(base_to_trial);
          goto done;
        }
      }
      for (UF_long t = 0; t < len; ++t) {
        const UF_long p = solver->col_ptr[oldcol] + t;
        col_buf[t].row = row_pos[solver->row_idx[p]];
        col_buf[t].base_pos = p;
        col_buf[t].value = values[p];
      }
      qsort(col_buf, (size_t)len, sizeof(*col_buf),
            kls_block_trial_entry_cmp);
      for (UF_long t = 0; t < len; ++t) {
        trial_row_idx[count] = col_buf[t].row;
        trial_values[count] = col_buf[t].value;
        base_to_trial[col_buf[t].base_pos] = count;
        count++;
      }
    }
    free(col_buf);
    trial_col_ptr[n] = count;
    /* compose the input-to-internal map */
    if (solver->input_to_csc != NULL) {
      for (UF_long q = 0; q < nnz; ++q) {
        trial_input[q] = base_to_trial[solver->input_to_csc[q]];
      }
    } else {
      for (UF_long q = 0; q < nnz; ++q) {
        trial_input[q] = base_to_trial[q];
      }
    }
    free(base_to_trial);
  }
  trial_row_perm = row_pos;
  row_pos = NULL;
  trial_col_perm = block_perm;
  block_perm = NULL;

  if (getenv("KLS_TRACE_PREDICTED") != NULL) {
    fprintf(stderr, "KLS block: match+build done %.2fs\n",
            kls_now_seconds() - start);
  }
  trial_symbolic = trilinos_klu_l_analyze_given(n, trial_col_ptr,
                                                trial_row_idx, NULL, NULL,
                                                &trial_common);
  if (getenv("KLS_TRACE_PREDICTED") != NULL) {
    fprintf(stderr, "KLS block: analyze done %.2fs\n",
            kls_now_seconds() - start);
  }
  if (trial_symbolic == NULL || trial_common.status < TRILINOS_KLU_OK) {
    if (getenv("KLS_TRACE_PREDICTED") != NULL) {
      fprintf(stderr, "KLS block: analyze failed status=%d\n",
              (int)trial_common.status);
    }
    goto done;
  }
  trial_common.scale = -1;
  /* The matched diagonal is the pivot policy: a low threshold keeps
     partial pivoting on it through elimination-phase value drift (at the
     default 1e-3 KLU wanders off and fill explodes ~7x); the solve probe
     and refinement police the numerics. */
  trial_common.tol =
    solver->common.tol < 1.0e-4 ? solver->common.tol : 1.0e-4;
  if (solver->symbolic != NULL && solver->symbolic->lnz > 0.0 &&
      solver->symbolic->unz > 0.0 && trial_symbolic->lnz > 0.0 &&
      trial_symbolic->unz > 0.0 &&
      trial_symbolic->lnz + trial_symbolic->unz >
        4.0 * (solver->symbolic->lnz + solver->symbolic->unz)) {
    /* The trial's own symbolic bound already dwarfs the baseline estimate;
       factoring it just to reject on exact fill wastes seconds.  Real fill
       runs below the bound, so keep a generous margin. */
    if (getenv("KLS_TRACE_PREDICTED") != NULL) {
      fprintf(stderr,
              "KLS block: trial skipped, est=%.3g > 4x base est=%.3g\n",
              trial_symbolic->lnz + trial_symbolic->unz,
              solver->symbolic->lnz + solver->symbolic->unz);
    }
    goto done;
  }

  {
    trilinos_klu_l_symbolic *old_symbolic = solver->symbolic;
    trilinos_klu_l_common old_common = solver->common;
    UF_long *old_col_ptr = solver->col_ptr;
    UF_long *old_row_idx = solver->row_idx;
    UF_long *old_input = solver->input_to_csc;
    double *old_values_arr = solver->values;

    solver->col_ptr = trial_col_ptr;
    solver->row_idx = trial_row_idx;
    solver->input_to_csc = trial_input;
    solver->row_perm = trial_row_perm;
    solver->user_col_perm = trial_col_perm;
    solver->values = trial_values;
    solver->symbolic = trial_symbolic;
    solver->common = trial_common;
    kls_invalidate_factor_etree_stats(solver);

    solver->block_trial_active = 1;
    int factored = 0;
    /* With a clean value-matched diagonal, real partial pivoting stays on
       it and the actual fill runs far below the predicted pattern's etree
       bound; every one of the cycle's refactorizations then pays the
       smaller pattern.  Weak or missing diagonals fall through to the
       predicted trial, whose restricted pivoting handles them. */
    {
      /* Threshold pivoting at the trial's tol=1e-4 stays on any diagonal
         at least that strong; a handful of weak columns just pivot
         off-diagonal locally, so only a weak *population* forces the
         predicted trial. */
      const UF_long weak_cap = n / 1024u > 4u ? n / 1024u : 4u;
      UF_long weak = 0;
      for (UF_long k = 0; k < n && weak <= weak_cap; ++k) {
        double colmax = 0.0;
        double diag = -1.0;
        for (UF_long p = trial_col_ptr[k]; p < trial_col_ptr[k + 1u];
             ++p) {
          const double av = fabs(trial_values[p]);
          colmax = colmax < av ? av : colmax;
          if (trial_row_idx[p] == k) {
            diag = av;
          }
        }
        if (diag < 0.0 || diag < 1.0e-4 * colmax) {
          weak++;
        }
      }
      if (weak <= weak_cap) {
        const double plain_start = kls_now_seconds();
        trilinos_klu_l_numeric *plain =
          trilinos_klu_l_factor(trial_col_ptr, trial_row_idx, trial_values,
                                solver->symbolic, &solver->common);
        if (plain != NULL && solver->common.status >= TRILINOS_KLU_OK &&
            solver->common.status != TRILINOS_KLU_SINGULAR &&
            plain->lnz + plain->unz <= 6 * (UF_long)nnz) {
          solver->numeric = plain;
          solver->numeric_is_predicted = 0;
          factored = 1;
          if (getenv("KLS_TRACE_PREDICTED") != NULL) {
            fprintf(stderr,
                    "KLS block: plain KLU trial lnz=%ld unz=%ld %.2fs\n",
                    (long)plain->lnz, (long)plain->unz,
                    kls_now_seconds() - plain_start);
          }
        } else {
          if (getenv("KLS_TRACE_PREDICTED") != NULL) {
            fprintf(stderr,
                    "KLS block: plain KLU trial rejected status=%d "
                    "fill=%ld\n",
                    (int)solver->common.status,
                    plain != NULL ? (long)(plain->lnz + plain->unz) : -1L);
          }
          if (plain != NULL) {
            trilinos_klu_l_free_numeric(&plain, &solver->common);
          }
          solver->common.status = TRILINOS_KLU_OK;
        }
      }
    }
    if (!factored) {
      factored =
        kls_predicted_pattern_first_factor(solver, trial_values, elapsed);
    }
    solver->block_trial_active = 0;
    if (getenv("KLS_TRACE_PREDICTED") != NULL) {
      fprintf(stderr, "KLS block: trial factor done %.2fs\n",
              kls_now_seconds() - start);
    }

    if (factored && solver->numeric != NULL) {
      /* Acceptance guard: the baseline ordering's symbolic estimate is the
         fill the cycle would have paid without the block trial; a trial
         numeric denser than that estimate loses every refactorization and
         must not replace it. */
      const double base_est = old_symbolic != NULL &&
                              old_symbolic->lnz > 0.0 &&
                              old_symbolic->unz > 0.0
                                ? old_symbolic->lnz + old_symbolic->unz
                                : -1.0;
      const double trial_fill =
        (double)(solver->numeric->lnz + solver->numeric->unz);
      if (base_est > 0.0 && trial_fill > base_est) {
        if (getenv("KLS_TRACE_PREDICTED") != NULL) {
          fprintf(stderr,
                  "KLS block: adoption rejected fill=%.3g > base est=%.3g\n",
                  trial_fill, base_est);
        }
        trilinos_klu_l_free_numeric(&solver->numeric, &solver->common);
        solver->numeric = NULL;
        kls_numeric_replaced_invalidate(solver);
        solver->numeric_is_predicted = 0;
        factored = 0;
      }
    }
    if (factored && solver->numeric != NULL) {
      solver->exact_matching_selected = 1;
      solver->spral_matching_selected = 1;
      solver->stats.selected_ordering = KLS_ORDERING_NATURAL;
      solver->stats.selected_orientation = solver->orientation;
      solver->stats.nblocks = (int64_t)solver->symbolic->nblocks;
      solver->stats.max_block = (int64_t)solver->symbolic->maxblock;
      solver->stats.structural_rank =
        (int64_t)solver->symbolic->structural_rank;
      solver->stats.estimated_flops = solver->symbolic->est_flops;
      /* The analyze-time METIS race and its scale-trial workers read
         the pattern arrays freed below: join and discard them first,
         exactly as the pre-static adoption does.  TSOPF_FS_b39_c19:
         two scale-trial threads segfaulted in lsolve_symbolic on the
         freed base pattern (16T runs corrupted into timeouts). */
      kls_metis_race_abandon(solver);
      trilinos_klu_l_free_symbolic(&old_symbolic, &old_common);
      free(old_col_ptr);
      free(old_row_idx);
      free(old_input);
      free(old_values_arr);
      trial_col_ptr = NULL;
      trial_row_idx = NULL;
      trial_input = NULL;
      trial_values = NULL;
      trial_row_perm = NULL;
      trial_col_perm = NULL;
      trial_symbolic = NULL;
      accepted = 1;
      if (getenv("KLS_TRACE_PREDICTED") != NULL) {
        fprintf(stderr, "KLS block: ordering adopted\n");
      }
    } else {
      solver->col_ptr = old_col_ptr;
      solver->row_idx = old_row_idx;
      solver->input_to_csc = old_input;
      solver->row_perm = NULL;
      solver->user_col_perm = NULL;
      solver->values = old_values_arr;
      solver->symbolic = old_symbolic;
      solver->common = old_common;
      kls_invalidate_factor_etree_stats(solver);
    }
  }

done:
  if (!accepted) {
    if (trial_symbolic != NULL) {
      trilinos_klu_l_free_symbolic(&trial_symbolic, &trial_common);
    }
    free(trial_col_ptr);
    free(trial_row_idx);
    free(trial_input);
    free(trial_values);
    free(trial_row_perm);
    free(trial_col_perm);
  }
  free(row_pos);
  free(block_perm);
  free(block_comp);
  *elapsed += kls_now_seconds() - start;
}
#endif

/* Rescue path: suffix-close the pattern, rebuild the derived structures on
   it, and refactor with pattern-restricted pivoting.  The closed fill must
   stay near the ordering estimate, or the pivoted KLU factorization is the
   better numeric and the rescue declines. */
static int kls_predicted_pivot_rescue(kls_solver *solver,
                                      const double *numeric_values,
                                      double *elapsed) {
  if ((getenv("KLS_ENABLE_PIVOT_FILL") == NULL &&
       !solver->block_trial_active) ||
      solver->numeric == NULL || solver->symbolic == NULL) {
    return 0;
  }
  const double est_total =
    solver->block_trial_active
      ? -1.0
      : solver->symbolic->lnz + solver->symbolic->unz;
  if (est_total > 0.0 &&
      (double)(solver->numeric->lnz + solver->numeric->unz) >
        1.2 * est_total) {
    return 0;
  }
  if (solver->snode_run_end == NULL) {
    free_snode_panels(solver);
    solver->snode_numeric_pre_sorted = 1;
    kls_maybe_prepare_snode_panels(solver, elapsed);
  kls_snb_maybe_accept(solver, numeric_values, elapsed);
  (void)kls_i32_solve_ready(solver);
  }
  if (solver->block_trial_active) {
    /* Block trials swap into structurally absent slots by design: the
       full closure is their legalization.  The probe-retry path instead
       fills first with vetoes and opens only the veto-hot spans - the
       full closure explodes on spiked patterns (rajat25: 11x est). */
    if (!kls_predicted_suffix_close_numeric(solver)) {
      return 0;
    }
    if (est_total > 0.0 &&
        (double)(solver->numeric->lnz + solver->numeric->unz) >
          1.6 * est_total) {
      if (getenv("KLS_TRACE_PREDICTED") != NULL) {
        fprintf(stderr,
                "KLS predicted: closure fill reject %.3e vs est %.3e\n",
                (double)(solver->numeric->lnz + solver->numeric->unz),
                est_total);
      }
      return 0;
    }
    if (getenv("KLS_TRACE_PREDICTED") != NULL) {
      fprintf(stderr, "KLS predicted: closed fill=%ld+%ld\n",
              (long)solver->numeric->lnz, (long)solver->numeric->unz);
    }
  }
  /* The fill must pivot within the very spans the closure guaranteed, so
     the pre-closure run partition stays live through it; the wider runs of
     the closed pattern are re-detected afterwards for the refactor path. */
  free_refactor_lu_pointer_cache(solver);
  free_refactor_map(solver);
  free_refactor_schedule(solver);
  unsigned char *span_hot =
    (unsigned char *)calloc((size_t)solver->n, sizeof(*span_hot));
  int filled =
    kls_predicted_pivoting_fill_hot(solver, numeric_values, span_hot);
  if (filled && span_hot != NULL) {
    UF_long weak_pivot = KLS_KLU_EMPTY;
    UF_long weak_col = KLS_KLU_EMPTY;
    int any_hot = 0;
    for (UF_long k = 0; k < solver->n && !any_hot; ++k) {
      any_hot = span_hot[k];
    }
    if (any_hot &&
        !kls_numeric_pivots_pass_threshold(solver, &weak_pivot, &weak_col)) {
      /* Weak pivots survive behind structural vetoes: open the veto-hot
         spans fully - the minimal closure that legalizes their swaps -
         and refill once. */
      if (getenv("KLS_TRACE_PREDICTED") != NULL) {
        UF_long hot_count = 0;
        for (UF_long k = 0; k < solver->n; ++k) {
          hot_count += span_hot[k];
        }
        fprintf(stderr,
                "KLS predicted: opening %ld veto-hot spans fully\n",
                (long)hot_count);
      }
      free_refactor_lu_pointer_cache(solver);
      free_refactor_map(solver);
      free_refactor_schedule(solver);
      if (kls_predicted_suffix_close_numeric_mode(solver, span_hot)) {
        if (getenv("KLS_TRACE_PREDICTED") != NULL) {
          fprintf(stderr, "KLS predicted: reclosed fill=%ld+%ld\n",
                  (long)solver->numeric->lnz, (long)solver->numeric->unz);
        }
        if (est_total > 0.0 &&
            (double)(solver->numeric->lnz + solver->numeric->unz) >
              1.6 * est_total) {
          if (getenv("KLS_TRACE_PREDICTED") != NULL) {
            fprintf(stderr,
                    "KLS predicted: reclose fill reject %.3e vs est %.3e\n",
                    (double)(solver->numeric->lnz + solver->numeric->unz),
                    est_total);
          }
          filled = 0;
        } else {
          filled = kls_predicted_pivoting_fill(solver, numeric_values);
        }
      }
    }
  }
  free(span_hot);
  /* The refill's own post-preps rebuilt the refactor map, and snode panel
     preparation declines while position-based retained structure exists;
     drop the map first so the closed pattern gets its run detection, then
     rebuild the map against the sorted storage. */
  free_refactor_map(solver);
  free_refactor_schedule(solver);
  free_snode_panels(solver);
  solver->snode_numeric_pre_sorted = 1;
  kls_maybe_prepare_snode_panels(solver, elapsed);
  kls_snb_maybe_accept(solver, numeric_values, elapsed);
  (void)kls_i32_solve_ready(solver);
  maybe_prepare_refactor_map(solver, elapsed);
  maybe_prepare_refactor_schedule(solver, elapsed);
  return filled;
}

/* ------------------------------------------------------------------ */
/* i32-stream solve (hard core 2, brick 1).  klu_l_solve walks 16     */
/* bytes per factor entry (UF_long index + double); the row indices   */
/* never change across refactorizations, so a solver-lifetime int32   */
/* mirror cuts the solve's index traffic in half while values are     */
/* still read straight from the packed numeric.  Semantics replicate  */
/* trilinos_klu_l_solve for one right-hand side.                      */
/* ------------------------------------------------------------------ */

/* ---------------------------------------------------------------------
   Subtree-partitioned parallel triangular solve (hard core C, brick 2).

   Built once per i32-stream build for the largest BTF block.  The block
   etree is derived from the streams themselves (parent = min over the
   L(:,j) off-diagonal rows and the U row constraints) and every stream
   entry is verified against it with Euler-tour intervals: row i of
   L(:,k) must be an ancestor of k, row i of U(:,k) a descendant of k.
   Any violation, unsorted stream, or unbalanced partition declines to
   the serial i32 sweep, so the structural assumptions are checked, not
   trusted.  Columns are partitioned into per-thread bins of disjoint
   subtrees; scatter targets outside the own chunk are provably in the
   shared ancestor "top" set, accumulated privately, and merged before
   the serial top sweeps.  Backward U needs no accumulation at all: its
   scatter targets are descendants, which stay inside the own chunk.  */
#define KLS_PTS_MAX_THREADS 16

typedef struct kls_pts_stream_run {
  int32_t row;
  uint32_t length;
} kls_pts_stream_run;

typedef struct kls_pts_s {
  UF_long block;                    /* BTF block this partition covers */
  UF_long k1;
  UF_long nk;
  int nthreads;
  int32_t *tcols;                   /* per-bin ascending column lists */
  int64_t tcols_off[KLS_PTS_MAX_THREADS + 1];
  int32_t *top_cols;                /* ascending shared-ancestor columns */
  int64_t ntop;
  int32_t *top_map;                 /* nk: compact top slot or -1 */
  int32_t *lsplit;                  /* nk: in-chunk L stream prefix */
  double *top_acc;                  /* nthreads * ntop, kept zeroed */
  double *xwork;                    /* nthreads * nk, lazy, kept zeroed */
  double *refactor_top_x;           /* ntop * nk two-stage workspaces */
  int solve_ok;                     /* tight top gate for the solve */
  int solve_direct_verified;        /* tight plan may skip timing */
  int refactor_ok;                  /* flop-weighted top gate */
  double refactor_top_flops;        /* realized work in shared top */
  double refactor_total_flops;      /* realized work in full block */
  int two_stage_refactor_decision;  /* 0 trial, 1 two-stage, -1 serial */
  int two_stage_refactor_samples[2];
  double two_stage_refactor_min[2];
  int top_run_scatter;              /* contiguous top stream walk */
  uint32_t *top_l_run_ptr;          /* ntop+1 offsets into top_l_runs */
  uint32_t *top_u_run_ptr;          /* ntop+1 offsets into top_u_runs */
  kls_pts_stream_run *top_l_runs;
  kls_pts_stream_run *top_u_runs;
  int solve_decision;               /* 0 trial, 1 PTS, -1 serial i32 */
  int solve_probe_samples[2];       /* [0] serial, [1] PTS */
  double solve_probe_min[2];        /* warm minimum for each solve arm */
  int egraph_terminal_decision;     /* 0 unknown, 1 eligible, -1 rejected */
  UF_long egraph_terminal_col;
} kls_pts;

typedef struct kls_pts_touch_job {
  double *x;
  size_t count;
  int tid;
  int nthreads;
} kls_pts_touch_job;
