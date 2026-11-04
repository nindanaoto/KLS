/* Exact input-slot plan for a fixed-pattern elimination forest. The caller
 * has already validated the single-block map32 layout. No values, RHSs,
 * readiness or accuracy verdicts are cached. Map/owner invalidation releases
 * the plan, including when an in-place numeric layout is replaced. */
static int kls_prepare_forest_map(kls_solver *s) {
  if (s->forest_map_state != 0) return s->forest_map_state > 0;
  s->forest_map_state = -1;
#if FLT_EVAL_METHOD == 0
  const trilinos_klu_l_numeric *num = s->numeric;
  const double *lu = (const double *)num->LUbx[0];
  if (s->n <= 1 || s->n > (UF_long)INT32_MAX ||
      s->nnz > (UF_long)(SIZE_MAX / sizeof(int32_t))) return 0;
  /* Most factors decline without allocating or examining their input map. */
  for (UF_long k = 0; k < s->n; ++k)
    if (num->Llen[k] > 1u) return 0;
  int32_t *slots = malloc((size_t)s->nnz * sizeof(*slots));
  int32_t *row_slot = malloc((size_t)s->n * sizeof(*row_slot));
  if (!slots || !row_slot) { free(slots); free(row_slot); return 0; }
  for (UF_long i = 0; i < s->n; ++i) row_slot[i] = -3;
  UF_long updates = 0;
  for (UF_long k = 0; k < s->n; ++k) {
    const UF_long first = s->refactor_col_ptr[k];
    const UF_long last = s->refactor_col_ptr[k + 1u];
    const UF_long count = num->Ulen[k];
    if (count >= s->n || last - first != count + num->Llen[k] + 1u)
      goto decline;
    const UF_long *ui = count ? (const UF_long *)(lu + num->Uip[k]) : NULL;
    row_slot[k] = -1;
    if (num->Llen[k]) {
      const UF_long lower = ((const UF_long *)(lu + num->Lip[k]))[0];
      if (lower <= k || lower >= s->n) goto decline;
      row_slot[lower] = -2;
    }
    for (UF_long up = 0; up < count; ++up) {
      const UF_long j = ui[up];
      if (j >= k || row_slot[j] != -3 || num->Llen[j] != 1u ||
          ((const UF_long *)(lu + num->Lip[j]))[0] != k) goto decline;
      row_slot[j] = (int32_t)up;
    }
    for (UF_long p = first; p < last; ++p) {
      const int32_t row = s->refactor_row_idx32[p];
      if (row < 0 || (UF_long)row >= s->n || row_slot[row] == -3)
        goto decline;
      slots[p] = row_slot[row];
      row_slot[row] = -3; /* Each original entry must occur exactly once. */
    }
    updates += count;
  }
  if (!updates) goto decline;
  free(row_slot);
  s->forest_input_slot = slots;
  s->forest_map_state = 1;
  return 1;
decline:
  free(slots); free(row_slot);
#endif
  return 0;
}

/* Direct stores publish U inputs earlier than the general scatter does.
 * Retain the incumbent if traps could observe that intermediate state. */
static int kls_forest_masked_environment(void) {
#if defined(__x86_64__) && FLT_EVAL_METHOD == 0
  unsigned short control;
  __asm__ volatile("fnstcw %0" : "=m"(control));
  return (control & 63u) == 63u &&
    (_mm_getcsr() & 0x1f80u) == 0x1f80u;
#else
  return 0;
#endif
}
