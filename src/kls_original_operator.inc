/* Private integration preparation: caller-frame CSC from the retained maps
 * and CURRENT caller-order values. Never undo rounded scaling to recover A.
 * Public numeric updates publish this snapshot transactionally. */
typedef struct private_original_operator {
  size_t n, nnz;
  int64_t *p, *rows;
  double *values;
} private_original_operator;

static void private_original_operator_free(private_original_operator *a) {
  if (!a) return;
  free(a->p); free(a->rows); free(a->values);
  memset(a, 0, sizeof(*a));
}

static int private_original_operator_build(const kls_solver *s,
    const double *caller_values, private_original_operator *out) {
  if (!s || !out || !caller_values || s->n <= 0 || s->nnz < 0 ||
      !s->col_ptr || (s->nnz && !s->row_idx) ||
      (s->orientation != KLS_ORIENTATION_NORMAL &&
       s->orientation != KLS_ORIENTATION_TRANSPOSE) ||
      (uintmax_t)s->n > SIZE_MAX / sizeof(int64_t) - 1 ||
      (uintmax_t)s->nnz > SIZE_MAX / sizeof(int64_t) ||
      (uintmax_t)s->nnz > SIZE_MAX / sizeof(double))
    return KLS_ERR_INVALID_ARGUMENT;
  const size_t n = (size_t)s->n, nnz = (size_t)s->nnz;
  private_original_operator result = {n, nnz, NULL, NULL, NULL};
  int status = KLS_ERR_OUT_OF_MEMORY;
  int64_t *inverse_row = malloc(n * sizeof(*inverse_row));
  int64_t *cursor = malloc(n * sizeof(*cursor));
  int64_t *source = malloc((nnz ? nnz : 1) * sizeof(*source));
  result.p = calloc(n + 1, sizeof(*result.p));
  result.rows = malloc((nnz ? nnz : 1) * sizeof(*result.rows));
  result.values = malloc((nnz ? nnz : 1) * sizeof(*result.values));
  if (!inverse_row || !cursor || !source || !result.p ||
      !result.rows || !result.values) goto done;
  status = KLS_ERR_INVALID_ARGUMENT;
  if (s->col_ptr[0] != 0 || s->col_ptr[n] != s->nnz) goto done;
  for (size_t j = 0; j < n; ++j) {
    inverse_row[j] = -1;
    cursor[j] = -1;
    if (s->col_ptr[j] < 0 || s->col_ptr[j] > s->col_ptr[j + 1] ||
        s->col_ptr[j + 1] > s->nnz) goto done;
  }
  for (size_t j = 0; j < n; ++j) {
    const UF_long row = s->row_perm ? s->row_perm[j] : (UF_long)j;
    const UF_long col = s->user_col_perm ? s->user_col_perm[j] : (UF_long)j;
    if (row < 0 || (uintmax_t)row >= n || inverse_row[row] != -1 ||
        col < 0 || (uintmax_t)col >= n || cursor[col] != -1) goto done;
    inverse_row[row] = (int64_t)j;
    cursor[col] = (int64_t)j;
  }
  for (size_t k = 0; k < nnz; ++k) source[k] = -1;
  for (size_t k = 0; k < nnz; ++k) {
    const UF_long target = s->input_to_csc ? s->input_to_csc[k] : (UF_long)k;
    if (target < 0 || (uintmax_t)target >= nnz || source[target] != -1)
      goto done;
    source[target] = (int64_t)k;
  }
  const int transposed = s->orientation == KLS_ORIENTATION_TRANSPOSE;
  for (size_t j = 0; j < n; ++j) {
    const int64_t col = s->user_col_perm ? s->user_col_perm[j] : (int64_t)j;
    for (UF_long k = s->col_ptr[j]; k < s->col_ptr[j + 1]; ++k) {
      if (s->row_idx[k] < 0 || (uintmax_t)s->row_idx[k] >= n) goto done;
      const int64_t row = inverse_row[s->row_idx[k]];
      ++result.p[(transposed ? row : col) + 1];
    }
  }
  for (size_t j = 0; j < n; ++j) {
    result.p[j + 1] += result.p[j];
    cursor[j] = result.p[j];
  }
  for (size_t j = 0; j < n; ++j) {
    const int64_t col = s->user_col_perm ? s->user_col_perm[j] : (int64_t)j;
    for (UF_long k = s->col_ptr[j]; k < s->col_ptr[j + 1]; ++k) {
      const int64_t row = inverse_row[s->row_idx[k]];
      const int64_t at = cursor[transposed ? row : col]++;
      result.rows[at] = transposed ? col : row;
      if(s->accuracy_policy==KLS_ACCURACY_APPLICATION_MANAGED &&
         !isfinite(caller_values[source[k]])) goto done;
      memcpy(result.values + at, caller_values + source[k], sizeof(double));
    }
  }
  private_original_operator_free(out);
  *out = result;
  memset(&result, 0, sizeof(result));
  status = KLS_OK;
done:
  free(inverse_row); free(cursor); free(source);
  private_original_operator_free(&result);
  return status;
}
