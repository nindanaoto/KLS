/* Private integration preparation: caller-frame CSC from the retained maps
 * and CURRENT caller-order values. Never undo rounded scaling to recover A.
 * Public numeric updates publish this snapshot transactionally. */
typedef struct private_original_operator {
  size_t n, nnz;
  int64_t *p, *rows;
  double *values;
  /* Optional exact structural snapshot and caller-index gather map. Values
   * and certificates are never cached across numeric updates. */
  UF_long *layout;
  int64_t *value_source;
  unsigned map_flags;
  int orientation;
} private_original_operator;

static void private_original_operator_free(private_original_operator *a) {
  if (!a) return;
  free(a->p); free(a->rows); free(a->values);
  free(a->layout); free(a->value_source);
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
  const unsigned flags=(s->row_perm?1u:0u)|(s->user_col_perm?2u:0u)|
                       (s->input_to_csc?4u:0u);
  if(out->layout && out->value_source && out->n==n && out->nnz==nnz &&
     out->orientation==(int)s->orientation && out->map_flags==flags) {
    const UF_long *saved=out->layout;
    int same=!memcmp(saved,s->col_ptr,(n+1)*sizeof(*saved));saved+=n+1;
    if(nnz) same &= !memcmp(saved,s->row_idx,nnz*sizeof(*saved));
    saved+=nnz;
    if(s->row_perm) same &= !memcmp(saved,s->row_perm,n*sizeof(*saved));
    saved+=n;
    if(s->user_col_perm) same &= !memcmp(saved,s->user_col_perm,n*sizeof(*saved));
    saved+=n;
    if(s->input_to_csc && nnz) same &= !memcmp(saved,s->input_to_csc,nnz*sizeof(*saved));
    if(same) {
      double *values=malloc((nnz?nnz:1)*sizeof(*values));
      if(!values) return KLS_ERR_OUT_OF_MEMORY;
      for(size_t k=0;k<nnz;++k) {
        const double *value=caller_values+out->value_source[k];
        if(s->accuracy_policy==KLS_ACCURACY_APPLICATION_MANAGED && !isfinite(*value)) {
          free(values);return KLS_ERR_INVALID_ARGUMENT;
        }
        memcpy(values+k,value,sizeof(*values));
      }
      free(out->values);out->values=values;
      return KLS_OK;
    }
  }
  private_original_operator result = {.n=n,.nnz=nnz};
  int status = KLS_ERR_OUT_OF_MEMORY;
  int64_t *inverse_row = malloc(n * sizeof(*inverse_row));
  int64_t *cursor = malloc(n * sizeof(*cursor));
  int64_t *source = malloc((nnz ? nnz : 1) * sizeof(*source));
  result.p = calloc(n + 1, sizeof(*result.p));
  result.rows = malloc((nnz ? nnz : 1) * sizeof(*result.rows));
  result.values = malloc((nnz ? nnz : 1) * sizeof(*result.values));
  if (!inverse_row || !cursor || !source || !result.p ||
      !result.rows || !result.values) goto done;
  /* Optional acceleration: allocation failure retains the full rebuild path.
   * Exact content comparisons also detect in-place permutation/map changes. */
  const size_t capacity=SIZE_MAX/sizeof(UF_long);
  if(n<=(capacity-1)/3 && nnz<=(capacity-1-3*n)/2) {
    result.layout=malloc((3*n+2*nnz+1)*sizeof(*result.layout));
    result.value_source=malloc((nnz?nnz:1)*sizeof(*result.value_source));
    if(!result.layout || !result.value_source) {
      free(result.layout);free(result.value_source);
      result.layout=NULL;result.value_source=NULL;
    }
  }
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
      if(result.value_source) result.value_source[at]=source[k];
    }
  }
  if(result.layout) {
    UF_long *saved=result.layout;
    memcpy(saved,s->col_ptr,(n+1)*sizeof(*saved));saved+=n+1;
    if(nnz) memcpy(saved,s->row_idx,nnz*sizeof(*saved));
    saved+=nnz;
    if(s->row_perm) memcpy(saved,s->row_perm,n*sizeof(*saved));
    saved+=n;
    if(s->user_col_perm) memcpy(saved,s->user_col_perm,n*sizeof(*saved));
    saved+=n;
    if(s->input_to_csc && nnz) memcpy(saved,s->input_to_csc,nnz*sizeof(*saved));
    result.map_flags=flags;result.orientation=s->orientation;
  }
  /* Positions in a certificate row plan refer to the caller CSC layout.
   * Validate that exact layout before retaining a plan across value updates;
   * pointer identity is insufficient because snapshots are replaced. */
  if(out->n!=n || out->nnz!=nnz ||
     memcmp(out->p,result.p,(n+1)*sizeof(*result.p)) ||
     memcmp(out->rows,result.rows,nnz*sizeof(*result.rows)))
    kls_componentwise_plan_reset(s->componentwise_plan);
  private_original_operator_free(out);
  *out = result;
  memset(&result, 0, sizeof(result));
  status = KLS_OK;
done:
  free(inverse_row); free(cursor); free(source);
  private_original_operator_free(&result);
  return status;
}
