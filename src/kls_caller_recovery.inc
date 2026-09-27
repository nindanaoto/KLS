/* Cold proposals use the exact owned caller matrix. Transposition is a sparse
 * coordinate operation, never an inverse of rounded row/column scaling.
 * Search limits are the existing private lattice policy, not acceptance rules. */
static int private_caller_recover(kls_solver *s,const private_original_operator *a,
    const double *b, double *x, int transpose, double deadline) {
  if (s->private_cold_trial || !a || !a->n || a->n > INT_MAX-64 ||
      !(deadline > kls_now_seconds())) return 0;
  const int64_t *p = a->p, *rows = a->rows;
  const double *values = a->values;
  int64_t *tp = NULL, *tr = NULL, *next = NULL;
  double *tv = NULL;
  int accepted = 0;
  if (transpose) {
    tp = calloc(a->n+1,sizeof(*tp));
    next = malloc(a->n*sizeof(*next));
    tr = malloc((a->nnz ? a->nnz : 1)*sizeof(*tr));
    tv = malloc((a->nnz ? a->nnz : 1)*sizeof(*tv));
    if (!tp || !next || !tr || !tv) goto done;
    for (size_t k=0;k<a->nnz;++k) ++tp[a->rows[k]+1];
    for (size_t j=0;j<a->n;++j) { tp[j+1]+=tp[j];next[j]=tp[j]; }
    for (size_t j=0;j<a->n;++j)
      for (int64_t k=a->p[j];k<a->p[j+1];++k) {
        int64_t at=next[a->rows[k]]++;
        tr[at]=(int64_t)j;tv[at]=a->values[k];
      }
    p=tp;rows=tr;values=tv;
  }
  {
    int proposals=0;
    accepted=private_fresh_proposal(s,a->n,p,rows,values,b,x,deadline);
    if(!accepted) accepted=kls_accuracy_lattice_schedule((int)a->n,p,rows,values,b,x,
      KLS_COLD_LATTICE_FIRST_WIDTH,KLS_COLD_LATTICE_LAST_WIDTH,
      KLS_COLD_LATTICE_SWEEPS,KLS_COLD_LATTICE_ROW_LIMIT,
      KLS_COLD_LATTICE_REGULARIZER,deadline,&proposals)==1;
  }
done:
  free(tp);free(tr);free(next);free(tv);
  return accepted;
}
