static int read_mtx(const char *path, I *n_out, I **cp_out, I **ri_out, double **vx_out) {
  FILE *f = fopen(path, "r");
  if (!f) return 0;
  char line[1024];
  if (!fgets(line, sizeof line, f)) return 0;
  int symmetric = strstr(line, "symmetric") != NULL;
  while (fgets(line, sizeof line, f) && line[0] == '%') {}
  long nr, nc, nnz;
  if (sscanf(line, "%ld %ld %ld", &nr, &nc, &nnz) != 3 || nr != nc) return 0;
  long cap = symmetric ? 2 * nnz : nnz;
  I *rr = malloc((size_t)cap * sizeof(I));
  I *cc = malloc((size_t)cap * sizeof(I));
  double *vv = malloc((size_t)cap * sizeof(double));
  long k = 0;
  for (long e = 0; e < nnz; ++e) {
    long r, c; double v;
    if (fscanf(f, "%ld %ld %lf", &r, &c, &v) != 3) { fclose(f); return 0; }
    rr[k] = r - 1; cc[k] = c - 1; vv[k] = v; k++;
    if (symmetric && r != c) { rr[k] = c - 1; cc[k] = r - 1; vv[k] = v; k++; }
  }
  fclose(f);
  nnz = k;
  I *cp = calloc((size_t)nr + 1, sizeof(I));
  for (long e = 0; e < nnz; ++e) cp[cc[e] + 1]++;
  for (long j = 0; j < nr; ++j) cp[j + 1] += cp[j];
  I *ri = malloc((size_t)nnz * sizeof(I));
  double *vx = malloc((size_t)nnz * sizeof(double));
  I *w = malloc((size_t)nr * sizeof(I));
  memcpy(w, cp, (size_t)nr * sizeof(I));
  for (long e = 0; e < nnz; ++e) { I p = w[cc[e]]++; ri[p] = rr[e]; vx[p] = vv[e]; }
  free(rr); free(cc); free(vv); free(w);
  for (long j = 0; j < nr; ++j) {
    I lo = cp[j], hi = cp[j + 1];
    for (I a = lo + 1; a < hi; ++a) {
      I r0 = ri[a]; double v = vx[a]; I b = a - 1;
      while (b >= lo && ri[b] > r0) { ri[b+1]=ri[b]; vx[b+1]=vx[b]; b--; }
      ri[b+1]=r0; vx[b+1]=v;
    }
  }
  *n_out = nr; *cp_out = cp; *ri_out = ri; *vx_out = vx;
  return 1;
}

