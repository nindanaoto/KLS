#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "spral_scaling.h"

int main(int argc, char **argv) {
  if (argc != 3) { fprintf(stderr, "usage: %s in.mtx out.mtx\n", argv[0]); return 1; }
  FILE *f = fopen(argv[1], "r");
  if (!f) { perror("open"); return 1; }
  char line[256];
  long n = 0, m = 0, nnz = 0;
  while (fgets(line, sizeof line, f)) {
    if (line[0] == '%') continue;
    sscanf(line, "%ld %ld %ld", &n, &m, &nnz);
    break;
  }
  long *ri = malloc(nnz * sizeof(long));
  long *ci = malloc(nnz * sizeof(long));
  double *vx = malloc(nnz * sizeof(double));
  char **raw = malloc(nnz * sizeof(char *));
  for (long p = 0; p < nnz; ++p) {
    long i, j; char valbuf[64];
    if (fscanf(f, "%ld %ld %63s", &i, &j, valbuf) != 3) { fprintf(stderr, "read fail %ld\n", p); return 1; }
    ri[p] = i - 1; ci[p] = j - 1; vx[p] = atof(valbuf); raw[p] = strdup(valbuf);
  }
  fclose(f);
  /* CSC build */
  int64_t *cp = calloc(n + 1, sizeof(int64_t));
  for (long p = 0; p < nnz; ++p) cp[ci[p] + 1]++;
  for (long j = 0; j < n; ++j) cp[j + 1] += cp[j];
  int *rows = malloc(nnz * sizeof(int));
  double *vals = malloc(nnz * sizeof(double));
  long *pos = malloc(nnz * sizeof(long));
  int64_t *cur = malloc((n + 1) * sizeof(int64_t));
  memcpy(cur, cp, (n + 1) * sizeof(int64_t));
  for (long p = 0; p < nnz; ++p) {
    int64_t d = cur[ci[p]]++;
    rows[d] = (int)ri[p]; vals[d] = vx[p]; pos[d] = p;
  }
  double *rs = malloc(n * sizeof(double));
  double *cs = malloc(n * sizeof(double));
  int *match = malloc(n * sizeof(int));
  struct spral_scaling_hungarian_options opt;
  struct spral_scaling_hungarian_inform inf;
  spral_scaling_hungarian_default_options(&opt);
  opt.array_base = 0;
  opt.scale_if_singular = false;
  spral_scaling_hungarian_unsym_long((int)n, (int)n, cp, rows, vals, rs, cs, match, &opt, &inf);
  if (inf.flag != 0 && inf.flag != -2) { fprintf(stderr, "spral flag %d\n", inf.flag); return 1; }
  /* match[row] = col; row perm: row goes to position match[row] */
  long *inv_row = malloc(n * sizeof(long));
  long matched = 0;
  for (long r = 0; r < n; ++r) inv_row[r] = -1;
  for (long r = 0; r < n; ++r) {
    if (match[r] >= 0 && match[r] < n) { inv_row[r] = match[r]; matched++; }
  }
  fprintf(stderr, "matched %ld/%ld\n", matched, n);
  if (matched != n) {
    /* fill unmatched rows into unmatched positions */
    char *taken = calloc(n, 1);
    for (long r = 0; r < n; ++r) if (inv_row[r] >= 0) taken[inv_row[r]] = 1;
    long next = 0;
    for (long r = 0; r < n; ++r) {
      if (inv_row[r] >= 0) continue;
      while (taken[next]) next++;
      inv_row[r] = next; taken[next] = 1;
    }
    free(taken);
  }
  FILE *g = fopen(argv[2], "w");
  fprintf(g, "%%%%MatrixMarket matrix coordinate real general\n%ld %ld %ld\n", n, n, nnz);
  for (long p = 0; p < nnz; ++p) {
    fprintf(g, "%ld %ld %s\n", inv_row[ri[p]] + 1, ci[p] + 1, raw[p]);
  }
  fclose(g);
  fprintf(stderr, "done\n");
  return 0;
}
