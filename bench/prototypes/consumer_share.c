/* Measure: within each strict consumer run, what fraction of producer-batch
 * flops is shared by ALL members (GEMM-able) vs member-private. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "trilinos_klu_decl.h"
typedef long I;
static inline void get_ptr(double *lu, const I *xip, const I *xlen, I k,
                           I **xi, double **xx, I *len) {
  double *xp = lu + xip[k];
  *len = xlen[k];
  *xi = (I *)xp;
  *xx = xp + *len;
}

#include "read_mtx.inc"
int main(int argc, char **argv) {
  I n, *cp, *ri; double *vx;
  if (!read_mtx(argv[1], &n, &cp, &ri, &vx)) { fprintf(stderr, "read fail\n"); return 1; }
  trilinos_klu_l_common common;
  trilinos_klu_l_defaults(&common);
  common.scale = 0;
  trilinos_klu_l_symbolic *sym = trilinos_klu_l_analyze(n, cp, ri, &common);
  trilinos_klu_l_numeric *num = trilinos_klu_l_factor(cp, ri, vx, sym, &common);
  if (!num) return 1;
  trilinos_klu_l_sort(sym, num, &common);
  const I *R = sym->R;
  double shared_fl = 0, priv_fl = 0, scalar_fl = 0;
  double gemm_fl = 0;
  for (I b = 0; b < sym->nblocks; ++b) {
    I k1 = R[b], k2 = R[b + 1], nk = k2 - k1;
    if (nk < 2) continue;
    I *Lip = num->Lip + k1; I *Llen = num->Llen + k1;
    I *Uip = num->Uip + k1; I *Ulen = num->Ulen + k1;
    double *LU = (double *)num->LUbx[b];
    /* strict run detection (consumer side) */
    I *run_end = calloc((size_t)nk, sizeof(I));
    I start = 0;
    for (I k = 0; k < nk; ++k) {
      int ext = 0;
      if (k + 1 < nk && Llen[k + 1] == Llen[k] - 1 && Llen[k] >= 1) {
        I *li, *li2; double *lx, *lx2; I l1, l2;
        get_ptr(LU, Lip, Llen, k, &li, &lx, &l1);
        get_ptr(LU, Lip, Llen, k + 1, &li2, &lx2, &l2);
        if (li[0] == k + 1 && memcmp(li + 1, li2, (size_t)l2 * sizeof(I)) == 0) ext = 1;
      }
      if (!ext) {
        if (k > start) for (I c = start; c <= k; ++c) run_end[c] = k + 1;
        start = k + 1;
      }
    }
    /* per consumer run: batches of first member vs the rest */
    I g = 0;
    while (g < nk) {
      I e = run_end[g];
      I wc = (e > g + 1) ? e - g : 1;
      /* count producer batches per member; a batch (j,t) is "shared" if the
         SAME j appears in every member's U list with t reaching >= same end */
      /* simple approximation: intersect the first member's batch set with
         each other member's batch starts */
      if (wc == 1) {
        I *ui; double *ux; I ul;
        get_ptr(LU, Uip, Ulen, g, &ui, &ux, &ul);
        for (I p = 0; p < ul; ++p) {
          I *li; double *lx; I ll;
          get_ptr(LU, Lip, Llen, ui[p], &li, &lx, &ll);
          scalar_fl += 2.0 * ll;
        }
        g++;
        continue;
      }
      /* mark first member's producers */
      static I *mark = NULL; static I marked_n = 0;
      if (mark == NULL) { mark = calloc((size_t)n, sizeof(I)); marked_n = 1; }
      I *ui0; double *ux0; I ul0;
      get_ptr(LU, Uip, Ulen, g, &ui0, &ux0, &ul0);
      for (I p = 0; p < ul0; ++p) mark[ui0[p]] = g + 1;
      for (I m = g; m < e; ++m) {
        I *ui; double *ux; I ul;
        get_ptr(LU, Uip, Ulen, m, &ui, &ux, &ul);
        for (I p = 0; p < ul; ++p) {
          I j = ui[p];
          if (j >= g) continue; /* in-run triangle */
          I *li; double *lx; I ll;
          get_ptr(LU, Lip, Llen, j, &li, &lx, &ll);
          if (mark[j] == g + 1) shared_fl += 2.0 * ll;
          else priv_fl += 2.0 * ll;
        }
      }
      g = e;
    }
    free(run_end);
  }
  double tot = shared_fl + priv_fl + scalar_fl;
  printf("%s: shared(GEMMable)=%.3e (%.1f%%) private=%.3e (%.1f%%) scalar=%.3e (%.1f%%)\n",
         argv[1], shared_fl, 100*shared_fl/tot, priv_fl, 100*priv_fl/tot, scalar_fl, 100*scalar_fl/tot);
  return 0;
}
