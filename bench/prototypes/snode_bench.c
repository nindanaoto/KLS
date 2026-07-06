/* Prototype: sorted-storage column-supernode refactor vs klu_l_refactor.
 * Factor once (scale=0), klu_l_sort, detect strict nested-tail runs, then
 * benchmark a hand-rolled left-looking refactor with panel batching against
 * the stock KLU refactor.  Verifies Udiag equality. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "trilinos_klu_decl.h"

typedef long I;
#define SN_MAX_BATCH 64
#define SN_MIN_BATCH 3

static double now_sec(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

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

/* GET_POINTER equivalent for UF_long build: indices then values */
static inline void get_ptr(double *lu, const I *xip, const I *xlen, I k,
                           I **xi, double **xx, I *len) {
  double *xp = lu + xip[k];
  *len = xlen[k];
  *xi = (I *)xp;
  *xx = xp + *len; /* UNITS(Int,len) == len units when Int is 8 bytes */
}

/* run panel registry: per run, a row-major [tlen x w] mirror of the run
 * columns' shared-tail values. run_of[gcol] = run index or -1. */
typedef struct {
  I gs, ge;      /* global col range [gs, ge) */
  I tlen;        /* shared tail length (pattern of col ge-1) */
  double *P;     /* [tlen * w] row-major, w = ge - gs */
} run_panel;
static run_panel *g_panels = NULL;
static I g_npanels = 0;
static I *g_run_of = NULL;

static void build_panel_registry(I n, trilinos_klu_l_symbolic *sym,
                                 trilinos_klu_l_numeric *num,
                                 const I *run_end) {
  g_run_of = malloc((size_t)n * sizeof(I));
  for (I i = 0; i < n; ++i) g_run_of[i] = -1;
  I cap = 1024;
  g_panels = malloc((size_t)cap * sizeof(run_panel));
  g_npanels = 0;
  const I *R = sym->R;
  for (I b = 0; b < sym->nblocks; ++b) {
    I k1 = R[b], k2 = R[b + 1];
    I g = k1;
    while (g < k2) {
      I e = run_end[g];
      if (e <= g + 1) { g++; continue; }
      I *Lip = num->Lip + k1;
      I *Llen = num->Llen + k1;
      double *LU = (double *)num->LUbx[b];
      I *li; double *lx; I ll;
      get_ptr(LU, Lip, Llen, e - 1 - k1, &li, &lx, &ll);
      I w = e - g;
      if (g_npanels == cap) { cap *= 2; g_panels = realloc(g_panels, (size_t)cap * sizeof(run_panel)); }
      run_panel *rp = &g_panels[g_npanels];
      rp->gs = g; rp->ge = e; rp->tlen = ll;
      rp->P = malloc((size_t)ll * (size_t)w * sizeof(double));
      for (I c = g; c < e; ++c) g_run_of[c] = g_npanels;
      g_npanels++;
      g = e;
    }
  }
}

/* custom refactor over sorted numeric with panel batching.
 * run_end[gcol] = exclusive end of the strict run containing gcol, or 0. */
static int my_refactor(I n, const I *Ap, const I *Ai, const double *Az,
                       trilinos_klu_l_symbolic *sym, trilinos_klu_l_numeric *num,
                       const I *run_end, double *X, int use_panels,
                       double *panel_flops, double *scalar_flops) {
  const I *Q = sym->Q;
  const I *R = sym->R;
  const I *Pinv = num->Pinv;
  const I *Offp = num->Offp;
  double *Offx = (double *)num->Offx;
  double *Udiag = (double *)num->Udiag;
  I nblocks = sym->nblocks;
  I poff = 0;
  I nzoff = num->nzoff;
  double xs[SN_MAX_BATCH];
  I *li_arr[SN_MAX_BATCH];
  double *lx_arr[SN_MAX_BATCH];

  for (I block = 0; block < nblocks; ++block) {
    I k1 = R[block], k2 = R[block + 1], nk = k2 - k1;
    if (nk == 1) {
      I oldcol = Q[k1];
      double s = 0.0;
      for (I p = Ap[oldcol]; p < Ap[oldcol + 1]; ++p) {
        I newrow = Pinv[Ai[p]] - k1;
        if (newrow < 0 && poff < nzoff) Offx[poff++] = Az[p];
        else s = Az[p];
      }
      Udiag[k1] = s;
      if (s == 0.0) return 0;
      continue;
    }
    I *Lip = num->Lip + k1;
    I *Llen = num->Llen + k1;
    I *Uip = num->Uip + k1;
    I *Ulen = num->Ulen + k1;
    double *LU = (double *)num->LUbx[block];
    for (I k = 0; k < nk; ++k) {
      I oldcol = Q[k + k1];
      for (I p = Ap[oldcol]; p < Ap[oldcol + 1]; ++p) {
        I newrow = Pinv[Ai[p]] - k1;
        if (newrow < 0 && poff < nzoff) Offx[poff++] = Az[p];
        else X[newrow] = Az[p];
      }
      I *Ui; double *Ux; I ulen;
      get_ptr(LU, Uip, Ulen, k, &Ui, &Ux, &ulen);
      I up = 0;
      while (up < ulen) {
        I j = Ui[up];
        if (use_panels) {
          I g = j + k1;
          I end = run_end[g];
          if (end > g + 1) {
            I tmax = end - g;
            if (tmax > SN_MAX_BATCH) tmax = SN_MAX_BATCH;
            if (tmax > ulen - up) tmax = ulen - up;
            I t = 1;
            while (t < tmax && Ui[up + t] == j + t) t++;
            if (t >= SN_MIN_BATCH) {
              for (I i = 0; i < t; ++i) {
                I dummy_len;
                get_ptr(LU, Lip, Llen, j + i, &li_arr[i], &lx_arr[i], &dummy_len);
                xs[i] = X[j + i];
                X[j + i] = 0.0;
              }
              /* in-panel unit-lower updates + U writes */
              for (I i = 0; i < t; ++i) {
                double u = xs[i];
                Ux[up + i] = u;
                const double *lxi = lx_arr[i];
                for (I r0 = 0; r0 < t - 1 - i; ++r0)
                  xs[i + 1 + r0] -= lxi[r0] * u;
              }
              /* shared extended tail: pattern of col j+t-1.
                 Chunked rows, column-inner contiguous streams so the
                 compiler vectorizes and the accumulators stay independent. */
              I *tli; double *tlx; I tlen;
              get_ptr(LU, Lip, Llen, j + t - 1, &tli, &tlx, &tlen);
              const I gj = j + k1;
              if (use_panels == 2 && gj + t == run_end[gj] &&
                  g_run_of[gj] >= 0) {
                const run_panel *rp = &g_panels[g_run_of[gj]];
                const I w = rp->ge - rp->gs;
                const I c0 = gj - rp->gs;
                const double *P = rp->P + c0;
                I p = 0;
                for (; p + 4 <= tlen; p += 4) {
                  const double *r0 = P + (size_t)p * (size_t)w;
                  const double *r1 = r0 + w;
                  const double *r2 = r1 + w;
                  const double *r3 = r2 + w;
                  double a0 = 0, a1 = 0, a2 = 0, a3 = 0;
                  for (I i = 0; i < t; ++i) {
                    const double xv = xs[i];
                    a0 += r0[i] * xv;
                    a1 += r1[i] * xv;
                    a2 += r2[i] * xv;
                    a3 += r3[i] * xv;
                  }
                  X[tli[p]] -= a0;
                  X[tli[p + 1]] -= a1;
                  X[tli[p + 2]] -= a2;
                  X[tli[p + 3]] -= a3;
                }
                for (; p < tlen; ++p) {
                  const double *row = P + (size_t)p * (size_t)w;
                  double acc = 0.0;
                  for (I i = 0; i < t; ++i) acc += row[i] * xs[i];
                  X[tli[p]] -= acc;
                }
              } else {
              for (I p0 = 0; p0 < tlen; p0 += 32) {
                I pc = tlen - p0 < 32 ? tlen - p0 : 32;
                double acc[32];
                {
                  const double *src = tlx + p0;
                  const double u = xs[t - 1];
                  for (I p = 0; p < pc; ++p) acc[p] = src[p] * u;
                }
                for (I i = 0; i < t - 1; ++i) {
                  const double *src = lx_arr[i] + (t - 1 - i) + p0;
                  const double u = xs[i];
                  for (I p = 0; p < pc; ++p) acc[p] += src[p] * u;
                }
                for (I p = 0; p < pc; ++p) X[tli[p0 + p]] -= acc[p];
              }
              }
              *panel_flops += 2.0 * (double)t * (double)tlen +
                              (double)t * (double)(t - 1);
              up += t;
              continue;
            }
          }
        }
        double ujk = X[j];
        X[j] = 0.0;
        Ux[up] = ujk;
        I *Li; double *Lx; I llen;
        get_ptr(LU, Lip, Llen, j, &Li, &Lx, &llen);
        if (ujk != 0.0) {
          for (I p = 0; p < llen; ++p) X[Li[p]] -= Lx[p] * ujk;
          *scalar_flops += 2.0 * (double)llen;
        }
        up++;
      }
      double ukk = X[k];
      X[k] = 0.0;
      if (ukk == 0.0) return 0;
      Udiag[k + k1] = ukk;
      I *Li; double *Lx; I llen;
      get_ptr(LU, Lip, Llen, k, &Li, &Lx, &llen);
      for (I p = 0; p < llen; ++p) {
        I i = Li[p];
        Lx[p] = X[i] / ukk;
        X[i] = 0.0;
      }
      if (use_panels == 2) {
        const I gk = k + k1;
        const I ridx = g_run_of != NULL ? g_run_of[gk] : -1;
        if (ridx >= 0) {
          const run_panel *rp = &g_panels[ridx];
          const I w = rp->ge - rp->gs;
          const I c = gk - rp->gs;
          double *P = rp->P + c;
          const double *src = Lx + (llen - rp->tlen);
          for (I p = 0; p < rp->tlen; ++p) P[(size_t)p * (size_t)w] = src[p];
        }
      }
    }
  }
  return 1;
}

int main(int argc, char **argv) {
  I n, *cp, *ri; double *vx;
  int reps = argc > 2 ? atoi(argv[2]) : 20;
  if (!read_mtx(argv[1], &n, &cp, &ri, &vx)) { fprintf(stderr, "read failed\n"); return 1; }
  trilinos_klu_l_common common;
  trilinos_klu_l_defaults(&common);
  common.scale = 0;
  trilinos_klu_l_symbolic *sym = trilinos_klu_l_analyze(n, cp, ri, &common);
  trilinos_klu_l_numeric *num = trilinos_klu_l_factor(cp, ri, vx, sym, &common);
  if (!num) { fprintf(stderr, "factor failed\n"); return 1; }
  trilinos_klu_l_flops(sym, num, &common);
  printf("n=%ld lnz=%ld unz=%ld flops=%.3e nblocks=%ld\n", (long)n,
         (long)num->lnz, (long)num->unz, common.flops, (long)sym->nblocks);
  if (!trilinos_klu_l_sort(sym, num, &common)) { fprintf(stderr, "sort failed\n"); return 1; }

  /* detect strict runs on sorted storage */
  I *run_end = calloc((size_t)n, sizeof(I));
  {
    const I *R = sym->R;
    I runs = 0, cols_in = 0;
    double entries_in = 0, entries_total = 0;
    for (I b = 0; b < sym->nblocks; ++b) {
      I k1 = R[b], k2 = R[b + 1], nk = k2 - k1;
      if (nk < 2) continue;
      I *Lip = num->Lip + k1;
      I *Llen = num->Llen + k1;
      double *LU = (double *)num->LUbx[b];
      I start = 0;
      for (I k = 0; k < nk; ++k) {
        entries_total += (double)Llen[k];
        int extends = 0;
        if (k + 1 < nk && Llen[k + 1] == Llen[k] - 1 && Llen[k] >= 1) {
          I *li, *li2; double *lx, *lx2; I l1, l2;
          get_ptr(LU, Lip, Llen, k, &li, &lx, &l1);
          get_ptr(LU, Lip, Llen, k + 1, &li2, &lx2, &l2);
          if (li[0] == k + 1 &&
              memcmp(li + 1, li2, (size_t)l2 * sizeof(I)) == 0)
            extends = 1;
        }
        if (!extends) {
          if (k > start) {
            runs++;
            for (I c = start; c <= k; ++c) {
              run_end[k1 + c] = k1 + k + 1;
              I *li; double *lx; I ll;
              get_ptr(LU, Lip, Llen, c, &li, &lx, &ll);
              entries_in += (double)ll;
            }
            cols_in += k - start + 1;
          }
          start = k + 1;
        }
      }
    }
    printf("runs=%ld cols_in=%ld/%ld entries_in=%.0f/%.0f (%.1f%%)\n",
           (long)runs, (long)cols_in, (long)n, entries_in, entries_total,
           100.0 * entries_in / (entries_total > 0 ? entries_total : 1));
  }

  double *X = calloc((size_t)sym->maxblock, sizeof(double));
  double pf = 0, sf = 0;

  /* verify: run both once, compare Udiag */
  double *ud_ref = malloc((size_t)n * sizeof(double));
  trilinos_klu_l_refactor(cp, ri, vx, sym, num, &common);
  memcpy(ud_ref, num->Udiag, (size_t)n * sizeof(double));
  if (!my_refactor(n, cp, ri, vx, sym, num, run_end, X, 1, &pf, &sf)) {
    fprintf(stderr, "my_refactor singular/fail\n");
    return 1;
  }
  double maxrel = 0;
  double *ud = (double *)num->Udiag;
  for (I i = 0; i < n; ++i) {
    double d = fabs(ud[i] - ud_ref[i]);
    double s = fabs(ud_ref[i]) > 1e-300 ? fabs(ud_ref[i]) : 1.0;
    if (d / s > maxrel) maxrel = d / s;
  }
  printf("udiag max rel diff vs klu: %.2e  panel_flops=%.3e scalar_flops=%.3e (panel share %.1f%%)\n",
         maxrel, pf, sf, 100.0 * pf / (pf + sf > 0 ? pf + sf : 1));

  double t0, t1;
  t0 = now_sec();
  for (int r = 0; r < reps; ++r) trilinos_klu_l_refactor(cp, ri, vx, sym, num, &common);
  t1 = now_sec();
  double klu_ms = (t1 - t0) * 1000.0 / reps;
  t0 = now_sec();
  for (int r = 0; r < reps; ++r) my_refactor(n, cp, ri, vx, sym, num, run_end, X, 0, &pf, &sf);
  t1 = now_sec();
  double scalar_ms = (t1 - t0) * 1000.0 / reps;
  t0 = now_sec();
  for (int r = 0; r < reps; ++r) my_refactor(n, cp, ri, vx, sym, num, run_end, X, 1, &pf, &sf);
  t1 = now_sec();
  double panel_ms = (t1 - t0) * 1000.0 / reps;
  build_panel_registry(n, sym, num, run_end);
  /* verify mode 2 too */
  my_refactor(n, cp, ri, vx, sym, num, run_end, X, 2, &pf, &sf);
  {
    double mr = 0; double *ud2 = (double *)num->Udiag;
    for (I i = 0; i < n; ++i) {
      double d = fabs(ud2[i] - ud_ref[i]);
      double s = fabs(ud_ref[i]) > 1e-300 ? fabs(ud_ref[i]) : 1.0;
      if (d / s > mr) mr = d / s;
    }
    printf("mode2 udiag max rel diff: %.2e (npanels=%ld)\n", mr, (long)g_npanels);
  }
  t0 = now_sec();
  for (int r = 0; r < reps; ++r) my_refactor(n, cp, ri, vx, sym, num, run_end, X, 2, &pf, &sf);
  t1 = now_sec();
  double dense_ms = (t1 - t0) * 1000.0 / reps;
  printf("refactor ms/iter: klu=%.2f  mine-scalar=%.2f  mine-panel=%.2f  mine-densepanel=%.2f (dense vs panel %.2fx, vs klu %.2fx)\n",
         klu_ms, scalar_ms, panel_ms, dense_ms, panel_ms / dense_ms, klu_ms / dense_ms);
  return 0;
}
