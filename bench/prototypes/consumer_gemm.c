/* Prototype: consumer-run supernodal left-looking refactor.
 * For each strict consumer run [g,e): walk the UNION of member U patterns
 * ascending; at each producer j read u_m = X_m[j] (structural zeros are
 * naturally 0.0), zero it, store to member U storage where patterned, and
 * update all member workspaces from L(:,j) in one pass (row value reused
 * wc times).  In-run producers finalize member pivots in order first.
 * Compare against the per-column scalar/batched refactor. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "trilinos_klu_decl.h"
typedef long I;
#define SN_MAX_BATCH 64
#define SN_MIN_BATCH 3
#define WC_CAP 32
static double now_sec(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + 1e-9 * ts.tv_nsec;
}
static inline void get_ptr(double *lu, const I *xip, const I *xlen, I k,
                           I **xi, double **xx, I *len) {
  double *xp = lu + xip[k];
  *len = xlen[k];
  *xi = (I *)xp;
  *xx = xp + *len;
}
#include "read_mtx.inc"

/* per consumer run precomputed union: sorted producer list (external only,
 * j < g), and for each member a cursor-mapped store position array is not
 * needed: member U patterns are sorted so walk cursors suffice. */
typedef struct {
  I g, e;          /* local col range in block */
  I ulen;          /* union external producer count */
  I *uj;           /* union producer local indices, ascending */
} crun;

static I *g_rpos = NULL; /* local row -> compact panel row, per active run */

static int refactor_runs(I n, const I *Ap, const I *Ai, const double *Az,
                         trilinos_klu_l_symbolic *sym,
                         trilinos_klu_l_numeric *num,
                         const I *run_end_all, double *W, I ldw,
                         crun **cruns_by_block, I *ncruns_by_block,
                         double *X) {
  const I *Q = sym->Q;
  const I *R = sym->R;
  const I *Pinv = num->Pinv;
  const I *Offp = num->Offp;
  double *Offx = (double *)num->Offx;
  double *Udiag = (double *)num->Udiag;
  I nblocks = sym->nblocks;
  I poff = 0;
  I nzoff = num->nzoff;
  (void)Offp;
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
    crun *cr = cruns_by_block[block];
    I ncr = ncruns_by_block[block];
    I cri = 0;
    I k = 0;
    while (k < nk) {
      if (cri < ncr && cr[cri].g == k) {
        /* joint consumer-run processing, compact-row panel:
           rows = union producers ++ [g..e) ++ tail(pattern of col g). */
        const crun *c = &cr[cri];
        const I g = c->g, e = c->e, wc = e - g;
        I nrows = 0;
        for (I p = 0; p < c->ulen; ++p) {
          g_rpos[c->uj[p]] = nrows++;
        }
        for (I m = g; m < e; ++m) {
          g_rpos[m] = nrows++;
        }
        {
          I *li; double *lx; I ll;
          get_ptr(LU, Lip, Llen, g, &li, &lx, &ll);
          for (I q = 0; q < ll; ++q) {
            if (li[q] >= e) {
              g_rpos[li[q]] = nrows++;
            }
          }
        }
        /* W reused as compact panel [nrows][wc], row-major */
        memset(W, 0, (size_t)nrows * (size_t)wc * sizeof(double));
        for (I m = 0; m < wc; ++m) {
          I oldcol = Q[g + m + k1];
          for (I p = Ap[oldcol]; p < Ap[oldcol + 1]; ++p) {
            I newrow = Pinv[Ai[p]] - k1;
            if (newrow < 0 && poff < nzoff) Offx[poff++] = Az[p];
            else W[(size_t)g_rpos[newrow] * wc + m] = Az[p];
          }
        }
        /* external union producers (all < g), with producer-run batching */
        const I *uj = c->uj;
        const I ul = c->ulen;
        double us[SN_MAX_BATCH][WC_CAP];
        I up = 0;
        while (up < ul) {
          I j = uj[up];
          I gj = j + k1;
          I t = 1;
          I re = run_end_all[gj];
          if (re > gj + 1) {
            I tmax = re - gj;
            if (tmax > SN_MAX_BATCH) tmax = SN_MAX_BATCH;
            if (tmax > ul - up) tmax = ul - up;
            while (t < tmax && uj[up + t] == j + t) t++;
          }
          /* collect u values for the t x wc tile */
          for (I i = 0; i < t; ++i) {
            double *row = W + (size_t)g_rpos[j + i] * wc;
            for (I m = 0; m < wc; ++m) {
              us[i][m] = row[m];
              row[m] = 0.0;
            }
          }
          /* in-batch unit-lower propagation (producer triangle) */
          for (I i = 0; i < t; ++i) {
            I *li; double *lx; I ll;
            get_ptr(LU, Lip, Llen, j + i, &li, &lx, &ll);
            for (I r0 = 0; r0 + i + 1 < t; ++r0) {
              const double lv = lx[r0];
              for (I m = 0; m < wc; ++m) us[i + 1 + r0][m] -= lv * us[i][m];
            }
          }
          /* shared tail: per member, the proven chunked column-stream loop;
             producer columns stay cache-resident across members so memory
             reads amortize over the whole consumer run */
          I *tli; double *tlx; I tlen;
          get_ptr(LU, Lip, Llen, j + t - 1, &tli, &tlx, &tlen);
          {
            const double *lx_arr[SN_MAX_BATCH];
            for (I i = 0; i < t - 1; ++i) {
              I *li0; double *lx0; I ll0;
              get_ptr(LU, Lip, Llen, j + i, &li0, &lx0, &ll0);
              lx_arr[i] = lx0 + (t - 1 - i);
            }
            /* outer-product microkernel: per tail row, vector FMA across
               the wc consumers */
            for (I p = 0; p < tlen; ++p) {
              double *wrow = W + (size_t)g_rpos[tli[p]] * wc;
              {
                const double lv = tlx[p];
                const double *u = us[t - 1];
                for (I m = 0; m < wc; ++m) wrow[m] -= lv * u[m];
              }
              for (I i = 0; i < t - 1; ++i) {
                const double lv = lx_arr[i][p];
                if (lv == 0.0) continue;
                const double *u = us[i];
                for (I m = 0; m < wc; ++m) wrow[m] -= lv * u[m];
              }
            }
          }
          /* store u values into member U columns where patterned */
          for (I m = 0; m < wc; ++m) {
            I *ui; double *ux; I uln;
            get_ptr(LU, Uip, Ulen, g + m, &ui, &ux, &uln);
            /* member patterns sorted: binary search once per batch */
            I lo = 0, hi = uln;
            while (lo < hi) { I mid = (lo + hi) / 2; if (ui[mid] < j) lo = mid + 1; else hi = mid; }
            for (I i = 0; i < t && lo < uln; ++i) {
              if (ui[lo] == j + i) { ux[lo] = us[i][m]; lo++; }
            }
          }
          up += t;
        }
        /* in-run: finalize members in order; consume within run through
           the compact panel (column m of W) */
        for (I m = 0; m < wc; ++m) {
          I *ui; double *ux; I uln;
          get_ptr(LU, Uip, Ulen, g + m, &ui, &ux, &uln);
          for (I p = 0; p < uln; ++p) {
            I j = ui[p];
            if (j < g) continue;
            double *urow = W + (size_t)g_rpos[j] * wc;
            double u = urow[m];
            urow[m] = 0.0;
            ux[p] = u;
            if (u != 0.0) {
              I *li; double *lx; I ll;
              get_ptr(LU, Lip, Llen, j, &li, &lx, &ll);
              for (I q = 0; q < ll; ++q) {
                W[(size_t)g_rpos[li[q]] * wc + m] -= lx[q] * u;
              }
            }
          }
          double *drow = W + (size_t)g_rpos[g + m] * wc;
          double ukk = drow[m];
          drow[m] = 0.0;
          if (ukk == 0.0) return 0;
          Udiag[g + m + k1] = ukk;
          I *li; double *lx; I ll;
          get_ptr(LU, Lip, Llen, g + m, &li, &lx, &ll);
          for (I q = 0; q < ll; ++q) {
            double *rrow = W + (size_t)g_rpos[li[q]] * wc;
            lx[q] = rrow[m] / ukk;
            rrow[m] = 0.0;
          }
        }
        k = e;
        cri++;
        continue;
      }
      /* plain scalar column */
      I oldcol = Q[k + k1];
      for (I p = Ap[oldcol]; p < Ap[oldcol + 1]; ++p) {
        I newrow = Pinv[Ai[p]] - k1;
        if (newrow < 0 && poff < nzoff) Offx[poff++] = Az[p];
        else X[newrow] = Az[p];
      }
      I *Ui; double *Ux; I ulen;
      get_ptr(LU, Uip, Ulen, k, &Ui, &Ux, &ulen);
      for (I p = 0; p < ulen; ++p) {
        I j = Ui[p];
        double ujk = X[j];
        X[j] = 0.0;
        Ux[p] = ujk;
        if (ujk != 0.0) {
          I *Li; double *Lx; I llen;
          get_ptr(LU, Lip, Llen, j, &Li, &Lx, &llen);
          for (I q = 0; q < llen; ++q) X[Li[q]] -= Lx[q] * ujk;
        }
      }
      double ukk = X[k];
      X[k] = 0.0;
      if (ukk == 0.0) return 0;
      Udiag[k + k1] = ukk;
      I *Li; double *Lx; I llen;
      get_ptr(LU, Lip, Llen, k, &Li, &Lx, &llen);
      for (I p = 0; p < llen; ++p) {
        I r = Li[p];
        Lx[p] = X[r] / ukk;
        X[r] = 0.0;
      }
      k++;
    }
  }
  return 1;
}

int main(int argc, char **argv) {
  I n, *cp, *ri; double *vx;
  int reps = argc > 2 ? atoi(argv[2]) : 10;
  if (!read_mtx(argv[1], &n, &cp, &ri, &vx)) { fprintf(stderr, "read fail\n"); return 1; }
  trilinos_klu_l_common common;
  trilinos_klu_l_defaults(&common);
  common.scale = 0;
  trilinos_klu_l_symbolic *sym = trilinos_klu_l_analyze(n, cp, ri, &common);
  trilinos_klu_l_numeric *num = trilinos_klu_l_factor(cp, ri, vx, sym, &common);
  if (!num) { fprintf(stderr, "factor fail\n"); return 1; }
  trilinos_klu_l_sort(sym, num, &common);
  const I *R = sym->R;
  /* run detection (global cols) */
  I *run_end = calloc((size_t)n, sizeof(I));
  for (I b = 0; b < sym->nblocks; ++b) {
    I k1 = R[b], k2 = R[b + 1], nk = k2 - k1;
    if (nk < 2) continue;
    I *Lip = num->Lip + k1; I *Llen = num->Llen + k1;
    double *LU = (double *)num->LUbx[b];
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
        if (k > start) for (I c = start; c <= k; ++c) run_end[k1 + c] = k1 + k + 1;
        start = k + 1;
      }
    }
  }
  /* consumer runs + union producer lists per block (width cap WC_CAP) */
  crun **cbb = calloc((size_t)sym->nblocks, sizeof(crun *));
  I *ncb = calloc((size_t)sym->nblocks, sizeof(I));
  I *seen = calloc((size_t)n, sizeof(I));
  I stamp = 0;
  I *tmp = malloc((size_t)n * sizeof(I));
  for (I b = 0; b < sym->nblocks; ++b) {
    I k1 = R[b], k2 = R[b + 1], nk = k2 - k1;
    if (nk < 2) continue;
    I cap = 64, cnt = 0;
    crun *arr = malloc((size_t)cap * sizeof(crun));
    I *Uip = num->Uip + k1; I *Ulen = num->Ulen + k1;
    double *LU = (double *)num->LUbx[b];
    I g = 0;
    while (g < nk) {
      I e = run_end[k1 + g];
      if (e <= k1 + g + 1) { g++; continue; }
      e -= k1;
      /* split wide runs into WC_CAP chunks */
      for (I gs = g; gs < e; gs += WC_CAP) {
        I ge = gs + WC_CAP < e ? gs + WC_CAP : e;
        stamp++;
        I un = 0;
        for (I m = gs; m < ge; ++m) {
          I *ui; double *ux; I ul;
          get_ptr(LU, Uip, Ulen, m, &ui, &ux, &ul);
          for (I p = 0; p < ul; ++p) {
            I j = ui[p];
            if (j >= gs) break; /* in-run (sorted) */
            if (seen[j] != stamp) { seen[j] = stamp; tmp[un++] = j; }
          }
        }
        /* sort union (small; insertion) */
        for (I a = 1; a < un; ++a) {
          I v = tmp[a]; I bq = a;
          while (bq > 0 && tmp[bq - 1] > v) { tmp[bq] = tmp[bq - 1]; bq--; }
          tmp[bq] = v;
        }
        if (cnt == cap) { cap *= 2; arr = realloc(arr, (size_t)cap * sizeof(crun)); }
        arr[cnt].g = gs; arr[cnt].e = ge; arr[cnt].ulen = un;
        arr[cnt].uj = malloc((size_t)(un > 0 ? un : 1) * sizeof(I));
        memcpy(arr[cnt].uj, tmp, (size_t)un * sizeof(I));
        cnt++;
      }
      g = e;
    }
    cbb[b] = arr; ncb[b] = cnt;
  }
  double *X = calloc((size_t)sym->maxblock, sizeof(double));
  double *W = calloc((size_t)sym->maxblock * WC_CAP, sizeof(double));
  g_rpos = malloc((size_t)sym->maxblock * sizeof(I));
  /* verify */
  double *ud_ref = malloc((size_t)n * sizeof(double));
  trilinos_klu_l_refactor(cp, ri, vx, sym, num, &common);
  memcpy(ud_ref, num->Udiag, (size_t)n * sizeof(double));
  if (!refactor_runs(n, cp, ri, vx, sym, num, run_end, W, sym->maxblock, cbb, ncb, X)) {
    fprintf(stderr, "joint refactor fail\n"); return 1;
  }
  double mr = 0; double *ud = (double *)num->Udiag;
  for (I i = 0; i < n; ++i) {
    double d = fabs(ud[i] - ud_ref[i]);
    double s = fabs(ud_ref[i]) > 1e-300 ? fabs(ud_ref[i]) : 1.0;
    if (d / s > mr) mr = d / s;
  }
  printf("joint udiag max rel diff: %.2e\n", mr);
  double t0 = now_sec();
  for (int r = 0; r < reps; ++r) trilinos_klu_l_refactor(cp, ri, vx, sym, num, &common);
  double t1 = now_sec();
  double klu_ms = (t1 - t0) * 1000.0 / reps;
  t0 = now_sec();
  for (int r = 0; r < reps; ++r) refactor_runs(n, cp, ri, vx, sym, num, run_end, W, sym->maxblock, cbb, ncb, X);
  t1 = now_sec();
  double joint_ms = (t1 - t0) * 1000.0 / reps;
  printf("refactor ms/iter: klu=%.2f joint=%.2f (%.2fx vs klu)\n", klu_ms, joint_ms, klu_ms / joint_ms);
  return 0;
}
