/* supernodal_bench: prototype harness for the supernodal dense-panel
 * refactor storage rebuild (task #12).
 *
 * Reconstruction of the lost scratchpad prototype (rounds 1-2), kept in
 * git this time.  Three refactor kernels run interleaved on the same
 * analyzed+factored pattern:
 *
 *   klu    - vendored trilinos_klu_l_refactor (scattered column storage)
 *   landed - replica of the landed KLS serial refactor kernel: precomputed
 *            A-scatter program, i32 L-index mirrors, prefetching scalar
 *            scatter-subtract, and the strict-nesting per-column batch
 *            kernel (kls_snode_batch_consume_cached analogue)
 *   snb    - supernodal dense-panel storage: amalgamated supernodes
 *            (greedy zeta-padding, wmax), column-major dense panels as
 *            primary storage, per-edge UB gather + unit-lower TRSM + GEMM
 *            into row-major temp + row scatter, fused no-pivot getrf/TRSM,
 *            hybrid scalar path for narrow producers (w < WNARROW).
 *
 * Usage: supernodal_bench matrix.mtx [passes] [zeta] [wnarrow] [wmax]
 *   env: SNB_SCALE (default 2), SNB_BTF (default 1), SNB_ORDERING (0=amd),
 *        SNB_TRACE=1 for waste counters and supernode histograms.
 *
 * Both non-klu kernels are verified against trilinos_klu_l_refactor to
 * roundoff before timing.
 */

#define _POSIX_C_SOURCE 200809L

#include "trilinos_klu_decl.h"

#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <time.h>

/* ------------------------------------------------------------------ */
/* small utilities                                                     */
/* ------------------------------------------------------------------ */

static double now_seconds(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1.0e-9 * (double)ts.tv_nsec;
}

static void *xmalloc(size_t bytes) {
  void *p = malloc(bytes ? bytes : 1u);
  if (p == NULL) {
    fprintf(stderr, "out of memory (%zu bytes)\n", bytes);
    exit(EXIT_FAILURE);
  }
  return p;
}

static void *xcalloc(size_t count, size_t size) {
  void *p = calloc(count ? count : 1u, size);
  if (p == NULL) {
    fprintf(stderr, "out of memory (%zu x %zu)\n", count, size);
    exit(EXIT_FAILURE);
  }
  return p;
}

/* ------------------------------------------------------------------ */
/* MatrixMarket reader (same conventions as klu_width_compare)         */
/* ------------------------------------------------------------------ */

typedef struct triplet {
  int64_t row;
  int64_t col;
  double value;
} triplet;

typedef struct matrix {
  int64_t n;
  int64_t nnz;
  int64_t *col_ptr;
  int64_t *row_idx;
  double *values;
} matrix;

static int cmp_triplet(const void *lhs, const void *rhs) {
  const triplet *a = (const triplet *)lhs;
  const triplet *b = (const triplet *)rhs;
  if (a->col != b->col) return (a->col < b->col) ? -1 : 1;
  if (a->row != b->row) return (a->row < b->row) ? -1 : 1;
  return 0;
}

static int starts_with_ci(const char *s, const char *prefix) {
  while (*prefix != '\0') {
    if (tolower((unsigned char)*s) != tolower((unsigned char)*prefix)) {
      return 0;
    }
    ++s;
    ++prefix;
  }
  return 1;
}

static int contains_token_ci(const char *s, const char *needle) {
  const size_t n = strlen(needle);
  for (; *s != '\0'; ++s) {
    if (strncasecmp(s, needle, n) == 0) return 1;
  }
  return 0;
}

static int read_matrix_market(const char *path, matrix *out) {
  FILE *fp = fopen(path, "r");
  if (fp == NULL) {
    fprintf(stderr, "cannot open %s: %s\n", path, strerror(errno));
    return 0;
  }
  char line[4096];
  if (fgets(line, sizeof(line), fp) == NULL ||
      !starts_with_ci(line, "%%MatrixMarket matrix coordinate")) {
    fprintf(stderr, "%s is not a coordinate MatrixMarket file\n", path);
    fclose(fp);
    return 0;
  }
  const int is_pattern = contains_token_ci(line, "pattern");
  const int is_complex = contains_token_ci(line, "complex");
  const int is_symmetric = contains_token_ci(line, "symmetric");
  const int is_skew = contains_token_ci(line, "skew-symmetric");
  if (is_complex) {
    fprintf(stderr, "%s is complex; real only\n", path);
    fclose(fp);
    return 0;
  }
  do {
    if (fgets(line, sizeof(line), fp) == NULL) {
      fclose(fp);
      return 0;
    }
  } while (line[0] == '%');
  int64_t rows = 0, cols = 0, entries = 0;
  if (sscanf(line, "%" SCNd64 " %" SCNd64 " %" SCNd64, &rows, &cols,
             &entries) != 3 ||
      rows <= 0 || rows != cols) {
    fprintf(stderr, "%s must be square\n", path);
    fclose(fp);
    return 0;
  }
  triplet *items = (triplet *)xmalloc(
    (size_t)(2 * entries + 1) * sizeof(*items));
  int64_t count = 0;
  for (int64_t k = 0; k < entries; ++k) {
    if (fgets(line, sizeof(line), fp) == NULL) {
      fclose(fp);
      free(items);
      return 0;
    }
    int64_t row = 0, col = 0;
    double value = 1.0;
    const int fields =
      is_pattern ? sscanf(line, "%" SCNd64 " %" SCNd64, &row, &col)
                 : sscanf(line, "%" SCNd64 " %" SCNd64 " %lf", &row, &col,
                          &value);
    if (fields < 2 || row <= 0 || row > rows || col <= 0 || col > cols) {
      fclose(fp);
      free(items);
      return 0;
    }
    --row;
    --col;
    items[count].row = row;
    items[count].col = col;
    items[count].value = value;
    ++count;
    if ((is_symmetric || is_skew) && row != col) {
      items[count].row = col;
      items[count].col = row;
      items[count].value = is_skew ? -value : value;
      ++count;
    }
  }
  fclose(fp);
  qsort(items, (size_t)count, sizeof(*items), cmp_triplet);
  int64_t unique = 0;
  for (int64_t i = 0; i < count;) {
    int64_t j = i + 1;
    double sum = items[i].value;
    while (j < count && items[j].row == items[i].row &&
           items[j].col == items[i].col) {
      sum += items[j].value;
      ++j;
    }
    items[unique] = items[i];
    items[unique].value = sum;
    ++unique;
    i = j;
  }
  out->n = rows;
  out->nnz = unique;
  out->col_ptr = (int64_t *)xcalloc((size_t)rows + 1u, sizeof(int64_t));
  out->row_idx = (int64_t *)xmalloc((size_t)unique * sizeof(int64_t));
  out->values = (double *)xmalloc((size_t)unique * sizeof(double));
  for (int64_t p = 0; p < unique; ++p) out->col_ptr[items[p].col + 1]++;
  for (int64_t c = 0; c < rows; ++c) out->col_ptr[c + 1] += out->col_ptr[c];
  int64_t *next = (int64_t *)xmalloc((size_t)rows * sizeof(int64_t));
  memcpy(next, out->col_ptr, (size_t)rows * sizeof(int64_t));
  for (int64_t p = 0; p < unique; ++p) {
    const int64_t dst = next[items[p].col]++;
    out->row_idx[dst] = items[p].row;
    out->values[dst] = items[p].value;
  }
  free(next);
  free(items);
  return 1;
}

/* ------------------------------------------------------------------ */
/* KLU packed-LU access (mirrors kls_klu_get_pointer)                  */
/* ------------------------------------------------------------------ */

static size_t units_for_indices(UF_long length) {
  const size_t bytes = (size_t)length * sizeof(UF_long);
  return (bytes + sizeof(double) - 1u) / sizeof(double);
}

static void lu_get(double *lu, const UF_long *offsets, const UF_long *lengths,
                   UF_long k, UF_long **indices_out, double **values_out,
                   UF_long *length_out) {
  const UF_long length = lengths[k];
  double *base = lu + offsets[k];
  *indices_out = (UF_long *)base;
  *values_out = base + units_for_indices(length);
  *length_out = length;
}

/* ------------------------------------------------------------------ */
/* shared: precomputed A-scatter program                               */
/*                                                                     */
/* Replays klu_l_refactor's input walk: for each block column (global  */
/* order), each A entry goes to the off-diagonal store or into the     */
/* working vector/panel.  src = A value index; scale by 1/Rs[oldrow]   */
/* when scaling is on (division kept to match KLU roundoff exactly).   */
/* ------------------------------------------------------------------ */

typedef struct scatter_entry {
  int64_t dst;  /* kernel-specific destination slot */
  int64_t src;  /* index into A values */
  int64_t oldrow;
} scatter_entry;

typedef struct scatter_program {
  scatter_entry *off;  /* off-diagonal entries: dst = Offx slot */
  int64_t off_count;
  scatter_entry *blk;  /* in-block entries: dst = kernel slot */
  int64_t blk_count;
  int64_t *col_off_ptr;  /* per global col: range in off[] */
  int64_t *col_blk_ptr;  /* per global col: range in blk[] */
} scatter_program;

static void scatter_program_free(scatter_program *sp) {
  free(sp->off);
  free(sp->blk);
  free(sp->col_off_ptr);
  free(sp->col_blk_ptr);
}

/* ------------------------------------------------------------------ */
/* landed-kernel replica                                               */
/* ------------------------------------------------------------------ */

#define SNODE_MIN_BATCH 3
#define SNODE_MAX_BATCH 64
#define SNODE_TAIL_CHUNK 32
#define SNODE_MIN_BATCH_WORK 192

typedef struct landed_state {
  matrix *a;
  trilinos_klu_l_symbolic *symbolic;
  trilinos_klu_l_numeric *numeric; /* pattern reference (not written) */
  double **lu_copy;                /* per block: private LUbx value copy */
  double *udiag;
  double *offx;
  double *rs;      /* recomputed per refactor when scale > 0 */
  double *x;       /* maxblock workspace */
  int32_t **li32;  /* per global col: 32-bit L row mirror */
  int32_t **ui32;  /* per global col: 32-bit U row mirror */
  UF_long *snode_run_end; /* strict-nesting runs, global cols */
  scatter_program sp;     /* blk dst = block-local row (needs k1 bias) */
  int scale;
  int status;
} landed_state;

static void landed_scatter_subtract_i32(double *restrict x,
                                        const int32_t *restrict rows,
                                        const double *restrict values,
                                        UF_long length, double scale) {
  if (scale == 0.0) return;
  UF_long p = 0;
  for (; p + 7 < length; p += 8) {
    if (p + 24 < length) {
      __builtin_prefetch(&x[rows[p + 16]], 1, 1);
      __builtin_prefetch(&x[rows[p + 24]], 1, 1);
    }
    x[rows[p]] -= values[p] * scale;
    x[rows[p + 1]] -= values[p + 1] * scale;
    x[rows[p + 2]] -= values[p + 2] * scale;
    x[rows[p + 3]] -= values[p + 3] * scale;
    x[rows[p + 4]] -= values[p + 4] * scale;
    x[rows[p + 5]] -= values[p + 5] * scale;
    x[rows[p + 6]] -= values[p + 6] * scale;
    x[rows[p + 7]] -= values[p + 7] * scale;
  }
  for (; p < length; ++p) {
    x[rows[p]] -= values[p] * scale;
  }
}

/* Replica of kls_snode_batch_consume_cached on the private LU copy.
 * Returns the number of producers consumed (0 = no batch here). */
static UF_long landed_batch_consume(const landed_state *st, UF_long k1,
                                    double *lu_vals, const UF_long *lip,
                                    const UF_long *llen, const int32_t *ui32,
                                    double *ux, UF_long ucol_len, UF_long up,
                                    double *restrict x) {
  const UF_long j = (UF_long)ui32[up];
  const UF_long run_end = st->snode_run_end[k1 + j];
  if (run_end <= k1 + j + 1) return 0;
  UF_long tmax = run_end - (k1 + j);
  if (tmax > SNODE_MAX_BATCH) tmax = SNODE_MAX_BATCH;
  if (tmax > ucol_len - up) tmax = ucol_len - up;
  UF_long t = 1;
  while (t < tmax && (UF_long)ui32[up + t] == j + t) t++;
  if (t < SNODE_MIN_BATCH) return 0;
  UF_long *tli;
  double *tlx;
  UF_long tlen;
  lu_get(lu_vals, lip, llen, j + t - 1, &tli, &tlx, &tlen);
  if (t * tlen < SNODE_MIN_BATCH_WORK) return 0;
  double xs[SNODE_MAX_BATCH];
  const double *lx_arr[SNODE_MAX_BATCH];
  for (UF_long i = 0; i < t; ++i) {
    UF_long *li_i;
    double *lx_i;
    UF_long len_i;
    lu_get(lu_vals, lip, llen, j + i, &li_i, &lx_i, &len_i);
    lx_arr[i] = lx_i;
    xs[i] = x[j + i];
  }
  for (UF_long i = 0; i < t; ++i) x[j + i] = 0.0;
  for (UF_long i = 0; i < t; ++i) {
    const double u = xs[i];
    ux[up + i] = u;
    const double *lxi = lx_arr[i];
    for (UF_long r0 = 0; r0 + i + 1 < t; ++r0) {
      xs[i + 1 + r0] -= lxi[r0] * u;
    }
  }
  const double *tlx_last = lx_arr[t - 1];
  for (UF_long p0 = 0; p0 < tlen; p0 += SNODE_TAIL_CHUNK) {
    const UF_long pc =
      tlen - p0 < SNODE_TAIL_CHUNK ? tlen - p0 : SNODE_TAIL_CHUNK;
    double acc[SNODE_TAIL_CHUNK];
    {
      const double *src = tlx_last + p0;
      const double u = xs[t - 1];
      for (UF_long p = 0; p < pc; ++p) acc[p] = src[p] * u;
    }
    for (UF_long i = 0; i + 1 < t; ++i) {
      const double *src = lx_arr[i] + (t - 1 - i) + p0;
      const double u = xs[i];
      for (UF_long p = 0; p < pc; ++p) acc[p] += src[p] * u;
    }
    for (UF_long p = 0; p < pc; ++p) x[tli[p0 + p]] -= acc[p];
  }
  return t;
}

static int landed_refactor(landed_state *st) {
  const trilinos_klu_l_symbolic *symbolic = st->symbolic;
  const trilinos_klu_l_numeric *numeric = st->numeric;
  const matrix *a = st->a;
  const UF_long nblocks = symbolic->nblocks;
  double *x = st->x;
  st->status = 0;

  if (st->scale > 0) {
    /* recompute row scales from A, matching TRILINOS_KLU_scale(2,...) */
    double *rs = st->rs;
    const int64_t n = a->n;
    for (int64_t i = 0; i < n; ++i) rs[i] = 0.0;
    for (int64_t c = 0; c < n; ++c) {
      for (int64_t p = a->col_ptr[c]; p < a->col_ptr[c + 1]; ++p) {
        const double v = fabs(a->values[p]);
        if (st->scale == 1) {
          rs[a->row_idx[p]] += v;
        } else if (v > rs[a->row_idx[p]]) {
          rs[a->row_idx[p]] = v;
        }
      }
    }
    for (int64_t i = 0; i < n; ++i) {
      if (rs[i] == 0.0) rs[i] = 1.0;
    }
  }
  const double *rs = st->rs;
  const double *az = a->values;
  const scatter_program *sp = &st->sp;

  for (UF_long block = 0; block < nblocks; ++block) {
    const UF_long k1 = symbolic->R[block];
    const UF_long k2 = symbolic->R[block + 1];
    const UF_long nk = k2 - k1;
    if (nk == 1) {
      double s = 0.0;
      for (int64_t p = sp->col_off_ptr[k1]; p < sp->col_off_ptr[k1 + 1];
           ++p) {
        st->offx[sp->off[p].dst] =
          st->scale > 0 ? az[sp->off[p].src] / rs[sp->off[p].oldrow]
                        : az[sp->off[p].src];
      }
      for (int64_t p = sp->col_blk_ptr[k1]; p < sp->col_blk_ptr[k1 + 1];
           ++p) {
        s = st->scale > 0 ? az[sp->blk[p].src] / rs[sp->blk[p].oldrow]
                          : az[sp->blk[p].src];
      }
      st->udiag[k1] = s;
      continue;
    }
    const UF_long *lip = numeric->Lip + k1;
    const UF_long *llen = numeric->Llen + k1;
    const UF_long *uip = numeric->Uip + k1;
    const UF_long *ulen = numeric->Ulen + k1;
    double *lu_vals = st->lu_copy[block];
    for (UF_long k = 0; k < nk; ++k) {
      const UF_long gk = k1 + k;
      for (int64_t p = sp->col_off_ptr[gk]; p < sp->col_off_ptr[gk + 1];
           ++p) {
        st->offx[sp->off[p].dst] =
          st->scale > 0 ? az[sp->off[p].src] / rs[sp->off[p].oldrow]
                        : az[sp->off[p].src];
      }
      for (int64_t p = sp->col_blk_ptr[gk]; p < sp->col_blk_ptr[gk + 1];
           ++p) {
        x[sp->blk[p].dst] =
          st->scale > 0 ? az[sp->blk[p].src] / rs[sp->blk[p].oldrow]
                        : az[sp->blk[p].src];
      }
      UF_long *ui;
      double *ux;
      UF_long ucol_len;
      lu_get(lu_vals, uip, ulen, k, &ui, &ux, &ucol_len);
      const int32_t *ui32 = st->ui32[gk];
      UF_long up = 0;
      while (up < ucol_len) {
        if (st->snode_run_end != NULL) {
          const UF_long consumed = landed_batch_consume(
            st, k1, lu_vals, lip, llen, ui32, ux, ucol_len, up, x);
          if (consumed != 0) {
            up += consumed;
            continue;
          }
        }
        const UF_long j = (UF_long)ui32[up];
        const double ujk = x[j];
        x[j] = 0.0;
        ux[up] = ujk;
        if (ujk != 0.0) {
          UF_long *li;
          double *lx;
          UF_long lcol_len;
          lu_get(lu_vals, lip, llen, j, &li, &lx, &lcol_len);
          landed_scatter_subtract_i32(x, st->li32[k1 + j], lx, lcol_len,
                                      ujk);
        }
        up++;
      }
      const double ukk = x[k];
      x[k] = 0.0;
      if (ukk == 0.0) st->status = 1; /* singular; keep going like KLU */
      st->udiag[gk] = ukk;
      UF_long *li;
      double *lx;
      UF_long lcol_len;
      lu_get(lu_vals, lip, llen, k, &li, &lx, &lcol_len);
      for (UF_long p = 0; p < lcol_len; ++p) {
        const UF_long i = li[p];
        lx[p] = x[i] / ukk;
        x[i] = 0.0;
      }
    }
  }
  return 1;
}

/* ------------------------------------------------------------------ */
/* supernodal dense-panel kernel                                       */
/* ------------------------------------------------------------------ */

typedef struct snb_edge {
  int32_t producer;   /* supernode id within the block */
  int32_t u_start;    /* target-panel local row of this edge's U rows */
  int32_t sub_count;  /* # producer columns touched by the target */
  int64_t sub_off;    /* offset into block subcols[] */
  int64_t dst_off;    /* offset into block edge_dsts[] (producer below
                         rows -> target panel rows, -1 = miss) */
} snb_edge;

typedef struct snb_snode {
  int64_t start;      /* first column (block-local) */
  int32_t width;
  int32_t below;      /* # below-diagonal union rows */
  int64_t panel;      /* offset into block panel storage */
  int32_t height;     /* uw + width + below (panel ld) */
  int32_t uw;         /* widened U rows */
  int64_t rows_off;   /* offset into block row-list storage (below rows) */
  int64_t edges_off;  /* offset into block edge storage */
  int32_t edge_count;
} snb_snode;

typedef struct snb_block {
  UF_long k1, nk;
  int64_t snode_count;
  snb_snode *snodes;
  int32_t *col2snode;          /* block-local col -> snode id */
  int32_t *below_rows;         /* concatenated below-row lists */
  snb_edge *edges;             /* concatenated edge lists */
  int32_t *subcols;            /* concatenated producer-local col lists */
  int32_t *edge_dsts;          /* precomputed per-edge scatter targets */
  double *panels;              /* all panels, column-major per snode */
  int64_t panel_doubles;
  int32_t *map;                /* block-local row -> slot (setup/verify) */
} snb_block;

typedef struct snb_state {
  matrix *a;
  trilinos_klu_l_symbolic *symbolic;
  trilinos_klu_l_numeric *numeric;
  snb_block *blocks;
  UF_long nblocks;
  double *udiag;
  double *offx;
  double *rs;
  scatter_program sp; /* blk dst = encoded (block, snode, panel slot) */
  int scale;
  int status;
  int wnarrow;
  /* scratch */
  double *gemm_temp;   /* max below x wmax */
  double *ub;          /* wmax x wmax, column-major, then row-major copy */
  /* waste counters (round 3 diagnostics) */
  int64_t scatter_hit;
  int64_t scatter_miss;
  int64_t gemm_flops;
  int64_t scalar_flops;
  /* SNB_PROFILE=1 rdtsc phase attribution */
  int profile;
  uint64_t phase_cycles[5]; /* scatter, narrow, trsm, gemm, getrf */
} snb_state;

static inline uint64_t snb_tsc(void) {
#if defined(__x86_64__)
  unsigned lo, hi;
  __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
  return ((uint64_t)hi << 32) | lo;
#else
  return 0;
#endif
}

static int cmp_uf(const void *x, const void *y) {
  const UF_long a = *(const UF_long *)x;
  const UF_long b = *(const UF_long *)y;
  return a < b ? -1 : (a > b ? 1 : 0);
}

static void snb_map_set(const snb_block *blk, const snb_snode *sn);
static void snb_map_clear(const snb_block *blk, const snb_snode *sn);

/* Greedy consecutive amalgamation with zeta padding budget.
 * Patterns must be sorted (klu_l_sort done).  Returns snode count. */
static int64_t snb_amalgamate(const UF_long *lip, const UF_long *llen,
                              double *lu, UF_long nk, double zeta, int wmax,
                              int32_t *col2snode, int64_t *starts_out,
                              UF_long *union_rows, UF_long *merge_tmp) {
  int64_t count = 0;
  UF_long c = 0;
  while (c < nk) {
    UF_long *li;
    double *lx;
    UF_long len;
    lu_get(lu, lip, llen, c, &li, &lx, &len);
    UF_long ulen_cur = len;
    memcpy(union_rows, li, (size_t)len * sizeof(UF_long));
    int64_t true_entries = len;
    UF_long w = 1;
    while (c + w < nk && w < (UF_long)wmax) {
      UF_long *li2;
      double *lx2;
      UF_long len2;
      lu_get(lu, lip, llen, c + w, &li2, &lx2, &len2);
      /* merge candidate pattern into union (rows beyond the new diag) */
      UF_long i = 0, j = 0, m = 0;
      while (i < ulen_cur && j < len2) {
        UF_long ri = union_rows[i];
        UF_long rj = li2[j];
        UF_long r = ri < rj ? ri : rj;
        if (ri == r) i++;
        if (rj == r) j++;
        merge_tmp[m++] = r;
      }
      while (i < ulen_cur) merge_tmp[m++] = union_rows[i++];
      while (j < len2) merge_tmp[m++] = li2[j++];
      /* dense size with the merge: below rows are union entries beyond
       * c+w (the new supernode end is c+w+1; diag block grows) */
      UF_long below = 0;
      for (UF_long t = 0; t < m; ++t) {
        if (merge_tmp[t] > c + w) below++;
      }
      const UF_long nw = w + 1;
      const double dense_entries =
        (double)nw * (double)below + (double)nw * (double)(nw - 1) / 2.0;
      const double truth = (double)(true_entries + len2);
      if (dense_entries > (1.0 + zeta) * truth) break;
      memcpy(union_rows, merge_tmp, (size_t)m * sizeof(UF_long));
      ulen_cur = m;
      true_entries += len2;
      w++;
    }
    starts_out[count] = c;
    for (UF_long i = 0; i < w; ++i) col2snode[c + i] = (int32_t)count;
    count++;
    c += w;
  }
  return count;
}

static void snb_free(snb_state *st) {
  if (st->blocks != NULL) {
    for (UF_long b = 0; b < st->nblocks; ++b) {
      snb_block *blk = &st->blocks[b];
      free(blk->snodes);
      free(blk->col2snode);
      free(blk->below_rows);
      free(blk->edges);
      free(blk->subcols);
      free(blk->edge_dsts);
      free(blk->panels);
      free(blk->map);
    }
    free(st->blocks);
  }
  free(st->gemm_temp);
  free(st->ub);
  scatter_program_free(&st->sp);
}

/* Build supernodal structures for every block.  Encoded panel slot for
 * the A-scatter program: (block b already implied by column) dst =
 * snode_id * 2^32 + local_panel_index is overkill; we store panel
 * global offset directly since panels live in one array per block. */
static void snb_setup(snb_state *st, double zeta, int wmax) {
  const trilinos_klu_l_symbolic *symbolic = st->symbolic;
  const trilinos_klu_l_numeric *numeric = st->numeric;
  const UF_long nblocks = symbolic->nblocks;
  st->nblocks = nblocks;
  st->blocks = (snb_block *)xcalloc((size_t)nblocks, sizeof(snb_block));
  UF_long maxblock = symbolic->maxblock;
  UF_long *union_rows =
    (UF_long *)xmalloc((size_t)maxblock * sizeof(UF_long));
  UF_long *merge_tmp =
    (UF_long *)xmalloc((size_t)maxblock * sizeof(UF_long));
  int64_t *starts = (int64_t *)xmalloc((size_t)maxblock * sizeof(int64_t));
  int64_t max_below = 0;

  for (UF_long block = 0; block < nblocks; ++block) {
    const UF_long k1 = symbolic->R[block];
    const UF_long k2 = symbolic->R[block + 1];
    const UF_long nk = k2 - k1;
    snb_block *blk = &st->blocks[block];
    blk->k1 = k1;
    blk->nk = nk;
    if (nk < 2) continue;
    const UF_long *lip = numeric->Lip + k1;
    const UF_long *llen = numeric->Llen + k1;
    const UF_long *uip = numeric->Uip + k1;
    const UF_long *ulen = numeric->Ulen + k1;
    double *lu = (double *)numeric->LUbx[block];

    blk->col2snode = (int32_t *)xmalloc((size_t)nk * sizeof(int32_t));
    const int64_t ns = snb_amalgamate(lip, llen, lu, nk, zeta, wmax,
                                      blk->col2snode, starts, union_rows,
                                      merge_tmp);
    blk->snode_count = ns;
    blk->snodes = (snb_snode *)xcalloc((size_t)ns, sizeof(snb_snode));
    for (int64_t s = 0; s < ns; ++s) {
      blk->snodes[s].start = starts[s];
      const int64_t end =
        s + 1 < ns ? starts[s + 1] : (int64_t)nk;
      blk->snodes[s].width = (int32_t)(end - starts[s]);
    }
    /* below-row unions */
    int64_t rows_total = 0;
    for (int64_t s = 0; s < ns; ++s) {
      snb_snode *sn = &blk->snodes[s];
      const int64_t end = sn->start + sn->width;
      UF_long ulen_cur = 0;
      for (int64_t c = sn->start; c < end; ++c) {
        UF_long *li;
        double *lx;
        UF_long len;
        lu_get(lu, lip, llen, c, &li, &lx, &len);
        UF_long i = 0, j = 0, m = 0;
        while (i < ulen_cur && j < len) {
          UF_long ri = union_rows[i];
          UF_long rj = li[j];
          UF_long r = ri < rj ? ri : rj;
          if (ri == r) i++;
          if (rj == r) j++;
          merge_tmp[m++] = r;
        }
        while (i < ulen_cur) merge_tmp[m++] = union_rows[i++];
        while (j < len) merge_tmp[m++] = li[j++];
        memcpy(union_rows, merge_tmp, (size_t)m * sizeof(UF_long));
        ulen_cur = m;
      }
      int64_t below = 0;
      for (UF_long t = 0; t < ulen_cur; ++t) {
        if (union_rows[t] >= (UF_long)end) below++;
      }
      sn->below = (int32_t)below;
      sn->rows_off = rows_total;
      rows_total += below;
      if (below > max_below) max_below = below;
    }
    blk->below_rows = (int32_t *)xmalloc((size_t)rows_total *
                                         sizeof(int32_t));
    for (int64_t s = 0; s < ns; ++s) {
      snb_snode *sn = &blk->snodes[s];
      const int64_t end = sn->start + sn->width;
      UF_long ulen_cur = 0;
      for (int64_t c = sn->start; c < end; ++c) {
        UF_long *li;
        double *lx;
        UF_long len;
        lu_get(lu, lip, llen, c, &li, &lx, &len);
        UF_long i = 0, j = 0, m = 0;
        while (i < ulen_cur && j < len) {
          UF_long ri = union_rows[i];
          UF_long rj = li[j];
          UF_long r = ri < rj ? ri : rj;
          if (ri == r) i++;
          if (rj == r) j++;
          merge_tmp[m++] = r;
        }
        while (i < ulen_cur) merge_tmp[m++] = union_rows[i++];
        while (j < len) merge_tmp[m++] = li[j++];
        memcpy(union_rows, merge_tmp, (size_t)m * sizeof(UF_long));
        ulen_cur = m;
      }
      int64_t below = 0;
      for (UF_long t = 0; t < ulen_cur; ++t) {
        if (union_rows[t] >= (UF_long)end) {
          blk->below_rows[sn->rows_off + below] = (int32_t)union_rows[t];
          below++;
        }
      }
    }
    /* edges: per target, the exact union of member columns' U rows,
     * grouped by producer supernode (rows sorted, so groups come out in
     * ascending producer order).  Each edge keeps the producer-local
     * column list; untouched producer columns carry exact zeros and are
     * excluded from the TRSM/GEMM. */
    int64_t edges_total = 0;
    int64_t subs_total = 0;
    unsigned char *inset = (unsigned char *)xcalloc((size_t)nk, 1u);
    UF_long *urows = union_rows; /* reuse block-sized scratch */
    for (int64_t s = 0; s < ns; ++s) {
      snb_snode *sn = &blk->snodes[s];
      const int64_t end = sn->start + sn->width;
      int64_t m = 0;
      for (int64_t c = sn->start; c < end; ++c) {
        UF_long *ui;
        double *ux;
        UF_long len;
        lu_get(lu, uip, ulen, c, &ui, &ux, &len);
        for (UF_long p = 0; p < len; ++p) {
          const UF_long j = ui[p];
          if (j >= (UF_long)sn->start) continue; /* inside the snode */
          if (!inset[j]) {
            inset[j] = 1;
            urows[m++] = j;
          }
        }
      }
      for (int64_t t = 0; t < m; ++t) inset[urows[t]] = 0;
      /* count edges = distinct producers among rows */
      if (m > 1) {
        /* sort rows ascending (they arrive per-column sorted; the
         * cross-column union is not).  qsort is setup-only cost. */
        int cmp_needed = 0;
        for (int64_t t = 1; t < m; ++t) {
          if (urows[t - 1] > urows[t]) {
            cmp_needed = 1;
            break;
          }
        }
        if (cmp_needed) {
          /* simple insertion sort is fine for short lists; fall back to
           * qsort for long ones */
          if (m <= 64) {
            for (int64_t i2 = 1; i2 < m; ++i2) {
              UF_long v = urows[i2];
              int64_t j2 = i2 - 1;
              while (j2 >= 0 && urows[j2] > v) {
                urows[j2 + 1] = urows[j2];
                j2--;
              }
              urows[j2 + 1] = v;
            }
          } else {
            qsort(urows, (size_t)m, sizeof(UF_long), cmp_uf);
          }
        }
      }
      int64_t cnt = 0;
      int32_t last_prod = -1;
      for (int64_t t = 0; t < m; ++t) {
        const int32_t prod = blk->col2snode[urows[t]];
        if (prod != last_prod) {
          cnt++;
          last_prod = prod;
        }
      }
      sn->edge_count = (int32_t)cnt;
      sn->edges_off = edges_total;
      edges_total += cnt;
      subs_total += m;
      sn->uw = (int32_t)m;
    }
    blk->edges = (snb_edge *)xmalloc((size_t)edges_total * sizeof(snb_edge));
    blk->subcols =
      (int32_t *)xmalloc((size_t)subs_total * sizeof(int32_t));
    int64_t panel_doubles = 0;
    int64_t sub_cursor = 0;
    for (int64_t s = 0; s < ns; ++s) {
      snb_snode *sn = &blk->snodes[s];
      const int64_t end = sn->start + sn->width;
      int64_t m = 0;
      for (int64_t c = sn->start; c < end; ++c) {
        UF_long *ui;
        double *ux;
        UF_long len;
        lu_get(lu, uip, ulen, c, &ui, &ux, &len);
        for (UF_long p = 0; p < len; ++p) {
          const UF_long j = ui[p];
          if (j >= (UF_long)sn->start) continue;
          if (!inset[j]) {
            inset[j] = 1;
            urows[m++] = j;
          }
        }
      }
      for (int64_t t = 0; t < m; ++t) inset[urows[t]] = 0;
      if (m > 1) {
        if (m <= 64) {
          for (int64_t i2 = 1; i2 < m; ++i2) {
            UF_long v = urows[i2];
            int64_t j2 = i2 - 1;
            while (j2 >= 0 && urows[j2] > v) {
              urows[j2 + 1] = urows[j2];
              j2--;
            }
            urows[j2 + 1] = v;
          }
        } else {
          qsort(urows, (size_t)m, sizeof(UF_long), cmp_uf);
        }
      }
      snb_edge *edges = blk->edges + sn->edges_off;
      int64_t cnt = -1;
      int32_t last_prod = -1;
      for (int64_t t = 0; t < m; ++t) {
        const int32_t prod = blk->col2snode[urows[t]];
        if (prod != last_prod) {
          cnt++;
          edges[cnt].producer = prod;
          edges[cnt].u_start = (int32_t)t;
          edges[cnt].sub_count = 0;
          edges[cnt].sub_off = sub_cursor + t;
          last_prod = prod;
        }
        edges[cnt].sub_count++;
        blk->subcols[sub_cursor + t] =
          (int32_t)(urows[t] - blk->snodes[prod].start);
      }
      sub_cursor += m;
      sn->height = sn->uw + sn->width + sn->below;
      sn->panel = panel_doubles;
      panel_doubles += (int64_t)sn->height * sn->width;
    }
    free(inset);
    blk->panel_doubles = panel_doubles;
    blk->panels = (double *)xmalloc((size_t)panel_doubles * sizeof(double));
    memset(blk->panels, 0, (size_t)panel_doubles * sizeof(double));
    blk->map = (int32_t *)xmalloc((size_t)nk * sizeof(int32_t));
    memset(blk->map, 0xff, (size_t)nk * sizeof(int32_t));
    /* precomputed per-edge scatter destinations (round 3: the refactor
     * map analogue) -- replaces map set/clear + random lookups in the
     * refactor hot loop with a sequential read */
    int64_t dsts_total = 0;
    for (int64_t s = 0; s < ns; ++s) {
      const snb_snode *sn = &blk->snodes[s];
      snb_edge *edges = blk->edges + sn->edges_off;
      for (int32_t e = 0; e < sn->edge_count; ++e) {
        edges[e].dst_off = dsts_total;
        dsts_total += blk->snodes[edges[e].producer].below;
      }
    }
    blk->edge_dsts =
      (int32_t *)xmalloc((size_t)(dsts_total > 0 ? dsts_total : 1) *
                         sizeof(int32_t));
    for (int64_t s = 0; s < ns; ++s) {
      const snb_snode *sn = &blk->snodes[s];
      snb_map_set(blk, sn);
      const snb_edge *edges = blk->edges + sn->edges_off;
      for (int32_t e = 0; e < sn->edge_count; ++e) {
        const snb_snode *ps = &blk->snodes[edges[e].producer];
        const int32_t *prows = blk->below_rows + ps->rows_off;
        int32_t *dsts = blk->edge_dsts + edges[e].dst_off;
        for (int32_t r = 0; r < ps->below; ++r) {
          dsts[r] = blk->map[prows[r]];
        }
      }
      snb_map_clear(blk, sn);
    }
  }
  free(union_rows);
  free(merge_tmp);
  free(starts);
  st->gemm_temp = (double *)xmalloc((size_t)(max_below > 0 ? max_below : 1) *
                                    (size_t)wmax * sizeof(double));
  st->ub = (double *)xmalloc(2u * (size_t)wmax * (size_t)wmax *
                             sizeof(double));
}

/* map maintenance: set the target's rows, returns nothing.  Slots:
 * [0, uw) producer columns of touched producers, [uw, uw+w) own diag
 * rows, [uw+w, H) below rows.  -1 means "not in this target" (scatter
 * skip). */
static void snb_map_set(const snb_block *blk, const snb_snode *sn) {
  int32_t *map = blk->map;
  const snb_edge *edges = blk->edges + sn->edges_off;
  for (int32_t e = 0; e < sn->edge_count; ++e) {
    const int64_t pstart = blk->snodes[edges[e].producer].start;
    const int32_t *sub = blk->subcols + edges[e].sub_off;
    const int32_t base = edges[e].u_start;
    for (int32_t i = 0; i < edges[e].sub_count; ++i) {
      map[pstart + sub[i]] = base + i;
    }
  }
  const int32_t uw = sn->uw;
  for (int32_t i = 0; i < sn->width; ++i) {
    map[sn->start + i] = uw + i;
  }
  const int32_t *rows = blk->below_rows + sn->rows_off;
  const int32_t base = uw + sn->width;
  for (int32_t i = 0; i < sn->below; ++i) {
    map[rows[i]] = base + i;
  }
}

static void snb_map_clear(const snb_block *blk, const snb_snode *sn) {
  int32_t *map = blk->map;
  const snb_edge *edges = blk->edges + sn->edges_off;
  for (int32_t e = 0; e < sn->edge_count; ++e) {
    const int64_t pstart = blk->snodes[edges[e].producer].start;
    const int32_t *sub = blk->subcols + edges[e].sub_off;
    for (int32_t i = 0; i < edges[e].sub_count; ++i) {
      map[pstart + sub[i]] = -1;
    }
  }
  for (int32_t i = 0; i < sn->width; ++i) {
    map[sn->start + i] = -1;
  }
  const int32_t *rows = blk->below_rows + sn->rows_off;
  for (int32_t i = 0; i < sn->below; ++i) {
    map[rows[i]] = -1;
  }
}

static int snb_refactor(snb_state *st, int trace) {
  const matrix *a = st->a;
  st->status = 0;
  st->scatter_hit = 0;
  st->scatter_miss = 0;

  if (st->scale > 0) {
    double *rs = st->rs;
    const int64_t n = a->n;
    for (int64_t i = 0; i < n; ++i) rs[i] = 0.0;
    for (int64_t c = 0; c < n; ++c) {
      for (int64_t p = a->col_ptr[c]; p < a->col_ptr[c + 1]; ++p) {
        const double v = fabs(a->values[p]);
        if (st->scale == 1) {
          rs[a->row_idx[p]] += v;
        } else if (v > rs[a->row_idx[p]]) {
          rs[a->row_idx[p]] = v;
        }
      }
    }
    for (int64_t i = 0; i < n; ++i) {
      if (rs[i] == 0.0) rs[i] = 1.0;
    }
  }
  const double *rs = st->rs;
  const double *az = a->values;
  const scatter_program *sp = &st->sp;
  const int wnarrow = st->wnarrow;

  for (UF_long block = 0; block < st->nblocks; ++block) {
    snb_block *blk = &st->blocks[block];
    const UF_long k1 = blk->k1;
    const UF_long nk = blk->nk;
    if (nk == 1) {
      double s = 0.0;
      for (int64_t p = sp->col_off_ptr[k1]; p < sp->col_off_ptr[k1 + 1];
           ++p) {
        st->offx[sp->off[p].dst] =
          st->scale > 0 ? az[sp->off[p].src] / rs[sp->off[p].oldrow]
                        : az[sp->off[p].src];
      }
      for (int64_t p = sp->col_blk_ptr[k1]; p < sp->col_blk_ptr[k1 + 1];
           ++p) {
        s = st->scale > 0 ? az[sp->blk[p].src] / rs[sp->blk[p].oldrow]
                          : az[sp->blk[p].src];
      }
      st->udiag[k1] = s;
      continue;
    }
    double *panels = blk->panels;
    const int prof = st->profile;
    for (int64_t s = 0; s < blk->snode_count; ++s) {
      const snb_snode *sn = &blk->snodes[s];
      const int32_t w = sn->width;
      const int32_t H = sn->height;
      const int32_t uw = sn->uw;
      double *W = panels + sn->panel;
      uint64_t t0 = prof ? snb_tsc() : 0;
      memset(W, 0, (size_t)H * (size_t)w * sizeof(double));
      /* scatter A columns (program dst = local panel index incl. col) */
      for (int32_t c = 0; c < w; ++c) {
        const UF_long gk = k1 + sn->start + c;
        for (int64_t p = sp->col_off_ptr[gk]; p < sp->col_off_ptr[gk + 1];
             ++p) {
          st->offx[sp->off[p].dst] =
            st->scale > 0 ? az[sp->off[p].src] / rs[sp->off[p].oldrow]
                          : az[sp->off[p].src];
        }
        for (int64_t p = sp->col_blk_ptr[gk]; p < sp->col_blk_ptr[gk + 1];
             ++p) {
          W[sp->blk[p].dst] =
            st->scale > 0 ? az[sp->blk[p].src] / rs[sp->blk[p].oldrow]
                          : az[sp->blk[p].src];
        }
      }
      if (prof) {
        const uint64_t t1 = snb_tsc();
        st->phase_cycles[0] += t1 - t0;
        t0 = t1;
      }
      /* edges, ascending producer */
      const snb_edge *edges = blk->edges + sn->edges_off;
      for (int32_t e = 0; e < sn->edge_count; ++e) {
        const snb_snode *ps = &blk->snodes[edges[e].producer];
        const int32_t pw = ps->width;
        const int32_t ph = ps->below;
        const int32_t sc = edges[e].sub_count;
        const int32_t *sub = blk->subcols + edges[e].sub_off;
        const int32_t *dsts = blk->edge_dsts + edges[e].dst_off;
        const int32_t ustart = edges[e].u_start;
        const double *Lp = panels + ps->panel;   /* producer panel */
        const int32_t pld = ps->height;
        const int32_t puw = ps->uw;
        const double *Ldiag = Lp + puw;          /* pw x pw diag block */
        const double *Lbelow = Lp + puw + pw;    /* ph x pw below block */
        if (sc < wnarrow) {
          /* hybrid scalar path: register u, immediate scatter, u==0
           * skip; precomputed dsts replace the map (round 3) */
          for (int32_t c = 0; c < w; ++c) {
            double *Wc = W + (size_t)c * H;
            for (int32_t i = 0; i < sc; ++i) {
              const double u = Wc[ustart + i];
              if (u == 0.0) continue;
              const double *lcol = Ldiag + (size_t)sub[i] * pld;
              for (int32_t i2 = i + 1; i2 < sc; ++i2) {
                Wc[ustart + i2] -= lcol[sub[i2]] * u;
              }
              const double *lb = Lbelow + (size_t)sub[i] * pld;
              for (int32_t r = 0; r < ph; ++r) {
                const int32_t dst = dsts[r];
                if (dst >= 0) {
                  Wc[dst] -= lb[r] * u;
                }
              }
            }
          }
          if (prof) {
            const uint64_t t1 = snb_tsc();
            st->phase_cycles[1] += t1 - t0;
            t0 = t1;
          }
          continue;
        }
        /* gather UB (sc x w), TRSM restricted to the touched producer
         * columns (untouched ones carry exact zeros), writeback */
        double *B = st->ub;                    /* column-major sc x w */
        for (int32_t c = 0; c < w; ++c) {
          const double *Wc = W + (size_t)c * H + ustart;
          double *Bc = B + (size_t)c * sc;
          for (int32_t i = 0; i < sc; ++i) Bc[i] = Wc[i];
        }
        for (int32_t i = 0; i < sc; ++i) {
          const double *lcol = Ldiag + (size_t)sub[i] * pld;
          for (int32_t c = 0; c < w; ++c) {
            double *Bc = B + (size_t)c * sc;
            const double u = Bc[i];
            if (u == 0.0) continue;
            for (int32_t i2 = i + 1; i2 < sc; ++i2) {
              Bc[i2] -= lcol[sub[i2]] * u;
            }
          }
        }
        for (int32_t c = 0; c < w; ++c) {
          double *Wc = W + (size_t)c * H + ustart;
          const double *Bc = B + (size_t)c * sc;
          for (int32_t i = 0; i < sc; ++i) Wc[i] = Bc[i];
        }
        if (prof) {
          const uint64_t t1 = snb_tsc();
          st->phase_cycles[2] += t1 - t0;
          t0 = t1;
        }
        if (ph == 0) continue;
        /* block path: GEMM into row-major temp (long contiguous streams,
         * independent FMAs -- measured faster than a fused accumulate),
         * then row scatter through the precomputed dsts */
        double *Brow = st->ub + (size_t)sc * w; /* row-major mirror */
        for (int32_t c = 0; c < w; ++c) {
          const double *Bc = B + (size_t)c * sc;
          for (int32_t i = 0; i < sc; ++i) {
            Brow[(size_t)i * w + c] = Bc[i];
          }
        }
        double *T = st->gemm_temp;
        {
          const double *lb = Lbelow + (size_t)sub[0] * pld;
          const double *br = Brow;
          for (int32_t r = 0; r < ph; ++r) {
            const double l = lb[r];
            double *tr = T + (size_t)r * w;
            for (int32_t c = 0; c < w; ++c) tr[c] = l * br[c];
          }
          for (int32_t i = 1; i < sc; ++i) {
            const double *lbi = Lbelow + (size_t)sub[i] * pld;
            const double *bri = Brow + (size_t)i * w;
            for (int32_t r = 0; r < ph; ++r) {
              const double l = lbi[r];
              double *tr = T + (size_t)r * w;
              for (int32_t c = 0; c < w; ++c) tr[c] += l * bri[c];
            }
          }
        }
        for (int32_t r = 0; r < ph; ++r) {
          const int32_t dst = dsts[r];
          const double *tr = T + (size_t)r * w;
          if (dst >= 0) {
            double *Wd = W + dst;
            for (int32_t c = 0; c < w; ++c) {
              Wd[(size_t)c * H] -= tr[c];
            }
            if (trace) st->scatter_hit++;
          } else if (trace) {
            st->scatter_miss++;
          }
        }
        if (prof) {
          const uint64_t t1 = snb_tsc();
          st->phase_cycles[3] += t1 - t0;
          t0 = t1;
        }
      }
      /* fused no-pivot getrf/TRSM on rows [uw, H) */
      {
        double *diag = W + uw;
        for (int32_t i = 0; i < w; ++i) {
          double *ci = diag + (size_t)i * H;
          const double pivot = ci[i];
          st->udiag[k1 + sn->start + i] = pivot;
          if (pivot == 0.0) st->status = 1;
          const int32_t hrem = H - uw - i - 1;
          double *sub = ci + i + 1;
          /* division (not reciprocal) to match KLU's roundoff */
          for (int32_t r = 0; r < hrem; ++r) sub[r] /= pivot;
          for (int32_t j = i + 1; j < w; ++j) {
            double *cj = diag + (size_t)j * H;
            const double u = cj[i];
            if (u == 0.0) continue;
            double *dstv = cj + i + 1;
            for (int32_t r = 0; r < hrem; ++r) {
              dstv[r] -= sub[r] * u;
            }
          }
        }
      }
      if (prof) {
        const uint64_t t1 = snb_tsc();
        st->phase_cycles[4] += t1 - t0;
      }
    }
  }
  return 1;
}

/* ------------------------------------------------------------------ */
/* setup shared pieces                                                 */
/* ------------------------------------------------------------------ */

/* Build the A-scatter program.  dst semantics per kernel:
 *  landed: block-local row (X index)
 *  snb:    local panel index (row + col*H) of the OWNING target panel
 * Both kernels share off-diagonal entries (Offx slots in KLU order). */
static void build_scatter_programs(const matrix *a,
                                   const trilinos_klu_l_symbolic *symbolic,
                                   const trilinos_klu_l_numeric *numeric,
                                   const snb_state *snb,
                                   scatter_program *landed_sp,
                                   scatter_program *snb_sp) {
  const UF_long n = symbolic->n;
  const UF_long nblocks = symbolic->nblocks;
  const UF_long *pinv = numeric->Pinv;
  const UF_long nzoff = symbolic->nzoff;

  int64_t off_count = 0;
  int64_t blk_count = 0;
  for (UF_long block = 0; block < nblocks; ++block) {
    const UF_long k1 = symbolic->R[block];
    const UF_long k2 = symbolic->R[block + 1];
    for (UF_long k = k1; k < k2; ++k) {
      const UF_long oldcol = symbolic->Q[k];
      for (int64_t p = a->col_ptr[oldcol]; p < a->col_ptr[oldcol + 1];
           ++p) {
        const UF_long newrow = pinv[a->row_idx[p]] - k1;
        if (newrow < 0 && off_count < nzoff) {
          off_count++;
        } else {
          blk_count++;
        }
      }
    }
  }
  landed_sp->off =
    (scatter_entry *)xmalloc((size_t)off_count * sizeof(scatter_entry));
  landed_sp->blk =
    (scatter_entry *)xmalloc((size_t)blk_count * sizeof(scatter_entry));
  landed_sp->col_off_ptr = (int64_t *)xcalloc((size_t)n + 1, sizeof(int64_t));
  landed_sp->col_blk_ptr = (int64_t *)xcalloc((size_t)n + 1, sizeof(int64_t));
  landed_sp->off_count = off_count;
  landed_sp->blk_count = blk_count;
  snb_sp->off =
    (scatter_entry *)xmalloc((size_t)off_count * sizeof(scatter_entry));
  snb_sp->blk =
    (scatter_entry *)xmalloc((size_t)blk_count * sizeof(scatter_entry));
  snb_sp->col_off_ptr = (int64_t *)xcalloc((size_t)n + 1, sizeof(int64_t));
  snb_sp->col_blk_ptr = (int64_t *)xcalloc((size_t)n + 1, sizeof(int64_t));
  snb_sp->off_count = off_count;
  snb_sp->blk_count = blk_count;

  int64_t poff = 0;
  int64_t pblk = 0;
  for (UF_long block = 0; block < nblocks; ++block) {
    const UF_long k1 = symbolic->R[block];
    const UF_long k2 = symbolic->R[block + 1];
    const snb_block *blk = &snb->blocks[block];
    for (UF_long k = k1; k < k2; ++k) {
      landed_sp->col_off_ptr[k] = poff;
      landed_sp->col_blk_ptr[k] = pblk;
      const UF_long oldcol = symbolic->Q[k];
      for (int64_t p = a->col_ptr[oldcol]; p < a->col_ptr[oldcol + 1];
           ++p) {
        const UF_long oldrow = a->row_idx[p];
        const UF_long newrow = pinv[oldrow] - k1;
        if (newrow < 0 && poff < nzoff) {
          landed_sp->off[poff].dst = poff;
          landed_sp->off[poff].src = p;
          landed_sp->off[poff].oldrow = oldrow;
          snb_sp->off[poff] = landed_sp->off[poff];
          poff++;
        } else {
          landed_sp->blk[pblk].dst = newrow;
          landed_sp->blk[pblk].src = p;
          landed_sp->blk[pblk].oldrow = oldrow;
          /* snb: panel slot of the owning target */
          int64_t dst = -1;
          if (blk->nk >= 2) {
            const int32_t s = blk->col2snode[k - k1];
            const snb_snode *sn = &blk->snodes[s];
            const int32_t c = (int32_t)(k - k1 - sn->start);
            /* local row within the panel */
            int32_t lr = -1;
            if (newrow >= sn->start &&
                newrow < sn->start + sn->width) {
              lr = sn->uw + (int32_t)(newrow - sn->start);
            } else if (newrow >= sn->start + sn->width) {
              /* below rows: binary search the row list */
              const int32_t *rows = blk->below_rows + sn->rows_off;
              int32_t lo = 0, hi = sn->below - 1;
              while (lo <= hi) {
                const int32_t mid = (lo + hi) >> 1;
                if (rows[mid] < (int32_t)newrow) {
                  lo = mid + 1;
                } else if (rows[mid] > (int32_t)newrow) {
                  hi = mid - 1;
                } else {
                  lr = sn->uw + sn->width + mid;
                  break;
                }
              }
            } else {
              /* U region: find the producer edge, then the row within
               * its touched-column list */
              const int32_t prod = blk->col2snode[newrow];
              const snb_edge *edges = blk->edges + sn->edges_off;
              for (int32_t e2 = 0; e2 < sn->edge_count; ++e2) {
                if (edges[e2].producer == prod) {
                  const int32_t *sub = blk->subcols + edges[e2].sub_off;
                  const int32_t want =
                    (int32_t)(newrow - blk->snodes[prod].start);
                  int32_t lo = 0, hi = edges[e2].sub_count - 1;
                  while (lo <= hi) {
                    const int32_t mid = (lo + hi) >> 1;
                    if (sub[mid] < want) {
                      lo = mid + 1;
                    } else if (sub[mid] > want) {
                      hi = mid - 1;
                    } else {
                      lr = edges[e2].u_start + mid;
                      break;
                    }
                  }
                  break;
                }
              }
            }
            if (lr < 0) {
              fprintf(stderr,
                      "internal error: A entry outside panel structure "
                      "(block %ld col %ld row %ld)\n",
                      (long)block, (long)k, (long)newrow);
              exit(EXIT_FAILURE);
            }
            dst = (int64_t)lr + (int64_t)c * sn->height;
          } else {
            dst = 0; /* singleton: value only */
          }
          snb_sp->blk[pblk].dst = dst;
          snb_sp->blk[pblk].src = p;
          snb_sp->blk[pblk].oldrow = oldrow;
          pblk++;
        }
      }
    }
  }
  landed_sp->col_off_ptr[n] = poff;
  landed_sp->col_blk_ptr[n] = pblk;
  snb_sp->col_off_ptr[n] = poff;
  snb_sp->col_blk_ptr[n] = pblk;
  memcpy(snb_sp->col_off_ptr, landed_sp->col_off_ptr,
         ((size_t)n + 1) * sizeof(int64_t));
  memcpy(snb_sp->col_blk_ptr, landed_sp->col_blk_ptr,
         ((size_t)n + 1) * sizeof(int64_t));
}

/* ------------------------------------------------------------------ */
/* verification                                                        */
/* ------------------------------------------------------------------ */

/* Componentwise reconstruction residual: per block column,
 *   r = scatter(A/Rs) - sum_{j in Upat(k)} L(:,j)*U(j,k) - L(:,k)*ukk
 * with unit L diagonals expanded, |every term| accumulated alongside.
 * max |r|/acc is O(eps * column ops) for ANY correct kernel regardless
 * of summation order, so it verifies factors on ill-conditioned inputs
 * where direct value comparison is chaotic (mc2depi, rajat25). */
static double reconstruction_error(const matrix *a,
                                   const trilinos_klu_l_symbolic *symbolic,
                                   const trilinos_klu_l_numeric *numeric,
                                   double *const *lu_by_block,
                                   const double *udiag, const double *rs,
                                   int scale, const scatter_program *sp) {
  const double *az = a->values;
  double *x = (double *)xcalloc((size_t)symbolic->maxblock, sizeof(double));
  double *xa = (double *)xcalloc((size_t)symbolic->maxblock, sizeof(double));
  double worst = 0.0;
  for (UF_long block = 0; block < symbolic->nblocks; ++block) {
    const UF_long k1 = symbolic->R[block];
    const UF_long k2 = symbolic->R[block + 1];
    const UF_long nk = k2 - k1;
    if (nk < 2) {
      double s = 0.0;
      double sa = 0.0;
      for (int64_t p = sp->col_blk_ptr[k1]; p < sp->col_blk_ptr[k1 + 1];
           ++p) {
        const double v = scale > 0 ? az[sp->blk[p].src] / rs[sp->blk[p].oldrow]
                                   : az[sp->blk[p].src];
        s = v;
        sa = fabs(v);
      }
      const double err = fabs(s - udiag[k1]);
      const double den = sa > fabs(udiag[k1]) ? sa : fabs(udiag[k1]);
      if (den > 0.0 && err / den > worst) worst = err / den;
      continue;
    }
    const UF_long *lip = numeric->Lip + k1;
    const UF_long *llen = numeric->Llen + k1;
    const UF_long *uip = numeric->Uip + k1;
    const UF_long *ulen = numeric->Ulen + k1;
    double *lu = lu_by_block[block];
    for (UF_long k = 0; k < nk; ++k) {
      const UF_long gk = k1 + k;
      for (int64_t p = sp->col_blk_ptr[gk]; p < sp->col_blk_ptr[gk + 1];
           ++p) {
        const double v = scale > 0 ? az[sp->blk[p].src] / rs[sp->blk[p].oldrow]
                                   : az[sp->blk[p].src];
        x[sp->blk[p].dst] += v;
        xa[sp->blk[p].dst] += fabs(v);
      }
      UF_long *ui;
      double *ux;
      UF_long uclen;
      lu_get(lu, uip, ulen, k, &ui, &ux, &uclen);
      for (UF_long up = 0; up < uclen; ++up) {
        const UF_long j = ui[up];
        const double u = ux[up];
        x[j] -= u;
        xa[j] += fabs(u);
        UF_long *li;
        double *lx;
        UF_long lclen;
        lu_get(lu, lip, llen, j, &li, &lx, &lclen);
        for (UF_long p = 0; p < lclen; ++p) {
          const double t = lx[p] * u;
          x[li[p]] -= t;
          xa[li[p]] += fabs(t);
        }
      }
      const double ukk = udiag[gk];
      x[k] -= ukk;
      xa[k] += fabs(ukk);
      UF_long *li;
      double *lx;
      UF_long lclen;
      lu_get(lu, lip, llen, k, &li, &lx, &lclen);
      for (UF_long p = 0; p < lclen; ++p) {
        const double t = lx[p] * ukk;
        x[li[p]] -= t;
        xa[li[p]] += fabs(t);
      }
      /* score and reset the touched rows */
      for (int64_t p = sp->col_blk_ptr[gk]; p < sp->col_blk_ptr[gk + 1];
           ++p) {
        const int64_t r = sp->blk[p].dst;
        if (xa[r] > 0.0 && fabs(x[r]) / xa[r] > worst) {
          worst = fabs(x[r]) / xa[r];
        }
        x[r] = 0.0;
        xa[r] = 0.0;
      }
      for (UF_long up = 0; up < uclen; ++up) {
        const UF_long j = ui[up];
        if (xa[j] > 0.0 && fabs(x[j]) / xa[j] > worst) {
          worst = fabs(x[j]) / xa[j];
        }
        x[j] = 0.0;
        xa[j] = 0.0;
        UF_long *li2;
        double *lx2;
        UF_long lclen2;
        lu_get(lu, lip, llen, j, &li2, &lx2, &lclen2);
        for (UF_long p = 0; p < lclen2; ++p) {
          const UF_long r = li2[p];
          if (xa[r] > 0.0 && fabs(x[r]) / xa[r] > worst) {
            worst = fabs(x[r]) / xa[r];
          }
          x[r] = 0.0;
          xa[r] = 0.0;
        }
      }
      if (xa[k] > 0.0 && fabs(x[k]) / xa[k] > worst) {
        worst = fabs(x[k]) / xa[k];
      }
      x[k] = 0.0;
      xa[k] = 0.0;
      for (UF_long p = 0; p < lclen; ++p) {
        const UF_long r = li[p];
        if (xa[r] > 0.0 && fabs(x[r]) / xa[r] > worst) {
          worst = fabs(x[r]) / xa[r];
        }
        x[r] = 0.0;
        xa[r] = 0.0;
      }
    }
  }
  free(x);
  free(xa);
  return worst;
}

/* export snb panel values into a klu-layout value copy so the
 * reconstruction checker can walk them uniformly */
static void snb_export(const snb_state *st, double *const *lu_out) {
  const trilinos_klu_l_symbolic *symbolic = st->symbolic;
  const trilinos_klu_l_numeric *numeric = st->numeric;
  for (UF_long block = 0; block < symbolic->nblocks; ++block) {
    const snb_block *blk = &st->blocks[block];
    const UF_long k1 = blk->k1;
    const UF_long nk = blk->nk;
    if (nk < 2) continue;
    const UF_long *lip = numeric->Lip + k1;
    const UF_long *llen = numeric->Llen + k1;
    const UF_long *uip = numeric->Uip + k1;
    const UF_long *ulen = numeric->Ulen + k1;
    double *lu = lu_out[block];
    for (UF_long k = 0; k < nk; ++k) {
      const int32_t s = blk->col2snode[k];
      const snb_snode *sn = &blk->snodes[s];
      const int32_t c = (int32_t)(k - sn->start);
      const double *W = blk->panels + sn->panel + (size_t)c * sn->height;
      UF_long *ri;
      double *rv;
      UF_long rl;
      lu_get(lu, lip, llen, k, &ri, &rv, &rl);
      for (UF_long p = 0; p < rl; ++p) {
        const UF_long row = ri[p];
        if (row < (UF_long)(sn->start + sn->width)) {
          rv[p] = W[sn->uw + (row - sn->start)];
        } else {
          const int32_t *rows = blk->below_rows + sn->rows_off;
          int32_t lo = 0, hi = sn->below - 1;
          while (lo <= hi) {
            const int32_t mid = (lo + hi) >> 1;
            if (rows[mid] < (int32_t)row) {
              lo = mid + 1;
            } else if (rows[mid] > (int32_t)row) {
              hi = mid - 1;
            } else {
              rv[p] = W[sn->uw + sn->width + mid];
              break;
            }
          }
        }
      }
      lu_get(lu, uip, ulen, k, &ri, &rv, &rl);
      for (UF_long p = 0; p < rl; ++p) {
        const UF_long row = ri[p];
        if (row >= (UF_long)sn->start) {
          rv[p] = W[sn->uw + (row - sn->start)];
        } else {
          const int32_t prod = blk->col2snode[row];
          const snb_edge *edges = blk->edges + sn->edges_off;
          for (int32_t e2 = 0; e2 < sn->edge_count; ++e2) {
            if (edges[e2].producer == prod) {
              const int32_t *sub = blk->subcols + edges[e2].sub_off;
              const int32_t want =
                (int32_t)(row - blk->snodes[prod].start);
              int32_t lo = 0, hi = edges[e2].sub_count - 1;
              while (lo <= hi) {
                const int32_t mid = (lo + hi) >> 1;
                if (sub[mid] < want) {
                  lo = mid + 1;
                } else if (sub[mid] > want) {
                  hi = mid - 1;
                } else {
                  rv[p] = W[edges[e2].u_start + mid];
                  break;
                }
              }
              break;
            }
          }
        }
      }
    }
  }
}

typedef struct verify_result {
  double max_rel;      /* max per-column normwise relative difference */
  double udiag_max_rel;
  double off_max_rel;
} verify_result;

static double rel_diff(double a, double b, double scale) {
  const double d = fabs(a - b);
  const double m = fabs(b) > scale ? fabs(b) : scale;
  return m > 0.0 ? d / m : d;
}

/* compare landed L/U values (in lu_copy) against the reference numeric */
static void verify_landed(const landed_state *st, verify_result *out) {
  const trilinos_klu_l_symbolic *symbolic = st->symbolic;
  const trilinos_klu_l_numeric *numeric = st->numeric;
  out->max_rel = 0.0;
  out->udiag_max_rel = 0.0;
  out->off_max_rel = 0.0;
  for (UF_long block = 0; block < symbolic->nblocks; ++block) {
    const UF_long k1 = symbolic->R[block];
    const UF_long k2 = symbolic->R[block + 1];
    const UF_long nk = k2 - k1;
    if (nk < 2) continue;
    const UF_long *lip = numeric->Lip + k1;
    const UF_long *llen = numeric->Llen + k1;
    const UF_long *uip = numeric->Uip + k1;
    const UF_long *ulen = numeric->Ulen + k1;
    double *lu_ref = (double *)numeric->LUbx[block];
    double *lu_mine = st->lu_copy[block];
    for (UF_long k = 0; k < nk; ++k) {
      UF_long *ri, *mi;
      double *rv, *mv;
      UF_long rl, ml;
      lu_get(lu_ref, lip, llen, k, &ri, &rv, &rl);
      lu_get(lu_mine, lip, llen, k, &mi, &mv, &ml);
      double cn = 0.0;
      for (UF_long p = 0; p < rl; ++p) {
        if (fabs(rv[p]) > cn) cn = fabs(rv[p]);
      }
      for (UF_long p = 0; p < rl; ++p) {
        const double r = rel_diff(mv[p], rv[p], cn);
        if (r > out->max_rel) out->max_rel = r;
      }
      lu_get(lu_ref, uip, ulen, k, &ri, &rv, &rl);
      lu_get(lu_mine, uip, ulen, k, &mi, &mv, &ml);
      cn = 0.0;
      for (UF_long p = 0; p < rl; ++p) {
        if (fabs(rv[p]) > cn) cn = fabs(rv[p]);
      }
      for (UF_long p = 0; p < rl; ++p) {
        const double r = rel_diff(mv[p], rv[p], cn);
        if (r > out->max_rel) out->max_rel = r;
      }
    }
  }
  const double *udiag_ref = (const double *)numeric->Udiag;
  for (UF_long k = 0; k < symbolic->n; ++k) {
    const double r = rel_diff(st->udiag[k], udiag_ref[k], 0.0);
    if (r > out->udiag_max_rel) out->udiag_max_rel = r;
  }
  const double *off_ref = (const double *)numeric->Offx;
  for (UF_long p = 0; p < numeric->nzoff; ++p) {
    const double r = rel_diff(st->offx[p], off_ref[p], 0.0);
    if (r > out->off_max_rel) out->off_max_rel = r;
  }
}

/* compare supernodal panels against the reference numeric */
static void verify_snb(const snb_state *st, verify_result *out) {
  const trilinos_klu_l_symbolic *symbolic = st->symbolic;
  const trilinos_klu_l_numeric *numeric = st->numeric;
  out->max_rel = 0.0;
  out->udiag_max_rel = 0.0;
  out->off_max_rel = 0.0;
  for (UF_long block = 0; block < symbolic->nblocks; ++block) {
    const snb_block *blk = &st->blocks[block];
    const UF_long k1 = blk->k1;
    const UF_long nk = blk->nk;
    if (nk < 2) continue;
    const UF_long *lip = numeric->Lip + k1;
    const UF_long *llen = numeric->Llen + k1;
    const UF_long *uip = numeric->Uip + k1;
    const UF_long *ulen = numeric->Ulen + k1;
    double *lu_ref = (double *)numeric->LUbx[block];
    for (UF_long k = 0; k < nk; ++k) {
      const int32_t s = blk->col2snode[k];
      const snb_snode *sn = &blk->snodes[s];
      const int32_t c = (int32_t)(k - sn->start);
      const double *W = blk->panels + sn->panel + (size_t)c * sn->height;
      UF_long *ri;
      double *rv;
      UF_long rl;
      /* L column */
      lu_get(lu_ref, lip, llen, k, &ri, &rv, &rl);
      double cn = 0.0;
      for (UF_long p = 0; p < rl; ++p) {
        if (fabs(rv[p]) > cn) cn = fabs(rv[p]);
      }
      for (UF_long p = 0; p < rl; ++p) {
        const UF_long row = ri[p];
        double mine;
        if (row < (UF_long)(sn->start + sn->width)) {
          mine = W[sn->uw + (row - sn->start)];
        } else {
          const int32_t *rows = blk->below_rows + sn->rows_off;
          int32_t lo = 0, hi = sn->below - 1, at = -1;
          while (lo <= hi) {
            const int32_t mid = (lo + hi) >> 1;
            if (rows[mid] < (int32_t)row) {
              lo = mid + 1;
            } else if (rows[mid] > (int32_t)row) {
              hi = mid - 1;
            } else {
              at = mid;
              break;
            }
          }
          if (at < 0) {
            fprintf(stderr, "verify: L row missing from panel\n");
            exit(EXIT_FAILURE);
          }
          mine = W[sn->uw + sn->width + at];
        }
        const double r = rel_diff(mine, rv[p], cn);
        if (r > out->max_rel) out->max_rel = r;
      }
      /* U column */
      lu_get(lu_ref, uip, ulen, k, &ri, &rv, &rl);
      cn = 0.0;
      for (UF_long p = 0; p < rl; ++p) {
        if (fabs(rv[p]) > cn) cn = fabs(rv[p]);
      }
      for (UF_long p = 0; p < rl; ++p) {
        const UF_long row = ri[p];
        double mine = 0.0;
        if (row >= (UF_long)sn->start) {
          mine = W[sn->uw + (row - sn->start)];
        } else {
          const int32_t prod = blk->col2snode[row];
          const snb_edge *edges = blk->edges + sn->edges_off;
          int found = 0;
          for (int32_t e2 = 0; e2 < sn->edge_count; ++e2) {
            if (edges[e2].producer == prod) {
              const int32_t *sub = blk->subcols + edges[e2].sub_off;
              const int32_t want =
                (int32_t)(row - blk->snodes[prod].start);
              int32_t lo = 0, hi = edges[e2].sub_count - 1;
              while (lo <= hi) {
                const int32_t mid = (lo + hi) >> 1;
                if (sub[mid] < want) {
                  lo = mid + 1;
                } else if (sub[mid] > want) {
                  hi = mid - 1;
                } else {
                  mine = W[edges[e2].u_start + mid];
                  found = 1;
                  break;
                }
              }
              break;
            }
          }
          if (!found) {
            fprintf(stderr, "verify: U row %ld of col %ld not in edges\n",
                    (long)row, (long)k);
            exit(EXIT_FAILURE);
          }
        }
        const double r = rel_diff(mine, rv[p], cn);
        if (r > out->max_rel) out->max_rel = r;
      }
    }
  }
  const double *udiag_ref = (const double *)numeric->Udiag;
  for (UF_long k = 0; k < symbolic->n; ++k) {
    const double r = rel_diff(st->udiag[k], udiag_ref[k], 0.0);
    if (r > out->udiag_max_rel) out->udiag_max_rel = r;
  }
  const double *off_ref = (const double *)numeric->Offx;
  for (UF_long p = 0; p < numeric->nzoff; ++p) {
    const double r = rel_diff(st->offx[p], off_ref[p], 0.0);
    if (r > out->off_max_rel) out->off_max_rel = r;
  }
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

int main(int argc, char **argv) {
  if (argc < 2) {
    fprintf(stderr,
            "Usage: %s matrix.mtx [passes=5] [zeta=0.2] [wnarrow=4] "
            "[wmax=32]\n",
            argv[0]);
    return EXIT_FAILURE;
  }
  const char *path = argv[1];
  const int passes = argc > 2 ? atoi(argv[2]) : 5;
  const double zeta = argc > 3 ? atof(argv[3]) : 0.2;
  const int wnarrow = argc > 4 ? atoi(argv[4]) : 4;
  const int wmax = argc > 5 ? atoi(argv[5]) : 32;
  const int postorder = argc > 6 ? atoi(argv[6]) : 0;
  const int scale = getenv("SNB_SCALE") ? atoi(getenv("SNB_SCALE")) : 2;
  const int btf = getenv("SNB_BTF") ? atoi(getenv("SNB_BTF")) : 1;
  const int ordering =
    getenv("SNB_ORDERING") ? atoi(getenv("SNB_ORDERING")) : 0;
  const int trace = getenv("SNB_TRACE") != NULL;

  matrix a;
  memset(&a, 0, sizeof(a));
  if (!read_matrix_market(path, &a)) return EXIT_FAILURE;

  UF_long *ap = (UF_long *)xmalloc(((size_t)a.n + 1) * sizeof(UF_long));
  UF_long *ai = (UF_long *)xmalloc((size_t)a.nnz * sizeof(UF_long));
  for (int64_t i = 0; i <= a.n; ++i) ap[i] = a.col_ptr[i];
  for (int64_t i = 0; i < a.nnz; ++i) ai[i] = a.row_idx[i];

  trilinos_klu_l_common common;
  trilinos_klu_l_defaults(&common);
  common.ordering = ordering;
  common.btf = btf;
  common.scale = scale;

  const double t_analyze = now_seconds();
  trilinos_klu_l_symbolic *symbolic =
    trilinos_klu_l_analyze(a.n, ap, ai, &common);
  if (symbolic == NULL) {
    fprintf(stderr, "analyze failed (%ld)\n", (long)common.status);
    return EXIT_FAILURE;
  }
  const double t_factor = now_seconds();
  trilinos_klu_l_numeric *numeric =
    trilinos_klu_l_factor(ap, ai, a.values, symbolic, &common);
  if (numeric == NULL) {
    fprintf(stderr, "factor failed (%ld)\n", (long)common.status);
    return EXIT_FAILURE;
  }
  const double t_sort = now_seconds();
  if (!trilinos_klu_l_sort(symbolic, numeric, &common)) {
    fprintf(stderr, "sort failed\n");
    return EXIT_FAILURE;
  }

  /* Round 4: etree-postorder amalgamation.  Recompute the ordering so
   * that every L-etree child is adjacent to its parent, then re-analyze
   * with the composed permutations and re-factor.  Widens supernodes
   * from the same fill pattern (nested patterns become consecutive). */
  if (postorder) {
    const UF_long n = a.n;
    UF_long *puser = (UF_long *)xmalloc((size_t)n * sizeof(UF_long));
    UF_long *quser = (UF_long *)xmalloc((size_t)n * sizeof(UF_long));
    UF_long *parent = (UF_long *)xmalloc((size_t)n * sizeof(UF_long));
    UF_long *first_child = (UF_long *)xmalloc((size_t)n * sizeof(UF_long));
    UF_long *next_sibling = (UF_long *)xmalloc((size_t)n * sizeof(UF_long));
    UF_long *stack = (UF_long *)xmalloc((size_t)n * sizeof(UF_long));
    UF_long *perm = (UF_long *)xmalloc((size_t)n * sizeof(UF_long));
    const int64_t old_lnz = numeric->lnz;
    const UF_long old_nblocks = symbolic->nblocks;
    for (UF_long block = 0; block < symbolic->nblocks; ++block) {
      const UF_long k1 = symbolic->R[block];
      const UF_long k2 = symbolic->R[block + 1];
      const UF_long nk = k2 - k1;
      if (nk < 2 || numeric->LUbx[block] == NULL) {
        for (UF_long k = k1; k < k2; ++k) {
          puser[k] = numeric->Pnum[k];
          quser[k] = symbolic->Q[k];
        }
        continue;
      }
      const UF_long *lip = numeric->Lip + k1;
      const UF_long *llen = numeric->Llen + k1;
      double *lu = (double *)numeric->LUbx[block];
      for (UF_long k = 0; k < nk; ++k) {
        UF_long *li;
        double *lx;
        UF_long len;
        lu_get(lu, lip, llen, k, &li, &lx, &len);
        parent[k] = len > 0 ? li[0] : -1; /* sorted: first below row */
        first_child[k] = -1;
        next_sibling[k] = -1;
      }
      /* child lists in increasing order (walk columns backwards) */
      for (UF_long k = nk; k-- > 0;) {
        const UF_long p = parent[k];
        if (p >= 0) {
          next_sibling[k] = first_child[p];
          first_child[p] = k;
        }
      }
      /* iterative DFS postorder over the forest, roots ascending */
      UF_long out = 0;
      for (UF_long root = 0; root < nk; ++root) {
        if (parent[root] >= 0) continue;
        UF_long top = 0;
        stack[top] = root;
        while (top != (UF_long)-1) {
          const UF_long v = stack[top];
          const UF_long child = first_child[v];
          if (child >= 0) {
            first_child[v] = next_sibling[child]; /* consume */
            stack[++top] = child;
          } else {
            perm[out++] = v;
            top--;
          }
        }
      }
      if (out != nk) {
        fprintf(stderr, "postorder: block %ld covered %ld of %ld cols\n",
                (long)block, (long)out, (long)nk);
        return EXIT_FAILURE;
      }
      for (UF_long i = 0; i < nk; ++i) {
        puser[k1 + i] = numeric->Pnum[k1 + perm[i]];
        quser[k1 + i] = symbolic->Q[k1 + perm[i]];
      }
    }
    trilinos_klu_l_free_numeric(&numeric, &common);
    trilinos_klu_l_free_symbolic(&symbolic, &common);
    symbolic = trilinos_klu_l_analyze_given(a.n, ap, ai, puser, quser,
                                            &common);
    if (symbolic == NULL) {
      fprintf(stderr, "postorder analyze_given failed (%ld)\n",
              (long)common.status);
      return EXIT_FAILURE;
    }
    numeric = trilinos_klu_l_factor(ap, ai, a.values, symbolic, &common);
    if (numeric == NULL) {
      fprintf(stderr, "postorder factor failed (%ld)\n",
              (long)common.status);
      return EXIT_FAILURE;
    }
    if (!trilinos_klu_l_sort(symbolic, numeric, &common)) {
      fprintf(stderr, "postorder sort failed\n");
      return EXIT_FAILURE;
    }
    printf("# postorder: nblocks %ld -> %ld, lnz %ld -> %ld (%+.1f%%)\n",
           (long)old_nblocks, (long)symbolic->nblocks, (long)old_lnz,
           (long)numeric->lnz,
           100.0 * ((double)numeric->lnz / (double)old_lnz - 1.0));
    free(puser);
    free(quser);
    free(parent);
    free(first_child);
    free(next_sibling);
    free(stack);
    free(perm);
  }
  const double t_setup = now_seconds();

  printf("# %s n=%ld nnz=%ld nblocks=%ld maxblock=%ld lnz=%ld unz=%ld\n",
         path, (long)a.n, (long)a.nnz, (long)symbolic->nblocks,
         (long)symbolic->maxblock, (long)numeric->lnz, (long)numeric->unz);
  printf("# analyze %.3fs factor %.3fs sort %.3fs\n", t_factor - t_analyze,
         t_sort - t_factor, t_setup - t_sort);

  /* ---------------- landed replica setup ---------------- */
  landed_state landed;
  memset(&landed, 0, sizeof(landed));
  landed.a = &a;
  landed.symbolic = symbolic;
  landed.numeric = numeric;
  landed.scale = scale;
  landed.udiag = (double *)xcalloc((size_t)a.n, sizeof(double));
  landed.offx = (double *)xcalloc(
    (size_t)(symbolic->nzoff > 0 ? symbolic->nzoff : 1), sizeof(double));
  landed.rs = (double *)xcalloc((size_t)a.n, sizeof(double));
  landed.x = (double *)xcalloc((size_t)symbolic->maxblock, sizeof(double));
  landed.lu_copy = (double **)xcalloc((size_t)symbolic->nblocks,
                                      sizeof(double *));
  landed.li32 = (int32_t **)xcalloc((size_t)a.n, sizeof(int32_t *));
  landed.ui32 = (int32_t **)xcalloc((size_t)a.n, sizeof(int32_t *));
  landed.snode_run_end = (UF_long *)xcalloc((size_t)a.n, sizeof(UF_long));
  {
    UF_long covered = 0;
    for (UF_long block = 0; block < symbolic->nblocks; ++block) {
      const UF_long k1 = symbolic->R[block];
      const UF_long k2 = symbolic->R[block + 1];
      const UF_long nk = k2 - k1;
      /* singleton blocks have no LU storage (LUbx NULL, LUsize unset) */
      if (nk < 2 || numeric->LUbx[block] == NULL) continue;
      const size_t bytes = numeric->LUsize[block] * sizeof(double);
      landed.lu_copy[block] = (double *)xmalloc(bytes ? bytes : 1);
      memcpy(landed.lu_copy[block], numeric->LUbx[block], bytes);
      const UF_long *lip = numeric->Lip + k1;
      const UF_long *llen = numeric->Llen + k1;
      const UF_long *uip = numeric->Uip + k1;
      const UF_long *ulen = numeric->Ulen + k1;
      double *lu = landed.lu_copy[block];
      for (UF_long k = 0; k < nk; ++k) {
        UF_long *idx;
        double *val;
        UF_long len;
        lu_get(lu, lip, llen, k, &idx, &val, &len);
        int32_t *m = (int32_t *)xmalloc((size_t)len * sizeof(int32_t));
        for (UF_long p = 0; p < len; ++p) m[p] = (int32_t)idx[p];
        landed.li32[k1 + k] = m;
        lu_get(lu, uip, ulen, k, &idx, &val, &len);
        m = (int32_t *)xmalloc((size_t)len * sizeof(int32_t));
        for (UF_long p = 0; p < len; ++p) m[p] = (int32_t)idx[p];
        landed.ui32[k1 + k] = m;
      }
      /* strict-nesting runs */
      UF_long start_col = 0;
      for (UF_long k = 0; k < nk; ++k) {
        int extends = 0;
        if (k + 1 < nk && llen[k] >= 1 && llen[k + 1] == llen[k] - 1) {
          UF_long *li, *li2;
          double *lx, *lx2;
          UF_long l1, l2;
          lu_get(lu, lip, llen, k, &li, &lx, &l1);
          lu_get(lu, lip, llen, k + 1, &li2, &lx2, &l2);
          if (li[0] == k + 1 &&
              memcmp(li + 1, li2, (size_t)l2 * sizeof(UF_long)) == 0) {
            extends = 1;
          }
        }
        if (!extends) {
          if (k > start_col) {
            for (UF_long c2 = start_col; c2 <= k; ++c2) {
              landed.snode_run_end[k1 + c2] = k1 + k + 1;
            }
            covered += k - start_col + 1;
          }
          start_col = k + 1;
        }
      }
    }
    if (trace) {
      printf("# landed nesting-run coverage: %ld cols\n", (long)covered);
    }
  }

  /* ---------------- supernodal setup ---------------- */
  snb_state snb;
  memset(&snb, 0, sizeof(snb));
  snb.a = &a;
  snb.symbolic = symbolic;
  snb.numeric = numeric;
  snb.scale = scale;
  snb.wnarrow = wnarrow;
  snb.profile = getenv("SNB_PROFILE") != NULL;
  snb.udiag = (double *)xcalloc((size_t)a.n, sizeof(double));
  snb.offx = (double *)xcalloc(
    (size_t)(symbolic->nzoff > 0 ? symbolic->nzoff : 1), sizeof(double));
  snb.rs = (double *)xcalloc((size_t)a.n, sizeof(double));
  const double t_snb0 = now_seconds();
  snb_setup(&snb, zeta, wmax);
  const double t_snb1 = now_seconds();

  build_scatter_programs(&a, symbolic, numeric, &snb, &landed.sp, &snb.sp);

  /* supernode stats */
  {
    int64_t snodes = 0, cols = 0, panel_doubles = 0, edges = 0;
    int64_t true_entries = numeric->lnz + numeric->unz;
    for (UF_long b = 0; b < symbolic->nblocks; ++b) {
      const snb_block *blk = &snb.blocks[b];
      if (blk->nk < 2) continue;
      snodes += blk->snode_count;
      cols += blk->nk;
      panel_doubles += blk->panel_doubles;
      for (int64_t s = 0; s < blk->snode_count; ++s) {
        edges += blk->snodes[s].edge_count;
      }
    }
    printf("# snb: %ld snodes over %ld cols (avg w %.2f), panel %.1f MB, "
           "%.2fx true LU, %ld edges, setup %.3fs\n",
           (long)snodes, (long)cols,
           snodes > 0 ? (double)cols / (double)snodes : 0.0,
           (double)panel_doubles * 8.0 / 1048576.0,
           true_entries > 0 ? (double)panel_doubles / (double)true_entries
                            : 0.0,
           (long)edges, t_snb1 - t_snb0);
  }

  /* ---------------- verification ---------------- */
  /* Direct value comparison is chaotic on ill-conditioned inputs
   * (mc2depi/rajat25 churn O(1) between ANY two kernels), so the gate is
   * a componentwise reconstruction residual, which is summation-order
   * independent.  Direct diffs are still printed as information. */
  if (!trilinos_klu_l_refactor(ap, ai, a.values, symbolic, numeric,
                               &common)) {
    fprintf(stderr, "klu refactor failed\n");
    return EXIT_FAILURE;
  }
  landed_refactor(&landed);
  snb_refactor(&snb, 1);
  verify_result vl, vs;
  verify_landed(&landed, &vl);
  verify_snb(&snb, &vs);
  printf("# direct diff landed: LU %.3g udiag %.3g off %.3g\n", vl.max_rel,
         vl.udiag_max_rel, vl.off_max_rel);
  printf("# direct diff snb:    LU %.3g udiag %.3g off %.3g\n", vs.max_rel,
         vs.udiag_max_rel, vs.off_max_rel);
  if (trace) {
    printf("# snb scatter hit %ld miss %ld (%.1f%% wasted)\n",
           (long)snb.scatter_hit, (long)snb.scatter_miss,
           snb.scatter_hit + snb.scatter_miss > 0
             ? 100.0 * (double)snb.scatter_miss /
                 (double)(snb.scatter_hit + snb.scatter_miss)
             : 0.0);
  }
  double err_klu, err_landed, err_snb;
  {
    /* klu's Rs is permuted to pivot order after refactor; the checker
     * wants old-row indexing */
    double *rs_unperm = (double *)xcalloc((size_t)a.n, sizeof(double));
    if (scale > 0 && numeric->Rs != NULL) {
      for (UF_long r = 0; r < a.n; ++r) {
        rs_unperm[r] = numeric->Rs[numeric->Pinv[r]];
      }
    }
    double **klu_lu = (double **)xcalloc((size_t)symbolic->nblocks,
                                         sizeof(double *));
    double **snb_lu = (double **)xcalloc((size_t)symbolic->nblocks,
                                         sizeof(double *));
    for (UF_long b = 0; b < symbolic->nblocks; ++b) {
      klu_lu[b] = (double *)numeric->LUbx[b];
      if (landed.lu_copy[b] != NULL) {
        const size_t bytes = numeric->LUsize[b] * sizeof(double);
        snb_lu[b] = (double *)xmalloc(bytes ? bytes : 1);
        memcpy(snb_lu[b], numeric->LUbx[b], bytes);
      }
    }
    snb_export(&snb, snb_lu);
    err_klu = reconstruction_error(&a, symbolic, numeric, klu_lu,
                                   (const double *)numeric->Udiag,
                                   rs_unperm, scale, &landed.sp);
    err_landed = reconstruction_error(&a, symbolic, numeric,
                                      landed.lu_copy, landed.udiag,
                                      landed.rs, scale, &landed.sp);
    err_snb = reconstruction_error(&a, symbolic, numeric, snb_lu,
                                   snb.udiag, snb.rs, scale, &landed.sp);
    for (UF_long b = 0; b < symbolic->nblocks; ++b) free(snb_lu[b]);
    free(snb_lu);
    free(klu_lu);
    free(rs_unperm);
  }
  printf("# residual klu %.3g landed %.3g snb %.3g\n", err_klu, err_landed,
         err_snb);
  const double tol = 1e-9;
  if (err_landed > tol || err_snb > tol) {
    fprintf(stderr, "VERIFICATION FAILED (residual tol %g, klu floor %g)\n",
            tol, err_klu);
    return EXIT_FAILURE;
  }

  /* ---------------- timing (interleaved) ---------------- */
  double best_klu = 1e300, best_landed = 1e300, best_snb = 1e300;
  printf("# pass  klu        landed     snb        landed/snb\n");
  for (int pass = 0; pass < passes; ++pass) {
    double t0 = now_seconds();
    trilinos_klu_l_refactor(ap, ai, a.values, symbolic, numeric, &common);
    double t1 = now_seconds();
    landed_refactor(&landed);
    double t2 = now_seconds();
    snb_refactor(&snb, 0);
    double t3 = now_seconds();
    const double dk = t1 - t0, dl = t2 - t1, ds = t3 - t2;
    if (dk < best_klu) best_klu = dk;
    if (dl < best_landed) best_landed = dl;
    if (ds < best_snb) best_snb = ds;
    printf("  %2d    %.6f   %.6f   %.6f   %.3f\n", pass, dk, dl, ds,
           ds > 0 ? dl / ds : 0.0);
  }
  printf("RESULT %s klu %.6f landed %.6f snb %.6f speedup_vs_landed %.3f "
         "speedup_vs_klu %.3f\n",
         path, best_klu, best_landed, best_snb,
         best_snb > 0 ? best_landed / best_snb : 0.0,
         best_snb > 0 ? best_klu / best_snb : 0.0);
  if (snb.profile) {
    uint64_t total = 0;
    for (int i = 0; i < 5; ++i) total += snb.phase_cycles[i];
    if (total > 0) {
      printf("# snb phases: scatter %.1f%% narrow %.1f%% trsm %.1f%% "
             "gemm %.1f%% getrf %.1f%%\n",
             100.0 * (double)snb.phase_cycles[0] / (double)total,
             100.0 * (double)snb.phase_cycles[1] / (double)total,
             100.0 * (double)snb.phase_cycles[2] / (double)total,
             100.0 * (double)snb.phase_cycles[3] / (double)total,
             100.0 * (double)snb.phase_cycles[4] / (double)total);
    }
  }

  /* cleanup */
  for (UF_long k = 0; k < a.n; ++k) {
    free(landed.li32[k]);
    free(landed.ui32[k]);
  }
  free(landed.li32);
  free(landed.ui32);
  for (UF_long b = 0; b < symbolic->nblocks; ++b) free(landed.lu_copy[b]);
  free(landed.lu_copy);
  free(landed.udiag);
  free(landed.offx);
  free(landed.rs);
  free(landed.x);
  free(landed.snode_run_end);
  scatter_program_free(&landed.sp);
  free(snb.udiag);
  free(snb.offx);
  free(snb.rs);
  snb_free(&snb);
  trilinos_klu_l_free_numeric(&numeric, &common);
  trilinos_klu_l_free_symbolic(&symbolic, &common);
  free(ap);
  free(ai);
  free(a.col_ptr);
  free(a.row_idx);
  free(a.values);
  return EXIT_SUCCESS;
}
