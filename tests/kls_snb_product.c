/* Deterministic private-kernel coverage, independent of policy timing and
   circuit matrices. The test executable exports no additional library API. */
#include "../src/kls.c"

int main(void) {
  const int rows[] = {0, 1, 7, 8, 9, 15, 16, 17, 63, 64, 65};
  const int widths[] = {1, 3, 4, 5, 7, 8, 9, 31, 32};
  const int counts[] = {1, 2, 7, 8, 17, 32};
  double L[73 * 37], B[32 * 32], scratch[65 * 32 + 2];
  int32_t sub[32];
  for (int p = 0; p < 73 * 37; ++p) L[p] = ((p * 17) % 29 - 14) / 11.0;
  for (int p = 0; p < 32 * 32; ++p) B[p] = ((p * 13) % 23 - 11) / 7.0;
  for (int i = 0; i < 32; ++i) sub[i] = (i * 7) % 37;
  for (unsigned a = 0; a < sizeof(rows) / sizeof(rows[0]); ++a)
    for (unsigned b = 0; b < sizeof(widths) / sizeof(widths[0]); ++b)
      for (unsigned d = 0; d < sizeof(counts) / sizeof(counts[0]); ++d) {
        const int nr = rows[a], w = widths[b], sc = counts[d];
        double *T = scratch + 1;
        scratch[0] = T[nr * w] = 12345.0;
        kls_snb_dense_product(L, 73, sub, sc, B, w, 3, 3 + nr, T);
        if (scratch[0] != 12345.0 || T[nr * w] != 12345.0) return 1;
        for (int c = 0; c < w; ++c)
          for (int r = 0; r < nr; ++r) {
            long double expected = 0.0;
            for (int i = 0; i < sc; ++i)
              expected += (long double)L[sub[i] * 73 + 3 + r] * B[c * sc + i];
            if (!isfinite(T[c * nr + r]) ||
                fabsl(T[c * nr + r] - expected) > 1e-12L) {
              fprintf(stderr, "SNB product mismatch rows=%d width=%d count=%d\n", nr, w, sc);
              return 1;
            }
          }
      }
  return 0;
}
