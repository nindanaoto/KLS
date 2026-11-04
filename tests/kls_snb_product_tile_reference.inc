static __attribute__((noinline)) void reference_snb_dense_product(
  const double *L, int32_t ld, const int32_t *sub, int32_t sc,
  const double *B, int32_t bld, int32_t w, int32_t r0, int32_t r1, double *T) {
  const int32_t nr = r1 - r0;
  int32_t c = 0;
#if defined(__AVX512F__)
  for (; c + 3 < w; c += 4) {
    const double *b0 = B + (size_t)c * bld;
    const double *b1 = b0 + bld, *b2 = b1 + bld, *b3 = b2 + bld;
    int32_t r = r0;
    for (; r + 7 < r1; r += 8) {
      __m512d l = _mm512_loadu_pd(L + (size_t)sub[0] * ld + r);
      __m512d t0 = _mm512_mul_pd(l, _mm512_set1_pd(b0[0]));
      __m512d t1 = _mm512_mul_pd(l, _mm512_set1_pd(b1[0]));
      __m512d t2 = _mm512_mul_pd(l, _mm512_set1_pd(b2[0]));
      __m512d t3 = _mm512_mul_pd(l, _mm512_set1_pd(b3[0]));
      for (int32_t i = 1; i < sc; ++i) {
        l = _mm512_loadu_pd(L + (size_t)sub[i] * ld + r);
        t0 = _mm512_add_pd(t0, _mm512_mul_pd(l, _mm512_set1_pd(b0[i])));
        t1 = _mm512_add_pd(t1, _mm512_mul_pd(l, _mm512_set1_pd(b1[i])));
        t2 = _mm512_add_pd(t2, _mm512_mul_pd(l, _mm512_set1_pd(b2[i])));
        t3 = _mm512_add_pd(t3, _mm512_mul_pd(l, _mm512_set1_pd(b3[i])));
      }
      double *out = T + (size_t)c * nr + r - r0;
      _mm512_storeu_pd(out, t0);
      _mm512_storeu_pd(out + nr, t1);
      _mm512_storeu_pd(out + 2u * nr, t2);
      _mm512_storeu_pd(out + 3u * nr, t3);
    }
    if (r < r1) {
      const __mmask8 mask = (__mmask8)((1u << (r1-r))-1u);
      __m512d l = _mm512_maskz_loadu_pd(mask,L+(size_t)sub[0]*ld+r);
      __m512d t0 = _mm512_maskz_mul_pd(mask,l,_mm512_set1_pd(b0[0]));
      __m512d t1 = _mm512_maskz_mul_pd(mask,l,_mm512_set1_pd(b1[0]));
      __m512d t2 = _mm512_maskz_mul_pd(mask,l,_mm512_set1_pd(b2[0]));
      __m512d t3 = _mm512_maskz_mul_pd(mask,l,_mm512_set1_pd(b3[0]));
      for (int32_t i = 1; i < sc; ++i) {
        l = _mm512_maskz_loadu_pd(mask,L+(size_t)sub[i]*ld+r);
        t0 = _mm512_add_pd(t0,_mm512_maskz_mul_pd(mask,l,_mm512_set1_pd(b0[i])));
        t1 = _mm512_add_pd(t1,_mm512_maskz_mul_pd(mask,l,_mm512_set1_pd(b1[i])));
        t2 = _mm512_add_pd(t2,_mm512_maskz_mul_pd(mask,l,_mm512_set1_pd(b2[i])));
        t3 = _mm512_add_pd(t3,_mm512_maskz_mul_pd(mask,l,_mm512_set1_pd(b3[i])));
      }
      double *out = T + (size_t)c * nr + r - r0;
      _mm512_mask_storeu_pd(out,mask,t0);
      _mm512_mask_storeu_pd(out+nr,mask,t1);
      _mm512_mask_storeu_pd(out+2u*nr,mask,t2);
      _mm512_mask_storeu_pd(out+3u*nr,mask,t3);
    }
  }
#endif
  for (; c < w; ++c) {
    const double *b = B + (size_t)c * bld;
#if defined(__AVX512F__)
    int32_t r = r0;
    for (; r + 7 < r1; r += 8) {
      __m512d t = _mm512_mul_pd(
        _mm512_loadu_pd(L + (size_t)sub[0] * ld + r), _mm512_set1_pd(b[0]));
      for (int32_t i = 1; i < sc; ++i)
        t = _mm512_add_pd(t, _mm512_mul_pd(
          _mm512_loadu_pd(L + (size_t)sub[i] * ld + r), _mm512_set1_pd(b[i])));
      _mm512_storeu_pd(T + (size_t)c * nr + r - r0, t);
    }
    if (r < r1) {
      const __mmask8 mask = (__mmask8)((1u << (r1-r))-1u);
      __m512d t = _mm512_maskz_mul_pd(mask,
        _mm512_maskz_loadu_pd(mask,L+(size_t)sub[0]*ld+r),_mm512_set1_pd(b[0]));
      for (int32_t i = 1; i < sc; ++i)
        t = _mm512_add_pd(t,_mm512_maskz_mul_pd(mask,
          _mm512_maskz_loadu_pd(mask,L+(size_t)sub[i]*ld+r),_mm512_set1_pd(b[i])));
      _mm512_mask_storeu_pd(T+(size_t)c*nr+r-r0,mask,t);
    }
#else
    /* Portable column update: contiguous L and T streams allow the compiler
       to use its available SIMD width without a row-major B mirror. */
    double *out = T + (size_t)c * nr;
    const double *l = L + (size_t)sub[0] * ld + r0;
#pragma omp simd
    for (int32_t row = 0; row < nr; ++row) out[row] = l[row] * b[0];
    for (int32_t i = 1; i < sc; ++i) {
      l = L + (size_t)sub[i] * ld + r0;
#pragma omp simd
      for (int32_t row = 0; row < nr; ++row) out[row] += l[row] * b[i];
    }
#endif
  }
}
