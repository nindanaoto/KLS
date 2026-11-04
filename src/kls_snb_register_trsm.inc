/* Solve sixteen producer rows for four independent target columns.
 * Keep current B in registers, preserving each column's producer order and
 * the incumbent skip for zero multipliers. */
#if defined(__AVX512F__) && defined(__AVX512VL__)
#if defined(__GNUC__) || defined(__clang__)
#define KLS_SNB_REGISTER_NOINLINE __attribute__((noinline,hot))
#elif defined(_MSC_VER)
#define KLS_SNB_REGISTER_NOINLINE __declspec(noinline)
#else
#define KLS_SNB_REGISTER_NOINLINE
#endif
static KLS_SNB_REGISTER_NOINLINE void kls_snb_register_trsm_16(
    double *restrict B,int H,const double *restrict L,int ld,int start) {
#ifdef KLS_SNB_REGISTER_TEST_EVENT
  KLS_SNB_REGISTER_TEST_EVENT(16);
#endif
  const __m256d c0_0=_mm256_loadu_pd(B+0*(size_t)H+0);
  const __m256d c0_1=_mm256_loadu_pd(B+1*(size_t)H+0);
  const __m256d c0_2=_mm256_loadu_pd(B+2*(size_t)H+0);
  const __m256d c0_3=_mm256_loadu_pd(B+3*(size_t)H+0);
  const __m256d a0=_mm256_unpacklo_pd(c0_0,c0_1);
  const __m256d b0=_mm256_unpackhi_pd(c0_0,c0_1);
  const __m256d c0=_mm256_unpacklo_pd(c0_2,c0_3);
  const __m256d d0=_mm256_unpackhi_pd(c0_2,c0_3);
  __m256d r0=_mm256_permute2f128_pd(a0,c0,0x20);
  __m256d r1=_mm256_permute2f128_pd(b0,d0,0x20);
  __m256d r2=_mm256_permute2f128_pd(a0,c0,0x31);
  __m256d r3=_mm256_permute2f128_pd(b0,d0,0x31);
  const __m256d c4_0=_mm256_loadu_pd(B+0*(size_t)H+4);
  const __m256d c4_1=_mm256_loadu_pd(B+1*(size_t)H+4);
  const __m256d c4_2=_mm256_loadu_pd(B+2*(size_t)H+4);
  const __m256d c4_3=_mm256_loadu_pd(B+3*(size_t)H+4);
  const __m256d a4=_mm256_unpacklo_pd(c4_0,c4_1);
  const __m256d b4=_mm256_unpackhi_pd(c4_0,c4_1);
  const __m256d c4=_mm256_unpacklo_pd(c4_2,c4_3);
  const __m256d d4=_mm256_unpackhi_pd(c4_2,c4_3);
  __m256d r4=_mm256_permute2f128_pd(a4,c4,0x20);
  __m256d r5=_mm256_permute2f128_pd(b4,d4,0x20);
  __m256d r6=_mm256_permute2f128_pd(a4,c4,0x31);
  __m256d r7=_mm256_permute2f128_pd(b4,d4,0x31);
  const __m256d c8_0=_mm256_loadu_pd(B+0*(size_t)H+8);
  const __m256d c8_1=_mm256_loadu_pd(B+1*(size_t)H+8);
  const __m256d c8_2=_mm256_loadu_pd(B+2*(size_t)H+8);
  const __m256d c8_3=_mm256_loadu_pd(B+3*(size_t)H+8);
  const __m256d a8=_mm256_unpacklo_pd(c8_0,c8_1);
  const __m256d b8=_mm256_unpackhi_pd(c8_0,c8_1);
  const __m256d c8=_mm256_unpacklo_pd(c8_2,c8_3);
  const __m256d d8=_mm256_unpackhi_pd(c8_2,c8_3);
  __m256d r8=_mm256_permute2f128_pd(a8,c8,0x20);
  __m256d r9=_mm256_permute2f128_pd(b8,d8,0x20);
  __m256d r10=_mm256_permute2f128_pd(a8,c8,0x31);
  __m256d r11=_mm256_permute2f128_pd(b8,d8,0x31);
  const __m256d c12_0=_mm256_loadu_pd(B+0*(size_t)H+12);
  const __m256d c12_1=_mm256_loadu_pd(B+1*(size_t)H+12);
  const __m256d c12_2=_mm256_loadu_pd(B+2*(size_t)H+12);
  const __m256d c12_3=_mm256_loadu_pd(B+3*(size_t)H+12);
  const __m256d a12=_mm256_unpacklo_pd(c12_0,c12_1);
  const __m256d b12=_mm256_unpackhi_pd(c12_0,c12_1);
  const __m256d c12=_mm256_unpacklo_pd(c12_2,c12_3);
  const __m256d d12=_mm256_unpackhi_pd(c12_2,c12_3);
  __m256d r12=_mm256_permute2f128_pd(a12,c12,0x20);
  __m256d r13=_mm256_permute2f128_pd(b12,d12,0x20);
  __m256d r14=_mm256_permute2f128_pd(a12,c12,0x31);
  __m256d r15=_mm256_permute2f128_pd(b12,d12,0x31);
  const __mmask8 m0=_mm256_cmp_pd_mask(r0,_mm256_setzero_pd(),_CMP_NEQ_UQ);
  if(m0) {
    const double *li=L+(size_t)(start+0)*ld+start;
    r1=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[1]),r0,r1,m0);
    r2=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[2]),r0,r2,m0);
    r3=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[3]),r0,r3,m0);
    r4=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[4]),r0,r4,m0);
    r5=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[5]),r0,r5,m0);
    r6=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[6]),r0,r6,m0);
    r7=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[7]),r0,r7,m0);
    r8=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[8]),r0,r8,m0);
    r9=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[9]),r0,r9,m0);
    r10=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[10]),r0,r10,m0);
    r11=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[11]),r0,r11,m0);
    r12=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[12]),r0,r12,m0);
    r13=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[13]),r0,r13,m0);
    r14=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[14]),r0,r14,m0);
    r15=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[15]),r0,r15,m0);
  }
  const __mmask8 m1=_mm256_cmp_pd_mask(r1,_mm256_setzero_pd(),_CMP_NEQ_UQ);
  if(m1) {
    const double *li=L+(size_t)(start+1)*ld+start;
    r2=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[2]),r1,r2,m1);
    r3=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[3]),r1,r3,m1);
    r4=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[4]),r1,r4,m1);
    r5=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[5]),r1,r5,m1);
    r6=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[6]),r1,r6,m1);
    r7=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[7]),r1,r7,m1);
    r8=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[8]),r1,r8,m1);
    r9=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[9]),r1,r9,m1);
    r10=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[10]),r1,r10,m1);
    r11=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[11]),r1,r11,m1);
    r12=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[12]),r1,r12,m1);
    r13=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[13]),r1,r13,m1);
    r14=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[14]),r1,r14,m1);
    r15=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[15]),r1,r15,m1);
  }
  const __mmask8 m2=_mm256_cmp_pd_mask(r2,_mm256_setzero_pd(),_CMP_NEQ_UQ);
  if(m2) {
    const double *li=L+(size_t)(start+2)*ld+start;
    r3=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[3]),r2,r3,m2);
    r4=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[4]),r2,r4,m2);
    r5=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[5]),r2,r5,m2);
    r6=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[6]),r2,r6,m2);
    r7=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[7]),r2,r7,m2);
    r8=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[8]),r2,r8,m2);
    r9=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[9]),r2,r9,m2);
    r10=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[10]),r2,r10,m2);
    r11=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[11]),r2,r11,m2);
    r12=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[12]),r2,r12,m2);
    r13=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[13]),r2,r13,m2);
    r14=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[14]),r2,r14,m2);
    r15=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[15]),r2,r15,m2);
  }
  const __mmask8 m3=_mm256_cmp_pd_mask(r3,_mm256_setzero_pd(),_CMP_NEQ_UQ);
  if(m3) {
    const double *li=L+(size_t)(start+3)*ld+start;
    r4=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[4]),r3,r4,m3);
    r5=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[5]),r3,r5,m3);
    r6=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[6]),r3,r6,m3);
    r7=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[7]),r3,r7,m3);
    r8=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[8]),r3,r8,m3);
    r9=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[9]),r3,r9,m3);
    r10=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[10]),r3,r10,m3);
    r11=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[11]),r3,r11,m3);
    r12=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[12]),r3,r12,m3);
    r13=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[13]),r3,r13,m3);
    r14=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[14]),r3,r14,m3);
    r15=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[15]),r3,r15,m3);
  }
  const __mmask8 m4=_mm256_cmp_pd_mask(r4,_mm256_setzero_pd(),_CMP_NEQ_UQ);
  if(m4) {
    const double *li=L+(size_t)(start+4)*ld+start;
    r5=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[5]),r4,r5,m4);
    r6=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[6]),r4,r6,m4);
    r7=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[7]),r4,r7,m4);
    r8=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[8]),r4,r8,m4);
    r9=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[9]),r4,r9,m4);
    r10=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[10]),r4,r10,m4);
    r11=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[11]),r4,r11,m4);
    r12=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[12]),r4,r12,m4);
    r13=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[13]),r4,r13,m4);
    r14=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[14]),r4,r14,m4);
    r15=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[15]),r4,r15,m4);
  }
  const __mmask8 m5=_mm256_cmp_pd_mask(r5,_mm256_setzero_pd(),_CMP_NEQ_UQ);
  if(m5) {
    const double *li=L+(size_t)(start+5)*ld+start;
    r6=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[6]),r5,r6,m5);
    r7=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[7]),r5,r7,m5);
    r8=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[8]),r5,r8,m5);
    r9=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[9]),r5,r9,m5);
    r10=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[10]),r5,r10,m5);
    r11=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[11]),r5,r11,m5);
    r12=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[12]),r5,r12,m5);
    r13=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[13]),r5,r13,m5);
    r14=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[14]),r5,r14,m5);
    r15=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[15]),r5,r15,m5);
  }
  const __mmask8 m6=_mm256_cmp_pd_mask(r6,_mm256_setzero_pd(),_CMP_NEQ_UQ);
  if(m6) {
    const double *li=L+(size_t)(start+6)*ld+start;
    r7=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[7]),r6,r7,m6);
    r8=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[8]),r6,r8,m6);
    r9=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[9]),r6,r9,m6);
    r10=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[10]),r6,r10,m6);
    r11=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[11]),r6,r11,m6);
    r12=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[12]),r6,r12,m6);
    r13=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[13]),r6,r13,m6);
    r14=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[14]),r6,r14,m6);
    r15=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[15]),r6,r15,m6);
  }
  const __mmask8 m7=_mm256_cmp_pd_mask(r7,_mm256_setzero_pd(),_CMP_NEQ_UQ);
  if(m7) {
    const double *li=L+(size_t)(start+7)*ld+start;
    r8=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[8]),r7,r8,m7);
    r9=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[9]),r7,r9,m7);
    r10=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[10]),r7,r10,m7);
    r11=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[11]),r7,r11,m7);
    r12=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[12]),r7,r12,m7);
    r13=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[13]),r7,r13,m7);
    r14=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[14]),r7,r14,m7);
    r15=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[15]),r7,r15,m7);
  }
  const __mmask8 m8=_mm256_cmp_pd_mask(r8,_mm256_setzero_pd(),_CMP_NEQ_UQ);
  if(m8) {
    const double *li=L+(size_t)(start+8)*ld+start;
    r9=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[9]),r8,r9,m8);
    r10=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[10]),r8,r10,m8);
    r11=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[11]),r8,r11,m8);
    r12=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[12]),r8,r12,m8);
    r13=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[13]),r8,r13,m8);
    r14=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[14]),r8,r14,m8);
    r15=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[15]),r8,r15,m8);
  }
  const __mmask8 m9=_mm256_cmp_pd_mask(r9,_mm256_setzero_pd(),_CMP_NEQ_UQ);
  if(m9) {
    const double *li=L+(size_t)(start+9)*ld+start;
    r10=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[10]),r9,r10,m9);
    r11=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[11]),r9,r11,m9);
    r12=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[12]),r9,r12,m9);
    r13=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[13]),r9,r13,m9);
    r14=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[14]),r9,r14,m9);
    r15=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[15]),r9,r15,m9);
  }
  const __mmask8 m10=_mm256_cmp_pd_mask(r10,_mm256_setzero_pd(),_CMP_NEQ_UQ);
  if(m10) {
    const double *li=L+(size_t)(start+10)*ld+start;
    r11=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[11]),r10,r11,m10);
    r12=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[12]),r10,r12,m10);
    r13=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[13]),r10,r13,m10);
    r14=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[14]),r10,r14,m10);
    r15=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[15]),r10,r15,m10);
  }
  const __mmask8 m11=_mm256_cmp_pd_mask(r11,_mm256_setzero_pd(),_CMP_NEQ_UQ);
  if(m11) {
    const double *li=L+(size_t)(start+11)*ld+start;
    r12=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[12]),r11,r12,m11);
    r13=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[13]),r11,r13,m11);
    r14=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[14]),r11,r14,m11);
    r15=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[15]),r11,r15,m11);
  }
  const __mmask8 m12=_mm256_cmp_pd_mask(r12,_mm256_setzero_pd(),_CMP_NEQ_UQ);
  if(m12) {
    const double *li=L+(size_t)(start+12)*ld+start;
    r13=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[13]),r12,r13,m12);
    r14=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[14]),r12,r14,m12);
    r15=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[15]),r12,r15,m12);
  }
  const __mmask8 m13=_mm256_cmp_pd_mask(r13,_mm256_setzero_pd(),_CMP_NEQ_UQ);
  if(m13) {
    const double *li=L+(size_t)(start+13)*ld+start;
    r14=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[14]),r13,r14,m13);
    r15=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[15]),r13,r15,m13);
  }
  const __mmask8 m14=_mm256_cmp_pd_mask(r14,_mm256_setzero_pd(),_CMP_NEQ_UQ);
  if(m14) {
    const double *li=L+(size_t)(start+14)*ld+start;
    r15=_mm256_mask3_fnmadd_pd(_mm256_set1_pd(li[15]),r14,r15,m14);
  }
  const __m256d oa0=_mm256_unpacklo_pd(r0,r1);
  const __m256d ob0=_mm256_unpackhi_pd(r0,r1);
  const __m256d oc0=_mm256_unpacklo_pd(r2,r3);
  const __m256d od0=_mm256_unpackhi_pd(r2,r3);
  _mm256_storeu_pd(B+0*(size_t)H+0,_mm256_permute2f128_pd(oa0,oc0,0x20));
  _mm256_storeu_pd(B+1*(size_t)H+0,_mm256_permute2f128_pd(ob0,od0,0x20));
  _mm256_storeu_pd(B+2*(size_t)H+0,_mm256_permute2f128_pd(oa0,oc0,0x31));
  _mm256_storeu_pd(B+3*(size_t)H+0,_mm256_permute2f128_pd(ob0,od0,0x31));
  const __m256d oa4=_mm256_unpacklo_pd(r4,r5);
  const __m256d ob4=_mm256_unpackhi_pd(r4,r5);
  const __m256d oc4=_mm256_unpacklo_pd(r6,r7);
  const __m256d od4=_mm256_unpackhi_pd(r6,r7);
  _mm256_storeu_pd(B+0*(size_t)H+4,_mm256_permute2f128_pd(oa4,oc4,0x20));
  _mm256_storeu_pd(B+1*(size_t)H+4,_mm256_permute2f128_pd(ob4,od4,0x20));
  _mm256_storeu_pd(B+2*(size_t)H+4,_mm256_permute2f128_pd(oa4,oc4,0x31));
  _mm256_storeu_pd(B+3*(size_t)H+4,_mm256_permute2f128_pd(ob4,od4,0x31));
  const __m256d oa8=_mm256_unpacklo_pd(r8,r9);
  const __m256d ob8=_mm256_unpackhi_pd(r8,r9);
  const __m256d oc8=_mm256_unpacklo_pd(r10,r11);
  const __m256d od8=_mm256_unpackhi_pd(r10,r11);
  _mm256_storeu_pd(B+0*(size_t)H+8,_mm256_permute2f128_pd(oa8,oc8,0x20));
  _mm256_storeu_pd(B+1*(size_t)H+8,_mm256_permute2f128_pd(ob8,od8,0x20));
  _mm256_storeu_pd(B+2*(size_t)H+8,_mm256_permute2f128_pd(oa8,oc8,0x31));
  _mm256_storeu_pd(B+3*(size_t)H+8,_mm256_permute2f128_pd(ob8,od8,0x31));
  const __m256d oa12=_mm256_unpacklo_pd(r12,r13);
  const __m256d ob12=_mm256_unpackhi_pd(r12,r13);
  const __m256d oc12=_mm256_unpacklo_pd(r14,r15);
  const __m256d od12=_mm256_unpackhi_pd(r14,r15);
  _mm256_storeu_pd(B+0*(size_t)H+12,_mm256_permute2f128_pd(oa12,oc12,0x20));
  _mm256_storeu_pd(B+1*(size_t)H+12,_mm256_permute2f128_pd(ob12,od12,0x20));
  _mm256_storeu_pd(B+2*(size_t)H+12,_mm256_permute2f128_pd(oa12,oc12,0x31));
  _mm256_storeu_pd(B+3*(size_t)H+12,_mm256_permute2f128_pd(ob12,od12,0x31));
}
#undef KLS_SNB_REGISTER_NOINLINE
#endif
