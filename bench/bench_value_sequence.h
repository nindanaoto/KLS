#ifndef KLS_BENCH_VALUE_SEQUENCE_H
#define KLS_BENCH_VALUE_SEQUENCE_H

#include <stdint.h>
#include <string.h>

typedef enum bench_refactor_value_mode {
  BENCH_REFACTOR_VALUES_UNCHANGED = 0,
  BENCH_REFACTOR_VALUES_RANK_PRESERVING = 1
} bench_refactor_value_mode;

static inline int bench_parse_refactor_value_mode(
    const char *name, bench_refactor_value_mode *mode_out) {
  if (name == NULL || mode_out == NULL) return 0;
  if (strcmp(name, "unchanged") == 0) {
    *mode_out = BENCH_REFACTOR_VALUES_UNCHANGED;
    return 1;
  }
  if (strcmp(name, "rank-preserving") == 0) {
    *mode_out = BENCH_REFACTOR_VALUES_RANK_PRESERVING;
    return 1;
  }
  return 0;
}

static inline const char *bench_refactor_value_mode_name(
    bench_refactor_value_mode mode) {
  return mode == BENCH_REFACTOR_VALUES_RANK_PRESERVING
           ? "rank-preserving" : "unchanged";
}

static inline uint64_t bench_value_mix64(uint64_t value) {
  value += UINT64_C(0x9e3779b97f4a7c15);
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

static inline double bench_axis_multiplier(uint64_t index,
                                           uint64_t generation,
                                           uint64_t salt,
                                           double amplitude) {
  const uint64_t key =
    (index + UINT64_C(1)) * UINT64_C(0xd6e8feb86659fd93) ^
    (generation + UINT64_C(1)) * UINT64_C(0xa0761d6478bd642f) ^ salt;
  const uint64_t bucket = bench_value_mix64(key) >> 48;
  const double signed_unit =
    ((double)bucket - 32767.5) * (1.0 / 32767.5);
  return 1.0 + amplitude * signed_unit;
}

/*
 * Define generation g as A_g = D_row(g) A_0 D_col(g).  For amplitude < 1
 * both diagonal matrices are nonsingular, so every generated matrix has
 * exactly the same rank as A_0.  The sequence depends only on coordinates
 * and generation, making CSC and CSR harnesses bitwise consistent without
 * using any matrix- or corpus-specific tuning.
 */
static inline double bench_refactor_value(double base_value,
                                          uint64_t row,
                                          uint64_t col,
                                          uint64_t generation,
                                          double amplitude) {
  const double row_scale = bench_axis_multiplier(
    row, generation, UINT64_C(0xe7037ed1a0b428db), amplitude);
  const double col_scale = bench_axis_multiplier(
    col, generation, UINT64_C(0x8ebc6af09c88c6e3), amplitude);
  return (base_value * row_scale) * col_scale;
}

#endif
