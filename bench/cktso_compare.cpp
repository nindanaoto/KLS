#include "cktso.h"
#include "bench_value_sequence.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <string>
#include <vector>

struct Entry {
  int row;
  int col;
  double value;
};

struct Matrix {
  int n = 0;
  std::vector<int> col_ptr;
  std::vector<int> row_idx;
  std::vector<double> values;
};

static bool read_matrix_market(const char *path, Matrix &a) {
  FILE *fp = std::fopen(path, "r");
  if (fp == nullptr) {
    std::perror(path);
    return false;
  }

  char line[4096];
  if (std::fgets(line, sizeof(line), fp) == nullptr ||
      std::strncmp(line, "%%MatrixMarket matrix coordinate", 32) != 0) {
    std::fprintf(stderr, "%s is not a coordinate MatrixMarket file\n", path);
    std::fclose(fp);
    return false;
  }
  const bool symmetric = std::strstr(line, "symmetric") != nullptr;
  const bool pattern = std::strstr(line, "pattern") != nullptr;
  const bool complex = std::strstr(line, "complex") != nullptr;
  if (complex) {
    std::fprintf(stderr, "complex matrices are not supported by this compare tool\n");
    std::fclose(fp);
    return false;
  }

  do {
    if (std::fgets(line, sizeof(line), fp) == nullptr) {
      std::fclose(fp);
      return false;
    }
  } while (line[0] == '%');

  long long rows = 0;
  long long cols = 0;
  long long entries = 0;
  if (std::sscanf(line, "%lld %lld %lld", &rows, &cols, &entries) != 3 ||
      rows <= 0 || rows != cols || rows > INT_MAX || entries < 0) {
    std::fprintf(stderr, "%s must be a square matrix fitting 32-bit CKTSO\n", path);
    std::fclose(fp);
    return false;
  }

  std::vector<Entry> items;
  items.reserve(static_cast<size_t>(entries) * (symmetric ? 2u : 1u));
  for (long long k = 0; k < entries; ++k) {
    if (std::fgets(line, sizeof(line), fp) == nullptr) {
      std::fclose(fp);
      return false;
    }
    int row = 0;
    int col = 0;
    double value = 1.0;
    const int parsed = pattern
      ? std::sscanf(line, "%d %d", &row, &col)
      : std::sscanf(line, "%d %d %lf", &row, &col, &value);
    if (parsed < 2 || row <= 0 || col <= 0 || row > rows || col > cols) {
      std::fprintf(stderr, "invalid MatrixMarket entry near %lld\n", k + 1);
      std::fclose(fp);
      return false;
    }
    --row;
    --col;
    items.push_back({row, col, value});
    if (symmetric && row != col) {
      items.push_back({col, row, value});
    }
  }
  std::fclose(fp);

  std::sort(items.begin(), items.end(), [](const Entry &lhs, const Entry &rhs) {
    return lhs.col == rhs.col ? lhs.row < rhs.row : lhs.col < rhs.col;
  });

  std::vector<Entry> unique;
  unique.reserve(items.size());
  for (size_t i = 0; i < items.size();) {
    size_t j = i + 1;
    double sum = items[i].value;
    while (j < items.size() && items[j].row == items[i].row && items[j].col == items[i].col) {
      sum += items[j].value;
      ++j;
    }
    if (sum != 0.0) {
      unique.push_back({items[i].row, items[i].col, sum});
    }
    i = j;
  }

  a.n = static_cast<int>(rows);
  a.col_ptr.assign(static_cast<size_t>(a.n) + 1u, 0);
  a.row_idx.resize(unique.size());
  a.values.resize(unique.size());
  for (const Entry &e : unique) {
    a.col_ptr[static_cast<size_t>(e.col) + 1u]++;
  }
  for (int col = 0; col < a.n; ++col) {
    a.col_ptr[static_cast<size_t>(col) + 1u] += a.col_ptr[static_cast<size_t>(col)];
  }
  std::vector<int> next = a.col_ptr;
  for (const Entry &e : unique) {
    const int dst = next[static_cast<size_t>(e.col)]++;
    a.row_idx[static_cast<size_t>(dst)] = e.row;
    a.values[static_cast<size_t>(dst)] = e.value;
  }
  return true;
}

static void matvec(const Matrix &a, const std::vector<double> &x, std::vector<double> &y) {
  std::fill(y.begin(), y.end(), 0.0);
  for (int col = 0; col < a.n; ++col) {
    for (int p = a.col_ptr[static_cast<size_t>(col)];
         p < a.col_ptr[static_cast<size_t>(col) + 1u]; ++p) {
      y[static_cast<size_t>(a.row_idx[static_cast<size_t>(p)])] +=
        a.values[static_cast<size_t>(p)] * x[static_cast<size_t>(col)];
    }
  }
}

static void make_refactor_values(Matrix &a,
                                 const std::vector<double> &base_values,
                                 bench_refactor_value_mode mode,
                                 uint64_t generation, double amplitude) {
  for (int col = 0; col < a.n; ++col) {
    for (int p = a.col_ptr[static_cast<size_t>(col)];
         p < a.col_ptr[static_cast<size_t>(col) + 1u]; ++p) {
      a.values[static_cast<size_t>(p)] = bench_generated_refactor_value(
        mode, base_values[static_cast<size_t>(p)],
        static_cast<uint64_t>(a.row_idx[static_cast<size_t>(p)]),
        static_cast<uint64_t>(col), static_cast<uint64_t>(a.n), generation,
        amplitude);
    }
  }
}

static double residual(const Matrix &a,
                       const std::vector<double> &x,
                       const std::vector<double> &b,
                       double *relative_out) {
  std::vector<double> ax(static_cast<size_t>(a.n), 0.0);
  matvec(a, x, ax);
  double r2 = 0.0;
  double b2 = 0.0;
  for (int i = 0; i < a.n; ++i) {
    const double r = ax[static_cast<size_t>(i)] - b[static_cast<size_t>(i)];
    r2 += r * r;
    b2 += b[static_cast<size_t>(i)] * b[static_cast<size_t>(i)];
  }
  const double rnorm = std::sqrt(r2);
  if (relative_out != nullptr) {
    *relative_out = b2 > 0.0 ? rnorm / std::sqrt(b2) : rnorm;
  }
  return rnorm;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "Usage: %s <matrix.mtx> [threads] [repeat] [refactor-repeat] [factor-repeat] [unchanged|rank-preserving|entrywise|localized-entrywise] [amplitude]\n", argv[0]);
    return EXIT_FAILURE;
  }
  const int threads = argc > 2 ? std::atoi(argv[2]) : 16;
  const int repeat = argc > 3 ? std::atoi(argv[3]) : 5;
  const int refactor_repeat = argc > 4 ? std::atoi(argv[4]) : repeat;
  const int factor_repeat = argc > 5 ? std::atoi(argv[5]) : repeat;
  bench_refactor_value_mode refactor_value_mode =
    BENCH_REFACTOR_VALUES_UNCHANGED;
  if (argc > 6 &&
      !bench_parse_refactor_value_mode(argv[6], &refactor_value_mode)) {
    std::fprintf(stderr, "unknown refactor value mode: %s\n", argv[6]);
    return EXIT_FAILURE;
  }
  const double refactor_value_amplitude = argc > 7 ? std::atof(argv[7]) : 1.0e-3;
  if (threads <= 0 || repeat <= 0 || refactor_repeat < 0 ||
      factor_repeat < 0 || !std::isfinite(refactor_value_amplitude) ||
      refactor_value_amplitude < 0.0 || refactor_value_amplitude >= 1.0 ||
      (refactor_value_mode != BENCH_REFACTOR_VALUES_UNCHANGED &&
       refactor_value_amplitude <= 0.0)) {
    std::fprintf(stderr, "threads/repeat arguments must be positive\n");
    return EXIT_FAILURE;
  }

  Matrix a;
  if (!read_matrix_market(argv[1], a)) {
    return EXIT_FAILURE;
  }
  const std::vector<double> base_values = a.values;
  std::vector<double> x_true(static_cast<size_t>(a.n), 1.0);
  for (int i = 0; i < a.n; ++i) {
    x_true[static_cast<size_t>(i)] = 1.0 + static_cast<double>(i % 17) * 0.01;
  }
  std::vector<double> b(static_cast<size_t>(a.n), 0.0);
  std::vector<double> x(static_cast<size_t>(a.n), 0.0);
  matvec(a, x_true, b);

  ICktSo inst = nullptr;
  int *iparm = nullptr;
  const long long *oparm = nullptr;
  int ret = CKTSO_CreateSolver(&inst, &iparm, &oparm);
  if (ret < 0) {
    std::fprintf(stderr, "CKTSO_CreateSolver failed: %d\n", ret);
    return EXIT_FAILURE;
  }
  iparm[0] = 1;
  ret = CKTSO_Analyze(inst, false, a.n, a.col_ptr.data(), a.row_idx.data(),
                      a.values.data(), threads);
  if (ret < 0) {
    std::fprintf(stderr, "CKTSO_Analyze failed: %d\n", ret);
    CKTSO_DestroySolver(inst);
    return EXIT_FAILURE;
  }
  const long long analysis_us = oparm[0];

  ret = CKTSO_Factorize(inst, a.values.data(), true);
  if (ret < 0) {
    std::fprintf(stderr, "CKTSO initial Factorize failed: %d\n", ret);
    CKTSO_DestroySolver(inst);
    return EXIT_FAILURE;
  }
  const long long initial_factor_us = oparm[1];
  long long factor_flops = -1;
  long long solve_flops = -1;
  long long factor_mem = -1;
  long long solve_mem = -1;
  (void)CKTSO_Statistics(inst, &factor_flops, &solve_flops,
                         &factor_mem, &solve_mem, true, -1, false);

  long long factor_total = 0;
  long long refactor_total = 0;
  long long refactor_first = 0;
  long long solve_total = 0;
  for (int i = 0; i < factor_repeat; ++i) {
    ret = CKTSO_Factorize(inst, a.values.data(), true);
    if (ret < 0) break;
    factor_total += oparm[1];
  }
  for (int i = 0; i < refactor_repeat && ret >= 0; ++i) {
    if (refactor_value_mode != BENCH_REFACTOR_VALUES_UNCHANGED) {
      make_refactor_values(a, base_values, refactor_value_mode,
                           static_cast<uint64_t>(i) + 1u,
                           refactor_value_amplitude);
      matvec(a, x_true, b);
    }
    ret = CKTSO_Refactorize(inst, a.values.data());
    if (ret < 0) break;
    refactor_total += oparm[1];
    if (i == 0) refactor_first = oparm[1];
    /* Match a Newton loop and give adaptive solve paths the same one-solve
       opportunity after every sampled refactor as KLS. */
    ret = CKTSO_Solve(inst, b.data(), x.data(), false, true);
  }
  for (int i = 0; i < repeat && ret >= 0; ++i) {
    ret = CKTSO_Solve(inst, b.data(), x.data(), false, true);
    if (ret < 0) break;
    solve_total += oparm[2];
  }
  if (ret < 0) {
    std::fprintf(stderr, "CKTSO run failed: %d\n", ret);
    CKTSO_DestroySolver(inst);
    return EXIT_FAILURE;
  }

  const double factor_us_avg = factor_repeat > 0
    ? static_cast<double>(factor_total) / static_cast<double>(factor_repeat)
    : static_cast<double>(initial_factor_us);
  const double refactor_us_avg = refactor_repeat > 0
    ? static_cast<double>(refactor_total) / static_cast<double>(refactor_repeat)
    : 0.0;
  const double refactor_steady_us_avg = refactor_repeat > 1
    ? static_cast<double>(refactor_total - refactor_first) /
        static_cast<double>(refactor_repeat - 1)
    : static_cast<double>(refactor_first);
  const double solve_us_avg = static_cast<double>(solve_total) / static_cast<double>(repeat);
  const double spice_cycle_seconds =
    1.0e-6 * (static_cast<double>(analysis_us + initial_factor_us) +
              solve_us_avg + static_cast<double>(refactor_first) +
              solve_us_avg + 98.0 * (refactor_steady_us_avg + solve_us_avg));
  double relative_residual = 0.0;
  const double residual_l2 = residual(a, x, b, &relative_residual);

  std::printf("{\"matrix\":\"%s\",\"n\":%d,\"nnz\":%d,"
              "\"threads\":%d,\"repeat\":%d,\"factor_repeat\":%d,"
              "\"refactor_repeat\":%d,"
              "\"refactor_value_mode\":\"%s\","
              "\"refactor_value_amplitude\":%.9g,"
              "\"analysis_us\":%lld,\"initial_factor_us\":%lld,"
              "\"factor_us_avg\":%.9g,\"refactor_us_avg\":%.9g,"
              "\"solve_us_avg\":%.9g,"
              "\"analysis_seconds\":%.9g,\"initial_factor_seconds\":%.9g,"
              "\"factor_seconds_avg\":%.9g,\"refactor_seconds_avg\":%.9g,"
              "\"solve_seconds_avg\":%.9g,\"spice_cycle_seconds\":%.9g,"
              "\"refactor_first_seconds\":%.9g,\"refactor_steady_seconds_avg\":%.9g,"
              "\"residual_l2\":%.9g,\"relative_residual_l2\":%.9g,"
              "\"nnz_l\":%lld,\"nnz_u\":%lld,"
              "\"selected_ordering\":%lld,\"supernodes\":%lld,"
              "\"offdiagonal_pivots\":%lld,"
              "\"predicted_nnz\":%lld,\"predicted_flops\":%lld,"
              "\"factor_flops\":%lld,\"solve_flops\":%lld,"
              "\"factor_memory_access_bytes\":%lld,"
              "\"solve_memory_access_bytes\":%lld,"
              "\"memory_bytes\":%lld,\"memory_peak_bytes\":%lld}\n",
              argv[1], a.n, a.col_ptr[static_cast<size_t>(a.n)], threads,
              repeat, factor_repeat, refactor_repeat,
              bench_refactor_value_mode_name(refactor_value_mode),
              refactor_value_mode != BENCH_REFACTOR_VALUES_UNCHANGED
                ? refactor_value_amplitude : 0.0,
              analysis_us,
              initial_factor_us,
              factor_us_avg, refactor_us_avg, solve_us_avg,
              1.0e-6 * static_cast<double>(analysis_us),
              1.0e-6 * static_cast<double>(initial_factor_us),
              1.0e-6 * factor_us_avg, 1.0e-6 * refactor_us_avg,
              1.0e-6 * solve_us_avg, spice_cycle_seconds,
              1.0e-6 * static_cast<double>(refactor_first),
              1.0e-6 * refactor_steady_us_avg,
              residual_l2, relative_residual, oparm[5], oparm[6],
              oparm[8], oparm[7], oparm[4], oparm[16], oparm[17],
              factor_flops, solve_flops, factor_mem, solve_mem,
              oparm[12], oparm[13]);
  CKTSO_DestroySolver(inst);
  return EXIT_SUCCESS;
}
