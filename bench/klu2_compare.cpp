#include "klu2_defaults.hpp"
#include "klu2_analyze.hpp"
#include "klu2_factor.hpp"

template <typename Entry, typename Int>
static Int kls_klu2_scale_for_refactor(Int scale,
                                       Int n,
                                       Int Ap[],
                                       Int Ai[],
                                       double Ax[],
                                       double Rs[],
                                       decltype(nullptr),
                                       klu_common<Entry, Int> *common) {
  return klu_scale<Entry, Int>(
    scale, n, Ap, Ai, reinterpret_cast<Entry *>(Ax), Rs,
    static_cast<Int *>(nullptr), common);
}

#pragma push_macro("NULL")
#undef NULL
#define NULL nullptr
#undef KLU_scale
#define KLU_scale kls_klu2_scale_for_refactor
#undef KLU_refactor
#define KLU_refactor kls_klu2_refactor_impl
#include "klu2_refactor.hpp"
#undef KLU_refactor
#define KLU_refactor klu_refactor
#undef KLU_scale
#define KLU_scale klu_scale
#pragma pop_macro("NULL")

#include "klu2_solve.hpp"
#include "klu2_diagnostics.hpp"
#include "klu2_free_symbolic.hpp"
#include "klu2_free_numeric.hpp"

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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

static double now_seconds() {
  using clock = std::chrono::steady_clock;
  return std::chrono::duration<double>(clock::now().time_since_epoch()).count();
}

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
  const bool skew = std::strstr(line, "skew-symmetric") != nullptr;
  const bool pattern = std::strstr(line, "pattern") != nullptr;
  const bool complex = std::strstr(line, "complex") != nullptr;
  if (complex) {
    std::fprintf(stderr, "complex matrices are not supported by klu2_compare\n");
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
    std::fprintf(stderr, "%s must be a square matrix fitting 32-bit KLU2\n", path);
    std::fclose(fp);
    return false;
  }

  std::vector<Entry> items;
  items.reserve(static_cast<size_t>(entries) * ((symmetric || skew) ? 2u : 1u));
  for (long long k = 0; k < entries; ++k) {
    if (std::fgets(line, sizeof(line), fp) == nullptr) {
      std::fprintf(stderr, "%s ended before all entries were read\n", path);
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
    if ((symmetric || skew) && row != col) {
      items.push_back({col, row, skew ? -value : value});
    }
  }
  std::fclose(fp);

  std::sort(items.begin(), items.end(), [](const Entry &lhs, const Entry &rhs) {
    return lhs.col == rhs.col ? lhs.row < rhs.row : lhs.col < rhs.col;
  });

  std::vector<Entry> unique;
  unique.reserve(items.size());
  for (size_t i = 0; i < items.size();) {
    size_t j = i + 1u;
    double sum = items[i].value;
    while (j < items.size() && items[j].row == items[i].row &&
           items[j].col == items[i].col) {
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
    a.col_ptr[static_cast<size_t>(col) + 1u] +=
      a.col_ptr[static_cast<size_t>(col)];
  }
  std::vector<int> next = a.col_ptr;
  for (const Entry &e : unique) {
    const int dst = next[static_cast<size_t>(e.col)]++;
    a.row_idx[static_cast<size_t>(dst)] = e.row;
    a.values[static_cast<size_t>(dst)] = e.value;
  }
  return true;
}

static void matvec(const Matrix &a,
                   const std::vector<double> &x,
                   std::vector<double> &y) {
  std::fill(y.begin(), y.end(), 0.0);
  for (int col = 0; col < a.n; ++col) {
    for (int p = a.col_ptr[static_cast<size_t>(col)];
         p < a.col_ptr[static_cast<size_t>(col) + 1u]; ++p) {
      y[static_cast<size_t>(a.row_idx[static_cast<size_t>(p)])] +=
        a.values[static_cast<size_t>(p)] * x[static_cast<size_t>(col)];
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
    std::fprintf(stderr, "Usage: %s <matrix.mtx> [repeat] [refactor-repeat]\n", argv[0]);
    return EXIT_FAILURE;
  }
  const int repeat = argc > 2 ? std::atoi(argv[2]) : 5;
  const int refactor_repeat = argc > 3 ? std::atoi(argv[3]) : repeat;
  if (repeat <= 0 || refactor_repeat < 0) {
    std::fprintf(stderr, "repeat arguments must be valid\n");
    return EXIT_FAILURE;
  }

  Matrix a;
  if (!read_matrix_market(argv[1], a)) {
    return EXIT_FAILURE;
  }

  std::vector<double> x_true(static_cast<size_t>(a.n), 1.0);
  for (int i = 0; i < a.n; ++i) {
    x_true[static_cast<size_t>(i)] = 1.0 + static_cast<double>(i % 17) * 0.01;
  }
  std::vector<double> b(static_cast<size_t>(a.n), 0.0);
  std::vector<double> x(static_cast<size_t>(a.n), 0.0);
  matvec(a, x_true, b);

  klu_common<double, int> common;
  klu_defaults<double, int>(&common);

  const double analyze_start = now_seconds();
  klu_symbolic<double, int> *symbolic =
    klu_analyze<double, int>(a.n, a.col_ptr.data(), a.row_idx.data(), &common);
  const double analysis_seconds = now_seconds() - analyze_start;
  if (symbolic == nullptr || common.status < KLU_OK) {
    std::fprintf(stderr, "KLU2 analyze failed: %d\n", common.status);
    return EXIT_FAILURE;
  }

  const double initial_start = now_seconds();
  klu_numeric<double, int> *numeric =
    klu_factor<double, int>(a.col_ptr.data(), a.row_idx.data(), a.values.data(),
                            symbolic, &common);
  const double initial_factor_seconds = now_seconds() - initial_start;
  if (numeric == nullptr || common.status < KLU_OK) {
    std::fprintf(stderr, "KLU2 initial factor failed: %d\n", common.status);
    klu_free_symbolic<double, int>(&symbolic, &common);
    return EXIT_FAILURE;
  }

  double factor_total = 0.0;
  for (int i = 0; i < repeat; ++i) {
    klu_free_numeric<double, int>(&numeric, &common);
    const double start = now_seconds();
    numeric =
      klu_factor<double, int>(a.col_ptr.data(), a.row_idx.data(), a.values.data(),
                              symbolic, &common);
    factor_total += now_seconds() - start;
    if (numeric == nullptr || common.status < KLU_OK) {
      std::fprintf(stderr, "KLU2 factor failed: %d\n", common.status);
      klu_free_symbolic<double, int>(&symbolic, &common);
      return EXIT_FAILURE;
    }
  }

  double refactor_total = 0.0;
  for (int i = 0; i < refactor_repeat; ++i) {
    const double start = now_seconds();
    const int ok =
      kls_klu2_refactor_impl<double, int>(
        a.col_ptr.data(), a.row_idx.data(), a.values.data(), symbolic, numeric,
        &common);
    refactor_total += now_seconds() - start;
    if (!ok || common.status < KLU_OK) {
      std::fprintf(stderr, "KLU2 refactor failed: %d\n", common.status);
      klu_free_numeric<double, int>(&numeric, &common);
      klu_free_symbolic<double, int>(&symbolic, &common);
      return EXIT_FAILURE;
    }
  }

  double solve_total = 0.0;
  for (int i = 0; i < repeat; ++i) {
    x = b;
    const double start = now_seconds();
    const int ok =
      klu_solve<double, int>(symbolic, numeric, a.n, 1, x.data(), &common);
    solve_total += now_seconds() - start;
    if (!ok || common.status < KLU_OK) {
      std::fprintf(stderr, "KLU2 solve failed: %d\n", common.status);
      klu_free_numeric<double, int>(&numeric, &common);
      klu_free_symbolic<double, int>(&symbolic, &common);
      return EXIT_FAILURE;
    }
  }

  (void)klu_flops<double, int>(symbolic, numeric, &common);
  (void)klu_rcond<double, int>(symbolic, numeric, &common);
  const double factor_seconds_avg = factor_total / static_cast<double>(repeat);
  const double refactor_seconds_avg = refactor_repeat > 0
    ? refactor_total / static_cast<double>(refactor_repeat)
    : 0.0;
  const double solve_seconds_avg = solve_total / static_cast<double>(repeat);
  const double spice_cycle_seconds =
    analysis_seconds + initial_factor_seconds + solve_seconds_avg +
    99.0 * (refactor_seconds_avg + solve_seconds_avg);
  double relative_residual = 0.0;
  const double residual_l2 = residual(a, x, b, &relative_residual);

  std::printf("{\"matrix\":\"%s\",\"n\":%d,\"nnz\":%d,"
              "\"repeat\":%d,\"refactor_repeat\":%d,"
              "\"analysis_seconds\":%.9g,\"initial_factor_seconds\":%.9g,"
              "\"factor_seconds_avg\":%.9g,\"refactor_seconds_avg\":%.9g,"
              "\"solve_seconds_avg\":%.9g,\"spice_cycle_seconds\":%.9g,"
              "\"residual_l2\":%.9g,\"relative_residual_l2\":%.9g,"
              "\"nblocks\":%d,\"max_block\":%d,\"structural_rank\":%d,"
              "\"offdiag_pivots\":%d,\"nnz_l\":%d,\"nnz_u\":%d,"
              "\"factor_flops\":%.9g,\"rcond\":%.9g,"
              "\"memory_bytes\":%zu,\"memory_peak_bytes\":%zu}\n",
              argv[1], a.n, a.col_ptr[static_cast<size_t>(a.n)],
              repeat, refactor_repeat,
              analysis_seconds, initial_factor_seconds,
              factor_seconds_avg, refactor_seconds_avg, solve_seconds_avg,
              spice_cycle_seconds, residual_l2, relative_residual,
              symbolic->nblocks, symbolic->maxblock, symbolic->structural_rank,
              common.noffdiag, numeric->lnz, numeric->unz, common.flops,
              common.rcond, common.memusage, common.mempeak);

  klu_free_numeric<double, int>(&numeric, &common);
  klu_free_symbolic<double, int>(&symbolic, &common);
  return EXIT_SUCCESS;
}
