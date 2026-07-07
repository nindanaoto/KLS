#include "subtree_lu.h"

#include <algorithm>
#include <cmath>
#include <chrono>
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
    std::fprintf(stderr, "Usage: %s <matrix.mtx> [threads] [repeat] [refactor-repeat]\n", argv[0]);
    return EXIT_FAILURE;
  }
  const int threads = argc > 2 ? std::atoi(argv[2]) : 16;
  const int repeat = argc > 3 ? std::atoi(argv[3]) : 5;
  const int refactor_repeat = argc > 4 ? std::atoi(argv[4]) : repeat;
  if (threads <= 0 || repeat <= 0 || refactor_repeat < 0) {
    std::fprintf(stderr, "threads/repeat arguments must be positive\n");
    return EXIT_FAILURE;
  }

  Matrix a;
  if (!read_matrix_market(argv[1], a)) {
    return EXIT_FAILURE;
  }
  std::vector<double> x_true(static_cast<size_t>(a.n), 1.0);
  std::vector<double> b(static_cast<size_t>(a.n), 0.0);
  std::vector<double> x(static_cast<size_t>(a.n), 0.0);
  matvec(a, x_true, b);

  subtree_lu::SubtreeLU<int, double> solver;
  long long *parm = solver.parm;
  int ret = solver.analyze(a.n, a.col_ptr.data(), a.row_idx.data(),
                           a.values.data(), threads, nullptr, nullptr);
  if (ret != subtree_lu::E_OK) {
    std::fprintf(stderr, "SubtreeLU analyze failed: %d\n", ret);
    return EXIT_FAILURE;
  }
  const long long analysis_us = parm[subtree_lu::O_ANALYZE_TIME];

  ret = solver.factorize(a.values.data());
  if (ret != subtree_lu::E_OK) {
    std::fprintf(stderr, "SubtreeLU initial factorize failed: %d\n", ret);
    return EXIT_FAILURE;
  }
  const long long initial_factor_us = parm[subtree_lu::O_FACTORIZE_TIME];

  long long factor_total = 0;
  long long refactor_total = 0;
  long long solve_total = 0;
  long long refactor_wall_total = 0;
  long long solve_wall_total = 0;
  for (int i = 0; i < repeat; ++i) {
    ret = solver.factorize(a.values.data());
    if (ret != subtree_lu::E_OK) break;
    factor_total += parm[subtree_lu::O_FACTORIZE_TIME];
  }
  for (int i = 0; i < refactor_repeat && ret == subtree_lu::E_OK; ++i) {
    const auto w0 = std::chrono::steady_clock::now();
    ret = solver.refactorize(a.values.data());
    const auto w1 = std::chrono::steady_clock::now();
    if (ret != subtree_lu::E_OK) break;
    refactor_total += parm[subtree_lu::O_FACTORIZE_TIME];
    refactor_wall_total +=
      std::chrono::duration_cast<std::chrono::microseconds>(w1 - w0).count();
  }
  for (int i = 0; i < repeat && ret == subtree_lu::E_OK; ++i) {
    const auto w0 = std::chrono::steady_clock::now();
    ret = solver.solve(b.data(), x.data());
    const auto w1 = std::chrono::steady_clock::now();
    if (ret != subtree_lu::E_OK) break;
    solve_total += parm[subtree_lu::O_SOLVE_TIME];
    solve_wall_total +=
      std::chrono::duration_cast<std::chrono::microseconds>(w1 - w0).count();
  }
  if (ret != subtree_lu::E_OK) {
    std::fprintf(stderr, "SubtreeLU run failed: %d\n", ret);
    return EXIT_FAILURE;
  }

  const double factor_us_avg = static_cast<double>(factor_total) / static_cast<double>(repeat);
  const double refactor_us_avg = refactor_repeat > 0
    ? static_cast<double>(refactor_total) / static_cast<double>(refactor_repeat)
    : 0.0;
  const double solve_us_avg = static_cast<double>(solve_total) / static_cast<double>(repeat);
  const double spice_cycle_seconds =
    1.0e-6 * (static_cast<double>(analysis_us + initial_factor_us) +
              solve_us_avg + 99.0 * (refactor_us_avg + solve_us_avg));
  std::fprintf(stderr,
               "subtreelu wall: refactor %.1fus (reported %.1fus) "
               "solve %.1fus (reported %.1fus)\n",
               refactor_repeat > 0
                 ? (double)refactor_wall_total / refactor_repeat : 0.0,
               refactor_us_avg,
               (double)solve_wall_total / repeat, solve_us_avg);
  double relative_residual = 0.0;
  const double residual_l2 = residual(a, x, b, &relative_residual);

  std::printf("{\"matrix\":\"%s\",\"n\":%d,\"nnz\":%d,"
              "\"threads\":%d,\"repeat\":%d,\"refactor_repeat\":%d,"
              "\"analysis_us\":%lld,\"initial_factor_us\":%lld,"
              "\"factor_us_avg\":%.9g,\"refactor_us_avg\":%.9g,"
              "\"solve_us_avg\":%.9g,"
              "\"analysis_seconds\":%.9g,\"initial_factor_seconds\":%.9g,"
              "\"factor_seconds_avg\":%.9g,\"refactor_seconds_avg\":%.9g,"
              "\"solve_seconds_avg\":%.9g,\"spice_cycle_seconds\":%.9g,"
              "\"residual_l2\":%.9g,\"relative_residual_l2\":%.9g,"
              "\"nnz_l\":%lld,\"nnz_u\":%lld,"
              "\"nsupernodes\":%lld,\"factorize_flops\":%lld}\n",
              argv[1], a.n, a.col_ptr[static_cast<size_t>(a.n)], threads,
              repeat, refactor_repeat, analysis_us, initial_factor_us,
              factor_us_avg, refactor_us_avg, solve_us_avg,
              1.0e-6 * static_cast<double>(analysis_us),
              1.0e-6 * static_cast<double>(initial_factor_us),
              1.0e-6 * factor_us_avg, 1.0e-6 * refactor_us_avg,
              1.0e-6 * solve_us_avg, spice_cycle_seconds,
              residual_l2, relative_residual,
              parm[subtree_lu::O_LNNZ], parm[subtree_lu::O_UNNZ],
              parm[subtree_lu::O_NSUPERNODES],
              parm[subtree_lu::O_FACTORIZE_FLOPS]);
  return EXIT_SUCCESS;
}
