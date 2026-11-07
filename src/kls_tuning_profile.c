#define _GNU_SOURCE
#define _POSIX_C_SOURCE 200809L
#include "kls_tuning_profile.h"
#include "kls/kls.h"
#include "kls_tuning.inc"

#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#ifdef __linux__
#include <sched.h>
#endif
#if defined(__x86_64__) || defined(__i386__)
#include <cpuid.h>
#endif

#ifndef KLS_TUNING_BUILD_FINGERPRINT
#define KLS_TUNING_BUILD_FINGERPRINT "unknown-build"
#endif

typedef enum kls_tuning_kind { KLS_TUNE_UNSIGNED, KLS_TUNE_U64,
                               KLS_TUNE_DOUBLE } kls_tuning_kind;
typedef struct kls_tuning_field {
  const char *name;
  size_t offset;
  kls_tuning_kind kind;
  double minimum, maximum;
} kls_tuning_field;

#define TF_U(name, member, lo, hi) \
  {name, offsetof(kls_tuning_values, member), KLS_TUNE_UNSIGNED, lo, hi}
#define TF_Q(name, member, lo, hi) \
  {name, offsetof(kls_tuning_values, member), KLS_TUNE_U64, lo, hi}
#define TF_D(name, member, lo, hi) \
  {name, offsetof(kls_tuning_values, member), KLS_TUNE_DOUBLE, lo, hi}

static const kls_tuning_field kls_tuning_fields[] = {
  TF_U("row_batch_min_rows", row_batch_min_rows, 2, 64),
  TF_U("row_batch_max_rows", row_batch_max_rows, 2, 128),
  TF_U("row_small_sort_max", row_small_sort_max, 8, 256),
  TF_D("row_dense_min_work", row_dense_min_work, 64, 1048576),
  TF_D("row_compact_panel_min_work", row_compact_panel_min_work, 256, 16777216),
  TF_D("row_compact_panel_min_work_per_entry", row_compact_panel_min_work_per_entry, 1, 64),
  TF_D("row_native_panel_min_work", row_native_panel_min_work, 256, 16777216),
  TF_U("row_blocked_trailing_min_rows", row_blocked_trailing_min_rows, 2, 128),
  TF_U("row_blocked_trailing_min_cols", row_blocked_trailing_min_cols, 2, 256),
  TF_D("row_batch_snode_min_work", row_batch_snode_min_work, 256, 16777216),
  TF_D("row_batch_snode_min_work_per_entry", row_batch_snode_min_work_per_entry, 1, 64),
  TF_U("dense_help_tail_slice", dense_help_tail_slice, 16, 4096),
  TF_U("btf_scalar_min_rows", btf_scalar_min_rows, 2, 4096),
  TF_U("btf_scalar_max_rows", btf_scalar_max_rows, 16, 65536),
  TF_Q("solve_dense_tail_min_nnz", solve_dense_tail_min_nnz, 4096, 1.0e10),
  TF_D("solve_dense_tail_min_fraction", solve_dense_tail_min_fraction, 0.25, 0.95),
  TF_U("solve_trapezoid_slices", solve_trapezoid_slices, 2, 32),
  TF_Q("solve_parallel_rect_min_nnz", solve_parallel_rect_min_nnz, 4096, 1.0e10),
  TF_Q("solve_min_entries_per_sync", solve_min_entries_per_sync, 4096, 1.0e10),
  TF_U("snode_min_batch", snode_min_batch, 2, KLS_SNODE_MAX_BATCH),
  TF_U("snode_min_batch_work", snode_min_batch_work, 16, 65536),
  TF_D("egraph_min_flops_per_thread", egraph_min_flops_per_thread, 1000, 1.0e9),
  TF_Q("egraph_min_size", egraph_min_size, 128, 100000000),
  TF_U("egraph_pool_spin_iters", egraph_pool_spin_iters, 0, 2000000),
  TF_U("egraph_worker_spin_iters", egraph_worker_spin_iters, 0, 2000000),
  TF_Q("snb_coop_min_doubles", snb_coop_min_doubles, 1024, 100000000)
};

static size_t kls_tuning_field_count(void) {
  return sizeof(kls_tuning_fields) / sizeof(kls_tuning_fields[0]);
}

typedef struct kls_cost_field { const char *name; size_t offset; } kls_cost_field;
#define CF(name, member) {name, offsetof(kls_snb_cost_model, member)}
static const kls_cost_field kls_cost_fields[] = {
  CF("snb_narrow_edge", narrow[0]), CF("snb_narrow_trsm", narrow[1]),
  CF("snb_narrow_scatter", narrow[2]), CF("snb_dense_edge", dense[0]),
  CF("snb_dense_trsm", dense[1]), CF("snb_dense_vectors", dense[2]),
  CF("snb_dense_loads", dense[3]), CF("snb_dense_tails", dense[4]),
  CF("snb_dense_scatter", dense[5]), CF("snb_panel_fixed", panel[0]),
  CF("snb_panel_divisions", panel[1]), CF("snb_panel_updates", panel[2]),
  CF("snb_panel_entries", panel[3])
};
#undef CF
#define KLS_COST_FIELD_COUNT (sizeof(kls_cost_fields)/sizeof(kls_cost_fields[0]))

static int kls_cost_model_valid(const kls_snb_cost_model *m) {
  if (m == NULL || !m->enabled) return 1;
  double total = 0.0;
  for (size_t i=0;i<KLS_COST_FIELD_COUNT;++i) {
    const double value=*(const double *)((const unsigned char *)m+kls_cost_fields[i].offset);
    if (!isfinite(value) || value < 0.0) return 0;
    total += value;
  }
  return isfinite(total) && m->enabled == 1 && m->narrow[0]+m->narrow[1]+m->narrow[2] > 0.0 &&
    m->dense[0]+m->dense[1]+m->dense[2]+m->dense[3]+m->dense[4]+m->dense[5] > 0.0 &&
    m->panel[0]+m->panel[1]+m->panel[2]+m->panel[3] > 0.0;
}

static int kls_set_cost_field(kls_snb_cost_model *m, const char *key,
                              const char *value, unsigned *seen) {
  for (size_t i=0;i<KLS_COST_FIELD_COUNT;++i) {
    if (strcmp(key,kls_cost_fields[i].name)) continue;
    char *end=NULL; errno=0;
    const double parsed=strtod(value,&end);
    if ((*seen & (1u<<i)) || errno || end==value || *end ||
        !isfinite(parsed) || parsed<0.0) return 0;
    *(double *)((unsigned char *)m+kls_cost_fields[i].offset)=parsed;
    *seen |= 1u<<i;
    return 1;
  }
  return -1;
}

static int kls_tuning_changed_fields(const kls_tuning_values *values) {
  kls_tuning_values defaults;
  kls_tuning_defaults(&defaults);
  int count = 0;
  for (size_t i = 0; i < kls_tuning_field_count(); ++i) {
    const kls_tuning_field *f = &kls_tuning_fields[i];
    const unsigned char *a = (const unsigned char *)values + f->offset;
    const unsigned char *b = (const unsigned char *)&defaults + f->offset;
    if (f->kind == KLS_TUNE_DOUBLE) count += *(const double *)a != *(const double *)b;
    else if (f->kind == KLS_TUNE_UNSIGNED) count += *(const unsigned *)a != *(const unsigned *)b;
    else count += *(const uint64_t *)a != *(const uint64_t *)b;
  }
  return count;
}

void kls_tuning_defaults(kls_tuning_values *v) {
  if (v == NULL) return;
  *v = (kls_tuning_values){
    KLS_ROW_REFACTOR_BATCH_MIN_ROWS, KLS_ROW_REFACTOR_BATCH_MAX_ROWS,
    KLS_ROW_REFACTOR_SMALL_SORT_MAX, KLS_ROW_REFACTOR_DENSE_MIN_WORK,
    KLS_ROW_REFACTOR_COMPACT_PANEL_MIN_WORK,
    KLS_ROW_REFACTOR_COMPACT_PANEL_MIN_WORK_PER_ENTRY,
    KLS_ROW_REFACTOR_NATIVE_PANEL_AUTO_MIN_WORK,
    KLS_ROW_REFACTOR_BLOCKED_TRAILING_MIN_ROWS,
    KLS_ROW_REFACTOR_BLOCKED_TRAILING_MIN_COLS,
    KLS_ROW_REFACTOR_BATCH_SUPERNODE_MIN_WORK,
    KLS_ROW_REFACTOR_BATCH_SUPERNODE_MIN_WORK_PER_ENTRY,
    KLS_DENSE_HELP_TAIL_SLICE, KLS_BTF_SCALAR_RUN_EXEC_MIN_ROWS,
    KLS_BTF_SCALAR_RUN_EXEC_MAX_ROWS,
    KLS_ROW_SOLVE_DENSE_TAIL_MIN_NNZ, KLS_ROW_SOLVE_DENSE_TAIL_MIN_FRACTION,
    KLS_ROW_SOLVE_TRAPEZOID_SLICES, KLS_ROW_SOLVE_PARALLEL_RECT_MIN_NNZ,
    KLS_ROW_SOLVE_PARALLEL_MIN_ENTRIES_PER_SYNC,
    KLS_SNODE_DEFAULT_MIN_BATCH, KLS_SNODE_DEFAULT_MIN_BATCH_WORK,
    KLS_EGRAPH_REFACTOR_MIN_FLOPS_PER_THREAD, KLS_EGRAPH_REFACTOR_MIN_SIZE,
    KLS_EGRAPH_POOL_SPIN_ITERS, KLS_EGRAPH_WORKER_SPIN_ITERS,
    KLS_SNB_COOP_MIN_DOUBLES
  };
}

const char *kls_tuning_build_fingerprint(void) {
  return KLS_TUNING_BUILD_FINGERPRINT;
}

static void kls_cpu_fingerprint(char *out, size_t cap) {
#if defined(__x86_64__) || defined(__i386__)
  unsigned eax = 0, ebx = 0, ecx = 0, edx = 0;
  char vendor[13] = {0};
  __get_cpuid(0, &eax, &ebx, &ecx, &edx);
  memcpy(vendor, &ebx, 4); memcpy(vendor + 4, &edx, 4);
  memcpy(vendor + 8, &ecx, 4);
  unsigned family = 0, model = 0, stepping = 0, features = 0;
  unsigned ext_features = 0;
  if (__get_cpuid(1, &eax, &ebx, &ecx, &edx)) {
    stepping = eax & 15u;
    family = (eax >> 8u) & 15u;
    model = (eax >> 4u) & 15u;
    if (family == 15u) family += (eax >> 20u) & 255u;
    if (family == 6u || family == 15u) model += ((eax >> 16u) & 15u) << 4u;
    features = ecx;
  }
  if (__get_cpuid_max(0, NULL) >= 7u)
    __cpuid_count(7, 0, eax, ext_features, ecx, edx);
  unsigned usable_isa = 0u;
#if defined(__GNUC__) || defined(__clang__)
  __builtin_cpu_init();
  if (__builtin_cpu_supports("avx2")) usable_isa |= 1u;
  if (__builtin_cpu_supports("avx512f")) usable_isa |= 2u;
  if (__builtin_cpu_supports("fma")) usable_isa |= 4u;
#endif
  long l2 = sysconf(_SC_LEVEL2_CACHE_SIZE);
  long l3 = sysconf(_SC_LEVEL3_CACHE_SIZE);
  snprintf(out, cap, "%s-%u-%u-%u-%08x-%08x-u%x-l2%ld-l3%ld", vendor, family,
           model, stepping, features, ext_features, usable_isa,
           l2 > 0 ? l2 : 0, l3 > 0 ? l3 : 0);
#else
  snprintf(out, cap, "unsupported");
#endif
}

static uint64_t kls_hash_update(uint64_t h, const char *s, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    h ^= (unsigned char)s[i];
    h *= UINT64_C(1099511628211);
  }
  return h;
}

static int kls_set_field(kls_tuning_values *v, const char *key,
                         const char *text, unsigned char *seen) {
  for (size_t i = 0; i < kls_tuning_field_count(); ++i) {
    const kls_tuning_field *f = &kls_tuning_fields[i];
    if (strcmp(key, f->name) != 0) continue;
    if (seen[i]) return 0;
    char *end = NULL;
    errno = 0;
    double parsed = strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0' || !isfinite(parsed) ||
        parsed < f->minimum || parsed > f->maximum) return 0;
    unsigned char *base = (unsigned char *)v + f->offset;
    if (f->kind == KLS_TUNE_DOUBLE) *(double *)base = parsed;
    else if (f->kind == KLS_TUNE_UNSIGNED) {
      if (floor(parsed) != parsed) return 0;
      *(unsigned *)base = (unsigned)parsed;
    } else {
      if (floor(parsed) != parsed) return 0;
      *(uint64_t *)base = (uint64_t)parsed;
    }
    seen[i] = 1;
    return 1;
  }
  return -1;
}

static int kls_tuning_consistent(const kls_tuning_values *v) {
  return v->row_batch_min_rows <= v->row_batch_max_rows &&
         v->btf_scalar_min_rows <= v->btf_scalar_max_rows &&
         v->solve_trapezoid_slices >= 2u;
}

static const char *const kls_policy_names[KLS_POLICY_COST_COUNT]={
  "policy_scan","policy_pair","policy_sort","policy_flop","policy_entry",
  "policy_solve_entry","policy_copy_entry","policy_dispatch",
  "policy_factor_flop","policy_factor_entry","policy_factor_dispatch",
  "policy_first_refactor_entry"};

uint64_t kls_tuning_policy_cache_bytes(void) {
#ifdef __linux__
  cpu_set_t pending; CPU_ZERO(&pending);
  if (sched_getaffinity(0, sizeof pending, &pending)) return 0;
  uint64_t minimum = UINT64_MAX;
  for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) if (CPU_ISSET(cpu, &pending)) {
    char path[160], line[512];
    snprintf(path, sizeof path, "/sys/devices/system/cpu/cpu%d/cache/index3/size", cpu);
    FILE *f = fopen(path, "r"); if (!f) return 0;
    int ok = fgets(line, sizeof line, f) != NULL; fclose(f); if (!ok) return 0;
    char *end = NULL; errno = 0;
    unsigned long long size = strtoull(line, &end, 10);
    if (errno || !size || (*end != 'K' && *end != 'M')) return 0;
    const uint64_t scale = *end == 'K' ? 1024u : 1024u * 1024u;
    if (size > UINT64_MAX / scale) return 0;
    const uint64_t bytes = size * scale;
    if (bytes < minimum) minimum = bytes;
    CPU_CLR(cpu, &pending);
    /* Members of one shared cache have the same capacity. Read each cache
       once per call, while still respecting a changed caller affinity. No
       learned policy or affinity-sensitive result is cached globally. */
    snprintf(path, sizeof path,
      "/sys/devices/system/cpu/cpu%d/cache/index3/shared_cpu_list", cpu);
    f = fopen(path, "r"); if (!f) continue;
    ok = fgets(line, sizeof line, f) != NULL; fclose(f);
    if (!ok || !strchr(line, '\n')) continue;
    cpu_set_t shared; CPU_ZERO(&shared);
    const char *at = line;
    while (*at && *at != '\n') {
      errno = 0; long first = strtol(at, &end, 10);
      if (errno || end == at || first < 0 || first >= CPU_SETSIZE) break;
      at = end; long last = first;
      if (*at == '-') {
        ++at; errno = 0; last = strtol(at, &end, 10);
        if (errno || end == at || last < first || last >= CPU_SETSIZE) break;
        at = end;
      }
      if (*at != ',' && *at != '\n' && *at) break;
      for (long member = first; member <= last; ++member) CPU_SET((int)member, &shared);
      if (*at == ',') ++at;
    }
    if (*at == '\n' || !*at)
      for (int member = 0; member < CPU_SETSIZE; ++member)
        if (CPU_ISSET(member, &shared)) CPU_CLR(member, &pending);
  }
  return minimum == UINT64_MAX ? 0 : minimum;
#else
  return 0;
#endif
}

static int kls_policy_model_valid(const kls_policy_cost_model *p) {
  if(!p || !p->enabled)return 1;
  if(p->enabled!=1 || !p->cache_bytes ||
      p->cache_bytes!=kls_tuning_policy_cache_bytes())return 0;
  for(int i=0;i<KLS_POLICY_COST_COUNT;++i)
    if(!isfinite(p->seconds[i]) || p->seconds[i]<0)return 0;
  return p->seconds[KLS_POLICY_SCAN]+p->seconds[KLS_POLICY_PAIR]>0 &&
    p->seconds[KLS_POLICY_FLOP]+p->seconds[KLS_POLICY_ENTRY]>0 &&
    p->seconds[KLS_POLICY_FACTOR_FLOP]+p->seconds[KLS_POLICY_FACTOR_ENTRY]>0 &&
    p->seconds[KLS_POLICY_SOLVE_ENTRY]>0 && p->seconds[KLS_POLICY_COPY_ENTRY]>0 &&
    p->seconds[KLS_POLICY_SORT]>0;
}

int kls_tuning_load_host_profile(const char *path, int threads,
                                 kls_tuning_values *values,
                                 kls_tuning_metadata *metadata) {
  if (path == NULL || values == NULL || metadata == NULL || threads <= 0)
    return KLS_ERR_INVALID_ARGUMENT;
  FILE *file = fopen(path, "r");
  if (file == NULL) return KLS_ERR_TUNING_PROFILE;
  kls_tuning_values candidate;
  kls_snb_cost_model cost = {0};
  unsigned cost_seen = 0;
  int cost_version_seen = 0;
  kls_policy_cost_model policy={0};
  unsigned policy_seen=0;
  int policy_version_seen=0,policy_cache_seen=0;
  kls_tuning_defaults(&candidate);
  unsigned char seen[sizeof(kls_tuning_fields) / sizeof(kls_tuning_fields[0])] = {0};
  int format_seen = 0, build_seen = 0, cpu_seen = 0, threads_seen = 0;
  int seed_seen = 0;
  int checksum_seen = 0;
  uint64_t hash = UINT64_C(1469598103934665603), expected_hash = 0;
  char host_cpu[192]; kls_cpu_fingerprint(host_cpu, sizeof(host_cpu));
  char line[512];
  while (fgets(line, sizeof(line), file) != NULL) {
    if (checksum_seen) goto fail;
    size_t len = strlen(line);
    if (len == sizeof(line) - 1u && line[len - 1u] != '\n') goto fail;
    if (line[0] == '#' || line[0] == '\n') continue;
    if (strncmp(line, "checksum=", 9) == 0) {
      char tail = 0;
      if (checksum_seen || sscanf(line + 9, "%" SCNx64 " %c", &expected_hash, &tail) != 1)
        goto fail;
      checksum_seen = 1; continue;
    }
    hash = kls_hash_update(hash, line, len);
    if (len > 0 && line[len - 1u] == '\n') line[--len] = '\0';
    char *equals = strchr(line, '=');
    if (equals == NULL || equals == line) goto fail;
    *equals = '\0'; const char *value = equals + 1;
    if (strcmp(line, "format") == 0) {
      if (format_seen++ || strcmp(value, "kls-tuning-v1") != 0) goto fail;
    } else if (strcmp(line, "build") == 0) {
      if (build_seen++ || strcmp(value, kls_tuning_build_fingerprint()) != 0) goto fail;
    } else if (strcmp(line, "cpu") == 0) {
      if (cpu_seen++ || strcmp(value, host_cpu) != 0) goto fail;
    } else if (strcmp(line, "threads") == 0) {
      char *end = NULL; long n = strtol(value, &end, 10);
      if (threads_seen++ || *value == '\0' || *end != '\0' || n != threads) goto fail;
    } else if (strcmp(line, "generator_seed") == 0) {
      char *end = NULL; (void)strtoull(value, &end, 10);
      if (seed_seen++ || *value == '\0' || *end != '\0') goto fail;
    } else if (strcmp(line, "snb_cost_version") == 0) {
      if (cost_version_seen++ || strcmp(value,"1") || threads != 1) goto fail;
      cost.enabled=1;
    } else if(strcmp(line,"policy_cost_version")==0) {
      if(policy_version_seen++ || strcmp(value,"1"))goto fail;
      policy.enabled=1;
    } else if(strcmp(line,"policy_cache_bytes")==0) {
      char *end=NULL;errno=0;
      unsigned long long bytes=strtoull(value,&end,10);
      if(policy_cache_seen++ || errno || !*value || *value=='-' || *end || !bytes)goto fail;
      policy.cache_bytes=bytes;
    } else {
      int result = kls_set_field(&candidate, line, value, seen);
      if (result < 0) result = kls_set_cost_field(&cost, line, value, &cost_seen);
      if(result<0)for(int i=0;i<KLS_POLICY_COST_COUNT;++i)
        if(!strcmp(line,kls_policy_names[i])) {
          char *end=NULL;errno=0;double v=strtod(value,&end);
          if((policy_seen&(1u<<i)) || errno || end==value || *end ||
             !isfinite(v) || v<0)goto fail;
          policy.seconds[i]=v;policy_seen|=1u<<i;result=1;break;
        }
      if (result <= 0) goto fail;
    }
  }
  fclose(file);
  for (size_t i = 0; i < kls_tuning_field_count(); ++i)
    if (!seen[i]) return KLS_ERR_TUNING_PROFILE;
  if (!format_seen || !build_seen || !cpu_seen || !threads_seen || !seed_seen ||
      !checksum_seen || hash != expected_hash || !kls_tuning_consistent(&candidate) ||
      (cost_version_seen ? cost_seen != (1u<<KLS_COST_FIELD_COUNT)-1u : cost_seen != 0) ||
      (policy_version_seen ? (!policy_cache_seen || policy_seen!=(1u<<KLS_POLICY_COST_COUNT)-1u)
                           : (policy_cache_seen || policy_seen)) ||
      !kls_policy_model_valid(&policy) ||
      !kls_cost_model_valid(&cost))
    return KLS_ERR_TUNING_PROFILE;
  *values = candidate;
  *metadata = (kls_tuning_metadata){.profile_id=hash,
    .field_count=kls_tuning_changed_fields(&candidate)+
      (cost.enabled?(int)KLS_COST_FIELD_COUNT:0)+
      (policy.enabled?KLS_POLICY_COST_COUNT:0),
    .threads=threads,.snb_cost=cost,.policy_cost=policy,
    .kernel_overrides=kls_tuning_changed_fields(&candidate)+(cost.enabled!=0)};
  return KLS_OK;
fail:
  fclose(file);
  return KLS_ERR_TUNING_PROFILE;
}

static int kls_emit(FILE *f, uint64_t *hash, const char *line) {
  size_t n = strlen(line);
  if (fputs(line, f) == EOF) return 0;
  *hash = kls_hash_update(*hash, line, n);
  return 1;
}

int kls_tuning_write_host_profile(const char *path, int threads,
                                  const kls_tuning_values *values,
                                  unsigned long long seed, int force,
                                  kls_tuning_metadata *metadata) {
  return kls_tuning_write_host_profile_model(path,threads,values,NULL,seed,force,metadata);
}

static int kls_tuning_write_profile(const char *path, int threads,
                                  const kls_tuning_values *values,
                                  const kls_snb_cost_model *model,
                                  const kls_policy_cost_model *policy,
                                  unsigned long long seed, int force,
                                  kls_tuning_metadata *metadata) {
  if (path == NULL || values == NULL || threads <= 0 ||
      !kls_tuning_consistent(values) || !kls_cost_model_valid(model) ||
      !kls_policy_model_valid(policy) ||
      (model != NULL && model->enabled && threads != 1)) return KLS_ERR_INVALID_ARGUMENT;
  if (!force && access(path, F_OK) == 0) return KLS_ERR_TUNING_PROFILE;
  char temp[4096];
  if (snprintf(temp, sizeof(temp), "%s.tmp.%ld", path, (long)getpid()) >= (int)sizeof(temp))
    return KLS_ERR_INVALID_ARGUMENT;
  int fd = open(temp, O_WRONLY | O_CREAT | O_EXCL, 0600);
  if (fd < 0) return KLS_ERR_TUNING_PROFILE;
  FILE *f = fdopen(fd, "w");
  if (f == NULL) { close(fd); unlink(temp); return KLS_ERR_TUNING_PROFILE; }
  uint64_t hash = UINT64_C(1469598103934665603);
  char line[512], cpu[192]; kls_cpu_fingerprint(cpu, sizeof(cpu));
#define EMIT(...) do { snprintf(line, sizeof(line), __VA_ARGS__); \
  if (!kls_emit(f, &hash, line)) goto write_fail; } while (0)
  EMIT("format=kls-tuning-v1\n");
  EMIT("build=%s\n", kls_tuning_build_fingerprint());
  EMIT("cpu=%s\n", cpu);
  EMIT("threads=%d\n", threads);
  EMIT("generator_seed=%llu\n", seed);
  for (size_t i = 0; i < kls_tuning_field_count(); ++i) {
    const kls_tuning_field *field = &kls_tuning_fields[i];
    const unsigned char *base = (const unsigned char *)values + field->offset;
    if (field->kind == KLS_TUNE_DOUBLE) EMIT("%s=%.17g\n", field->name, *(const double *)base);
    else if (field->kind == KLS_TUNE_UNSIGNED) EMIT("%s=%u\n", field->name, *(const unsigned *)base);
    else EMIT("%s=%" PRIu64 "\n", field->name, *(const uint64_t *)base);
  }
  if (model != NULL && model->enabled) {
    EMIT("snb_cost_version=1\n");
    for (size_t i=0;i<KLS_COST_FIELD_COUNT;++i) {
      const double value=*(const double *)((const unsigned char *)model+kls_cost_fields[i].offset);
      EMIT("%s=%.17g\n",kls_cost_fields[i].name,value);
    }
  }
  if(policy && policy->enabled) {
    EMIT("policy_cost_version=1\n");
    EMIT("policy_cache_bytes=%" PRIu64 "\n",policy->cache_bytes);
    for(int i=0;i<KLS_POLICY_COST_COUNT;++i)
      EMIT("%s=%.17g\n",kls_policy_names[i],policy->seconds[i]);
  }
  if (fprintf(f, "checksum=%016" PRIx64 "\n", hash) < 0 ||
      fflush(f) != 0 || fsync(fileno(f)) != 0 || fclose(f) != 0) {
    f = NULL; goto write_fail;
  }
  f = NULL;
  if (force) {
    if (rename(temp, path) != 0) goto write_fail;
  } else {
    if (link(temp, path) != 0 || unlink(temp) != 0) goto write_fail;
  }
  if (metadata != NULL)
    *metadata = (kls_tuning_metadata){.profile_id=hash,
      .field_count=kls_tuning_changed_fields(values)+
        (model&&model->enabled?(int)KLS_COST_FIELD_COUNT:0)+
        (policy&&policy->enabled?KLS_POLICY_COST_COUNT:0),
      .threads=threads,.snb_cost=model?*model:(kls_snb_cost_model){0},
      .policy_cost=policy?*policy:(kls_policy_cost_model){0},
      .kernel_overrides=kls_tuning_changed_fields(values)+(model&&model->enabled)};
  return KLS_OK;
write_fail:
  if (f != NULL) fclose(f);
  unlink(temp);
  return KLS_ERR_TUNING_PROFILE;
#undef EMIT
}

int kls_tuning_write_host_profile_model(const char *path,int threads,
    const kls_tuning_values *values,const kls_snb_cost_model *model,
    unsigned long long seed,int force,kls_tuning_metadata *metadata) {
  return kls_tuning_write_profile(path,threads,values,model,NULL,seed,force,metadata);
}

int kls_tuning_write_policy_profile(const char *path,int threads,
    const kls_policy_cost_model *model,unsigned long long seed,int force) {
  if(!model || !model->enabled)return KLS_ERR_INVALID_ARGUMENT;
  kls_tuning_values values;kls_tuning_defaults(&values);
  return kls_tuning_write_profile(path,threads,&values,NULL,model,seed,force,NULL);
}
