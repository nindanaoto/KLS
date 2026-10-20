#define _POSIX_C_SOURCE 200809L
#include "kls/kls.h"
#include "kls_componentwise.h"
#include <stdio.h>
#include <time.h>


#include "trilinos_klu_decl.h"
#include <stdatomic.h>
#include <stdint.h>
enum { OUTSIDE, ANALYZE, FACTOR, REFACTOR, SOLVE, COMBINED, CERTIFICATE, CONTEXTS };
static _Thread_local unsigned phase_context;
static _Atomic uint64_t klu_ns[CONTEXTS][4], klu_calls[CONTEXTS][4];
static uint64_t cert_parent_ns[CONTEXTS];
static double combined_solve_stats_seconds;
static uint64_t tick(void) {
  struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);
  return (uint64_t)t.tv_sec*UINT64_C(1000000000)+(uint64_t)t.tv_nsec;
}
static void klu_record(unsigned ctx,unsigned op,uint64_t start) {
  uint64_t elapsed=tick()-start;
  atomic_fetch_add_explicit(&klu_ns[ctx][op],elapsed,memory_order_relaxed);
  atomic_fetch_add_explicit(&klu_calls[ctx][op],1,memory_order_relaxed);
}
extern UF_long __real_trilinos_klu_l_refactor(UF_long*,UF_long*,double*,trilinos_klu_l_symbolic*,trilinos_klu_l_numeric*,trilinos_klu_l_common*);
UF_long __wrap_trilinos_klu_l_refactor(UF_long*p,UF_long*i,double*a,trilinos_klu_l_symbolic*s,trilinos_klu_l_numeric*n,trilinos_klu_l_common*c) {
  unsigned ctx=phase_context;uint64_t t=tick();UF_long r=__real_trilinos_klu_l_refactor(p,i,a,s,n,c);klu_record(ctx,0,t);return r;
}
extern trilinos_klu_l_numeric *__real_trilinos_klu_l_factor(UF_long*,UF_long*,double*,trilinos_klu_l_symbolic*,trilinos_klu_l_common*);
trilinos_klu_l_numeric *__wrap_trilinos_klu_l_factor(UF_long*p,UF_long*i,double*a,trilinos_klu_l_symbolic*s,trilinos_klu_l_common*c) {
  unsigned ctx=phase_context;uint64_t t=tick();trilinos_klu_l_numeric*r=__real_trilinos_klu_l_factor(p,i,a,s,c);klu_record(ctx,1,t);return r;
}
#define WRAP_SOLVE(NAME,OP) \
extern UF_long __real_##NAME(trilinos_klu_l_symbolic*,trilinos_klu_l_numeric*,UF_long,UF_long,double*,trilinos_klu_l_common*); \
UF_long __wrap_##NAME(trilinos_klu_l_symbolic*s,trilinos_klu_l_numeric*n,UF_long ld,UF_long nrhs,double*b,trilinos_klu_l_common*c) { \
 unsigned ctx=phase_context;uint64_t t=tick();UF_long r=__real_##NAME(s,n,ld,nrhs,b,c);klu_record(ctx,OP,t);return r; }
WRAP_SOLVE(trilinos_klu_l_solve,2)
WRAP_SOLVE(trilinos_klu_l_tsolve,3)
__attribute__((destructor)) static void phase_report(void) {
  const char*names[]={"outside","analyze","factor","refactor","solve","combined","certificate"};
  const char*ops[]={"refactor","factor","solve","tsolve"};
  fprintf(stderr,"klu_phase_profile combined_solve_stats_seconds=%.17g",combined_solve_stats_seconds);
  for(unsigned c=0;c<CONTEXTS;++c) {
    fprintf(stderr," %s_certificate_ns=%llu",names[c],(unsigned long long)cert_parent_ns[c]);
    for(unsigned op=0;op<4;++op)fprintf(stderr," %s_klu_%s_calls=%llu %s_klu_%s_ns=%llu",names[c],ops[op],(unsigned long long)atomic_load(&klu_calls[c][op]),names[c],ops[op],(unsigned long long)atomic_load(&klu_ns[c][op]));
  }
  fputc('\n',stderr);
}

/* Diagnostic for a serialized benchmark owner, not production telemetry. */
static size_t solve_calls,certificate_calls,certificate_misses,certificate_errors;
static size_t certificate_legacy_calls,certificate_verdict_calls;
static double solve_seconds,certificate_seconds,standalone_solve_wall_seconds;
static size_t analyze_calls,factor_calls,refactor_calls;
static double analyze_seconds,factor_seconds,refactor_seconds;
static size_t refactor_solve_calls,stats_errors;
static double refactor_solve_wall_seconds;
static size_t public_refactor_calls,refactor_path_calls[KLS_REFACTOR_PATH_UNCHANGED+1];
static double first_public_refactor_wall_seconds;
/* These path totals cover the entire public call. Combined calls include
 * their solve/certificate time; they are not isolated factor kernel times. */
static double refactor_path_wall_seconds[KLS_REFACTOR_PATH_UNCHANGED+1];
static size_t certificate_min_n,certificate_max_n;
static int64_t certificate_min_nnz,certificate_max_nnz;
static double now(void) {
  struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);
  return (double)t.tv_sec+1e-9*(double)t.tv_nsec;
}
/* Read a fixed prefix only: no retained-structure census in these getters.
 * The timer has already stopped before collecting diagnostic metadata. */
static int refactor_stats(kls_solver *s,double elapsed,kls_stats *stats) {
  if(!public_refactor_calls++) first_public_refactor_wall_seconds=elapsed;
  *stats=(kls_stats){0};
  stats->struct_size=offsetof(kls_stats,last_refactor_path)+sizeof(stats->last_refactor_path);
  if(kls_get_stats(s,stats)!=KLS_OK) {++stats_errors;return 0;}
  const int path=(int)stats->last_refactor_path;
  if(path<KLS_REFACTOR_PATH_NONE || path>KLS_REFACTOR_PATH_UNCHANGED) {
    ++stats_errors;return 0;
  }
  ++refactor_path_calls[path];refactor_path_wall_seconds[path]+=elapsed;
  return 1;
}
extern int __real_kls_analyze_csr(kls_solver *,kls_index_type,int64_t,
  const void *,const void *,int,const kls_options *);
int __wrap_kls_analyze_csr(kls_solver *s,kls_index_type index,int64_t n,
    const void *p,const void *cols,int base,const kls_options *options) {
  unsigned saved_context=phase_context;phase_context=ANALYZE;
  double start=now();
  int status=__real_kls_analyze_csr(s,index,n,p,cols,base,options);
  phase_context=saved_context;
  analyze_seconds+=now()-start;++analyze_calls;return status;
}
extern int __real_kls_factor(kls_solver *,const double *);
int __wrap_kls_factor(kls_solver *s,const double *a) {
  unsigned saved_context=phase_context;phase_context=FACTOR;
  double start=now();int status=__real_kls_factor(s,a);
  phase_context=saved_context;
  factor_seconds+=now()-start;++factor_calls;return status;
}
extern int __real_kls_refactor(kls_solver *,const double *);
int __wrap_kls_refactor(kls_solver *s,const double *a) {
  unsigned saved_context=phase_context;phase_context=REFACTOR;
  double start=now();int status=__real_kls_refactor(s,a);
  phase_context=saved_context;
  const double elapsed=now()-start;
  refactor_seconds+=elapsed;++refactor_calls;
  kls_stats stats;refactor_stats(s,elapsed,&stats);return status;
}
extern int __real_kls_solve(kls_solver *,int64_t,const double *,int64_t,double *,int64_t);
int __wrap_kls_solve(kls_solver *s,int64_t nrhs,const double *b,int64_t ldb,
                    double *x,int64_t ldx) {
  unsigned saved_context=phase_context;phase_context=SOLVE;
  double start=now();int status=__real_kls_solve(s,nrhs,b,ldb,x,ldx);
  phase_context=saved_context;
  const double elapsed=now()-start;
  solve_seconds+=elapsed;standalone_solve_wall_seconds+=elapsed;
  ++solve_calls;return status;
}
extern int __real_kls_refactor_solve(kls_solver *,const double *,int64_t,
  const double *,int64_t,double *,int64_t);
int __wrap_kls_refactor_solve(kls_solver *s,const double *a,int64_t nrhs,
  const double *b,int64_t ldb,double *x,int64_t ldx) {
  unsigned saved_context=phase_context;phase_context=COMBINED;
  double start=now();
  int status=__real_kls_refactor_solve(s,a,nrhs,b,ldb,x,ldx);
  phase_context=saved_context;
  const double elapsed=now()-start;
  refactor_solve_wall_seconds+=elapsed;++refactor_solve_calls;
  kls_stats stats;
  if(refactor_stats(s,elapsed,&stats)) {solve_seconds+=stats.solve_seconds;combined_solve_stats_seconds+=stats.solve_seconds;}
  ++solve_calls;return status;
}
extern int __real_kls_componentwise_certify_correction(kls_componentwise_plan **,
  int,size_t,const int64_t *,const int64_t *,const double *,const double *,
  const double *,int,long double *,kls_componentwise_executor,void *,double **);
int __wrap_kls_componentwise_certify_correction(kls_componentwise_plan **plan,
  int threads,size_t n,const int64_t *p,const int64_t *rows,const double *a,
  const double *b,const double *x,int transpose,long double *upper,
  kls_componentwise_executor execute,void *ctx,double **correction) {
  unsigned saved_context=phase_context;phase_context=CERTIFICATE;
  uint64_t cert_start=tick();
  double start=now();
  int status=__real_kls_componentwise_certify_correction(plan,threads,n,p,rows,
      a,b,x,transpose,upper,execute,ctx,correction);
  phase_context=saved_context;cert_parent_ns[saved_context]+=tick()-cert_start;
  certificate_seconds+=now()-start;
  const int64_t nnz=p?p[n]:-1;
  if(!certificate_calls) {
    certificate_min_n=certificate_max_n=n;
    certificate_min_nnz=certificate_max_nnz=nnz;
  } else {
    if(n<certificate_min_n) certificate_min_n=n;
    if(n>certificate_max_n) certificate_max_n=n;
    if(nnz<certificate_min_nnz) certificate_min_nnz=nnz;
    if(nnz>certificate_max_nnz) certificate_max_nnz=nnz;
  }
  ++certificate_legacy_calls;
  ++certificate_calls;
  certificate_misses+=status==0;certificate_errors+=status<0;return status;
}
#ifdef KLS_PROFILE_VERDICT_INTERFACE
extern int __real_kls_componentwise_certify_verdict_correction(kls_componentwise_plan **,
  int,size_t,const int64_t *,const int64_t *,const double *,const double *,
  const double *,int,long double *,kls_componentwise_executor,void *,double **);
int __wrap_kls_componentwise_certify_verdict_correction(kls_componentwise_plan **plan,
  int threads,size_t n,const int64_t *p,const int64_t *rows,const double *a,
  const double *b,const double *x,int transpose,long double *upper,
  kls_componentwise_executor execute,void *ctx,double **correction) {
  unsigned saved_context=phase_context;phase_context=CERTIFICATE;
  uint64_t cert_start=tick();
  double start=now();
  int status=__real_kls_componentwise_certify_verdict_correction(plan,threads,n,p,rows,
      a,b,x,transpose,upper,execute,ctx,correction);
  phase_context=saved_context;cert_parent_ns[saved_context]+=tick()-cert_start;
  certificate_seconds+=now()-start;
  const int64_t nnz=p?p[n]:-1;
  if(!certificate_calls) {
    certificate_min_n=certificate_max_n=n;
    certificate_min_nnz=certificate_max_nnz=nnz;
  } else {
    if(n<certificate_min_n) certificate_min_n=n;
    if(n>certificate_max_n) certificate_max_n=n;
    if(nnz<certificate_min_nnz) certificate_min_nnz=nnz;
    if(nnz>certificate_max_nnz) certificate_max_nnz=nnz;
  }
  ++certificate_verdict_calls;
  ++certificate_calls;
  certificate_misses+=status==0;certificate_errors+=status<0;return status;
}
#endif

__attribute__((destructor)) static void report(void) {
  const double api_wall_seconds=analyze_seconds+factor_seconds+refactor_seconds+
      standalone_solve_wall_seconds+refactor_solve_wall_seconds;
  fprintf(stderr,"solve_profile solve_calls=%zu solve_seconds=%.17g "
      "stats_errors=%zu certificate_calls=%zu certificate_misses=%zu certificate_errors=%zu "
      "certificate_seconds=%.17g refactor_solve_calls=%zu "
      "refactor_solve_wall_seconds=%.17g certificate_min_n=%zu certificate_max_n=%zu "
      "certificate_min_nnz=%lld certificate_max_nnz=%lld "
      "analyze_calls=%zu analyze_seconds=%.17g factor_calls=%zu factor_seconds=%.17g "
      "refactor_calls=%zu refactor_seconds=%.17g standalone_solve_wall_seconds=%.17g "
      "api_wall_seconds=%.17g public_refactor_calls=%zu "
      "first_public_refactor_wall_seconds=%.17g",
      solve_calls,solve_seconds,stats_errors,certificate_calls,
      certificate_misses,certificate_errors,certificate_seconds,refactor_solve_calls,
      refactor_solve_wall_seconds,certificate_min_n,certificate_max_n,
      (long long)certificate_min_nnz,(long long)certificate_max_nnz,
      analyze_calls,analyze_seconds,factor_calls,factor_seconds,
      refactor_calls,refactor_seconds,standalone_solve_wall_seconds,api_wall_seconds,
      public_refactor_calls,first_public_refactor_wall_seconds);
  for(int path=KLS_REFACTOR_PATH_NONE;path<=KLS_REFACTOR_PATH_UNCHANGED;++path)
    fprintf(stderr," refactor_path_%d_calls=%zu refactor_path_%d_wall_seconds=%.17g",
        path,refactor_path_calls[path],path,refactor_path_wall_seconds[path]);
  fprintf(stderr," certificate_legacy_calls=%zu certificate_verdict_calls=%zu",
      certificate_legacy_calls,certificate_verdict_calls);
  fputc('\n',stderr);
}
