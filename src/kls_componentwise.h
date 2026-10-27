#ifndef KLS_COMPONENTWISE_H
#define KLS_COMPONENTWISE_H
#include <stddef.h>
#include <stdint.h>
/* Contract constants, not CPU tuning parameters. */
#define KLS_COMPONENTWISE_TOLERANCE_DENOMINATOR 100000000
#define KLS_COMPONENTWISE_MAX_CORRECTIONS 3
/* Private optional row plan. The owner must reset it before changing the CSC
 * pattern, including entry order. Values may change on every call. A thread
 * budget of one uses serial rows and never launches timing-trial workers. */
typedef struct kls_componentwise_plan kls_componentwise_plan;
/* Synchronous private executor: invoke work(job,tid) exactly once for each
 * tid in [0,threads), and join every invocation before returning, even on
 * failure. Return one only when the whole dispatch completed. */
typedef int (*kls_componentwise_executor)(void *context,
    void (*work)(void *,int),void *job,int threads);
void kls_componentwise_plan_reset(kls_componentwise_plan *plan);
void kls_componentwise_plan_destroy(kls_componentwise_plan *plan);
int kls_componentwise_certify_planned(kls_componentwise_plan **plan,int threads,
    size_t n,const int64_t *p,const int64_t *rows,const double *a,
    const double *b,const double *x,int transpose,long double *upper);
int kls_componentwise_certify_executor(kls_componentwise_plan **plan,int threads,
    size_t n,const int64_t *p,const int64_t *rows,const double *a,
    const double *b,const double *x,int transpose,long double *upper,
    kls_componentwise_executor executor,void *context);
/* Lazily allocate 2*n doubles at *correction when an interval pass is needed.
 * On an uncertified return with finite *upper, its first n values contain the
 * freshly computed correction RHS. The caller owns/frees the workspace, and
 * may reuse it only at the same n. Each call checks the current values/RHS/x. */
int kls_componentwise_certify_correction(kls_componentwise_plan **plan,int threads,
    size_t n,const int64_t *p,const int64_t *rows,const double *a,
    const double *b,const double *x,int transpose,long double *upper,
    kls_componentwise_executor executor,void *context,double **correction);
/* 1 certified, 0 not certified, -1 allocation failure. Residual, when
 * requested, is an ordinary double correction RHS, never a certificate. */
int kls_componentwise_certify(size_t n, const int64_t *p,
    const int64_t *rows, const double *a, const double *b, const double *x,
    int transpose, double *residual, long double *upper);
/* Private verdict/correction interface. All current rows must satisfy the
 * same componentwise tolerance. The successful bound is not part of this
 * interface's output contract. On rejection, failure_upper and correction
 * have exactly the existing correction-interface semantics. Legacy detailed
 * interfaces above retain their accepted-bound contract. */
int kls_componentwise_certify_verdict_correction(kls_componentwise_plan **plan,int threads,
    size_t n,const int64_t *p,const int64_t *rows,const double *a,
    const double *b,const double *x,int transpose,long double *failure_upper,
    kls_componentwise_executor executor,void *context,double **correction);
#endif
