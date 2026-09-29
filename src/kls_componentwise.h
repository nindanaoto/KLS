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
void kls_componentwise_plan_reset(kls_componentwise_plan *plan);
void kls_componentwise_plan_destroy(kls_componentwise_plan *plan);
int kls_componentwise_certify_planned(kls_componentwise_plan **plan,int threads,
    size_t n,const int64_t *p,const int64_t *rows,const double *a,
    const double *b,const double *x,int transpose,long double *upper);
/* 1 certified, 0 not certified, -1 allocation failure. Residual, when
 * requested, is an ordinary double correction RHS, never a certificate. */
int kls_componentwise_certify(size_t n, const int64_t *p,
    const int64_t *rows, const double *a, const double *b, const double *x,
    int transpose, double *residual, long double *upper);
#endif
