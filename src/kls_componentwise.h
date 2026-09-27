#ifndef KLS_COMPONENTWISE_H
#define KLS_COMPONENTWISE_H
#include <stddef.h>
#include <stdint.h>
/* Contract constants, not CPU tuning parameters. */
#define KLS_COMPONENTWISE_TOLERANCE_DENOMINATOR 100000000
#define KLS_COMPONENTWISE_MAX_CORRECTIONS 3
/* 1 certified, 0 not certified, -1 allocation failure. Residual, when
 * requested, is an ordinary double correction RHS, never a certificate. */
int kls_componentwise_certify(size_t n, const int64_t *p,
    const int64_t *rows, const double *a, const double *b, const double *x,
    int transpose, double *residual, long double *upper);
#endif
