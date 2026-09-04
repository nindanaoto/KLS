/* Internal numeric factorization and refactor engine implementation.
 * Included by kls.c at the original definition point to preserve a single
 * translation unit and its optimization behavior.
 */

#include "kls_numeric_prepare_repair.inc"
#include "kls_row_refactor.inc"
#include "kls_egraph_refactor.inc"
#include "kls_numeric_schedule.inc"
#include "kls_first_factor.inc"
#include "kls_pts.inc"
#include "kls_compact_solve.inc"
