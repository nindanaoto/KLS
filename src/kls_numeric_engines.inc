/* Internal numeric factorization and refactor engine implementation.
 * Included by kls.c at the original definition point to preserve a single
 * translation unit and its optimization behavior.
 *
 * Include order is significant.  Preparation and repair define the shared
 * metadata builders; row and EGraph define executors; scheduling connects
 * those executors; first-factor and PTS consume all three; compact solve is
 * last because it can reuse every retained representation.  Do not reorder
 * these files merely to resolve a declaration—place a narrow forward
 * declaration at the ownership boundary instead.
 */

#include "kls_numeric_prepare_repair.inc"
#include "kls_row_refactor.inc"
#include "kls_egraph_refactor.inc"
#include "kls_numeric_schedule.inc"
#include "kls_first_factor.inc"
#include "kls_pts.inc"
#include "kls_compact_solve.inc"
