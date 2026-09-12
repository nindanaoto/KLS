/* KLS inspects native factor metadata for every column, including singleton
 * BTF blocks. Poison allocations so untouched fields fail deterministically. */
#include "trilinos_klu_decl.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void *poison_malloc(size_t size) {
  void *p = malloc(size);
  if (p != NULL) memset(p, 0xa5, size);
  return p;
}

int main(void) {
  UF_long ap[] = {0, 1, 3, 5, 6};
  UF_long ai[] = {0, 1, 2, 1, 2, 3};
  double ax[] = {2, 4, 1, 1, 3, 5};
  trilinos_klu_l_common common;
  if (!trilinos_klu_l_defaults(&common)) return 1;
  common.btf = 1;
  trilinos_klu_l_symbolic *s = trilinos_klu_l_analyze(4, ap, ai, &common);
  if (s == NULL) return 1;
  common.malloc_memory = poison_malloc;
  trilinos_klu_l_numeric *n = trilinos_klu_l_factor(ap, ai, ax, s, &common);
  int ok = n != NULL;
  int singletons = 0;
  if (n != NULL) {
    for (UF_long block = 0; block < s->nblocks; ++block) {
      UF_long col = s->R[block];
      if (s->R[block + 1] - col != 1) continue;
      ++singletons;
      if (n->Lip[col] != 0 || n->Uip[col] != 0 || n->Llen[col] != 0 ||
          n->Ulen[col] != 0 || n->LUsize[block] != 0 || n->LUbx[block] != NULL) {
        fprintf(stderr, "singleton column %ld has unset factor metadata\n", (long)col);
        ok = 0;
      }
    }
    for (int update = 0; update < 2; ++update) {
      ax[0] = 2 + update;
      double rhs[] = {2 + update, 5, 4, 5};
      if ((update && !trilinos_klu_l_refactor(ap, ai, ax, s, n, &common)) ||
          !trilinos_klu_l_solve(s, n, 4, 1, rhs, &common)) ok = 0;
      for (int row = 0; row < 4; ++row)
        if (!isfinite(rhs[row]) || fabs(rhs[row] - 1) > 1e-12) ok = 0;
    }
  }
  if (singletons != 2) ok = 0;
  trilinos_klu_l_free_numeric(&n, &common);
  trilinos_klu_l_free_symbolic(&s, &common);
  return ok ? 0 : 1;
}
