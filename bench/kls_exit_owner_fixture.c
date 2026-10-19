/* Standalone diagnostic contract fixture; no production solver is linked. */
#define _GNU_SOURCE
#include "kls/kls.h"
#include <errno.h>
#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <xmmintrin.h>
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"exit-owner fixture line %d: %s\n",__LINE__,#c); abort(); } } while (0)
struct kls_solver { int token; };
static kls_solver *expected;
static int expected_errno, expected_round, expected_flags, mode, getter_calls, real_calls, real_effect;
static unsigned expected_csr, real_csr;
static unsigned short expected_cw;
void __wrap_kls_destroy(kls_solver *);
int kls_get_stats(const kls_solver *solver, kls_stats *stats) {
  CHECK(solver == expected && solver->token == 19);
  CHECK(stats->struct_size == offsetof(kls_stats,compact_solve_index_bytes)+sizeof(stats->compact_solve_index_bytes));
  ++getter_calls;
  memset(stats,0,sizeof(*stats));stats->struct_size=sizeof(*stats);
  stats->n=65;stats->nnz=433;stats->nnz_l=193;stats->nnz_u=193;
  stats->last_refactor_path=KLS_REFACTOR_PATH_MAPPED;
  stats->compact_solve_index_bytes=2;
  /* The observer must mask traps and discard these getter/IO side effects. */
  CHECK(feraiseexcept(FE_INVALID)==0);
  CHECK(fesetround(FE_DOWNWARD)==0);
  _mm_setcsr(_mm_getcsr() ^ ((1u<<15)|(1u<<6)));
  errno=999;
  return mode ? KLS_ERR_INVALID_ARGUMENT : KLS_OK;
}
void __real_kls_destroy(kls_solver *solver) {
  const int seen_errno=errno;
  unsigned short cw;__asm__ volatile("fnstcw %0":"=m"(cw));
  CHECK(solver==expected && seen_errno==expected_errno);
  CHECK(_mm_getcsr()==expected_csr && cw==expected_cw);
  CHECK(fegetround()==expected_round && fetestexcept(FE_ALL_EXCEPT)==expected_flags);
  CHECK(getter_calls==(solver!=NULL));
  ++real_calls;
  if(real_effect) CHECK(feraiseexcept(FE_DIVBYZERO)==0);
  real_csr=_mm_getcsr();
  errno=71;
}
int main(void) {
  fenv_t original;CHECK(fegetenv(&original)==0);
  const unsigned original_csr=_mm_getcsr();
  const int rounds[]={FE_TONEAREST,FE_DOWNWARD,FE_UPWARD,FE_TOWARDZERO};
  struct kls_solver solver={19};unsigned contexts=0;
  for(int rm=0;rm<4;++rm)for(int denorm=0;denorm<4;++denorm)
    for(int flags=0;flags<2;++flags)for(int failure=0;failure<2;++failure)
      for(int present=0;present<2;++present) {
        CHECK(fesetenv(&original)==0);CHECK(fesetround(rounds[rm])==0);
        CHECK(feclearexcept(FE_ALL_EXCEPT)==0);
        if(flags) CHECK(feraiseexcept(FE_INEXACT)==0);
        unsigned csr=_mm_getcsr()&~((1u<<15)|(1u<<6));
        if(denorm&1)csr|=1u<<15;
        if(denorm&2)csr|=1u<<6;
        _mm_setcsr(csr);
        expected=present?&solver:NULL;mode=failure;getter_calls=real_calls=0;real_effect=1;
        expected_csr=_mm_getcsr();expected_round=fegetround();expected_flags=fetestexcept(FE_ALL_EXCEPT);
        __asm__ volatile("fnstcw %0":"=m"(expected_cw));expected_errno=449;errno=expected_errno;
        __wrap_kls_destroy(expected);
        CHECK(errno==71 && real_calls==1 && getter_calls==present);
        CHECK(fegetround()==expected_round && _mm_getcsr()==real_csr);
        CHECK(fetestexcept(FE_ALL_EXCEPT)==(expected_flags|FE_DIVBYZERO));
        ++contexts;
      }
  /* Unmasked invalid traps are restored before the real destructor, while
   * the diagnostic getter itself raises invalid safely under feholdexcept. */
  CHECK(fesetenv(&original)==0);CHECK(feclearexcept(FE_ALL_EXCEPT)==0);
  CHECK(feenableexcept(FE_INVALID)>=0);
  expected=&solver;mode=0;getter_calls=real_calls=0;real_effect=0;
  expected_csr=_mm_getcsr();expected_round=fegetround();expected_flags=0;
  __asm__ volatile("fnstcw %0":"=m"(expected_cw));expected_errno=449;errno=expected_errno;
  __wrap_kls_destroy(expected);CHECK(errno==71 && real_calls==1 && getter_calls==1);
  CHECK(_mm_getcsr()==expected_csr && fegetround()==expected_round);
  CHECK(fesetenv(&original)==0);_mm_setcsr(original_csr);
  printf("PASS exit-owner pointer, order, status, errno and FP contract: %u contexts plus unmasked trap\n",contexts);
  return 0;
}
