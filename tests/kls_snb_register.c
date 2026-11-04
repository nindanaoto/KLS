/* Complete edge arithmetic and FP-status oracle, with actual activation. */
static _Thread_local unsigned long long register_calls;
#define KLS_SNB_REGISTER_TEST_EVENT(width) (++register_calls)
#include "../src/kls.c"
#undef NDEBUG
#include <assert.h>
#include "kls_snb_register_reference.inc"
#include "kls_snb_register_checks.inc"
int main(void) {
  register_trsm_checks();
#if defined(__AVX512F__) && defined(__AVX512VL__) && FLT_EVAL_METHOD==0
  assert(register_calls>0);
#endif
  return 0;
}
