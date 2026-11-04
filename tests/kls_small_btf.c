/* Include the private implementation to compare actual native payloads and FP
 * flags against the preserved scalar factor and vendored solve, per ordering. */
static _Thread_local unsigned long long compact_btf_factor_calls,compact_btf_solve_calls;
static _Thread_local unsigned long long direct_btf_solve_calls;
#define KLS_COMPACT_BTF_FACTOR_TEST_EVENT() (++compact_btf_factor_calls)
#define KLS_COMPACT_BTF_SOLVE_TEST_EVENT() (++compact_btf_solve_calls)
#define KLS_DIRECT_BTF_SOLVE_TEST_EVENT(s) (++direct_btf_solve_calls)
static _Thread_local unsigned long long small_btf_proposals,small_btf_selected;
#define KLS_SMALL_BTF_AUTO_TEST_EVENT(c,trial,score,accepted) ((void)(c),(void)(trial),(void)(score),++small_btf_proposals,small_btf_selected+=(accepted))
static _Thread_local int small_btf_fail_pattern;
#define KLS_SMALL_BTF_PATTERN_MALLOC(n) (small_btf_fail_pattern ? (small_btf_fail_pattern=0,(void*)0) : malloc(n))
#include "../src/kls.c"
#include "kls_small_btf_oracles.inc"
#define KLS_SCALED_BTF_MAPPED_TEST 1
#include "kls_small_btf_native.inc"
#include "kls_small_btf_selection.inc"
int main(void) {if(native_btf_checks())return 1;return small_btf_selection_checks();}
