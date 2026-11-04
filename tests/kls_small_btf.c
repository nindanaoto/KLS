static _Thread_local unsigned long long caller_local_frames,caller_heap_frames;
#define KLS_CALLER_LOCAL_FRAME_TEST_EVENT(local) ((local)?++caller_local_frames:++caller_heap_frames)
static _Thread_local unsigned long long compact_btf_descriptor_calls;
static _Thread_local int descriptor_allocation_failure;
#define KLS_COMPACT_BTF_DESCRIPTOR_TEST_EVENT() (++compact_btf_descriptor_calls)
#define KLS_BTF_DESCRIPTOR_MALLOC(n) (descriptor_allocation_failure ? (descriptor_allocation_failure=0,(void*)0) : malloc(n))
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
static _Thread_local unsigned long long singleton_calls;
#define KLS_COMPACT_BTF_SINGLETON_TEST_EVENT() (++singleton_calls)
#include "../src/kls.c"
#include "kls_small_btf_oracles.inc"
#include "kls_btf_descriptor_checks.inc"
#include "kls_small_frame_checks.inc"
#define KLS_SCALED_BTF_MAPPED_TEST 1
#include "kls_small_btf_native.inc"
#include "kls_small_btf_selection.inc"
int main(void) {if(!small_frame_checks())return 1;if(native_btf_checks())return 1;int result=small_btf_selection_checks();if(!singleton_calls){fprintf(stderr,"singleton branch was not exercised\n");return 1;}if(!compact_btf_descriptor_calls)return 1;printf("BTF_DESCRIPTOR_ORACLES activation=%llu ownership=1 allocation_fallback=1 geometry=1\n",compact_btf_descriptor_calls);return result;}
