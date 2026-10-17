#define _POSIX_C_SOURCE 200809L
#define KLS_COARSE_COUNT 3
static const char *const kls_coarse_names[] = {"owner", "child", "worker"};
static const char *const kls_coarse_roles[] = {"caller", "caller", "worker"};
#include "../scripts/coarse_profile.inc"
#include <assert.h>
#include <fenv.h>
#include <pthread.h>
#if defined(__x86_64__) || defined(__i386__)
#include <xmmintrin.h>
#endif

static int child(int early) {
  KLS_COARSE_SCOPE(1);
  void *p = __wrap_malloc(25);
  void *q = __wrap_calloc(3, 8);
  assert(p && q);
  p = __wrap_realloc(p, 80);
  assert(p);
  void *r = __wrap_aligned_alloc(64, 64), *s = NULL;
  assert(r && __wrap_posix_memalign(&s, 64, 64) == 0 && s);
  __wrap_free(p); __wrap_free(q); __wrap_free(r); __wrap_free(s);
  errno = EDOM;
  if (early) return 7;
  return 9;
}

static void *worker(void *unused) {
  (void)unused;
  KLS_COARSE_SCOPE(2);
  for (int i = 0; i < 100; ++i) {
    assert(child(i & 1) == ((i & 1) ? 7 : 9));
    assert(errno == EDOM);
  }
  return NULL;
}

static void owner(void) {
  KLS_COARSE_SCOPE(0);
  fenv_t original;
  assert(fegetenv(&original) == 0);
  const int modes[] = {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO};
  for (unsigned m = 0; m < 4; ++m) {
    assert(fesetround(modes[m]) == 0);
    for (unsigned flags = 0; flags < 64; ++flags) {
#if defined(__x86_64__) || defined(__i386__)
      unsigned control = (_mm_getcsr() & ~63u) | flags;
      _mm_setcsr(control);
      unsigned short sw_before, cw_before, sw_after, cw_after;
      __asm__ volatile("fnstsw %0" : "=am"(sw_before));
      __asm__ volatile("fnstcw %0" : "=m"(cw_before));
#endif
      assert(child(flags & 1) == ((flags & 1) ? 7 : 9));
      assert(errno == EDOM && fegetround() == modes[m]);
#if defined(__x86_64__) || defined(__i386__)
      __asm__ volatile("fnstsw %0" : "=am"(sw_after));
      __asm__ volatile("fnstcw %0" : "=m"(cw_after));
      assert(_mm_getcsr() == control && sw_before == sw_after && cw_before == cw_after);
#endif
    }
  }
  assert(fesetenv(&original) == 0);
  pthread_t threads[2];
  for (int i = 0; i < 2; ++i) assert(pthread_create(threads + i, NULL, worker, NULL) == 0);
  for (int i = 0; i < 2; ++i) assert(pthread_join(threads[i], NULL) == 0);
}

int main(void) {
  owner();
  assert(kls_coarse_top == NULL && atomic_load(&kls_coarse_errors) == 0);
  assert(atomic_load(&kls_coarse_data[0].calls) == 1);
  assert(atomic_load(&kls_coarse_data[1].calls) == 456);
  assert(atomic_load(&kls_coarse_data[1].threads) == 3);
  assert(atomic_load(&kls_coarse_data[1].allocation_calls) == 5 * 456);
  assert(atomic_load(&kls_coarse_data[1].requested_bytes) == 257 * 456);
  assert(atomic_load(&kls_coarse_data[1].free_calls) == 4 * 456);
  assert(atomic_load(&kls_coarse_data[1].allocation_failures) == 0);
  assert(atomic_load(&kls_coarse_data[2].calls) == 2);
  assert(atomic_load(&kls_coarse_data[2].threads) == 2);
  for (unsigned i = 0; i < KLS_COARSE_COUNT; ++i) {
    assert(atomic_load(&kls_coarse_data[i].inclusive_ns) > 0);
    assert(atomic_load(&kls_coarse_data[i].exclusive_ns) <=
           atomic_load(&kls_coarse_data[i].inclusive_ns));
  }
  puts("coarse profile controls passed: nesting, early returns, allocation counts, threads and FP state");
  return 0;
}
