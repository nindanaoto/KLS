/* Diagnostic builds only. Integer clocks preserve the caller's FP state.
 * Exclusive time subtracts nested scopes on this thread, never other workers.
 * Allocation bytes are requests at wrapped link boundaries, not live memory.
 * Define KLS_COARSE_COUNT, kls_coarse_names and kls_coarse_roles before use. */
#include <errno.h>
#include <stdint.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

_Static_assert(KLS_COARSE_COUNT <= 64, "thread-local scope membership mask");
typedef struct kls_coarse_scope {
  struct kls_coarse_scope *parent;
  uint64_t start, children;
  unsigned slot;
  int in_api, api_root;
} kls_coarse_scope;
typedef struct {
  _Atomic uint64_t calls, threads, inclusive_ns, exclusive_ns;
  _Atomic uint64_t api_calls, api_inclusive_ns, api_exclusive_ns;
  _Atomic uint64_t api_root_calls, api_root_inclusive_ns;
  _Atomic uint64_t allocation_calls, requested_bytes, allocation_failures;
  _Atomic uint64_t free_calls;
} kls_coarse_totals;
static kls_coarse_totals kls_coarse_data[KLS_COARSE_COUNT];
static _Atomic uint64_t kls_coarse_errors;
static _Thread_local kls_coarse_scope *kls_coarse_top;
static _Thread_local uint64_t kls_coarse_seen;

static void kls_coarse_add(_Atomic uint64_t *counter, uint64_t n) {
  uint64_t before = atomic_fetch_add_explicit(counter, n, memory_order_relaxed);
  if (before > UINT64_MAX - n)
    atomic_fetch_add_explicit(&kls_coarse_errors, 1, memory_order_relaxed);
}

static uint64_t kls_coarse_now(void) {
  struct timespec t = {0};
  int saved = errno;
  int status = clock_gettime(CLOCK_MONOTONIC, &t);
  errno = saved;
  if (status || t.tv_sec < 0 || (uint64_t)t.tv_sec > UINT64_MAX / 1000000000) {
    kls_coarse_add(&kls_coarse_errors, 1);
    return 0;
  }
  return (uint64_t)t.tv_sec * UINT64_C(1000000000) + (uint64_t)t.tv_nsec;
}

static void kls_coarse_enter(kls_coarse_scope *s, unsigned slot) {
  *s = (kls_coarse_scope){.parent = kls_coarse_top,
                          .start = kls_coarse_now(), .slot = slot,
                          .in_api = (kls_coarse_top && kls_coarse_top->in_api) ||
                                    kls_coarse_roles[slot][0] == 'a',
                          .api_root = !(kls_coarse_top && kls_coarse_top->in_api) &&
                                      kls_coarse_roles[slot][0] == 'a'};
  kls_coarse_top = s;
  kls_coarse_add(&kls_coarse_data[slot].calls, 1);
  if (s->in_api) kls_coarse_add(&kls_coarse_data[slot].api_calls, 1);
  if (s->api_root) kls_coarse_add(&kls_coarse_data[slot].api_root_calls, 1);
  uint64_t mask = UINT64_C(1) << slot;
  if (!(kls_coarse_seen & mask)) {
    kls_coarse_seen |= mask;
    kls_coarse_add(&kls_coarse_data[slot].threads, 1);
  }
}

static void kls_coarse_leave(kls_coarse_scope *s) {
  uint64_t end = kls_coarse_now();
  if (kls_coarse_top != s || !s->start || end < s->start ||
      end - s->start < s->children) {
    kls_coarse_add(&kls_coarse_errors, 1);
  } else {
    uint64_t elapsed = end - s->start;
    kls_coarse_add(&kls_coarse_data[s->slot].inclusive_ns, elapsed);
    kls_coarse_add(&kls_coarse_data[s->slot].exclusive_ns,
                   elapsed - s->children);
    if (s->in_api) {
      kls_coarse_add(&kls_coarse_data[s->slot].api_inclusive_ns, elapsed);
      kls_coarse_add(&kls_coarse_data[s->slot].api_exclusive_ns,
                     elapsed - s->children);
    }
    if (s->api_root)
      kls_coarse_add(&kls_coarse_data[s->slot].api_root_inclusive_ns, elapsed);
    if (s->parent) {
      if (s->parent->children > UINT64_MAX - elapsed)
        kls_coarse_add(&kls_coarse_errors, 1);
      else s->parent->children += elapsed;
    }
  }
  kls_coarse_top = s->parent;
}

#define KLS_COARSE_SCOPE(slot) \
  kls_coarse_scope kls_coarse_local \
    __attribute__((cleanup(kls_coarse_leave))); \
  kls_coarse_enter(&kls_coarse_local, (slot))

static void kls_coarse_allocation(size_t bytes, int failed) {
  if (!kls_coarse_top) return;
  kls_coarse_totals *d = &kls_coarse_data[kls_coarse_top->slot];
  kls_coarse_add(&d->allocation_calls, 1);
  kls_coarse_add(&d->requested_bytes, bytes);
  kls_coarse_add(&d->allocation_failures, failed != 0);
}

void *__real_malloc(size_t);
void *__real_calloc(size_t, size_t);
void *__real_realloc(void *, size_t);
void __real_free(void *);
void *__real_aligned_alloc(size_t, size_t);
int __real_posix_memalign(void **, size_t, size_t);

void *__wrap_malloc(size_t n) {
  void *p = __real_malloc(n);
  kls_coarse_allocation(n, p == NULL && n != 0);
  return p;
}
void *__wrap_calloc(size_t n, size_t width) {
  void *p = __real_calloc(n, width);
  if (width && n > SIZE_MAX / width) {
    kls_coarse_add(&kls_coarse_errors, 1);
    kls_coarse_allocation(0, 1);
  } else kls_coarse_allocation(n * width, p == NULL && n && width);
  return p;
}
void *__wrap_realloc(void *old, size_t n) {
  void *p = __real_realloc(old, n);
  kls_coarse_allocation(n, p == NULL && n != 0);
  return p;
}
void __wrap_free(void *p) {
  if (kls_coarse_top)
    kls_coarse_add(&kls_coarse_data[kls_coarse_top->slot].free_calls, 1);
  __real_free(p);
}
void *__wrap_aligned_alloc(size_t alignment, size_t n) {
  void *p = __real_aligned_alloc(alignment, n);
  kls_coarse_allocation(n, p == NULL && n != 0);
  return p;
}
int __wrap_posix_memalign(void **p, size_t alignment, size_t n) {
  int status = __real_posix_memalign(p, alignment, n);
  kls_coarse_allocation(n, status != 0);
  return status;
}

__attribute__((destructor)) static void kls_coarse_report(void) {
  for (unsigned i = 0; i < KLS_COARSE_COUNT; ++i) {
    kls_coarse_totals *d = &kls_coarse_data[i];
    fprintf(stderr, "kls_coarse slot=%u name=%s role=%s calls=%llu threads=%llu "
      "inclusive_ns=%llu exclusive_ns=%llu allocation_calls=%llu "
      "requested_bytes=%llu allocation_failures=%llu free_calls=%llu "
      "api_calls=%llu api_inclusive_ns=%llu api_exclusive_ns=%llu "
      "api_root_calls=%llu api_root_inclusive_ns=%llu\n", i,
      kls_coarse_names[i], kls_coarse_roles[i],
      (unsigned long long)atomic_load(&d->calls),
      (unsigned long long)atomic_load(&d->threads),
      (unsigned long long)atomic_load(&d->inclusive_ns),
      (unsigned long long)atomic_load(&d->exclusive_ns),
      (unsigned long long)atomic_load(&d->allocation_calls),
      (unsigned long long)atomic_load(&d->requested_bytes),
      (unsigned long long)atomic_load(&d->allocation_failures),
      (unsigned long long)atomic_load(&d->free_calls),
      (unsigned long long)atomic_load(&d->api_calls),
      (unsigned long long)atomic_load(&d->api_inclusive_ns),
      (unsigned long long)atomic_load(&d->api_exclusive_ns),
      (unsigned long long)atomic_load(&d->api_root_calls),
      (unsigned long long)atomic_load(&d->api_root_inclusive_ns));
  }
  fprintf(stderr, "kls_coarse_summary errors=%llu\n",
    (unsigned long long)atomic_load(&kls_coarse_errors));
}
