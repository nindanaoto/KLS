/* Exercise simultaneous first feature queries without a prior warmup. */
#include "../src/kls.c"
#undef NDEBUG
#include <assert.h>

#if KLS_HAVE_AVX512_KERNELS
enum { FEATURE_THREADS = 8 };
static pthread_barrier_t start;
static int observed[FEATURE_THREADS];

static void *query(void *arg) {
  const size_t index=(size_t)(uintptr_t)arg;
  pthread_barrier_wait(&start);
  int value=kls_avx512_scatter_enabled();
  for(int i=0;i<1024;++i) assert(kls_avx512_scatter_enabled()==value);
  observed[index]=value;
  return NULL;
}
#endif

int main(void) {
#if KLS_HAVE_AVX512_KERNELS
  pthread_t workers[FEATURE_THREADS];
  assert(pthread_barrier_init(&start,NULL,FEATURE_THREADS+1)==0);
  for(size_t i=0;i<FEATURE_THREADS;++i)
    assert(pthread_create(&workers[i],NULL,query,(void *)(uintptr_t)i)==0);
  pthread_barrier_wait(&start);
  for(size_t i=0;i<FEATURE_THREADS;++i)assert(pthread_join(workers[i],NULL)==0);
  const int expected=kls_avx512_scatter_enabled();
  for(size_t i=0;i<FEATURE_THREADS;++i)assert(observed[i]==expected);
  assert(pthread_barrier_destroy(&start)==0);
#endif
  return 0;
}
