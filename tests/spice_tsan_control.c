#define _GNU_SOURCE
#include <stdio.h>
#include <pthread.h>
__attribute__((target_clones("default", "arch=x86-64-v3", "arch=x86-64-v4"), noinline))
static int value(int x) { return x + 1; }
#ifdef RACE
static volatile int shared;
static pthread_barrier_t barrier;
static void *writer(void *unused) {
  (void)unused;
  pthread_barrier_wait(&barrier);
  for(int i=0;i<10000;++i) shared=i;
  return NULL;
}
#endif
int main(void) {
  fputs("CONTROL_MAIN_ENTERED\n",stderr);
#ifdef RACE
  pthread_t a,b;
  pthread_barrier_init(&barrier,NULL,3);
  pthread_create(&a,NULL,writer,NULL);
  pthread_create(&b,NULL,writer,NULL);
  pthread_barrier_wait(&barrier);
  pthread_join(a,NULL);pthread_join(b,NULL);
  pthread_barrier_destroy(&barrier);
#endif
  return value(0)!=1;
}
