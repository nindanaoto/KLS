#define _POSIX_C_SOURCE 200809L
#include "trilinos_klu_internal.h"
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

/* Each Common has a genuinely different allocator family. Prefix tags prove
 * that AMD workspaces are released by their own allocator even concurrently. */
typedef union allocation_header {
  max_align_t alignment;
  struct {size_t bytes;unsigned family;} value;
} allocation_header;
static _Atomic long live[2],calls[2];
static _Atomic int bad;
static _Thread_local int fail_after=-1;
static void *owned_allocate(unsigned family,size_t bytes) {
  if(fail_after==0)return NULL;
  if(fail_after>0)--fail_after;
  if(bytes>SIZE_MAX-sizeof(allocation_header))return NULL;
  allocation_header *p=malloc(sizeof(*p)+bytes);if(!p)return NULL;
  p->value.bytes=bytes;p->value.family=family;
  atomic_fetch_add(live+family,1);atomic_fetch_add(calls+family,1);return p+1;
}
static void owned_release(unsigned family,void *memory) {
  if(!memory)return;
  allocation_header *p=(allocation_header*)memory-1;
  if(p->value.family!=family)atomic_store(&bad,1);
  atomic_fetch_sub(live+p->value.family,1);free(p);
}
static void *owned_resize(unsigned family,void *memory,size_t bytes) {
  if(!memory)return owned_allocate(family,bytes);
  if(fail_after==0)return NULL;
  if(fail_after>0)--fail_after;
  if(bytes>SIZE_MAX-sizeof(allocation_header))return NULL;
  allocation_header *p=(allocation_header*)memory-1;
  if(p->value.family!=family)atomic_store(&bad,1);
  p=realloc(p,sizeof(*p)+bytes);if(!p)return NULL;p->value.bytes=bytes;return p+1;
}
#define FAMILY(ID) \
static void *allocate##ID(size_t n){return owned_allocate(ID,n);} \
static void release##ID(void *p){owned_release(ID,p);} \
static void *resize##ID(void *p,size_t n){return owned_resize(ID,p,n);} \
static void *zero##ID(size_t n,size_t z){if(z&&n>SIZE_MAX/z)return NULL;void *p=owned_allocate(ID,n*z);if(p)memset(p,0,n*z);return p;}
FAMILY(0)
FAMILY(1)
static void initialize(TRILINOS_KLU_common *c,unsigned family) {
  TRILINOS_KLU_defaults(c);c->ordering=0;
  c->malloc_memory=family?allocate1:allocate0;c->free_memory=family?release1:release0;
  c->realloc_memory=family?resize1:resize0;c->calloc_memory=family?zero1:zero0;
}
static TRILINOS_KLU_symbolic *analyze(TRILINOS_KLU_common *c) {
  Int p[6]={0,5,10,15,20,25},row[25];
  for(int i=0;i<25;++i)row[i]=i%5;
  return TRILINOS_KLU_analyze(5,p,row,c);
}
static void *worker(void *arg) {
  unsigned family=(unsigned)(uintptr_t)arg;
  for(int i=0;i<64;++i) {
    TRILINOS_KLU_common c;initialize(&c,family);
    TRILINOS_KLU_symbolic *s=analyze(&c);
    if(!s || c.status!=TRILINOS_KLU_OK){atomic_store(&bad,1);return NULL;}
    TRILINOS_KLU_free_symbolic(&s,&c);
  }return NULL;
}
int main(void) {
  void *(*original_malloc)(size_t)=trilinos_amd_malloc;
  void (*original_free)(void*)=trilinos_amd_free;
  void *(*original_calloc)(size_t,size_t)=trilinos_amd_calloc;
  void *(*original_realloc)(void*,size_t)=trilinos_amd_realloc;
  pthread_t threads[4];
  for(uintptr_t i=0;i<4;++i)if(pthread_create(threads+i,NULL,worker,(void*)(i&1u)))return 1;
  for(int i=0;i<4;++i)if(pthread_join(threads[i],NULL))return 1;
  int failed_sites=0;
  for(int site=0;site<64;++site) {
    TRILINOS_KLU_common c;initialize(&c,0);fail_after=site;
    TRILINOS_KLU_symbolic *s=analyze(&c);fail_after=-1;
    if(s){TRILINOS_KLU_free_symbolic(&s,&c);break;}
    if(c.status!=TRILINOS_KLU_OUT_OF_MEMORY)return 1;
    ++failed_sites;
  }
  if(atomic_load(&bad)||atomic_load(live)||atomic_load(live+1)||
     !atomic_load(calls)||!atomic_load(calls+1)||failed_sites<5||failed_sites>=64||
     trilinos_amd_malloc!=original_malloc||trilinos_amd_free!=original_free||
     trilinos_amd_calloc!=original_calloc||trilinos_amd_realloc!=original_realloc)return 1;
  printf("PASS KLU allocator ownership: four concurrent owners, two tagged allocator families, %d OOM sites, public AMD hooks unchanged\n",failed_sites);
  return 0;
}
