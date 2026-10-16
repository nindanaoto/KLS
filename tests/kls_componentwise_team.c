#undef NDEBUG
#include <assert.h>
#include <errno.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static int allocation,allocation_fail,live,initialization,initialization_fail;
static int creation,creation_fail;
static void *fault_calloc(size_t n,size_t size) {
  if(++allocation==allocation_fail) return NULL;
  void *p=calloc(n,size);if(p) ++live;return p;
}
static void fault_free(void *p) {if(p) --live;free(p);}
static int fault_mutex_init(pthread_mutex_t *p,const pthread_mutexattr_t *a) {
  if(++initialization==initialization_fail) return EAGAIN;
  return pthread_mutex_init(p,a);
}
static int fault_cond_init(pthread_cond_t *p,const pthread_condattr_t *a) {
  if(++initialization==initialization_fail) return EAGAIN;
  return pthread_cond_init(p,a);
}
static int fault_create(pthread_t *p,const pthread_attr_t *a,
    void *(*work)(void *),void *job) {
  if(++creation==creation_fail) return EAGAIN;
  return pthread_create(p,a,work,job);
}
#define calloc fault_calloc
#define free fault_free
#define pthread_mutex_init fault_mutex_init
#define pthread_cond_init fault_cond_init
#define pthread_create fault_create
#include "../src/kls_componentwise_team.inc"
#undef calloc
#undef free
#undef pthread_mutex_init
#undef pthread_cond_init
#undef pthread_create

typedef struct {
  int count,hits[8];
  kls_componentwise_team *team;
} job;
static void count_member(void *opaque,int tid) {
  job *j=opaque;assert(tid>=0 && tid<j->count);++j->hits[tid];
  /* A reentrant dispatch must fail even if the other members finish first. */
  if(tid==0) {
    pthread_mutex_lock(&j->team->lock);
    while(j->team->pending) pthread_cond_wait(&j->team->done_cond,&j->team->lock);
    pthread_mutex_unlock(&j->team->lock);
    assert(!kls_componentwise_team_execute(j->team,count_member,j,j->count));
  }
}
static void exercise(int count) {
  kls_componentwise_team *team=kls_componentwise_team_create(count);assert(team);
  job j={.count=count,.team=team};
  assert(!kls_componentwise_team_execute(team,count_member,&j,count+1));
  assert(!kls_componentwise_team_execute(team,NULL,&j,count));
  for(int repeat=0;repeat<64;++repeat) {
    assert(kls_componentwise_team_execute(team,count_member,&j,count));
    assert(!team->pending && !team->work && !team->job);
    for(int tid=0;tid<count;++tid) assert(j.hits[tid]==repeat+1);
  }
  kls_componentwise_team_destroy(team);assert(!live);
}
int main(void) {
  for(int fail=1;fail<=3;++fail) {
    allocation=initialization=creation=0;allocation_fail=fail;
    assert(!kls_componentwise_team_create(4));assert(!live);
  }
  allocation_fail=0;
  for(int fail=1;fail<=3;++fail) {
    allocation=initialization=creation=0;initialization_fail=fail;
    assert(!kls_componentwise_team_create(4));assert(!live);
  }
  initialization_fail=0;
  for(int fail=1;fail<=3;++fail) {
    allocation=initialization=creation=0;creation_fail=fail;
    assert(!kls_componentwise_team_create(4));assert(!live);
  }
  creation_fail=0;
  assert(!kls_componentwise_team_create(0));
  assert(!kls_componentwise_team_create(1));
  assert(!kls_componentwise_team_execute(NULL,count_member,NULL,2));
  exercise(2);exercise(3);exercise(8);
  puts("PASS sleeping certificate team: allocation/init/thread faults, repeated jobs, reentrancy, joined cleanup");
}
