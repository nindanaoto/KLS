/* Optional certificate team for a solver without a numerical worker pool.
 * The caller is member zero; the remaining members sleep between jobs.
 * Dispatch, destruction and solver calls on one owner must be serialized.
 * Only the current callback/job is retained, until every member has joined. */
typedef struct kls_componentwise_team kls_componentwise_team;
typedef struct {
  kls_componentwise_team *team;
  int tid;
} kls_componentwise_team_member;
struct kls_componentwise_team {
  pthread_mutex_t lock;
  pthread_cond_t work_cond,done_cond;
  pthread_t *threads;
  kls_componentwise_team_member *members;
  int count,created,pending,shutdown;
  int lock_initialized,work_initialized,done_initialized;
  unsigned long generation;
  void (*work)(void *,int);
  void *job;
};

static void *kls_componentwise_team_worker(void *opaque) {
  kls_componentwise_team_member *member=opaque;
  kls_componentwise_team *team=member->team;
  unsigned long seen=0;
  pthread_mutex_lock(&team->lock);
  for(;;) {
    while(!team->shutdown && team->generation==seen)
      pthread_cond_wait(&team->work_cond,&team->lock);
    if(team->shutdown) break;
    seen=team->generation;
    void (*work)(void *,int)=team->work;
    void *job=team->job;
    pthread_mutex_unlock(&team->lock);
    work(job,member->tid);
    pthread_mutex_lock(&team->lock);
    if(--team->pending==0) pthread_cond_signal(&team->done_cond);
  }
  pthread_mutex_unlock(&team->lock);
  return NULL;
}

static void kls_componentwise_team_destroy(kls_componentwise_team *team) {
  if(!team) return;
  if(team->created) {
    pthread_mutex_lock(&team->lock);
    team->shutdown=1;
    pthread_cond_broadcast(&team->work_cond);
    pthread_mutex_unlock(&team->lock);
    for(int i=0;i<team->created;++i) pthread_join(team->threads[i],NULL);
  }
  if(team->done_initialized) pthread_cond_destroy(&team->done_cond);
  if(team->work_initialized) pthread_cond_destroy(&team->work_cond);
  if(team->lock_initialized) pthread_mutex_destroy(&team->lock);
  free(team->members);free(team->threads);free(team);
}

static kls_componentwise_team *kls_componentwise_team_create(int count) {
  if(count<2) return NULL;
  const size_t members=(size_t)(count-1);
  if(members>SIZE_MAX/sizeof(pthread_t) ||
     members>SIZE_MAX/sizeof(kls_componentwise_team_member)) return NULL;
  kls_componentwise_team *team=calloc(1,sizeof(*team));
  if(!team) return NULL;
  team->count=count;
  if(pthread_mutex_init(&team->lock,NULL)) goto failed;
  team->lock_initialized=1;
  if(pthread_cond_init(&team->work_cond,NULL)) goto failed;
  team->work_initialized=1;
  if(pthread_cond_init(&team->done_cond,NULL)) goto failed;
  team->done_initialized=1;
  team->threads=calloc(members,sizeof(*team->threads));
  team->members=calloc(members,sizeof(*team->members));
  if(!team->threads || !team->members) goto failed;
  for(int i=0;i<count-1;++i) {
    team->members[i].team=team;team->members[i].tid=i+1;
    if(pthread_create(team->threads+i,NULL,kls_componentwise_team_worker,
                      team->members+i)) goto failed;
    ++team->created;
  }
  return team;
failed:
  kls_componentwise_team_destroy(team);return NULL;
}

static int kls_componentwise_team_execute(kls_componentwise_team *team,
    void (*work)(void *,int),void *job,int count) {
  if(!team || !work || count!=team->count) return 0;
  pthread_mutex_lock(&team->lock);
  if(team->shutdown || team->work) {
    pthread_mutex_unlock(&team->lock);return 0;
  }
  team->work=work;team->job=job;team->pending=count-1;
  ++team->generation;
  pthread_cond_broadcast(&team->work_cond);
  pthread_mutex_unlock(&team->lock);
  work(job,0);
  pthread_mutex_lock(&team->lock);
  while(team->pending) pthread_cond_wait(&team->done_cond,&team->lock);
  team->work=NULL;team->job=NULL;
  pthread_mutex_unlock(&team->lock);
  return 1;
}
