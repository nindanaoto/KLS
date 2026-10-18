#define _GNU_SOURCE
#include <dlfcn.h>
#include <errno.h>
#include <limits.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Diagnostic only. All original arguments and return codes are delegated
   unchanged. File/registry failures disable recording, never change the API.
   Integer/byte operations preserve coefficients and floating-point flags.
   Public calls on one instance must not overlap, as in the Xyce adapter.
   Separate instances may be recorded concurrently under the file mutex. */
#define ANALYZE "_ZN10subtree_lu9SubtreeLUIidE7analyzeEiPKiS3_PKdiPiS6_"
#define FACTOR "_ZN10subtree_lu9SubtreeLUIidE9factorizeEPKd"
#define REFACTOR "_ZN10subtree_lu9SubtreeLUIidE11refactorizeEPKd"
#define SOLVE "_ZN10subtree_lu9SubtreeLUIidE5solveEPKdPd"
#define CSR_MAGIC UINT64_C(0x53544c5543535231)
#define STREAM_MAGIC UINT64_C(0x53544c5553545231)
_Static_assert(sizeof(int)==4 && sizeof(double)==8,"Native int32/binary64 ABI");
enum {ST_ANALYZE=1,ST_FACTOR=2,ST_REFACTOR=3,ST_SOLVE=4};
enum {ENTER=1,LEAVE=2,NONE=0,DOUBLES=1,INTS=2};
typedef struct context {
 void *self;int n,nnz;uint64_t generation,sequence;FILE *stream;int failed;
 struct context *next;
} context;
typedef struct {context *c;uint64_t generation,sequence,flags;int n,nnz;} ticket;
static context *contexts;
static pthread_mutex_t lock=PTHREAD_MUTEX_INITIALIZER;
static uint64_t generation_counter;
static void event(unsigned stage,unsigned edge,int status,uint64_t sequence) {
 int saved=errno;
 fprintf(stderr,"SUBTREELU_STREAM stage=%u edge=%u status=%d sequence=%llu\n",
   stage,edge,status,(unsigned long long)sequence);fflush(stderr);errno=saved;
}
static void recording_error(context *c,const char *what) {
 if(c)c->failed=1;
 fprintf(stderr,"SUBTREELU_STREAM_RECORDING_ERROR %s\n",what);fflush(stderr);
}
static void *resolve(const char *symbol) {
 int saved=errno;void *p=dlsym(RTLD_NEXT,symbol);
 if(!p){fprintf(stderr,"SubtreeLU callstream cannot delegate original symbol: %s\n",dlerror());_Exit(125);}
 errno=saved;return p;
}
static context *lookup(void *self,int create) {
 for(context *c=contexts;c;c=c->next)if(c->self==self)return c;
 if(!create)return NULL;
 context *c=calloc(1,sizeof(*c));if(!c){recording_error(NULL,"registry allocation");return NULL;}
 c->self=self;c->next=contexts;contexts=c;return c;
}
static int bytes(FILE *f,const void *p,size_t width,size_t count) {
 return count==0 || (p && fwrite(p,width,count,f)==count);
}
static void record_locked(ticket t,unsigned stage,unsigned edge,int status,
    unsigned kind,const void *a,size_t na,const void *b,size_t nb) {
 context *c=t.c;
 if(!c || c->generation!=t.generation || !c->stream || c->failed)return;
 uint64_t header[8]={stage,edge,t.sequence,(uint64_t)(int64_t)status,
   kind,(uint64_t)na+(uint64_t)nb,t.flags,0};
 size_t width=kind==INTS?sizeof(int):sizeof(double);
 if(!bytes(c->stream,header,sizeof(header),1) ||
    !bytes(c->stream,a,width,na) || !bytes(c->stream,b,width,nb) ||
    fflush(c->stream)!=0)recording_error(c,"call record write/flush");
}
static int capture_initial(const char *folder,int n,int nnz,const int *ap,
    const int *ai,const double *ax,int threads,int *rp,int *cp) {
 char path[PATH_MAX];int length=snprintf(path,sizeof(path),"%s/initial.csr",folder);
 if(length<0 || (size_t)length>=sizeof(path))return 0;
 FILE *f=fopen(path,"wb");if(!f)return 0;
 const uint64_t flags=(ax?1u:0u)|(rp?2u:0u)|(cp?4u:0u)|(ai?8u:0u)|16u;
 uint64_t header[8]={CSR_MAGIC,2,(uint64_t)n,(uint64_t)nnz,(uint64_t)threads,flags,0,0};
 int ok=bytes(f,header,sizeof(header),1) && bytes(f,ap,sizeof(int),(size_t)n+1) &&
   bytes(f,ai,sizeof(int),(size_t)nnz) && (!ax || bytes(f,ax,sizeof(double),(size_t)nnz)) &&
   (!rp || bytes(f,rp,sizeof(int),(size_t)n)) && (!cp || bytes(f,cp,sizeof(int),(size_t)n));
 if(fclose(f)!=0)ok=0;
 return ok;
}
static ticket begin_analyze(void *self,int n,const int *ap,const int *ai,
    const double *ax,int threads,int *rp,int *cp) {
 int saved=errno;ticket t={0};pthread_mutex_lock(&lock);
 context *c=lookup(self,1);
 if(c) {
  if(c->stream){fclose(c->stream);c->stream=NULL;}
  c->generation=++generation_counter;c->sequence=0;c->failed=0;c->n=n;c->nnz=0;
  t=(ticket){c,c->generation,++c->sequence,(ax?1u:0u)|(rp?2u:0u)|(cp?4u:0u),n,0};
  const char *base=getenv("KLS_SUBTREELU_CALLSTREAM_DIR");
  if(!base)recording_error(c,"missing capture directory");
  else if(n<=0 || !ap || ap[0]!=0 || ap[n]<0 || (ap[n]>0 && !ai))
   recording_error(c,"unsupported invalid CSR arguments; original API still delegated");
  else {
   c->nnz=ap[n];t.nnz=c->nnz;int valid=1;
   for(int i=0;i<n;++i)if(ap[i]<0 || ap[i]>ap[i+1]){valid=0;break;}
   for(int i=0;valid && i<c->nnz;++i)if(ai[i]<0 || ai[i]>=n){valid=0;break;}
   char folder[PATH_MAX],path[PATH_MAX];
   int length=snprintf(folder,sizeof(folder),"%s/generation-%ld-%llu-XXXXXX",base,
      (long)getpid(),(unsigned long long)c->generation);
   if(!valid || length<0 || (size_t)length>=sizeof(folder) || !mkdtemp(folder))
    recording_error(c,"CSR validation/capture directory");
   else if(!capture_initial(folder,n,c->nnz,ap,ai,ax,threads,rp,cp))
    recording_error(c,"initial CSR write");
   else {
    length=snprintf(path,sizeof(path),"%s/calls.bin",folder);
    if(length<0 || (size_t)length>=sizeof(path))recording_error(c,"stream path");
    else if(!(c->stream=fopen(path,"wb")))recording_error(c,"stream open");
    else {
     uint64_t header[8]={STREAM_MAGIC,1,(uint64_t)n,(uint64_t)c->nnz,
       (uint64_t)threads,((uint64_t)sizeof(int)<<32)|sizeof(double),c->generation,0};
     if(!bytes(c->stream,header,sizeof(header),1) || fflush(c->stream)!=0)
      recording_error(c,"stream header write/flush");
    }
   }
  }
  record_locked(t,ST_ANALYZE,ENTER,0,NONE,NULL,0,NULL,0);
 }
 pthread_mutex_unlock(&lock);event(ST_ANALYZE,ENTER,0,t.sequence);errno=saved;return t;
}
static ticket begin_numeric(void *self,unsigned stage,const double *values,
    const double *rhs,double *x) {
 int saved=errno;ticket t={0};pthread_mutex_lock(&lock);context *c=lookup(self,0);
 if(c) {
  uint64_t flags=stage==ST_SOLVE?((rhs?1u:0u)|(x?2u:0u)|((rhs && rhs==x)?4u:0u)):(values?1u:0u);
  t=(ticket){c,c->generation,++c->sequence,flags,c->n,c->nnz};
  if(stage==ST_SOLVE)record_locked(t,stage,ENTER,0,DOUBLES,rhs,rhs?(size_t)c->n:0,NULL,0);
  else record_locked(t,stage,ENTER,0,DOUBLES,values,values?(size_t)c->nnz:0,NULL,0);
 } else recording_error(NULL,"numeric call has no captured analyze generation");
 pthread_mutex_unlock(&lock);event(stage,ENTER,0,t.sequence);errno=saved;return t;
}
static void finish(ticket t,unsigned stage,int status,int *rp,int *cp,double *x) {
 int saved=errno;pthread_mutex_lock(&lock);
 if(t.c && t.c->generation==t.generation) {
  if(stage==ST_ANALYZE && status==0)
   record_locked(t,stage,LEAVE,status,INTS,rp,rp?(size_t)t.n:0,cp,cp?(size_t)t.n:0);
  else if(stage==ST_SOLVE && status==0)
   record_locked(t,stage,LEAVE,status,DOUBLES,x,x?(size_t)t.n:0,NULL,0);
  else record_locked(t,stage,LEAVE,status,NONE,NULL,0,NULL,0);
 } else if(t.c)recording_error(t.c,"overlapping/reentrant analyze generations");
 pthread_mutex_unlock(&lock);event(stage,LEAVE,status,t.sequence);errno=saved;
}
int trace_analyze(void *,int,const int *,const int *,const double *,int,int *,int *) __asm__(ANALYZE);
int trace_analyze(void *self,int n,const int *ap,const int *ai,const double *ax,int threads,int *rp,int *cp) {
 typedef int(*fn)(void *,int,const int *,const int *,const double *,int,int *,int *);
 fn next=(fn)resolve(ANALYZE);ticket t=begin_analyze(self,n,ap,ai,ax,threads,rp,cp);
 int status=next(self,n,ap,ai,ax,threads,rp,cp);finish(t,ST_ANALYZE,status,rp,cp,NULL);return status;
}
int trace_factor(void *,const double *) __asm__(FACTOR);
int trace_factor(void *self,const double *ax) {
 typedef int(*fn)(void *,const double *);fn next=(fn)resolve(FACTOR);
 ticket t=begin_numeric(self,ST_FACTOR,ax,NULL,NULL);int status=next(self,ax);
 finish(t,ST_FACTOR,status,NULL,NULL,NULL);return status;
}
int trace_refactor(void *,const double *) __asm__(REFACTOR);
int trace_refactor(void *self,const double *ax) {
 typedef int(*fn)(void *,const double *);fn next=(fn)resolve(REFACTOR);
 ticket t=begin_numeric(self,ST_REFACTOR,ax,NULL,NULL);int status=next(self,ax);
 finish(t,ST_REFACTOR,status,NULL,NULL,NULL);return status;
}
int trace_solve(void *,const double *,double *) __asm__(SOLVE);
int trace_solve(void *self,const double *b,double *x) {
 typedef int(*fn)(void *,const double *,double *);fn next=(fn)resolve(SOLVE);
 ticket t=begin_numeric(self,ST_SOLVE,NULL,b,x);int status=next(self,b,x);
 finish(t,ST_SOLVE,status,NULL,NULL,x);return status;
}
__attribute__((destructor)) static void close_streams(void) {
 pthread_mutex_lock(&lock);
 while(contexts){context *c=contexts;contexts=c->next;if(c->stream && fclose(c->stream)!=0)recording_error(c,"stream close");free(c);}
 pthread_mutex_unlock(&lock);
}
