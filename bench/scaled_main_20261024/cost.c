#define _POSIX_C_SOURCE 200809L
#include "kls/kls.h"
#include "kls_componentwise.h"
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/resource.h>
#include <math.h>
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static long rss(void){struct rusage u;getrusage(RUSAGE_SELF,&u);return u.ru_maxrss;}
static int cmp(const void*a,const void*b){double x=*(const double*)a,y=*(const double*)b;return (x>y)-(x<y);}
int main(int argc,char**argv){
 int healthy=argc==4&&!strcmp(argv[1],"healthy");if(!healthy && (argc!=3||strcmp(argv[1],"frame")))return 2;
 int64_t n,z,*p,*j;double *a,*b,*bt=NULL,*x;kls_options o;kls_default_options(&o);
 if(healthy){n=atoll(argv[2]);o.threads=atoi(argv[3]);z=3*n-2;p=malloc((n+1)*8);j=malloc(z*8);a=malloc(z*8);b=calloc(n,8);bt=calloc(n,8);int64_t at=0;
  for(int64_t r=0;r<n;++r){p[r]=at;if(r){j[at]=r-1;a[at++]=-1;}j[at]=r;a[at++]=4;if(r+1<n){j[at]=r+1;a[at++]=-.5;}}p[n]=at;
  for(int64_t r=0;r<n;++r)for(int64_t k=p[r];k<p[r+1];++k){b[r]+=a[k];bt[j[k]]+=a[k];}
 }else{FILE*f=fopen(argv[2],"rb");uint64_t h[4];if(!f||fread(h,sizeof(h),1,f)!=1||h[0]!=UINT64_C(0x4b4c534449414731)||h[3]!=1||fread(&o,sizeof(o),1,f)!=1)return 3;n=h[1];z=h[2];p=malloc((n+1)*8);j=malloc(z*8);a=malloc(z*8);b=malloc(n*8);if(fread(p,8,n+1,f)!=(size_t)n+1||fread(j,8,z,f)!=(size_t)z||fread(a,8,z,f)!=(size_t)z||fread(b,8,n,f)!=(size_t)n)return 4;fclose(f);}
 x=calloc(n,8);kls_solver*s=NULL;int st=kls_create(&s);if(!st)st=kls_analyze_csr(s,KLS_INDEX_INT64,n,p,j,0,&o);if(!st)st=kls_factor(s,a);if(st){printf("{\"setup_status\":%d}\n",st);return 0;}
 const long before=rss();double start=now();st=kls_solve(s,1,b,0,x,0);double cold=now()-start;long double bound;int cert=kls_componentwise_certify(n,p,j,a,b,x,1,NULL,&bound);long after=rss();int initial=st,initial_cert=cert;long double initial_bound=bound;
 int count=healthy?(n<=128?1024:n<=4096?256:64):3;double *times=malloc(count*sizeof(double));int successes=st==KLS_OK&&cert==1;double fresh=0;
 for(int i=0;i<count;++i){if(!healthy&&i==count-1){st=kls_refactor(s,a);if(st)return 5;}
  int tr=healthy&&(i%2);const double*rhs=tr?bt:b;start=now();st=healthy&&i%4==2?kls_refactor_solve(s,a,1,rhs,0,x,0):tr?kls_solve_transpose(s,1,rhs,0,x,0):kls_solve(s,1,rhs,0,x,0);times[i]=now()-start;
  if(!healthy||i==count-1){cert=kls_componentwise_certify(n,p,j,a,rhs,x,!tr,NULL,&bound);if(st==KLS_OK&&cert==1)++successes;}else if(st==KLS_OK)++successes;
  if(!healthy&&i==count-1)fresh=times[i];
 }
 double mean=0;for(int i=0;i<count;++i)mean+=times[i]/count;
 qsort(times,count,sizeof(double),cmp);kls_stats stats={0};stats.struct_size=sizeof(stats);kls_get_stats(s,&stats);
 printf("{\"setup_status\":0,\"n\":%ld,\"nnz\":%ld,\"threads\":%d,\"initial_status\":%d,\"initial_certificate\":%d,\"initial_bound\":%.21Lg,\"cold_seconds\":%.17g,\"mean_repeated_seconds\":%.17g,\"median_repeated_seconds\":%.17g,\"new_epoch_seconds\":%.17g,\"peak_before_kb\":%ld,\"peak_after_first_kb\":%ld,\"peak_final_kb\":%ld,\"successes\":%d,\"calls\":%d,\"incumbent_scale\":%d}\n",n,z,o.threads,initial,initial_cert,initial_bound,cold,mean,times[count/2],fresh,before,after,rss(),successes,count+1,stats.selected_scale);
 kls_destroy(s);free(p);free(j);free(a);free(b);free(bt);free(x);free(times);return 0;
}
