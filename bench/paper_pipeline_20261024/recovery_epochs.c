#define _GNU_SOURCE
#define _POSIX_C_SOURCE 200809L
#include KLS_CORE_SOURCE
#include <sys/resource.h>
static double elapsed(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
int main(int argc,char**argv){
 if(argc!=3)return 2;FILE*f=fopen(argv[1],"rb");uint64_t h[4];kls_options o;
 if(!f||fread(h,sizeof(h),1,f)!=1||fread(&o,sizeof(o),1,f)!=1||h[0]!=UINT64_C(0x4b4c534449414731)||h[3]!=1)return 3;
 int64_t n=h[1],z=h[2];int family=atoi(argv[2]);
 int64_t*p=malloc((n+1)*8),*j=malloc(z*8);double*a=malloc(z*8),*b=malloc(n*8),*v=malloc(z*8),*rhs=malloc(n*8),*x=calloc(n,8);
 if(fread(p,8,n+1,f)!=(size_t)n+1||fread(j,8,z,f)!=(size_t)z||fread(a,8,z,f)!=(size_t)z||fread(b,8,n,f)!=(size_t)n)return 4;fclose(f);
 kls_solver*s=NULL;kls_create(&s);int st=kls_analyze_csr(s,KLS_INDEX_INT64,n,p,j,0,&o);if(!st)st=kls_factor(s,a);if(st)return 5;
 double total=0,cold=0;int valid=0,updates=8;long double worst=0;
 for(int epoch=0;epoch<=updates;++epoch){
  for(int64_t r=0;r<n;++r){int power=epoch?((r+epoch)%3)-1:0;long double sum=0;
   rhs[r]=family==0?ldexp(b[r],power):b[r];
   for(int64_t q=p[r];q<p[r+1];++q){
    v[q]=family==0?ldexp(a[q],power):family==1?ldexp(a[q],epoch?((j[q]+epoch)%3)-1:0):a[q]*(epoch?1+1e-4*((q+epoch)%7-3):1);
    if(family==2 && epoch)sum+=(long double)v[q]*(1+(j[q]%5)*.125);
   }if(family==2&&epoch)rhs[r]=(double)sum;
  }
  double start=elapsed();
  if(!epoch)st=kls_solve(s,1,rhs,0,x,0);
  else if(epoch%2)st=kls_refactor_solve(s,v,1,rhs,0,x,0);
  else{st=kls_refactor(s,v);if(!st)st=kls_solve(s,1,rhs,0,x,0);}
  double seconds=elapsed()-start;if(epoch)total+=seconds;else cold=seconds;
  long double bound=INFINITY;int cert=kls_componentwise_certify(n,p,j,v,rhs,x,1,NULL,&bound);
  if(st==KLS_OK&&cert==1)++valid;if(isfinite(bound)&&bound>worst)worst=bound;
  printf("epoch index=%d status=%d certificate=%d seconds=%.17g bound=%.21Lg\n",epoch,st,cert,seconds,bound);
 }
 struct rusage ru;getrusage(RUSAGE_SELF,&ru);uint64_t attempts=0,accepted=0,builds=0;
#ifdef PAPER_RECOVERY_REUSE
 attempts=s->private_scaled_recovery_updates;accepted=s->private_scaled_recovery_accepts;builds=s->private_scaled_recovery_builds;
#endif
 printf("{\"n\":%ld,\"nnz\":%ld,\"family\":%d,\"valid\":%d,\"epochs\":%d,\"cold_seconds\":%.17g,\"update_seconds\":%.17g,\"mean_update_seconds\":%.17g,\"peak_kb\":%ld,\"reuse_attempts\":%lu,\"reuse_accepts\":%lu,\"fresh_builds\":%lu,\"worst_bound\":%.21Lg}\n",n,z,family,valid,updates+1,cold,total,total/updates,ru.ru_maxrss,attempts,accepted,builds,worst);
 kls_destroy(s);free(p);free(j);free(a);free(b);free(v);free(rhs);free(x);return 0;
}
