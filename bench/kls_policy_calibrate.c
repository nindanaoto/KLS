/* Synthetic-only calibration of AUTO setup and lifecycle costs. */
#define _POSIX_C_SOURCE 200809L
#include "kls/kls.h"
#include "kls_tuning_profile.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

typedef struct {int64_t n,nz,*p,*i;uint64_t seed;double *a,*base,*b,*x;} policy_matrix;
typedef struct {double x[KLS_POLICY_COST_COUNT],y;} policy_sample;
static double policy_clock(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static int policy_index_cmp(const void *a,const void *b){int64_t x=*(const int64_t*)a,y=*(const int64_t*)b;return (x>y)-(x<y);}
static void policy_free(policy_matrix *m){free(m->p);free(m->i);free(m->a);free(m->base);free(m->b);free(m->x);memset(m,0,sizeof *m);}
static int policy_generate(policy_matrix *m,int n,int family,uint64_t seed){
  memset(m,0,sizeof *m);m->n=n;m->seed=seed;
  m->p=calloc((size_t)n+1,sizeof *m->p);m->i=malloc((size_t)n*20*sizeof *m->i);
  m->a=malloc((size_t)n*20*sizeof *m->a);m->base=malloc((size_t)n*20*sizeof *m->base);m->b=calloc((size_t)n,sizeof *m->b);m->x=malloc((size_t)n*sizeof *m->x);
  if(!m->p||!m->i||!m->a||!m->base||!m->b||!m->x){policy_free(m);return 0;}
  int side=(int)sqrt(n);double *rows=calloc((size_t)n,sizeof *rows);if(!rows){policy_free(m);return 0;}
  for(int c=0;c<n;++c){
    int64_t v[20];int count=0;v[count++]=c;
    if(family==0){
      if(c%side)v[count++]=c-1;
      if(c%side!=side-1 && c+1<n)v[count++]=c+1;
      if(c>=side)v[count++]=c-side;
      if(c+side<n)v[count++]=c+side;
    }else if(family==1){
      for(int j=1;j<=3+(n%7);++j)v[count++]=(c*(2*j+1)+j+(int)(seed%29))%n;
    }else{
      int first=(c/32)*32,last=first+32<n?first+32:n;
      for(int j=1;j<=2;++j){if(c-j>=first)v[count++]=c-j;if(c+j<last)v[count++]=c+j;}
      if(c%32==0 && c+32<n)v[count++]=c+32;
    }
    qsort(v,(size_t)count,sizeof *v,policy_index_cmp);
    m->p[c]=m->nz;
    for(int j=0;j<count;++j){if(j&&v[j]==v[j-1])continue;
      int64_t r=v[j];m->i[m->nz]=r;m->a[m->nz]=r==c?1:-.01*(1+(r+c)%5);rows[r]+=r==c?0:fabs(m->a[m->nz]);++m->nz;}
  }
  m->p[n]=m->nz;
  for(int c=0;c<n;++c)for(int64_t p=m->p[c];p<m->p[c+1];++p){int r=(int)m->i[p];if(r==c)m->a[p]=1+2*rows[r];m->b[r]+=m->a[p]*(.25*((c+seed)%11)-1.25);}
  memcpy(m->base,m->a,(size_t)m->nz*sizeof *m->a);free(rows);return 1;
}
static int policy_measure(policy_matrix *m,int threads,kls_ordering order,
    const char *profile,int cycles,policy_sample *samples,int *count,double *total){
  kls_options o;kls_default_options(&o);o.threads=threads;o.ordering=order;
  o.orientation=KLS_ORIENTATION_NORMAL;o.expected_refactorizations=cycles-1;o.expected_solves=cycles;o.tuning_profile_path=profile;
  memcpy(m->a,m->base,(size_t)m->nz*sizeof *m->a);
  memset(m->b,0,(size_t)m->n*sizeof *m->b);
  for(int64_t c=0;c<m->n;++c)for(int64_t p=m->p[c];p<m->p[c+1];++p)m->b[m->i[p]]+=m->a[p]*(.25*((c+m->seed)%11)-1.25);
  kls_solver *s=NULL;double start=policy_clock();int status=kls_create(&s);
  if(status==KLS_OK)status=kls_analyze_csc(s,KLS_INDEX_INT64,m->n,m->p,m->i,0,&o);
  double analysis=policy_clock()-start;
  double begin=policy_clock();if(status==KLS_OK)status=kls_factor(s,m->a);double factor=policy_clock()-begin;
  double solve=0,warm_refactors=0,first_refactor=0;*total=analysis+factor;
  for(int k=0;k<cycles&&status==KLS_OK;++k){
    if(k){
      memset(m->b,0,(size_t)m->n*sizeof *m->b);
      for(int64_t c=0;c<m->n;++c)for(int64_t p=m->p[c];p<m->p[c+1];++p){
        m->a[p]=m->base[p]*(1+1e-4*(double)((p*17+k*13)%29-14));
        m->b[m->i[p]]+=m->a[p]*(.25*((c+m->seed)%11)-1.25);
      }
      begin=policy_clock();status=kls_refactor(s,m->a);
      double ref_seconds=policy_clock()-begin;*total+=ref_seconds;
      if(k==1)first_refactor=ref_seconds;else warm_refactors+=ref_seconds;
    }
    begin=policy_clock();if(status==KLS_OK)status=kls_solve(s,1,m->b,m->n,m->x,m->n);double duration=policy_clock()-begin;solve+=duration;*total+=duration;
    for(int64_t r=0;r<m->n&&status==KLS_OK;++r)if(!isfinite(m->x[r]) || fabs(m->x[r]-(.25*((r+m->seed)%11)-1.25))>1e-9)status=KLS_ERR_SOLVE_FAILED;
  }
  if(samples&&status==KLS_OK){
    kls_stats st={.struct_size=sizeof st};status=kls_get_stats(s,&st);
    double pairs=0;for(int64_t c=0;c<m->n;++c){double d=(double)(m->p[c+1]-m->p[c]);pairs+=d*d;}
    policy_sample *a=&samples[(*count)++];memset(a,0,sizeof *a);a->x[KLS_POLICY_SCAN]=(double)(m->n+m->nz);a->x[KLS_POLICY_PAIR]=pairs;a->y=analysis;
    a=&samples[(*count)++];memset(a,0,sizeof *a);
    a->x[KLS_POLICY_FACTOR_FLOP]=fmax(0,st.factor_flops);
    a->x[KLS_POLICY_FACTOR_ENTRY]=(double)(st.nnz_l+st.nnz_u);
    a->x[KLS_POLICY_FACTOR_DISPATCH]=(double)st.nblocks;
    a->x[KLS_POLICY_COPY_ENTRY]=(double)m->nz;a->y=factor;
    a=&samples[(*count)++];memset(a,0,sizeof *a);
    a->x[KLS_POLICY_FLOP]=fmax(0,st.factor_flops);
    a->x[KLS_POLICY_ENTRY]=(double)(st.nnz_l+st.nnz_u);
    a->x[KLS_POLICY_DISPATCH]=(double)st.nblocks;
    a->x[KLS_POLICY_COPY_ENTRY]=(double)m->nz;
    a->y=cycles>2?warm_refactors/(cycles-2):first_refactor;
    policy_sample steady=*a;
    a=&samples[(*count)++];*a=steady;
    a->x[KLS_POLICY_FIRST_REFACTOR_ENTRY]=(double)(st.nnz_l+st.nnz_u);
    a->y=first_refactor;
    a=&samples[(*count)++];memset(a,0,sizeof *a);a->x[KLS_POLICY_SOLVE_ENTRY]=(double)(st.nnz_l+st.nnz_u);a->x[KLS_POLICY_SCAN]=(double)m->nz;a->x[KLS_POLICY_DISPATCH]=(double)st.nblocks;a->y=solve/cycles;
  }
  kls_destroy(s);return status==KLS_OK;
}
static double policy_median(double *a,int n){for(int i=1;i<n;++i){double v=a[i];int j=i;while(j&&a[j-1]>v){a[j]=a[j-1];--j;}a[j]=v;}return n%2?a[n/2]:.5*(a[n/2-1]+a[n/2]);}
/* Fixed-count projected coordinate descent; all coefficients are nonnegative. */
static void policy_fit(policy_sample *a,int n,kls_policy_cost_model *m){
  double gram[KLS_POLICY_COST_COUNT][KLS_POLICY_COST_COUNT]={{0}},rhs[KLS_POLICY_COST_COUNT]={0},scale[KLS_POLICY_COST_COUNT],v[KLS_POLICY_COST_COUNT]={0};
  for(int r=0;r<n;++r)if(a[r].y>0&&isfinite(a[r].y))for(int i=0;i<KLS_POLICY_COST_COUNT;++i){double x=a[r].x[i]/a[r].y;rhs[i]+=x;for(int j=0;j<KLS_POLICY_COST_COUNT;++j)gram[i][j]+=x*a[r].x[j]/a[r].y;}
  for(int i=0;i<KLS_POLICY_COST_COUNT;++i)scale[i]=fmax(1,sqrt(gram[i][i]));
  for(int i=0;i<KLS_POLICY_COST_COUNT;++i){rhs[i]/=scale[i];for(int j=0;j<KLS_POLICY_COST_COUNT;++j)gram[i][j]/=scale[i]*scale[j];}
  for(int it=0;it<2000;++it)for(int i=0;i<KLS_POLICY_COST_COUNT;++i)if(gram[i][i]>0){double r=rhs[i];for(int j=0;j<KLS_POLICY_COST_COUNT;++j)r-=gram[i][j]*v[j];v[i]=fmax(0,v[i]+r/gram[i][i]);}
  for(int i=0;i<KLS_POLICY_COST_COUNT;++i)m->seconds[i]=fmax(0,v[i]/scale[i]);
}
int kls_policy_calibrate_model(const char *output,int threads,uint64_t seed,int force,int quick){
  policy_sample samples[48];int count=0;double elapsed;
  for(int f=0;f<3;++f)for(int size=0;size<3;++size){
    policy_matrix a;int n=(quick?32:128)<<(2*size);
    if(!policy_generate(&a,n,f,seed))return 0;
    int ok=policy_measure(&a,threads,KLS_ORDERING_AMD,NULL,33,samples,&count,&elapsed);policy_free(&a);if(!ok)return 0;
  }
  kls_policy_cost_model model={.enabled=1,.cache_bytes=kls_tuning_policy_cache_bytes()};
  /* Sorting and copies are calibrated independently of solver arithmetic. */
  enum{N=8192};int64_t *a=malloc(N*sizeof *a),*b=malloc(N*sizeof *b);if(!a||!b){free(a);free(b);return 0;}
  double copy[5],sort[5];volatile int64_t checksum=0;
  for(int k=0;k<5;++k){for(int j=0;j<N;++j)a[j]=(j*7919+(int)(seed%97))%N;
    double t=policy_clock();for(int r=0;r<64;++r){memcpy(b,a,N*sizeof *b);
#if defined(__GNUC__) || defined(__clang__)
      __asm__ __volatile__("" : : "r"(b) : "memory");
#endif
      checksum+=b[r];}copy[k]=(policy_clock()-t)/(64*N);
    t=policy_clock();qsort(a,N,sizeof *a,policy_index_cmp);sort[k]=(policy_clock()-t)/(N*log2(N));}
  policy_sample *row=&samples[count++];memset(row,0,sizeof *row);
  row->x[KLS_POLICY_COPY_ENTRY]=64*N;row->y=policy_median(copy,5)*64*N;
  row=&samples[count++];memset(row,0,sizeof *row);
  row->x[KLS_POLICY_SORT]=N*log2(N);row->y=policy_median(sort,5)*N*log2(N);
  policy_fit(samples,count,&model);free(a);free(b);(void)checksum;
  char path[4096];if(snprintf(path,sizeof path,"%s.candidate.%ld",output,(long)getpid())>=(int)sizeof path)return 0;
  if(kls_tuning_write_policy_profile(path,threads,&model,seed,0)!=KLS_OK){
    fprintf(stderr,"policy calibration declined: invalid or incomplete cost fit; no profile published\n");return 0;
  }
  double logs=0;int ok=1,accepted=1;
  for(int f=0;f<3&&ok;++f){policy_matrix a;if(!policy_generate(&a,quick?96:768,f,seed+1)){ok=0;break;}
    double ratios[8];
    for(int pass=0;pass<8&&ok;++pass){double t[2];for(int j=0;j<2&&ok;++j){int side=(pass+j)%2;ok=policy_measure(&a,threads,KLS_ORDERING_AUTO,side?path:NULL,quick?2:100,NULL,NULL,&t[side]);}if(ok)ratios[pass]=t[1]/t[0];}
    if(ok){double ratio=policy_median(ratios,8);fprintf(stderr,"policy synthetic validation family %d ratio %.6g\n",f,ratio);accepted &= isfinite(ratio)&&ratio>0&&ratio<=1.02;logs+=log(ratio);}policy_free(&a);
  }
  unlink(path);
  if(!ok || !accepted || exp(logs/3)>1){fprintf(stderr,"policy calibration declined: synthetic lifecycle gate failed; no profile published\n");return 0;}
  return kls_tuning_write_policy_profile(output,threads,&model,seed,force)==KLS_OK;
}
