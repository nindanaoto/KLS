/* Private C orchestration for bounded local lattice recovery.
 * No matrix identity, stored answer or external optimization package is used.
 * Return 1 only for a certified vector, 0 for an exhausted search, -1 for
 * invalid input/workspace failure. Failed setup preserves the caller's x. */
#define _POSIX_C_SOURCE 200809L
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <time.h>
#include "kls_tuning.inc"

extern int kls_accuracy_lattice_propose(int,int,const double*,const double*,double*,int);
extern int kls_accuracy_certify(size_t,const int64_t*,const int64_t*,const double*,
                        const double*,const double*,int,double*);
typedef struct { long double score;int64_t index; } ranked;
static int ascending(const void *aa,const void *bb) {
  const ranked *a=aa,*b=bb;
  if(a->score<b->score) return -1;
  if(a->score>b->score) return 1;
  return (a->index>b->index)-(a->index<b->index);
}
static int descending(const void *a,const void *b) { return -ascending(a,b); }
static double seconds_now(void) {
  struct timespec ts;if(clock_gettime(CLOCK_MONOTONIC,&ts)) return INFINITY;
  return (double)ts.tv_sec+1e-9*(double)ts.tv_nsec;
}
static long double input_ld(double x) { volatile long double value=x;return value; }
static long double spacing_ld(double x) {
  long double value=fabsl(input_ld(x));int exponent=0;
  if(value>0) (void)frexpl(value,&exponent);
  exponent=value>0?exponent-DBL_MANT_DIG:DBL_MIN_EXP-DBL_MANT_DIG;
  if(exponent<DBL_MIN_EXP-DBL_MANT_DIG) exponent=DBL_MIN_EXP-DBL_MANT_DIG;
  return ldexpl(1.L,exponent);
}
/* Offline candidate: permit useful moves in tiny/zero components without
 * enumerating astronomical numbers of their native ULPs. The artificial
 * power-of-two grid is no coarser than an equal share of the existing error
 * budget divided by the column norm. Only the final certificate accepts x. */
static long double search_spacing(double x,long double norm,long double budget,int n) {
  long double step=spacing_ld(x);
  if(norm>0) {
    int exponent=0;
    long double allowance=budget/((long double)n*norm);
    if(!(allowance>0) || !isfinite(allowance)) return step;
    (void)frexpl(allowance,&exponent);
    long double grid=ldexpl(1.L,exponent-1);
    if(grid>step) step=grid;
  }
  return step;
}
static void residual_csc(int n,const int64_t *p,const int64_t *rows,const double *a,
                         const double *b,const double *x,long double budget,long double *r) {
  for(int j=0;j<n;++j) r[j]=input_ld(b[j]);
  for(int j=0;j<n;++j) {
    long double value=input_ld(x[j]);
    for(int64_t k=p[j];k<p[j+1];++k) r[rows[k]]-=input_ld(a[k])*value;
  }
  for(int j=0;j<n;++j) r[j]/=budget;
}

static int lattice_recover(int n,const int64_t *p,const int64_t *rows,const double *a,
                    const double *b,double *x,int width,int sweeps,int row_limit,
                    double regularizer,double seconds,int *proposals_out,
                    int *no_progress_out) {
  if(no_progress_out) *no_progress_out=0;
  if(n<=0 || n>INT_MAX-64 || width<=0 || width>64 || sweeps<=0 || row_limit<=0 ||
     !p || !rows || !a || !b || !x || p[0]!=0 || p[n]<0 ||
     !isfinite(regularizer) || regularizer<=0 || !isfinite(seconds) || seconds<=0 ||
     (size_t)n>SIZE_MAX/64/sizeof(double)-64 ||
     (uint64_t)p[n]>SIZE_MAX/sizeof(int64_t)) return -1;
  const size_t nz=(size_t)p[n],count=(size_t)n;
  /* Setup and the initial certificate consume this sweep's remaining
   * allowance too; widening must not hide repeated setup from the deadline. */
  const double began=seconds_now();
  if(!isfinite(began)) return -1;
  int status=-1,proposals=0;
  int64_t *rp=NULL,*rc=NULL,*cursor=NULL,*map=NULL,*touched=NULL;
  double *rv=NULL,*work=NULL,*best=NULL,*dense=NULL,*target=NULL;
  long double *res=NULL,*trial=NULL,*norms=NULL;
  unsigned char *marked=NULL,*rowmarked=NULL;
  ranked *ranking=NULL,*candidates=NULL;
  rp=calloc(count+1,sizeof(*rp));rc=malloc((nz?nz:1)*sizeof(*rc));
  cursor=malloc(count*sizeof(*cursor));rv=malloc((nz?nz:1)*sizeof(*rv));
  map=malloc(count*sizeof(*map));touched=malloc(count*sizeof(*touched));
  work=malloc(count*sizeof(*work));best=malloc(count*sizeof(*best));
  dense=malloc((count+(size_t)width)*(size_t)width*sizeof(*dense));
  target=malloc((count+(size_t)width)*sizeof(*target));
  res=malloc(count*sizeof(*res));trial=malloc(count*sizeof(*trial));
  norms=calloc(count,sizeof(*norms));marked=malloc(count);rowmarked=malloc(count);
  ranking=malloc(count*sizeof(*ranking));candidates=malloc(count*sizeof(*candidates));
  if(!rp||!rc||!cursor||!rv||!map||!touched||!work||!best||!dense||!target||
     !res||!trial||!norms||!marked||!rowmarked||!ranking||!candidates) goto done;
  long double bn=0;
  for(int j=0;j<n;++j) {
    if(p[j]<0 || p[j]>p[j+1] || p[j+1]>p[n] ||
       !isfinite(input_ld(b[j])) || !isfinite(input_ld(x[j]))) goto done;
    long double value=input_ld(b[j]);bn+=value*value;
    for(int64_t k=p[j];k<p[j+1];++k) {
      if(rows[k]<0 || rows[k]>=n || !isfinite(input_ld(a[k]))) goto done;
      ++rp[rows[k]+1];value=input_ld(a[k]);norms[j]+=value*value;
    }
    norms[j]=sqrtl(norms[j]);
  }
  for(int j=0;j<n;++j) { rp[j+1]+=rp[j];cursor[j]=rp[j]; }
  for(int j=0;j<n;++j) for(int64_t k=p[j];k<p[j+1];++k) {
    int64_t dest=cursor[rows[k]]++;rc[dest]=j;rv[dest]=a[k];
  }
  memcpy(work,x,count*sizeof(*x));memcpy(best,x,count*sizeof(*x));
  long double budget=1e-8L*(bn>0?sqrtl(bn):1);
  double best_ratio=INFINITY;
  status=kls_accuracy_certify(count,p,rows,a,b,best,0,&best_ratio)?1:0;
  if(status) goto finished;
  for(int sweep=0;sweep<sweeps;++sweep) {
    residual_csc(n,p,rows,a,b,work,budget,res);
    for(int j=0;j<n;++j) ranking[j]=(ranked){fabsl(res[j]),j};
    qsort(ranking,count,sizeof(*ranking),descending);
    int changed=0,completed=1,limit=row_limit<n?row_limit:n;
    for(int visit=0;visit<limit;++visit) {
      if(seconds_now()-began>seconds) { completed=0;break; }
      int64_t row=ranking[visit].index,cols[64];int nc=0,np=0;
      memset(marked,0,count);memset(rowmarked,0,count);
      for(int64_t k=rp[row];k<rp[row+1];++k) {
        int64_t col=rc[k];if(marked[col]) continue;
        marked[col]=1;candidates[np++]=(ranked){fabsl(input_ld(rv[k])*search_spacing(work[col],norms[col],budget,n)),col};
      }
      qsort(candidates,(size_t)np,sizeof(*candidates),ascending);
      memset(marked,0,count);
      for(int j=np>width?np-width:0;j<np;++j) { cols[nc++]=candidates[j].index;marked[candidates[j].index]=1; }
      if(!nc) continue;
      if(nc<width) {
        for(int j=0;j<nc;++j) for(int64_t k=p[cols[j]];k<p[cols[j]+1];++k) rowmarked[rows[k]]=1;
        np=0;
        for(int r=0;r<n;++r) if(rowmarked[r]) for(int64_t k=rp[r];k<rp[r+1];++k) {
          int64_t col=rc[k];if(marked[col]) continue;
          marked[col]=2;long double dot=0;
          for(int64_t t=p[col];t<p[col+1];++t) dot+=input_ld(a[t])*res[rows[t]];
          candidates[np++]=(ranked){norms[col]>0?fabsl(dot)/norms[col]:0,col};
        }
        qsort(candidates,(size_t)np,sizeof(*candidates),ascending);
        int slots=width-nc;
        for(int j=np>slots?np-slots:0;j<np;++j) cols[nc++]=candidates[j].index;
      }
      memset(rowmarked,0,count);
      for(int j=0;j<nc;++j) for(int64_t k=p[cols[j]];k<p[cols[j]+1];++k) rowmarked[rows[k]]=1;
      int nr=0;
      for(int r=0;r<n;++r) if(rowmarked[r]) { map[r]=nr;touched[nr++]=r; }
      int m=nr+nc;
      memset(dense,0,(size_t)m*nc*sizeof(*dense));memset(target,0,(size_t)m*sizeof(*target));
      long double scales[64];double delta[64],replacement[64];
      for(int r=0;r<nr;++r) target[r]=(double)res[touched[r]];
      for(int j=0;j<nc;++j) {
        scales[j]=search_spacing(work[cols[j]],norms[cols[j]],budget,n);
        for(int64_t k=p[cols[j]];k<p[cols[j]+1];++k)
          dense[(size_t)map[rows[k]]*nc+j]+=(double)(input_ld(a[k])*scales[j]/budget);
        dense[(size_t)(nr+j)*nc+j]=regularizer;
      }
      int proposed=kls_accuracy_lattice_propose(m,nc,dense,target,delta,KLS_COLD_LATTICE_MAX_STEPS);
      if(proposed<0) { completed=0;continue; }
      if(!proposed) continue;
      ++proposals;memcpy(trial,res,count*sizeof(*trial));int finite=1;
      for(int j=0;j<nc;++j) {
        replacement[j]=(double)(input_ld(work[cols[j]])+scales[j]*input_ld(delta[j]));
        if(!isfinite(input_ld(replacement[j]))) { finite=0;break; }
        long double change=input_ld(replacement[j])-input_ld(work[cols[j]]);
        for(int64_t k=p[cols[j]];k<p[cols[j]+1];++k) trial[rows[k]]-=input_ld(a[k])*change/budget;
      }
      if(!finite) continue;
      long double before=0,after=0;
      for(int j=0;j<n;++j) { before+=res[j]*res[j];after+=trial[j]*trial[j]; }
      if(after<before) {
        for(int j=0;j<nc;++j) work[cols[j]]=replacement[j];
        memcpy(res,trial,count*sizeof(*res));++changed;
        /* The cheap residual only triggers a certificate attempt. It cannot
         * accept an answer. Avoid waiting until the entire sweep/deadline
         * when the current representable vector already meets the contract. */
        if(after<=1) {
          double ratio=INFINITY;
          if(kls_accuracy_certify(count,p,rows,a,b,work,0,&ratio)) {
            memcpy(best,work,count*sizeof(*best));status=1;goto finished;
          }
        }
      }
    }
    double ratio=INFINITY;int passed=kls_accuracy_certify(count,p,rows,a,b,work,0,&ratio);
    if(ratio<best_ratio || passed) { memcpy(best,work,count*sizeof(*best));best_ratio=ratio; }
    if(passed) { status=1;break; }
    /* Only a fully visited sweep with no accepted internal changes is a
     * fixed point. An unchanged published best vector alone is insufficient:
     * work may have changed, or a deadline may have truncated the search. */
    if(no_progress_out && completed && !changed) *no_progress_out=1;
    if(!changed || seconds_now()-began>seconds) break;
  }
finished:
  memcpy(x,best,count*sizeof(*x));
done:
  if(proposals_out) *proposals_out=proposals;
  free(rp);free(rc);free(cursor);free(map);free(touched);free(rv);free(work);free(best);
  free(dense);free(target);free(res);free(trial);free(norms);free(marked);free(rowmarked);
  free(ranking);free(candidates);return status;
}

int kls_accuracy_lattice_recover(int n,const int64_t *p,const int64_t *rows,
    const double *a,const double *b,double *x,int width,int sweeps,int row_limit,
    double regularizer,double seconds,int *proposals_out) {
  return lattice_recover(n,p,rows,a,b,x,width,sweeps,row_limit,regularizer,
                         seconds,proposals_out,NULL);
}

/* The scheduler needs a stronger signal than equality of the best answer.
 * Keep the multi-sweep entry point's existing result convention unchanged. */
int kls_accuracy_lattice_sweep(int n,const int64_t *p,const int64_t *rows,
    const double *a,const double *b,double *x,int width,int row_limit,
    double regularizer,double seconds,int *proposals_out,int *no_progress_out) {
  return lattice_recover(n,p,rows,a,b,x,width,1,row_limit,regularizer,
                         seconds,proposals_out,no_progress_out);
}
