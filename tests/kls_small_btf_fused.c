/* Application-managed small-BTF fusion: the combined call must return the
   answer of a separate refactor and solve, engage only in that mode, and
   leave no stale answer behind after a failed update. */
static int fused_events,lane_events;
#define KLS_SMALL_BTF_FUSED_TEST_EVENT() (++fused_events)
#define KLS_FAST_LANE_TEST_EVENT() (++lane_events)
#include "../src/kls.c"
#undef NDEBUG
#include <assert.h>

/* One tridiagonal block coupled to a chain of singleton blocks. */
static int build(int n,int *ptr,int *idx,double *a) {
  const int m=n-n/4;int nz=0;
  for(int i=0;i<n;++i) {
    ptr[i]=nz;
    if(i<m) {
      if(i>0){idx[nz]=i-1;a[nz++]=-1;}
      idx[nz]=i;a[nz++]=4;
      if(i+1<m){idx[nz]=i+1;a[nz++]=-1.5;}
      if(m<n){idx[nz]=m+i%(n-m);a[nz++]=0.5;}
    } else {
      idx[nz]=i;a[nz++]=3;
      if(i+1<n){idx[nz]=i+1;a[nz++]=0.25;}
    }
  }
  ptr[n]=nz;return nz;
}
static kls_solver *make(int n,const int *ptr,const int *idx,const double *a,int app) {
  kls_solver *s=NULL;kls_options o;kls_default_options(&o);o.threads=1;
  assert(kls_create(&s)==KLS_OK);
  if(app)assert(kls_set_accuracy_policy(s,KLS_ACCURACY_APPLICATION_MANAGED)==KLS_OK);
  assert(kls_analyze_csr(s,KLS_INDEX_INT32,n,ptr,idx,0,&o)==KLS_OK);
  assert(kls_factor(s,a)==KLS_OK);
  return s;
}
static double residual(int n,const int *ptr,const int *idx,const double *a,
    const double *b,const double *x) {
  double worst=0,scale=0;
  for(int i=0;i<n;++i) {
    double r=b[i],row=fabs(b[i]);
    for(int p=ptr[i];p<ptr[i+1];++p){r-=a[p]*x[idx[p]];row+=fabs(a[p]*x[idx[p]]);}
    if(fabs(r)>worst)worst=fabs(r);
    if(row>scale)scale=row;
  }
  return worst/scale;
}
int main(void) {
  static const int sizes[]={14,37,60};
  for(size_t c=0;c<sizeof sizes/sizeof *sizes;++c) {
    const int n=sizes[c];
    int *ptr=malloc((size_t)(n+1)*sizeof *ptr),*idx=malloc((size_t)5*n*sizeof *idx);
    double *a=malloc((size_t)5*n*sizeof *a),*b=malloc((size_t)n*sizeof *b);
    double *x=malloc((size_t)n*sizeof *x),*y=malloc((size_t)n*sizeof *y);
    assert(ptr && idx && a && b && x && y);
    const int nz=build(n,ptr,idx,a);
    kls_solver *fused=make(n,ptr,idx,a,1),*split=make(n,ptr,idx,a,1);
    kls_solver *certified=make(n,ptr,idx,a,0);
    fused_events=lane_events=0;
    int replays_after_full=0;
    for(int k=0;k<400;++k) {
      for(int p=0;p<nz;++p)a[p]*=1+1e-3*((p+k)%5-2);
      for(int i=0;i<n;++i)b[i]=1+(i+k)%7;
      const int before=fused_events,lane_before=lane_events;
      assert(kls_refactor_solve(certified,a,1,b,0,y,0)==KLS_OK);
      assert(fused_events==before);             /* certified mode never fuses */
      assert(kls_refactor(split,a)==KLS_OK && kls_solve(split,1,b,0,y,0)==KLS_OK);
      assert(fused_events==before);             /* nor do separate calls */
      assert(lane_events==lane_before);         /* and neither replays */
      if(k%97==96) {
        /* Separate calls on the replaying solver switch replay off; the
           following combined call must be a full one. */
        if(k&1)assert(kls_solve(fused,1,b,0,x,0)==KLS_OK);
        else assert(kls_refactor(fused,a)==KLS_OK);
        assert(!fused->fast_lane_armed);
        const int events=lane_events;
        assert(kls_refactor_solve(fused,a,1,b,0,x,0)==KLS_OK);
        assert(lane_events==events);
      }
      if(fused->fast_lane_armed &&
         fused->fast_lane_budget==KLS_FAST_LANE_BUDGET)++replays_after_full;
      if(k%3==2) {                              /* aliased right-hand side */
        memcpy(x,b,(size_t)n*sizeof *x);
        assert(kls_refactor_solve(fused,a,1,x,0,x,0)==KLS_OK);
      } else assert(kls_refactor_solve(fused,a,1,b,0,x,0)==KLS_OK);
      double scale=0;
      for(int i=0;i<n;++i)if(fabs(y[i])>scale)scale=fabs(y[i]);
      for(int i=0;i<n;++i)assert(fabs(x[i]-y[i])<=1e-12*scale);
      assert(residual(n,ptr,idx,a,b,x)<=1e-12);
    }
    assert(fused_events>0 && lane_events>0);
    /* Replay is interrupted by a full call more than once. */
    assert(replays_after_full>2);
    /* While replay is armed, every irregular call fails exactly as a full
       call does: nothing is published and the next good update is solved
       from its own values. */
    for(int kind=0;kind<3;++kind) {
      for(int k=0;k<2*KLS_FAST_LANE_STREAK && !fused->fast_lane_armed;++k)
        assert(kls_refactor_solve(fused,a,1,b,0,x,0)==KLS_OK);
      assert(fused->fast_lane_armed && fused->fast_lane_budget>0);
      const double kept=a[ptr[n-1]],rhs=b[0];
      memcpy(y,x,(size_t)n*sizeof *y);
      if(kind==0)a[ptr[n-1]]=0;                 /* last singleton's pivot */
      else if(kind==1)a[ptr[n-1]]=NAN;
      else b[0]=INFINITY;
      assert(kls_refactor_solve(fused,a,1,b,0,x,0)!=KLS_OK);
      assert(!memcmp(x,y,(size_t)n*sizeof *x));
      assert(!fused->fast_lane_armed);
      a[ptr[n-1]]=kept;b[0]=rhs;
      assert(kls_refactor_solve(fused,a,1,b,0,x,0)==KLS_OK);
      assert(residual(n,ptr,idx,a,b,x)<=1e-12);
      assert(kls_refactor(split,a)==KLS_OK && kls_solve(split,1,b,0,y,0)==KLS_OK);
      for(int i=0;i<n;++i)assert(fabs(x[i]-y[i])<=1e-12*(1+fabs(y[i])));
    }
    kls_destroy(fused);kls_destroy(split);kls_destroy(certified);
    free(ptr);free(idx);free(a);free(b);free(x);free(y);
  }
  return 0;
}
