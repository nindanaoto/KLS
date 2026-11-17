/* Settled-solve replay: separate solves on an unchanged small system must
   return what the full solve returns, replay only for an application-managed
   caller on the compact walk, publish nothing on an irregular right-hand
   side, and keep taking a full solve at the budget. */
static int lane_events,combined_events;
#define KLS_SOLVE_LANE_TEST_EVENT() (++lane_events)
#define KLS_FAST_LANE_TEST_EVENT() (++combined_events)
#include "../src/kls.c"
#undef NDEBUG
#include <assert.h>

/* `blocks` tridiagonal blocks of `m` rows chained by one upper entry each;
   one block is a single strongly connected system. */
static int build(int blocks,int m,int *ptr,int *idx,double *a) {
  const int n=blocks*m;int nz=0;
  for(int i=0;i<n;++i) {
    const int local=i%m;
    ptr[i]=nz;
    if(local>0){idx[nz]=i-1;a[nz++]=-1;}
    idx[nz]=i;a[nz++]=4+0.125*(i%3);
    if(local+1<m){idx[nz]=i+1;a[nz++]=-1.5;}
    else if(i+1<n){idx[nz]=i+1;a[nz++]=0.5;}
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
  static const int shapes[][2]={{1,2},{1,5},{1,12},{3,4},{7,5}};
  for(size_t c=0;c<sizeof shapes/sizeof *shapes;++c) {
    const int n=shapes[c][0]*shapes[c][1];
    int *ptr=malloc((size_t)(n+1)*sizeof *ptr),*idx=malloc((size_t)3*n*sizeof *idx);
    double *a=malloc((size_t)3*n*sizeof *a),*b=malloc((size_t)n*sizeof *b);
    double *x=malloc((size_t)n*sizeof *x),*y=malloc((size_t)n*sizeof *y);
    assert(ptr && idx && a && b && x && y);
    const int nz=build(shapes[c][0],shapes[c][1],ptr,idx,a);
    kls_solver *lane=make(n,ptr,idx,a,1),*full=make(n,ptr,idx,a,1);
    kls_solver *certified=make(n,ptr,idx,a,0);
    lane_events=0;
    int replays=0,full_after_budget=0;
    for(int k=0;k<600;++k) {
      for(int i=0;i<n;++i)b[i]=1+(i+k)%7;
      if(k%40==0) {
        /* A new matrix arrives through the combined call on every solver. */
        for(int p=0;p<nz;++p)a[p]*=1+1e-3*((p+k)%5-2);
        assert(kls_refactor_solve(lane,a,1,b,0,x,0)==KLS_OK);
        assert(kls_refactor_solve(full,a,1,b,0,y,0)==KLS_OK);
        assert(kls_refactor_solve(certified,a,1,b,0,y,0)==KLS_OK);
      }
      int events=lane_events;
      assert(kls_solve(certified,1,b,0,y,0)==KLS_OK);
      assert(lane_events==events && !certified->solve_lane_armed);
      full->solve_lane_armed=0;full->solve_lane_streak=0;  /* always a full solve */
      assert(kls_solve(full,1,b,0,y,0)==KLS_OK);
      assert(lane_events==events);
      const int armed=lane->solve_lane_armed;
      const unsigned budget=lane->solve_lane_budget;
      if(k%3==2) {                              /* aliased right-hand side */
        memcpy(x,b,(size_t)n*sizeof *x);
        assert(kls_solve(lane,1,x,0,x,0)==KLS_OK);
      } else assert(kls_solve(lane,1,b,0,x,0)==KLS_OK);
      if(lane_events!=events) {assert(armed && budget>0);++replays;}
      else if(armed) {assert(!budget);++full_after_budget;}
      for(int i=0;i<n;++i)assert(x[i]==y[i]);   /* the same walk, the same bits */
      assert(residual(n,ptr,idx,a,b,x)<=1e-12);
    }
    assert(replays>400 && full_after_budget>0);
    /* An irregular right-hand side is declined before anything is published
       and switches replay off. */
    while(!lane->solve_lane_armed)assert(kls_solve(lane,1,b,0,x,0)==KLS_OK);
    for(int i=0;i<n;++i)y[i]=91;
    b[n/2]=NAN;
    const int events=lane_events;
    assert(kls_solve(lane,1,b,0,y,0)!=KLS_OK && lane_events==events);
    assert(!lane->solve_lane_armed);
    for(int i=0;i<n;++i)assert(y[i]==91);
    /* A separate refactor switches it off as well. */
    b[n/2]=1;
    while(!lane->solve_lane_armed)assert(kls_solve(lane,1,b,0,x,0)==KLS_OK);
    assert(kls_refactor(lane,a)==KLS_OK && !lane->solve_lane_armed);
    /* One new matrix followed by a few solves, repeated: both replays stay
       in use, and every answer matches separate full calls. */
    lane_events=combined_events=0;
    for(int k=0;k<400;++k) {
      for(int p=0;p<nz;++p)a[p]*=1+1e-3*((p+k)%5-2);
      for(int i=0;i<n;++i)b[i]=1+(i+k)%7;
      assert(kls_refactor_solve(lane,a,1,b,0,x,0)==KLS_OK);
      assert(kls_refactor(full,a)==KLS_OK && kls_solve(full,1,b,0,y,0)==KLS_OK);
      double scale=0;
      for(int i=0;i<n;++i)if(fabs(y[i])>scale)scale=fabs(y[i]);
      for(int i=0;i<n;++i)assert(fabs(x[i]-y[i])<=1e-12*scale);
      for(int j=0;j<4;++j) {
        for(int i=0;i<n;++i)b[i]=2+(i+k+j)%5;
        assert(kls_solve(lane,1,b,0,x,0)==KLS_OK);
        assert(kls_solve(full,1,b,0,y,0)==KLS_OK);
        for(int i=0;i<n;++i)assert(fabs(x[i]-y[i])<=1e-12*scale+1e-300);
        assert(residual(n,ptr,idx,a,b,x)<=1e-12);
      }
    }
    assert(lane_events>1200 && combined_events>300);
    kls_destroy(lane);kls_destroy(full);kls_destroy(certified);
    free(ptr);free(idx);free(a);free(b);free(x);free(y);
  }
  return 0;
}
