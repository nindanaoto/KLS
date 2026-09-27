#include "../src/kls.c"
extern int kls_accuracy_certify(size_t,const int64_t*,const int64_t*,const double*,const double*,const double*,int,double*);

/* A=[1 2^100;0 1], b=[1,1]. A residual <=sqrt(2)*1e-8
 * forces x2 near1 and x1 near-2^100. Both terms in row1 are then
 * multiples of2^47, so their sum cannot approximate1 to that accuracy.
 * This is an intentional, mathematically infeasible binary64 contract,
 * not a recovery-performance benchmark. Rejection is the correct result. */
int main(void) {
  int failures=0;
  int64_t p[]={0,1,3},i[]={0,0,1};
  for(int frame=0;frame<3;++frame)
  for(int configuration=0;configuration<4;++configuration)
  for(int btf=0;btf<2;++btf)
  for(int hard=0;hard<2;++hard)
  for(int transpose=0;transpose<2;++transpose)
  for(int mode=0;mode<4;++mode) {
    int threads=configuration==1?8:1;
    double a[]={1,hard?0x1p100:2,1};
    double b[]={1,1,0,1,1,0},x[]={0,0,0,0,0,0};
    int nrhs=(mode&2)?2:1,inplace=mode&1;
    kls_options o;kls_default_options(&o);o.threads=threads;o.use_btf=btf;
    if(frame==1) o.scale=2;
    if(frame==2) o.orientation=KLS_ORIENTATION_TRANSPOSE;
    if(configuration==2) o.backend=KLS_BACKEND_SERIAL;
    if(configuration==3) o.backend=KLS_BACKEND_KLS;
    kls_solver *s=NULL;int status=kls_create(&s);
    if(status==KLS_OK) status=kls_set_accuracy_policy(s,KLS_ACCURACY_STRICT_RHS_L2);
    if(status==KLS_OK) status=kls_analyze_csc(s,KLS_INDEX_INT64,2,p,i,0,&o);
    for(int epoch=0;epoch<3;++epoch) {
    memset(x,0,sizeof x);if(inplace) memcpy(x,b,sizeof x);
    /* Changed values followed by an identical-value refactor exercise both
     * snapshot renewal and exact-repeat early exits. */
    if(epoch==1) a[1]*=2;
    if(epoch==0) { if(status==KLS_OK) status=kls_factor(s,a); }
    else status=kls_refactor(s,a);
    if(status==KLS_OK) status=transpose
      ? kls_solve_transpose(s,nrhs,inplace?x:b,3,x,3)
      : kls_solve(s,nrhs,inplace?x:b,3,x,3);
    int certified=1;
    for(int rhs=0;rhs<nrhs;++rhs) {
      double ratio=INFINITY;
      if(!kls_accuracy_certify(2,p,i,a,b+3*rhs,x+3*rhs,transpose,&ratio)) certified=0;
    }
    int valid=hard ? status==KLS_ERR_SOLVE_FAILED : status==KLS_OK && certified;
    if(!valid || x[2]!=0 || x[5]!=0) ++failures;
    printf("frame=%d configuration=%d threads=%d btf=%d hard=%d transpose=%d mode=%d epoch=%d status=%d certified=%d valid=%d\n",
           frame,configuration,threads,btf,hard,transpose,mode,epoch,status,certified,valid);
    fflush(stdout);
    }
    kls_destroy(s);
  }
  return failures!=0;
}
