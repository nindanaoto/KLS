#include <kls/kls.h>
#include "../src/kls_componentwise.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do { if(!(c)) {fprintf(stderr,"line %d: %s\n",__LINE__,#c);return 1;} } while(0)
int main(void) {
  int64_t p[]={0,1,3},r[]={0,0,1};
  double a[]={1,0x1p100,1},b[]={1,1},x[]={-0x1p100,1};
  long double bound;
  CHECK(kls_componentwise_certify(2,p,r,a,b,x,0,NULL,&bound)==1);
  CHECK(bound<1e-8L);
  double bad[]={0,1};
  CHECK(kls_componentwise_certify(2,p,r,a,b,bad,0,NULL,&bound)==0);
  int64_t dp[]={0,3},dr[]={0,0,0};
  double da[]={0x1p54,1,-0x1p54},db[]={1},dx[]={0};
  CHECK(kls_componentwise_certify(1,dp,dr,da,db,dx,0,NULL,&bound)==0);
  int64_t zp[]={0,1},zr[]={0};double za[]={1},zb[]={0},zx[]={0};
  CHECK(kls_componentwise_certify(1,zp,zr,za,zb,zx,0,NULL,&bound)==1);
  zx[0]=NAN;CHECK(kls_componentwise_certify(1,zp,zr,za,zb,zx,0,NULL,&bound)==0);
  for(int trans=0;trans<2;++trans) for(int btf=0;btf<2;++btf)
  for(int scale=0;scale<3;scale+=2) for(int csr=0;csr<2;++csr) {
    kls_solver *s=NULL;kls_accuracy_policy policy;
    CHECK(kls_create(&s)==KLS_OK);
    CHECK(kls_get_accuracy_policy(s,&policy)==KLS_OK && policy==KLS_ACCURACY_COMPONENTWISE_BACKWARD_ERROR);
    CHECK(kls_set_accuracy_policy(s,KLS_ACCURACY_COMPONENTWISE_BACKWARD_ERROR)==KLS_OK);
    kls_options o;kls_default_options(&o);o.use_btf=btf;o.scale=scale;
    o.ordering=KLS_ORDERING_NATURAL;o.orientation=KLS_ORIENTATION_NORMAL;
    int64_t cp[]={0,2,3},cr[]={0,1,1};
    CHECK((csr?kls_analyze_csr(s,KLS_INDEX_INT64,2,cp,cr,0,&o):
               kls_analyze_csc(s,KLS_INDEX_INT64,2,p,r,0,&o))==KLS_OK);
    CHECK(kls_set_accuracy_policy(s,KLS_ACCURACY_STRICT_RHS_L2)==KLS_ERR_INVALID_ARGUMENT);
    CHECK(kls_factor(s,a)==KLS_OK);
    double answer[]={1,1,99,1,1,99};
    CHECK((trans?kls_solve_transpose(s,2,answer,3,answer,3):
                 kls_solve(s,2,answer,3,answer,3))==KLS_OK);
    CHECK(answer[2]==99 && answer[5]==99);
    CHECK(kls_componentwise_certify(2,p,r,a,b,answer,trans,NULL,&bound)==1);
    CHECK(kls_refactor_solve(s,a,1,b,0,answer,0)==KLS_OK);
    CHECK(kls_componentwise_certify(2,p,r,a,b,answer,0,NULL,&bound)==1);
    double fail[]={NAN,1},sentinel[]={41,42};
    CHECK(kls_solve(s,1,fail,0,sentinel,0)!=KLS_OK);
    CHECK(sentinel[0]==41 && sentinel[1]==42);
    kls_destroy(s);
  }
  return 0;
}
