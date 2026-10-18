#include <subtree_lu.h>
#include <cfenv>
#include <cmath>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <xmmintrin.h>

// Compare every public status, input/output byte and immediate FP/errno state
// between unhooked original-library calls and the diagnostic recorder.
static void bytes(const char *name,const void *p,size_t n) {
 std::printf(" %s=",name);
 for(size_t i=0;i<n;++i)std::printf("%02x",static_cast<const unsigned char *>(p)[i]);
}
static void state(const char *edge,int stage,int status=0) {
 const int error=errno,exceptions=std::fetestexcept(FE_ALL_EXCEPT),round=std::fegetround();
 const unsigned mxcsr=_mm_getcsr();unsigned short cw;asm volatile("fnstcw %0":"=m"(cw));
 std::printf("STATE %s stage=%d status=%d errno=%d flags=%d round=%d mxcsr=%u cw=%u\n",edge,stage,status,error,exceptions,round,mxcsr,unsigned(cw));
 errno=error;
}
int main(int argc,char **argv) {
 if(argc!=3)return 90;
 const int threads=std::atoi(argv[1]),permutations=std::atoi(argv[2]);
 if((threads!=1 && threads!=8)||(permutations!=0 && permutations!=1))return 91;
 int ap[]={0,3,6,9,12},ai[]={0,1,3,0,1,2,1,2,3,0,2,3};
 double ax[]={8,1,-.25,-1,7,.5,1,6,-.5,.25,1,9};
 int rp[]={0,1,2,3},cp[]={0,1,2,3};
 subtree_lu::SubtreeLU<int,double> solver;
 solver.parm[subtree_lu::I_PIVTOL]=1000;solver.parm[subtree_lu::I_MEM_FACTOR]=150;
 errno=EDOM;std::feraiseexcept(FE_INEXACT);
 std::printf("INPUT analyze");bytes("ap",ap,sizeof(ap));bytes("ai",ai,sizeof(ai));bytes("ax",ax,sizeof(ax));if(permutations){bytes("rp",rp,sizeof(rp));bytes("cp",cp,sizeof(cp));}std::puts("");
 state("enter",1);int status=solver.analyze(4,ap,ai,ax,threads,permutations?rp:nullptr,permutations?cp:nullptr);state("leave",1,status);
 if(status){return threads==8 && permutations && status==-1?0:1;}
 if(permutations){std::printf("OUTPUT analyze");bytes("rp",rp,sizeof(rp));bytes("cp",cp,sizeof(cp));std::puts("");}
 const double wanted[]={1,-2,3,-4};
 for(int generation=0;generation<4;++generation) {
  if(generation==1){ax[0]+=2;ax[4]+=2;ax[7]+=2;ax[11]+=2;}
  double good[12];std::memcpy(good,ax,sizeof(ax));if(generation==2)std::memset(ax,0,sizeof(ax));
  const int stage=generation==1 || generation==2?3:2;
  std::printf("INPUT numeric stage=%d",stage);bytes("ax",ax,sizeof(ax));std::puts("");
  state("enter",stage);status=stage==3?solver.refactorize(ax):solver.factorize(ax);state("leave",stage,status);
  if(generation==2){if(status!=-3)return 2;std::memcpy(ax,good,sizeof(ax));continue;}
  if(status)return 3;
  for(int alias=0;alias<2;++alias) {
   double rhs[4]={},out[4];
   for(int i=0;i<4;++i)for(int p=ap[i];p<ap[i+1];++p)rhs[i]+=ax[p]*wanted[ai[p]];
   std::printf("INPUT solve alias=%d",alias);bytes("b",rhs,sizeof(rhs));std::puts("");
   state("enter",4);status=solver.solve(rhs,alias?rhs:out);state("leave",4,status);
   if(status)return 4;
   std::printf("OUTPUT solve");bytes("x",alias?rhs:out,sizeof(out));std::puts("");
   for(int i=0;i<4;++i)if(!std::isfinite((alias?rhs:out)[i]) || std::fabs((alias?rhs:out)[i]-wanted[i])>1e-12)return 5;
  }
 }
 return 0;
}
