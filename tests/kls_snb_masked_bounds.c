#define _GNU_SOURCE
#include "../src/kls.c"
#include <fenv.h>
#include <sys/mman.h>
#include <unistd.h>

/* Put the last legal element immediately before an inaccessible page.
   Every masked tail length must neither read nor write inactive lanes. */
int main(void) {
  const size_t page=(size_t)sysconf(_SC_PAGESIZE);
  const size_t span=2*page;
  void *maps[3];
  for(int i=0;i<3;++i) {
    maps[i]=mmap(NULL,span+page,PROT_READ|PROT_WRITE,
      MAP_PRIVATE|MAP_ANONYMOUS,-1,0);
    if(maps[i]==MAP_FAILED || mprotect((char *)maps[i]+span,page,PROT_NONE))return 1;
  }
  const int widths[]={1,3,4,5,7,8,31,32};
  const int counts[]={1,8,17,32};
  int32_t sub[32];for(int i=0;i<32;++i)sub[i]=i;
  for(int nr=1;nr<=7;++nr)
    for(unsigned wi=0;wi<sizeof(widths)/sizeof(*widths);++wi)
      for(unsigned si=0;si<sizeof(counts)/sizeof(*counts);++si)
        for(int infinite=0;infinite<2;++infinite) {
          const int w=widths[wi],sc=counts[si];
          if((size_t)sc*w*sizeof(double)>span)return 2;
          double *L=(double *)((char *)maps[0]+span)-sc*nr;
          double *B=(double *)((char *)maps[1]+span)-sc*w;
          double *T=(double *)((char *)maps[2]+span)-nr*w;
          for(int i=0;i<sc*nr;++i)L[i]=1.0;
          for(int i=0;i<sc*w;++i)B[i]=infinite?INFINITY:2.0;
          feclearexcept(FE_ALL_EXCEPT);
          kls_snb_dense_product(L,nr,sub,sc,B,sc,w,0,nr,T);
          if(fetestexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW))return 3;
          for(int i=0;i<nr*w;++i)
            if(infinite?T[i]!=INFINITY:T[i]!=2.0*sc)return 4;
        }
  for(int i=0;i<3;++i)if(munmap(maps[i],span+page))return 5;
  return 0;
}
