#define _GNU_SOURCE
#define main existing_product_main
#include "kls_snb_product.c"
#undef main
#include <fenv.h>
#include <sys/mman.h>
#include <unistd.h>

#if defined(__AVX512F__)
/* Independent scalar oracle for the explicitly fused contiguous kernel.
   Portable builds retain the existing scalar expression and reference. */
static void reference_trsm_fused(double *B, const double *L, int sc, int w) {
  for(int i=0;i<sc;++i)
    for(int c=0;c<w;++c) {
      const double u=B[c*sc+i];
      if(u==0.0)continue;
      for(int j=i+1;j<sc;++j)
        B[c*sc+j]=fma(-L[i*sc+j],u,B[c*sc+j]);
    }
}
#endif

/* Compare against a scalar reference, ending both input panels at guard pages.
   Force the dense path for every triangular width, including all tail sizes.
   Compare output bits and exception flags for finite and nonfinite inputs. */
int main(void) {
  const long page_query = sysconf(_SC_PAGESIZE);
  if (page_query <= 0) return 1;
  const size_t page = (size_t)page_query;
  const size_t span = ((32u*32u*sizeof(double)+page-1)/page)*page;
  char *map = mmap(NULL, 2*span+2*page, PROT_READ|PROT_WRITE,
                   MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
  if (map == MAP_FAILED) return 1;
  if (mprotect(map+span,page,PROT_NONE) ||
      mprotect(map+2*span+page,page,PROT_NONE)) return 2;
  const int widths[] = {1,4,7,32};
  int32_t sub[32];
  double expected[2048], scratch[1];
#if !defined(__AVX512F__)
  double ub[1024];
#endif
  for(int sc=1;sc<=32;++sc)for(unsigned wi=0;wi<4;++wi)
    for(int pattern=0;pattern<5;++pattern) {
      const int w=widths[wi];
      double *L=(double *)(map+span)-sc*sc;
      double *B=(double *)(map+2*span+page)-sc*w;
      kls_snb_snode nodes[2]={{0}};
      nodes[0].width=sc; nodes[0].height=sc;
      nodes[1].width=w; nodes[1].height=sc;
      kls_snb_edge edge={0}; edge.sub_count=sc;
      kls_snb_block blk={0}; blk.snodes=nodes; blk.subcols=sub;
      int32_t dummy_dst=0; blk.edge_dsts=&dummy_dst;
      for(int i=0;i<sc;++i)sub[i]=i;
      for(int i=0;i<sc*sc;++i)L[i]=expected[i]=(i%7+1)*0.125;
      for(int c=0;c<w;++c)for(int i=0;i<sc;++i) {
        double value=(c+i+1)*0.25;
        if(pattern==1)value=i%2?0.0:-0.0;
        if(pattern==2&&i==0)value=INFINITY;
        if(pattern==3&&i==0)value=NAN;
        if(pattern==4&&i==0)value=-INFINITY;
        B[c*sc+i]=expected[1024+c*sc+i]=value;
      }
      blk.panels=expected; nodes[1].panel=1024;
      feclearexcept(FE_ALL_EXCEPT);
#if defined(__AVX512F__)
      reference_trsm_fused(expected+1024,expected,sc,w);
#else
      reference_edge(&blk,&nodes[1],&edge,1,ub,scratch);
#endif
      const int flags=fetestexcept(FE_ALL_EXCEPT);
      blk.panels=L; nodes[1].panel=(B-L);
      feclearexcept(FE_ALL_EXCEPT);
      kls_snb_update_edge(&blk,&nodes[1],&edge,1,scratch);
      const int actual_flags=fetestexcept(FE_ALL_EXCEPT);
      if(flags!=actual_flags ||
         memcmp(B,expected+1024,(size_t)sc*w*sizeof(double))) {
        fprintf(stderr,"TRSM contract sc=%d w=%d pattern=%d\n",sc,w,pattern);
        fprintf(stderr,"flags expected=%x actual=%x\n",flags,actual_flags);
        for(int j=0;j<sc*w;++j)
          if(memcmp(B+j,expected+1024+j,sizeof(double))) {
            fprintf(stderr,"first mismatch index=%d expected=%a actual=%a\n",
                    j,expected[1024+j],B[j]);
            break;
          }
        return 3;
      }
    }
  return munmap(map,2*span+2*page)!=0;
}
