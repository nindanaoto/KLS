/* Compare the larger register tile to the frozen incumbent kernel, including
 * status flags, storage guards, arbitrary producer order and all row tails. */
#include "../src/kls.c"
#include "../src/kls_snb_row8.inc"
#include <fenv.h>
#undef NDEBUG
#include <assert.h>
#include "kls_snb_product_tile_reference.inc"

#if defined(__x86_64__) && FLT_EVAL_METHOD==0
static void compare(int nr,int width,int count,int offset,int pattern) {
  double l[141*37],b[66*32],actual[129*32+2],expected[129*32+2];
  int32_t sub[32];
  for(int i=0;i<32;++i) sub[i]=(i*7)%37;
  for(int i=0;i<141*37;++i) {
    l[i]=((i*17)%29-14)/11.0;
    if(pattern==1) l[i]=i%3?0x1p-1022:DBL_TRUE_MIN;
    if(pattern==2) l[i]=i%3?0x1p512:-0x1p512;
    if(pattern==3) l[i]=i%3?INFINITY:0;
    if(pattern==4) l[i]=i%3?NAN:-0.0;
  }
  for(int i=0;i<66*32;++i) {
    b[i]=((i*13)%23-11)/7.0;
    if(pattern==1) b[i]=i%3?0.5:0x1p-512;
    if(pattern==2) b[i]=i%3?0x1p512:-0x1p512;
    if(pattern==3) b[i]=i%3?0:INFINITY;
  }
  for(int i=0;i<129*32+2;++i) actual[i]=expected[i]=12345;
  unsigned before,reference_flags,actual_flags;
  __asm__ volatile("stmxcsr %0":"=m"(before)::"memory");
  before&=~0x3fu;
  __asm__ volatile("ldmxcsr %0"::"m"(before):"memory");
  reference_snb_dense_product(l,141,sub,count,b,66,width,offset,offset+nr,expected+1);
  __asm__ volatile("stmxcsr %0":"=m"(reference_flags)::"memory");
  __asm__ volatile("ldmxcsr %0"::"m"(before):"memory");
  kls_snb_dense_product(l,141,sub,count,b,66,width,offset,offset+nr,actual+1);
  __asm__ volatile("stmxcsr %0":"=m"(actual_flags)::"memory");
  assert(actual_flags==reference_flags);
  for(int i=0;i<129*32+2;++i) {
    if(isnan(expected[i])) assert(isnan(actual[i]));
    else assert(!memcmp(actual+i,expected+i,sizeof(double)));
  }
}
#endif
#include "kls_snb_narrow_checks.inc"
#include "kls_snb_row8_checks.inc"

int main(void) {
#if defined(__x86_64__) && FLT_EVAL_METHOD==0
  fenv_t saved;assert(!fegetenv(&saved));
  unsigned original;__asm__ volatile("stmxcsr %0":"=m"(original)::"memory");
  const int rows[]={0,1,7,8,9,15,16,17,23,24,25,31,32,33,63,64,65,127,128,129};
  const int widths[]={1,3,4,5,7,8,9,31,32},counts[]={1,2,7,8,17,32};
  assert(!fesetround(FE_TONEAREST));
#if defined(__AVX512F__) && defined(__AVX512DQ__)
  row8_checks();
#endif
  narrow_edge_checks();
  for(unsigned r=0;r<sizeof(rows)/sizeof(rows[0]);++r)
    for(unsigned w=0;w<sizeof(widths)/sizeof(widths[0]);++w)
      for(unsigned c=0;c<sizeof(counts)/sizeof(counts[0]);++c)
        for(int offset=0;offset<=3;offset+=3) compare(rows[r],widths[w],counts[c],offset,0);
  const int rounding[]={FE_TONEAREST,FE_UPWARD,FE_DOWNWARD,FE_TOWARDZERO};
  for(unsigned r=0;r<4;++r) {
    assert(!fesetround(rounding[r]));
    unsigned control;__asm__ volatile("stmxcsr %0":"=m"(control)::"memory");
    for(unsigned mode=0;mode<4;++mode) {
      unsigned selected=(control|0x1f80u)&~(0x3fu|0x8040u);
      selected|=(mode&1?1u<<15:0)|(mode&2?1u<<6:0);
      __asm__ volatile("ldmxcsr %0"::"m"(selected):"memory");
      for(int pattern=0;pattern<8;++pattern)for(int width=4;width<=9;++width)
        for(int count=1;count<=7;count+=2)for(int holes=0;holes<3;++holes)
          compare_narrow_edge(width,33,count,pattern,holes);
      for(int pattern=0;pattern<5;++pattern)
        for(int nr=16;nr<=33;++nr)
          for(int width=4;width<=9;++width) compare(nr,width,17,3,pattern);
    }
  }
  assert(!fesetenv(&saved));
  __asm__ volatile("ldmxcsr %0"::"m"(original):"memory");
#endif
  return 0;
}
