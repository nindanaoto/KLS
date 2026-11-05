/* Direct bitwise, scratch-bound and environment oracles for the repair move. */
#undef NDEBUG
#include <assert.h>
#include <fenv.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(__x86_64__) || defined(__i386__)
#include <xmmintrin.h>
#endif
#include "../src/kls_snb_inplace_panel.inc"
static void exercise(int height) {
  const size_t count=(size_t)height*16,bits=count/8;
  uint64_t *data=malloc((count+2)*sizeof(*data)),*reference=malloc(count*sizeof(*reference));
  unsigned char *scratch=malloc(bits+2);assert(data&&reference&&scratch);
  static const uint64_t patterns[]={0,UINT64_C(0x8000000000000000),1,
    UINT64_C(0x0010000000000000),UINT64_C(0x7ff0000000000000),
    UINT64_C(0xfff0000000000000),UINT64_C(0x7ff8000000000123),
    UINT64_C(0x7ff0000000000456),UINT64_C(0x3ff0000000000000)};
  data[0]=data[count+1]=UINT64_C(0x23456789abcdef01);
  for(size_t i=0;i<count;++i) {
    data[i+1]=patterns[i%(sizeof(patterns)/sizeof(patterns[0]))]^(i<<10);
    size_t column=i/(size_t)height,row=i-column*(size_t)height;
    reference[(column/8)*8*(size_t)height+row*8+column%8]=data[i+1];
  }
  memset(scratch,0xa5,bits+2);
  assert(!kls_snb_column_panel_to_rows_inplace((double *)(data+1),height,scratch+1,bits-1));
  assert(!kls_snb_column_panel_to_rows_inplace((double *)(data+1),height,data+1,bits));
  assert(!kls_snb_column_panel_to_rows_inplace((double *)(data+1),height,(unsigned char *)(data+1)+count*sizeof(*data)-1,bits));
  if(bits>sizeof(*data))assert(!kls_snb_column_panel_to_rows_inplace((double *)(data+1),height,data,bits));
  for(size_t i=0;i<count;++i)assert(data[i+1]==(patterns[i%(sizeof(patterns)/sizeof(patterns[0]))]^(i<<10)));
  for(size_t i=0;i<bits+2;++i)assert(scratch[i]==0xa5);
  const int flags=fetestexcept(FE_ALL_EXCEPT),round=fegetround();
#if defined(__x86_64__) || defined(__i386__)
  const unsigned csr=_mm_getcsr();
#endif
  assert(kls_snb_column_panel_to_rows_inplace((double *)(data+1),height,scratch+1,bits));
  assert(!memcmp(data+1,reference,count*sizeof(*reference)));
  assert(data[0]==UINT64_C(0x23456789abcdef01)&&data[count+1]==data[0]);
  assert(scratch[0]==0xa5&&scratch[bits+1]==0xa5);
  assert(fetestexcept(FE_ALL_EXCEPT)==flags&&fegetround()==round);
#if defined(__x86_64__) || defined(__i386__)
  assert(_mm_getcsr()==csr);
#endif
  free(data);free(reference);free(scratch);
}
int main(void) {
  fenv_t saved;assert(!fegetenv(&saved));
  double dummy[16]={0};unsigned char work[2]={0};
  assert(!kls_snb_column_panel_to_rows_inplace(NULL,1,work,sizeof(work)));
  assert(!kls_snb_column_panel_to_rows_inplace(dummy,1,NULL,sizeof(work)));
  assert(!kls_snb_column_panel_to_rows_inplace(dummy,0,work,sizeof(work)));
  assert(!kls_snb_column_panel_to_rows_inplace(dummy,-1,work,sizeof(work)));
  assert(!kls_snb_column_panel_to_rows_inplace(dummy,INT32_MAX,work,0));
  const int rounds[]={FE_TONEAREST,FE_DOWNWARD,FE_UPWARD,FE_TOWARDZERO};
  const int heights[]={1,2,3,7,8,15,16,17,31,32,124,1967,4096,4097};
  for(unsigned r=0;r<4;++r)for(unsigned mode=0;mode<4;++mode) {
    assert(!fesetround(rounds[r]));
#if defined(__x86_64__) || defined(__i386__)
    unsigned csr=_mm_getcsr();csr=(csr&~0x8040u)|((mode&1)?0x8000u:0)|((mode&2)?0x40u:0);_mm_setcsr(csr);
#endif
    assert(!feclearexcept(FE_ALL_EXCEPT));assert(!feraiseexcept(FE_INEXACT));
    for(unsigned h=0;h<sizeof(heights)/sizeof(heights[0]);++h)exercise(heights[h]);
  }
  assert(!fesetenv(&saved));
  puts("INPLACE_PANEL_ORACLES heights=14 environments=16 guarded=1 nan_bits=1 flags=1");
  return 0;
}
