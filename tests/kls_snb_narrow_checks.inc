/* Complete-edge oracle: original producer and target-column order. */
#include "kls_snb_narrow_reference.inc"
#if defined(__x86_64__) && FLT_EVAL_METHOD==0
static unsigned narrow_widths,narrow_rows,narrow_zero_masks;
static void compare_narrow_edge(int width,int rows,int count,int pattern,int holes) {
  double actual[4098],expected[4098],scratch[65*32];int32_t sub[7],dest[65];
  kls_snb_snode nodes[2]={{0}};
  nodes[0].width=7;nodes[0].height=73;nodes[0].below=rows;
  nodes[1].width=width;nodes[1].height=96;nodes[1].panel=1025;
  kls_snb_edge edge={0};edge.sub_count=count;edge.u_start=1;
  kls_snb_block block={0};block.snodes=nodes;block.subcols=sub;block.edge_dsts=dest;
  const uint64_t exceptional[]={1,UINT64_C(0x0010000000000000),
    UINT64_C(0x7fefffffffffffff),UINT64_C(0x7ff0000000000000),
    UINT64_C(0x7ff8000000000042),UINT64_C(0x7ff0000000000042),
    UINT64_C(0x8000000000000000)};
  for(int k=0;k<4098;++k) {
    actual[k]=.00390625*(1+k%13);
    if(pattern>0&&pattern<8)memcpy(actual+k,exceptional+(k+pattern-1)%7,sizeof(double));
  }
  for(int i=0;i<count;++i)sub[i]=pattern==25?2*i+1:i;
  for(int r=0;r<rows;++r)dest[r]=holes==2||(holes==1&&r%3==0)?-1:31+r;
  if(pattern==26)for(int r=0;r<rows;++r)if(dest[r]>=0)dest[r]=16+2*r;
  if(pattern>=8 && pattern<24) {
    const unsigned mask=(unsigned)(pattern-8);narrow_zero_masks|=1u<<mask;
    for(int c=0;c<width;++c)for(int i=0;i<count;++i)
      actual[1025+c*96+1+i]=(mask&(1u<<(c%4)))?.125:-0.0;
  }
  if(pattern==24)for(int i=0;i<count;++i)for(int r=0;r<rows;++r)
    memcpy(actual+(size_t)i*73+7+r,exceptional+(i+r)%7,sizeof(double));
  memcpy(expected,actual,sizeof(actual));
  unsigned control,flags,observed;__asm__ volatile("stmxcsr %0":"=m"(control)::"memory");
  control&=~0x3fu;__asm__ volatile("ldmxcsr %0"::"m"(control):"memory");
  block.panels=expected;reference_narrow_edge(&block,nodes+1,&edge,8,scratch);
  __asm__ volatile("stmxcsr %0":"=m"(flags)::"memory");
  __asm__ volatile("ldmxcsr %0"::"m"(control):"memory");
  const int saved_owned=kls_snb_owned_masked_environment;
  kls_snb_owned_masked_environment=fegetround()==FE_TONEAREST;
  block.panels=actual;kls_snb_update_edge(&block,nodes+1,&edge,8,scratch);
  kls_snb_owned_masked_environment=saved_owned;
  __asm__ volatile("stmxcsr %0":"=m"(observed)::"memory");
  if(flags!=observed)fprintf(stderr,"NARROW_FLAGS width=%d rows=%d count=%d pattern=%d holes=%d round=%d expected=%x actual=%x\n",width,rows,count,pattern,holes,fegetround(),flags,observed);
  assert(flags==observed);
  for(int k=0;k<4098;++k)assert((isnan(expected[k])&&isnan(actual[k]))||!memcmp(expected+k,actual+k,sizeof(double)));
  narrow_widths|=1u<<(width%4);narrow_rows|=1u<<(rows%8);
  __asm__ volatile("ldmxcsr %0"::"m"(control):"memory");
}
static void narrow_edge_checks(void) {
  const int rows[]={0,1,2,3,4,5,6,7,8,9,15,16,17,31,32,33,63,64,65};
  for(int width=1;width<=32;++width)for(int count=1;count<=7;++count)
    for(unsigned r=0;r<sizeof(rows)/sizeof(rows[0]);++r)
      for(int holes=0;holes<3;++holes)compare_narrow_edge(width,rows[r],count,0,holes);
  assert(narrow_widths==15&&narrow_rows==255);
  for(int holes=0;holes<3;++holes)for(int count=1;count<=3;++count) {
    compare_narrow_edge(16,33,count,25,holes); /* Gapped producer subset. */
    compare_narrow_edge(16,33,count,26,holes); /* Gapped target destinations. */
  }
  /* Every possible zero/nonzero combination across the four target columns. */
  for(int mask=0;mask<16;++mask)for(int count=1;count<=7;++count)
    compare_narrow_edge(9,33,count,8+mask,1);
  assert(narrow_zero_masks==65535);
  /* Exceptional values, all SSE rounding directions, gradual underflow
     and FTZ/DAZ must match the original scalar edge and its full flags. */
  fenv_t saved;assert(!feholdexcept(&saved));
  const int modes[]={FE_TONEAREST,FE_UPWARD,FE_DOWNWARD,FE_TOWARDZERO};
  for(unsigned mode=0;mode<4;++mode)for(unsigned flush=0;flush<2;++flush) {
    assert(!fesetround(modes[mode]));
    _mm_setcsr((_mm_getcsr()&~((1u<<15)|(1u<<6)))|(flush?((1u<<15)|(1u<<6)):0));
    for(int pattern=1;pattern<=7;++pattern)for(int count=1;count<=7;++count)
      for(int holes=0;holes<3;++holes) {
        compare_narrow_edge(9,33,count,pattern,holes);
        compare_narrow_edge(8,16,count,pattern,holes);
      }
    for(int width=1;width<=3;++width)for(int count=1;count<=7;++count)
      for(int pattern=1;pattern<=7;++pattern)for(int holes=0;holes<3;++holes) {
        compare_narrow_edge(width,33,count,pattern,holes);
        compare_narrow_edge(width,7,count,pattern,holes);
        compare_narrow_edge(width,2,count,24,holes);
        compare_narrow_edge(width,33,count,24,holes);
      }
    for(int mask=0;mask<16;++mask)for(int holes=0;holes<3;++holes)
      compare_narrow_edge(9,33,7,8+mask,holes);
  }
  assert(!fesetenv(&saved));
  unsigned control;__asm__ volatile("stmxcsr %0":"=m"(control)::"memory");
  double dummy[4*96]={0},lower[73]={0};int32_t subset[1]={0},dest[33];
  for(int r=0;r<33;++r)dest[r]=31+r;
  assert(kls_snb_narrow_columns(dummy,96,4,lower,lower,73,subset,1,dest,33,1)==4);
  const unsigned unmasked=control&~(1u<<7);
  __asm__ volatile("ldmxcsr %0"::"m"(unmasked):"memory");
  assert(kls_snb_narrow_columns(dummy,96,4,lower,lower,73,subset,1,dest,33,1)==0);
  __asm__ volatile("ldmxcsr %0"::"m"(control):"memory");
}
#endif
