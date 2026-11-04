#if defined(__AVX512F__) && defined(__AVX512DQ__) && FLT_EVAL_METHOD==0
static unsigned row8_successes;
static void row8_compare(int rows,int count,int holes,int gapped,int pattern) {
  enum {PH=160,H=288,PANEL=5200,SIZE=10000};
  double original[SIZE],expected[SIZE],actual[SIZE],R[16*H],T[16*97];int32_t sub[32],dest[97];
  for(int k=0;k<SIZE;++k)original[k]=.015625*((k*7)%17-8);
  for(int i=0;i<count;++i)sub[i]=gapped?2*i:i;
  for(int r=0;r<rows;++r)dest[r]=holes==2 || (holes==1 && r%3==0)?-1:64+(holes==3?2*r:r);
  if(pattern==1)for(int c=0;c<32;++c)for(int r=0;r<rows;++r)original[c*PH+32+r]=DBL_TRUE_MIN;
  if(pattern==2)for(int c=0;c<16;++c)for(int i=0;i<count;++i)original[PANEL+c*H+4+i]=(c%2)?NAN:INFINITY;
  if(pattern==3)for(int c=0;c<32;++c)for(int r=0;r<32;++r)original[c*PH+r]=NAN;
  if(pattern==4)for(int c=0;c<32;++c)for(int r=0;r<rows;++r)original[c*PH+32+r]=(r%2)?INFINITY:NAN;
  if(pattern==5 || pattern==6)for(int c=0;c<16;++c)for(int i=0;i<count;++i)
    if(pattern==5 || c%2==0)original[PANEL+c*H+4+i]=-0.;
  memcpy(expected,original,sizeof original);memcpy(actual,original,sizeof original);
  kls_snb_snode nodes[2]={{0}};nodes[0].width=32;nodes[0].height=PH;nodes[0].below=rows;
  nodes[1].width=16;nodes[1].height=H;nodes[1].panel=PANEL;
  kls_snb_edge edge={0};edge.sub_count=count;edge.u_start=4;
  kls_snb_block block={0};block.snodes=nodes;block.subcols=sub;block.edge_dsts=dest;
  unsigned control=_mm_getcsr()&~63u;_mm_setcsr(control);
  block.panels=expected;reference_narrow_edge(&block,nodes+1,&edge,8,T);
  const unsigned expected_flags=_mm_getcsr();_mm_setcsr(control);
  for(int c=0;c<16;++c)for(int r=0;r<H;++r)R[(c/8)*8*H+r*8+c%8]=original[PANEL+c*H+r];
  block.panels=actual;int success=0;
  if((control&0x7f80u)==0x1f80u) {
    success=kls_snb_row8_edge(&block,nodes+1,&edge,8,R,T);
    if(_mm_getcsr()&31u)success=0;
    for(int r=0;r<16*H && success;++r)if(!kls_snb_row8_scalar(R[r]))success=0;
  }
  if(success) {
    for(int c=0;c<16;++c)for(int r=0;r<H;++r)actual[PANEL+c*H+r]=R[(c/8)*8*H+r*8+c%8];
    ++row8_successes;
  } else {
    _mm_setcsr(control);reference_narrow_edge(&block,nodes+1,&edge,8,T);
  }
  const unsigned observed=_mm_getcsr();
  if(expected_flags!=observed)fprintf(stderr,"ROW8 flags rows=%d count=%d holes=%d gap=%d pattern=%d success=%d expected=%x actual=%x\n",rows,count,holes,gapped,pattern,success,expected_flags,observed);
  assert(expected_flags==observed);
  if(memcmp(actual,expected,sizeof actual)) {
    for(int k=0;k<SIZE;++k)if(memcmp(actual+k,expected+k,sizeof(double))) {
      fprintf(stderr,"ROW8 value rows=%d count=%d holes=%d gap=%d pattern=%d success=%d index=%d expected=%a actual=%a\n",rows,count,holes,gapped,pattern,success,k,expected[k],actual[k]);break;
    }
    assert(!memcmp(actual,expected,sizeof actual));
  }
  if(!pattern && (control&0x7f80u)==0x1f80u)assert(success);
  _mm_setcsr(control);
}
static void row8_checks(void) {
  fenv_t saved;assert(!feholdexcept(&saved));assert(!fesetround(FE_TONEAREST));
  double before[64],after[64];__m512d tile[8];
  for(int i=0;i<64;++i){uint64_t bits=UINT64_C(0x7ff0000000000000)^(UINT64_C(0x123456789abcdef)*i);memcpy(before+i,&bits,8);}
  unsigned flags=_mm_getcsr();
  for(int i=0;i<8;++i)tile[i]=_mm512_loadu_pd(before+8*i);
  kls_snb_transpose8(tile);for(int i=0;i<8;++i)_mm512_storeu_pd(after+8*i,tile[i]);
  for(int i=0;i<8;++i)for(int j=0;j<8;++j)assert(!memcmp(after+8*i+j,before+8*j+i,8));
  kls_snb_transpose8(tile);for(int i=0;i<8;++i)_mm512_storeu_pd(after+8*i,tile[i]);
  assert(!memcmp(before,after,sizeof before) && flags==_mm_getcsr());
  const int sizes[]={0,1,2,7,8,9,31,32,33,65,97},counts[]={1,2,7,8,16,17,32};
  for(unsigned r=0;r<sizeof sizes/sizeof *sizes;++r)for(unsigned c=0;c<sizeof counts/sizeof *counts;++c)
    for(int holes=0;holes<4;++holes)for(int gap=0;gap<2;++gap)
      if(!gap || counts[c]<=16)row8_compare(sizes[r],counts[c],holes,gap,0);
  const int modes[]={FE_TONEAREST,FE_UPWARD,FE_DOWNWARD,FE_TOWARDZERO};
  for(int mode=0;mode<4;++mode)for(int flush=0;flush<2;++flush) {
    assert(!fesetround(modes[mode]));unsigned control=_mm_getcsr();
    _mm_setcsr((control&~((1u<<15)|(1u<<6)))|(flush?((1u<<15)|(1u<<6)):0));
    for(int pattern=1;pattern<=6;++pattern)for(int count=1;count<=16;count*=4)
      for(int holes=0;holes<4;++holes)row8_compare(33,count,holes,0,pattern);
  }
  assert(row8_successes>100);assert(!fesetenv(&saved));
}
#endif
