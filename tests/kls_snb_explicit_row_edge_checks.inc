/* Check row producers against the qualified column executor, with an
 * explicit external destination. Both the ordinary panel and producer must
 * remain unchanged. Finite random mantissas expose contraction differences. */
#ifdef KLS_SNB_HAVE_ROW_PANEL
static uint64_t explicit_edge_cases;
static double explicit_edge_value(unsigned index,int pattern) {
  uint64_t word=((uint64_t)index+17)*UINT64_C(0x9e3779b97f4a7c15);
  word^=word>>29;word*=UINT64_C(0xbf58476d1ce4e5b9);word^=word>>31;
  uint64_t bits=(word&UINT64_C(0x800fffffffffffff))|
    ((uint64_t)(1003+(word>>57)%41)<<52);
  if(pattern==1 && index%5==0)bits=word&UINT64_C(0x8000000000000000);
  if(pattern==2) {
    const uint64_t small[]={1,UINT64_C(0x8000000000000001),
      UINT64_C(0x0010000000000000),UINT64_C(0x8010000000000000)};
    bits=small[index%4];
  }
  double value;memcpy(&value,&bits,8);return value;
}
static void explicit_row_edge_case(int width,int count,int rows,int gap,
    int holes,int pattern,int math) {
  enum{PH=192,H=384,BASE=4096,SIZE=BASE+32*H};
  double *expected=malloc(SIZE*sizeof(double)),*actual=malloc(SIZE*sizeof(double));
  double *before=malloc(SIZE*sizeof(double)),*target=malloc((32*H+2)*sizeof(double));
  double scratch[128*32];int32_t sub[16],dest[95];
  assert(expected&&actual&&before&&target);
  kls_snb_snode nodes[2]={{0}};nodes[0].width=16;nodes[0].height=PH;
  nodes[0].uw=3;nodes[0].below=rows;nodes[1].width=width;
  nodes[1].height=H;nodes[1].panel=BASE;
  kls_snb_edge edge={0};edge.sub_count=count;edge.u_start=4;
  kls_snb_block block={0};block.snodes=nodes;block.subcols=sub;block.edge_dsts=dest;
  for(unsigned i=0;i<SIZE;++i)expected[i]=explicit_edge_value(i,pattern);
  memcpy(actual,expected,SIZE*sizeof(double));
  for(int c=0;c<16;++c)for(int r=0;r<PH;++r)
    memcpy(actual+kls_rpanel_at(PH,c,r),expected+c*PH+r,8);
  memcpy(before,actual,SIZE*sizeof(double));
  for(int i=0;i<count;++i)sub[i]=gap?2*i:i;
  for(int r=0;r<rows;++r)dest[r]=holes==2||(holes==1&&r%3==0)?-1:80+2*r;
  uint64_t guard=UINT64_C(0x123456789abcdef0);
  memcpy(target,&guard,8);memcpy(target+width*H+1,&guard,8);
  memcpy(target+1,expected+BASE,(size_t)width*H*sizeof(double));
  assert(!feclearexcept(FE_ALL_EXCEPT));unsigned control=_mm_getcsr()&~63u;
  _mm_setcsr(control);block.panels=expected;
  kls_snb_update_edge_to(&block,nodes+1,&edge,8,scratch,expected+BASE);
  const unsigned expected_csr=_mm_getcsr();const int expected_flags=fetestexcept(FE_ALL_EXCEPT);
  assert(!feclearexcept(FE_ALL_EXCEPT));_mm_setcsr(control);block.panels=actual;
  kls_rpanel_small_target_edge_to(&block,nodes+1,&edge,8,math,scratch,target+1);
  const unsigned actual_csr=_mm_getcsr();const int actual_flags=fetestexcept(FE_ALL_EXCEPT);
  if(expected_csr!=actual_csr||expected_flags!=actual_flags||
      memcmp(target+1,expected+BASE,(size_t)width*H*sizeof(double))) {
    fprintf(stderr,"EXPLICIT_ROW_EDGE width=%d count=%d rows=%d gap=%d holes=%d pattern=%d expected_flags=%x/%x actual_flags=%x/%x\n",
      width,count,rows,gap,holes,pattern,expected_csr,expected_flags,actual_csr,actual_flags);abort();
  }
  assert(!memcmp(before,actual,SIZE*sizeof(double)));
  assert(!memcmp(target,&guard,8)&&!memcmp(target+width*H+1,&guard,8));
  ++explicit_edge_cases;_mm_setcsr(control);
  free(expected);free(actual);free(before);free(target);
}
static void explicit_row_edge_checks(void) {
  fenv_t saved;assert(!feholdexcept(&saved));assert(!fesetround(FE_TONEAREST));
  const unsigned original=_mm_getcsr();explicit_edge_cases=0;
  const int saved_owned=kls_snb_owned_masked_environment;kls_snb_owned_masked_environment=1;
  const int widths[]={1,3,4,8,16,17,32},counts[]={1,2,3,7,8,15,16};
  const int rows[]={0,1,7,8,15,16,17,31,32,95};
  for(unsigned mode=0;mode<4;++mode) {
    _mm_setcsr((original&~(0x8040u|63u))|((mode&1)?0x8000u:0)|((mode&2)?0x40u:0));
    int math=kls_snb_row_product_math();assert(math>=0);
    for(unsigned w=0;w<7;++w)for(unsigned c=0;c<7;++c)for(unsigned r=0;r<10;++r)
      for(int gap=0;gap<2;++gap)if((counts[c]-1)*(gap?2:1)<16)
        for(int holes=0;holes<3;++holes)for(int pattern=0;pattern<3;++pattern)
          explicit_row_edge_case(widths[w],counts[c],rows[r],gap,holes,pattern,math);
  }
  assert(explicit_edge_cases==30240);
  fprintf(stderr,"EXPLICIT_ROW_EDGE_ORACLES cases=%llu ftz_daz=4 external_target=1 producer_immutable=1 flags=1\n",
    (unsigned long long)explicit_edge_cases);
  kls_snb_owned_masked_environment=saved_owned;assert(!fesetenv(&saved));
}
#else
static void explicit_row_edge_checks(void) {}
#endif
