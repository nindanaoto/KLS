#include "kls_snb_narrow_reference.inc"

#if defined(KLS_SNB_HAVE_ROW_PANEL)
static unsigned row_normal_cases;
static int row_normal_math;
static char row_normal_context[192];
static double row_normal_value(int i,int pattern) {
  static const uint64_t bits[]={1,UINT64_C(0x0010000000000000),
    UINT64_C(0x7fefffffffffffff),UINT64_C(0x7ff0000000000000),
    UINT64_C(0x7ff8000000000042),UINT64_C(0x7ff0000000000042),
    UINT64_C(0x8000000000000000)};
  double x=(i%29-14)*0x1p-8;
  if(pattern>0&&pattern<=7)memcpy(&x,bits+(i+pattern-1)%7,8);
  return x;
}
static void row_normal_equal(const double *a,const double *b,int n,const char *kind) {
  for(int i=0;i<n;++i)if(!(isnan(a[i])&&isnan(b[i]))&&memcmp(a+i,b+i,8)) {
    fprintf(stderr,"%s %s value index=%d expected=%a actual=%a\n",row_normal_context,kind,i,a[i],b[i]);abort();
  }
}
static void row_normal_edge(int count,int rows,int pattern,int holes,int row_producer,int gap) {
  snprintf(row_normal_context,sizeof(row_normal_context),"edge count=%d rows=%d pattern=%d holes=%d layout=%d gap=%d",count,rows,pattern,holes,row_producer,gap);
  enum{H=288,PH=160,BASE=6200,SIZE=13000};
  double expected[SIZE],actual[SIZE],R[16*H],scratch[128*16],producer[32*PH];
  int32_t subset[32],dest[128];kls_snb_snode nodes[2]={{0}};
  nodes[0].width=row_producer?16:32;nodes[0].uw=3;
  nodes[0].height=PH;nodes[0].below=rows;
  nodes[1].width=16;nodes[1].height=H;nodes[1].panel=BASE;
  kls_snb_edge edge={0};edge.sub_count=count;edge.u_start=4;
  kls_snb_block block={0};block.snodes=nodes;block.subcols=subset;block.edge_dsts=dest;
  for(int i=0;i<SIZE;++i)expected[i]=row_normal_value(i,pattern);
  for(int i=0;i<count;++i)subset[i]=gap?2*i:i;
  for(int r=0;r<rows;++r)dest[r]=holes==2||(holes==1&&r%3==0)?-1:64+2*r;
  if(pattern>=8)for(int c=0;c<16;++c)for(int i=0;i<count;++i)
    expected[BASE+c*H+4+i]=((pattern-8)&(1u<<(c%4)))?.125:-0.;
  memcpy(actual,expected,sizeof(actual));
  for(int c=0;c<16;++c)for(int r=0;r<H;++r)R[kls_rpanel_at(H,c,r)]=expected[BASE+c*H+r];
  if(row_producer)for(int c=0;c<16;++c)for(int r=0;r<PH;++r)
    actual[kls_rpanel_at(PH,c,r)]=expected[c*PH+r];
  memcpy(producer,actual,sizeof(producer));
  unsigned control=_mm_getcsr()&~63u;_mm_setcsr(control);
  block.panels=expected;reference_narrow_edge(&block,nodes+1,&edge,8,scratch);
  unsigned expected_flags=_mm_getcsr();_mm_setcsr(control);
  block.panels=actual;assert(kls_rpanel_edge(&block,nodes+1,&edge,8,R,row_producer,row_normal_math));
  unsigned actual_flags=_mm_getcsr();
  if(expected_flags!=actual_flags) {
    fprintf(stderr,"edge flags count=%d rows=%d pattern=%d holes=%d layout=%d gap=%d expected=%x actual=%x\n",count,rows,pattern,holes,row_producer,gap,expected_flags,actual_flags);abort();
  }
  assert(!memcmp(producer,actual,sizeof(producer)));
  for(int c=0;c<16;++c)for(int r=0;r<H;++r)actual[BASE+c*H+r]=R[kls_rpanel_at(H,c,r)];
  row_normal_equal(expected+BASE,actual+BASE,16*H,"edge");++row_normal_cases;_mm_setcsr(control);
}
static __attribute__((noinline)) int row_normal_reference_getrf(double *W,int H,int uw,double *D) {
  int singular=-1;
  for(int i=0;i<16;++i) {
    double pivot=W[(size_t)i*H+uw+i];D[i]=pivot;
    if(pivot==0.&&singular<0)singular=i;
    for(int r=uw+i+1;r<H;++r)W[(size_t)i*H+r]/=pivot;
    for(int c=i+1;c<16;++c) {
      double u=W[(size_t)c*H+uw+i];if(u==0.)continue;
      for(int r=uw+i+1;r<H;++r)W[(size_t)c*H+r]-=W[(size_t)i*H+r]*u;
    }
  }
  return singular;
}
static void row_normal_getrf(int H,int pattern,int zero) {
  snprintf(row_normal_context,sizeof(row_normal_context),"getrf H=%d pattern=%d zero=%d",H,pattern,zero);
  double W[80*16],R[80*16],actual[80*16],D[16],RD[16];int uw=3;
  for(int c=0;c<16;++c)for(int r=0;r<H;++r) {
    double x=row_normal_value(c*H+r,pattern);
    if(!pattern&&r==uw+c)x=2.;
    if(zero>=0&&c==zero&&r==uw+c)x=0.;
    W[c*H+r]=R[kls_rpanel_at(H,c,r)]=x;
  }
  unsigned control=_mm_getcsr()&~63u;_mm_setcsr(control);
  int es=row_normal_reference_getrf(W,H,uw,D);unsigned ef=_mm_getcsr();_mm_setcsr(control);
  int as=kls_rpanel_getrf(R,H,uw,RD);unsigned af=_mm_getcsr();
  if(ef!=af||es!=as) {
    fprintf(stderr,"getrf H=%d pattern=%d zero=%d expected=%x/%d actual=%x/%d\n",H,pattern,zero,ef,es,af,as);abort();
  }
  for(int c=0;c<16;++c)for(int r=0;r<H;++r)actual[c*H+r]=R[kls_rpanel_at(H,c,r)];
  row_normal_equal(W,actual,16*H,"getrf");row_normal_equal(D,RD,16,"diagonal");++row_normal_cases;_mm_setcsr(control);
}
#endif
#if defined(KLS_SNB_HAVE_ROW_PANEL)
static void row_normal_small(int width,int rows,int count,int holes,int gap) {
  snprintf(row_normal_context,sizeof(row_normal_context),"small width=%d rows=%d count=%d holes=%d gap=%d",width,rows,count,holes,gap);
  enum{SIZE=18000,H=288,PH=160,BASE=6200};
  double expected[SIZE],actual[SIZE],producer[16*PH],scratch[128*32];
  int32_t sub[16],dst[95];kls_snb_snode nodes[2]={{0}};
  nodes[0].width=16;nodes[0].height=PH;nodes[0].uw=3;nodes[0].below=rows;
  nodes[1].width=width;nodes[1].height=H;nodes[1].panel=BASE;
  kls_snb_edge edge={0};edge.sub_count=count;edge.u_start=4;
  kls_snb_block block={0};block.snodes=nodes;block.subcols=sub;block.edge_dsts=dst;
  for(int k=0;k<SIZE;++k)expected[k]=row_normal_value(k,0);
  for(int i=0;i<count;++i)sub[i]=gap?2*i:i;
  for(int r=0;r<rows;++r)dst[r]=holes==2||(holes==1&&r%3==0)?-1:64+2*r;
  memcpy(actual,expected,sizeof(actual));
  for(int c=0;c<16;++c)for(int r=0;r<PH;++r)actual[kls_rpanel_at(PH,c,r)]=expected[c*PH+r];
  memcpy(producer,actual,sizeof producer);unsigned control=_mm_getcsr()&~63u;_mm_setcsr(control);
  block.panels=expected;reference_narrow_edge(&block,nodes+1,&edge,8,scratch);unsigned ef=_mm_getcsr();_mm_setcsr(control);
  block.panels=actual;kls_rpanel_small_target_edge(&block,nodes+1,&edge,8,row_normal_math,scratch);unsigned af=_mm_getcsr();
  if(ef!=af){fprintf(stderr,"small row-producer flags width=%d rows=%d count=%d expected=%x actual=%x\n",width,rows,count,ef,af);abort();}
  assert(!memcmp(producer,actual,sizeof producer));
  row_normal_equal(expected+BASE,actual+BASE,width*H,"small row producer");++row_normal_cases;_mm_setcsr(control);
}
#endif
static void kls_snb_row_normal_checks_mode(int ftz) {
#if defined(KLS_SNB_HAVE_ROW_PANEL)
  unsigned original=_mm_getcsr();fenv_t saved;assert(!feholdexcept(&saved));
  assert(!fesetround(FE_TONEAREST));_mm_setcsr(((_mm_getcsr()|0x1f80u)&~(0xe040u|63u))|((ftz&1?0x8000u:0)|(ftz&2?0x40u:0)));
  row_normal_cases=0;
  unsigned probe_csr=_mm_getcsr();unsigned short before_probe,after_probe;
  __asm__ volatile("fnstsw %0":"=am"(before_probe));
  row_normal_math=kls_snb_row_product_math();assert(row_normal_math>=0);
  __asm__ volatile("fnstsw %0":"=am"(after_probe));
  assert(_mm_getcsr()==probe_csr&&(before_probe&63u)==(after_probe&63u));
  for(int H=1;H<=65;++H) {
    double source[16*65],column[16*65];
    for(int i=0;i<16*H;++i){uint64_t bits=UINT64_C(0x7ff0000000000000)^(UINT64_C(0x123456789abcdef)*(unsigned)i);memcpy(source+i,&bits,8);}
    unsigned flags=_mm_getcsr();kls_rpanel_to_columns(source,H,column);
    assert(flags==_mm_getcsr());
    double roundtrip[16*65];kls_rpanel_from_columns(column,H,roundtrip);
    assert(flags==_mm_getcsr());assert(!memcmp(source,roundtrip,(size_t)16*H*sizeof(double)));
    for(int c=0;c<16;++c)for(int r=0;r<H;++r)assert(!memcmp(column+(size_t)c*H+r,source+kls_rpanel_at(H,c,r),8));
  }
  const int counts[]={1,2,7,8,16,17,32},rows[]={0,1,2,3,4,7,16,33,95};
  for(int layout=0;layout<2;++layout)for(unsigned c=0;c<sizeof(counts)/sizeof(*counts);++c)
    for(int gap=0;gap<2;++gap)if((counts[c]-1)*(gap?2:1)<(layout?16:32))
      for(unsigned r=0;r<sizeof(rows)/sizeof(*rows);++r)
        for(int holes=0;holes<3;++holes) {
          row_normal_edge(counts[c],rows[r],0,holes,layout,gap);
          for(int mask=8;mask<24;++mask)row_normal_edge(counts[c],rows[r],mask,holes,layout,gap);
        }
  const int widths[]={1,2,3,4,7,8,15,17,31,32};
  const int small_counts[]={1,2,7,8,16};
  for(unsigned w=0;w<sizeof(widths)/sizeof(*widths);++w)
    for(unsigned c=0;c<sizeof(small_counts)/sizeof(*small_counts);++c)
      for(unsigned r=0;r<sizeof(rows)/sizeof(*rows);++r)
        for(int gap=0;gap<2;++gap)if((small_counts[c]-1)*(gap?2:1)<16)
          for(int holes=0;holes<3;++holes)row_normal_small(widths[w],rows[r],small_counts[c],holes,gap);
  for(int H=19;H<=40;++H)row_normal_getrf(H,0,-1);
  assert(row_normal_cases==12091);
  fprintf(stderr,"ROW_NORMAL_EXACT_CASES %u contraction=%d\n",row_normal_cases,row_normal_math);
  assert(!fesetenv(&saved));_mm_setcsr(original);
#endif
}

static void kls_snb_row_normal_checks(void){for(int mode=0;mode<4;++mode)kls_snb_row_normal_checks_mode(mode);}
