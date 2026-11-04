/* Complete-edge comparison against the exact pre-prototype function body. */
static void register_trsm_checks(void) {
#if defined(__x86_64__) && FLT_EVAL_METHOD==0
  enum {SIZE=5120};double actual[SIZE],expected[SIZE],scratch[65*32];
  int32_t sub[16],dest[65];
  const int counts[]={7,8,15,16},widths[]={3,4,5,8,9},rows[]={0,1,7,8,9,15,16,17,33,65};
  const uint64_t exceptional[]={1,UINT64_C(0x0010000000000000),
    UINT64_C(0x7fefffffffffffff),UINT64_C(0x7ff0000000000000),
    UINT64_C(0x7ff8000000000042),UINT64_C(0x7ff0000000000042),UINT64_C(0x8000000000000000)};
  fenv_t saved;assert(!fegetenv(&saved));assert(!fesetround(FE_TONEAREST));unsigned original=_mm_getcsr();
  unsigned long long cases=0;
  for(int mode=0;mode<4;++mode)
    for(unsigned ww=0;ww<5;++ww)for(unsigned nn=0;nn<4;++nn)
      for(int shape=0;shape<3;++shape)for(unsigned rr=0;rr<10;++rr)
        for(int pattern=0;pattern<28;++pattern) {
          const int w=widths[ww],sc=counts[nn],ph=rows[rr];
          kls_snb_snode nodes[2]={{0}};
          nodes[0].width=32;nodes[0].height=104;nodes[0].uw=3;nodes[0].below=ph;
          nodes[1].width=w;nodes[1].height=112;nodes[1].uw=32;nodes[1].panel=3584;
          kls_snb_edge edge={0};edge.sub_count=sc;edge.u_start=1;
          kls_snb_block block={0};block.snodes=nodes;block.subcols=sub;block.edge_dsts=dest;
          for(int i=0;i<sc;++i)sub[i]=shape==2?2*i:i+(shape?3:0);
          for(int r=0;r<ph;++r)dest[r]=pattern%3==1 && r%3==0?-1:31+r;
          for(int i=0;i<SIZE;++i) {
            actual[i]=.00390625*(1+i%13);
            if(pattern>0 && pattern<8)memcpy(actual+i,exceptional+(i+pattern-1)%7,8);
            if(pattern==27)actual[i]=INFINITY;
          }
          if(pattern>=8 && pattern<24) {
            const int mask=pattern-8;
            for(int c=0;c<w;++c)for(int i=0;i<sc;++i)
              actual[3584+c*112+1+i]=mask&(1<<(c%4))?.125:-0.0;
          }
          if(pattern>=24)for(int c=0;c<w;++c)for(int i=0;i<sc;++i) {
            double v=0.;
            if(pattern==24 && i==sc-1)memcpy(&v,exceptional+5,8);
            if(pattern==25 && i==sc-1)v=DBL_TRUE_MIN;
            if(pattern==26 && i==sc-1)v=DBL_MIN;
            actual[3584+c*112+1+i]=v;
          }
          memcpy(expected,actual,sizeof actual);
          unsigned control=(original|0x1f80u)&~(0x3fu|0x8040u);
          control|=(mode&1?0x8000u:0)|(mode&2?0x40u:0);
          _mm_setcsr(control);feclearexcept(FE_ALL_EXCEPT);
          block.panels=expected;reference_register_edge(&block,nodes+1,&edge,8,scratch);
          unsigned flags=_mm_getcsr();int exceptions=fetestexcept(FE_ALL_EXCEPT);
          _mm_setcsr(control);feclearexcept(FE_ALL_EXCEPT);
          block.panels=actual;kls_snb_update_edge(&block,nodes+1,&edge,8,scratch);
          if(_mm_getcsr()!=flags || fetestexcept(FE_ALL_EXCEPT)!=exceptions) {
            fprintf(stderr,"register FP mismatch mode=%d w=%d sc=%d shape=%d ph=%d pattern=%d expected=%x actual=%x\n",mode,w,sc,shape,ph,pattern,flags,_mm_getcsr());abort();
          }
          for(int i=0;i<SIZE;++i)if(!isnan(expected[i]) || !isnan(actual[i]))
            if(memcmp(expected+i,actual+i,8)) {
              fprintf(stderr,"register value mismatch w=%d sc=%d shape=%d ph=%d pattern=%d at=%d\n",w,sc,shape,ph,pattern,i);abort();
            }
          ++cases;
        }
  assert(cases==67200);assert(!fesetenv(&saved));
  fprintf(stderr,"REGISTER_TRSM_ORACLE_CASES %llu\n",cases);
#endif
}
