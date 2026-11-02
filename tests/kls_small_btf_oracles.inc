static int compact_btf_factors_equal(kls_solver *a,kls_solver *b) {
  if(a->n!=b->n || a->symbolic->nblocks!=b->symbolic->nblocks)return 0;
  if(memcmp(a->numeric->Udiag,b->numeric->Udiag,(size_t)a->n*sizeof(double)))return 0;
  if(memcmp(a->numeric->Offx,b->numeric->Offx,(size_t)a->numeric->Offp[a->n]*sizeof(double)))return 0;
  for(UF_long block=0;block<a->symbolic->nblocks;++block) {
    const UF_long first=a->symbolic->R[block],last=a->symbolic->R[block+1];
    if(last-first<2)continue;
    for(UF_long k=0;k<last-first;++k)for(int upper=0;upper<2;++upper) {
      UF_long *ai,*bi,an,bn;double *av,*bv;
      kls_klu_get_pointer(a->numeric->LUbx[block],
        (upper?a->numeric->Uip:a->numeric->Lip)+first,
        (upper?a->numeric->Ulen:a->numeric->Llen)+first,k,&ai,&av,&an);
      kls_klu_get_pointer(b->numeric->LUbx[block],
        (upper?b->numeric->Uip:b->numeric->Lip)+first,
        (upper?b->numeric->Ulen:b->numeric->Llen)+first,k,&bi,&bv,&bn);
      if(an!=bn || memcmp(ai,bi,(size_t)an*sizeof(*ai)) ||
         memcmp(av,bv,(size_t)an*sizeof(*av)))return 0;
    }
  }
  return 1;
}

static int compact_btf_solve_oracle(kls_solver *s) {
#if defined(__x86_64__)
  const size_t n=(size_t)s->n;int ok=1;
  double *vectors=malloc(3*n*sizeof(double));if(!vectors)return 0;
  double *rhs=vectors,*got=rhs+n,*expected=got+n;
  fenv_t environment;fegetenv(&environment);
  const unsigned original_mxcsr=_mm_getcsr();
  const int rounds[]={FE_TONEAREST,FE_DOWNWARD,FE_UPWARD,FE_TOWARDZERO};
  for(int rounding=0;rounding<4 && ok;++rounding)
    for(unsigned flush=0;flush<4 && ok;++flush)for(int seed=0;seed<32 && ok;++seed) {
      for(size_t i=0;i<n;++i) {
        uint64_t bits=(i+seed)%5==0?((seed&1)?UINT64_C(0x8000000000000000):0):
          (i+seed)%5==1?(uint64_t)(seed+1)*(i+1):
          (i+seed)%5==2?(UINT64_C(0x8000000000000000)|((uint64_t)(seed+3)*(i+5))):
          (i+seed)%5==3?(UINT64_C(0x3fd0000000000000)|((uint64_t)(seed+1)*(i+3)<<32)):
            (UINT64_C(0xbff8000000000000)|((uint64_t)(seed+3)*(i+5)<<29));
        memcpy(rhs+i,&bits,sizeof(bits));
      }
      fesetround(rounds[rounding]);feclearexcept(FE_ALL_EXCEPT);
      unsigned control=_mm_getcsr();
      control=(control&~(0x8040u|0x3fu))|((flush&1)?0x8000u:0)|((flush&2)?0x40u:0);
      memcpy(got,rhs,n*sizeof(double));_mm_setcsr(control);
      if(rounding!=0) {
        const unsigned before=_mm_getcsr();
        if(kls_compact_btf_transpose_prechecked(s,got) ||
           memcmp(got,rhs,n*sizeof(double)) || _mm_getcsr()!=before)ok=0;
      }
      const int used=kls_compact_btf_transpose_solve(s,got);
      ok=used==(rounding==0);
      if(!used) {
        ok=ok && !memcmp(got,rhs,n*sizeof(double)) && (_mm_getcsr()&0x3fu)==0;
        ok=ok && trilinos_klu_l_tsolve(s->symbolic,s->numeric,(UF_long)n,1,got,&s->common);
      }
      const unsigned got_flags=_mm_getcsr()&0x3fu;
      memcpy(expected,rhs,n*sizeof(double));_mm_setcsr(control);
      ok=ok && trilinos_klu_l_tsolve(s->symbolic,s->numeric,(UF_long)n,1,expected,&s->common);
      const unsigned expected_flags=_mm_getcsr()&0x3fu;
      if(memcmp(got,expected,n*sizeof(double)) || got_flags!=expected_flags) {
        fprintf(stderr,"compact BTF oracle mismatch round=%d flush=%u seed=%d flags=%u/%u\n",
          rounding,flush,seed,got_flags,expected_flags);ok=0;
      }
    }
  /* Mixed x87/SSE rounding and an unmasked trap also decline read-only. */
  fesetround(FE_TONEAREST);feclearexcept(FE_ALL_EXCEPT);
  const unsigned nearest=_mm_getcsr()&~0x3fu;
  for(int mode=0;mode<2 && ok;++mode) {
    const unsigned control=mode?(nearest&~0x80u):((nearest&~0x6000u)|0x2000u);
    memcpy(got,rhs,n*sizeof(double));_mm_setcsr(control);
    const unsigned prechecked_control=_mm_getcsr();
    if(kls_compact_btf_transpose_prechecked(s,got) ||
       memcmp(got,rhs,n*sizeof(double)) || _mm_getcsr()!=prechecked_control) {
      ok=0;_mm_setcsr(nearest);break;
    }

    const int used=kls_compact_btf_transpose_solve(s,got);
    ok=used==0 && !memcmp(got,rhs,n*sizeof(double)) && (_mm_getcsr()&0x3fu)==0;
    _mm_setcsr(nearest);
  }
  fesetenv(&environment);_mm_setcsr(original_mxcsr);free(vectors);
  if(!ok)return 0;
  /* A missing optional mirror declines before writing caller storage. */
  uint16_t *saved=s->i16solve_u;s->i16solve_u=NULL;
  double *out=malloc(n*sizeof(double));if(!out){s->i16solve_u=saved;return 0;}
  for(size_t i=0;i<n;++i)out[i]=123.;
  ok=kls_compact_btf_transpose_solve(s,out)==0;
  for(size_t i=0;i<n;++i)ok=ok && out[i]==123.;
  s->i16solve_u=saved;free(out);return ok;
#else
  (void)s;return 1;
#endif
}
