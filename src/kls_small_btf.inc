#include <fenv.h>
#ifdef __linux__
static uint64_t kls_execution_domain_llc_bytes(void);
#endif

/* Reuse the existing immutable i16 cache. Numerical coefficients remain in
   the ordinary KLU factor and are loaded anew on every solve. */
static int kls_small_btf_environment(void) {
#if defined(__x86_64__) && FLT_EVAL_METHOD==0
  unsigned short control;__asm__ volatile("fnstcw %0":"=m"(control));
  return (control&63u)==63u && (_mm_getcsr()&0x1f80u)==0x1f80u &&
    fegetround()==FE_TONEAREST &&
    (_mm_getcsr()&0x6000u)==0;
#else
  return 0;
#endif
}

static int kls_small_btf_cache_ready(const kls_solver *s) {
  return s->i32solve_state>0 && s->i16solve_l && s->i16solve_u &&
    s->i16solve_loff && s->i16solve_uoff && s->i16solve_lx && s->i16solve_ux &&
    s->i16solve_q && s->i16solve_pnum && s->i16solve_r &&
    s->i16solve_offp && s->i16solve_offi;
}

static int kls_small_btf_cache_capable(kls_solver *s) {
  /* An ineligible symbolic frame must not traverse numerical ownership state
     on every ordinary solve/refactor. This is rejection only; positive frames
     still pass the complete ownership, factor-count, budget and FP gates. */
  if(!s || s->n<=1 || s->n>4096 || s->nnz<=0 || s->nnz>UINT16_MAX ||
     s->orientation!=KLS_ORIENTATION_TRANSPOSE)return 0;
  if(!s->symbolic || !s->numeric ||
     s->options.backend!=KLS_BACKEND_AUTO || s->stats.tuning_profile_active ||
     s->accuracy_policy!=KLS_ACCURACY_COMPONENTWISE_BACKWARD_ERROR ||
     s->common.scale>0 || s->numeric->Rs || !s->numeric->Offp ||
     !s->numeric->Udiag || !s->numeric->Xwork || !s->numeric->LUbx ||
     s->snb || s->snb_group ||
     kls_pts_refactor_ready(s) || s->dense_tail_cols || s->padded_run_of ||
     s->row_perm || s->user_col_perm || s->row_scale || s->col_scale ||
     s->row_refactor_values_ready || s->row_refactor_values_dirty ||
     s->refactor_l_packed_valid || s->refactor_u_packed_valid ||
     s->numeric_is_predicted || s->numeric_needs_refinement || s->tight_tol_refine ||
     s->pivot_nudge_count || s->common.kls_perturb_count) return 0;
  /* Only the immutable geometry/budget verdict persists. Numerical ownership
     above and the executor's current floating environment are never cached. */
  const int policy=atomic_load_explicit(&s->compact_btf_policy,memory_order_relaxed);
  if(policy)return policy>0;
  if(s->symbolic->nblocks<=1 || s->symbolic->maxblock<=1 ||
     s->numeric->lnz<0 || s->numeric->unz<0 ||
     s->numeric->lnz>UINT16_MAX || s->numeric->unz>UINT16_MAX ||
     s->numeric->Offp[s->n]<0 || s->numeric->Offp[s->n]>UINT16_MAX) {
    atomic_store_explicit(&s->compact_btf_policy,-1,memory_order_relaxed);return 0;
  }
  if(!kls_small_btf_environment())return 0;
#ifdef __linux__
  /* Include simultaneous i32/i16 mirrors, pointer tables, and native input
     maps. A complete pre-existing cache and map require no additional memory. */
  const uint64_t bytes=256+160*(uint64_t)s->n+16*(uint64_t)s->nnz+
    8*((uint64_t)s->numeric->lnz+(uint64_t)s->numeric->unz)+
    4*(uint64_t)s->numeric->Offp[s->n]+8*(uint64_t)s->symbolic->nblocks;
  const uint64_t llc=kls_execution_domain_llc_bytes();
  const uint64_t fallback_bytes=kls_small_btf_fallback_bytes(&s->small_btf_fallback);
  const int allowed=llc && fallback_bytes<=UINT64_MAX-bytes &&
    bytes+fallback_bytes<=llc/32;
  atomic_store_explicit(&s->compact_btf_policy,allowed?1:-1,memory_order_relaxed);
  return allowed;
#else
  return 0;
#endif
}

/* The linked long-index one-RHS KLU walk uses separate multiply/subtract,
   without zero skips. Preserve that native column order and division.
   There is one solve of the actual RHS, followed by the existing refinement
   and original-operator certificate. No reciprocal substitution or fusion. */
#if defined(__GNUC__) && !defined(__clang__)
__attribute__((noinline,hot,aligned(64),optimize("fp-contract=off")))
#elif defined(__clang__)
__attribute__((noinline,hot,aligned(64)))
#endif
/* The settled caller has just proved current numeric/cache capability.
   Execution still checks the current FP environment, including mixed controls. */
static UF_long kls_compact_btf_transpose_prechecked(kls_solver *s,double *out) {
#if defined(__clang__)
#pragma clang fp contract(off)
#endif
  if(!kls_small_btf_environment())return 0; /* current execution environment */
  const UF_long n=s->n;
  double *work=s->numeric->Xwork;
  const double *diagonal=s->numeric->Udiag,*offx=s->numeric->Offx;
  const uint16_t *q=s->i16solve_q,*pnum=s->i16solve_pnum;
  const uint16_t *offp=s->i16solve_offp,*offi=s->i16solve_offi;
  const uint16_t *up=s->i16solve_uoff,*ui=s->i16solve_u;
  const uint16_t *lp=s->i16solve_loff,*li=s->i16solve_l;
  for(UF_long k=0;k<n;++k)work[k]=out[q[k]];
  for(UF_long block=0;block<s->symbolic->nblocks;++block) {
    const UF_long first=s->i16solve_r[block],last=s->i16solve_r[block+1];
    double *local=work+first;
    for(UF_long k=first;k<last;++k) {
      double value=work[k];
      for(uint32_t z=offp[k];z<offp[k+1];++z)value-=offx[z]*work[offi[z]];
      work[k]=value;
    }
    for(UF_long k=first;k<last;++k) {
      double value=work[k];const double *ux=s->i16solve_ux[k];
      const uint32_t begin=up[k],end=up[k+1];
      for(uint32_t z=begin;z<end;++z) {
        const double previous=local[ui[z]];
        value-=ux[z-begin]*previous;
      }
      work[k]=value/diagonal[k];
    }
    for(UF_long k=last;k-->first;) {
      double value=work[k];const double *lx=s->i16solve_lx[k];
      const uint32_t begin=lp[k],end=lp[k+1];
      for(uint32_t z=begin;z<end;++z) {
        const double previous=local[li[z]];
        value-=lx[z-begin]*previous;
      }
      work[k]=value;
    }
  }
  for(UF_long k=0;k<n;++k)out[pnum[k]]=work[k];
  s->common.status=TRILINOS_KLU_OK;
#ifdef KLS_COMPACT_BTF_SOLVE_TEST_EVENT
  KLS_COMPACT_BTF_SOLVE_TEST_EVENT();
#endif
#ifdef KLS_COMPACT_BTF_SOLVE_DIAGNOSTIC_EVENT
  KLS_COMPACT_BTF_SOLVE_DIAGNOSTIC_EVENT(s);
#endif
  return 1;
}

/* Preserve the complete admission and optional cache preparation for every
   raw/ordinary caller. No solve or answer is cached by either entry. */
static UF_long kls_compact_btf_transpose_solve(kls_solver *s,double *out) {
  if(!s || s->n<=1 || s->n>4096 || s->nnz<=0 || s->nnz>UINT16_MAX ||
     s->orientation!=KLS_ORIENTATION_TRANSPOSE)return 0;
  if(!kls_small_btf_environment() || !kls_small_btf_cache_capable(s) ||
     !kls_i32_solve_ready(s) || !kls_small_btf_cache_ready(s))return 0;
  return kls_compact_btf_transpose_prechecked(s,out);
}
