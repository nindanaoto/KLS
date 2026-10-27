/* Private verdict/correction path only. Existing detailed-bound interfaces
 * retain their original helpers and results. Every row still computes its
 * complete residual/magnitude/error enclosure from current A,b,x. This
 * helper returns a tolerance witness, not the old tight quotient. No
 * verdict or coefficient value is retained across calls. The caller consumes only a verdict on success; detailed public bounds
 * still use the original range. */
static int componentwise_binary64_budget_row_verdict_witness(
    const kls_componentwise_plan *plan,size_t row,const int64_t *p,
    const int64_t *rows,const double *a,const double *b,const double *x,
    int transpose,long double tolerance,componentwise_tolerance64_cache *cache,
    double *upper) {
  double residual=b[row],magnitude=fabs(residual);
  if(!isfinite(residual) || !isfinite(x[row])) return 0;
  const size_t first=transpose?(size_t)p[row]:plan->rowptr[row];
  const size_t last=transpose?(size_t)p[row+1]:plan->rowptr[row+1];
  for(size_t k=first;k<last;++k) {
    const size_t at=transpose?k:(size_t)plan->positions[k];
    const size_t col=transpose?(size_t)rows[k]:(size_t)plan->columns[k];
    const double product=componentwise_binary64_mul(a[at],x[col]);
    residual=componentwise_binary64_sub(residual,product);
    magnitude=componentwise_binary64_add(magnitude,fabs(product));
  }
  /* Any nonfinite operand/product or overflowing accumulation propagates
   * into these final results. A zero rounded magnitude may conceal underflow
   * in a product; the zero case requires an exact-product check below. */
  if(!isfinite(residual) || !isfinite(magnitude)) return 0;
  if(magnitude==0) {
    /* Rounded zero products may conceal underflow. Recheck only this case
     * to establish that each exact product is zero, retaining the ordinary
     * zero-vector fast path without adding tests to nonzero row loops. */
    for(size_t k=first;k<last;++k) {
      const size_t at=transpose?k:(size_t)plan->positions[k];
      const size_t col=transpose?(size_t)rows[k]:(size_t)plan->columns[k];
      if(a[at]!=0 && x[col]!=0) return 0;
    }
    *upper=0;return residual==0;
  }
  const double coefficient=(double)(int64_t)((last-first)+1)*0x1p-50;
  if(!(coefficient<=0.25)) return 0;
  /* count <= 2^48-1. count * 2^-1074 is exactly representable. The same
   * error bounds both the true residual and the true denominator's error.
   * The rounded relative product is enclosed by an adjacent upward value. */
  const double relative_error=componentwise_positive_up(
      componentwise_binary64_mul(coefficient,magnitude));
  /* Every adjacent positive binary64 step is at least 2^-1074. Advancing
   * the encoding by the exact integer count therefore adds at least the
   * required count*eta budget without arithmetic on subnormal operands.
   * The result may be more conservative than adding count*eta directly. */
  uint64_t error_bits;memcpy(&error_bits,&relative_error,sizeof(error_bits));
  const uint64_t steps=(uint64_t)(last-first);
  if(error_bits>=UINT64_C(0x7ff0000000000000) ||
     steps>=UINT64_C(0x7ff0000000000000)-error_bits) return 0;
  error_bits+=steps;
  double error;memcpy(&error,&error_bits,sizeof(error));
  const double numerator=componentwise_positive_up(
      componentwise_binary64_add(fabs(residual),error));
  if(!(magnitude>error) || !isfinite(numerator)) return 0;
  const double denominator=componentwise_positive_down(
      componentwise_binary64_sub(magnitude,error));
  if(!(denominator>0)) return -1;
  /* Initialize only at the original compiled first-use point. The same
   * masked conversion result and its sticky flags persist for this range. */
  if(!cache->ready) {
    cache->value=componentwise_positive_down((double)tolerance);
    cache->ready=1;
  }
  const double tolerance64=cache->value;
  if(!(tolerance64>0) || !isfinite(tolerance64)) return 0;
  const double threshold=componentwise_positive_down(
      componentwise_binary64_mul(tolerance64,denominator));
  /* Affirmative ordered witness: threshold <= tolerance64*denominator,
   * so numerator <= threshold proves the same row error tolerance without
   * computing an unused tight quotient. Unordered comparisons decline. */
  if(!(numerator<=threshold)) return 0;
  *upper=tolerance64;return 1;
}

static int componentwise_plan_range_verdict_exact_domain(const kls_componentwise_plan *plan,
    size_t first,size_t last,const int64_t *p,const int64_t *rows,
    const double *a,const double *b,const double *x,int transpose,
    int binary64_permitted,long double tolerance,long double *upper) {
  unsigned short control;__asm__ volatile("fnstcw %0":"=m"(control));
  const unsigned csr=_mm_getcsr();
  if(!binary64_permitted || (control&63u)!=63u ||
     (csr&0x1f80u)!=0x1f80u || (csr&((1u<<15)|(1u<<6)|(3u<<13))) ||
     (control&0x0f00u)!=0x0300u || first==last)
    return componentwise_plan_range_original(plan,first,last,p,rows,a,b,x,
        transpose,binary64_permitted,tolerance,upper);
  componentwise_tolerance64_cache cache={0};
  double binary64_maximum=0;
  long double extended_maximum=0;
  for(size_t row=first;row<last;++row) {
    double bound;
    const int status=componentwise_binary64_budget_row_verdict_witness(
        plan,row,p,rows,a,b,x,transpose,tolerance,&cache,&bound);
    /* Original compiled cast placement can vary before a failed denominator
     * guard. Replay its unchanged range under the same masked environment. */
    if(status<0) return componentwise_plan_range_original(plan,first,last,p,
        rows,a,b,x,transpose,binary64_permitted,tolerance,upper);
    if(status>0) {
      if(bound>binary64_maximum) binary64_maximum=bound;
    } else {
      long double extended_bound;
      if(!componentwise_plan_row(plan,row,p,rows,a,b,x,transpose,0,tolerance,
          &extended_bound)) return 0;
      if(extended_bound>extended_maximum) extended_maximum=extended_bound;
    }
  }
  const long double binary64_extended=(long double)binary64_maximum;
  *upper=extended_maximum>binary64_extended?extended_maximum:binary64_extended;
  return 1;
}
