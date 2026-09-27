/* Conservative accumulated-error fast gate. A miss leaves the existing
 * interval and binary128 certificates authoritative. Requires the outward
 * arithmetic helpers from kls_strict_interval.inc. */
static int gamma_csc_certificate(size_t n,const int64_t*p,const int64_t*i,const double*a,
             const double*b,const double*x,int t,
             long double *upper_out,long double *budget_out) {
#if LDBL_MANT_DIG != 64 || LDBL_MAX_EXP != 16384 || !defined(__x86_64__)
  (void)n;(void)p;(void)i;(void)a;(void)b;(void)x;(void)t;
  (void)upper_out;(void)budget_out;
  return 0;
#else
  unsigned short cw;__asm__ volatile("fnstcw %0":"=m"(cw));
  if((cw&0x0f00u)!=0x0300u || fegetround()!=FE_TONEAREST || !p||!b||!x ||
     n>SIZE_MAX/(2*sizeof(long double)) || p[0]!=0 || p[n]<0 ||
     (p[n]&&(!a||!i))) return 0;
  long double*r=malloc((n?2*n:2)*sizeof(*r));
  uint64_t*degree=calloc(n?n:1,sizeof(*degree));
  if(!r||!degree){free(r);free(degree);return 0;}
  long double*s=r+n,bn=0,upper=0;int zero=1,ok=0;
  for(size_t j=0;j<n;++j){
    volatile long double bj=b[j],xj=x[j];
    if(!isfinite(bj)||!isfinite(xj))goto done;
    r[j]=bj;s[j]=fabsl(bj);zero&=bj==0;
    long double sq=interval_outward(bj*bj,-INFINITY);
    if(sq<0)sq=0;
    bn=interval_outward(bn+sq,-INFINITY);if(bn<0)bn=0;
  }
  for(size_t j=0;j<n;++j){
    if(p[j]<0||p[j]>p[j+1]||p[j+1]>p[n])goto done;
    for(int64_t k=p[j];k<p[j+1];++k){
      if(i[k]<0||(uint64_t)i[k]>=n)goto done;
      size_t dest=t?j:(size_t)i[k],src=t?(size_t)i[k]:j;
      volatile long double av=a[k],xv=x[src];
      /* The initial full-vector pass already validated x, including entries
       * unused by A. Retain the extended load for subnormal inputs. */
      if(!isfinite(av))goto done;
      long double v=av*xv;r[dest]-=v;s[dest]+=fabsl(v);++degree[dest];
    }
  }
  for(size_t j=0;j<n;++j){
    /* Products of binary64 values cannot underflow/overflow x87 extended.
     * At most 2*k rounded operations affect the residual or magnitude sum.
     * For c=4*(k+1)*epsilon <=1/4, c dominates
     * gamma_(2k)/(1-gamma_(2k)), including magnitude-sum error.
     * Outward rounding also encloses the evaluation of the bound itself. */
    long double c=4*((long double)degree[j]+1)*LDBL_EPSILON;
    if(!(c<=0.25L))goto done;
    long double error=interval_outward(c*s[j],INFINITY);
    long double bound=interval_outward(fabsl(r[j])+error,INFINITY);
    upper=interval_outward(upper+interval_outward(bound*bound,INFINITY),INFINITY);
  }
  long double tol=interval_outward(1e-8L,-INFINITY);
  long double budget=interval_outward(interval_outward(tol*tol,-INFINITY)*(zero?1:bn),-INFINITY);
  ok=isfinite(upper)&&budget>0&&upper<=budget;
  if(ok){if(upper_out)*upper_out=upper;if(budget_out)*budget_out=budget;}
done:free(r);free(degree);return ok;
#endif
}
