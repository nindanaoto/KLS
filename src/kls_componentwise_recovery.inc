/* At most two scaled proposals after an observed componentwise failure. The incumbent
 * and its caller-frame snapshot stay intact even if analyze/factor/solve fails.
 * A proposal never recursively proposes another factor or changes pivot rules. */
static int private_componentwise_scaled_recover(kls_solver *s,int transpose,
    int64_t nrhs,const double *b,double *x,int64_t ldx) {
  const private_original_operator *a=s->private_original;
  if(s->private_cold_trial || s->options.scale!=-2 || s->common.scale>0 ||
     !s->private_original_valid || !a || a->n!=(size_t)s->n)
    return KLS_ERR_SOLVE_FAILED;
  /* Reuse a factor only after the incumbent fails on this actual RHS. A
   * certificate for a previous RHS must never displace a usable incumbent.
   * The cached factor's public solve performs its own bounded correction. */
  if(s->private_scaled_recovery) {
    int status=transpose
      ?kls_solve_transpose(s->private_scaled_recovery,nrhs,b,0,x,ldx)
      :kls_solve(s->private_scaled_recovery,nrhs,b,0,x,ldx);
    if(status==KLS_OK)
      status=private_componentwise_finish(s,transpose,nrhs,b,x,ldx,status);
    return status;
  }
  const int scales[]={2,1};
  for(size_t attempt=0;attempt<sizeof(scales)/sizeof(scales[0]);++attempt) {
    kls_solver *trial=NULL;
    int status=kls_create(&trial);
    if(status!=KLS_OK) return status;
    trial->private_cold_trial=1;
    kls_options options=s->options;
    options.scale=scales[attempt];
    options.tuning_profile_path=NULL;
    status=kls_analyze_csc(trial,KLS_INDEX_INT64,(int64_t)a->n,
        a->p,a->rows,0,&options);
    if(status==KLS_OK) {
      /* The parent owns already loaded tuning; its options retain no borrowed
       * profile path. Preserve those numeric controls in the isolated factor. */
      trial->tuning=s->tuning;
      trial->snb_cost=s->snb_cost;
      status=kls_factor(trial,a->values);
    }
    if(status==KLS_OK)
      status=transpose?kls_solve_transpose(trial,nrhs,b,0,x,ldx)
                      :kls_solve(trial,nrhs,b,0,x,ldx);
    /* Independently check every RHS against the parent's current original
     * operator before retaining the factor, including CSR/permuted inputs. */
    if(status==KLS_OK)
      status=private_componentwise_finish(s,transpose,nrhs,b,x,ldx,status);
    if(status==KLS_OK) {s->private_scaled_recovery=trial;return KLS_OK;}
    kls_destroy(trial);
    if(status==KLS_ERR_OUT_OF_MEMORY) return status;
  }
  return KLS_ERR_SOLVE_FAILED;
}
