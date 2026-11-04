/* Immutable input positions and native coefficient offsets for the proved
   forest's remaining columns. Numerical values are never retained here. */
static void kls_prepare_forest_remaining_plan(kls_solver *s,
    kls_forest_leaf_plan *p,uint64_t llc) {
  if(!llc)return;
  const int32_t *rest=p->indices+6*p->count;
  const trilinos_klu_l_numeric *numeric=s->numeric;
  uintmax_t updates=0;
  for(UF_long at=0;at<p->remaining;++at)updates+=numeric->Ulen[rest[at]];
  const uintmax_t entries=9*(uintmax_t)p->remaining+3*updates;
  const long double bytes=sizeof(*p)+((long double)p->n+5*p->count+
      s->nnz+entries)*sizeof(int32_t);
  if(entries>SIZE_MAX/sizeof(int32_t) || entries>INT32_MAX || bytes>llc/32)return;
  int32_t *plan=malloc((size_t)(entries?entries:1)*sizeof(*plan));
  if(!plan)return;
  int32_t *deps=plan+9*p->remaining;
  UF_long next=0;
  const double *lu=numeric->LUbx[0];
  for(UF_long at=0;at<p->remaining;++at) {
    const UF_long k=(UF_long)rest[at],count=numeric->Ulen[k];
    const UF_long units=kls_klu_units_for_indices(count);
    const UF_long q=s->symbolic->Q[k];
    if(count>INT32_MAX || next>INT32_MAX-count || q<0 || q>=s->n ||
       (count && (numeric->Uip[k]<0 || numeric->Uip[k]>INT32_MAX-units ||
        numeric->Uip[k]+units>INT32_MAX-count)))goto decline;
    int32_t *row=plan+9*at;
    row[0]=(int32_t)k;row[1]=row[2]=row[3]=row[4]=-1;
    row[5]=(int32_t)next;row[6]=(int32_t)count;row[7]=(int32_t)q;
    row[8]=count?(int32_t)(numeric->Uip[k]+units):0;
    if(numeric->Llen[k]) {
      const UF_long parent=((const UF_long *)(lu+numeric->Lip[k]))[0];
      const UF_long lower_units=kls_klu_units_for_indices(1);
      if(numeric->Llen[k]!=1 || parent<=k || parent>=s->n ||
         numeric->Lip[k]<0 || numeric->Lip[k]>INT32_MAX-lower_units)goto decline;
      row[3]=(int32_t)parent;row[4]=(int32_t)(numeric->Lip[k]+lower_units);
    }
    const UF_long *ui=count?(const UF_long *)(lu+numeric->Uip[k]):NULL;
    for(UF_long up=0;up<count;++up) {
      const UF_long j=ui[up],lower_units=kls_klu_units_for_indices(1);
      if(j<0 || j>=k || numeric->Llen[j]!=1 || numeric->Lip[j]<0 ||
         numeric->Lip[j]>INT32_MAX-lower_units ||
         ((const UF_long *)(lu+numeric->Lip[j]))[0]!=k)goto decline;
      deps[3*(next+up)]=(int32_t)j;
      deps[3*(next+up)+1]=-1;
      deps[3*(next+up)+2]=(int32_t)(numeric->Lip[j]+lower_units);
    }
    for(UF_long pos=s->refactor_col_ptr[k];pos<s->refactor_col_ptr[k+1];++pos) {
      const int32_t slot=s->forest_input_slot[pos],input=s->refactor_input_pos32[pos];
      if(input<0 || (UF_long)input>=s->nnz)goto decline;
      const UF_long row_index=(UF_long)s->refactor_row_idx32[pos];
      if(slot>=0) {
        if((UF_long)slot>=count || deps[3*(next+slot)+1]>=0 ||
           row_index!=(UF_long)deps[3*(next+slot)])goto decline;
        deps[3*(next+slot)+1]=input;
      } else if(slot==-1 && row[1]<0 && row_index==k)row[1]=input;
      else if(slot==-2 && row[2]<0 && row_index==(UF_long)row[3])row[2]=input;
      else goto decline;
    }
    if(row[1]<0 || (row[3]>=0 && row[2]<0))goto decline;
    for(UF_long up=0;up<count;++up)if(deps[3*(next+up)+1]<0)goto decline;
    next+=count;
  }
  if((uintmax_t)next!=updates)goto decline;
  p->remaining_plan=plan;return;
decline:
  free(plan);
}

static KLS_ALWAYS_INLINE double kls_forest_planned_input(const double *values,
    int32_t input,int32_t row,const double *scale,const double *rs) {
  const double value=values[input];
  return !scale?value:(isfinite(scale[row])?value*scale[row]:value/rs[row]);
}

/* All dependencies retain their ascending native order. Each producer has
   been recomputed in this call before its L coefficient or forward value is
   read; every U value is published to the ordinary native factor stream. */
static int kls_forest_remaining_planned(kls_solver *s,const double *values,
    const double *rhs,const double *scale,double *y,double *diag_recip) {
#ifdef KLS_FOREST_REMAINING_TEST_EVENT
  KLS_FOREST_REMAINING_TEST_EVENT();
#endif
  const kls_forest_leaf_plan *p=s->forest_leaves;
  const int32_t *plan=p->remaining_plan,*deps=plan+9*p->remaining;
  const trilinos_klu_l_numeric *numeric=s->numeric;
  double *lu=numeric->LUbx[0],*udiag=numeric->Udiag;
  const double *rs=numeric->Rs;
  for(UF_long at=0;at<p->remaining;++at) {
    const int32_t *row=plan+9*at;
    const int32_t k=row[0],parent=row[3];
    double diagonal=kls_forest_planned_input(values,row[1],k,scale,rs);
    const double lower=parent<0?0:kls_forest_planned_input(values,row[2],parent,scale,rs);
    double yk=rhs[row[7]];
    double *ux=lu+row[8];
    for(int32_t up=0;up<row[6];++up) {
      const int32_t *dep=deps+3*((UF_long)row[5]+up);
      const int32_t j=dep[0];
      const double u=kls_forest_planned_input(values,dep[1],j,scale,rs);
      ux[up]=u;
#if defined(__FMA__)
      yk=fma(-u,y[j],yk);
#else
      yk-=u*y[j];
#endif
      if(u!=0)diagonal-=lu[dep[2]]*u;
    }
    if(diagonal==0) {
      s->common.status=TRILINOS_KLU_SINGULAR;
      if(s->common.numerical_rank==KLS_KLU_EMPTY) {
        s->common.numerical_rank=k;s->common.singular_col=row[7];
      }
      if(s->common.halt_if_singular)return 0;
    }
    udiag[k]=diagonal;const double inverse=1.0/diagonal;
    y[k]=isfinite(inverse)?yk*inverse:yk/diagonal;
    if(diag_recip)diag_recip[k]=inverse;
    if(parent>=0)lu[row[4]]=isfinite(inverse)?lower*inverse:lower/diagonal;
  }
  return 1;
}
