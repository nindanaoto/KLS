static int btf_descriptor_geometry(const kls_solver *s) {
  const kls_btf_l_descriptor *d=s->compact_btf_l_descriptors;
  if(!d)return 0;
  for(UF_long k=0;k<s->n;++k) {
    if(d[k].values!=s->i16solve_lx[k] || d[k].offset!=s->i16solve_loff[k] ||
       d[k].count!=s->i16solve_loff[k+1]-s->i16solve_loff[k])return 0;
    for(unsigned j=0;j<3 && j<d[k].count;++j)
      if(d[k].short_rows[j]!=s->i16solve_l[d[k].offset+j])return 0;
  }
  return 1;
}
static int btf_descriptor_checks(kls_solver *s,kls_solver *reference,double *a) {
  if(!btf_descriptor_geometry(s))return 0;
  kls_invalidate_i32_solve(s);
  if(s->compact_btf_l_descriptors)return 0;
  descriptor_allocation_failure=1;
  if(!kls_i32_solve_ready(s) || descriptor_allocation_failure ||
     !kls_small_btf_cache_ready(s) || s->compact_btf_l_descriptors)return 0;
  memcpy(reference->numeric->Xwork,s->numeric->Xwork,
    (size_t)s->symbolic->maxblock*sizeof(double));
  const unsigned long long before=compact_btf_descriptor_calls;
  if(kls_native_btf_map32_refactor(reference,a)!=1 ||
     kls_compact_btf_map32_refactor(s,a)!=1 ||
     compact_btf_descriptor_calls!=before || !compact_btf_factors_equal(s,reference))return 0;
  kls_invalidate_i32_solve(s);
  if(s->compact_btf_l_descriptors || !kls_i32_solve_ready(s) ||
     !btf_descriptor_geometry(s))return 0;
  return 1;
}
