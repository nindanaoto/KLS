#define _GNU_SOURCE
#define _POSIX_C_SOURCE 200809L
#include <execinfo.h>
#define KLS_TEST_FIRST_LIMIT
#include KLS_CORE_SOURCE
static int64_t n,z,*p,*j;
static double *a,*b,*x;
int main(int argc,char**argv){
 if(argc!=4)return 2;FILE*f=fopen(argv[1],"rb");uint64_t h[4];kls_options o;
 if(!f||fread(h,sizeof(h),1,f)!=1||fread(&o,sizeof(o),1,f)!=1||h[0]!=UINT64_C(0x4b4c534449414731))return 3;
 n=h[1];z=h[2];p=malloc((n+1)*8);j=malloc(z*8);a=malloc(z*8);b=calloc(n,8);x=calloc(n,8);
 if(fread(p,8,n+1,f)!=(size_t)n+1||fread(j,8,z,f)!=(size_t)z||fread(a,8,z,f)!=(size_t)z)return 4;fclose(f);
 for(int64_t r=0;r<n;++r){long double v=0;for(int64_t k=p[r];k<p[r+1];++k)v+=(long double)a[k]*(1+(j[k]%7)*.125);b[r]=(double)v;}
 o.threads=atoi(argv[2]);int matched=atoi(argv[3]);kls_solver*s=NULL;kls_create(&s);
 int status=kls_analyze_csr(s,KLS_INDEX_INT64,n,p,j,0,&o);
 UF_long *perm=NULL,*cp=NULL,*ri=NULL,*map=NULL,mc=0;double *rv=NULL,*cs=NULL,*rs=NULL; kls_solver*t=NULL;
 if(status==0 && matched){
  double*values=NULL;status=prepare_numeric_values(s,a,&values);
  if(!status)status=build_spral_hungarian_row_match_scaling(n,z,s->col_ptr,s->row_idx,values,&perm,&mc,&rs,&cs);
  if(!status && mc!=(UF_long)n)status=KLS_ERR_UNSUPPORTED;
  if(!status)status=build_sorted_row_permuted_pattern(n,z,s->col_ptr,s->row_idx,values,perm,NULL,&cp,&ri,&rv,&map);
  if(!status){kls_options options=o;options.orientation=KLS_ORIENTATION_NORMAL;options.static_pivoting=0;kls_create(&t);status=kls_analyze_csc(t,KLS_INDEX_INT64,n,cp,ri,0,&options);}
 }
 kls_first_test_limit=21850128;
 if(!status)status=kls_factor(matched?t:s,matched?rv:a);
 kls_solver*fact=matched?t:s;
 if(!status && matched){double*rhs=malloc(n*8),*result=malloc(n*8);
  if(s->orientation==KLS_ORIENTATION_NORMAL){for(int64_t r=0;r<n;++r)rhs[perm[r]]=b[r];status=kls_solve(t,1,rhs,0,x,0);}
  else{status=kls_solve_transpose(t,1,b,0,result,0);if(!status)for(int64_t r=0;r<n;++r)x[r]=result[perm[r]];}
  free(rhs);free(result);
 }else if(!status)status=kls_solve(s,1,b,0,x,0);
 long double bound=INFINITY;int cert=kls_componentwise_certify(n,p,j,a,b,x,1,NULL,&bound);
 printf("probe threads=%d matched=%d orientation=%s status=%d certificate=%d bound=%.21Lg path=%s peak=%zu limit=%zu reason=%d row=%ld match_rows=%ld\n",o.threads,matched,kls_orientation_name(s->orientation),status,cert,bound,fact?kls_factor_path_name(fact->stats.last_factor_path):"none",fact?fact->subtree_last_peak:0,fact?fact->subtree_last_limit:0,fact?fact->subtree_last_failure:0,fact?fact->subtree_last_row:0,mc);
 kls_destroy(t);kls_destroy(s);free(perm);free(cp);free(ri);free(rv);free(map);free(rs);free(cs);free(p);free(j);free(a);free(b);free(x);
 return status!=KLS_OK || cert!=1;
}
