/* Private floating LLL proposal generator, never an accuracy certificate. */
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* A changed basis column can affect only itself and later columns. Preserve
 * the untouched prefix; retain the original two-pass arithmetic in the suffix. */
static int orthogonalize(int m,int n,int first,int end,const long double *basis,
                        long double *orth,long double *mu,long double *norm) {
  for(int i=first;i<end;++i) {
    memset(mu+(size_t)i*n,0,(size_t)n*sizeof(*mu));
    long double *v=orth+(size_t)i*m;
    memcpy(v,basis+(size_t)i*m,(size_t)m*sizeof(*v));
    for(int pass=0;pass<2;++pass) for(int j=0;j<i;++j) {
      long double dot=0;
      for(int r=0;r<m;++r) dot+=v[r]*orth[(size_t)j*m+r];
      if(!(norm[j]>0) || !isfinite(dot)) return 0;
      long double q=dot/norm[j];mu[(size_t)i*n+j]+=q;
      for(int r=0;r<m;++r) v[r]-=q*orth[(size_t)j*m+r];
    }
    norm[i]=0;
    for(int r=0;r<m;++r) norm[i]+=v[r]*v[r];
    if(!(norm[i]>0) || !isfinite(norm[i])) return 0;
  }
  return 1;
}

int kls_accuracy_lattice_propose(int m,int n,const double *matrix,const double *target,
                    double *answer,int max_steps) {
  if(m<n || n<=0 || n>64 || max_steps<=0 || !matrix || !target || !answer ||
     (size_t)m>SIZE_MAX/(size_t)n/sizeof(long double)) return 0;
  long double *basis=calloc((size_t)m*n,sizeof(*basis));
  long double *orth=calloc((size_t)m*n,sizeof(*orth));
  long double *mu=calloc((size_t)n*n,sizeof(*mu));
  long double *transform=calloc((size_t)n*n,sizeof(*transform));
  long double *norm=calloc((size_t)n,sizeof(*norm));
  long double *work=calloc((size_t)m,sizeof(*work));
  long double *result=calloc((size_t)n,sizeof(*result));
  int ok=0;
  /* Distinguish transient workspace failure from deterministic no-proposal
   * outcomes, so the scheduler cannot mistake allocation failure for a fixed
   * point. The answer is untouched on either unsuccessful outcome. */
  if(!basis||!orth||!mu||!transform||!norm||!work||!result) { ok=-1;goto done; }
  for(int r=0;r<m;++r) {
    if(!isfinite(target[r])) goto done;
    work[r]=target[r];
    for(int j=0;j<n;++j) {
      double value=matrix[(size_t)r*n+j];
      if(!isfinite(value)) goto done;
      basis[(size_t)j*m+r]=value;
    }
  }
  for(int j=0;j<n;++j) transform[(size_t)j*n+j]=1;
  if(!orthogonalize(m,n,0,n,basis,orth,mu,norm)) goto done;
  /* Earliest changed column since the last rebuild; n means none. Analytic
   * geometry updates below do not change the original periodic refresh cadence. */
  int dirty=n;
  int k=1,steps=0,swaps=0;
  while(k<n) {
    if(++steps>max_steps) goto done;
    for(int j=k-1;j>=0;--j) {
      long double q=roundl(mu[(size_t)k*n+j]);
      if(!isfinite(q) || fabsl(q)>0x1p52L) goto done;
      if(q==0) continue;
      if(k<dirty) dirty=k;
      for(int r=0;r<m;++r) basis[(size_t)k*m+r]-=q*basis[(size_t)j*m+r];
      for(int r=0;r<n;++r) {
        transform[(size_t)k*n+r]-=q*transform[(size_t)j*n+r];
        if(fabsl(transform[(size_t)k*n+r])>0x1p52L) goto done;
      }
      for(int r=0;r<j;++r) mu[(size_t)k*n+r]-=q*mu[(size_t)j*n+r];
      mu[(size_t)k*n+j]-=q;
    }
    long double q=mu[(size_t)k*n+k-1];
    if(norm[k]>=(.75L-q*q)*norm[k-1]) ++k;
    else {
      if(k-1<dirty) dirty=k-1;
      long double alpha=norm[k]+q*q*norm[k-1];
      if(!(alpha>0) || !isfinite(alpha)) goto done;
      long double replacement=q*norm[k-1]/alpha;
      norm[k]=norm[k]*norm[k-1]/alpha;norm[k-1]=alpha;
      for(int r=0;r<k-1;++r) {
        long double swap=mu[(size_t)k*n+r];
        mu[(size_t)k*n+r]=mu[(size_t)(k-1)*n+r];mu[(size_t)(k-1)*n+r]=swap;
      }
      mu[(size_t)k*n+k-1]=replacement;
      for(int r=k+1;r<n;++r) {
        long double old=mu[(size_t)r*n+k];
        mu[(size_t)r*n+k]=mu[(size_t)r*n+k-1]-q*old;
        mu[(size_t)r*n+k-1]=old+replacement*mu[(size_t)r*n+k];
      }
      for(int r=0;r<m;++r) {
        long double swap=basis[(size_t)k*m+r];
        basis[(size_t)k*m+r]=basis[(size_t)(k-1)*m+r];basis[(size_t)(k-1)*m+r]=swap;
      }
      for(int r=0;r<n;++r) {
        long double swap=transform[(size_t)k*n+r];
        transform[(size_t)k*n+r]=transform[(size_t)(k-1)*n+r];transform[(size_t)(k-1)*n+r]=swap;
      }
      /* Periodically rebuild the approximate geometry; it never certifies x. */
      if(++swaps%n==0) {
        if(!orthogonalize(m,n,dirty,n,basis,orth,mu,norm)) goto done;
        dirty=n;
      }
      if(k>1) --k;
    }
  }
  if(!orthogonalize(m,n,dirty,n,basis,orth,mu,norm)) goto done;
  for(int j=n-1;j>=0;--j) {
    long double dot=0;
    for(int r=0;r<m;++r) dot+=work[r]*orth[(size_t)j*m+r];
    long double z=roundl(dot/norm[j]);
    if(!isfinite(z) || fabsl(z)>0x1p52L) goto done;
    for(int r=0;r<m;++r) work[r]-=z*basis[(size_t)j*m+r];
    for(int r=0;r<n;++r) result[r]+=z*transform[(size_t)j*n+r];
  }
  for(int j=0;j<n;++j) if(!isfinite(result[j]) || fabsl(result[j])>0x1p52L) goto done;
  for(int j=0;j<n;++j) answer[j]=(double)result[j];
  ok=1;
done:
  free(basis);free(orth);free(mu);free(transform);free(norm);free(work);free(result);
  return ok;
}
