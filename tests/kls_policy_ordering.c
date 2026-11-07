/* Verify that a declined consultation publishes the actual AMD allocation. */
#define _GNU_SOURCE
static int block_probes, amd_calls;
#define KLS_POLICY_ANALYZE_EVENT(order) (amd_calls += ((order) == KLS_ORDERING_AMD))
#define KLS_POLICY_BLOCK_EVENT() (++block_probes)
static void record_policy_seed(const void *p);
static void record_policy_reuse(const void *p);
#define KLS_POLICY_SEED_EVENT(p) record_policy_seed(p)
#define KLS_POLICY_REUSE_EVENT(p) record_policy_reuse(p)
#include "../src/kls.c"
#undef NDEBUG
#include <assert.h>
static int seed_count,reuse_count;
static const void *policy_seed;
static void record_policy_seed(const void *p){policy_seed=p;++seed_count;}
static void record_policy_reuse(const void *p){assert(policy_seed==p);++reuse_count;}

static UF_long policy_gcd(UF_long a,UF_long b){while(b){UF_long r=a%b;a=b;b=r;}return a;}
static void boundary_and_relabel_checks(void) {
  const UF_long sizes[]={32767,32768,32769,131071,131072,131073};
  for(size_t which=0;which<sizeof sizes/sizeof sizes[0];++which)
    for(int relabel=0;relabel<2;++relabel) {
      const UF_long n=sizes[which],nz=3*n-2;
      UF_long *ptr=calloc((size_t)n+1,sizeof *ptr),*idx=malloc((size_t)nz*sizeof *idx);
      UF_long *permutation=malloc((size_t)n*sizeof *permutation);
      double *a=malloc((size_t)nz*sizeof *a),*b=calloc((size_t)n,sizeof *b),*x=malloc((size_t)n*sizeof *x);
      assert(ptr&&idx&&permutation&&a&&b&&x);
      UF_long multiplier=17;while(policy_gcd(multiplier,n)!=1)multiplier+=2;
      for(UF_long c=0;c<n;++c)permutation[c]=relabel?(multiplier*c+11)%n:c;
      for(UF_long c=0;c<n;++c)ptr[permutation[c]+1]=1+(c>0)+(c+1<n);
      for(UF_long c=0;c<n;++c)ptr[c+1]+=ptr[c];
      for(UF_long c=0;c<n;++c){
        UF_long rows[3],count=0;double values[3];
        if(c){rows[count]=permutation[c-1];values[count++]=-1;}
        rows[count]=permutation[c];values[count++]=4;
        if(c+1<n){rows[count]=permutation[c+1];values[count++]=-1;}
        for(UF_long i=1;i<count;++i){UF_long row=rows[i],j=i;double value=values[i];
          while(j&&rows[j-1]>row){rows[j]=rows[j-1];values[j]=values[j-1];--j;}rows[j]=row;values[j]=value;}
        for(UF_long i=0;i<count;++i){UF_long z=ptr[permutation[c]]+i;idx[z]=rows[i];a[z]=values[i];}
      }
      kls_solver *solver=NULL;assert(kls_create(&solver)==KLS_OK);
      kls_options options;kls_default_options(&options);options.orientation=KLS_ORIENTATION_NORMAL;
      seed_count=reuse_count=block_probes=amd_calls=0;
      assert(kls_analyze_csc(solver,KLS_INDEX_INT64,n,ptr,idx,0,&options)==KLS_OK);
      assert(seed_count==1 && reuse_count==1 && block_probes==1 && amd_calls==1);
      assert(kls_factor(solver,a)==KLS_OK);
      for(UF_long z=0;z<nz;++z)a[z]*=1.03125;
      for(UF_long c=0;c<n;++c)for(UF_long z=ptr[permutation[c]];z<ptr[permutation[c]+1];++z)
        b[idx[z]]+=a[z]*(.25*(double)(c%9)-1);
      assert(kls_refactor_solve(solver,a,1,b,0,x,0)==KLS_OK);
      for(UF_long c=0;c<n;++c)assert(fabs(x[permutation[c]]-(.25*(double)(c%9)-1))<1e-10);
      kls_destroy(solver);free(ptr);free(idx);free(permutation);free(a);free(b);free(x);
    }
}

static void concrete_block_checks(void) {
  const UF_long sizes[] = {301, 9999, 10000, 10001};
  for (size_t trial = 0; trial < sizeof sizes / sizeof sizes[0]; ++trial) {
    const UF_long n = sizes[trial], interior = n - 1;
    UF_long *ptr = calloc((size_t)n + 1, sizeof(*ptr));
    UF_long *idx = malloc((size_t)5 * n * sizeof(*idx));
    double *a = malloc((size_t)5 * n * sizeof(*a));
    assert(ptr && idx && a);
    UF_long nz = 0;
    for (UF_long j = 0; j < n; ++j) {
      ptr[j] = nz;
      if (j == interior) {
        for (UF_long i = 0; i < interior; i += 5) {
          idx[nz] = i; a[nz++] = -0.125;
        }
      } else {
        const UF_long block = j * 3 / interior;
        if (j && (j - 1) * 3 / interior == block) {
          idx[nz] = j - 1; a[nz++] = -0.125;
        }
        idx[nz] = j; a[nz++] = 4;
        if (j + 1 < interior && (j + 1) * 3 / interior == block) {
          idx[nz] = j + 1; a[nz++] = -0.125;
        }
        if (j % 5 == 0) { idx[nz] = interior; a[nz++] = -0.125; }
        continue;
      }
      idx[nz] = j; a[nz++] = 4;
    }
    ptr[n] = nz;
    UF_long *perm = NULL, *comp = NULL;
    assert(kls_build_block_structured_order(n, ptr, idx, &perm, &comp));
    unsigned char *seen = calloc((size_t)n, 1);
    assert(seen && perm[n - 1] == interior);
    for (UF_long k = 0; k < n; ++k) {
      assert(perm[k] < n && !seen[perm[k]]); seen[perm[k]] = 1;
    }
    free(seen); free(perm); free(comp);
    kls_solver *solver = NULL; assert(kls_create(&solver) == KLS_OK);
    kls_options o; kls_default_options(&o); o.orientation = KLS_ORIENTATION_NORMAL;
    o.expected_refactorizations = 100; o.expected_solves = 101;
    seed_count = reuse_count = block_probes = 0;
    assert(kls_analyze_csc(solver, KLS_INDEX_INT64, n, ptr, idx, 0, &o) == KLS_OK);
#ifdef KLS_HAVE_SPRAL_SCALING
    assert(seed_count == 1);
    if (n >= 10000) {
      assert(block_probes == 1 && solver->block_order_checked && solver->block_order_perm);
    } else {
      assert(block_probes == 0 && !solver->block_order_perm);
    }
#endif
    assert(kls_factor(solver, a) == KLS_OK);
    assert(block_probes <= 1);
    kls_destroy(solver); free(ptr); free(idx); free(a);
  }
}

static void matched_representation_checks(void) {
  const UF_long n = 16385;
  UF_long *ptr = malloc(((size_t)n + 1) * sizeof *ptr);
  UF_long *idx = malloc((size_t)2 * n * sizeof *idx);
  double *a = malloc((size_t)2 * n * sizeof *a);
  double *b = calloc((size_t)n, sizeof *b), *x = malloc((size_t)n * sizeof *x);
  assert(ptr && idx && a && b && x);
  for (UF_long c = 0; c < n; ++c) {
    ptr[c] = 2 * c;
    UF_long first = (c + 1) % n, second = (c + 2) % n;
    idx[2 * c] = first < second ? first : second;
    idx[2 * c + 1] = first < second ? second : first;
    a[2 * c] = first < second ? 4 : -0.125;
    a[2 * c + 1] = first < second ? -0.125 : 4;
  }
  ptr[n] = 2 * n;
  double fill[2];
  for (int btf = 0; btf < 2; ++btf) {
    kls_solver *reference = NULL; assert(kls_create(&reference) == KLS_OK);
    kls_options o; kls_default_options(&o); o.ordering = KLS_ORDERING_AMD;
    o.orientation = KLS_ORIENTATION_NORMAL; o.static_pivoting = 0; o.use_btf = btf;
    assert(kls_analyze_csc(reference, KLS_INDEX_INT64, n, ptr, idx, 0, &o) == KLS_OK);
    fill[btf] = symbolic_score(reference->symbolic);
    kls_destroy(reference);
  }
  assert(fill[1] < fill[0]);
  kls_solver *solver = NULL; assert(kls_create(&solver) == KLS_OK);
  kls_options o; kls_default_options(&o); o.orientation = KLS_ORIENTATION_NORMAL;
  o.static_pivoting = 0;
  seed_count = reuse_count = block_probes = amd_calls = 0;
  assert(kls_analyze_csc(solver, KLS_INDEX_INT64, n, ptr, idx, 0, &o) == KLS_OK);
  assert(solver->symbolic->do_btf && symbolic_score(solver->symbolic) == fill[1]);
  assert(amd_calls == 2 && seed_count == 1 && reuse_count == 1);
  assert(kls_factor(solver, a) == KLS_OK);
  for (UF_long c = 0; c < n; ++c) for (UF_long p = ptr[c]; p < ptr[c + 1]; ++p)
    b[idx[p]] += a[p] * (.25 * (c % 11) - 1.25);
  assert(kls_solve(solver, 1, b, 0, x, 0) == KLS_OK);
  for (UF_long c = 0; c < n; ++c) assert(fabs(x[c] - (.25 * (c % 11) - 1.25)) < 1e-10);
  kls_destroy(solver); free(ptr); free(idx); free(a); free(b); free(x);
}


static void one_component_checks(void) {
  enum { N = 257, NZ = 3 * N };
  UF_long ptr[N + 1], idx[NZ]; double a[NZ], b[N], x[N];
  for (UF_long c = 0; c < N; ++c) {
    ptr[c] = 3 * c;
    UF_long rows[3] = {(c + 1) % N, (c + 2) % N, (c + N - 1) % N};
    double values[3] = {4, -0.125, -0.125};
    for (int i = 1; i < 3; ++i) {
      UF_long row = rows[i]; double value = values[i]; int j = i;
      while (j && rows[j - 1] > row) { rows[j] = rows[j - 1]; values[j] = values[j - 1]; --j; }
      rows[j] = row; values[j] = value;
    }
    for (int i = 0; i < 3; ++i) { idx[3 * c + i] = rows[i]; a[3 * c + i] = values[i]; }
  }
  ptr[N] = NZ;
  for (int explicit_order = 0; explicit_order < 2; ++explicit_order) {
    kls_solver *solver = NULL; assert(kls_create(&solver) == KLS_OK);
    kls_options o; kls_default_options(&o); o.orientation = KLS_ORIENTATION_NORMAL;
    if (explicit_order) o.ordering = KLS_ORDERING_AMD;
    assert(kls_analyze_csc(solver, KLS_INDEX_INT64, N, ptr, idx, 0, &o) == KLS_OK);
    assert(solver->symbolic->nblocks == 1);
    assert(solver->symbolic->do_btf);
    assert(kls_factor(solver, a) == KLS_OK);
    for (int epoch = 0; epoch < 3; ++epoch) {
      memset(b, 0, sizeof b);
      for (UF_long c = 0; c < N; ++c) for (UF_long p = ptr[c]; p < ptr[c + 1]; ++p)
        b[idx[p]] += a[p] * (.25 * ((c + epoch) % 11) - 1.25);
      assert(kls_refactor_solve(solver, a, 1, b, 0, x, 0) == KLS_OK);
      for (UF_long c = 0; c < N; ++c) assert(fabs(x[c] - (.25 * ((c + epoch) % 11) - 1.25)) < 1e-10);
    }
    kls_destroy(solver);
  }
  /* Symmetric odd cycle with a missing structural diagonal: symmetry alone
     must not discard a useful BTF matching permutation. */
  for (UF_long c = 0; c < N; ++c) {
    ptr[c] = 2 * c;
    UF_long left = (c + N - 1) % N, right = (c + 1) % N;
    idx[2 * c] = left < right ? left : right;
    idx[2 * c + 1] = left < right ? right : left;
    a[2 * c] = a[2 * c + 1] = 1;
  }
  ptr[N] = 2 * N;
  for (int explicit_order = 0; explicit_order < 2; ++explicit_order) {
    kls_solver *solver = NULL; assert(kls_create(&solver) == KLS_OK);
    kls_options o; kls_default_options(&o); o.orientation = KLS_ORIENTATION_NORMAL;
    if (explicit_order) o.ordering = KLS_ORDERING_AMD;
    assert(kls_analyze_csc(solver, KLS_INDEX_INT64, N, ptr, idx, 0, &o) == KLS_OK);
    assert(solver->symbolic->nblocks == 1);
    assert(solver->symbolic->do_btf);
    assert(kls_factor(solver, a) == KLS_OK);
    memset(b, 0, sizeof b);
    for (UF_long c = 0; c < N; ++c) for (UF_long p = ptr[c]; p < ptr[c + 1]; ++p)
      b[idx[p]] += a[p] * (.25 * (c % 11) - 1.25);
    assert(kls_solve(solver, 1, b, 0, x, 0) == KLS_OK);
    for (UF_long c = 0; c < N; ++c) assert(fabs(x[c] - (.25 * (c % 11) - 1.25)) < 1e-10);
    kls_destroy(solver);
  }
  kls_solver *solver = NULL; assert(kls_create(&solver) == KLS_OK);
  kls_options o; kls_default_options(&o); o.ordering = KLS_ORDERING_NATURAL; o.use_btf = 0;
  o.orientation = KLS_ORIENTATION_NORMAL;
  assert(kls_analyze_csc(solver, KLS_INDEX_INT64, N, ptr, idx, 0, &o) == KLS_OK);
  assert(solver->symbolic->est_flops == TRILINOS_KLU_EMPTY);
  kls_stats st = {.struct_size = sizeof st}; assert(kls_get_stats(solver, &st) == KLS_OK);
  assert(st.estimated_flops == TRILINOS_KLU_EMPTY);
  kls_destroy(solver);
}

static void model_only_profile_checks(void) {
  char path[]="/tmp/kls-policy-ordering-XXXXXX";
  int fd=mkstemp(path);assert(fd>=0);close(fd);unlink(path);
  kls_policy_context p;kls_policy_init(&p,NULL);p.cost.enabled=1;
  p.cost.cache_bytes=kls_tuning_policy_cache_bytes();
  if(!p.cost.cache_bytes)return;
  assert(kls_tuning_write_policy_profile(path,1,&p.cost,123,0)==KLS_OK);
  kls_options o;kls_default_options(&o);o.orientation=KLS_ORIENTATION_NORMAL;
  const UF_long ptr[]={0,2,5,7},idx[]={0,1,0,1,2,1,2};
  for(int legacy=0;legacy<2;++legacy){
    kls_solver *s=NULL;assert(kls_create(&s)==KLS_OK);
    if(legacy){kls_tuning_values v;kls_tuning_defaults(&v);
      assert(kls_tuning_write_host_profile(path,1,&v,123,1,NULL)==KLS_OK);}
    o.tuning_profile_path=path;
    assert(kls_analyze_csc(s,KLS_INDEX_INT64,3,ptr,idx,0,&o)==KLS_OK);
    assert(s->stats.tuning_profile_active);
    assert(kls_kernel_tuning_overridden(s)==legacy);
    kls_destroy(s);
  }
  unlink(path);
}

int main(void) {
  enum{N=129,NZ=3*N-2};UF_long ptr[N+1],idx[NZ];double a[NZ],b[N],x[N];
  int z=0;for(int c=0;c<N;++c){ptr[c]=z;if(c){idx[z]=c-1;a[z++]=-1;}idx[z]=c;a[z++]=4;if(c+1<N){idx[z]=c+1;a[z++]=-1;}}ptr[N]=z;
  for(int format=0;format<2;++format)for(int threads=1;threads<=2;++threads){
    kls_solver *s=NULL;assert(kls_create(&s)==KLS_OK);
    kls_options o;kls_default_options(&o);o.threads=threads;o.orientation=KLS_ORIENTATION_NORMAL;
    seed_count=reuse_count=block_probes=amd_calls=0;
    int status=format?kls_analyze_csr(s,KLS_INDEX_INT64,N,ptr,idx,0,&o):kls_analyze_csc(s,KLS_INDEX_INT64,N,ptr,idx,0,&o);
    assert(status==KLS_OK);
    assert(seed_count>=1 && reuse_count>=1 && block_probes==0 && amd_calls==1);
    assert(s->stats.selected_ordering==KLS_ORDERING_AMD);
    assert(kls_factor(s,a)==KLS_OK);
    assert(block_probes == 0);
    for(int epoch=0;epoch<3;++epoch){
      memset(b,0,sizeof b);for(int c=0;c<N;++c)for(UF_long p=ptr[c];p<ptr[c+1];++p)b[idx[p]]+=a[p]*(.25*((c+epoch)%9)-1);
      assert(kls_refactor_solve(s,a,1,b,0,x,0)==KLS_OK);
      for(int c=0;c<N;++c)assert(fabs(x[c]-(.25*((c+epoch)%9)-1))<1e-10);
    }
    kls_destroy(s);
  }
  boundary_and_relabel_checks();
  concrete_block_checks();
  matched_representation_checks();
  one_component_checks();
  model_only_profile_checks();
  puts("policy profiles, ordering and public lifecycle checks passed");return 0;
}
