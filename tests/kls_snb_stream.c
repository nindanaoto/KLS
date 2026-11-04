static int fail_worker=-1;
#define KLS_SNB_STREAM_TEST_ENV_FAILURE(tid) ((tid)==fail_worker)
#include "../src/kls.c"
#include "../src/kls_snb_direct_product.inc"
#include "../src/kls_snb_slice.inc"
#undef NDEBUG
#include <assert.h>
#include "kls_snb_row_normal.inc"
static size_t capture(kls_solver *s,double *out) {
  size_t at=0;trilinos_klu_l_numeric *n=s->numeric;
  for(UF_long k=0;k<s->n;++k) {
    out[at++]=((double *)n->Udiag)[k];
    for(int upper=0;upper<2;++upper) {
      UF_long *ri,len;double *v;
      kls_klu_get_pointer(n->LUbx[0],upper?n->Uip:n->Lip,
        upper?n->Ulen:n->Llen,k,&ri,&v,&len);
      for(UF_long p=0;p<len;++p)out[at++]=v[p];
    }
  }
  return at;
}
static void private_closure(void) {
 kls_snb_state *st=calloc(1,sizeof(*st));assert(st);
 st->nblocks=1;st->nworkers=4;st->slevel_count=1;
 st->blocks=calloc(1,sizeof(*st->blocks));
 st->slevel_ptr=calloc(2,sizeof(*st->slevel_ptr));
 st->slevel_snode=calloc(8,sizeof(*st->slevel_snode));
 assert(st->blocks && st->slevel_ptr && st->slevel_snode);
 kls_snb_block *b=st->blocks;b->snode_count=8;
 b->snodes=calloc(8,sizeof(*b->snodes));b->edges=calloc(5,sizeof(*b->edges));
 assert(b->snodes && b->edges);st->slevel_ptr[1]=8;
 for(int j=0;j<8;++j){b->snodes[j].width=1;b->snodes[j].height=1;st->slevel_snode[j]=j;}
 const int pred[5]={0,1,4,5,6};for(int e=0;e<5;++e)b->edges[e].producer=pred[e];
 b->snodes[4].edges_off=0;b->snodes[4].edge_count=1;
 b->snodes[5].edges_off=1;b->snodes[5].edge_count=1;
 b->snodes[6].edges_off=2;b->snodes[6].edge_count=2;
 b->snodes[7].edges_off=4;b->snodes[7].edge_count=1;
 assert(kls_snb_build_private(st,0,SIZE_MAX));
 for(int j=0;j<4;++j)assert(st->private_owner[j]==j);
 assert(st->private_owner[4]==0 && st->private_owner[5]==1);
 assert(st->private_owner[6]==-1 && st->private_owner[7]==-1);
 assert(st->private_ptr[4]==6 && st->private_pipeline_count==2);
 assert(st->private_pipeline[0]==6 && st->private_pipeline[1]==7);
 kls_snb_destroy(st);
}
static void direct_product_oracle(void) {
#if defined(__AVX512F__)
 double l[128*64],b[64*12],t[40*12],expected[128*12],actual[128*12];
 int32_t sub[32],dst[40];
 for(int i=0;i<128*64;++i)l[i]=.03125*((i*7)%19-9);
 for(int i=0;i<64*12;++i)b[i]=.0625*((i*11)%17-8);
 const int counts[]={1,2,7,16,17,32};
 for(int scase=0;scase<6;++scase)for(int w=1;w<=9;++w)
 for(int rows=1;rows<=33;++rows)for(int gap=0;gap<2;++gap) {
  const int sc=counts[scase];
  for(int i=0;i<sc;++i)sub[i]=gap?2*i:i;
  for(int i=0;i<rows;++i)dst[i]=gap?(i%5==0?-1:3*i):i;
  for(int i=0;i<128*12;++i)expected[i]=actual[i]=.125*(i%7-3);
  kls_snb_dense_product(l,128,sub,sc,b,64,w,0,rows,t);
  for(int c=0;c<w;++c)for(int r=0;r<rows;++r)
   if(dst[r]>=0)expected[c*128+dst[r]]-=t[c*rows+r];
  assert(kls_snb_direct_product(l,128,sub,sc,b,64,w,rows,dst,actual,128));
  assert(!memcmp(actual,expected,sizeof actual));
 }
#endif
}
int main(void) {
  kls_snb_row_normal_checks();
  direct_product_oracle();
  private_closure();
  enum {N=129};int64_t ap[N+1],ai[N*N];double a[N*N],b[N],x[N];
  for(int c=0;c<N;++c){ap[c]=c*N;for(int r=0;r<N;++r)ai[c*N+r]=r;}ap[N]=N*N;
  double oracle[2*N*N+N],actual[2*N*N+N];
  kls_options o;kls_default_options(&o);o.threads=4;o.use_btf=0;o.scale=0;
  o.ordering=KLS_ORDERING_NATURAL;o.orientation=KLS_ORIENTATION_NORMAL;
  kls_solver *s=NULL;assert(kls_create(&s)==KLS_OK);
  assert(kls_analyze_csc(s,KLS_INDEX_INT64,N,ap,ai,0,&o)==KLS_OK);
  for(int c=0;c<N;++c)for(int r=0;r<N;++r)a[c*N+r]=r==c?N+1.:.001;
  assert(kls_factor(s,a)==KLS_OK);
  assert(kls_parallel_lu_sort(s));s->snode_prepared=1;s->snb_declined=0;
  assert(kls_snb_prepare(s));assert(s->snb->nworkers>=4);
  assert(kls_snb_build_private(s->snb,0,SIZE_MAX));
  assert(s->snb->private_threads==4);
  assert(ensure_egraph_refactor_pool(s,4));
  for(int epoch=0;epoch<12;++epoch) {
    s->snb->private_threads=(epoch&1)?0:4;
    memset(b,0,sizeof b);
    for(int c=0;c<N;++c)for(int r=0;r<N;++r) {
      a[c*N+r]=r==c?N+1.+epoch:.001*(1+(r+c+epoch)%7);
      b[r]+=a[c*N+r]*(1.+.125*((c+epoch)%7));
    }
    assert(kls_refactor(s,a)==KLS_OK);
    double *v=NULL;assert(prepare_numeric_values(s,a,&v)==KLS_OK);
    kls_egraph_refactor_shared *shared=&s->egraph_pool->shared;
    assert(!pthread_mutex_lock(&shared->lock));
    assert(kls_snb_stream_refactor(s,v,4)==-1);
    assert(!pthread_mutex_unlock(&shared->lock));
    assert(!fesetround(FE_UPWARD));assert(kls_snb_stream_refactor(s,v,4)==-1);
    assert(!fesetround(FE_TONEAREST));
    fail_worker=1;assert(kls_snb_stream_refactor(s,v,4)==0);
    fail_worker=-1;
    assert(!s->componentwise_work && !s->componentwise_work_context);
    assert(kls_snb_owned_masked_environment==0);
    assert(kls_snb_refactor(s,v)==1);
    const size_t count=capture(s,oracle);
    assert(kls_snb_refactor_prelude(s,v));
    for(int64_t j=0;j<s->snb->blocks[0].snode_count;++j) {
      const int width=s->snb->blocks[0].snodes[j].width;
      for(int begin=0;begin<width;begin+=4)
        kls_snb_process_target_slice(s,s->snb,&s->snb->blocks[0],j,v,NULL,-1,
          s->snb->gemm_w[0],0,begin,begin+4<width?begin+4:width);
      kls_snb_finalize_target(s,s->snb,&s->snb->blocks[0],j);
    }
    assert(kls_snb_refactor_epilogue(s));
    assert(capture(s,actual)==count && !memcmp(actual,oracle,count*sizeof(double)));

    assert(kls_snb_stream_refactor(s,v,4)==1);
    assert(capture(s,actual)==count);
    assert(!memcmp(actual,oracle,count*sizeof(double)));
    assert(!s->componentwise_work && !s->componentwise_work_context);
    assert(kls_snb_owned_masked_environment==0);
    assert(kls_solve(s,1,b,0,x,0)==KLS_OK);
    for(int c=0;c<N;++c)assert(fabs(x[c]-(1.+.125*((c+epoch)%7)))<1e-11);
  }
  kls_destroy(s);puts("PASS owned SNB stream current values/RHS, busy pool, FP decline and joins");
  return 0;
}
