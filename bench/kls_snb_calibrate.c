/* Private executable-side calibration of the exact compiled kernels.
   Like the private kernel tests, this TU does not add a public solver API. */
#include "../src/kls.c"
#if defined(__GNUC__) || defined(__clang__)
#define COST_NOINLINE __attribute__((noinline))
#else
#define COST_NOINLINE
#endif

typedef struct kls_cost_sample {
  double x[6], y;
  int kind, validation;
} kls_cost_sample;

static double cost_median5(double *v) {
  for(int i=1;i<5;++i){double x=v[i];int j=i;while(j&&v[j-1]>x){v[j]=v[j-1];--j;}v[j]=x;}
  return v[2];
}

static void cost_fit(const kls_cost_sample *rows,int count,int kind,int nf,double *c) {
  double gram[6][6]={{0}},rhs[6]={0},scale[6],v[6]={0};
  for(int r=0;r<count;++r)if(rows[r].kind==kind&&!rows[r].validation) {
    /* Relative-error weighting prevents large panels dominating panel setup. */
    const double weight=kind==2?1.0/rows[r].y:1.0;
    for(int i=0;i<nf;++i) {
      double xi=rows[r].x[i]*weight;
      rhs[i]+=xi*rows[r].y*weight;
      for(int j=0;j<nf;++j)gram[i][j]+=xi*rows[r].x[j]*weight;
    }
  }
  for(int i=0;i<nf;++i)scale[i]=fmax(sqrt(gram[i][i]),1.0);
  for(int i=0;i<nf;++i){rhs[i]/=scale[i];for(int j=0;j<nf;++j)gram[i][j]/=scale[i]*scale[j];}
  for(int it=0;it<2000;++it)for(int i=0;i<nf;++i)if(gram[i][i]>0.0) {
    double residual=rhs[i];for(int j=0;j<nf;++j)residual-=gram[i][j]*v[j];
    v[i]=fmax(0.0,v[i]+residual/gram[i][i]);
  }
  for(int i=0;i<nf;++i)c[i]=v[i]/scale[i];
}

static COST_NOINLINE void cost_edge(kls_snb_block *b,kls_snb_snode *sn,
    kls_snb_edge *e,int dense,double *scratch) {
  kls_snb_update_edge(b,sn,e,dense?0:33,scratch);
}

/* Arithmetic reference only, never used for timing or fitted costs. */
static void cost_panel_reference(double *p,int w,int h) {
  for(int i=0;i<w;++i) {
    double *ci=p+(size_t)i*h;const double pivot=ci[i];
    for(int r=i+1;r<h;++r)ci[r]/=pivot;
    for(int j=i+1;j<w;++j) {
      double *cj=p+(size_t)j*h;const double u=cj[i];
      if(u==0.0)continue;
      for(int r=i+1;r<h;++r)cj[r]-=ci[r]*u;
    }
  }
}

static int cost_panel_sample(int w,int ph,double *seconds,double *observed) {
  const int h=w+ph,n=h>1?h:2,nz=w*h+n-w;
  int *ptr=malloc(((size_t)n+1)*sizeof(int)),*idx=malloc((size_t)nz*sizeof(int));
  double *a=malloc((size_t)nz*sizeof(double));
  double *reference=malloc((size_t)w*h*sizeof(double));
  kls_solver *s=NULL;kls_snb_state *st=NULL;
  kls_snb_partition partition={0};
  int ok=ptr&&idx&&a&&reference;
  if(!ok)goto done;
  int at=0;
  for(int c=0;c<n;++c) {
    ptr[c]=at;
    if(c<w)for(int r=0;r<h;++r) {
      idx[at]=r;a[at++]=r==c?8.0:((r+c*3)%17+1)*0.001;
    } else {idx[at]=c;a[at++]=8.0;}
  }
  ptr[n]=at;
  kls_options options;kls_default_options(&options);
  options.threads=1;options.backend=KLS_BACKEND_SERIAL;
  options.ordering=KLS_ORDERING_NATURAL;options.use_btf=0;
  options.orientation=KLS_ORIENTATION_NORMAL;
  if(kls_create(&s)||kls_analyze_csc(s,KLS_INDEX_INT32,n,ptr,idx,0,&options)||
     kls_factor(s,a)||s->symbolic->nblocks!=1){ok=0;goto done;}
  s->snode_prepared=1;
  partition.starts=malloc((size_t)n*sizeof(int64_t));
  partition.counts=malloc(sizeof(int64_t));
  if(!partition.starts||!partition.counts||!kls_parallel_lu_sort(s)){ok=0;goto done;}
  partition.counts[0]=1+n-w;partition.starts[0]=0;
  for(int i=1;i<partition.counts[0];++i)partition.starts[i]=w+i-1;
  st=kls_snb_build(s,&partition,0.0,SIZE_MAX);
  if(!st){ok=0;goto done;}
  kls_snb_block *block=&st->blocks[0];
  const kls_snb_snode *sn=&block->snodes[0];
  if(sn->width!=w||sn->height!=h||sn->uw||sn->edge_count){ok=0;goto done;}
  memcpy(reference,a,(size_t)w*h*sizeof(double));
  cost_panel_reference(reference,w,h);
  kls_snb_process_target(s,st,block,0,a,NULL,0,st->gemm_w[0],0);
  for(int i=0;i<w*h;++i)if(!isfinite(block->panels[sn->panel+i])||
      fabs(block->panels[sn->panel+i]-reference[i])>1e-12){ok=0;goto done;}
  double samples[5];
  for(int rep=0;rep<5;++rep) {
    const double start=kls_now_seconds();
    for(int it=0;it<64;++it)
      kls_snb_process_target(s,st,block,0,a,NULL,0,st->gemm_w[0],0);
    samples[rep]=(kls_now_seconds()-start)/64.0;
  }
  *seconds=cost_median5(samples);*observed=block->panels[sn->panel+w*h-1];
done:
  free(partition.starts);free(partition.counts);kls_snb_destroy(st);kls_destroy(s);
  free(ptr);free(idx);free(a);free(reference);return ok;
}

/* Layout is independent of subset width. A full producer and a singleton
   have no genuinely gapped alternative; do not duplicate those samples. */
static int cost_subcols(int producer_width,int sc,int gapped,int *sub) {
  if(sc<1||sc>producer_width||(gapped&&(sc==1||sc==producer_width)))return 0;
  for(int i=0;i<sc;++i)
    sub[i]=gapped?(int)((int64_t)i*(producer_width-1)/(sc-1)):i;
  return 1;
}

int kls_snb_calibrate_model(kls_snb_cost_model *model,uint64_t seed) {
  const int widths[2][11]={{1,2,3,4,7,8,9,16,24,31,32},{5,6,10,12,15,17,20,23,25,29,30}};
  const int heights[2][9]={{0,1,7,8,9,17,33,65,129},{2,5,6,10,15,31,63,127,255}};
  const int panel_heights[]={0,1,7,8,9,17,33,65,129,255};
  int capacity=32*10;
  for(int set=0;set<2;++set)for(int si=0;si<11;++si)
    capacity+=11*9*2*2*((widths[set][si]>1&&widths[set][si]<32)?2:1);
  kls_cost_sample *rows=calloc((size_t)capacity,sizeof(*rows));
  double *base=malloc(25000*sizeof(double)),*panels=malloc(25000*sizeof(double));
  double *answer=malloc(25000*sizeof(double)),*scratch=malloc(8192*sizeof(double));
  int count=0,ok=rows&&base&&panels&&answer&&scratch;
  volatile double checksum=0.0;
  memset(model,0,sizeof(*model));
  if(!ok)goto done;
  for(int set=0;set<2;++set)for(int wi=0;wi<11;++wi)
   for(int si=0;si<11;++si)for(int hi=0;hi<9;++hi)for(int holes=0;holes<2;++holes)
    for(int gapped=0;gapped<2;++gapped) {
    const int w=widths[set][wi],sc=widths[set][si],ph=heights[set][hi];
    const int h=sc+w+ph,start=32*(32+ph),length=start+w*h;
    int32_t sub[32],dst[255];kls_snb_snode sn[2]={{0}};
    if(!cost_subcols(32,sc,gapped,sub))continue;
    kls_snb_edge edge={0};kls_snb_block block={0};
    sn[0].width=32;sn[0].below=ph;sn[0].height=32+ph;
    sn[1].start=32;sn[1].width=w;sn[1].uw=sc;sn[1].below=ph;
    sn[1].height=h;sn[1].panel=start;sn[1].edge_count=1;
    edge.sub_count=sc;block.edges=&edge;block.snodes=sn;block.panels=panels;block.subcols=sub;block.edge_dsts=dst;
    /* Both row unions are sorted, as in prepared SNB schedules. The target
       includes its complete diagonal block; holes represent absent rows. */
    for(int r=0;r<ph;++r)dst[r]=holes&&r%3==0?-1:sc+w+r;
    for(int i=0;i<length;++i)base[i]=((i*(17+2*set)+(int)(seed%29))%29-14)/10000.0;
    memcpy(panels,base,(size_t)length*sizeof(double));cost_edge(&block,&sn[1],&edge,0,scratch);
    memcpy(answer,panels,(size_t)length*sizeof(double));
    memcpy(panels,base,(size_t)length*sizeof(double));cost_edge(&block,&sn[1],&edge,1,scratch);
    for(int i=0;i<length;++i)if(!isfinite(panels[i])||fabs(panels[i]-answer[i])>1e-12){ok=0;goto done;}
    double times[2][5];
    for(int rep=0;rep<5;++rep)for(int arm=0;arm<2;++arm) {
      const int dense=(arm+rep)%2;
      memcpy(panels,base,(size_t)length*sizeof(double));cost_edge(&block,&sn[1],&edge,dense,scratch);
      const double t=kls_now_seconds();
      for(int it=0;it<64;++it)cost_edge(&block,&sn[1],&edge,dense,scratch);
      times[dense][rep]=(kls_now_seconds()-t)/64.0;checksum+=panels[length-1];
    }
    const double tr=(double)w*sc*(sc-1)/2;
    rows[count++]=(kls_cost_sample){{1,tr,(double)w*sc*ph,0,0,0},cost_median5(times[0]),0,set};
    rows[count++]=(kls_cost_sample){{1,tr,(double)w*(ph/8)*sc,
      (double)(w/4+w%4)*(ph/8)*sc,(double)w*(ph%8)*sc,(double)w*ph},cost_median5(times[1]),1,set};
  }
  for(int w=1;w<=32;++w)for(int hi=0;hi<10;++hi) {
    const int ph=panel_heights[hi],h=w+ph;
    double seconds=0.0,observed=0.0,updates=0.0;
    for(int i=0;i<w;++i)updates+=(double)(w-i-1)*(w+ph-i-1);
    if(!cost_panel_sample(w,ph,&seconds,&observed)){ok=0;goto done;}
    checksum+=observed;
    /* As for edges, train the endpoints and validate intervening widths;
       extrapolation outside the calibrated kernel envelope is not required. */
    int validation=1;
    for(int i=0;i<11;++i)if(w==widths[0][i])validation=0;
    rows[count++]=(kls_cost_sample){{1,(double)w*ph+(double)w*(w-1)/2,updates,(double)w*h,0,0},seconds,2,validation};
  }
  if(count!=capacity){ok=0;goto done;}
  for(int i=0;i<count;++i)if(!isfinite(rows[i].y)||rows[i].y<=0.0){ok=0;goto done;}
  cost_fit(rows,count,0,3,model->narrow);cost_fit(rows,count,1,6,model->dense);
  cost_fit(rows,count,2,4,model->panel);
  /* Fit training only. Reject a model whose held-out shapes exceed the
     worst training relative error; no circuit data informs coefficients. */
  for(int kind=0;kind<3;++kind) {
    const int nf=kind==0?3:kind==1?6:4;
    const double *c=kind==0?model->narrow:kind==1?model->dense:model->panel;
    double errors[2]={0,0};
    for(int r=0;r<count;++r)if(rows[r].kind==kind) {
      double predicted=0.0;for(int j=0;j<nf;++j)predicted+=c[j]*rows[r].x[j];
      /* NNLS may assign zero intercept to a genuinely empty edge (one
         producer column, no below rows). Zero arithmetic may have zero
         modeled cost; every nonempty kernel still needs a positive price.
         Its measured call overhead remains in the validation error below. */
      int empty_edge=kind!=2;
      for(int j=1;j<nf;++j)if(rows[r].x[j]!=0.0)empty_edge=0;
      if(!isfinite(predicted)||predicted<0.0||(predicted==0.0&&!empty_edge)) {
        fprintf(stderr,"SNB invalid prediction: kind=%d row=%d cost=%.9g\n",kind,r,predicted);
        ok=0;goto done;
      }
      double error=fabs(predicted-rows[r].y)/rows[r].y;
      if(error>errors[rows[r].validation])errors[rows[r].validation]=error;
    }
    fprintf(stderr,"SNB generated cost %d: train-max-error=%.6g validation-max-error=%.6g\n",kind,errors[0],errors[1]);
    if(errors[1]>errors[0])ok=0;
  }
done:
  fprintf(stderr,"SNB generated calibration: samples=%d valid=%d checksum=%.9g\n",count,ok,(double)checksum);
  free(rows);free(base);free(panels);free(answer);free(scratch);
  if(ok)model->enabled=1;else memset(model,0,sizeof(*model));
  return ok;
}
