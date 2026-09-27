/* CPU-profile-priced partition proposals. All values remain in the existing
   native factor; this search changes only an unpublished start array.
   External consumers are priced by changed incoming edges, not full rescans. */
typedef struct {
  int n, stamp, expired;
  double deadline;
  const kls_snb_cost_model *cost;
  double *lu;
  const UF_long *lip,*llen,*uip,*ulen;
  int *width,*below,*map,*mark,*pmark,*counts,*producers;
} kls_snb_group_model;
typedef struct { double time, bytes; } kls_snb_group_price;

static int kls_snb_model_stamp(kls_snb_group_model *m) {
  if(m->stamp==INT_MAX){memset(m->mark,0,(size_t)m->n*sizeof(int));memset(m->pmark,0,(size_t)m->n*sizeof(int));m->stamp=0;}
  return ++m->stamp;
}

static double kls_snb_model_edge(const kls_snb_cost_model *cost,int w,int sc,int ph) {
  const double tr=(double)w*sc*(sc-1)/2;
  if(sc<8)return cost->narrow[0]+cost->narrow[1]*tr+
    cost->narrow[2]*(double)w*sc*ph;
  return cost->dense[0]+cost->dense[1]*tr+
    cost->dense[2]*(double)w*(ph/8)*sc+
    cost->dense[3]*(double)(w/4+w%4)*(ph/8)*sc+
    cost->dense[4]*(double)w*(ph%8)*sc+
    cost->dense[5]*(double)w*ph;
}
static int kls_snb_model_live(kls_snb_group_model *m) {
  if(m->expired || (m->deadline>0.0 && kls_now_seconds()>=m->deadline)) {
    m->expired=1;return 0;
  }
  return 1;
}
static int kls_snb_model_below(kls_snb_group_model *m,int c,int w) {
  int count=0,stamp=kls_snb_model_stamp(m);
  for(int j=c;j<c+w;++j) {
    if(!kls_snb_model_live(m)) break;
    UF_long *ri,len;double *rv;
    kls_klu_get_pointer(m->lu,m->lip,m->llen,j,&ri,&rv,&len);
    for(UF_long p=0;p<len;++p)if(ri[p]>=c+w&&m->mark[ri[p]]!=stamp) {
      m->mark[ri[p]]=stamp;++count;
    }
  }
  return count;
}
/* All split candidates share their column ranges. Grow the left union while
   dropping rows that enter its diagonal block; grow the right union backward
   with its fixed lower boundary. Each L column is scanned at most twice. */
static int kls_snb_model_split_below(kls_snb_group_model *m,int c,int w,
    int *left,int *right) {
  int count=0,stamp=kls_snb_model_stamp(m);
  for(int j=c;j<c+w-1;++j) {
    if(!kls_snb_model_live(m))return 0;
    if(m->mark[j]==stamp)--count;
    UF_long *ri,len;double *rv;
    kls_klu_get_pointer(m->lu,m->lip,m->llen,j,&ri,&rv,&len);
    for(UF_long p=0;p<len;++p)if(ri[p]>j&&m->mark[ri[p]]!=stamp) {
      m->mark[ri[p]]=stamp;++count;
    }
    left[j-c+1]=count;
  }
  count=0;stamp=kls_snb_model_stamp(m);
  for(int j=c+w-1;j>c;--j) {
    if(!kls_snb_model_live(m))return 0;
    UF_long *ri,len;double *rv;
    kls_klu_get_pointer(m->lu,m->lip,m->llen,j,&ri,&rv,&len);
    for(UF_long p=0;p<len;++p)if(ri[p]>=c+w&&m->mark[ri[p]]!=stamp) {
      m->mark[ri[p]]=stamp;++count;
    }
    right[j-c]=count;
  }
  return 1;
}
static kls_snb_group_price kls_snb_model_price(kls_snb_group_model *m,int c) {
  int w=m->width[c],ph=m->below[c],upper=0,edges=0,stamp=kls_snb_model_stamp(m);
  kls_snb_group_price price={0,0};if(!w)return price;
  for(int j=c;j<c+w;++j) {
    if(!kls_snb_model_live(m)) break;
    UF_long *ri,len;double *rv;
    kls_klu_get_pointer(m->lu,m->uip,m->ulen,j,&ri,&rv,&len);
    for(UF_long p=0;p<len;++p)if(ri[p]<c&&m->mark[ri[p]]!=stamp) {
      int producer=m->map[ri[p]];
      m->mark[ri[p]]=stamp;++upper;
      if(m->pmark[producer]!=stamp) {
        m->pmark[producer]=stamp;m->counts[producer]=0;m->producers[edges++]=producer;
      }
      ++m->counts[producer];
    }
  }
  double updates=0;
  for(int i=0;i<w;++i)updates+=(double)(w-i-1)*(w+ph-i-1);
  double panel=(double)w*(upper+w+ph);
  price.time=m->cost->panel[0]+m->cost->panel[1]*((double)w*ph+(double)w*(w-1)/2)+
    m->cost->panel[2]*updates+m->cost->panel[3]*panel;
  price.bytes=panel*sizeof(double)+sizeof(kls_snb_snode)+sizeof(int32_t);
  for(int e=0;e<edges;++e) {
    int producer=m->producers[e],sc=m->counts[producer],below=m->below[producer];
    price.time+=kls_snb_model_edge(m->cost,w,sc,below);
    price.bytes+=sizeof(kls_snb_edge)+(double)(sc+below)*sizeof(int32_t);
  }
  return price;
}
/* Below-row unions depend only on the immutable column range. Restoring a
   trial or committing its winner need not scan the factor again. */
static void kls_snb_model_restore(kls_snb_group_model *m,int c,int width,int below) {
  m->width[c]=width;
  if(width) {
    m->below[c]=below;
    for(int i=c;i<c+width;++i)m->map[i]=c;
  }
}
static void kls_snb_model_set(kls_snb_group_model *m,int c,int width) {
  kls_snb_model_restore(m,c,width,width?kls_snb_model_below(m,c,width):0);
}
/* External consumers keep their panels and all other incoming edges.
   Price only the edges from the changing range, using unions built once
   from the reverse U graph. No consumer U scan is needed per proposal. */
static kls_snb_group_price kls_snb_model_local(kls_snb_group_model *m,int c,int end,
    const int *affected,int count,const uint64_t *mask) {
  kls_snb_group_price total={0,0};
  for(int g=c;g<end;g+=m->width[g]) {
    kls_snb_group_price p=kls_snb_model_price(m,g);
    total.time+=p.time;total.bytes+=p.bytes;
    uint64_t bits=((UINT64_C(1)<<m->width[g])-1)<<(g-c);
    for(int i=0;i<count;++i) {
      int consumer=affected[i];
      int sc=__builtin_popcountll(mask[consumer]&bits);
      if(sc) {
        total.time+=kls_snb_model_edge(m->cost,m->width[consumer],sc,m->below[g]);
        total.bytes+=sizeof(kls_snb_edge)+(double)(sc+m->below[g])*sizeof(int32_t);
      }
    }
  }
  return total;
}
static int64_t kls_snb_model_regroup(double *lu,const UF_long *lip,const UF_long *llen,
    const UF_long *uip,const UF_long *ulen,UF_long nk,int64_t *starts,int64_t ns,const kls_snb_cost_model *cost,
    double deadline,double *base_out,double *candidate_out) {
  *base_out=*candidate_out=0.0;
  /* Capacity checks only: the incumbent's existing admission policy and
     the caller's workspace/time budgets already bound practical searches. */
  if(nk<2 || nk>INT_MAX || ns<1 ||
      (uint64_t)nk>(SIZE_MAX-sizeof(int))/(12*sizeof(int)+sizeof(uint64_t)))return ns;
  const int n=(int)nk;
  kls_snb_group_model m={0};m.n=n;m.cost=cost;m.deadline=deadline;m.lu=lu;m.lip=lip;m.llen=llen;m.uip=uip;m.ulen=ulen;
  int *workspace=calloc((size_t)n*11+(size_t)n+1,sizeof(int));
  if(!workspace)return ns;
  m.width=workspace;m.below=m.width+n;m.map=m.below+n;m.mark=m.map+n;
  m.pmark=m.mark+n;m.counts=m.pmark+n;m.producers=m.counts+n;
  int *affected=m.producers+n,*seen=affected+n,*offset=seen+n,*cursor=offset+n+1;
  int *consumers=NULL;
  uint64_t *mask=NULL;
  for(int64_t s=0;s<ns;++s){int c=starts[s],end=s+1<ns?starts[s+1]:n;kls_snb_model_set(&m,c,end-c);}
  for(int c=0;c<n;++c){
    if(!kls_snb_model_live(&m))goto cleanup;
    UF_long *ri,len;double *rv;kls_klu_get_pointer(lu,uip,ulen,c,&ri,&rv,&len);
    for(UF_long p=0;p<len;++p)++offset[ri[p]+1];}
  for(int c=0;c<n;++c){if(offset[c+1]>INT_MAX-offset[c])goto cleanup;offset[c+1]+=offset[c];}
  consumers=malloc((size_t)(offset[n]?offset[n]:1)*sizeof(int));
  if(!consumers){free(workspace);return ns;}
  mask=calloc((size_t)n,sizeof(uint64_t));
  if(!mask){free(consumers);free(workspace);return ns;}
  memcpy(cursor,offset,(size_t)n*sizeof(int));
  for(int c=0;c<n;++c){
    if(!kls_snb_model_live(&m))goto cleanup;
    UF_long *ri,len;double *rv;kls_klu_get_pointer(lu,uip,ulen,c,&ri,&rv,&len);
    for(UF_long p=0;p<len;++p)consumers[cursor[ri[p]]++]=c;}
  double base_time=0,base_bytes=0;
  for(int c=0;c<n;c+=m.width[c]){kls_snb_group_price p=kls_snb_model_price(&m,c);base_time+=p.time;base_bytes+=p.bytes;}
  if(m.expired)goto cleanup;
  double current_time=base_time,current_bytes=base_bytes;int changed=0,token=0;
  /* Refine from the opposite end only if the original search budget remains. */
  for(int pass=0;pass<2;++pass)for(int c=pass?m.map[n-1]:0;c>=0&&c<n;) {
    if(!kls_snb_model_live(&m)) break;
    int w=m.width[c],next=c+w,nw=next<n?m.width[next]:0;
    const int old_below=m.below[c],next_below=nw?m.below[next]:0;
    int merge=nw&&w+nw<=32,extent=merge?w+nw:w;
    int count=0; ++token;
    /* Union producer columns for external consumers once per search site. */
    for(int r=c;r<c+extent;++r) {
      for(int p=offset[r];p<offset[r+1];++p){int g=m.map[consumers[p]];
        if(g<c+extent)continue;
        if(seen[g]!=token){seen[g]=token;affected[count++]=g;mask[g]=0;}
        mask[g]|=UINT64_C(1)<<(r-c);}
    }
    kls_snb_group_price before=kls_snb_model_local(&m,c,c+extent,affected,count,mask),best=before;
    int left[33],right[33]; /* Existing maximum panel width is 32. */
    if(!kls_snb_model_split_below(&m,c,w,left,right))break;
    int best_split=0,best_left=0,best_right=0;
    for(int split=1;split<w;++split) {
      kls_snb_model_restore(&m,c,split,left[split]);kls_snb_model_restore(&m,c+split,w-split,right[split]);
      kls_snb_group_price after=kls_snb_model_local(&m,c,c+extent,affected,count,mask);
      if(!m.expired && after.time<best.time&&current_bytes+after.bytes-before.bytes<=base_bytes){best=after;best_split=split;best_left=m.below[c];best_right=m.below[c+split];}
      m.width[c+split]=0;kls_snb_model_restore(&m,c,w,old_below);
      if(m.expired)break;
    }
    if(merge && !m.expired) {
      m.width[next]=0;kls_snb_model_set(&m,c,w+nw);
      kls_snb_group_price after=kls_snb_model_local(&m,c,c+extent,affected,count,mask);
      if(!m.expired && after.time<best.time&&current_bytes+after.bytes-before.bytes<=base_bytes){best=after;best_split=-1;best_left=m.below[c];}
      kls_snb_model_restore(&m,c,w,old_below);kls_snb_model_restore(&m,next,nw,next_below);
    }
    if(m.expired)break;
    if(best_split) {
      ++changed;current_time+=best.time-before.time;current_bytes+=best.bytes-before.bytes;
      if(best_split<0){m.width[next]=0;kls_snb_model_restore(&m,c,w+nw,best_left);}
      else{kls_snb_model_restore(&m,c,best_split,best_left);kls_snb_model_restore(&m,c+best_split,w-best_split,best_right);}
    }
    if(pass)c=c>0?m.map[c-1]:-1;
    else c+=best_split<0?w+nw:w;
  }
  /* Interrupted proposals were restored above. Every committed delta includes
     outgoing consumers, so a completed prefix needs no unbudgeted full scan.
     Independent builder/pricing tests check this additive accounting. */
  if(changed && isfinite(base_time) && isfinite(current_time) &&
     current_time<base_time && current_bytes<=base_bytes) {
    *base_out=base_time;*candidate_out=current_time;
    ns=0;
    for(int c=0;c<n;c+=m.width[c]){starts[ns++]=c;}
  }
cleanup:
  free(mask);free(consumers);free(workspace);return ns;
}
