/* Private queues contain only dependency-closed subgraphs. Any merge of
 * different owners, or any shared predecessor, keeps the node shared. */
static int kls_snb_build_private(kls_snb_state *st,double deadline,size_t limit) {
  const int threads=st->nworkers;
  if(st->nblocks!=1 || threads<2 || threads>64)return 1;
  const kls_snb_block *blk=&st->blocks[0];const int64_t n=blk->snode_count;
  if(n<=0 || n>INT32_MAX || !st->slevel_ptr ||
     st->slevel_ptr[st->slevel_count]!=n)return 1;
  st->private_owner=kls_snb_allocate(st,(size_t)n,sizeof(*st->private_owner),0,deadline,limit);
  st->private_nodes=kls_snb_allocate(st,(size_t)n,sizeof(*st->private_nodes),0,deadline,limit);
  st->private_pipeline=kls_snb_allocate(st,(size_t)n,sizeof(*st->private_pipeline),0,deadline,limit);
  if(!st->private_owner || !st->private_nodes || !st->private_pipeline)return 0;
  int64_t counts[64]={0};long double work[64]={0};
  for(int64_t j=0;j<n;++j) {
    if(deadline>0 && kls_now_seconds()>=deadline)return 0;
    const kls_snb_snode *node=&blk->snodes[j];int owner=-2;
    for(int e=0;e<node->edge_count;++e) {
      const int32_t pred=blk->edges[node->edges_off+e].producer;
      if(pred<0 || pred>=j)return 0;
      const int prior=st->private_owner[pred];
      if(prior<0 || (owner>=0 && owner!=prior)){owner=-1;break;}
      owner=prior;
    }
    if(owner==-2) {
      owner=0;for(int t=1;t<threads;++t)if(work[t]<work[owner])owner=t;
    }
    st->private_owner[j]=(int16_t)owner;
    if(owner>=0) {
      ++counts[owner];
      work[owner]+=(long double)node->width*node->width*node->height;
      for(int e=0;e<node->edge_count;++e) {
        const kls_snb_edge *edge=&blk->edges[node->edges_off+e];
        work[owner]+=2.L*node->width*edge->sub_count*blk->snodes[edge->producer].below;
      }
    }
  }
  for(int t=0;t<threads;++t)st->private_ptr[t+1]=st->private_ptr[t]+counts[t];
  int64_t cursor[64];memcpy(cursor,st->private_ptr,(size_t)threads*sizeof(*cursor));
  for(int64_t j=0;j<n;++j)if(st->private_owner[j]>=0)
    st->private_nodes[cursor[st->private_owner[j]]++]=(int32_t)j;
  for(int64_t at=0;at<n;++at) {
    const int32_t j=st->slevel_snode[at];
    if(j<0 || j>=n)return 0;
    if(st->private_owner[j]<0)st->private_pipeline[st->private_pipeline_count++]=j;
  }
  st->private_threads=threads;return 1;
}
