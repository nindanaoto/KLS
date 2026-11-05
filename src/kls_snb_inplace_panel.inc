/* Convert one unpublished width-16 column panel to the two eight-column
 * row streams using only caller-owned scratch bits. No floating operations.
 * The caller must finish all uses of its GEMM scratch before this call. */
#ifndef KLS_SNB_INPLACE_PANEL_INC
#define KLS_SNB_INPLACE_PANEL_INC
_Static_assert(sizeof(double)==sizeof(uint64_t),"64-bit panel words required");
static int kls_snb_column_panel_to_rows_inplace(double *panel,int32_t height,
    void *scratch,size_t scratch_bytes) {
  if(!panel || !scratch || height<=0 ||
      (size_t)height>SIZE_MAX/(16*sizeof(double)))return 0;
  const size_t count=(size_t)height*16,bytes=count*sizeof(double),bits=count/8;
  if(scratch_bytes<bits)return 0;
  const uintptr_t first=(uintptr_t)panel,work=(uintptr_t)scratch;
  if(first>UINTPTR_MAX-bytes || work>UINTPTR_MAX-bits ||
      (first<work+bits && work<first+bytes))return 0;
  unsigned char *visited=scratch;memset(visited,0,bits);
  for(size_t start=0;start<count;++start) {
    if(visited[start/8]&(unsigned char)(1u<<(start%8)))continue;
    uint64_t value;memcpy(&value,panel+start,sizeof(value));
    size_t at=start;
    do {
      const size_t column=at/(size_t)height,row=at-column*(size_t)height;
      const size_t next=(column/8)*8*(size_t)height+row*8+column%8;
      uint64_t displaced;memcpy(&displaced,panel+next,sizeof(displaced));
      memcpy(panel+next,&value,sizeof(value));
      visited[at/8]|=(unsigned char)(1u<<(at%8));
      value=displaced;at=next;
    } while(at!=start);
  }
  return 1;
}
#endif
