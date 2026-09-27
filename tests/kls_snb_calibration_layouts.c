#include "../bench/kls_snb_calibrate.c"

int main(void) {
  int sub[64];
  for(int width=1;width<=64;++width)for(int sc=1;sc<=width;++sc)
    for(int gap=0;gap<2;++gap) {
      const int possible=!gap||(sc>1&&sc<width);
      if(cost_subcols(width,sc,gap,sub)!=possible)return 1;
      if(!possible)continue;
      if(sub[0]<0||sub[sc-1]>=width)return 2;
      for(int i=1;i<sc;++i)if(sub[i]<=sub[i-1])return 3;
      if((sub[sc-1]-sub[0]!=sc-1)!=gap)return 4;
    }
  if(cost_subcols(32,0,0,sub)||cost_subcols(32,33,0,sub))return 5;
  return 0;
}
