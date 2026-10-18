// Diagnostic only. Never preload the recorder or run replay during competitive timing.
#include <subtree_lu.h>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

static_assert(sizeof(int)==4 && sizeof(double)==8, "native capture ABI");
using Header=std::array<uint64_t,8>;
static void require(bool p,const char *s) { if(!p) throw std::runtime_error(s); }
struct File {
 std::ifstream f;
 uint64_t size=0,pos=0;
 explicit File(const std::string& path): f(path,std::ios::binary) {
  require(bool(f),"cannot open input"); f.seekg(0,std::ios::end);
  auto end=f.tellg();require(end>=0,"cannot size input");size=uint64_t(end);seek(0);
 }
 void seek(uint64_t off) {
  require(off<=size && off<=uint64_t(std::numeric_limits<std::streamoff>::max()),"invalid offset");
  f.clear();f.seekg(static_cast<std::streamoff>(off));require(bool(f),"seek failed");pos=off;
 }
 void read(void *p,uint64_t bytes) {
  require(bytes<=size-pos && bytes<=uint64_t(std::numeric_limits<std::streamsize>::max()),"truncated input");
  if(bytes) { f.read(static_cast<char*>(p),static_cast<std::streamsize>(bytes));require(bool(f),"read failed"); }
  pos+=bytes;
 }
 Header header() { Header h{};read(h.data(),sizeof(h));return h; }
 template<class T> std::vector<T> array(uint64_t n) {
  require(n<=(size-pos)/sizeof(T) && n<=std::numeric_limits<size_t>::max()/sizeof(T),"array extent");
  std::vector<T> v(static_cast<size_t>(n));read(v.data(),n*sizeof(T));return v;
 }
};
struct Initial {
 int n,threads;uint64_t nnz,flags;
 std::vector<int> ap,ai,rp,cp;std::vector<double> ax;
};
static Initial initial(File& f) {
 auto h=f.header();
 require(h[0]==UINT64_C(0x53544c5543535231) && h[1]==2,"initial magic/version");
 require(h[2]>0 && h[2]<=INT32_MAX && h[3]<=INT32_MAX && h[4]>0 && h[4]<=INT32_MAX,"initial dimensions/threads");
 require(!(h[5]&~UINT64_C(31)) && (h[5]&16) && (h[3]==0 || (h[5]&8)) && h[6]==0 && h[7]==0,"initial flags");
 Initial a{int(h[2]),int(h[4]),h[3],h[5],{},{},{},{},{}};
 a.ap=f.array<int>(uint64_t(a.n)+1);a.ai=f.array<int>(a.nnz);
 if(a.flags&1)a.ax=f.array<double>(a.nnz);
 if(a.flags&2)a.rp=f.array<int>(a.n);
 if(a.flags&4)a.cp=f.array<int>(a.n);
 require(f.pos==f.size && a.ap.front()==0 && a.ap.back()==int(a.nnz),"CSR extent");
 for(int i=0;i<a.n;++i)require(a.ap[i]>=0 && a.ap[i]<=a.ap[i+1],"CSR pointers");
 for(int i:a.ai)require(i>=0 && i<a.n,"CSR index");
 return a;
}
static int status(uint64_t bits) {
 uint32_t low=static_cast<uint32_t>(bits);int s;std::memcpy(&s,&low,sizeof(s));
 require(bits==static_cast<uint64_t>(static_cast<int64_t>(s)),"noncanonical status bits");return s;
}
struct Record { Header h;uint64_t payload;int result; };
struct Call { Record enter,leave;bool returned=false; };
static std::vector<Call> preflight(File& f,const Initial& a) {
 auto h=f.header();
 require(h[0]==UINT64_C(0x53544c5553545231) && h[1]==1 && h[2]==uint64_t(a.n) && h[3]==a.nnz && h[4]==uint64_t(a.threads),"stream header mismatch");
 require(h[5]==((UINT64_C(4)<<32)|8) && h[6]>0 && h[7]==0,"stream ABI/generation");
 std::vector<Call> calls;bool open=false;uint64_t next=1;
 while(f.pos<f.size) {
  Record r{f.header(),0,0};auto& q=r.h;r.payload=f.pos;r.result=status(q[3]);
  require(q[0]>=1 && q[0]<=4 && (q[1]==1 || q[1]==2) && q[7]==0,"record stage/edge/reserved");
  require(q[4]<=2,"payload kind");
  const bool enter=q[1]==1;
  if(enter) {
   require(!open && q[2]==next && r.result==0,"enter order/status");
   require((next==1 && q[0]==1) || (next>1 && q[0]!=1),"analyze generation order");
  } else {
   require(open && !calls.empty(),"unmatched leave");
   const auto& e=calls.back().enter.h;
   require(q[0]==e[0] && q[2]==e[2] && q[6]==e[6],"leave identity/flags");
  }
  uint64_t kind=0,count=0;
  if(q[0]==1) {
   require(q[6]==(a.flags&7),"analyze pointer flags");
   if(!enter && r.result==0) {kind=2;count=uint64_t(a.n)*((q[6]&2?1:0)+(q[6]&4?1:0));}
  } else if(q[0]==2 || q[0]==3) {
   require(q[6]<=1,"numeric pointer flags");
   if(enter) {kind=1;count=q[6]?a.nnz:0;}
  } else {
   require(q[6]<=7 && (!(q[6]&4) || (q[6]&3)==3),"solve pointer flags");
   if(enter) {kind=1;count=(q[6]&1)?uint64_t(a.n):0;}
   else if(r.result==0) {kind=1;count=(q[6]&2)?uint64_t(a.n):0;}
  }
  require(q[4]==kind && q[5]==count,"payload type/count");
  uint64_t width=kind==2?4:8;
  require(count<=(f.size-f.pos)/width,"truncated payload");f.seek(f.pos+count*width);
  if(enter) {calls.push_back(Call{r,{},false});open=true;}
  else {calls.back().leave=r;calls.back().returned=true;open=false;++next;}
 }
 require(!calls.empty(),"empty stream");return calls;
}
template<class T> static std::vector<T> payload(File& f,const Record& r) {
 f.seek(r.payload);return f.array<T>(r.h[5]);
}
template<class T> static bool same(const std::vector<T>& a,const std::vector<T>& b) {
 return a.size()==b.size() && (a.empty() || std::memcmp(a.data(),b.data(),a.size()*sizeof(T))==0);
}
static void boundary(const Call& c,const char *edge,int s) {
 std::fprintf(stderr,"REPLAY sequence=%llu stage=%llu %s status=%d captured_return=%d\n",
  (unsigned long long)c.enter.h[2],(unsigned long long)c.enter.h[0],edge,s,c.returned?1:0);std::fflush(stderr);
}
int main(int argc,char **argv) {
 try {
  require(argc==3 || argc==4,"usage: replay GENERATION_DIR --validate-only|--execute-xyce [--include-unreturned]");
  bool validate=std::string(argv[2])=="--validate-only";
  require(validate || std::string(argv[2])=="--execute-xyce","explicit replay mode required");
  bool include=argc==4;
  require(!include || (!validate && std::string(argv[3])=="--include-unreturned"),"invalid unreturned-call option");
  File csr(std::string(argv[1])+"/initial.csr"),stream(std::string(argv[1])+"/calls.bin");
  auto a=initial(csr);auto calls=preflight(stream,a);
  std::fprintf(stderr,"REPLAY preflight n=%d nnz=%llu threads=%d calls=%zu trailing_enter=%d\n",a.n,
   (unsigned long long)a.nnz,a.threads,calls.size(),calls.back().returned?0:1);std::fflush(stderr);
  if(validate)return 0;
  // Explicit Xyce diagnostic preset; the format does not serialize parm[] or
  // constructor history. This is not a replay claim for arbitrary callers.
  subtree_lu::SubtreeLU<int,double> solver;
  solver.parm[subtree_lu::I_PIVTOL]=1000;solver.parm[subtree_lu::I_MEM_FACTOR]=150;
  int idummy=0;double ddummy=0;
  std::vector<double> values,rhs,out;
  for(const auto& c:calls) {
   if(!c.returned && !include) {
    std::fprintf(stderr,"REPLAY retained unreturned sequence=%llu; not invoked\n",(unsigned long long)c.enter.h[2]);return 4;
   }
   uint64_t stage=c.enter.h[0],flags=c.enter.h[6];int result=0;
   if(stage==2 || stage==3)values=payload<double>(stream,c.enter);
   if(stage==4) {rhs=payload<double>(stream,c.enter);out.assign(a.n,0);}
   boundary(c,"enter",0);
   if(stage==1)result=solver.analyze(a.n,a.ap.data(),(a.flags&8)?(a.ai.empty()?&idummy:a.ai.data()):nullptr,
    (flags&1)?(a.ax.empty()?&ddummy:a.ax.data()):nullptr,a.threads,
    (flags&2)?a.rp.data():nullptr,(flags&4)?a.cp.data():nullptr);
   else if(stage==2)result=solver.factorize((flags&1)?(values.empty()?&ddummy:values.data()):nullptr);
   else if(stage==3)result=solver.refactorize((flags&1)?(values.empty()?&ddummy:values.data()):nullptr);
   else result=solver.solve((flags&1)?rhs.data():nullptr,(flags&2)?((flags&4)?rhs.data():out.data()):nullptr);
   boundary(c,"leave",result);
   if(!c.returned) {
    std::fprintf(stderr,"REPLAY captured call had no return; fresh replay returned %d (not a captured status match)\n",result);return 5;
   }
   require(result==c.leave.result,"replayed status differs; replay stopped");
   if(stage==1 && result==0) {
    auto expected=payload<int>(stream,c.leave);std::vector<int> actual;
    if(flags&2)actual.insert(actual.end(),a.rp.begin(),a.rp.end());
    if(flags&4)actual.insert(actual.end(),a.cp.begin(),a.cp.end());
    require(same(actual,expected),"replayed output permutations differ; replay stopped");
   }
   if(stage==4 && result==0 && (flags&2)) {
    auto expected=payload<double>(stream,c.leave);
    require(same((flags&4)?rhs:out,expected),"replayed solution bits differ; replay stopped");
   }
  }
  std::fprintf(stderr,"REPLAY all captured returns and successful output bytes matched\n");return 0;
 } catch(const std::exception& e) {std::fprintf(stderr,"REPLAY error: %s\n",e.what());return 3;}
}
