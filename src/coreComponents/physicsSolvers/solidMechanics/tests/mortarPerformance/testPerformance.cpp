/* SPDX-License-Identifier: LGPL-2.1-only
 * Standalone production-kernel tests. Threaded collectives below are explicitly
 * an EMULATOR for scheduler/control-flow tests, not an MPI transport test.
 */
#include "../../kernels/MPMBoundaryMortarQR.hpp"
#include "../../kernels/MPMBoundaryMortarComponents.hpp"
#include "../../kernels/MPMMortarSearch.hpp"
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <mutex>
#include <random>
#include <thread>

namespace bm = geos::mpm::boundaryMortar;
using Clock = std::chrono::steady_clock;
int checks = 0;
void check(bool c, char const * message) { if(!c) throw std::runtime_error(message); ++checks; }
double maxdiff(std::vector<double> const & a,std::vector<double> const & b)
{ double e=0; if(a.size()!=b.size())return 1e100; for(std::size_t i=0;i<a.size();++i)e=std::max(e,std::abs(a[i]-b[i]));return e; }

struct Fixture
{
  std::vector<double> mass{1,2,3,4,5,0,2,3}, cutoff=std::vector<double>(8,1e-20);
  std::vector<bm::Row> rows{{{{0,1,.6},{2,3,.4}},1},{{{2,3,.3},{6,7,.7}},1}};
  bm::Result solve(int c, bm::MinimumNormSolver backend) const
  {
    std::vector<double> v{1.0+c,2.1-c,-.3,.4,1.3,-.2,.1,-.1};
    std::vector<bm::Relation> bc;
    if(c==0) bc={{5,0,-1,0}};
    if(c==1) bc={{4,bm::invalidIndex,1,.15},{5,1,-1,0}};
    if(c==2) bc={{4,6,-1,.2},{5,0,-1,.1}};
    bm::Options opt;
    return bm::project(mass,v,cutoff,bc,{},rows,opt,backend);
  }
};

struct Round
{
  int root=-2, arrived=0, departed=0;
  std::size_t bytes=0;
  std::vector<char> buffer;
  std::vector<int> integers;
};
struct World
{
  explicit World(int n):size(n){}
  int size;
  std::mutex mutex;
  std::condition_variable changed;
  std::map<int,Round> rounds;
  void wait(std::unique_lock<std::mutex> & lock, std::function<bool()> const & pred)
  { if(!changed.wait_for(lock,std::chrono::seconds(15),pred))throw std::runtime_error("emulated collective timeout"); }
};
struct ThreadComm
{
  World & world; int id; mutable int sequence=0;
  int rank() const { return id; } int size() const { return world.size; }
  template<class T> void broadcast(T * data,int count,int root) const
  {
    std::unique_lock<std::mutex> lock(world.mutex);
    int const key=sequence++; Round & r=world.rounds[key];
    if(r.arrived==0){r.root=root;r.bytes=std::size_t(count)*sizeof(T);r.buffer.resize(r.bytes);}
    if(r.root!=root||r.bytes!=std::size_t(count)*sizeof(T))throw std::runtime_error("collective order/type mismatch");
    if(id==root&&r.bytes)std::memcpy(r.buffer.data(),data,r.bytes);
    ++r.arrived;world.changed.notify_all();
    world.wait(lock,[&]{return r.arrived==size();});
    if(r.bytes)std::memcpy(data,r.buffer.data(),r.bytes);
    if(++r.departed==size())world.rounds.erase(key);
  }
  int allMinimum(int value) const
  {
    std::unique_lock<std::mutex> lock(world.mutex);
    int const key=sequence++; Round & r=world.rounds[key];
    if(r.arrived==0){r.root=-1;r.integers.resize(size());}
    if(r.root!=-1)throw std::runtime_error("reduction order mismatch");
    r.integers[id]=value;++r.arrived;world.changed.notify_all();
    world.wait(lock,[&]{return r.arrived==size();});
    int result=*std::min_element(r.integers.begin(),r.integers.end());
    if(++r.departed==size())world.rounds.erase(key);
    return result;
  }
};

void testScheduler(bm::MinimumNormSolver backend)
{
  Fixture f;
  std::array<bm::Result,3> reference;
  for(int c=0;c<3;++c)reference[c]=f.solve(c,backend);
  for(int n:{1,2,3,4,8})for(int dim:{2,3})for(int workers:{1,2,3})
  {
    for(int fail:{-1,1})
    {
      World world(n);std::atomic<int> started{0};
      std::vector<std::thread> threads;
      std::vector<std::string> errors(n);
      std::vector<std::array<bm::Result,3>> outputs(n);
      int const active=std::min({n,dim,workers});
      for(int rank=0;rank<n;++rank)threads.emplace_back([&,rank]{
        try
        {
          ThreadComm comm{world,rank};bool first=true;
          outputs[rank]=bm::solveComponents(comm,dim,workers,8,2,[&](int c){
            if(first)
            {
              first=false;++started;auto deadline=Clock::now()+std::chrono::seconds(5);
              while(started.load()<active)
              { if(Clock::now()>deadline)throw std::runtime_error("component solves did not start concurrently");std::this_thread::yield(); }
            }
            if(c==fail)throw std::runtime_error("injected factor failure");
            return f.solve(c,backend);
          });
        }catch(std::exception const & e){errors[rank]=e.what();}
      });
      for(auto & t:threads)t.join();
      check(started==active,"wrong number of solve workers");
      for(int rank=0;rank<n;++rank)
      {
        if(fail>=0)
        {
          check(errors[rank].find("injected factor failure")!=std::string::npos,"failure not propagated");
          check(errors[rank]==errors[0],"ranks received different failure messages");
        }
        else
        {
          check(errors[rank].empty(),errors[rank].c_str());
          for(int c=0;c<dim;++c)
          {
            auto const & r=outputs[rank][c];
            check(r.velocity==reference[c].velocity,"scheduler changed velocities");
            check(r.mortarImpulse==reference[c].mortarImpulse,"scheduler changed impulses");
            check(r.reactionImpulse==reference[c].reactionImpulse,"scheduler changed reactions");
            check(r.accepted==reference[c].accepted&&r.rank==reference[c].rank,"scheduler metadata mismatch");
          }
        }
      }
    }
  }
  // No rows/dofs: every rank still follows the same collective path.
  World emptyWorld(4);std::vector<std::thread> emptyThreads;std::atomic<int> done{0};
  for(int r=0;r<4;++r)emptyThreads.emplace_back([&,r]{ThreadComm c{emptyWorld,r};
    bm::solveComponents(c,3,3,0,0,[](int){return bm::Result{};});++done;});
  for(auto & t:emptyThreads) t.join();
  check(done==4,"empty scheduler failed");
}

void testSearch()
{
  using Key=std::array<int,3>;using Index=int;
  std::mt19937 rng(931);
  for(int trial=0;trial<2000;++trial)
  {
    std::map<Key,std::vector<Index>> bins;
    int const n=trial%2?500:0;
    for(int i=0;i<n;++i)bins[Key{{int(rng()%3),int(rng()%3),int(rng()%3)}}].push_back(i);
    std::array<Key,27> keys{};std::size_t count=rng()%28;
    for(std::size_t k=0;k<count;++k)keys[k]={{int(rng()%4),int(rng()%4),int(rng()%4)}};
    std::vector<int> expected,actual;
    for(std::size_t k=0;k<count;++k){auto it=bins.find(keys[k]);if(it!=bins.end())expected.insert(expected.end(),it->second.begin(),it->second.end());}
    std::sort(expected.begin(),expected.end());expected.erase(std::unique(expected.begin(),expected.end()),expected.end());
    geos::mpm::mortarSearch::gatherSortedMarkers(keys,count,bins,actual);
    check(actual==expected,"bin merge changed sorted candidate sequence");
    // The exact input to the unchanged distance comparator is identical, even
    // for periodic duplicate bins and ties; fit accumulation order is preserved.
  }
}

int qrDriver(int argc,char**argv)
{
  if(argc!=5)throw std::runtime_error("--qr backend input.bin output.bin");
  std::ifstream in(argv[3],std::ios::binary);
  unsigned long long dims[2];double tolerance;
  in.read(reinterpret_cast<char*>(dims),sizeof(dims));in.read(reinterpret_cast<char*>(&tolerance),sizeof(tolerance));
  std::size_t n=dims[0],m=dims[1];std::vector<double>a(n*m),b(m);
  in.read(reinterpret_cast<char*>(a.data()),a.size()*sizeof(double));in.read(reinterpret_cast<char*>(b.data()),b.size()*sizeof(double));
  if(!in)throw std::runtime_error("truncated QR fixture");
  auto begin=Clock::now();auto r=bm::minimumNormSolver(bm::selectQRBackend(argv[2]))(a,n,m,b,tolerance);
  double seconds=std::chrono::duration<double>(Clock::now()-begin).count();
  std::ofstream out(argv[4],std::ios::binary);unsigned long long rank=r.rank;
  out.write(reinterpret_cast<char*>(&rank),sizeof(rank));out.write(reinterpret_cast<char*>(&seconds),sizeof(seconds));
  out.write(reinterpret_cast<char*>(r.correction.data()),r.correction.size()*sizeof(double));
  out.write(reinterpret_cast<char*>(r.dual.data()),r.dual.size()*sizeof(double));
  if(!out) throw std::runtime_error("could not write QR result");
  return 0;
}

int main(int argc,char**argv)
{
  try
  {
    if(argc>1&&std::string(argv[1])=="--qr")return qrDriver(argc,argv);
    std::string backend=argc>1?argv[1]:"native";
    auto solver=bm::minimumNormSolver(bm::selectQRBackend(backend));
    testSearch();testScheduler(solver);
    check(bm::componentWorker(2,3,2,3)==0,"two-rank schedule");
    bool failed=false;try{bm::selectQRBackend("bogus");}catch(...){failed=true;}check(failed,"invalid backend not rejected");
    failed=false;try{bm::componentWorker(0,3,8,4);}catch(...){failed=true;}check(failed,"invalid worker limit");
    std::cout<<"PASS "<<checks<<" search/scheduler checks; backend="<<backend
             <<"; collectives=thread-emulated (NOT MPI transport)\n";
  }catch(std::exception const & e){std::cerr<<"FAIL: "<<e.what()<<"\n";return 1;}
}
