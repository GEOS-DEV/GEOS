/* SPDX-License-Identifier: LGPL-2.1-only */
#include "../../kernels/MPMMortarIndexedSearch.hpp"
#include "../../kernels/MPMMortarSearch.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <tuple>

using Key=std::array<long long,3>;
using Point=std::array<double,3>;
struct Marker {Point x,normal; long long id; int group;};
using Search=geos::mpm::mortarSearch::IndexedBins<Key,int>;
using Neighbor=Search::Neighbor;
int checks=0;
void check(bool ok,char const * message){++checks;if(!ok)throw std::runtime_error(message);}

struct Fixture
{
  int dims=3,count=3;
  double length=1.0;
  std::array<bool,3> periodic{{false,false,false}};
  std::vector<Marker> markers;
  std::map<Key,std::vector<int>> bins;
  Key bin(Point const & x)const
  {
    Key out{};
    for(int d=0;d<dims;++d)
    {
      auto i=static_cast<long long>(std::floor(x[d]*count/length));
      out[d]=periodic[d]?((i%count)+count)%count:std::max(0LL,std::min(static_cast<long long>(count-1),i));
    }
    return out;
  }
  Point canonical(Point x) const
  {
    for(int d=0;d<dims;++d)if(periodic[d])x[d]-=std::floor(x[d]/length)*length;
    return x;
  }
  void build(){bins.clear();for(std::size_t i=0;i<markers.size();++i)bins[bin(markers[i].x)].push_back(static_cast<int>(i));}
  std::pair<std::array<Key,27>,std::size_t> neighbors(Point const & x)const
  {
    Key b=bin(x); std::array<Key,27> out{};std::size_t n=0;
    for(int k=dims==3?-1:0;k<=(dims==3?1:0);++k)for(int j=-1;j<=1;++j)for(int i=-1;i<=1;++i)
    {
      Key key{{b[0]+i,b[1]+j,b[2]+k}};bool valid=true;
      for(int d=0;d<dims;++d)
      {
        if(periodic[d])key[d]=(key[d]%count+count)%count;
        else if(key[d]<0||key[d]>=count)valid=false;
      }
      if(valid)out[n++]=key;
    }
    return {out,n};
  }
  double distance(Point const & q,int id)const
  {
    double s=0;for(int d=0;d<dims;++d){double a=q[d]-markers[id].x[d];if(periodic[d])a-=std::round(a/length)*length;s+=a*a;}return s;
  }
  bool less(Neighbor const & a,Neighbor const & b)const
  {
    if(a.first!=b.first)return a.first<b.first;
    Marker const & l=markers[a.second], &r=markers[b.second];
    return std::tie(l.id,l.group,l.x[0],l.x[1],l.x[2])<std::tie(r.id,r.group,r.x[0],r.x[1],r.x[2]);
  }
  std::vector<Neighbor> reference(Point const & q,std::size_t k,std::size_t & visited) const
  {
    auto keys=neighbors(q);std::vector<int> ids;
    geos::mpm::mortarSearch::gatherSortedMarkers(keys.first,keys.second,bins,ids);
    visited=ids.size();std::vector<Neighbor> out;out.reserve(visited);
    for(int id:ids)out.push_back({distance(q,id),id});
    k=std::min(k,out.size());std::partial_sort(out.begin(),out.begin()+k,out.end(),[&](Neighbor const &a,Neighbor const &b){return less(a,b);});out.resize(k);return out;
  }
  // Same weighted tangent-plane calculation; compare exact accumulation order.
  std::array<double,4> fit(Point const & q,std::vector<Neighbor>const & ids)const
  {
    double wsum=0,phi=0;Point norm{};
    for(auto const & t:ids){double w=1.0/(t.first+0.0001);wsum+=w;double proj=0;
      for(int d=0;d<dims;++d){double a=q[d]-markers[t.second].x[d];if(periodic[d])a-=std::round(a/length)*length;proj+=a*markers[t.second].normal[d];norm[d]+=w*markers[t.second].normal[d];}
      phi+=w*proj;
    }
    return {{wsum==0?0:phi/wsum,norm[0],norm[1],norm[2]}};
  }
};

void tests()
{
  std::mt19937_64 rng(928172);std::uniform_real_distribution<double> uniform(-.7,1.5);
  std::size_t queries=0,replays=0,visits=0,evaluations=0;
  for(int scene=0;scene<160;++scene)
  {
    Fixture f;f.dims=scene%2?2:3;f.count=1+scene%6;
    f.length=scene%8==0?1.e-6:1.0;
    for(int d=0;d<f.dims;++d)f.periodic[d]=((scene>>(d+1))&1)!=0;
    int n=scene%13==0?0:(scene%17==0?7:400+scene*3);
    for(int i=0;i<n;++i){Marker m;for(int d=0;d<3;++d){m.x[d]=uniform(rng)*f.length;m.normal[d]=uniform(rng);}m.x=f.canonical(m.x);m.id=i;m.group=i%2;f.markers.push_back(m);}
    // Comparator-equivalent duplicates must be detected and safely replayed.
    if(scene%9==0&&n>0)for(int i=0;i<4;++i)f.markers.push_back(f.markers[std::size_t(i)%n]);
    f.build();Search tree(f.bins,[&](int id){return f.canonical(f.markers[id].x);},f.dims,f.periodic,Point{{f.length,f.length,f.length}});
    for(int query=0;query<40;++query)
    {
      Point q{{uniform(rng)*f.length,uniform(rng)*f.length,uniform(rng)*f.length}};q=f.canonical(q);
      if(query==0&&!f.markers.empty())q=f.markers[0].x;
      auto keys=f.neighbors(q);std::size_t const k=query%7==0?1000:12;
      std::size_t count=0;auto ref=f.reference(q,k,count);std::vector<Neighbor> out;Search::Statistics stats;
      bool replay=tree.query(keys.first,keys.second,q,k,[&](int id){return f.distance(q,id);},[&](Neighbor const&a,Neighbor const&b){return f.less(a,b);},out,stats);
      ++queries;visits+=count;evaluations+=stats.distances;
      if(replay){++replays;out=ref;}
      check(count==stats.candidates,"Candidate cardinality changed");
      check(out==ref,"Selected marker/distance sequence differs");
      check(f.fit(q,out)==f.fit(q,ref),"Fitted values changed");
    }
  }
  // Adversarial equal-distance shell, periodic seam, and neighbors outside R.
  Fixture f;f.count=1;
  for(int i=0;i<60;++i)f.markers.push_back(Marker{{{i%2?-1.:1.,0,0}},{{1,0,0}},60-i,0});
  f.build();Search tree(f.bins,[&](int id){return f.markers[id].x;},3,f.periodic,Point{{1,1,1}});
  Point q{{0,0,0}};auto keys=f.neighbors(q);std::size_t count;auto ref=f.reference(q,12,count);Search::Statistics stats;std::vector<Neighbor> out;
  bool replay=tree.query(keys.first,keys.second,q,12,[&](int i){return f.distance(q,i);},[&](Neighbor const&a,Neighbor const&b){return f.less(a,b);},out,stats);
  check(!replay && out==ref,"Equal-distance metadata ordering changed");
  f.markers[1]=f.markers[0];f.build();Search duplicate(f.bins,[&](int i){return f.markers[i].x;},3,f.periodic,Point{{1,1,1}});
  replay=duplicate.query(keys.first,keys.second,f.markers[0].x,60,[&](int i){return f.distance(f.markers[0].x,i);},[&](Neighbor const&a,Neighbor const&b){return f.less(a,b);},out,stats);
  check(replay,"Ambiguous comparator tie not replayed");
  std::cout<<"PASS checks="<<checks<<" queries="<<queries<<" equivalent-tie-replays="<<replays<<" legacy-candidates="<<visits<<" indexed-distances="<<evaluations<<"\n";
}

void benchmark()
{
  Fixture f; f.count=2;std::mt19937_64 rng(16243);std::uniform_real_distribution<double> u(0.0,1.0);
  // Dense symmetry-extended sphere: synthetic fixture, not the user's run data.
  for(int p=0;p<10690;++p)
  {
    double z=u(rng),t=1.5707963267948966*u(rng),rho=std::sqrt(1-z*z);
    Point x{{.22*rho*std::cos(t),.22*rho*std::sin(t),.22*z}};
    for(int image=0;image<8;++image)
    {Marker m;m.id=p;m.group=p%2;for(int d=0;d<3;++d){m.x[d]=x[d]*((image&(1<<d))?-1:1);m.normal[d]=m.x[d]/.22;}f.markers.push_back(m);}
  }
  f.build();std::vector<Point> query;
  for(int k=0;k<11;++k)for(int j=0;j<11;++j)for(int i=0;i<10;++i)query.push_back(Point{{i*.04,j*.04,k*.04}});
  using Clock=std::chrono::steady_clock;auto start=Clock::now();
  Search tree(f.bins,[&](int id){return f.markers[id].x;},3,f.periodic,Point{{1,1,1}});
  double build=std::chrono::duration<double>(Clock::now()-start).count();
  std::vector<std::vector<Neighbor>> refs;std::size_t eligible=0;start=Clock::now();
  for(auto const&q:query){std::size_t n;refs.push_back(f.reference(q,12,n));eligible+=n;}
  double legacy=std::chrono::duration<double>(Clock::now()-start).count();
  std::size_t eval=0,nodes=0,replays=0;start=Clock::now();
  for(std::size_t j=0;j<query.size();++j)
  {
    auto const&q=query[j];auto keys=f.neighbors(q);std::vector<Neighbor> out;Search::Statistics stats;
    bool replay=tree.query(keys.first,keys.second,q,12,[&](int id){return f.distance(q,id);},[&](Neighbor const&a,Neighbor const&b){return f.less(a,b);},out,stats);
    if(replay){++replays;std::size_t n;out=f.reference(q,12,n);stats.distances+=n;}
    check(out==refs[j],"Benchmark changed selected neighbors");eval+=stats.distances;nodes+=stats.treeNodes;
  }
  double indexed=std::chrono::duration<double>(Clock::now()-start).count();
  std::cout<<std::setprecision(9)<<"{\"scope\":\"synthetic symmetry-extended sphere search, not GEOS\",\"markers\":"<<f.markers.size()<<",\"queries\":"<<query.size()<<",\"legacy_candidates\":"<<eligible<<",\"indexed_distance_evaluations\":"<<eval<<",\"tree_nodes\":"<<nodes<<",\"tie_replays\":"<<replays<<",\"build_seconds\":"<<build<<",\"legacy_seconds\":"<<legacy<<",\"indexed_query_seconds\":"<<indexed<<",\"speedup_including_build\":"<<legacy/(indexed+build)<<"}\n";
}
int main(int argc,char**argv){try{if(argc>1&&std::string(argv[1])=="--benchmark")benchmark();else tests();}catch(std::exception const&e){std::cerr<<e.what()<<"\n";return 1;}}
