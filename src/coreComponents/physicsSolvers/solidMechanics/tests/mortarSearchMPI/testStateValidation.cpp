/* SPDX-License-Identifier: LGPL-2.1-only */
#include "../../kernels/MPMMortarStateValidation.hpp"
#include "../../kernels/MPMBoundaryMortar.hpp"
#include <iostream>
#include <map>
#include <random>
namespace ms=geos::mpm::mortarState;
namespace bm=geos::mpm::boundaryMortar;
int checks=0;
void check(bool ok,char const*m){++checks;if(!ok)throw std::runtime_error(m);}
int main()
{
 try
 {
  ms::Lattice lattice=ms::Lattice::make({{-.08,-.08,-.08}},{{1.04,1.04,1.04}},{{14,14,14}},{{false,false,false}});
  bm::Faces faces;faces.cells=lattice.cells;
  for(int d=0;d<3;++d){faces.kind[2*d]=1;faces.controlled[2*d][d]=true;faces.low[d]=0;faces.high[d]=.96;}
  // Exhaustive valid seed keys: closure CANNOT leave the valid global lattice.
  for(int k=0;k<=14;++k)for(int j=0;j<=14;++j)for(int i=0;i<=14;++i)
  {
    std::set<ms::Key> keys{{{i,j,k}}};bm::closeBoundaryDependencies(keys,faces);
    for(auto const & key:keys)check(lattice.contains(key),"Boundary closure left global lattice");
    auto first=keys;bm::closeBoundaryDependencies(keys,faces);check(keys==first,"Closure not idempotent");
  }
  faces.kind[1]=2;faces.controlled[1][0]=true;faces.kind[3]=2;faces.controlled[3][1]=true;faces.kind[5]=2;faces.controlled[5][2]=true;
  for(int k=0;k<=14;++k)for(int j=0;j<=14;++j)for(int i=0;i<=14;++i)
  {
    std::set<ms::Key> keys{{{i,j,k}}};bm::closeBoundaryDependencies(keys,faces);
    for(auto const & key:keys)check(lattice.contains(key),"Moving closure left global lattice");
  }
  check(!lattice.contains({{-1,1,1}})&&!lattice.contains({{15,1,1}}),"Invalid seed accepted");
  bool threw=false;try{ms::Lattice::make({{0,0,0}},{{1,1,1}},{{1,4,4}},{{false,false,false}});}catch(std::exception const&){threw=true;}
  check(threw,"Too-small nonperiodic grid accepted");
  // Global formula resolves rank-local spacing subtraction differences without
  // changing actual node state; displaced/nonuniform coordinates must be flagged.
  for(int i=0;i<=14;++i)
  {
    ms::Key key{{i,2,7}};auto x=lattice.point(key);
    auto noisy=x;for(int d=0;d<3;++d)noisy[d]=std::nextafter(x[d],2.0);
    check(lattice.key(noisy)==key,"Roundoff changed canonical key");
    check(lattice.snap(noisy)==x,"Roundoff changed canonical point");
    check(lattice.mismatch(noisy)<1.e-8,"Roundoff rejected");
    noisy[0]+=.01;check(lattice.mismatch(noisy)>1.e-8,"Off-lattice node was not detected");
  }
  // Emulate request filtering with Cartesian owners. A requested valid remote
  // node can be absent even though it exists on its owner: owner used its own
  // (different) request set. Union BEFORE filtering repairs that exact case.
  for(int ranks: {1,2,3,4,7,8})
  {
    std::vector<std::set<ms::Key>> request(ranks);
    ms::Key far{{12,4,4}};request[0].insert(far);request[0].insert({{0,2,2}});
    auto owner=[&](ms::Key const & key){return std::min(ranks-1,int(key[0]*ranks/15));};
    std::vector<ms::Key> announced;for(auto const&r:request)announced.insert(announced.end(),r.begin(),r.end());
    auto united=ms::requestUnion(announced);auto closed=united;bm::closeBoundaryDependencies(closed,faces);
    std::vector<ms::Key> keys(closed.begin(),closed.end());
    std::vector<int> oldPriority(keys.size()*2,-1),newPriority(oldPriority);
    for(std::size_t i=0;i<keys.size();++i)
    {
      int o=owner(keys[i]);
      if(request[o].count(keys[i]))oldPriority[2*i]=oldPriority[2*i+1]=2;
      if(closed.count(keys[i]))newPriority[2*i]=newPriority[2*i+1]=2;
    }
    check(ms::missingRecords(keys,2,{0,1},newPriority,united).empty(),"Global request union lost a valid state");
    if(owner(far)!=0)check(!ms::missingRecords(keys,2,{0,1},oldPriority,united).empty(),"Did not expose local-request filtering bug");
    // Real absent records stay fatal; a zero-mass inactive record remains present.
    newPriority[0]=-1;auto missing=ms::missingRecords(keys,2,{0,1},newPriority,united);
    check(missing.size()==1&&missing[0].field==0,"Absent record hidden");
    newPriority[0]=1;check(ms::missingRecords(keys,2,{0,1},newPriority,united).empty(),"Inactive ghost record treated as missing");
    auto deps=bm::makeBoundaryRelations(keys,2,0,faces);check(!deps.empty(),"Fixture lost boundary dependencies");
  }
  ms::Fingerprint a,b,c;for(int i=0;i<100;++i){a.integer(i);b.integer(i);c.integer(i);}a.scalar(0);b.scalar(-0.0);c.scalar(1.e-30);
  check(a.value==b.value&&a.value!=c.value,"Operator audit did not detect discrepancy");
  auto periodic=lattice;periodic.periodic[0]=true;
  check(periodic.key({{1.04,.08,.48}})[0]==0,"Periodic endpoint not canonical");
  check(!periodic.contains({{14,2,7}}),"Duplicate periodic high key accepted");
  std::cout<<"PASS state/closure/lattice checks="<<checks<<"; request layouts=1,2,3,4,7,8 (emulated ownership; no MPI transport)\n";
 }
 catch(std::exception const&e){std::cerr<<e.what()<<"\n";return 1;}
}
