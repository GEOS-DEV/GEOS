/* SPDX-License-Identifier: LGPL-2.1-only */
// Actual MPI transport test for request union/state completeness, NOT GEOS mesh migration.
#include "../../kernels/MPMMortarStateValidation.hpp"
#include "../../kernels/MPMBoundaryMortar.hpp"
#include <mpi.h>
#include <algorithm>
#include <iostream>
#include <type_traits>
namespace ms=geos::mpm::mortarState;
namespace bm=geos::mpm::boundaryMortar;
template<class T>std::vector<T> all(std::vector<T>const & local)
{
 static_assert(std::is_trivially_copyable<T>::value,"MPI records must be trivially copyable");
 int size=0;MPI_Comm_size(MPI_COMM_WORLD,&size);int bytes=int(local.size()*sizeof(T));
 std::vector<int>counts(size),offsets(size);MPI_Allgather(&bytes,1,MPI_INT,counts.data(),1,MPI_INT,MPI_COMM_WORLD);
 for(int i=1;i<size;++i)offsets[i]=offsets[i-1]+counts[i-1];
 std::vector<T>out(std::size_t(offsets.back()+counts.back())/sizeof(T));char dummy=0;
 MPI_Allgatherv(local.empty()?static_cast<void*>(&dummy):const_cast<T*>(local.data()),bytes,MPI_BYTE,
               out.empty()?static_cast<void*>(&dummy):out.data(),counts.data(),offsets.data(),MPI_BYTE,MPI_COMM_WORLD);
 return out;
}
struct Record{ms::Key key;int field;double mass;};
int main(int argc,char**argv)
{
 MPI_Init(&argc,&argv);int rank=0,size=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
 int failure=0;
 try
 {
  std::vector<ms::Key>request;if(rank==0)request={{{12,4,4}},{{0,2,2}}};
  auto original=ms::requestUnion(all(request));auto closed=original;
  bm::Faces faces;faces.cells={{14,14,14}};for(int d=0;d<3;++d){faces.kind[2*d]=1;faces.controlled[2*d][d]=true;}
  bm::closeBoundaryDependencies(closed,faces);
  std::vector<ms::Key>keys(closed.begin(),closed.end());
  auto owner=[&](ms::Key const & key){return std::min(size-1,int(key[0]*size/15));};
  std::vector<Record>local;
  for(auto const & key:keys)if(owner(key)==rank)for(int field=0;field<2;++field)local.push_back({key,field,field?0.0:1.0});
  auto gathered=all(local);
  auto audit=[&](std::vector<Record>const & records)
  {
   std::vector<int>priorities(keys.size()*2,-1);
   for(auto const &record:records){auto i=std::lower_bound(keys.begin(),keys.end(),record.key)-keys.begin();priorities[std::size_t(i)*2+record.field]=2;}
   return ms::missingRecords(keys,2,{0,1},priorities,original);
  };
  if(!audit(gathered).empty())throw std::runtime_error("Global request union lost owner state");
  ms::Key removed{{12,4,4}};
  local.erase(std::remove_if(local.begin(),local.end(),[&](Record const&r){return r.key==removed;}),local.end());
  auto missing=audit(all(local));
  if(missing.size()!=2||missing[0].key!=removed||missing[0].addedByClosure)throw std::runtime_error("Injected true absence was hidden");
  if(rank==0)std::cout<<"PASS actual MPI request/state fixture, ranks="<<size<<"; zero-mass records retained and injected absence detected\n";
 }
 catch(std::exception const & e){std::cerr<<"rank="<<rank<<" "<<e.what()<<"\n";failure=1;}
 int any=0;MPI_Allreduce(&failure,&any,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);MPI_Finalize();return any;
}
