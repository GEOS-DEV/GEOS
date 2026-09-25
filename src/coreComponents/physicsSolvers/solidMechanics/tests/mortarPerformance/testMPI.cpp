/* SPDX-License-Identifier: LGPL-2.1-only
 * Real-MPI transport test of the production component scheduler, not a GEOS
 * partition/migration test. Build with mpicxx and MPMBoundaryMortarQR.cpp.
 */
#include "../../kernels/MPMBoundaryMortarComponents.hpp"
#include "../../kernels/MPMBoundaryMortarQR.hpp"
#include <mpi.h>
#include <iostream>
#include <type_traits>

namespace bm=geos::mpm::boundaryMortar;
struct Comm
{
  int rank() const {int r;MPI_Comm_rank(MPI_COMM_WORLD,&r);return r;}
  int size() const {int n;MPI_Comm_size(MPI_COMM_WORLD,&n);return n;}
  int allMinimum(int value) const {int result;MPI_Allreduce(&value,&result,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);return result;}
  template<class T> void broadcast(T*data,int count,int root) const
  {
    static_assert(std::is_same<T,int>::value||std::is_same<T,double>::value||std::is_same<T,char>::value,
                  "unsupported test datatype");
    MPI_Datatype type=std::is_same<T,int>::value?MPI_INT:(std::is_same<T,double>::value?MPI_DOUBLE:MPI_CHAR);
    if(MPI_Bcast(data,count,type,root,MPI_COMM_WORLD)!=MPI_SUCCESS)throw std::runtime_error("MPI_Bcast failed");
  }
};

int main(int argc,char**argv)
{
  MPI_Init(&argc,&argv);
  Comm comm;
  try
  {
    auto backend=bm::minimumNormSolver(bm::selectQRBackend(argc>1?argv[1]:"native"));
    auto solve=[&](int c)
    {
      return bm::project({2,3,4,1},{1.0+c,-2.0+c,.5,-.5},{0,0,0,0}, {},{},
        {{{{0,1,.4},{2,3,.6}},1}}, bm::Options{}, backend);
    };
    int checks=0;
    for(int dimensions:{2,3})for(int workers:{1,2,3})
    {
      for(int fail:{-1,1})
      {
        bool failed=false;std::string message;
        std::array<bm::Result,3> result;
        std::array<int,3> calls{{0,0,0}};
        try
        {
          result=bm::solveComponents(comm,dimensions,workers,4,1,[&](int c){
            ++calls[c];if(c==fail)throw std::runtime_error("injected worker failure");return solve(c);});
        }
        catch(std::exception const & e){failed=true;message=e.what();}
        int success=fail<0?!failed:(failed&&message.find("injected worker failure")!=std::string::npos);
        if(!comm.allMinimum(success))throw std::runtime_error("collective failure propagation mismatch");
        for(int c=0;c<dimensions;++c)
        {
          int globalCalls=0;MPI_Allreduce(&calls[c],&globalCalls,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
          if(globalCalls!=1)throw std::runtime_error("component computed more or less than once");
          if(fail<0)
          {
            auto expected=solve(c);
            int same=result[c].velocity==expected.velocity && result[c].mortarImpulse==expected.mortarImpulse &&
                     result[c].reactionImpulse==expected.reactionImpulse && result[c].rank==expected.rank &&
                     result[c].accepted==expected.accepted;
            if(!comm.allMinimum(same))throw std::runtime_error("MPI result differs from serial reference");
          }
        }
        ++checks;
      }
      auto empty=bm::solveComponents(comm,dimensions,workers,0,0,[](int){return bm::Result{};});
      if(!empty[0].velocity.empty())throw std::runtime_error("empty result mismatch");
      ++checks;
    }
    if(comm.rank()==0)std::cout<<"PASS real MPI scheduler: ranks="<<comm.size()<<", cases="<<checks<<"\n";
  }
  catch(std::exception const & e)
  {
    std::cerr<<"rank "<<comm.rank()<<": "<<e.what()<<"\n";
    MPI_Abort(MPI_COMM_WORLD,1);
    return 1;
  }
  MPI_Finalize();return 0;
}
