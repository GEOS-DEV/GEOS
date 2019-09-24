/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


#include "cxx-utilities/src/tests/testUtils.hpp"
#include "finiteElement/FiniteElementShapeFunctionKernel.hpp"

#include "rajaInterface/GEOS_RAJA_Interface.hpp"

#include "gtest/gtest.h"

#include <chrono>

namespace
{
int global_argc;
char** global_argv;
}

namespace geosx
{
double randomOne()
{
return ( 1 + ( 2 * static_cast<double>(rand())/static_cast<double>(RAND_MAX) - 1 )/1e2);
}

template< typename POLICY >
void testKernelDriver()
{
  localIndex const numElems = atoi(global_argv[1]);
  int const seed = atoi(global_argv[2]);
  int const numLaunches = atoi(global_argv[3]);

  localIndex const numNodes = (numElems+1)*4;

  std::cout<<"numElems = "<<numElems<<std::endl;
  std::cout<<"numLaunches = "<<numLaunches<<std::endl;
  srand(seed);

  auto t0 = std::chrono::high_resolution_clock::now();


  array2d<localIndex> elemsToNodesArray;
  elemsToNodesArray.resize(numElems,8);
  arrayView2d<localIndex const> const & elemsToNodes = elemsToNodesArray;

  array1d<R1Tensor> xArray;
  xArray.resize(numNodes);
  arrayView1d<R1Tensor const> const & X = xArray;

  array2d<real64> resultArray;
  resultArray.resize(numElems,3);
  arrayView2d<real64> const & result = resultArray;

  constexpr  real64 pCoords[3][8] = { { -1,  1, -1,  1, -1,  1, -1,  1 },
                                      { -1, -1,  1,  1, -1, -1,  1,  1 },
                                      { -1, -1, -1, -1,  1,  1,  1,  1 } };

  for( localIndex k=0 ; k<numElems ; ++k )
  {
    for( localIndex a=0 ; a<8 ; ++a )
    {
      localIndex const nodeIndex = a + k*4;
      elemsToNodesArray[k][a]= nodeIndex;
      xArray[nodeIndex][0] = pCoords[0][a] * randomOne();
      xArray[nodeIndex][1] = pCoords[1][a] * randomOne();
      xArray[nodeIndex][2] = pCoords[2][a] * randomOne() + k ;
    }
  }

  auto t1 = std::chrono::high_resolution_clock::now();
  auto setupTime = std::chrono::duration_cast<std::chrono::microseconds>( t1 - t0 ).count();
  std::cout<<"setuptime = "<<setupTime/1e6<<std::endl;

  std::chrono::high_resolution_clock::time_point t2;
  for( int iter=0 ; iter<numLaunches ; ++iter )
  {
    if( iter==1 )
    {
      t2 = std::chrono::high_resolution_clock::now();
      auto moveTime = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
      std::cout<<"movetime = "<<moveTime/1e6<<std::endl;

    }
    RAJA::forall< POLICY >( RAJA::TypedRangeSegment< localIndex >( 0, numElems ),
                            GEOSX_HOST_DEVICE_LAMBDA ( localIndex const k )
    {

      for( localIndex q=0 ; q<8 ; ++q )
      {
        real64 dNdX_data[3][8];
        FiniteElementShapeKernel::shapeFunctionDerivatives( k,
                                                            q,
                                                            elemsToNodes,
                                                            X,
                                                            dNdX_data );

        for( localIndex a=0 ; a<8 ; ++a )
        {
          for( localIndex i=0 ; i<3 ; ++i )
          {
            result(k,i) += dNdX_data[i][a];
          }
        }
      }
    });
  }
  auto t3 = std::chrono::high_resolution_clock::now();

  auto runtime = std::chrono::duration_cast<std::chrono::microseconds>( t3 - t2 ).count();
  std::cout<<"runtime = "<<runtime/1e6<<std::endl;


  RAJA::ReduceSum< parallelHostReduce, double > sum(0);

  RAJA::forall< serialPolicy >( RAJA::TypedRangeSegment< localIndex >( 0, numElems ),
                                GEOSX_LAMBDA ( localIndex const k )
  {
    for( localIndex i=0 ; i<3 ; ++i )
    {
      sum += result(k,i);
    }
  });
  std::cout<<sum.get()<<std::endl;


}

TEST( FiniteElementShapeFunctions, testKernelHost )
{
  testKernelDriver< parallelHostPolicy >();
}

#ifdef USE_CUDA

CUDA_TEST( FiniteElementShapeFunctions, testKernelCuda )
{
  testKernelDriver<parallelDevicePolicy<1024> >();
}

#endif

}

int main( int argc, char* argv[] )
{
  logger::InitializeLogger();

  testing::InitGoogleTest( &argc, argv );

  global_argc = argc;
  global_argv = new char*[static_cast<unsigned int>( global_argc )];
  for( int i = 0 ; i < argc ; ++i )
  {
    global_argv[i] = argv[i];
  }

  int const result = RUN_ALL_TESTS();

  delete[] global_argv;

  logger::FinalizeLogger();

#ifdef USE_CHAI
  chai::ArrayManager::finalize();
#endif

  return result;
}
