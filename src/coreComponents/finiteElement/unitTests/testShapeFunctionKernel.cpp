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

#include "finiteElement/FiniteElementShapeFunctionKernel.hpp"
#include "managers/initialization.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

#include "gtest/gtest.h"

#include <chrono>

using namespace geosx;

static real64 inverse( real64 (& J)[3][3] )
{
  real64 scratch[3][3];
  scratch[0][0] = J[1][1]*J[2][2] - J[1][2]*J[2][1];
  scratch[1][0] = J[0][2]*J[2][1] - J[0][1]*J[2][2];
  scratch[2][0] = J[0][1]*J[1][2] - J[0][2]*J[1][1];
  scratch[0][1] = J[1][2]*J[2][0] - J[1][0]*J[2][2];
  scratch[1][1] = J[0][0]*J[2][2] - J[0][2]*J[2][0];
  scratch[2][1] = J[0][2]*J[1][0] - J[0][0]*J[1][2];
  scratch[0][2] = J[1][0]*J[2][1] - J[1][1]*J[2][0];
  scratch[1][2] = J[0][1]*J[2][0] - J[0][0]*J[2][1];
  scratch[2][2] = J[0][0]*J[1][1] - J[0][1]*J[1][0];
  real64 invDet = 1 / ( J[0][0] * scratch[0][0] + J[1][0] * scratch[1][0] + J[2][0] * scratch[2][0] );

  for( int i=0; i<3; ++i )
  {
    for( int j=0; j<3; ++j )
    {
      J[i][j] = scratch[j][i] * invDet;
    }
  }

  return invDet;
}


template< typename POLICY >
void testKernelDriver()
{
  constexpr int numNodes = 8;
  constexpr int numQuadraturePoints = 8;

  array1d< real64 > arrDetJ( numQuadraturePoints );
  array2d< real64 > arrN( numQuadraturePoints, numNodes );
  array3d< real64 > arrdNdX( numQuadraturePoints, numNodes, 3 );

  arrayView1d< real64 > const & viewDetJ = arrDetJ;
  arrayView2d< real64 > const & viewN = arrN;
  arrayView3d< real64 > const & viewdNdX = arrdNdX;

  constexpr real64 xCoords[numNodes][3] = {
    { -1.1, -1.3, -1.1 },
    {  1.3, -1.1, -1.2 },
    { -1.2, 1.1, -1.1 },
    {  1.1, 1.2, -1.3 },
    { -1.3, -1.2, 1.1 },
    {  1.1, -1.3, 1.2 },
    { -1.2, 1.2, 1.3 },
    {  1.2, 1.1, 1.1 }
  };

  forAll< POLICY >( 1,
                    [=] GEOSX_HOST_DEVICE ( localIndex const )
  {
    real64 N[numQuadraturePoints][numNodes] = {{0}};
    real64 dNdX[numQuadraturePoints][numNodes][3] = {{{0}}};

    for( localIndex q=0; q<numQuadraturePoints; ++q )
    {
      FiniteElementShapeKernel::shapeFunctionValues( q, N[q] );
      viewDetJ[q] = FiniteElementShapeKernel::shapeFunctionDerivatives( q,
                                                                        xCoords,
                                                                        dNdX[q] );

      for( localIndex a=0; a<numNodes; ++a )
      {
        viewN( q, a ) = N[q][a];
        for( int i = 0; i < 3; ++i )
        {
          viewdNdX( q, a, i ) = dNdX[q][a][i];
        }
      }
    }
  } );


  constexpr real64 pCoords[3][numNodes] = {
    { -1, 1, -1, 1, -1, 1, -1, 1 },
    { -1, -1, 1, 1, -1, -1, 1, 1 },
    { -1, -1, -1, -1, 1, 1, 1, 1 }
  };

  constexpr static real64 quadratureFactor = 1.0 / 1.732050807568877293528;

  forAll< serialPolicy >( 1,
                          [=] ( localIndex const )
  {
    for( localIndex q=0; q<numQuadraturePoints; ++q )
    {
      real64 const xi[3] = { quadratureFactor *pCoords[0][q],
                             quadratureFactor*pCoords[1][q],
                             quadratureFactor*pCoords[2][q] };

      for( localIndex a=0; a<numNodes; ++a )
      {
        real64 N = 0.125 * ( 1 + xi[ 0 ]*pCoords[ 0 ][ a ] ) *
                   ( 1 + xi[ 1 ]*pCoords[ 1 ][ a ] ) *
                   ( 1 + xi[ 2 ]*pCoords[ 2 ][ a ] );
        EXPECT_FLOAT_EQ( N, viewN[q][a] );
      }

      real64 J[3][3] = {{0}};
      for( localIndex a=0; a<numNodes; ++a )
      {
        real64 dNdXi[3] = { 0.125 * pCoords[ 0 ][ a ] *
                            ( 1 + xi[ 1 ] * pCoords[ 1 ][ a ] ) *
                            ( 1 + xi[ 2 ] * pCoords[ 2 ][ a ] ),
                            0.125 * ( 1 + xi[ 0 ] * pCoords[ 0 ][ a ] ) *
                            pCoords[ 1 ][ a ] *
                            ( 1 + xi[ 2 ] * pCoords[ 2 ][ a ] ),
                            0.125 * ( 1 + xi[ 0 ] * pCoords[ 0 ][ a ] ) *
                            ( 1 + xi[ 1 ] * pCoords[ 1 ][ a ] ) *
                            pCoords[ 2 ][ a ] };
        for( int i = 0; i < 3; ++i )
        {
          for( int j = 0; j < 3; ++j )
          {
            J[i][j] = J[i][j] + xCoords[a][i] * dNdXi[j];
          }
        }
      }
      real64 const detJ = 1/inverse( J );
      EXPECT_FLOAT_EQ( detJ, viewDetJ[q] );

      for( localIndex a=0; a<numNodes; ++a )
      {
        real64 dNdX[3] = {0};
        real64 dNdXi[3] = { 0.125 * pCoords[ 0 ][ a ] *
                            ( 1 + xi[ 1 ] * pCoords[ 1 ][ a ] ) *
                            ( 1 + xi[ 2 ] * pCoords[ 2 ][ a ] ),
                            0.125 * ( 1 + xi[ 0 ] * pCoords[ 0 ][ a ] ) *
                            pCoords[ 1 ][ a ] *
                            ( 1 + xi[ 2 ] * pCoords[ 2 ][ a ] ),
                            0.125 * ( 1 + xi[ 0 ] * pCoords[ 0 ][ a ] ) *
                            ( 1 + xi[ 1 ] * pCoords[ 1 ][ a ] ) *
                            pCoords[ 2 ][ a ] };

        for( int i = 0; i < 3; ++i )
        {
          for( int j = 0; j < 3; ++j )
          {
            dNdX[i] += dNdXi[j] * J[j][i];
          }
        }

        EXPECT_FLOAT_EQ( dNdX[0], viewdNdX[q][a][0] );
        EXPECT_FLOAT_EQ( dNdX[1], viewdNdX[q][a][1] );
        EXPECT_FLOAT_EQ( dNdX[2], viewdNdX[q][a][2] );

      }

    }
  } );
}


#ifdef USE_CUDA
TEST( FiniteElementShapeFunctions, testKernelCuda )
{
  testKernelDriver< geosx::parallelDevicePolicy< 32 > >();
}
#endif
TEST( FiniteElementShapeFunctions, testKernelHost )
{
  testKernelDriver< serialPolicy >();
}



using namespace geosx;
int main( int argc, char * argv[] )
{

  basicSetup( argc, argv, false );

  int const result = RUN_ALL_TESTS();

  basicCleanup();

  return result;
}
