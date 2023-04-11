/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file testH1_Hexahedron_Lagrange2_GaussLegendre3
 */

#include "common/GEOS_RAJA_Interface.hpp"

#include "gtest/gtest.h"

#include <chrono>

#include "finiteElement/elementFormulations/H1_Hexahedron_Lagrange2_GaussLegendre3.hpp"

using namespace geosx;
using namespace finiteElement;

template< typename POLICY >
void testKernelDriver()
{
  constexpr int numNodes = 27;
  constexpr int numQuadraturePoints = 27;

  array1d< real64 > arrDetJ( numQuadraturePoints );
  array2d< real64 > arrN( numQuadraturePoints, numNodes );
  array3d< real64 > arrdNdX( numQuadraturePoints, numNodes, 3 );

  arrayView1d< real64 > const & viewDetJ = arrDetJ;
  arrayView2d< real64 > const & viewN = arrN;
  arrayView3d< real64 > const & viewdNdX = arrdNdX;

  constexpr real64 xCoords[numNodes][3] = {
    { -1.1, -1.3, -1.1 },
    { 0.1, -1.2, -1.15 },
    {  1.3, -1.1, -1.2 },
    { -1.15, -0.1, -1.1 },
    { 0.025, -0.025, -1.175 },
    { 1.2, 0.05, -1.25 },
    { -1.2, 1.1, -1.1 },
    { -0.05, 1.15, -1.2 },
    {  1.1, 1.2, -1.3 },
    { -1.2, -1.25, 0.0 },
    { 0.0, -1.225, 0.0 },
    { 1.2, -1.2, 0.0 },
    { -1.2, -0.05, 0.05 },
    { -0.0125, -0.0375, 0.0 },
    { 1.175, -0.025, -0.05 },
    { -1.2, 1.15, 0.1 },
    { -0.025, 1.15, 0.0 },
    { 1.15, 1.15, -0.1 },
    { -1.3, -1.2, 1.1 },
    { -0.1, -1.25, 1.15 },
    {  1.1, -1.3, 1.2 },
    { -1.25, 0.0, 1.2 },
    { -0.05, -0.05, 1.175 },
    { 1.15, -0.1, 1.15 },
    { -1.2, 1.2, 1.3 },
    { 0.0, 1.15, 1.2 },
    {  1.2, 1.1, 1.1 }
  };

  constexpr real64 xCoordsLin[8][3] = {
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

    for( localIndex q=0; q<numQuadraturePoints; ++q )
    {
      real64 N[numNodes] = {0};
      H1_Hexahedron_Lagrange2_GaussLegendre3::calcN( q, N );
      for( localIndex a=0; a<numNodes; ++a )
      {
        viewN( q, a ) = N[a];
      }
    }
  } );

  forAll< POLICY >( 1,
                    [=] GEOSX_HOST_DEVICE ( localIndex const )
  {

    for( localIndex q=0; q<numQuadraturePoints; ++q )
    {
      real64 dNdX[numNodes][3] = {{0}};
      viewDetJ[q] = H1_Hexahedron_Lagrange2_GaussLegendre3::calcGradN( q,
                                                                       xCoords,
                                                                       dNdX );


      for( localIndex a=0; a<numNodes; ++a )
      {
        for( int i = 0; i < 3; ++i )
        {
          viewdNdX( q, a, i ) = dNdX[a][i];
        }
      }
    }
  } );


  constexpr real64 pCoords[3][numNodes] = {
    {-1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1},
    {-1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1},
  };

  constexpr real64 pCoordsLin[3][8] = {
    { -1, 1, -1, 1, -1, 1, -1, 1 },
    { -1, -1, 1, 1, -1, -1, 1, 1 },
    { -1, -1, -1, -1, 1, 1, 1, 1 }
  };

  constexpr static real64 quadratureFactor = 0.774596669241483;

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

        real64 N = 0.0;
        real64 xPoly = 0.0;
        real64 yPoly = 0.0;
        real64 zPoly = 0.0;

        switch((int) pCoords[0][a] )
        {
          case -1:
            xPoly += 0.5*xi[0]*(xi[0]-1.0);
            break;
          case 0:
            xPoly += 1.0-(xi[0]*xi[0]);
            break;
          case 1:
            xPoly = 0.5*xi[0]*(xi[0]+1.0);
            break;
        }

        switch((int) pCoords[1][a] )
        {
          case -1:
            yPoly += 0.5*xi[1]*(xi[1]-1.0);
            break;
          case 0:
            yPoly += 1.0-(xi[1]*xi[1]);
            break;
          case 1:
            yPoly = 0.5*xi[1]*(xi[1]+1.0);
            break;
        }

        switch((int) pCoords[2][a] )
        {
          case -1:
            zPoly += 0.5*xi[2]*(xi[2]-1.0);
            break;
          case 0:
            zPoly += 1.0-(xi[2]*xi[2]);
            break;
          case 1:
            zPoly = 0.5*xi[2]*(xi[2]+1.0);
            break;
        }

        N = xPoly*yPoly*zPoly;

        EXPECT_FLOAT_EQ( N, viewN[q][a] );
      }

      real64 J[3][3] = {{0}};
      for( localIndex a=0; a<8; ++a )
      {
        real64 dNdXi[3] = { 0.125 * pCoordsLin[ 0 ][ a ] *
                            ( 1 + xi[ 1 ] * pCoordsLin[ 1 ][ a ] ) *
                            ( 1 + xi[ 2 ] * pCoordsLin[ 2 ][ a ] ),
                            0.125 * ( 1 + xi[ 0 ] * pCoordsLin[ 0 ][ a ] ) *
                            pCoordsLin[ 1 ][ a ] *
                            ( 1 + xi[ 2 ] * pCoordsLin[ 2 ][ a ] ),
                            0.125 * ( 1 + xi[ 0 ] * pCoordsLin[ 0 ][ a ] ) *
                            ( 1 + xi[ 1 ] * pCoordsLin[ 1 ][ a ] ) *
                            pCoordsLin[ 2 ][ a ] };
        for( int i = 0; i < 3; ++i )
        {
          for( int j = 0; j < 3; ++j )
          {
            J[i][j] = J[i][j] + xCoordsLin[a][i] * dNdXi[j];
          }
        }
      }

      real64 const detJ = LvArray::tensorOps::invert< 3 >( J );
      EXPECT_FLOAT_EQ( detJ, viewDetJ[q] );

      for( localIndex a=0; a<numNodes; ++a )
      {

        real64 xPoly = 0.0;
        real64 xDeriv = 0.0;
        real64 yPoly = 0.0;
        real64 yDeriv = 0.0;
        real64 zPoly = 0.0;
        real64 zDeriv = 0.0;

        switch((int) pCoords[0][a] )
        {
          case -1:
            xPoly += 0.5*xi[0]*(xi[0]-1.0);
            xDeriv += xi[0] - 0.5;
            break;
          case 0:
            xPoly += 1.0-(xi[0]*xi[0]);
            xDeriv += -2.0*xi[0];
            break;
          case 1:
            xPoly += 0.5*xi[0]*(xi[0]+1.0);
            xDeriv += xi[0] + 0.5;
            break;
        }

        switch((int) pCoords[1][a] )
        {
          case -1:
            yPoly += 0.5*xi[1]*(xi[1]-1.0);
            yDeriv += xi[1] - 0.5;
            break;
          case 0:
            yPoly += 1.0-(xi[1]*xi[1]);
            yDeriv += -2.0*xi[1];
            break;
          case 1:
            yPoly += 0.5*xi[1]*(xi[1]+1.0);
            yDeriv += xi[1] + 0.5;
            break;
        }

        switch((int) pCoords[2][a] )
        {
          case -1:
            zPoly += 0.5*xi[2]*(xi[2]-1.0);
            zDeriv += xi[2] - 0.5;
            break;
          case 0:
            zPoly += 1.0-(xi[2]*xi[2]);
            zDeriv += -2.0*xi[2];
            break;
          case 1:
            zPoly += 0.5*xi[2]*(xi[2]+1.0);
            zDeriv += xi[2] + 0.5;
            break;
        }

        real64 dNdXi[3] = {xDeriv *yPoly *zPoly, xPoly*yDeriv*zPoly, xPoly*yPoly*zDeriv};

        real64 dNdX[3] = {0};

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
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
