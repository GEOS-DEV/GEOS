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
 * @file testH1_Tetrahedron_Lagrange2_Gauss4.cpp
 */

#include "finiteElement/elementFormulations/H1_Tetrahedron_Lagrange2_Gauss4.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include "gtest/gtest.h"

#include <chrono>

using namespace geosx;
using namespace finiteElement;

template< typename POLICY >
void testKernelDriver()
{
  constexpr int numNodes = 10;
  constexpr int numQuadraturePoints = 4;
  //constexpr real64 weight = 1.0 / 24.0;

  array1d< real64 > arrDetJxW( numQuadraturePoints );
  array2d< real64 > arrN( numQuadraturePoints, numNodes );
  array3d< real64 > arrdNdX( numQuadraturePoints, numNodes, 3 );

  arrayView1d< real64 > const & viewDetJxW = arrDetJxW;
  arrayView2d< real64 > const & viewN = arrN;
  arrayView3d< real64 > const & viewdNdX = arrdNdX;

  //constexpr real64 xCoordsLin[4][3] = {
  //  { 1.0, 0.0, 0.0 },
  //  { 2.0, 0.0, 0.0 },
  //  { 0.0, 4.0, -0.5 },
  //  { 1.0, 1.0, 2.0 }
  //};

  constexpr real64 xCoords[numNodes][3] = {
    { 1.0, 0.0, 0.0 },
    { 2.0, 0.0, 0.0 },
    { 0.0, 4.0, -0.5 },
    { 1.0, 1.0, 2.0 },
    { 1.5, 0.0, 0.0 },
    { 1.0, 2.0, -0.25},
    { 0.5, 2.0, -0.25 },
    { 1.0, 0.5, 1.0 },
    { 1.5, 0.5, 1.0 },
    { 0.5, 2.5, 0.75 }
  };

  forAll< POLICY >( 1,
                    [=] GEOSX_HOST_DEVICE ( localIndex const )
  {

    for( localIndex q=0; q<numQuadraturePoints; ++q )
    {
      real64 N[numNodes] = {0};
      H1_Tetrahedron_Lagrange2_Gauss4::calcN( q, N );
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
      viewDetJxW[q] = H1_Tetrahedron_Lagrange2_Gauss4::calcGradN( q,
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

  constexpr real64 N_ref[numNodes][numQuadraturePoints] = {
    { -0.099999999496880, -0.099999999496880, -0.099999999496880, 0.100000004528080},
    { 0.100000004528080, -0.099999999496880, -0.099999999496880, -0.099999999496880},
    { -0.099999999496880, 0.100000004528080, -0.099999999496880, -0.099999999496880},
    { -0.099999999496880, -0.099999999496880, 0.100000004528080, -0.099999999496880},
    { 0.323606796981280, 0.076393201006240, 0.076393201006240, 0.323606796981280},
    { 0.323606796981280, 0.323606796981280, 0.076393201006240, 0.076393201006240},
    { 0.076393201006240, 0.323606796981280, 0.076393201006240, 0.323606796981280},
    { 0.076393201006240, 0.076393201006240, 0.323606796981280, 0.323606796981280},
    { 0.323606796981280, 0.076393201006240, 0.323606796981280, 0.076393201006240},
    { 0.076393201006240, 0.323606796981280, 0.323606796981280, 0.076393201006240}
  };

  real64 const detJxW_ref = 0.354166666666667;

  constexpr real64 dNdx_ref[numQuadraturePoints][numNodes][3] = {
    { {0.447213600000000, 0.236760141176471, 0.105226729411765},
      {1.341640800000000, 0.315680188235294, -0.157840094117647},
      {0.0, -0.105226729411765, 0.052613364705882},
      {0.0, -0.026306682352941, -0.210453458823529},
      {-1.788854400000000, -1.109624800000000, -0.616008000000000},
      { 0.552786400000000, 0.681041694117647, -0.340520847058824},
      {-0.552786400000000, -0.162584235294118, -0.195101082352941},
      {-0.552786400000000, -0.260134776470588, 0.130067388235294},
      {0.552786400000000, 0.267810964705882, 1.036914917647059},
      {0.0, 0.162584235294118, 0.195101082352941} },
    { {0.447213600000000, 0.236760141176471, 0.105226729411765},
      {-0.447213600000000, -0.105226729411765, 0.052613364705882},
      {0.0, 0.315680188235294, -0.157840094117647},
      {0.0, -0.026306682352941, -0.210453458823529},
      {0.0, -0.162584235294118, -0.195101082352941},
      { 2.341640800000000, 0.681041694117647, -0.340520847058824},
      {-2.341640800000000, -1.109624800000000, -0.616008000000000},
      {-0.552786400000000, -0.260134776470588, 0.130067388235294},
      {0.552786400000000, 0.162584235294118, 0.195101082352941},
      {0.0, 0.267810964705882, 1.036914917647059} },
    { {0.447213600000000, 0.236760141176471, 0.105226729411765},
      {-0.447213600000000, -0.105226729411765, 0.052613364705882},
      {0.0, -0.105226729411765, 0.052613364705882},
      {0.0, 0.078920047058824, 0.631360376470588},
      {0.0, -0.162584235294118, -0.195101082352941},
      {0.552786400000000, 0.260134776470588, -0.130067388235294},
      {-0.552786400000000, -0.162584235294118, -0.195101082352941},
      {-2.341640800000000, -1.207175341176471, -0.290839529411765},
      {2.341640800000000, 0.583491152941176, -0.015352376470588},
      {0.0, 0.583491152941176, -0.015352376470588} },
    { {-1.341640800000000, -0.710280423529412, -0.315680188235294},
      {-0.447213600000000, -0.105226729411765, 0.052613364705882},
      {0.0, -0.105226729411765, 0.052613364705882},
      {0.0, -0.026306682352941, -0.210453458823529},
      {1.788854400000000, 0.258322682352941, -0.405554541176471},
      {0.552786400000000, 0.260134776470588, -0.130067388235294},
      {-0.552786400000000, 0.258322682352941, -0.405554541176471},
      {-0.552786400000000, -0.154908047058824, 0.971881223529412},
      {0.552786400000000, 0.162584235294118, 0.195101082352941},
      {0.0, 0.162584235294118, 0.195101082352941} }
  };


  forAll< serialPolicy >( 1,
                          [=] ( localIndex const )
  {
    for( localIndex q=0; q<numQuadraturePoints; ++q )
    {

      for( localIndex a=0; a<numNodes; ++a )
      {
        EXPECT_FLOAT_EQ( N_ref[a][q], viewN[q][a] );

        for( int j = 0; j < 3; ++j )
        {
          EXPECT_FLOAT_EQ( dNdx_ref[q][a][j], viewdNdX[q][a][j] );
        }
      }

      EXPECT_FLOAT_EQ( detJxW_ref, viewDetJxW[q] );

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
