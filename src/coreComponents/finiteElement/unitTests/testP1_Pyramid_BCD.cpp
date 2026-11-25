/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file testH1_Hexahedron_Lagrange1_GaussLegendre2
 */

#include "mainInterface/initialization.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include "gtest/gtest.h"

#include <chrono>

#include "finiteElement/elementFormulations/Pk_Pyramid_BCD.hpp"

using namespace geos;
using namespace finiteElement;

template< typename POLICY >
void testKernelDriver()
{
  constexpr int numNodes = 5;
  constexpr int numQuadraturePoints = 5;

  array2d< real64 > arrN( numQuadraturePoints, numNodes );
  array2d< real64 > NtestArray( numQuadraturePoints, numNodes );

  arrayView2d< real64 > const & viewN = arrN;
  arrayView2d< real64 > const & Ntest = NtestArray;

  array2d< real64 > VDM_inv;
  VDM_inv.resize( numNodes, numNodes );
  VDM_inv = Pk_Pyramid_BCD_impl< 1 >::computeVanderMondeMatrixInverse( );

  arrayView2d< real64 const > const & VDM_inv_view = VDM_inv;

  forAll< POLICY >( 1,
                    [=] GEOS_HOST_DEVICE ( localIndex const )
  {
    for( localIndex q=0; q<numQuadraturePoints; ++q )
    {
      real64 N[numNodes] = {0};
      Pk_Pyramid_BCD_impl< 1 >::calcN( q, VDM_inv_view, N );
      for( localIndex a=0; a<numNodes; ++a )
      {
        if( fabs( N[a] )<1e-9 )
        {
          viewN( q, a ) = 0;
        }
        else
        {
          viewN( q, a ) = N[a];
        }
        Ntest[q][a] = 0.0;
      }
      Ntest[q][q] = 1.0;
    }
  } );

  forAll< serialPolicy >( 1,
                          [=] ( localIndex const )
  {
    for( localIndex q=0; q<numQuadraturePoints; ++q )
    {
      for( localIndex a=0; a<numNodes; ++a )
      {
        EXPECT_FLOAT_EQ( Ntest[q][a], viewN[q][a] );
      }
    }
  } );

  //Test on mass matrix
  array2d< real64 > MtestArray( numNodes, numNodes );
  arrayView2d< real64 > const & Mtest = MtestArray;
  Pk_Pyramid_BCD_impl< 1 >::computeMassTerm( VDM_inv_view, [=] GEOS_HOST_DEVICE ( const localIndex i, const localIndex j, const real64 Mij )
  {
    Mtest[i][j] = Mij;
  } );

  real64 sum = 0.0;
  for( localIndex i=0; i<numNodes; ++i )
  {
    for( localIndex j=0; j<numNodes; ++j )
    {
      sum+= Mtest[i][j];
    }
  }

  EXPECT_FLOAT_EQ( sum, 1.33336 );

  //Test on gradient of shape functions
  //Coords of the support points
  real64 coords[numNodes][3] = { { -1.0, 1.0, 0.0 },
    { -1.0, -1.0, 0.0 },
    { 1.0, -1.0, 0.0 },
    { 1.0, 1.0, 0.0 },
    { 0.0, 0.0, 1.0 } };

  //Derivatives of the five shape functions in P1.

  array2d< real64 > gradPhi1testArray( numNodes, 3 );
  arrayView2d< real64 > const & gradPhi1test = gradPhi1testArray;
  array2d< real64 > gradPhi2testArray( numNodes, 3 );
  arrayView2d< real64 > const & gradPhi2test = gradPhi2testArray;
  array2d< real64 > gradPhi3testArray( numNodes, 3 );
  arrayView2d< real64 > const & gradPhi3test = gradPhi3testArray;
  array2d< real64 > gradPhi4testArray( numNodes, 3 );
  arrayView2d< real64 > const & gradPhi4test = gradPhi4testArray;
  array2d< real64 > gradPhi5testArray( numNodes, 3 );
  arrayView2d< real64 > const & gradPhi5test = gradPhi5testArray;

  real64 const epsilon = 1e-10; // Small value to avoid compairing floating point numbers directly
  for( localIndex i = 0; i < numNodes; ++i )
  {

    if( LvArray::math::abs( coords[i][2]-1.0 ) < epsilon )
    {
      gradPhi1test[i][0] = 0.0;
      gradPhi2test[i][0] = 0.0;
      gradPhi3test[i][0] = 0.0;
      gradPhi4test[i][0] = 0.0;
      gradPhi5test[i][0] = 0.0;
      gradPhi1test[i][1] = 0.0;
      gradPhi2test[i][1] = 0.0;
      gradPhi3test[i][1] = 0.0;
      gradPhi4test[i][1] = 0.0;
      gradPhi5test[i][1] = 0.0;
      gradPhi1test[i][2] = 0.0;
      gradPhi2test[i][2] = 0.0;
      gradPhi3test[i][2] = 0.0;
      gradPhi4test[i][2] = 0.0;
      gradPhi5test[i][2] = 0.0;
      continue;

    }
    real64 const fracz = 1.0 / (1-coords[i][2]);
    real64 const fracz2 = -1.0/ pow( 1.0 - coords[i][2], 2 );
    gradPhi1test[i][0] = 0.25*(-1-coords[i][1]*fracz);
    gradPhi2test[i][0] = 0.25*(-1+coords[i][1]*fracz);
    gradPhi3test[i][0] = 0.25*(1-coords[i][1]*fracz);
    gradPhi4test[i][0] = 0.25*(1+coords[i][1]*fracz);
    gradPhi5test[i][0] = 0.0;

    gradPhi1test[i][1] = 0.25*(1-coords[i][0]*fracz);
    gradPhi2test[i][1] = 0.25*(-1+coords[i][0]*fracz);
    gradPhi3test[i][1] = 0.25*(-1-coords[i][0]*fracz);
    gradPhi4test[i][1] = 0.25*(1+coords[i][0]*fracz);
    gradPhi5test[i][1] = 0.0;

    gradPhi1test[i][2] = 0.25*(-1+coords[i][0]*coords[i][1]*fracz2);
    gradPhi2test[i][2] = 0.25*(-1-coords[i][0]*coords[i][1]*fracz2);
    gradPhi3test[i][2] = 0.25*(-1+coords[i][0]*coords[i][1]*fracz2);
    gradPhi4test[i][2] = 0.25*(-1-coords[i][0]*coords[i][1]*fracz2);
    gradPhi5test[i][2] = 1.0;


  }

  real64 gradNtest[numNodes][3];
  for( localIndex i = 0; i < numNodes; ++i )
  {
    Pk_Pyramid_BCD_impl< 1 >::calcGradN( coords[i], VDM_inv_view, gradNtest );
    EXPECT_FLOAT_EQ( gradNtest[0][0], gradPhi1test[i][0] );
    EXPECT_FLOAT_EQ( gradNtest[1][0], gradPhi2test[i][0] );
    EXPECT_FLOAT_EQ( gradNtest[2][0], gradPhi3test[i][0] );
    EXPECT_FLOAT_EQ( gradNtest[3][0], gradPhi4test[i][0] );
    EXPECT_FLOAT_EQ( gradNtest[4][0], gradPhi5test[i][0] );

    EXPECT_FLOAT_EQ( gradNtest[0][1], gradPhi1test[i][1] );
    EXPECT_FLOAT_EQ( gradNtest[1][1], gradPhi2test[i][1] );
    EXPECT_FLOAT_EQ( gradNtest[2][1], gradPhi3test[i][1] );
    EXPECT_FLOAT_EQ( gradNtest[3][1], gradPhi4test[i][1] );
    EXPECT_FLOAT_EQ( gradNtest[4][1], gradPhi5test[i][1] );

    EXPECT_FLOAT_EQ( gradNtest[0][2], gradPhi1test[i][2] );
    EXPECT_FLOAT_EQ( gradNtest[1][2], gradPhi2test[i][2] );
    EXPECT_FLOAT_EQ( gradNtest[2][2], gradPhi3test[i][2] );
    EXPECT_FLOAT_EQ( gradNtest[3][2], gradPhi4test[i][2] );
    EXPECT_FLOAT_EQ( gradNtest[4][2], gradPhi5test[i][2] );
  }

  //Test on stiffness matrix
  array2d< real64 > RtestArray( numNodes, numNodes );
  arrayView2d< real64 > const & Rtest = RtestArray;
  Pk_Pyramid_BCD_impl< 1 >::computeStiffnessTerm( VDM_inv_view, [=] GEOS_HOST_DEVICE ( const localIndex i, const localIndex j, const real64 Rij )
  {
    Rtest[i][j] = Rij;  // Initialize the mass term
  } );

  sum = 0.0;
  for( localIndex i=0; i<numNodes; ++i )
  {
    for( localIndex j=0; j<numNodes; ++j )
    {
      sum+= Rtest[i][j];
    }
  }

  EXPECT_TRUE( LvArray::math::abs( sum ) < 1e-15 );

}
#ifdef GEOS_USE_DEVICE
TEST( FiniteElementShapeFunctions, testKernelCuda )
{
  testKernelDriver< geos::parallelDevicePolicy< > >();
}
#endif
TEST( FiniteElementShapeFunctions, testKernelHost )
{
  testKernelDriver< serialPolicy >();
}



using namespace geos;

int main( int argc, char * argv[] )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
