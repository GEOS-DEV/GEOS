/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file testQ5_Hexahedron_Lagrange_GaussLobatto
 */

#include "mainInterface/initialization.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include "gtest/gtest.h"

#include <chrono>

#include "finiteElement/elementFormulations/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"

using namespace geos;
using namespace finiteElement;

template< typename POLICY >
void testKernelDriver()
{
  constexpr int numNodes = 216;
  constexpr int numQuadraturePoints = 216;

  array1d< real64 > arrDetJ( numQuadraturePoints );
  array2d< real64 > arrN( numQuadraturePoints, numNodes );
  array3d< real64 > arrdNdX( numQuadraturePoints, numNodes, 3 );
  array3d< real64 > arrdNdXcheck( numQuadraturePoints, numNodes, 3 );
  array2d< real64 > Ntestarray( numQuadraturePoints, numNodes );

  arrayView1d< real64 > const & viewDetJ = arrDetJ;
  arrayView2d< real64 > const & viewN = arrN;
  arrayView3d< real64 > const & viewdNdX = arrdNdX;
  arrayView3d< real64 > const & viewdNdXcheck = arrdNdXcheck;
  arrayView2d< real64 > const & Ntest = Ntestarray;

  forAll< POLICY >( 1,
                    [=] GEOS_HOST_DEVICE ( localIndex const )
  {
    for( localIndex q=0; q<numQuadraturePoints; ++q )
    {
      real64 N[numNodes] = {0};
      Q5_Hexahedron_Lagrange_GaussLobatto::calcN( q, N );
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
      Ntest[q][q]=1.0;
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


  array2d< real64 > xCoordsData( numNodes, 3 );
  arrayView2d< real64 > const & xCoords = xCoordsData;
  xCoords[0][0]=-1.0;
  xCoords[1][0]=-sqrt( 1.0/21.0*(7.0+2.0*sqrt( 7.0 )));
  xCoords[2][0]=-sqrt( 1.0/21.0*(7.0-2.0*sqrt( 7.0 )));
  xCoords[3][0]=sqrt( 1.0/21.0*(7.0-2.0*sqrt( 7.0 )));
  xCoords[4][0]=sqrt( 1.0/21.0*(7.0+2.0*sqrt( 7.0 )));
  xCoords[5][0]=1.0;

  xCoords[0][1]=-1.0;
  xCoords[6][1]=-sqrt( 1.0/21.0*(7.0+2.0*sqrt( 7.0 )));
  xCoords[12][1]=-sqrt( 1.0/21.0*(7.0-2.0*sqrt( 7.0 )));
  xCoords[18][1]=sqrt( 1.0/21.0*(7.0-2.0*sqrt( 7.0 )));
  xCoords[24][1]=sqrt( 1.0/21.0*(7.0+2.0*sqrt( 7.0 )));
  xCoords[30][1]=1.0;

  xCoords[0][2]=-1.0;
  xCoords[36][2]=-sqrt( 1.0/21.0*(7.0+2.0*sqrt( 7.0 )));
  xCoords[72][2]=-sqrt( 1.0/21.0*(7.0-2.0*sqrt( 7.0 )));
  xCoords[108][2]=sqrt( 1.0/21.0*(7.0-2.0*sqrt( 7.0 )));
  xCoords[144][2]=sqrt( 1.0/21.0*(7.0+2.0*sqrt( 7.0 )));
  xCoords[180][2]=1.0;

  for( localIndex k=0; k<6; ++k )
  {
    for( localIndex j=0; j<6; ++j )
    {
      for( localIndex i=0; i<6; ++i )
      {
        xCoords[i+6*j+36*k][0]=xCoords[i][0];
        xCoords[i+6*j+36*k][1]=xCoords[6*j][1];
        xCoords[i+6*j+36*k][2]=xCoords[36*k][2];
      }
    }
  }

  array2d< real64 > gradNxtestArray( numNodes, numQuadraturePoints );
  array2d< real64 > gradNytestArray( numNodes, numQuadraturePoints );
  array2d< real64 > gradNztestArray( numNodes, numQuadraturePoints );
  arrayView2d< real64 > const & gradNxtest = gradNxtestArray;
  arrayView2d< real64 > const & gradNytest = gradNytestArray;
  arrayView2d< real64 > const & gradNztest = gradNztestArray;

  gradNxtest[0][0]=-15.0/2.0;
  gradNxtest[0][1]=1.0/6.0*(-2.0-sqrt( 7.0 )-sqrt( 21.0+6.0*sqrt( 7.0 )));
  gradNxtest[0][2]=1.0/6.0*(-2.0+sqrt( 7.0 )+sqrt( 21.0-6.0*sqrt( 7.0 )));
  gradNxtest[0][3]=1.0/6.0*(-2.0+sqrt( 7.0 )-sqrt( 21.0-6.0*sqrt( 7.0 )));
  gradNxtest[0][4]=-1.0/6.0*(2.0+sqrt( 7.0 )-sqrt( 21.0+6.0*sqrt( 7.0 )));
  gradNxtest[0][5]= -1.0/2.0;

  gradNxtest[1][0]=1.0/4.0*(7.0+4.0*sqrt( 7.0 )+sqrt( 343.0+70.0*sqrt( 7.0 )));
  gradNxtest[1][1]=0.0;
  gradNxtest[1][2]=-1.0/24.0*sqrt( 28.0-2.0*sqrt( 7.0 ))*(5.0+sqrt( 3.0 )-sqrt( 7.0 )+sqrt( 21.0 ));
  gradNxtest[1][3]=1.0/24.0*sqrt( 28.0-2.0*sqrt( 7.0 ))*(-5.0+sqrt( 3.0 )+sqrt( 7.0 )+sqrt( 21.0 ));
  gradNxtest[1][4]=-1.0/2.0*sqrt( 7-2.0*sqrt( 7.0 ));
  gradNxtest[1][5]=1.0/4.0*(-7.0-4.0*sqrt( 7.0 )+sqrt( 343.0+70.0*sqrt( 7.0 )));

  gradNxtest[2][0]=1.0/4.0*(7.0-4.0*sqrt( 7.0 )-sqrt( 343.0-70.0*sqrt( 7.0 )));
  gradNxtest[2][1]=1.0/12.0*sqrt( 7.0+sqrt( 7.0 )/2.0 )*(5.0-sqrt( 3.0 )+sqrt( 7.0 )+sqrt( 21.0 ));
  gradNxtest[2][2]=0.0;
  gradNxtest[2][3]=-1.0/2.0*sqrt( 7.0+2.0*sqrt( 7.0 ));
  gradNxtest[2][4]=1.0/12.0*sqrt( 7.0+sqrt( 7.0 )/2.0 )*(5.0+sqrt( 3.0 )+sqrt( 7.0 )-sqrt( 21.0 ));
  gradNxtest[2][5]=1.0/4.0*(-7.0+4.0*sqrt( 7.0 )-sqrt( 343.0-70.0*sqrt( 7.0 )));

  gradNxtest[3][0]=-gradNxtest[2][5];
  gradNxtest[3][1]=-gradNxtest[2][4];
  gradNxtest[3][2]=-gradNxtest[2][3];
  gradNxtest[3][3]=-gradNxtest[2][2];
  gradNxtest[3][4]=-gradNxtest[2][1];
  gradNxtest[3][5]=-gradNxtest[2][0];

  gradNxtest[4][0]=-gradNxtest[1][5];
  gradNxtest[4][1]=-gradNxtest[1][4];
  gradNxtest[4][2]=-gradNxtest[1][3];
  gradNxtest[4][3]=-gradNxtest[1][2];
  gradNxtest[4][4]=-gradNxtest[1][1];
  gradNxtest[4][5]=-gradNxtest[1][0];

  gradNxtest[5][0]=-gradNxtest[0][5];
  gradNxtest[5][1]=-gradNxtest[0][4];
  gradNxtest[5][2]=-gradNxtest[0][3];
  gradNxtest[5][3]=-gradNxtest[0][2];
  gradNxtest[5][4]=-gradNxtest[0][1];
  gradNxtest[5][5]=-gradNxtest[0][0];

  for( localIndex k=0; k<6; ++k )
  {
    for( localIndex j=0; j<6; ++j )
    {
      for( localIndex i=0; i<6; ++i )
      {
        for( localIndex l=0; l<6; ++l )
        {
          gradNxtest[i+6*j+36*k][l+6*j+36*k]=gradNxtest[i][l];
          gradNxtest[i+6*j+36*k][l+6*j+36*k]=gradNxtest[i][l];

          gradNytest[i+6*j+36*k][i+6*l+36*k]=gradNxtest[j][l];
          gradNytest[i+6*j+36*k][i+6*l+36*k]=gradNxtest[j][l];

          gradNztest[i+6*j+36*k][i+6*j+36*l]=gradNxtest[k][l];
          gradNztest[i+6*j+36*k][i+6*j+36*l]=gradNxtest[k][l];
        }
      }
    }
  }

  forAll< POLICY >( 1,
                    [=] GEOS_HOST_DEVICE ( localIndex const )
  {

    real64 xLocal[numNodes][3];

    for( localIndex a=0; a< numNodes; ++a )
    {
      for( localIndex i=0; i<3; ++i )
      {
        xLocal[a][i] = xCoords[a][i];
      }
    }

    for( localIndex q=0; q<numQuadraturePoints; ++q )
    {

      real64 dNdX[numNodes][3] = {{0}};
      real64 dNdXcheck[numNodes][3] = {{0}};
      // check the explicit calculation of gradient values
      viewDetJ[q] = Q5_Hexahedron_Lagrange_GaussLobatto::calcGradN( xLocal[ q ],
                                                                    xLocal,
                                                                    dNdXcheck );
      viewDetJ[q] = Q5_Hexahedron_Lagrange_GaussLobatto::calcGradN( q,
                                                                    xLocal,
                                                                    dNdX );

      for( localIndex a=0; a<numNodes; ++a )
      {
        for( int i = 0; i < 3; ++i )
        {
          if( fabs( dNdX[a][i] )<1e-9 )
          {
            viewdNdX( q, a, i ) = 0;
          }
          else
          {
            viewdNdX( q, a, i ) = dNdX[a][i];
          }
        }
      }
      for( localIndex a=0; a<numNodes; ++a )
      {
        for( int i = 0; i < 3; ++i )
        {
          if( fabs( dNdXcheck[a][i] )<1e-9 )
          {
            viewdNdXcheck( q, a, i ) = 0;
          }
          else
          {
            viewdNdXcheck( q, a, i ) = dNdXcheck[a][i];
          }
        }
      }
    }
  } );

  forAll< serialPolicy >( 1,
                          [=] ( localIndex const )
  {

    for( localIndex q=0; q<numQuadraturePoints; ++q )
    {
      for( localIndex a=0; a<numNodes; ++a )
      {
        for( int i = 0; i < 3; ++i )
        {
          EXPECT_NEAR( viewdNdXcheck( q, a, i ), viewdNdX( q, a, i ), 1e-9 );
        }
        EXPECT_FLOAT_EQ( gradNxtest[a][q], viewdNdX( q, a, 0 ) );
        EXPECT_FLOAT_EQ( gradNytest[a][q], viewdNdX( q, a, 1 ) );
        EXPECT_FLOAT_EQ( gradNztest[a][q], viewdNdX( q, a, 2 ) );
      }
    }
  } );

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
//int main( int argc, char * argv[] )
//{
//  testing::InitGoogleTest();
//
//  int const result = RUN_ALL_TESTS();
//
//  return result;
//}
//
int main( int argc, char * argv[] )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
