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
 * @file testComputationalGeometry.cpp
 */
#include "../utilities/ComputationalGeometry.hpp"
#include <gtest/gtest.h>
#include <numeric>

namespace geosx
{

TEST( testComputationalGeometry, checkCentroid3DPolygon )
{
  constexpr localIndex n = 12;
  real64 t = 2.0*M_PI/n;
  real64 s = (1.0 + sin( t ))/cos( t ) - 1.0;
  constexpr localIndex nodeDim = 3;
  constexpr localIndex numNodes = n*3;

  real64 e[3] = { -1.0e-6, 0.0, 1.0e-6 };

  array2d< real64, nodes::REFERENCE_POSITION_PERM > points;

  points.resize( numNodes, nodeDim );

  localIndex idx = 0;
  for( localIndex z = 0; z < 3; ++z )
  {
    for( localIndex i = 0; i < n; ++i )
    {
      real64 r = 25.0*(1.0 + (i%2)*s + e[z]);
      real64 x = r*cos( i*t );
      real64 y = r*sin( i*t );
      points( idx, 0 ) = x;
      points( idx, 1 ) = y;
      points( idx, 2 ) = 10.0*z;
      ++idx;
    }
  }

  constexpr real64 areaTolerance = 0.0001;

  for( localIndex ci = 0; ci < 3; ++ci )
  {
    array1d< localIndex > indices;
    for( localIndex i = 0; i < n; ++i )
    {
      indices.emplace_back( ci*n + i );
    }

    real64 faceCenter[ 3 ], faceNormal[ 3 ];
    computationalGeometry::centroid_3DPolygon( indices.toSliceConst(), points.toViewConst(), faceCenter, faceNormal, areaTolerance );

    real64 norm = LvArray::math::square( faceNormal[0] ) + LvArray::math::square( faceNormal[1] ) + LvArray::math::square( faceNormal[2] );
    EXPECT_EQ( abs( sqrt( norm ) - 1.0 ) < 0.0001, true );
  }
}

TEST( testComputationalGeometry, checkNormalsConsistency )
{
  real64 pArray[ 8 ][ 3 ] = {
    { 0.1084442213177681, 0.0, 0.0                },
    { 0.2075795382261276, 0.0, 0.0                },
    { 0.1058338209986687, 0.1052169278264046, 0.0                },
    { 0.2090811282396317, 0.100415863096714, 0.0                },
    { 0.1048644408583641, 0.0, 0.1068341732025146 },
    { 0.2063014656305313, 0.0, 0.1084721013903618 },
    { 0.1081591323018074, 0.1056273579597473, 0.1087869703769684 },
    { 0.2010060846805573, 0.1080376580357552, 0.1048369631171227 }
  };
  array2d< real64, nodes::REFERENCE_POSITION_PERM > points( 8, 3 );
  LvArray::tensorOps::copy< 8, 3 >( points, pArray );
  real64 centroid[3];
  LvArray::tensorOps::fill< 3 >( centroid, 0 );
  for( localIndex i = 0; i < 8; ++i )
  {
    LvArray::tensorOps::add< 3 >( centroid, points[ i ] );
  }
  LvArray::tensorOps::scale< 3 >( centroid, 1.0/8.0 );

  localIndex fArray[6][4] = {
    { 0, 2, 3, 1 },
    { 6, 4, 5, 7 },
    { 6, 7, 3, 2 },
    { 5, 4, 0, 1 },
    { 7, 5, 1, 3 },
    { 2, 0, 4, 6 }
  };
  array2d< localIndex > faces( 6, 4 );
  LvArray::tensorOps::copy< 6, 4 >( faces, fArray );

  real64 divergence[ 3 ]; // numerical divergences of (1,0,0), (0,1,0), (0,0,1)
  LvArray::tensorOps::fill< 3 >( divergence, 0.0 );
  for( localIndex f = 0; f < 6; ++f )
  {
    real64 faceCenter[ 3 ];
    real64 normal[ 3 ];
    real64 area = computationalGeometry::centroid_3DPolygon( faces[f], points, faceCenter, normal );
    real64 signTestVector[3];
    LvArray::tensorOps::copy< 3 >( signTestVector, faceCenter );
    LvArray::tensorOps::subtract< 3 >( signTestVector, centroid );
    real64 const dot = LvArray::tensorOps::AiBi< 3 >( signTestVector, normal );
    if( dot < 0 )
    {
      LvArray::tensorOps::scale< 3 >( normal, -1.0 );
    }
    std::cout << "normal: " << area*normal[0] << " " << area*normal[1] << " " << area*normal[2] << std::endl ;
    LvArray::tensorOps::scaledAdd< 3 >( divergence, normal, area );
    std::cout << "div: " << divergence[0] << " " << divergence[1] << " " << divergence[2] << std::endl;
  }
  real64 const norm = LvArray::tensorOps::l2Norm< 3 >(divergence);
  EXPECT_TRUE( norm < 10*LvArray::NumericLimits< real64 >::epsilon ) << norm;
}

TEST( testComputationalGeometry, checkNormalsConsistency_NEW )
{
  real64 pArray[ 8 ][ 3 ] = {
    { 0.1084442213177681, 0.0, 0.0                },
    { 0.2075795382261276, 0.0, 0.0                },
    { 0.1058338209986687, 0.1052169278264046, 0.0                },
    { 0.2090811282396317, 0.100415863096714, 0.0                },
    { 0.1048644408583641, 0.0, 0.1068341732025146 },
    { 0.2063014656305313, 0.0, 0.1084721013903618 },
    { 0.1081591323018074, 0.1056273579597473, 0.1087869703769684 },
    { 0.2010060846805573, 0.1080376580357552, 0.1048369631171227 }
  };
  array2d< real64, nodes::REFERENCE_POSITION_PERM > points( 8, 3 );
  LvArray::tensorOps::copy< 8, 3 >( points, pArray );
  real64 centroid[3];
  LvArray::tensorOps::fill< 3 >( centroid, 0 );
  for( localIndex i = 0; i < 8; ++i )
  {
    LvArray::tensorOps::add< 3 >( centroid, points[ i ] );
  }
  LvArray::tensorOps::scale< 3 >( centroid, 1.0/8.0 );

  localIndex fArray[6][4] = {
    { 0, 2, 3, 1 },
    { 6, 4, 5, 7 },
    { 6, 7, 3, 2 },
    { 5, 4, 0, 1 },
    { 7, 5, 1, 3 },
    { 2, 0, 4, 6 }
  };
  array2d< localIndex > faces( 6, 4 );
  LvArray::tensorOps::copy< 6, 4 >( faces, fArray );

  real64 divergence[ 3 ]; // numerical divergences of (1,0,0), (0,1,0), (0,0,1)
  LvArray::tensorOps::fill< 3 >( divergence, 0.0 );
  for( localIndex f = 0; f < 6; ++f )
  {
    real64 faceCenter[ 3 ];
    real64 normal[ 3 ];
    real64 area = computationalGeometry::centroid_3DPolygon_NEW( faces[f], points, faceCenter, normal );
    real64 signTestVector[3];
    LvArray::tensorOps::copy< 3 >( signTestVector, faceCenter );
    LvArray::tensorOps::subtract< 3 >( signTestVector, centroid );
    real64 const dot = LvArray::tensorOps::AiBi< 3 >( signTestVector, normal );
    if( dot < 0 )
    {
      LvArray::tensorOps::scale< 3 >( normal, -1.0 );
    }
    std::cout << "normal: " << area*normal[0] << " " << area*normal[1] << " " << area*normal[2] << std::endl ;
    LvArray::tensorOps::scaledAdd< 3 >( divergence, normal, area );
    std::cout << "div: " << divergence[0] << " " << divergence[1] << " " << divergence[2] << std::endl;
  }
  std::cout << std::endl;
  real64 const norm = LvArray::tensorOps::l2Norm< 3 >(divergence);
  EXPECT_TRUE( LvArray::tensorOps::l2Norm< 3 >(divergence) < 10*LvArray::NumericLimits< real64 >::epsilon ) << divergence[0] << " " << divergence[1] << " " << divergence[2];
}

} /* namespace geosx */
