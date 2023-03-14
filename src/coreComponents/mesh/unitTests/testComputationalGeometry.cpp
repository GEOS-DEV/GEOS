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
  constexpr int n = 12;
  real64 t = 2.0*M_PI/n;
  real64 s = (1.0 + sin( t ))/cos( t ) - 1.0;
  constexpr int nodeDim = 3;

  real64 e[3] = { -1.0e-6, 0.0, 1.0e-6 };

  array2d< real64 > points( n*3, nodeDim );

  int idx = 0;
  for( int z = 0; z < 3; ++z )
  {
    for( int i = 0; i < n; ++i )
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

  std::cout << idx << std::endl;

  constexpr real64 areaTolerance = 0.0001;

  for( int ci = 0; ci < 3; ++ci )
  {
    array1d< int > indices;
    for( int i = 0; i < n; ++i )
    {
      indices.emplace_back( ci*n + i );
    }

    GEOSX_LOG_VAR( indices );

    real64 faceCenter[ 3 ], faceNormal[ 3 ];
    real64 const faceArea = computationalGeometry::centroid_3DPolygon( indices, points, faceCenter, faceNormal, areaTolerance );

    real64 norm = LvArray::math::square( faceNormal[0] ) + LvArray::math::square( faceNormal[1] ) + LvArray::math::square( faceNormal[2] );
    GEOSX_LOG( GEOSX_FMT( "Cell: {} Normal: [{}, {}, {}], Area: {}", ci, faceNormal[0], faceNormal[1], faceNormal[2], faceArea ));
    std::cout << norm << " " << abs( sqrt( norm ) - 1.0 ) << std::endl;
    EXPECT_EQ( abs( sqrt( norm ) - 1.0 ) < 0.0001, true );
  }
}

} /* namespace geosx */
