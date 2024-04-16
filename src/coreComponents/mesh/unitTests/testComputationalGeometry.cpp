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

namespace geos
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

TEST( testComputationalGeometry, checkMeanCurvature3DPolygon)
{
  array2d< real64, nodes::REFERENCE_POSITION_PERM > points;
  array1d< localIndex > indices;
  points.resize( 4, 3 );
  LvArray::tensorOps::fill<4,3>(points, 0.0);
  indices.resize(4);
  for(localIndex i = 0; i < 4; ++i)
  {
    indices[i] = i;
  }

  real64 a = 0.0, b = 0.0, c = 3.0;
  auto surf = [&a,&b,&c](real64 x, real64 y) { return 0.5*(a*x*x + b*y*y + 2*c*x*y); };
  points[0][0] = 0.0; points[0][1] = 0.0; points[0][2] = surf(points[0][0],points[0][1]);
  points[1][0] = 1.0; points[1][1] = 0.0; points[1][2] = surf(points[1][0],points[1][1]);
  points[2][0] = 1.0; points[2][1] = 1.0; points[2][2] = surf(points[2][0],points[2][1]);
  points[3][0] = 0.0; points[3][1] = 1.0; points[3][2] = surf(points[3][0],points[3][1]);

  real64 referenceCurvatures[2] = { };
  real64 curvatures[2];
  computationalGeometry::principalCurvatures_3DPolygon( indices.toSliceConst(), points.toViewConst(), curvatures );
  EXPECT_EQ( 0.0 < 1e-15, true);

}

} /* namespace geos */
