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
 * @file testComputationalGeometry.cpp
 */
#include "../utilities/ComputationalGeometry.hpp"
#include <cmath>
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

TEST( testComputationalGeometry, checkHighAspectRatio)
{
  constexpr localIndex numNodes = 3;

  array2d< real64, nodes::REFERENCE_POSITION_PERM > points;
  points.resize(numNodes, 3);
  array1d< localIndex > indices(3);
  indices(0) = 0;
  indices(1) = 1;
  indices(2) = 2;

  //mimic a perfect triangle but with huge offset
  constexpr real64 xa = 1.e4, ya = 1.e-2, za = 1.e0;
  constexpr real64 xb = 1.e-1, yb = 1.e2, zb = -1.e-2;
  constexpr real64 TOL = 1e-6;

  points(0,0) = 1.e6;
  points(0,1) = 1.e6;
  points(0,2) = 1.e6;
  
  points(1,0) = 1.e6 + xa;
  points(1,1) = 1.e6 + ya;
  points(1,2) = 1.e6 + za;

  points(2,0) = 1.e6 + xb;
  points(2,1) = 1.e6 + yb;
  points(2,2) = 1.e6 + zb;

  real64 faceCenter[ 3 ], faceNormal[ 3 ];
  real64 faceCenterOld[ 3 ], faceNormalOld[ 3 ];
  auto faceArea = computationalGeometry::centroid_3DPolygon(indices.toSliceConst(), points.toViewConst(), faceCenter, faceNormal, TOL);

  constexpr real64 EXPECTED_AREA = 0.5*std::sqrt(LvArray::math::square(xa*yb-ya*xb) + 
                                                    LvArray::math::square(ya*zb-za*yb) +
                                                    LvArray::math::square(za*xb-xa*zb) );

  constexpr real64 EXPECTED_CENTER[3] = { 1.e6 + (xa+xb)/3 , 1.e6 + (ya+yb)/3. , 1.e6 + (za+zb)/3};

  constexpr real64 EXPECTED_NORMAL[3] = { (ya*zb-za*yb)/2/EXPECTED_AREA, (za*xb-xa*zb)/2/EXPECTED_AREA, (xa*yb-ya*xb)/2/EXPECTED_AREA};

  EXPECT_LT( abs(faceArea - EXPECTED_AREA), 1e-6); 
  EXPECT_LT( abs(faceNormal[0]-EXPECTED_NORMAL[0]), 1e-6);
  EXPECT_LT( abs(faceNormal[1]-EXPECTED_NORMAL[1]), 1e-6);
  EXPECT_LT( abs(faceNormal[2]-EXPECTED_NORMAL[2]), 1e-6);
  EXPECT_LT( abs(faceCenter[0]-EXPECTED_CENTER[0]), 1e-6);
  EXPECT_LT( abs(faceCenter[0]-EXPECTED_CENTER[0]), 1e-6);
  EXPECT_LT( abs(faceCenter[1]-EXPECTED_CENTER[1]), 1e-6);
  EXPECT_LT( abs(faceCenter[2]-EXPECTED_CENTER[2]), 1e-6);

  real64 norm = LvArray::math::square( faceNormal[0] ) + LvArray::math::square( faceNormal[1] ) + LvArray::math::square( faceNormal[2] );
  EXPECT_LT(abs( sqrt( norm ) - 1.0 ), 1e-6);

}

} /* namespace geos */
