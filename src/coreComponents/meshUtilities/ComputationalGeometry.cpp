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
 * @file ComputationalGeometry.cpp
 */

#include "ComputationalGeometry.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{


namespace computationalGeometry
{

constexpr real64 machinePrecision = std::numeric_limits< real64 >::epsilon();


//*************************************************************************************************
void FixNormalOrientation_3D( arraySlice1d< real64 > const normal )
{
  real64 const orientationTolerance = 1.e+1*machinePrecision;

  // Orient local normal in global sense.
  // First check: align with z direction
  if( normal[ 2 ] <= -orientationTolerance )
  {
    LvArray::tensorOps::scale< 3 >( normal, -1.0 );
  }
  else if( std::fabs( normal[ 2 ] ) < orientationTolerance )
  {
    // If needed, second check: align with y direction
    if( normal[ 1 ] <= -orientationTolerance )
    {
      LvArray::tensorOps::scale< 3 >( normal, -1.0 );
    }
    else if( std::fabs( normal[ 1 ] ) < orientationTolerance )
    {
      // If needed, third check: align with x direction
      if( normal[ 0 ] <= -orientationTolerance )
      {
        LvArray::tensorOps::scale< 3 >( normal, -1.0 );
      }
    }
  }
}

//*************************************************************************************************
void RotationMatrix_3D( arraySlice1d< real64 const > const normal,
                        arraySlice2d< real64 > const rotationMatrix )
{
  real64 m1[ 3 ] = { normal[ 2 ], 0.0, -normal[ 0 ] };
  real64 m2[ 3 ] = { 0.0, normal[ 2 ], -normal[ 1 ] };
  real64 const norm_m1 = LvArray::tensorOps::l2Norm< 3 >( m1 );
  real64 const norm_m2 = LvArray::tensorOps::l2Norm< 3 >( m2 );

  // If present, looks for a vector with 0 norm
  // Fix the uncertain case of norm_m1 very close to norm_m2
  if( norm_m1+1.e+2*machinePrecision > norm_m2 )
  {
    LvArray::tensorOps::crossProduct( m2, normal, m1 );
    LvArray::tensorOps::normalize< 3 >( m2 );
    LvArray::tensorOps::normalize< 3 >( m1 );
  }
  else
  {
    LvArray::tensorOps::crossProduct( m1, normal, m2 );
    LvArray::tensorOps::scale< 3 >( m1, -1 );
    LvArray::tensorOps::normalize< 3 >( m1 );
    LvArray::tensorOps::normalize< 3 >( m2 );
  }

  // Save everything in the standard form (3x3 rotation matrix)
  rotationMatrix( 0, 0 ) = normal[ 0 ];
  rotationMatrix( 1, 0 ) = normal[ 1 ];
  rotationMatrix( 2, 0 ) = normal[ 2 ];
  rotationMatrix( 0, 1 ) = m1[ 0 ];
  rotationMatrix( 1, 1 ) = m1[ 1 ];
  rotationMatrix( 2, 1 ) = m1[ 2 ];
  rotationMatrix( 0, 2 ) = m2[ 0 ];
  rotationMatrix( 1, 2 ) = m2[ 1 ];
  rotationMatrix( 2, 2 ) = m2[ 2 ];

  GEOSX_ERROR_IF( std::fabs( LvArray::tensorOps::determinant< 3 >( rotationMatrix ) - 1.0 ) > 1.e+1*machinePrecision,
                  "Rotation matrix with determinant different from +1.0" );

  return;
}

//*************************************************************************************************
template< typename T >
int sgn( T val )
{
  return (T( 0 ) < val) - (val < T( 0 ));
}

//*************************************************************************************************
bool IsPointInsidePolyhedron( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodeCoordinates,
                              array1d< array1d< localIndex > > const & faceNodeIndicies,
                              real64 const ( &point )[3],
                              real64 const areaTolerance )
{
  localIndex const numFaces = faceNodeIndicies.size( 0 );
  real64 faceCenter[ 3 ], faceNormal[ 3 ];
  int sign = 0;

  for( localIndex kf = 0; kf < numFaces; ++kf )
  {
    Centroid_3DPolygon( faceNodeIndicies[kf], nodeCoordinates, faceCenter, faceNormal, areaTolerance );

    LvArray::tensorOps::subtract< 3 >( faceCenter, point );
    int const s = sgn( LvArray::tensorOps::AiBi< 3 >( faceNormal, faceCenter ) );

    // all dot products should be non-negative (for outward normals) or non-positive (for inward normals)
    if( sign * s < 0 )
    {
      return false;
    }
    sign = s;
  }

  return true;
}

//*************************************************************************************************
GEOSX_HOST_DEVICE
real64 TetVolume( real64 const X[][3] )
{
  real64 X1_X0[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( X[1] );
  LvArray::tensorOps::subtract< 3 >( X1_X0, X[0] );

  real64 X2_X0[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( X[2] );
  LvArray::tensorOps::subtract< 3 >( X2_X0, X[0] );

  real64 X3_X0[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( X[3] );
  LvArray::tensorOps::subtract< 3 >( X3_X0, X[0] );

  real64 X2_X0crossX3_X0[ 3 ];
  LvArray::tensorOps::crossProduct( X2_X0crossX3_X0, X2_X0, X3_X0 );

  return std::fabs( LvArray::tensorOps::AiBi< 3 >( X1_X0, X2_X0crossX3_X0 ) / 6.0 );
}

//*************************************************************************************************
GEOSX_HOST_DEVICE
real64 WedgeVolume( real64 const X[][3] )
{
  real64 const tet1[4][3] = { LVARRAY_TENSOROPS_INIT_LOCAL_3( X[0] ),
                              LVARRAY_TENSOROPS_INIT_LOCAL_3( X[1] ),
                              LVARRAY_TENSOROPS_INIT_LOCAL_3( X[2] ),
                              LVARRAY_TENSOROPS_INIT_LOCAL_3( X[4] ) };

  real64 const tet2[4][3] = { LVARRAY_TENSOROPS_INIT_LOCAL_3( X[0] ),
                              LVARRAY_TENSOROPS_INIT_LOCAL_3( X[2] ),
                              LVARRAY_TENSOROPS_INIT_LOCAL_3( X[4] ),
                              LVARRAY_TENSOROPS_INIT_LOCAL_3( X[5] ) };

  real64 const tet3[4][3] = { LVARRAY_TENSOROPS_INIT_LOCAL_3( X[0] ),
                              LVARRAY_TENSOROPS_INIT_LOCAL_3( X[3] ),
                              LVARRAY_TENSOROPS_INIT_LOCAL_3( X[4] ),
                              LVARRAY_TENSOROPS_INIT_LOCAL_3( X[5] ) };

  return TetVolume( tet1 ) + TetVolume( tet2 ) + TetVolume( tet3 );
}

//*************************************************************************************************
GEOSX_HOST_DEVICE
real64 PyramidVolume( real64 const X[][3] )
{
  real64 const tet1[4][3] = { LVARRAY_TENSOROPS_INIT_LOCAL_3( X[0] ),
                              LVARRAY_TENSOROPS_INIT_LOCAL_3( X[1] ),
                              LVARRAY_TENSOROPS_INIT_LOCAL_3( X[2] ),
                              LVARRAY_TENSOROPS_INIT_LOCAL_3( X[4] ) };

  real64 const tet2[4][3] = { LVARRAY_TENSOROPS_INIT_LOCAL_3( X[0] ),
                              LVARRAY_TENSOROPS_INIT_LOCAL_3( X[2] ),
                              LVARRAY_TENSOROPS_INIT_LOCAL_3( X[3] ),
                              LVARRAY_TENSOROPS_INIT_LOCAL_3( X[4] ) };

  return TetVolume( tet1 ) + TetVolume( tet2 );
}

//*************************************************************************************************
void GetBoundingBox( localIndex elemIndex,
                     arrayView2d< localIndex const, cells::NODE_MAP_USD > const & pointIndices,
                     arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & pointCoordinates,
                     real64 ( & boxDims )[ 3 ] )
{
  localIndex constexpr dim = 3;

  // This holds the min coordinates of the set in each direction
  real64 minCoords[ 3 ] = { std::numeric_limits< double >::max(), std::numeric_limits< double >::max(), std::numeric_limits< double >::max() };

  // boxDims is used to hold the max coordinates.
  LvArray::tensorOps::fill< 3 >( boxDims, std::numeric_limits< double >::min() );

  // loop over all the vertices of the element to get the min and max coords
  for( localIndex a = 0; a < pointIndices.size( 1 ); ++a )
  {
    localIndex const id = pointIndices( elemIndex, a );
    for( localIndex d = 0; d < dim; ++d )
    {
      minCoords[ d ] = std::min( minCoords[ d ], pointCoordinates( id, d ) );
      boxDims[ d ] = std::max( boxDims[ d ], pointCoordinates( id, d ) );
    }
  }

  LvArray::tensorOps::subtract< 3 >( boxDims, minCoords );
}

}

} /* namespace geosx */
