/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ComputationalGeometry.hpp
 */

#ifndef GEOSX_MESHUTILITIES_COMPUTATIONALGEOMETRY_HPP_
#define GEOSX_MESHUTILITIES_COMPUTATIONALGEOMETRY_HPP_

#include "common/DataTypes.hpp"
#include "common/DataLayouts.hpp"

namespace geosx
{
namespace computationalGeometry
{

/**
 * Calculates the centroid of a convex 3D polygon as well as the normal
 * @param[in] pointIndices list of index references for the points array in
 * order (CW or CCW) about the polygon loop
 * @param [in] numPoints the number of points in the polygon
 * @param[in] points 3D point list
 * @param[out] center 3D center of the given ordered polygon point list
 * @param[out] normal Normal to the face
 * @return area of the convex 3D polygon
 */
real64 Centroid_3DPolygon( localIndex const * const pointsIndices,
                           localIndex const numPoints,
                           arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & points,
                           R1Tensor & center,
                           R1Tensor & normal,
                           real64 areaTolerance = 0.0 );

/**
 * @author settgast
 * Calculates the centroid of a convex 3D polygon as well as the normal
 * @param[in] pointIndices list of index references for the points array in
 * order (CW or CCW) about the polygon loop
 * @param[in] points 3D point list
 * @param[out] center 3D center of the given ordered polygon point list
 * @param[out] normal Normal to the face
 * @return area of the convex 3D polygon
 */
real64 Centroid_3DPolygon( arrayView1d<localIndex const> const & pointsIndices,
                           arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & points,
                           R1Tensor & center,
                           R1Tensor & normal,
                           real64 areaTolerance = 0.0 );

/**
 * Calculates the centroid of a convex 3D polygon as well as the normal
 * @param[in] pointIndices list of index references for the points array in
 * order (CW or CCW) about the polygon loop
 * @param[in] numPoints the numberof points in the polygon
 * @param[in] pointReferences 3D reference point list
 * @param[in] pointDisplacements 3D displacement list
 * @param[out] center 3D center of the given ordered polygon point list
 * @param[out] normal Normal to the face
 * @return area of the convex 3D polygon
 */
real64 Centroid_3DPolygon( localIndex const * const pointsIndices,
                           localIndex const numPoints,
                           arrayView1d<R1Tensor const> const & pointReferences,
                           arrayView1d<R1Tensor const> const & pointDisplacements,
                           R1Tensor & center,
                           R1Tensor & normal );

/**
 * @author settgast
 * Calculates the centroid of a convex 3D polygon as well as the normal
 * @param[in] pointIndices list of index references for the points array in
 * order (CW or CCW) about the polygon loop
 * @param[in] pointReferences 3D reference point list
 * @param[in] pointDisplacements 3D displacement list
 * @param[out] center 3D center of the given ordered polygon point list
 * @param[out] normal Normal to the face
 * @return area of the convex 3D polygon
 */
real64 Centroid_3DPolygon( arrayView1d<localIndex const> const & pointsIndices,
                           arrayView1d<R1Tensor const> const & pointReferences,
                           arrayView1d<R1Tensor const> const & pointDisplacements,
                           R1Tensor & center,
                           R1Tensor & normal );

/**
 * @brief Check if a point is inside a convex polyhedron (3D polygon)
 * @param[in] nodeCoordinates a global array of nodal coordinates
 * @param[in] faceNodeIndicies ordered lists of node indices for each face of the polyhedron
 * @param[in] point coordinates of the query point
 * @param[in] areaTolerance same as in Centroid_3DPolygon
 * @return whether the point is inside
 *
 * @note Face nodes must all be ordered the same way (i.e. CW or CCW),
 * resulting in all face normals pointing either outside or inside the polyhendron
 *
 * @note For faces with n>3 nodes that are non-planar, average normal is used
 */
bool IsPointInsidePolyhedron( arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & nodeCoordinates,
                              array1d<array1d<localIndex>> const & faceNodeIndicies,
                              R1Tensor const & point,
                              real64 const areaTolerance = 0.0 );

inline
GEOSX_HOST_DEVICE
real64 HexVolume( R1Tensor const * const X )
{
  R1Tensor X7_X1( X[7] );
  X7_X1 -= X[1];

  R1Tensor X6_X0( X[6] );
  X6_X0 -= X[0];

  R1Tensor X7_X2( X[7] );
  X7_X2 -= X[2];

  R1Tensor X3_X0( X[3] );
  X3_X0 -= X[0];

  R1Tensor X5_X0( X[5] );
  X5_X0 -= X[0];

  R1Tensor X7_X4( X[7] );
  X7_X4 -= X[4];

  R1Tensor X7_X1plusX6_X0( X7_X1 );
  X7_X1plusX6_X0 += X6_X0;

  R1Tensor X7_X2plusX5_X0( X7_X2 );
  X7_X2plusX5_X0 += X5_X0;

  R1Tensor X7_X4plusX3_X0( X7_X4 );
  X7_X4plusX3_X0 += X3_X0;

  return 1.0/12.0 * ( Dot( X7_X1plusX6_X0, Cross( X7_X2, X3_X0 ) ) +
                      Dot( X6_X0, Cross( X7_X2plusX5_X0, X7_X4 ) ) +
                      Dot( X7_X1, Cross( X5_X0, X7_X4plusX3_X0 ) ) );
}

real64 TetVolume( R1Tensor const * const points );

real64 WedgeVolume( R1Tensor const * const points );

real64 PyramidVolume( R1Tensor const * const points );

}
} /* namespace geosx */

#endif /* GEOSX_MESHUTILITIES_COMPUTATIONALGEOMETRY_HPP_ */
