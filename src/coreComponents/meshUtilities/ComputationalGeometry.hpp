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
#include "mesh/InterObjectRelation.hpp"

namespace geosx
{
namespace computationalGeometry
{


/**
 * @brief Calculate the intersection between a line and a plane.
 * @param[in] lineDir vector defining direction of the line
 * @param[in] linePoint one point of the line
 * @param[in] planeNormal normal to plane
 * @param[in] planeOrigin plane origin
 * @return the area of the convex 3D polygon
 */
R1Tensor LinePlaneIntersection( R1Tensor lineDir,
                                R1Tensor linePoint,
                                R1Tensor planeNormal,
                                R1Tensor planeOrigin );

/**
 * @brief Calculate the area of a polygon given the set of points defining it.
 * @param[in] points coordinates of the points
 * @param[in] numPoints number of points
 * @param[in] normal unit normal vector to the surface
 * @return the area of the polygon
 */
real64 ComputeSurfaceArea( array1d< R1Tensor > const & points,
                           localIndex const numPoints,
                           R1Tensor const & normal );

/**
 * @brief Order a set of points counter-clockwise.
 * @param[in] points coordinates of the points
 * @param[in] numPoints number of points
 * @param[in] normal unit normal vector to the surface
 * @return the reordered set of points
 */
array1d< R1Tensor > orderPointsCCW( array1d< R1Tensor > const & points,
                                    localIndex const numPoints,
                                    R1Tensor const & normal );
/**
 * @brief Calculate the centroid of a convex 3D polygon as well as the normal.
 * @param[in] pointsIndices list of index references for the points array in
 * order (CW or CCW) about the polygon loop
 * @param[in] numPoints the number of points in the polygon
 * @param[in] points 3D point list
 * @param[out] center 3D center of the given ordered polygon point list
 * @param[out] normal normal to the face
 * @param[in] areaTolerance tolerance used in the geometric computations
 * @return area of the convex 3D polygon
 * @details if area < - areaTolerance, this function will throw an error,
 *          and if (- areaTolerance <= area <= areaTolerance), the area is set to zero
 */
real64 Centroid_3DPolygon( localIndex const * const pointsIndices,
                           localIndex const numPoints,
                           arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & points,
                           R1Tensor & center,
                           R1Tensor & normal,
                           real64 areaTolerance = 0.0 );

/**
 * @brief Calculate the centroid of a convex 3D polygon as well as the normal and the rotation matrix.
 * @param[in] pointsIndices list of index references for the points array in
 * order (CW or CCW) about the polygon loop
 * @param[in] numPoints the number of points in the polygon
 * @param[in] points 3D point list
 * @param[out] center 3D center of the given ordered polygon point list
 * @param[out] normal normal to the face
 * @param[out] rotationMatrix rotation matrix for the face
 * @param[in] areaTolerance tolerance used in the geometric computations
 * @return area of the convex 3D polygon
 * @details if area < - areaTolerance, this function will throw an error,
 *          and if (- areaTolerance <= area <= areaTolerance), the area is set to zero
 */
real64 Centroid_3DPolygon( localIndex const * const pointsIndices,
                           localIndex const numPoints,
                           arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & points,
                           R1Tensor & center,
                           R1Tensor & normal,
                           R2Tensor & rotationMatrix,
                           real64 const areaTolerance = 0.0 );

/**
 * @brief Change the orientation of the input vector to be consistent in a global sense.
 * @param[inout] normal normal to the face
 */
void FixNormalOrientation_3D( R1Tensor & normal );

/**
 * @brief Calculate the rotation matrix for a face in the 3D space
 * @param[in] normal normal to the face
 * @param[out] rotationMatrix rotation matrix for the face
 */
void RotationMatrix_3D( R1Tensor const & normal,
                        R2Tensor & rotationMatrix );

/**
 * @brief Calculate the centroid of a convex 3D polygon as well as the normal.
 * @param[in] pointsIndices list of index references for the points array in
 * order (CW or CCW) about the polygon loop
 * @param[in] points 3D point list
 * @param[out] center 3D center of the given ordered polygon point list
 * @param[out] normal normal to the face
 * @param[in] areaTolerance tolerance used in the geometric computations
 * @return area of the convex 3D polygon
 */
real64 Centroid_3DPolygon( arrayView1d< localIndex const > const & pointsIndices,
                           arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & points,
                           R1Tensor & center,
                           R1Tensor & normal,
                           real64 areaTolerance = 0.0 );

/**
 * @brief Calculate the centroid of a convex 3D polygon as well as the normal.
 * @param[in] pointsIndices list of index references for the points array in
 * order (CW or CCW) about the polygon loop
 * @param[in] numPoints the numberof points in the polygon
 * @param[in] pointReferences 3D reference point list
 * @param[in] pointDisplacements 3D displacement list
 * @param[out] center 3D center of the given ordered polygon point list
 * @param[out] normal normal to the face
 * @return area of the convex 3D polygon
 */
real64 Centroid_3DPolygon( localIndex const * const pointsIndices,
                           localIndex const numPoints,
                           arrayView1d< R1Tensor const > const & pointReferences,
                           arrayView1d< R1Tensor const > const & pointDisplacements,
                           R1Tensor & center,
                           R1Tensor & normal );

/**
 * @brief Calculate the centroid of a convex 3D polygon as well as the normal.
 * @param[in] pointsIndices list of index references for the points array in
 * order (CW or CCW) about the polygon loop
 * @param[in] pointReferences 3D reference point list
 * @param[in] pointDisplacements 3D displacement list
 * @param[out] center 3D center of the given ordered polygon point list
 * @param[out] normal normal to the face
 * @return area of the convex 3D polygon
 */
real64 Centroid_3DPolygon( arrayView1d< localIndex const > const & pointsIndices,
                           arrayView1d< R1Tensor const > const & pointReferences,
                           arrayView1d< R1Tensor const > const & pointDisplacements,
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
bool IsPointInsidePolyhedron( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodeCoordinates,
                              array1d< array1d< localIndex > > const & faceNodeIndicies,
                              R1Tensor const & point,
                              real64 const areaTolerance = 0.0 );

/**
 * @brief Compute the dimensions of the bounding box containing the element
          defined here by the coordinates of its vertices
 * @param[in] elemIndex index of the element in pointIndices
 * @param[in] pointIndices the indices of the vertices in pointCoordinates
 * @param[in] pointCoordinates the vertices coordinates
 * @return an R1Tensor containing the dimensions of the box
 */
template< typename NODEMAP >
R1Tensor GetBoundingBox( localIndex elemIndex,
                         NODEMAP const & pointIndices,
                         arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & pointCoordinates );

/**
 * @brief Compute the volume of an hexahedron
 * @param[in] X vertices of the hexahedron
 * @return the volume of the hexahedron
 */
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

/**
 * @brief Compute the volume of an tetrahedron
 * @param[in] points vertices of the tetrahedron
 * @return the volume of the tetrahedron
 */
real64 TetVolume( R1Tensor const * const points );

/**
 * @brief Compute the volume of a wedge
 * @param[in] points vertices of the wedge
 * @return the volume of the wedge
 */
real64 WedgeVolume( R1Tensor const * const points );

/**
 * @brief Compute the volume of a pyramid
 * @param[in] points vertices of the pyramid
 * @return the volume of the pyramid
 */
real64 PyramidVolume( R1Tensor const * const points );

}

/// @cond DO_NOT_DOCUMENT

extern template R1Tensor computationalGeometry::GetBoundingBox( localIndex elemIndex,
                                                                InterObjectRelation< array2d< localIndex, RAJA::PERM_IJ > > const & pointIndices,
                                                                arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & pointCoordinates );
extern template R1Tensor computationalGeometry::GetBoundingBox( localIndex elemIndex,
                                                                InterObjectRelation< array2d< localIndex, RAJA::PERM_JI > > const & pointIndices,
                                                                arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & pointCoordinates );

/// @endcond

} /* namespace geosx */

#endif /* GEOSX_MESHUTILITIES_COMPUTATIONALGEOMETRY_HPP_ */
