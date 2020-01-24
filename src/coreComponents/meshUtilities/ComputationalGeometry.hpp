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
#include "mesh/InterObjectRelation.hpp"

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
                           arrayView1d<R1Tensor const> const & points,
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
                           arrayView1d<R1Tensor const> const & points,
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
bool IsPointInsidePolyhedron( arrayView1d<R1Tensor const> const & nodeCoordinates,
                              array1d<array1d<localIndex>> const & faceNodeIndicies,
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
template<typename NODEMAP>
R1Tensor GetBoundingBox( localIndex elemIndex,
                         NODEMAP const & pointIndices,
                         arrayView1d<R1Tensor const> const & pointCoordinates );

real64 HexVolume( R1Tensor const * const points );

real64 TetVolume( R1Tensor const * const points );

real64 WedgeVolume( R1Tensor const * const points );

real64 PyramidVolume( R1Tensor const * const points );
  
inline void VectorDifference( array1d< R1Tensor > const & X,
                              localIndex const index0,
                              localIndex const index1,
                              R1Tensor & vec )
{
  vec = X[index1];
  vec -= X[index0];
}

template< int N >
inline void VectorMean( array1d< R1Tensor > const & X,
                        arrayView1d<localIndex> const indices,
                        R1Tensor & vec )
{
  vec = 0;
  for( int a=0 ; a<N ; ++a )
  {
    vec += X[indices[a]];
  }
  vec /= N;
}


}


template<typename NODEMAP>
R1Tensor computationalGeometry::GetBoundingBox( localIndex elemIndex,
                                                NODEMAP const & pointIndices,
                                                arrayView1d<R1Tensor const> const & pointCoordinates )
{
  localIndex constexpr dim = 3;
  
  // these arrays will store the min and max coordinates of the elem in each direction
  R1Tensor minCoords(  1e99 );
  R1Tensor maxCoords( -1e99 );

  // loop over all the vertices of the element to get the min and max coords
  for (localIndex a = 0; a < pointIndices.size( 1 ); ++a)
  {
    localIndex const id = pointIndices( elemIndex, a );
    R1Tensor const coords = pointCoordinates[id];

    for (localIndex d = 0; d < dim; ++d)
    {  
      if (coords[d] < minCoords[d])
      {
        minCoords[d] = coords[d];
      }
      else if (coords[d] > maxCoords[d])
      {
        maxCoords[d] = coords[d];
      }
    }
  }

  // compute the dimensions of the bounding box
  R1Tensor box( 0 );
  for (localIndex d = 0; d < dim; ++d)
  {
    box[d] = maxCoords[d] - minCoords[d];
  }
 
  return box;
}


} /* namespace geosx */

#endif /* GEOSX_MESHUTILITIES_COMPUTATIONALGEOMETRY_HPP_ */
