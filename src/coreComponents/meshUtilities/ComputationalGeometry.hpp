/*
 * ComputationalGeometry.hpp
 *
 *  Created on: Jun 26, 2018
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MESHUTILITIES_COMPUTATIONALGEOMETRY_HPP_
#define SRC_COMPONENTS_CORE_SRC_MESHUTILITIES_COMPUTATIONALGEOMETRY_HPP_

#include "common/DataTypes.hpp"

namespace geosx
{
namespace computationalGeometry
{

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
real64 Centroid_3DPolygon(const localIndex_array& pointsIndices,
                         const array1d<R1Tensor>& points,
                         R1Tensor& center,
                         R1Tensor& normal );

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
real64 Centroid_3DPolygon(const localIndex_array& pointsIndices,
                         const array1d<R1Tensor>& pointReferences,
                         const array1d<R1Tensor>& pointDisplacements,
                         R1Tensor& center,
                         R1Tensor& normal );


real64 HexVolume( R1Tensor const * const points );

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
} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MESHUTILITIES_COMPUTATIONALGEOMETRY_HPP_ */
