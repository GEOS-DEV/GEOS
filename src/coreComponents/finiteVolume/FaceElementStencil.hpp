/*
 * FaceElementStencil.hpp
 *
 *  Created on: Jul 1, 2019
 *      Author: settgast
 */

#ifndef SRC_CORECOMPONENTS_FINITEVOLUME_FACEELEMENTSTENCIL_HPP_
#define SRC_CORECOMPONENTS_FINITEVOLUME_FACEELEMENTSTENCIL_HPP_

#include "common/DataTypes.hpp"

namespace geosx
{

class FaceElementStencil
{
public:
  using INDEX_TYPE = ArrayOfArrays<localIndex>;
  using WEIGHT_TYPE = ArrayOfArrays<real64>;

  using INDEX_VIEW_TYPE = ArrayOfArraysView<localIndex>;
  using INDEX_VIEW_CONST_TYPE = ArrayOfArraysView<localIndex const>;

  using WEIGHT_VIEW_TYPE = ArrayOfArraysView<real64>;
  using WEIGHT_VIEW_CONST_TYPE = ArrayOfArraysView<real64 const>;

  /**
   * @brief Number of points the flux is between (normally 2)
   */
  static localIndex constexpr NUM_POINT_IN_FLUX = 6;

  /**
   * @brief Maximum number of points in a stencil (required to use static arrays in kernels)
   */
  static localIndex constexpr MAX_STENCIL_SIZE = 6;



  FaceElementStencil();
  ~FaceElementStencil();


  void add( localIndex const numPts,
            localIndex  const * const elementRegionIndices,
            localIndex  const * const elementSubRegionIndices,
            localIndex  const * const elementIndices,
            real64 const * const weights,
            localIndex const connectorIndex );

  bool zero( localIndex const connectorIndex );

  INDEX_VIEW_CONST_TYPE const &  getElementRegionIndices() const { return m_elementRegionIndices; }
  INDEX_VIEW_CONST_TYPE const &  getElementSubRegionIndices() const { return m_elementSubRegionIndices; }
  INDEX_VIEW_CONST_TYPE const &  getElementIndices() const { return m_elementIndices; }
  WEIGHT_VIEW_CONST_TYPE const & getWeights() const { return m_weights; }

private:
  INDEX_TYPE  m_elementRegionIndices;
  INDEX_TYPE  m_elementSubRegionIndices;
  INDEX_TYPE  m_elementIndices;
  WEIGHT_TYPE m_weights;
  map<localIndex, localIndex> m_connectorIndices;

};

} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_FINITEVOLUME_FACEELEMENTSTENCIL_HPP_ */
