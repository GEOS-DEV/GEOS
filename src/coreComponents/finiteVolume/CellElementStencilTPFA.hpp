/*
 * CellElementStencil.hpp
 *
 *  Created on: Jul 10, 2019
 *      Author: settgast
 */

#ifndef SRC_CORECOMPONENTS_FINITEVOLUME_CELLELEMENTSTENCILTPFA_HPP_
#define SRC_CORECOMPONENTS_FINITEVOLUME_CELLELEMENTSTENCILTPFA_HPP_

#include "common/DataTypes.hpp"

namespace geosx
{

class CellElementStencilTPFA
{
public:

  using INDEX_TYPE = array2d<localIndex>;
  using WEIGHT_TYPE = array2d<real64>;

  using INDEX_VIEW_TYPE = arrayView2d<localIndex>;
  using INDEX_VIEW_CONST_TYPE = arrayView2d<localIndex const>;

  using WEIGHT_VIEW_TYPE = arrayView2d<real64>;
  using WEIGHT_VIEW_CONST_TYPE = arrayView2d<real64 const>;

  /**
   * @brief Number of points the flux is between (normally 2)
   */
  static localIndex constexpr NUM_POINT_IN_FLUX = 2;

  /**
   * @brief Maximum number of points in a stencil (required to use static arrays in kernels)
   */
  static localIndex constexpr MAX_STENCIL_SIZE = 2;

  CellElementStencilTPFA();
//  ~CellElementStencilTPFA();

  void reserve( localIndex const size );

  void add( localIndex const numPts,
            localIndex  const * const elementRegionIndices,
            localIndex  const * const elementSubRegionIndices,
            localIndex  const * const elementIndices,
            real64 const * const weights,
            localIndex const connectorIndex );

  bool zero( localIndex const connectorIndex );

  localIndex size() const { return m_elementRegionIndices.size(0); }

  constexpr localIndex stencilSize( localIndex index ) const { return MAX_STENCIL_SIZE; }

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

#endif /* SRC_CORECOMPONENTS_FINITEVOLUME_CELLELEMENTSTENCILTPFA_HPP_ */
