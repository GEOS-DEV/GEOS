/*
 * FaceElementStencil.hpp
 *
 *  Created on: Jul 1, 2019
 *      Author: settgast
 */

#ifndef SRC_CORECOMPONENTS_FINITEVOLUME_FACEELEMENTSTENCIL_HPP_
#define SRC_CORECOMPONENTS_FINITEVOLUME_FACEELEMENTSTENCIL_HPP_

#include "FluxStencil.hpp"

namespace geosx
{

class FaceElementStencil : public FluxStencil<CellDescriptor,real64>
{
public:
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


};

} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_FINITEVOLUME_FACEELEMENTSTENCIL_HPP_ */
