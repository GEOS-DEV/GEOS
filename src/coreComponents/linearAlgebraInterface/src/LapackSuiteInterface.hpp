/*
 * LapackInterface.hpp
 *
 *  Created on: Feb 4, 2019
 *      Author: castelletto1
 */

#ifndef CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_LAPACKSUITEINTERFACE_HPP_
#define CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_LAPACKSUITEINTERFACE_HPP_

#include "BlasMatrix.hpp"
#include "BlasVector.hpp"

namespace geosx
{

/**
 * \class LapackInterface
 * \brief This class holds aliases based on LAPACK and BLAS library.
 */

class LapackSuiteInterface
{
public:

  // Lapack matrix and vector wrappers
  using DenseMatrix = BlasMatrix;
  using Vector = BlasVector;

  //! @name Constructor/Destructor Methods
  //@{
  /**
   * @brief Empty constructor.
   */
  LapackSuiteInterface() = default;

  /**
   * @brief Destructor.
   *
   */
  ~LapackSuiteInterface() = default;
  //@}

};

} /* namespace geosx */



#endif /* CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_LAPACKSUITEINTERFACE_HPP_ */
