/*
 * LapackInterface.hpp
 *
 *  Created on: Feb 4, 2019
 *      Author: castelletto1
 */

#ifndef CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_LAPACKINTERFACE_HPP_
#define CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_LAPACKINTERFACE_HPP_

#include "LapackMatrix.hpp"
//#include "LapackVector.hpp"

namespace geosx
{

/**
 * \class LapackInterface
 * \brief This class holds aliases based on LAPACK and BLAS library.
 */

class LapackInterface
{
public:

  // Lapack matrix and vector wrappers
  using DenseMatrix = LapackMatrix;
//  using DenseVector = LapackVector;

  //! @name Constructor/Destructor Methods
  //@{
  /**
   * @brief Empty constructor.
   */
  LapackInterface() = default;

  /**
   * @brief Destructor.
   *
   */
  ~LapackInterface() = default;
  //@}

};

} /* namespace geosx */



#endif /* CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_LAPACKINTERFACE_HPP_ */
