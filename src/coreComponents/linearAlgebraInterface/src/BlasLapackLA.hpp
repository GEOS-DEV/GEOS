/*
 * BlasLapackLA.hpp
 *
 *  Created on: Apr 8, 2019
 *      Author: castelletto1
 */

#ifndef CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLASLAPACKLA_HPP_
#define CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLASLAPACKLA_HPP_

#include "common/DataTypes.hpp"
#include "Logger.hpp"

#include "cblas.h"

namespace geosx
{

class BlasLapackLA
{
  /**
   * \class BlasLapackLA
   * \brief This class contains a collection of BLAS and LAPACK linear
   *        algebra operations (dense) for GEOSX array1d and array2d
   */

public:

  //----------------------------------------------------------------------------
  //! @name Mathematical methods
  //@{

  /**
   * @brief Returns the 1-norm of the vector.
   */
  real64 vectorNorm1( array1d<real64> const & X ) const;

  /**
   * @brief Returns the two norm of the vector.
   */
  real64 vectorNorm2( array1d<real64> const & X ) const;

  /**
   * @brief Infinity-norm of the vector.
   */
  real64 vectorNormInf( array1d<real64> const & X ) const;

  /**
   * @brief Vector-Vector sum;
   * <tt>y</tt> = alpha*<tt>x</tt> + <tt>y</tt>.
   *
   * Computes (alpha*<tt>x</tt> + <tt>y</tt>) and overwrites the result on
   * <tt>y</tt>, with optional scaling.
   *
   * \param IN
   * <tt>x</tt> - GEOSX array1d.
   * \param [IN]
   * alpha - Optional scalar to multiply with <tt>Vec</tt>.
   *
   * \param INOUT
   * <tt>y</tt> - GEOSX array1d.
   *
   * @warning
   * Assumes that <tt>x</tt> and <tt>y</tt> have the same size.
   */
  void vectorVectorAdd( array1d<real64> const & X,
                        array1d<real64> & Y,
                        real64 const alpha = 1. );

  /**
   * @brief In-place scalar-vector product;
   * <tt>x</tt> = alpha*<tt>x<tt>
   *
   * \param IN
   * scalarThis - Scalar to multiply with \a this.
   */
  void vectorScale( array1d<real64> & X,
                    real64 alpha );

  /**
   * @brief Dot product of two vectors.
   *
   * \param IN
   * <tt>x</tt> - GEOSX array1d.
   * \param IN
   * <tt>y</tt> - GEOSX array1d.
   *
   */
  real64 vectorDot( array1d<real64> const & X,
                    array1d<real64> const & Y);

  /**
   * @brief Vector copy;
   * <tt>y</tt> = <tt>x<tt>
   *
   * \param IN
   * <tt>x</tt> - GEOSX array1d.
   *
   * \param INOUT
   * <tt>y</tt> - GEOSX array1d.
   *
   * @warning
   * Assumes that <tt>x</tt> and <tt>y</tt> have the same size.
   *
   */
  void vectorCopy( array1d<real64> const & X,
                   array1d<real64> & Y );
  //@}

  //----------------------------------------------------------------------------
  //! @name I/O methods
  //@{

  /**
   * @brief Print service method; defines behavior of ostream << operator.
   */
  void print(array1d<real64> const & X);

  //@}

};

}

#endif /* CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLASLAPACKLA_HPP_ */
