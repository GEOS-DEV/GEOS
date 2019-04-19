/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file BlasVector.hpp
 */

#ifndef CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLASVECTOR_HPP_
#define CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLASVECTOR_HPP_

#include "common/DataTypes.hpp"
#include "Logger.hpp"

#include <random>

#include "cblas.h"
//#include "lapacke.h"
//
#define __cplusplus_temp __cplusplus
#undef __cplusplus
#define __cplusplus 199711L
#include "lapacke.h"
#undef __cplusplus
#define __cplusplus __cplusplus_temp
#undef __cplusplus_temp

namespace geosx
{

/**
 * \class BlasVector
 * \brief This class creates and provides support for manipulating a BLAS-style
 *        column vector
 */

class BlasVector
{
public:

  //----------------------------------------------------------------------------
  //! @name Constructor/Destructor Methods
  //@{

  /**
   * @brief Default constructor.
   */
  BlasVector();

  /**
   * @brief Shaped constructor; defines a variable-sized vector.
   *
   * \param IN
   * length - vector size.
   *
   * @note Values are not initialized to 0.
   */
  BlasVector( localIndex length );

  /**
   * @brief Copy constructor.
   *
   * \param IN
   * src - BlasVector
   *
   */
  BlasVector( BlasVector const & vec );

  /**
   * @brief Matrix destructor.
   */
  ~BlasVector();

  //@}

  //----------------------------------------------------------------------------
  //! @name Shaping/sizing/permuting methods
  //@{

  /**
   * @brief Resize vector.
   *
   * \param IN
   * length - Vector length
   *
   * @note Values are not initialized to 0.
   */
  void resize( localIndex length );

  /**
   * @brief Reinitialize the vector.
   *
   * Sets all elements to zero.
   *
   */
  void zero();

  /**
   * @brief Set vector elements to uniformly distributed random
   *        entries in the interval [0.0, 1.0).
   *
   */
  void rand();

  /**
   * @brief Set vector elements to uniformly distributed random
   *        entries in the interval [\a rangeFrom, \a rangeTo)r.
   *
   * * \param IN
   * \a rangeFrom - Interval left bound
   * \a rangeTo - Interval right bound
   *
   */
  void rand(real64 const rangeFrom,
            real64 const rangeTo);

  /**
   * @brief Operator that sets all matrix entries equal to a constant value
   *
   * * \param IN
   * \a value - Constant value to be assigned to each entry of the vector
   *
   */
  BlasVector &operator=( double value );

  /**
   * @brief Rearranges the vector coefficient as specified by a permutation
   *        vector
   *
   * * \param IN
   * \a permutationVector - contains the permutation vector
   * \a forwardPermuation - optional argument
   *
   * @note Given an M*1 vector V and a permutation vector PERM
   *
   *       Forward permutation (forwardPermutation = true):
   *       V(PERM(I),*) is moved V(I,*) for I = 1,2,...,M.
   *
   *       Backward permutation (forwardPermutation = false):
   *       V(I,*) is moved V(PERM(I),*) for I = 1,2,...,M.
   *
   * @warning Permutation indeces start from 1.
   *
   */
  void permute( array1d<int> permutationVector,
                const bool forwardPermutation = true );
  //@}

  //----------------------------------------------------------------------------
  //! @name Mathematical methods
  //@{

  /**
   * @brief Returns the 1-norm of the vector.
   */
  real64 norm1() const;

  /**
   * @brief Returns the two norm of the vector.
   */
  real64 norm2() const;

  /**
   * @brief Infinity-norm of the vector.
   */
  real64 normInf() const;

  /**
   * @brief Vector-Vector sum;
   * \a this = scalarVec*<tt>Vec</tt> + \a this.
   *
   * Computes (scalarVec*<tt>Vec</tt> + \a this) and overwrites the result on the
   * current \a this, with optional scaling.
   *
   * \param IN
   * <tt>Vec</tt> - BlasVector vector.
   * \param [IN]
   * scalarVec - Optional scalar to multiply with <tt>Vec</tt>.
   *
   * @warning
   * Assumes that <tt>Vec</tt> and \a this have the same size.
   */
  void vectorAdd( BlasVector const & Vec,
                  real64 const scalarVec = 1. );

  /**
   * @brief In-place scalar-vector product;
   * \a this = scalarThis* \a this
   *
   * \param IN
   * scalarThis - Scalar to multiply with \a this.
   */
  void scale( real64 scalarThis );

  /**
   * @brief Dot product with the vector vec.
   *
   * \param vec EpetraVector to dot-product with.
   *
   */
  real64 dot( BlasVector const &vec );

  /**
   * @brief Update vector \a this as \a this = <tt>vec</tt>.
   *
   * @note The naming convention follows the BLAS library.
   *
   * \param x EpetraVector to copy.
   *
   */
  void copy( BlasVector const &vec );

  //@}

  //----------------------------------------------------------------------------
  //! @name Data Accessor methods
  //@{

  /**
   * @brief Coefficient access function.
   *
   * The parentheses operator returns the reference to the coefficient in
   * position (index).
   */
  inline real64 & operator ()( localIndex index );

  /**
   * @brief Constant coefficient access function.
   *
   * The parentheses operator returns the reference to the coefficient in
   * position (index).
   */
  inline real64 const & operator ()( localIndex index ) const;

  /**
   * @brief Returns vector size.
   */
  localIndex getSize() const;

  /**
   * @brief Returns a const pointer to the underlying
   *        array1d<real64> object.
   */
  array1d<real64> const * getValues() const;

  /**
   * @brief Returns a non-const pointer to the underlying
   *        array1d<real64> object.
   */
  array1d<real64>* getValues();

  //@}

  //@}

  //----------------------------------------------------------------------------
  //! @name I/O methods
  //@{

  /**
   * @brief Print service method; defines behavior of ostream << operator.
   */
  void print();

  //@}

protected:

  localIndex m_size = 0; ///< Number of vector entries (lenght)
  array1d<real64> m_values; ///< array1d storing vector entries

};

// Inline methods
inline real64 & BlasVector::operator ()( localIndex index )
{
  GEOS_ASSERT_MSG( 0 <= index &&
                   index <= m_size,
                   "Requested value out of bounds" );
  return m_values[index];
}

inline real64 const & BlasVector::operator ()( localIndex index) const
{
  GEOS_ASSERT_MSG( 0 <= index &&
                   index <= m_size,
                   "Requested value out of bounds" );
  return m_values[index];
}

} // namespace geosx


#endif /* CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLASVECTOR_HPP_ */
