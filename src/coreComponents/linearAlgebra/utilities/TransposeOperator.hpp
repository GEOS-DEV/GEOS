/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TransposeOperator.hpp
 */

#ifndef GEOS_LINEARALGEBRA_UTILITIES_TRANSPOSEOPERATOR_HPP_
#define GEOS_LINEARALGEBRA_UTILITIES_TRANSPOSEOPERATOR_HPP_

#include "linearAlgebra/common/LinearOperator.hpp"

namespace geos
{

/**
 * @brief Simple class that wraps a matrix and represents its transpose as a linear operator.
 * @tparam MATRIX the linear algebra matrix type
 */
template< typename MATRIX >
class TransposeOperator : public LinearOperator< typename MATRIX::Vector >
{
public:

  /// Alias for base type
  using Base = LinearOperator< typename MATRIX::Vector >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = MATRIX;

  /**
   * @brief Constructor.
   * @param mat the underlying matrix
   */
  explicit TransposeOperator( Matrix const & mat )
    : Base(),
    m_matrix( mat )
  { }

  /**
   * @brief Apply operator to a vector.
   * @param src input vector (x)
   * @param dst output vector (b)
   *
   * @warning @p src and @p dst cannot alias the same vector (some implementations may allow this).
   */
  virtual void apply( Vector const & src, Vector & dst ) const override
  {
    m_matrix.applyTranspose( src, dst );
  }

  /**
   * @brief @return the number of global rows.
   */
  virtual globalIndex numGlobalRows() const override
  {
    return m_matrix.numGlobalCols();
  }

  /**
   * @brief @return the number of global columns.
   */
  virtual globalIndex numGlobalCols() const override
  {
    return m_matrix.numGlobalRows();
  }

  /**
   * @brief @return the number of local rows.
   */
  virtual localIndex numLocalRows() const override
  {
    return m_matrix.numLocalCols();
  }

  /**
   * @brief @return the number of local columns.
   */
  virtual localIndex numLocalCols() const override
  {
    return m_matrix.numLocalRows();
  }

  /**
   * @brief Get the MPI communicator the matrix was created with
   * @return MPI communicator of the underlying matrix
   */
  virtual MPI_Comm comm() const override
  {
    return m_matrix.comm();
  }

private:

  Matrix const & m_matrix;
};

}

#endif //GEOS_LINEARALGEBRA_UTILITIES_TRANSPOSEOPERATOR_HPP_
