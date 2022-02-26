/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file LinearOperator.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_LINEAROPERATOR_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_LINEAROPERATOR_HPP_

#include "common/DataTypes.hpp"

namespace geos
{

/**
 * @brief Abstract base class for linear operators.
 *
 * @tparam VECTOR Type of vector this operator can be applied to
 */
template< typename VECTOR >
class LinearOperator
{
public:

  /// Alias for template parameter
  using Vector = VECTOR;

  /**
   * @brief Constructor.
   */
  LinearOperator() = default;

  /**
   * @brief Destructor.
   */
  virtual ~LinearOperator() = default;

  /**
   * @brief Apply operator to a vector, <tt>dst = this(src)</tt>.
   * @param src input vector
   * @param dst output vector
   *
   * @warning @p src and @p dst cannot alias the same vector (some implementations may allow this).
   */
  virtual void apply( Vector const & src, Vector & dst ) const = 0;

  /**
   * @brief Compute residual <tt>r = b - this(x)</tt>.
   *
   * @param x Input solution.
   * @param b Input right hand side.
   * @param r Output residual.
   *
   * @warning @p b and @p x may alias the same vector.
   *          @p r cannot alias any of the other two vectors (some implementations may allow this).
   */
  virtual void residual( Vector const & x, Vector const & b, Vector & r ) const
  {
    this->apply( x, r );
    r.axpby( 1.0, b, -1.0 );
  }

  /**
   * @brief @return the number of global rows.
   */
  virtual globalIndex numGlobalRows() const = 0;

  /**
   * @brief @return the number of global columns.
   */
  virtual globalIndex numGlobalCols() const = 0;

  /**
   * @brief @return the total number of nonzero entries in the operator.
   *
   * @note While "non-zero entries" is not a well-defined term for the abstract concept of a linear operator,
   *       in practice this value is needed to track runtime/memory complexity of the operator implementation.
   *       The function returns the number of nonzero floating-point values stored/used when applying to a vector.
   *       A matrix-free operator implemented without any additional FP storage may either return 0,
   *       or the (approximate) number of floating-point multiply operations performed by apply().
   */
  virtual globalIndex numGlobalNonzeros() const
  {
    return 0;
  }

  /**
   * @brief @return the number of local rows.
   */
  virtual localIndex numLocalRows() const = 0;

  /**
   * @brief @return the number of local columns.
   *
   * @note The use of term "local columns" refers not to physical partitioning of columns across ranks
   *       (as e.g. matrices are partitioned by rows and typically physically store all column entries),
   *       but to the partitioning of a compatible vector object that this operator can be applied to.
   */
  virtual localIndex numLocalCols() const = 0;

  /**
   * @brief @return the number of nonzero entries in the local portion of the operator.
   */
  virtual localIndex numLocalNonzeros() const
  {
    return 0;
  }

  /**
   * @brief Get the MPI communicator the matrix was created with
   * @return MPI communicator passed in @p create...()
   *
   * @note when build without MPI, may return anything
   *       (MPI_Comm will be a mock type defined in MpiWrapper)
   */
  virtual MPI_Comm comm() const = 0;
};

}

#endif //GEOS_LINEARALGEBRA_INTERFACES_LINEAROPERATOR_HPP_
