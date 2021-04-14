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
 * @file LinearOperator.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_LINEAROPERATOR_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_LINEAROPERATOR_HPP_

#include "common/DataTypes.hpp"

namespace geosx
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
   * @brief Constructor
   */
  LinearOperator() = default;

  /**
   * @brief Destructor
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
   * @brief Get the number of global rows.
   * @return Number of global rows in the operator.
   */
  virtual globalIndex numGlobalRows() const = 0;

  /**
   * @brief Get the number of global columns.
   * @return Number of global columns in the operator.
   */
  virtual globalIndex numGlobalCols() const = 0;

  /**
   * @brief Get the number of local rows.
   * @return Number of local rows in the operator.
   */
  virtual localIndex numLocalRows() const = 0;

  /**
   * @brief Get the number of local columns.
   * @return Number of local columns in the operator.
   *
   * @note The use of term "local columns" refers not to physical partitioning of columns across ranks
   *       (as e.g. matrices are partitioned by rows and typically physically store all column entries),
   *       but to the partitioning of a compatible vector object that this operator can be applied to.
   */
  virtual localIndex numLocalCols() const = 0;

  /**
   * @brief Get the MPI communicator the matrix was created with
   * @return MPI communicator passed in @p create...()
   *
   * @note when build without MPI, may return anything
   *       (MPI_Comm will be a mock type defined in MpiWrapper)
   */
  virtual MPI_Comm getComm() const = 0;
};

}

#endif //GEOSX_LINEARALGEBRA_INTERFACES_LINEAROPERATOR_HPP_
