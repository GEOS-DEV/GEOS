/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DiagonalOperator.hpp
 */

#ifndef GEOS_LINEARALGEBRA_UTILITIES_DIAGONALOPERATOR_HPP_
#define GEOS_LINEARALGEBRA_UTILITIES_DIAGONALOPERATOR_HPP_

#include "linearAlgebra/common/LinearOperator.hpp"

namespace geos
{

/**
 * @class DiagonalOperator
 * @brief Matrix-free diagonal operator, dst = diag * src.
 * @tparam VECTOR type of vector
 *
 * The operator keeps a reference to the vector of diagonal entries, which must outlive it.
 */
template< typename VECTOR >
class DiagonalOperator : public LinearOperator< VECTOR >
{
public:

  /// Alias for the base type
  using Base = LinearOperator< VECTOR >;

  /// Alias for the vector type
  using Vector = typename Base::Vector;

  /**
   * @brief Constructor.
   * @param diagonal the diagonal entries
   */
  explicit DiagonalOperator( Vector const & diagonal ):
    m_diagonal( diagonal )
  {}

  virtual void apply( Vector const & src, Vector & dst ) const override
  {
    dst.copy( src );
    dst.pointwiseProduct( m_diagonal );
  }

  virtual globalIndex numGlobalRows() const override { return m_diagonal.globalSize(); }

  virtual globalIndex numGlobalCols() const override { return m_diagonal.globalSize(); }

  virtual localIndex numLocalRows() const override { return m_diagonal.localSize(); }

  virtual localIndex numLocalCols() const override { return m_diagonal.localSize(); }

  virtual localIndex numLocalNonzeros() const override { return m_diagonal.localSize(); }

  virtual globalIndex numGlobalNonzeros() const override { return m_diagonal.globalSize(); }

  virtual MPI_Comm comm() const override { return m_diagonal.comm(); }

private:

  /// Diagonal entries
  Vector const & m_diagonal;
};

} // namespace geos

#endif /* GEOS_LINEARALGEBRA_UTILITIES_DIAGONALOPERATOR_HPP_ */
