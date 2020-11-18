/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file Arnoldi.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_SOLVERS_ARNOLDI_HPP_
#define GEOSX_LINEARALGEBRA_SOLVERS_ARNOLDI_HPP_

#include "linearAlgebra/interfaces/VectorBase.hpp"
#include "linearAlgebra/interfaces/BlasLapackLA.hpp"
#include "linearAlgebra/interfaces/LinearOperator.hpp"

namespace geosx
{

/**
 * @brief NormalOperator Simple class to apply the operator A^T * A to a vector
 */
template< typename MATRIX, typename VECTOR >
class NormalOperator : public LinearOperator< VECTOR >
{
public:

  /**
   * @brief Sets the matrix
   * @param matrix the matrix
   * @param comm the MPI communicator
   */
  void set( MATRIX const & matrix, MPI_Comm const comm )
  {
    m_matrix = &matrix;
    m_comm = comm;
  }

  /**
   * @brief Returns the global number of rows
   * @return the global number of rows
   */
  globalIndex numGlobalRows() const override
  {
    return m_matrix->numGlobalRows();
  }

  /**
   * @brief Returns the global number of columns
   * @return the global number of columns
   */
  globalIndex numGlobalCols() const override
  {
    return m_matrix->numGlobalCols();
  }

  /**
   * @brief Returns the local number of rows
   * @return the local number of rows
   */
  localIndex numLocalRows() const
  {
    return m_matrix->numLocalRows();
  }

  /**
   * @brief Returns the communicator
   * @return the communicator
   */
  MPI_Comm const & getComm() const
  {
    return m_comm;
  }

  /**
   * @brief Applies the matrix and its transpose to a vector
   * @param x the input vector
   * @param y the output vector
   */
  void apply( VECTOR const & x, VECTOR & y ) const override
  {
    m_matrix->gemv( 1.0, x, 0.0, y, false );
    m_matrix->gemv( 1.0, y, 0.0, y, true );
  }

private:

  /// the matrix object
  MATRIX const * m_matrix;

  /// the communicator object
  MPI_Comm m_comm;
};

/**
 * @brief Function implementing the Arnoldi scheme to compute the largest eigenvalue
 * @param op the operator whose largest eigenvalue is required
 * @param m the number of iterations (size of the Krylov subspace)
 * @return the largest eigenvalue
 */
template< typename Operator >
real64 ArnoldiLargestEigenvalue( Operator const & op, localIndex const m = 4 )
{
  using Vector = typename Operator::Vector;

  localIndex const numGlobalRows = LvArray::integerConversion< localIndex >( op.numGlobalRows() );
  localIndex const numLocalRows = op.numLocalRows();
  localIndex const mInternal = ( m > numGlobalRows ) ? numGlobalRows : m;

  // Initialize data structure (Hessenberg matrix and Krylov subspace)
  array2d< real64, MatrixLayout::ROW_MAJOR_PERM > H( mInternal+1, mInternal );
  array1d< Vector > V( mInternal+1 );

  // Initial unitary vector
  V[0].createWithLocalSize( numLocalRows, op.getComm() );
  V[0].set( 1.0 / sqrt( static_cast< real64 >( numGlobalRows ) ) );

  for( localIndex j = 0; j < mInternal; ++j )
  {
    // Apply operator
    V[j+1].createWithLocalSize( numLocalRows, op.getComm() );
    op.apply( V[j], V[j+1] );
    // Arnoldi process
    for( localIndex i = 0; i <= j; ++i )
    {
      H( i, j ) = V[i].dot( V[j+1] );
      V[j+1].axpy( -H( i, j ), V[i] );
    }
    H( j+1, j ) = V[j+1].norm2();
    V[j+1].scale( 1.0 / H( j+1, j ) );
  }

  // Disregard the last entry and make the matrix square
  // Note: this is ok since we are using ROW_MAJOR_PERM
  H.resize( mInternal, mInternal );

  // Compute the eigenvalues
  array1d< std::complex< real64 > > lambda( mInternal );
  BlasLapackLA::matrixEigenvalues( H, lambda );

  // Find the largest eigenvalues
  real64 lambdaMax = 0.0;
  for( localIndex i = 0; i < mInternal; ++i )
  {
    lambdaMax = ( std::abs( lambda[i] ) > lambdaMax ) ? std::abs( lambda[i] ) : lambdaMax;
  }

  return lambdaMax;
}

} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_SOLVERS_ARNOLDI_HPP_
