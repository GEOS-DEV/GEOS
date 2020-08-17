/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_LINEARALGEBRA_SOLVERS_PRECONDITIONERJACOBI_HPP_
#define GEOSX_LINEARALGEBRA_SOLVERS_PRECONDITIONERJACOBI_HPP_

#include "linearAlgebra/interfaces/LinearOperator.hpp"
#include "linearAlgebra/solvers/PreconditionerBase.hpp"

namespace geosx
{

/**
 * @brief Common interface for identity preconditioning operator
 * @tparam LAI linear algebra interface providing vectors, matrices and solvers
 */
template< typename LAI >
class PreconditionerJacobi : public PreconditionerBase< LAI >
{
public:

  /// Alias for base type
  using Base = PreconditionerBase< LAI >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename Base::Matrix;

  virtual ~PreconditionerJacobi() = default;

  /**
   * @brief Compute the preconditioner from a matrix.
   * @param mat the matrix to precondition.
   */
  virtual void compute( Matrix const & mat ) override
  {
    GEOSX_LAI_ASSERT( mat.ready() );
    m_diag = new Vector();
    m_diag->createWithLocalSize( mat.numLocalRows(), mat.getComm() );
    mat.extractDiagonal( *m_diag );
  }

  /**
   * @brief Compute the preconditioner from a matrix
   * @param mat the matrix to precondition
   * @param dofManager the Degree-of-Freedom manager associated with matrix
   */
  virtual void compute( Matrix const & mat,
                        DofManager const & dofManager ) override
  {
    GEOSX_UNUSED_VAR( dofManager );
    compute( mat );
  }

  /**
   * @brief Clean up the preconditioner setup.
   *
   * Releases memory used and allows the matrix to be deleted cleanly.
   * This method should be called before the matrix used to compute the preconditioner
   * goes out of scope or is re-created. Some implementations require the matrix
   * to outlive the preconditioner (for example, Trilinos/ML may crash the program if
   * deleted after the matrix).
   *
   * @note Should be properly overridden in derived classes, which may call this method.
   */
  virtual void clear() override
  {
    delete m_diag;
    m_diag = nullptr;
  }

  /**
   * @brief Get the number of global rows.
   * @return Number of global rows in the operator.
   */
  virtual globalIndex numGlobalRows() const override final
  {
    return m_diag->globalSize();
  }

  /**
   * @brief Get the number of global columns.
   * @return Number of global columns in the operator.
   */
  virtual globalIndex numGlobalCols() const override final
  {
    return m_diag->globalSize();
  }

  /**
   * @brief Apply operator to a vector.
   *
   * @param src Input vector (src).
   * @param dst Output vector (dst).
   */
  virtual void apply( Vector const & src,
                      Vector & dst ) const override
  {
    GEOSX_LAI_ASSERT_EQ( this->numGlobalRows(), dst.globalSize() );
    GEOSX_LAI_ASSERT_EQ( this->numGlobalCols(), src.globalSize() );
    GEOSX_LAI_ASSERT_EQ( src.localSize(), dst.localSize() );

    real64 const * const src_data = src.extractLocalVector();
    real64 * const dst_data = dst.extractLocalVector();
    real64 const * const diag_data = m_diag->extractLocalVector();
    dst.copy( src );
    for( localIndex i = 0; i < src.localSize(); ++i )
    {
      dst_data[i] = src_data[i] / diag_data[i];
    }
  }

private:

  /// The diagonal of the matrix
  Vector * m_diag;
};

}

#endif //GEOSX_LINEARALGEBRA_SOLVERS_PRECONDITIONERJACOBI_HPP_
