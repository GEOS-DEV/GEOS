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

#ifndef GEOSX_LINEARALGEBRA_SOLVERS_MATRIXFREEPRECONDITIONERTWOLEVEL_HPP_
#define GEOSX_LINEARALGEBRA_SOLVERS_MATRIXFREEPRECONDITIONERTWOLEVEL_HPP_

#include "linearAlgebra/common/LinearOperator.hpp"
#include "linearAlgebra/common/PreconditionerBase.hpp"

namespace geosx
{

struct twoLevelStructuredMesh;

/**
 * @brief Common interface for two-level preconditioning operator
 * @tparam LAI linear algebra interface providing vectors, matrices and solvers
 */
template< typename LAI >
class MatrixFreePreconditionerTwoLevel : public PreconditionerBase< LAI >
{
public:

  /// Alias for base type
  using Base = PreconditionerBase< LAI >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename Base::Matrix;

  /**
   * @brief Constructor.
   * @param numComp number of components in the field
   * @param precond the actual preconditioner to apply to filtered matrix (ownership transferred)
   */
  MatrixFreePreconditionerTwoLevel ( LinearOperator< Vector > const & operatorFine,
                                     Matrix const & matCoarse,
		                     Vector const & diagonalInverse );

  virtual ~MatrixFreePreconditionerTwoLevel() override;

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
     m_coarseSolver.reset();
   }

  /**
   * @brief Apply operator to a vector.
   *
   * @param src Input vector (src).
   * @param dst Output vector (dst).
   */
  virtual void apply( Vector const & src,
                      Vector & dst ) const override;

  virtual globalIndex numGlobalRows() const override
  {
    return m_fineOperator.numGlobalRows();
  }

  virtual globalIndex numGlobalCols() const override
  {
    return m_fineOperator.numGlobalCols();
  }

  virtual localIndex numLocalRows() const override
  {
    return m_fineOperator.numLocalRows();
  }

  virtual localIndex numLocalCols() const override
  {
    return m_fineOperator.numLocalCols();
  }

  virtual MPI_Comm comm() const override
  {
    return m_fineOperator.comm();
  }

private:

  // Nested mesh data structure
  std::unique_ptr< twoLevelStructuredMesh > m_mesh;

  // Fine scale operator
  LinearOperator< Vector > const & m_fineOperator;

  // Coarse operator
  Matrix const & m_coarseMat;

  // Coarse solver
  std::unique_ptr< PreconditionerBase< LAI > > m_coarseSolver;

  // The inverse of the diagonal of the fine operator (Jacobi smoother)
  Vector const & m_diagInv;

  // Temporary vectors
  mutable Vector fineTmp;
  mutable Vector fineTmp2;
  mutable Vector fineRhs;
  mutable Vector coarseTmp;
  mutable Vector coarseRhs;
};

}

#endif //GEOSX_LINEARALGEBRA_SOLVERS_MATRIXFREEPRECONDITIONERTWOLEVEL_HPP_
