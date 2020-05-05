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

/**
 * @file PreconditionerMultiphase.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_SOLVERS_MULTIPHASEFLOWPRECONDITIONER_HPP_
#define GEOSX_LINEARALGEBRA_SOLVERS_MULTIPHASEFLOWPRECONDITIONER_HPP_

#include <DofManager.hpp>
#include "common/DataTypes.hpp"
#include "linearAlgebra/solvers/PreconditionerBase.hpp"
#include "linearAlgebra/utilities/BlockOperatorWrapper.hpp"
#include "linearAlgebra/utilities/BlockVector.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

namespace geosx
{

/// Type of Schur complement approximation
enum class SchurApproximationType
{
  BLOCK_DIAGONAL,        // Quasi-IMPES-like reduction
  COLSUM_BLOCK_DIAGONAL  // True-IMPES-like reduction with colsum
};

template< typename LAI >
class MultiphasePreconditioner final : public PreconditionerBase< LAI >
{
public:

  using Base = PreconditionerBase< LAI >;
  using Vector = typename Base::Vector;
  using Matrix = typename Base::Matrix;

  MultiphasePreconditioner( SchurApproximationType schurType,
                            LinearSolverParameters parameters );

  virtual ~MultiphasePreconditioner() override;

  void setPrimaryDofComponents( DofManager const & dofManager,
                                std::vector< DofManager::SubComponent > primaryDofs );

  virtual void compute( Matrix const & mat,
                        DofManager const & dofManager ) override;

  virtual void apply( Vector const & src, Vector & dst ) const override;

  virtual void clear() override;

private:

  /// (Re-)create storage for intermediate vectors used in solution phase
  void initialize( DofManager const & dofManager );

  /// Create and populate sparsity pattern for Schur complement matrix
  void createReduction( DofManager const & dofManager );

  void computeColsumReduction( DofManager const & dofManager, Matrix const & PS, Matrix const & SS );

  /// Type of Schur complement reduction to use
  SchurApproximationType m_schurType;

  /// Dof components picked in the first stage
  std::vector< DofManager::SubComponent > m_dofComponentsPrimary;

  /// Remaining dof components
  std::vector< DofManager::SubComponent > m_dofComponentsOther;

  /// Pressure prolongation operator
  Matrix m_prolongationOpPrimary;

  /// Pressure restriction operator
  Matrix m_restrictionOpPrimary;

  /// Other variable prolongation operator
  Matrix m_prolongationOpOther;

  /// Other variable restriction operator
  Matrix m_restrictionOpOther;

  /// System reduction operator
  Matrix m_reductionOp;

  /// Colsum result for pressure
  array1d< Vector > m_colsumResultPrimary;

  /// Colsum "operator" for pressure
  array1d< Vector > m_colsumOpPrimary;

  /// Colsum result for other dofs
  array1d< Vector > m_colsumResultOther;

  /// Colsum "operator" for other dofs
  array1d< Vector > m_colsumOpOther;

  /// Contains @f$ A_{pp} - \tilde{A}_{ps} \tilde{A}_{ss}^{-1} A_{sp} @f$
  Matrix m_systemMatrixPrimary;

  /// Parameters for first stage preconditioner
  LinearSolverParameters m_parameters;

  /// Preconditioner for first stage system
  std::unique_ptr< PreconditionerBase< LAI > > m_solverPrimary;

  /// Pressure right-hand side vector
  mutable Vector m_rhsPrimary;

  /// Pressure solution vector
  mutable Vector m_solPrimary;
};

} //namespace geosx

#endif //GEOSX_LINEARALGEBRA_SOLVERS_MULTIPHASEFLOWPRECONDITIONER_HPP_
