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
 * @file ArnoldiEigenSolver.hpp
 */

#ifndef GEOS_LINEARALGEBRA_SOLVERS_ARNOLDIEIGENSOLVER_HPP_
#define GEOS_LINEARALGEBRA_SOLVERS_ARNOLDIEIGENSOLVER_HPP_

#include "linearAlgebra/solvers/EigenSolverBase.hpp"

namespace geos
{

/**
 * @class ArnoldiEigenSolver
 * @brief Shift-and-invert block Krylov-Schur eigensolver for K x = lambda M x.
 * @tparam VECTOR type of the vectors
 *
 * The solver builds an M-orthonormal Krylov basis of the operator T = (K - sigma M)^{-1} M, which is
 * self-adjoint for the M inner product, and extracts the Ritz pairs of largest modulus (the eigenvalues
 * lambda closest to the shift sigma). The basis is restarted with the Krylov-Schur (thick-restart) technique.
 * For symmetric pencils the Arnoldi process reduces to the Lanczos process; the implementation keeps the
 * full projected matrix and uses full M-reorthogonalization so it does not depend on this structure. This
 * is the same family of method as the implicitly restarted Arnoldi/Lanczos method of ARPACK in shift-invert
 * mode.
 *
 * Each iteration costs one application of the shift-and-invert operator per basis vector, which dominates
 * the cost. The accuracy of the operator application limits the attainable accuracy of the eigenpairs.
 */
template< typename VECTOR >
class ArnoldiEigenSolver : public GeneralizedEigenSolver< VECTOR >
{
public:

  /// Base type
  using Base = GeneralizedEigenSolver< VECTOR >;

  /// Type of the vectors
  using Vector = typename Base::Vector;

  /// Type of the problem description
  using Problem = typename Base::Problem;

  /// Type of the constraint space
  using ConstraintSpace = typename Base::ConstraintSpace;

  /**
   * @brief Constructor.
   * @param parameters solver parameters
   */
  explicit ArnoldiEigenSolver( EigenSolverParameters const & parameters ):
    Base( parameters )
  {}

  virtual EigenSolverResult solve( Problem const & problem,
                                   Vector const & prototype,
                                   stdVector< Vector > & modes ) const override;
};

} // namespace geos

#endif /* GEOS_LINEARALGEBRA_SOLVERS_ARNOLDIEIGENSOLVER_HPP_ */
