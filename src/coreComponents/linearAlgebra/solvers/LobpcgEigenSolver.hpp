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
 * @file LobpcgEigenSolver.hpp
 */

#ifndef GEOS_LINEARALGEBRA_SOLVERS_LOBPCGEIGENSOLVER_HPP_
#define GEOS_LINEARALGEBRA_SOLVERS_LOBPCGEIGENSOLVER_HPP_

#include "linearAlgebra/solvers/EigenSolverBase.hpp"

namespace geos
{

/**
 * @class LobpcgEigenSolver
 * @brief Locally optimal block preconditioned conjugate gradient eigensolver for K x = lambda M x.
 * @tparam VECTOR type of the vectors
 *
 * Computes the numEigenvalues smallest eigenvalues by minimizing the Rayleigh quotient on the span of the
 * current iterates, their preconditioned residuals and the previous search directions (Knyazev, SIAM J. Sci.
 * Comput. 23, 2001). The basis is orthonormalized in the M inner product with the SVQB algorithm, which also
 * removes numerically dependent directions. Converged vectors are soft-locked: they remain in the Rayleigh-Ritz
 * space but no longer generate residual and search directions.
 *
 * Unlike the Arnoldi solver, no linear system is solved to high accuracy: the method only needs a
 * preconditioner T that approximates (K - sigma M)^{-1} (one multigrid cycle is enough). In each iteration it
 * applies the preconditioner and M to the residual of every iterate that is not locked, because the convergence
 * test uses the preconditioned residual, and it applies K to the residual of every active vector. A locked
 * (converged) iterate is tested only every few iterations, and before the solver stops.
 * The shift only enters through the preconditioner, and the eigenvalues found are the smallest ones of the
 * pencil, so the shift should be at or below the smallest eigenvalue. K and M may be singular.
 *
 * The convergence test is on the M-norm of the preconditioned residual T ( K x - lambda M x ) of the M-normalized
 * vector, which estimates ||(K - sigma M)^{-1} r||_M, the error measure of the Arnoldi solver. Unlike the plain
 * residual, it is not limited by the rounding noise of K x for the eigenvalues close to zero.
 */
template< typename VECTOR >
class LobpcgEigenSolver : public GeneralizedEigenSolver< VECTOR >
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
  explicit LobpcgEigenSolver( EigenSolverParameters const & parameters ):
    Base( parameters )
  {}

  virtual EigenSolverResult solve( Problem const & problem,
                                   Vector const & prototype,
                                   std::vector< Vector > & modes ) const override;
};

} // namespace geos

#endif /* GEOS_LINEARALGEBRA_SOLVERS_LOBPCGEIGENSOLVER_HPP_ */
