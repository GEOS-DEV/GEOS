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
 * @file SolidMechanicsMixedVEM.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSOLIDMECHANICSMIXEDVEM_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSOLIDMECHANICSMIXEDVEM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief SolidMechanicsMixedVEM strategy.
 *
 * Hellinger-Reissner mixed VEM elasticity, traction moments then rigid motions,
 *
 *   A [sigma; u] = [M B^T ; -B 0] [sigma; u] = [0 ; -f],
 *
 * with N_sigma = 6 |F_h| and N_u = 6 |T_h|.
 *
 * Point markers 0-5 label the traction moments of T_h(f), 6-11 the RM(E) unknowns.
 *
 * One reduction with F = {sigma}, C = {u} and A_CC = 0: div Sigma_h(E) is fixed by the
 * traction moments, so every F unknown couples to C and a single level suffices.
 *
 * A_FF = blkdiag_f(D) + C_E, with D the equation (15) face Gram matrices and C_E element
 * local of rank at most 2 rank(Pi_E). Block Jacobi interpolation takes that face blocking
 * as the surrogate, Ahat_FF = blkdiag_f(A_FF); with restriction by injection
 *
 *   P = [W_p ; I], W_p = -Ahat_FF^{-1} B^T,
 *   A_C = A_CC + A_CF W_p = B Ahat_FF^{-1} B^T,
 *
 * symmetric positive definite, approximating the Schur complement S = B M^{-1} B^T.
 *
 * A_C is an interior penalty form: on interior faces (A_C u, u) = sum_f w_f |[u]_f|^2 with
 * w_f ~ 2 mu |f| / h_E, plus boundary traces, and cond(A_C) = O(h^-2).
 *
 * F-relaxation is Jacobi. The lambda degenerate direction is a constant hydrostatic stress,
 * which lies in ker B and is annihilated by the coarse correction, so a point smoother is
 * lambda uniform here.
 *
 * The coarse solve is one BoomerAMG V-cycle with Chebyshev relaxation and unknown based
 * coarsening on the six RM(E) functions. The near null space of A_C is the jump free
 * fields, which classical interpolation does not reproduce, so the coarse cycle is not
 * h-uniform and the iteration count grows slowly under refinement.
 *
 * The cycle is a fixed linear operator; a flexible outer Krylov method is not required.
 */
class SolidMechanicsMixedVEM : public MGRStrategyBase< 1 >
{
public:

  /// Number of Jacobi F-relaxation sweeps
  static constexpr HYPRE_Int numFRelaxSweeps = 3;

  /// Number of unknowns of RM(E) carried by the coarse operator
  static constexpr HYPRE_Int numCoarseFunctions = 6;

  /// Number of traction moments carried by each face
  static constexpr HYPRE_Int numFaceMoments = 6;

  /**
   * @brief Constructor.
   */
  explicit SolidMechanicsMixedVEM( arrayView1d< int const > const & )
    : MGRStrategyBase( 12 )
  {
    // Level 0: eliminate the traction moments, keep the rigid motions
    m_labels[0] = { 6, 7, 8, 9, 10, 11 };

    setupLabels();

    m_levelFRelaxType[0]         = MGRFRelaxationType::jacobi;
    m_levelFRelaxIters[0]        = numFRelaxSweeps;
    m_levelInterpType[0]         = MGRInterpolationType::blockJacobi;
    m_levelRestrictType[0]       = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]   = MGRCoarseGridMethod::galerkin;
  }

  /**
   * @brief Setup the MGR strategy.
   * @param mgrParams MGR configuration parameters
   * @param precond preconditioner wrapper
   * @param mgrData auxiliary MGR data
   */
  void setup( LinearSolverParameters::MGR const & mgrParams,
              HyprePrecWrapper & precond,
              HypreMGRData & mgrData )
  {
    GEOS_UNUSED_VAR( mgrParams );

    setReduction( precond, mgrData );

    // the equation (15) face Gram matrices are the block diagonal of A_FF
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetBlockJacobiBlockSize( precond.ptr, numFaceMoments ) );

    // one V-cycle on A_C, which carries six unknowns per element
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &mgrData.coarseSolver.ptr ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( mgrData.coarseSolver.ptr, 0.0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( mgrData.coarseSolver.ptr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumFunctions( mgrData.coarseSolver.ptr, numCoarseFunctions ) );
    // error operator I - p(A) A is partition independent, unlike hybrid Gauss-Seidel's rank local splitting
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( mgrData.coarseSolver.ptr, 16 ) );

    mgrData.coarseSolver.setup = HYPRE_BoomerAMGSetup;
    mgrData.coarseSolver.solve = HYPRE_BoomerAMGSolve;
    mgrData.coarseSolver.destroy = HYPRE_BoomerAMGDestroy;
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSOLIDMECHANICSMIXEDVEM_HPP_*/
