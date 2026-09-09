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
 * Hellinger-Reissner mixed virtual element elasticity, ordered as traction moments then
 * rigid motion moments,
 *
 *   A [sigma; u] = [M B^T ; B 0] [sigma; u] = [0 ; -f],
 *
 * of order N_sigma + N_u with N_sigma = 6 |F_h| and N_u = 6 |T_h|.
 *
 * dofLabel: 0-5  = the six traction moments of T_h(f)
 * dofLabel: 6-11 = the six rigid motion moments of RM(E)
 *
 * One reduction with F = {sigma}, C = {u} and A_CC = 0. The divergence of a discrete
 * stress is a function of the traction moments alone, so every stress unknown is touched
 * by the constraint and the reduction is single level.
 *
 * With the interpolation surrogate Ahat_FF = diag(A_FF) the prolongation is
 * P = [W_p ; I], W_p = -diag(M)^{-1} B^T, restriction is by injection, and the Galerkin
 * coarse operator collapses to
 *
 *   A_C = A_CC + A_CF W_p = -B diag(M)^{-1} B^T,
 *
 * the discrete approximation of the Schur complement S = -B M^{-1} B^T. S is symmetric
 * negative definite, cell based, six unknowns per element on a face neighbour stencil,
 * with cond(S) = O(h^-2), so one BoomerAMG V-cycle is an h-uniform coarse solver.
 *
 * A constant hydrostatic stress lies in ker B and carries a_h(sigma, sigma) that vanishes
 * as lambda grows, so the incompressible degeneracy sits in A_FF and never enters A_C.
 * Incompressibility robustness is therefore a property of the F-relaxation: a point
 * smoother there tracks sqrt(2 mu + 3 lambda), an inner multigrid solve does not.
 *
 * The sweep count must be odd.
 *
 * The preconditioner is one MGR cycle and contains an inner multigrid solve, so it is not
 * a fixed linear operator and the outer Krylov method must be the flexible variant.
 */
class SolidMechanicsMixedVEM : public MGRStrategyBase< 1 >
{
public:

  /// Number of F-relaxation sweeps, odd by construction
  static constexpr HYPRE_Int numFRelaxSweeps = 3;

  /// Number of unknowns of RM(E) carried by the coarse operator
  static constexpr HYPRE_Int numCoarseFunctions = 6;

  /**
   * @brief Constructor.
   */
  explicit SolidMechanicsMixedVEM( arrayView1d< int const > const & )
    : MGRStrategyBase( 12 )
  {
    static_assert( numFRelaxSweeps % 2 == 1, "MGR F-relaxation sweeps must be odd" );

    // Level 0: eliminate the traction moments, keep the rigid motions
    m_labels[0] = { 6, 7, 8, 9, 10, 11 };

    setupLabels();

    m_levelFRelaxType[0]         = MGRFRelaxationType::amgVCycle;
    m_levelFRelaxIters[0]        = numFRelaxSweeps;
    m_levelInterpType[0]         = MGRInterpolationType::jacobi;
    m_levelRestrictType[0]       = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]   = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[0] = MGRGlobalSmootherType::none;
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

    // one V-cycle on A_C, which carries six unknowns per element
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &mgrData.coarseSolver.ptr ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( mgrData.coarseSolver.ptr, 0.0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( mgrData.coarseSolver.ptr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumFunctions( mgrData.coarseSolver.ptr, numCoarseFunctions ) );
    // error operator I - p(A) A is partition independent, unlike hybrid Gauss-Seidel's rank local splitting
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( mgrData.coarseSolver.ptr, 16 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetChebyOrder( mgrData.coarseSolver.ptr, 2 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( mgrData.coarseSolver.ptr, 0 ) );

    mgrData.coarseSolver.setup = HYPRE_BoomerAMGSetup;
    mgrData.coarseSolver.solve = HYPRE_BoomerAMGSolve;
    mgrData.coarseSolver.destroy = HYPRE_BoomerAMGDestroy;
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSOLIDMECHANICSMIXEDVEM_HPP_*/
