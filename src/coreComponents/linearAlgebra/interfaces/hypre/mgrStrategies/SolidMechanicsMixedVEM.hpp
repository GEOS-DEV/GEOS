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
 * Mixed virtual element elasticity in saddle point form,
 *
 *   [ K  B^T ] [ sigma ]   [  g ]
 *   [ B   0  ] [   u   ] = [ -f ],
 *
 * with six face traction unknowns per face and six rigid body motion unknowns per cell.
 *
 * dofLabel: 0-5  = face stress, the six traction modes of T_h(f)
 * dofLabel: 6-11 = cell rigid body motion, the six modes of RM(E)
 *
 * Ingredients:
 * 1. F-points face stress (0-5), C-points rigid body motions (6-11)
 * 2. F-points smoother: l1-Jacobi, the stress block is symmetric positive definite
 * 3. C-points coarse-grid/Schur complement solver: boomer AMG on six functions per cell
 * 4. Global smoother: none
 *
 * The reduction follows the structure of K. Element by element the stiffness is
 *
 *   K_E = [P_E ; Chat_E]^T [ W  -I ] [P_E ; Chat_E]  +  kappa_E h_E blockdiag_f( M2_f ),
 *                          [ -I  0 ]
 *
 * because the consistency term and the two stabilization cross terms all factor through
 * the six rows of P_E: the first part has rank at most twelve per element, while the
 * stabilization Gram term is exactly block diagonal with one 6x6 block M2_f per face.
 *
 * So the natural pivot for eliminating the stress is the 6x6 face block, not a scalar. The
 * six traction modes of a face are two constant tangential, one constant normal, one in
 * plane rotational and two linear normal functionals; they carry different scalings and
 * M2_f couples the last two. A scalar row sum averages all six into one number and loses
 * exactly the part of K that is genuinely diagonal, whereas a block Jacobi pivot keeps it.
 *
 * With D = blockdiag_f(K) the Galerkin coarse operator is
 *
 *   S = 0 - (-B) D^{-1} B^T = B D^{-1} B^T,
 *
 * which is sparse on the cell adjacency graph exactly as the scalar version is, since a
 * block diagonal inverse is still block diagonal: the better approximation costs no
 * fill in. It leaves 6 N_c unknowns instead of 6 N_f + 6 N_c, and it is definite rather
 * than indefinite because the balance is assembled as -(div sigma, v) = (f, v).
 *
 * K is symmetric positive definite, so F-relaxation is a full multigrid V-cycle rather
 * than a point smoother; with a rank twelve dense part per element a point smoother
 * leaves far too much of the stress error behind.
 */
class SolidMechanicsMixedVEM : public MGRStrategyBase< 1 >
{
public:

  /**
   * @brief Constructor.
   */
  explicit SolidMechanicsMixedVEM( arrayView1d< int const > const & )
    : MGRStrategyBase( 12 )
  {
    // Level 0: eliminate the face stress, keep the rigid body motions
    m_labels[0] = { 6, 7, 8, 9, 10, 11 };

    setupLabels();

    // Level 0
    m_levelFRelaxType[0]         = MGRFRelaxationType::amgVCycle;
    m_levelFRelaxIters[0]        = 1;
    m_levelInterpType[0]         = MGRInterpolationType::l1jacobi;
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
    setReduction( precond, mgrData );

    // The reduced system carries six rigid body motions per cell, so multigrid treats it
    // as a system of six functions rather than a scalar problem.
    //
    // That alone is not enough. S is a discrete elasticity operator on the skeleton of
    // cells and its near null space is the global rigid body motions, which on element E
    // read (a + omega ^ x_E, omega). The three translations are constant per function and
    // an unknown based coarsening reproduces them, but the three rotations are linear in
    // the cell center and it cannot. Left out, the coarse correction never damps them and
    // the iteration count grows with the mesh. Passing them as interpolation vectors is
    // what makes the reduction mesh independent.
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &mgrData.coarseSolver.ptr ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( mgrData.coarseSolver.ptr, 0.0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( mgrData.coarseSolver.ptr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxRowSum( mgrData.coarseSolver.ptr, 1.0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetStrongThreshold( mgrData.coarseSolver.ptr, 0.6 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( mgrData.coarseSolver.ptr, 0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxOrder( mgrData.coarseSolver.ptr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumFunctions( mgrData.coarseSolver.ptr, 6 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetFilterFunctions( mgrData.coarseSolver.ptr, mgrParams.separateComponents ) );

    if( !mgrData.nearNullSpace.empty() )
    {
      HYPRE_Int const nodal                 = 4;
      HYPRE_Int const nodalDiag             = 1;
      HYPRE_Int const relaxCoarse           = 8;
      HYPRE_Int const interpVecVariant      = 2;
      HYPRE_Int const qMax                  = 4;
      HYPRE_Int const smoothInterpVectors   = 1;
      HYPRE_Int const interpRefine          = 1;

      GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNodal( mgrData.coarseSolver.ptr, nodal ) );
      GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNodalDiag( mgrData.coarseSolver.ptr, nodalDiag ) );
      GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleRelaxType( mgrData.coarseSolver.ptr, relaxCoarse, 3 ) );
      GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetInterpVecVariant( mgrData.coarseSolver.ptr, interpVecVariant ) );
      GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetInterpVecQMax( mgrData.coarseSolver.ptr, qMax ) );
      GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetSmoothInterpVectors( mgrData.coarseSolver.ptr, smoothInterpVectors ) );
      GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetInterpRefine( mgrData.coarseSolver.ptr, interpRefine ) );
      GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetInterpVectors( mgrData.coarseSolver.ptr,
                                                             mgrData.nearNullSpace.size(),
                                                             mgrData.nearNullSpace.data() ) );
    }

    mgrData.coarseSolver.setup = HYPRE_BoomerAMGSetup;
    mgrData.coarseSolver.solve = HYPRE_BoomerAMGSolve;
    mgrData.coarseSolver.destroy = HYPRE_BoomerAMGDestroy;
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSOLIDMECHANICSMIXEDVEM_HPP_*/
