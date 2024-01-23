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

/**
 * @file SinglePhasePoromechanicsEmbeddedFractures.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASEPOROMECHANICSEMBEDDEDFRACTURES_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASEPOROMECHANICSEMBEDDEDFRACTURES_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

#define BRANCH_MGR_FSOLVER

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief SinglePhasePoromechanicsEmbeddedFractures strategy.
 *
 * dofLabel: 0 = displacement, x-component
 * dofLabel: 1 = displacement, y-component
 * dofLabel: 2 = displacement, z-component
 * dofLabel: 3 = wn
 * dofLabel: 4 = wt1
 * dofLabel: 5 = wt2
 * dofLabel: 6 = pressure (cell elem + fracture elems)
 *
 * Ingredients:
 * 1. F-points displacement (0,1,2), C-points pressure (3)
 * 2. F-points smoother: BoomerAMG, single V-cycle
 * 3. C-points coarse-grid/Schur complement solver: BoomerAMG
 * 4. Global smoother: none
 */
class SinglePhasePoromechanicsEmbeddedFractures : public MGRStrategyBase< 2 >
{
public:

  /**
   * @brief Constructor.
   */
  explicit SinglePhasePoromechanicsEmbeddedFractures( arrayView1d< int const > const & )
    : MGRStrategyBase( 7 )
  {
    // we keep u and p
    m_labels[0].push_back( 0 );
    m_labels[0].push_back( 1 );
    m_labels[0].push_back( 2 );
    m_labels[0].push_back( 6 );
    // we keep p
    m_labels[1].push_back( 6 );

    setupLabels();

    // Level 0
    m_levelFRelaxType[0]         = MGRFRelaxationType::gsElimWPivoting; // gaussian elimination for the dispJump block
    m_levelFRelaxIters[0]        = 1;
    m_levelGlobalSmootherType[0] = MGRGlobalSmootherType::none;
    m_levelInterpType[0]         = MGRInterpolationType::blockJacobi;
    m_levelRestrictType[0]       = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]   = MGRCoarseGridMethod::galerkin;

    // Level 1
    m_levelFRelaxType[1]         = MGRFRelaxationType::amgVCycle;
    m_levelGlobalSmootherType[1] = MGRGlobalSmootherType::none;
    m_levelInterpType[1]         = MGRInterpolationType::jacobi;
    m_levelRestrictType[1]       = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1]   = MGRCoarseGridMethod::nonGalerkin;

  }

  /**
   * @brief Setup the MGR strategy.
   * @param precond preconditioner wrapper
   * @param mgrData auxiliary MGR data
   */
  void setup( LinearSolverParameters::MGR const &,
              HyprePrecWrapper & precond,
              HypreMGRData & mgrData )
  {
    setReduction( precond, mgrData );

    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetPMaxElmts( precond.ptr, 0 ));

    // Configure the BoomerAMG solver used as F-relaxation for the second level
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &mgrData.mechSolver.ptr ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( mgrData.mechSolver.ptr, 0.0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( mgrData.mechSolver.ptr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxRowSum( mgrData.mechSolver.ptr, 1.0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetStrongThreshold( mgrData.mechSolver.ptr, 0.6 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( mgrData.mechSolver.ptr, 0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumFunctions( mgrData.mechSolver.ptr, 3 ) );

#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCoarsenType( mgrData.mechSolver.ptr, hypre::getAMGCoarseningType( LinearSolverParameters::AMG::CoarseningType::PMIS ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( mgrData.mechSolver.ptr, hypre::getAMGRelaxationType( LinearSolverParameters::AMG::SmootherType::chebyshev ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( mgrData.mechSolver.ptr, 1 ) );
#else
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxOrder( mgrData.mechSolver.ptr, 1 ) );
#endif
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetFSolverAtLevel( 1, precond.ptr, mgrData.mechSolver.ptr ) );

    // Configure the BoomerAMG solver used as mgr coarse solver for the pressure reduced system
    setPressureAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASEPOROMECHANICSEMBEDDEDFRACTURES_HPP_*/
