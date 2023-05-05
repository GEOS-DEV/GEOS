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
 * @file HybridSinglePhasePoromechanics.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRHYBRIDSINGLEPHASEPOROMECHANICS_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRHYBRIDSINGLEPHASEPOROMECHANICS_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief HybridSinglePhasePoromechanics strategy.
 *
 * Labels description stored in point_marker_array
 *   - dofLabel: 0 = nodal displacement, x-component
 *   - dofLabel: 1 = nodal displacement, y-component
 *   - dofLabel: 2 = nodal displacement, z-component
 *   - dofLabel: 3 = cell pressure
 *   - dofLabel: 4 = face pressure
 *
 * 2-level MGR reduction:
 *   - 1st level: eliminate displacements (0,1,2)
 *   - 2nd level: eliminate cell pressure (3)
 *   - The coarse grid solved with boomer AMG
 *
 */
class HybridSinglePhasePoromechanics : public MGRStrategyBase< 2 >
{
public:

  /**
   * @brief Constructor.
   */
  explicit HybridSinglePhasePoromechanics( arrayView1d< int const > const & )
    : MGRStrategyBase( 5 )
  {
    // Level 0: eliminate displacement degrees of freedom
    m_labels[0].push_back( 3 );
    m_labels[0].push_back( 4 );
    // Level 1: eliminate cell pressure degrees of freedom
    m_labels[1].push_back( 4 );

    setupLabels();

    // Level 0
    m_levelFRelaxMethod[0]     = MGRFRelaxationMethod::amgVCycle;
    m_levelInterpType[0]       = MGRInterpolationType::jacobi;
    m_levelRestrictType[0]     = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0] = MGRCoarseGridMethod::nonGalerkin;

    // Level 1
    m_levelFRelaxMethod[1]     = MGRFRelaxationMethod::singleLevel; //default, i.e. Jacobi (to be confirmed), switch to l1-jacobi
    // when hypre allows
    m_levelInterpType[1]       = MGRInterpolationType::jacobi;
    m_levelRestrictType[1]     = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1] = MGRCoarseGridMethod::galerkin;
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
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetCpointsByPointMarkerArray( precond.ptr,
                                                                 m_numBlocks, numLevels,
                                                                 m_numLabels, m_ptrLabels,
                                                                 mgrData.pointMarkers.data() ) );

    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelFRelaxMethod( precond.ptr, toUnderlyingPtr( m_levelFRelaxMethod ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( precond.ptr, toUnderlyingPtr( m_levelInterpType ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelRestrictType( precond.ptr, toUnderlyingPtr( m_levelRestrictType ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseGridMethod( precond.ptr, toUnderlyingPtr( m_levelCoarseGridMethod ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( precond.ptr, 1 ));
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetPMaxElmts( precond.ptr, 0 ));

    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &mgrData.coarseSolver.ptr ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( mgrData.coarseSolver.ptr, 0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( mgrData.coarseSolver.ptr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( mgrData.coarseSolver.ptr, 0.0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxOrder( mgrData.coarseSolver.ptr, 1 ) );

    mgrData.coarseSolver.setup = HYPRE_BoomerAMGSetup;
    mgrData.coarseSolver.solve = HYPRE_BoomerAMGSolve;
    mgrData.coarseSolver.destroy = HYPRE_BoomerAMGDestroy;

    // Configure the BoomerAMG solver used as F-relaxation for the first level
    setMechanicsFSolver( precond, mgrData );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRHYBRIDSINGLEPHASEPOROMECHANICS_HPP_*/
