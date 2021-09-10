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
 * @file LagrangianContactMechanics.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRLAGRANGIACONTACTMECHANICS_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRLAGRANGIACONTACTMECHANICS_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geosx
{

namespace hypre
{

namespace mgr
{

/**
 * @brief LagrangianContactMechanics strategy
 *
 * Contact mechanics with face-centered lagrangian multipliers
 *
 * dofLabel: 0 = displacement, x-component
 * dofLabel: 1 = displacement, y-component
 * dofLabel: 2 = displacement, z-component
 * dofLabel: 3 = pressure
 *
 * Ingredients:
 * 1. F-points displacement (0,1,2), C-points pressure (3)
 * 2. F-points smoother: AMG, single V-cycle, separate displacemente components
 * 3. C-points coarse-grid/Schur complement solver: boomer AMG
 * 4. Global smoother: none
 *
 * @todo: (Sergey) I changed this from contiguous block setup to an equivalent marker array setup
 *        to make it a little easier (avoid passing DofManager, etc.). Will need to come up with a
 *        way to support both ways equally well. Maybe pass block offsets as constructor input.
 */
class LagrangianContactMechanics : public MGRStrategyBase< 1 >
{
public:

  /**
   * @brief Constructor.
   */
  explicit LagrangianContactMechanics( arrayView1d< int const > const & )
    : MGRStrategyBase( 2 )
  {
    // Level 0: all three displacements kept
    m_labels[0] = { 0, 1, 2 };
    setupLabels();

    m_levelFRelaxMethod[0] = 0; // Jacobi
    m_levelInterpType[0] = 2; // diagonal scaling (Jacobi)
    m_levelRestrictType[0] = 0; // injection
    m_levelCoarseGridMethod[0] = 0; // all

    m_numRestrictSweeps = 0;
    m_numInterpSweeps = 0;
    m_relaxType = 0;
    m_numRelaxSweeps = 0; // skip F-relaxation
    m_globalSmoothType = 0; // block Jacobi
    m_numGlobalSmoothSweeps = 1;
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
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCpointsByPointMarkerArray( precond.ptr,
                                                                  m_numBlocks, numLevels,
                                                                  m_numLabels, m_ptrLabels,
                                                                  mgrData.pointMarkers.data() ) );

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelRestrictType( precond.ptr, m_levelRestrictType ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNumRestrictSweeps( precond.ptr, m_numRestrictSweeps ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNumInterpSweeps( precond.ptr, m_numInterpSweeps ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNumRelaxSweeps( precond.ptr, m_numRelaxSweeps ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( precond.ptr, 1 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetRelaxType( precond.ptr, m_relaxType ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelFRelaxMethod( precond.ptr, m_levelFRelaxMethod ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( precond.ptr, m_levelInterpType ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseGridMethod( precond.ptr, m_levelCoarseGridMethod ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetGlobalsmoothType( precond.ptr, m_globalSmoothType ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetMaxGlobalsmoothIters( precond.ptr, m_numGlobalSmoothSweeps ) );

    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &mgrData.coarseSolver.ptr ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( mgrData.coarseSolver.ptr, 0 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( mgrData.coarseSolver.ptr, 1 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( mgrData.coarseSolver.ptr, 0.0 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxOrder( mgrData.coarseSolver.ptr, 1 ) );

    mgrData.coarseSolver.setup = HYPRE_BoomerAMGSetup;
    mgrData.coarseSolver.solve = HYPRE_BoomerAMGSolve;
    mgrData.coarseSolver.destroy = HYPRE_BoomerAMGDestroy;
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRLAGRANGIACONTACTMECHANICS_HPP_*/
