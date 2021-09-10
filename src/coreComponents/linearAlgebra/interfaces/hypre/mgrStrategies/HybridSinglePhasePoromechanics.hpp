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

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRHYBRIDSINGLEPHASEPOROMECHANICS_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRHYBRIDSINGLEPHASEPOROMECHANICS_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geosx
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

    m_levelFRelaxMethod[0] = 2;  // AMG V-cycle
    m_levelFRelaxMethod[1] = 18; // l1-Jacobi

    m_levelInterpType[0] = 2;       // diagonal scaling (Jacobi)
    m_levelCoarseGridMethod[0] = 1; // diagonal sparsification
    m_levelInterpType[1] = 2;       // diagonal scaling (Jacobi)
    m_levelCoarseGridMethod[1] = 0; // Galerkin coarse grid computation using RAP
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

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelFRelaxMethod( precond.ptr, m_levelFRelaxMethod ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( precond.ptr, 1 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetPMaxElmts( precond.ptr, 0 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( precond.ptr, m_levelInterpType ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseGridMethod( precond.ptr, m_levelCoarseGridMethod ) );

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

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRHYBRIDSINGLEPHASEPOROMECHANICS_HPP_*/
