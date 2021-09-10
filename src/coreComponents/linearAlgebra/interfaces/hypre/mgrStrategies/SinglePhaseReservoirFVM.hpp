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
 * @file SinglePhaseReservoirFVM.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASERESERVOIRFVM_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASERESERVOIRFVM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geosx
{

namespace hypre
{

namespace mgr
{

/**
 * @brief SinglePhaseReservoirFVM strategy.
 *
 * Labels description stored in point_marker_array
 *   dofLabel: 0 = reservoir pressure (numResLabels = 1)
 *   dofLabel: 1 = well pressure
 *   dofLabel: 2 = well rate (numWellLabels = 2)
 *
 * Ingredients
 *
 * 1. F-points cell-centered pressure, C-points well vars
 * 2. F-points smoother: boomer AMG
 * 3. C-points coarse-grid/Schur complement solver: direct solver
 * 4. Global smoother: none
 */
class SinglePhaseReservoirFVM : public MGRStrategyBase< 1 >
{
public:
  /**
   * @brief Constructor.
   */
  explicit SinglePhaseReservoirFVM( arrayView1d< int const > const & )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( 3 ) )
  {
    // Level 0: eliminate the cell-centered pressure
    m_labels[0].push_back( 1 );
    m_labels[0].push_back( 2 );
    setupLabels();

    m_levelInterpType[0] = 2; // diagonal scaling
    m_levelCoarseGridMethod[0] = 0; // Galerkin coarse grid computation using RAP

    m_fRelaxMethod = 2; // AMG V-cycle
    m_numGlobalSmoothSweeps = 0;
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

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetFRelaxMethod( precond.ptr, m_fRelaxMethod ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( precond.ptr, 1 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( precond.ptr, m_levelInterpType ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseGridMethod( precond.ptr, m_levelCoarseGridMethod ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetMaxGlobalsmoothIters( precond.ptr, m_numGlobalSmoothSweeps ) );

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRDirectSolverCreate( &mgrData.coarseSolver.ptr ) );

    mgrData.coarseSolver.setup = HYPRE_MGRDirectSolverSetup;
    mgrData.coarseSolver.solve = HYPRE_MGRDirectSolverSolve;
    mgrData.coarseSolver.destroy = HYPRE_MGRDirectSolverDestroy;
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASERESERVOIRFVM_HPP_*/
