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
 * @file SinglePhaseReservoirHybridFVM.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASERESERVOIRHYBRIDFVM_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASERESERVOIRHYBRIDFVM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geosx
{

namespace hypre
{

namespace mgr
{

/**
 * @brief SinglePhaseReservoirHybridFVM strategy.
 *
 * Labels description stored in point_marker_array
 *   dofLabel: 0 = cell-centered reservoir pressure
 *   dofLabel: 1 = face-centered reservoir pressure
 *   dofLabel: 2 = well pressure
 *   dofLabel: 3 = well rate
 *
 * Ingredients
 *
 * 2-level MGR reduction strategy
 * 1st level: eliminate the cell-centered reservoir pressure
 * 2nd level: eliminate the face-centered reservoir pressure (Lagrange multiplier)
 * The coarse grid is the well vars
 * The coarse grid solved with a direct solver
 */
class SinglePhaseReservoirHybridFVM : public MGRStrategyBase< 2 >
{
public:
  /**
   * @brief Constructor.
   */
  explicit SinglePhaseReservoirHybridFVM( arrayView1d< int const > const & )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( 4 ) )
  {
    // Level 0: eliminate the cell-centered pressure
    m_labels[0].push_back( 1 );
    m_labels[0].push_back( 2 );
    m_labels[0].push_back( 3 );
    // Level 1: eliminate the face-centered pressure
    m_labels[1].push_back( 2 );
    m_labels[1].push_back( 3 );

    setupLabels();

    // Level 0
    m_levelFRelaxMethod[0]     = hypre::MGRFRelaxationMethod::singleLevel; //default, i.e. Jacobi (to be confirmed)
    m_levelInterpType[0]       = hypre::MGRInterpolationType::jacobi;
    m_levelRestrictType[0]     = hypre::MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0] = hypre::MGRCoarseGridMethod::nonGalerkin;

    // Level 1
    m_levelFRelaxMethod[1]     = hypre::MGRFRelaxationMethod::amgVCycle;
    m_levelInterpType[1]       = hypre::MGRInterpolationType::jacobi;
    m_levelRestrictType[1]     = hypre::MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1] = hypre::MGRCoarseGridMethod::galerkin;

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

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelFRelaxMethod( precond.ptr, toUnderlyingPtr( m_levelFRelaxMethod ) ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( precond.ptr, toUnderlyingPtr( m_levelInterpType ) ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelRestrictType( precond.ptr, toUnderlyingPtr( m_levelRestrictType ) ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseGridMethod( precond.ptr, toUnderlyingPtr( m_levelCoarseGridMethod ) ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( precond.ptr, 1 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetMaxGlobalsmoothIters( precond.ptr, m_numGlobalSmoothSweeps ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetTruncateCoarseGridThreshold( precond.ptr, 1e-14 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetPMaxElmts( precond.ptr, 15 ));

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRDirectSolverCreate( &mgrData.coarseSolver.ptr ) );

    mgrData.coarseSolver.setup = HYPRE_MGRDirectSolverSetup;
    mgrData.coarseSolver.solve = HYPRE_MGRDirectSolverSolve;
    mgrData.coarseSolver.destroy = HYPRE_MGRDirectSolverDestroy;
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASERESERVOIRHYBRIDFVM_HPP_*/
