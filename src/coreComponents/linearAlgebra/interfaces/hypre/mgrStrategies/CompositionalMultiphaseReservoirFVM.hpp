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
 * @file CompositionalMultiphaseReservoirFVM.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASERESERVOIRFVM_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASERESERVOIRFVM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geosx
{

namespace hypre
{

namespace mgr
{

/**
 * @brief CompositionalMultiphaseReservoirFVM strategy.
 *
 * Labels description stored in point_marker_array
 *                0 = reservoir pressure
 *                1 = reservoir density
 *              ... = ... (reservoir densities)
 * numResLabels - 1 = reservoir density
 *     numResLabels = well pressure
 * numResLabels + 1 = well density
 *              ... = ... (well densities)
 * numResLabels + numWellLabels - 2 = well density
 * numResLabels + numWellLabels - 1 = well rate
 *
 * 3-level MGR reduction strategy which seems to work well for 2 components:
 *   - 1st level: eliminate the reservoir density associated with the volume constraint
 *   - 2nd level: eliminate the rest of the reservoir densities
 *   - 3rd level: eliminate the pressure
 *   - The coarse grid is the well block and solved with ILU(0)
 *
 * @todo:
 *   - Use block Jacobi for F-relaxation/interpolation of the reservoir densities (2nd level)
 */
class CompositionalMultiphaseReservoirFVM : public MGRStrategyBase< 3 >
{
public:

  /**
   * @brief Constructor.
   * @param numComponentsPerField array with number of components for each field
   */
  explicit CompositionalMultiphaseReservoirFVM( arrayView1d< int const > const & numComponentsPerField )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] + numComponentsPerField[1] ) )
  {
    HYPRE_Int const numResLabels = LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] );
    // Level 0: eliminate the last density of the reservoir block
    m_labels[0].resize( m_numBlocks - 1 );
    std::iota( m_labels[0].begin(), m_labels[0].begin() + numResLabels - 1, 0 );
    std::iota( m_labels[0].begin() + numResLabels - 1, m_labels[0].end(), numResLabels );
    // Level 1: eliminate remaining densities of the reservoir block
    m_labels[1].resize( m_numBlocks - numResLabels + 1 );
    m_labels[1][0] = 0;
    std::iota( m_labels[1].begin() + 1, m_labels[1].end(), numResLabels );
    // Level 2: eliminate reservoir pressure
    m_labels[2].resize( m_numBlocks - numResLabels );
    std::iota( m_labels[2].begin(), m_labels[2].end(), numResLabels );

    setupLabels();

    // Level 0
    m_levelFRelaxMethod[0]     = hypre::MGRFRelaxationMethod::singleLevel; //default, i.e. Jacobi (to be confirmed)
    m_levelInterpType[0]       = hypre::MGRInterpolationType::jacobi;
    m_levelRestrictType[0]     = hypre::MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0] = hypre::MGRCoarseGridMethod::nonGalerkin;

    // Level 1
    m_levelFRelaxMethod[1]     = hypre::MGRFRelaxationMethod::singleLevel; //default, i.e. Jacobi (to be confirmed)
    m_levelInterpType[1]       = hypre::MGRInterpolationType::jacobi;
    m_levelRestrictType[1]     = hypre::MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1] = hypre::MGRCoarseGridMethod::nonGalerkin;

    // Level 2
    m_levelFRelaxMethod[2]     = hypre::MGRFRelaxationMethod::amgVCycle;
    m_levelInterpType[2]       = hypre::MGRInterpolationType::jacobi;
    m_levelRestrictType[2]     = hypre::MGRRestrictionType::injection;
    m_levelCoarseGridMethod[2] = hypre::MGRCoarseGridMethod::galerkin;

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
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetTruncateCoarseGridThreshold( precond.ptr, 1e-14 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetPMaxElmts( precond.ptr, 15 ));
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

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASERESERVOIRFVM_HPP_*/
