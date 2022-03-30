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
 * @file CompositionalMultiphaseHybridFVM.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASEHYBRIDFVM_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASEHYBRIDFVM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geosx
{

namespace hypre
{

namespace mgr
{

/**
 * @brief
 *
 * Labels description stored in point_marker_array
 *                         0 = pressure
 *                         1 = density
 *                       ... = ... (densities)
 * numCellCenteredLabels - 1 = density
 * numLabels - 1             = face pressure
 *
 * 3-level MGR reduction strategy inspired from CompositionalMultiphaseReservoir:
 *   - 1st level: eliminate the density associated with the volume constraint
 *   - 2nd level: eliminate the rest of the densities
 *   - 3rd level: eliminate the cell-centered pressure
 *   - The coarse grid is the interface pressure system and is solved with BoomerAMG
 */
class CompositionalMultiphaseHybridFVM : public MGRStrategyBase< 3 >
{
public:
  /**
   * @brief Constructor.
   * @param numComponentsPerField array with number of components for each field
   */
  explicit CompositionalMultiphaseHybridFVM( arrayView1d< int const > const & numComponentsPerField )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] + numComponentsPerField[1] ) )
  {
    // Level 0: eliminate the last density of the cell-centered block
    m_labels[0].resize( m_numBlocks - 1 );
    std::iota( m_labels[0].begin(), m_labels[0].begin() + m_numBlocks - 2, 0 );
    m_labels[0][m_numBlocks - 2] = m_numBlocks - 1;
    // Level 1: eliminate remaining densities of the cell-centered block
    m_labels[1].resize( 2 );
    m_labels[1][0] = 0;
    m_labels[1][1] = m_numBlocks - 1;
    // Level 2: eliminate reservoir pressure
    m_labels[2].resize( 1 );
    m_labels[2][0] = m_numBlocks - 1;

    setupLabels();

    // Level 0
    m_levelFRelaxMethod[0]     = hypre::MGRFRelaxationMethod::singleLevel; //default, i.e. Jacobi (to be confirmed)
    m_levelInterpType[0]       = hypre::MGRInterpolationType::jacobi;
    m_levelRestrictType[0]     = hypre::MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0] = hypre::MGRCoarseGridMethod::galerkin;

    // Level 1
    m_levelFRelaxMethod[1]     = hypre::MGRFRelaxationMethod::singleLevel; //default, i.e. Jacobi (to be confirmed)
    m_levelInterpType[1]       = hypre::MGRInterpolationType::jacobi;
    m_levelRestrictType[1]     = hypre::MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1] = hypre::MGRCoarseGridMethod::galerkin;

    // Level 2
    m_levelFRelaxMethod[2]     = hypre::MGRFRelaxationMethod::singleLevel; //default, i.e. Jacobi (to be confirmed)
    m_levelInterpType[2]       = hypre::MGRInterpolationType::jacobi;
    m_levelRestrictType[2]     = hypre::MGRRestrictionType::injection;
    m_levelCoarseGridMethod[2] = hypre::MGRCoarseGridMethod::galerkin;

    m_globalSmoothType = hypre::MGRGlobalSmootherType::ilu0; // TODO switch to levelGlobalSmoother when hypre allows
    m_numGlobalSmoothSweeps = 0; // No global smoother
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
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetGlobalsmoothType( precond.ptr, toUnderlying( m_globalSmoothType ) ) );
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

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASEHYBRIDFVM_HPP_*/
