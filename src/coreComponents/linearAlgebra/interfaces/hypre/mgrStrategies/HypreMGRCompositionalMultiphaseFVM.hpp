/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HypreMGRCompositionalMultiphaseFVM.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASEFVM_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASEFVM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGRStrategies.hpp"

namespace geosx
{

namespace hypre
{

namespace mgr
{

/**
 * @brief CompositionalMultiphaseFVM strategy.
 *
 * Labels description stored in point_marker_array
 *               0 = pressure
 *               1 = density
 *             ... = densities
 *   numLabels - 1 = density
 *
 * 2-level MGR reduction strategy which seems to work well for 2 components:
 *   - 1st level: eliminate the reservoir density associated with the volume constraint
 *   - 2nd level: eliminate the pressure
 *   - The coarse grid solved with ILU(0).
 *
 * @todo:
 *   - Experiment with block Jacobi for F-relaxation/interpolation of the reservoir densities
 *   - Explore ways to reduce onto the pressure variable and use AMG for coarse-grid solve
 */
class CompositionalMultiphaseFVM : public MGRStrategyBase< 2 >
{
public:
  /**
   * @brief Constructor.
   * @param numComponentsPerField array with number of components for each field
   */
  explicit CompositionalMultiphaseFVM( arrayView1d< localIndex const > const & numComponentsPerField )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] ) )
  {
    // Level 0: eliminate last density which corresponds to the volume constraint equation
    m_labels[0].resize( m_numBlocks - 1 );
    std::iota( m_labels[0].begin(), m_labels[0].end(), 0 );
    // Level 1: eliminate pressure
    m_labels[1].resize( m_numBlocks - 2 );
    std::iota( m_labels[1].begin(), m_labels[1].end(), 1 );

    setupLabels();

    m_levelFRelaxMethod[0] = 0; // Jacobi
    m_levelFRelaxMethod[1] = 2; // AMG V-cycle

    m_globalSmoothType = 16; // ILU(0)
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
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelFRelaxMethod( precond.ptr, m_levelFRelaxMethod ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( precond.ptr, 1 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetGlobalsmoothType( precond.ptr, m_globalSmoothType ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetMaxGlobalsmoothIters( precond.ptr, m_numGlobalSmoothSweeps ) );

    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUCreate( &mgrData.coarseSolver.ptr ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetType( mgrData.coarseSolver.ptr, 0 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetLevelOfFill( mgrData.coarseSolver.ptr, 0 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetMaxIter( mgrData.coarseSolver.ptr, 1 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetTol( mgrData.coarseSolver.ptr, 0.0 ) );

    mgrData.coarseSolver.setup = HYPRE_ILUSetup;
    mgrData.coarseSolver.solve = HYPRE_ILUSolve;
    mgrData.coarseSolver.destroy = HYPRE_ILUDestroy;
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASEFVM_HPP_*/
