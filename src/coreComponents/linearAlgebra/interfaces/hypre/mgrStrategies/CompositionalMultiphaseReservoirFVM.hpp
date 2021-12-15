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

    // Level 0: eliminate the well block
    m_labels[0].resize( numResLabels );
    std::iota( m_labels[0].begin(), m_labels[0].end(), 0 );
    // Level 1: eliminate the last density of the reservoir block
    m_labels[1].resize( numResLabels - 1 );
    std::iota( m_labels[1].begin(), m_labels[1].end(), 0 );
    // Level 2: eliminate the rest of the densities
    m_labels[2].push_back( 0 );

    setupLabels();

    m_levelCoarseGridMethod[0] = 0; // Standard Galerkin
    m_levelCoarseGridMethod[1] = 0; // Standard Galerkin
    m_levelCoarseGridMethod[2] = 3; // Quasi-IMPES reduction

    m_levelInterpType[0] = 12; // Exact well block elimination with block-Jacobi
    m_levelInterpType[1] = 2; // Diagonal scaling (Jacobi)
    m_levelInterpType[2] = 0; // Injection

    // Set global smoothing type/iterations at each level
    // Use block-GS for the condensed system
    m_levelSmoothType[2] = 1;
    m_levelSmoothIters[2] = 1;
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

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( precond.ptr, 1 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetRelaxType( precond.ptr, 0 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNumRelaxSweeps( precond.ptr, 1 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( precond.ptr, m_levelInterpType ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelSmoothType( precond.ptr, m_levelSmoothType ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelSmoothIters( precond.ptr, m_levelSmoothIters ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseGridMethod( precond.ptr, m_levelCoarseGridMethod ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetTruncateCoarseGridThreshold( precond.ptr, 1e-20 )); // Low tolerance to remove only zeros
#ifdef GEOSX_USE_HYPRE_CUDA
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetRelaxType( precond.ptr, 18 )); // l1-Jacobi
#endif

    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &mgrData.coarseSolver.ptr ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( mgrData.coarseSolver.ptr, 0 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( mgrData.coarseSolver.ptr, 1 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetAggNumLevels( mgrData.coarseSolver.ptr, 1 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( mgrData.coarseSolver.ptr, 0.0 ) );
#ifdef GEOSX_USE_HYPRE_CUDA
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCoarsenType( mgrData.coarseSolver.ptr, 8 ) ); // PMIS
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( mgrData.coarseSolver.ptr, 18 ) ); // l1-Jacobi
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( mgrData.coarseSolver.ptr, 2 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxRowSum( mgrData.coarseSolver.ptr, 1.0 ) );
#else
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxOrder( mgrData.coarseSolver.ptr, 1 ) );
#endif

    mgrData.coarseSolver.setup = HYPRE_BoomerAMGSetup;
    mgrData.coarseSolver.solve = HYPRE_BoomerAMGSolve;
    mgrData.coarseSolver.destroy = HYPRE_BoomerAMGDestroy;
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASERESERVOIRFVM_HPP_*/
