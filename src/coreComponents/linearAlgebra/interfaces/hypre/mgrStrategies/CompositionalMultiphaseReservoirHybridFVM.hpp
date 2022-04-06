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
 * @file CompositionalMultiphaseReservoirHybridFVM.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASERESERVOIRHYBRIDFVM_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASERESERVOIRHYBRIDFVM_HPP_

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
 *          numResLabels - 1 = face pressure
 *              numResLabels = well pressure
 *                         1 = well density
 *                       ... = ... (densities)
 *             numLabels - 1 = well rate
 *
 * 4-level MGR reduction strategy inspired from CompositionalMultiphaseReservoir
 * 1st level: eliminate the density associated with the volume constraint
 * 2nd level: eliminate the rest of the densities
 * 3rd level: eliminate the cell-centered pressure
 * 4th level: eliminate the face-centered pressure
 * The coarse grid is the well system and is solved with a direct solver
 */
class CompositionalMultiphaseReservoirHybridFVM : public MGRStrategyBase< 4 >
{
public:
  /**
   * @brief Constructor.
   * @param numComponentsPerField array with number of components for each field
   */
  explicit CompositionalMultiphaseReservoirHybridFVM( arrayView1d< int const > const & numComponentsPerField )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] + numComponentsPerField[1] + numComponentsPerField[2] ) )
  {
    HYPRE_Int const numResCellCenteredLabels = LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] );
    HYPRE_Int const numResFaceCenteredLabels = LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[1] );
    HYPRE_Int const numResLabels = numResCellCenteredLabels + numResFaceCenteredLabels;

    // Level 0: eliminate the last density of the reservoir block
    m_labels[0].resize( m_numBlocks - 1 );
    std::iota( m_labels[0].begin(), m_labels[0].begin() + numResCellCenteredLabels - 1, 0 );
    std::iota( m_labels[0].begin() + numResCellCenteredLabels - 1, m_labels[0].end(), numResCellCenteredLabels );
    // Level 1: eliminate remaining densities of the reservoir block
    m_labels[1].resize( m_numBlocks - numResCellCenteredLabels + 1 );
    m_labels[1][0] = 0;
    std::iota( m_labels[1].begin() + 1, m_labels[1].end(), numResCellCenteredLabels );
    // Level 2: eliminate reservoir cell centered pressure
    m_labels[2].resize( m_numBlocks - numResCellCenteredLabels );
    std::iota( m_labels[2].begin(), m_labels[2].end(), numResCellCenteredLabels );
    // Level 3: eliminate reservoir face pressure
    m_labels[3].resize( m_numBlocks - numResLabels );
    std::iota( m_labels[3].begin(), m_labels[3].end(), numResLabels );

    setupLabels();

    m_levelFRelaxMethod[0] = 0; // Jacobi
    m_levelFRelaxMethod[1] = 0; // Jacobi
    m_levelFRelaxMethod[2] = 0; // Jacobi
    m_levelFRelaxMethod[3] = 2; // AMG V-cycle

    m_levelCoarseGridMethod[0] = 1;
    m_levelCoarseGridMethod[1] = 1;
    m_levelCoarseGridMethod[2] = 1;
    m_levelCoarseGridMethod[3] = 0;

    m_levelInterpType[0] = 2;
    m_levelInterpType[1] = 2;
    m_levelInterpType[2] = 2;
    m_levelInterpType[3] = 2;

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

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelFRelaxMethod( precond.ptr, m_levelFRelaxMethod ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( precond.ptr, 1 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetTruncateCoarseGridThreshold( precond.ptr, 1e-14 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetPMaxElmts( precond.ptr, 13 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( precond.ptr, m_levelInterpType ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseGridMethod( precond.ptr, m_levelCoarseGridMethod ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetMaxGlobalSmoothIters( precond.ptr, m_numGlobalSmoothSweeps ) );

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRDirectSolverCreate( &mgrData.coarseSolver.ptr ) );

    mgrData.coarseSolver.setup = HYPRE_MGRDirectSolverSetup;
    mgrData.coarseSolver.solve = HYPRE_MGRDirectSolverSolve;
    mgrData.coarseSolver.destroy = HYPRE_MGRDirectSolverDestroy;
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASERESERVOIRHYBRIDFVM_HPP_*/
