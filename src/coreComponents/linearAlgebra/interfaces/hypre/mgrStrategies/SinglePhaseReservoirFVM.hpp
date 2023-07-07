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

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASERESERVOIRFVM_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASERESERVOIRFVM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
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
 * 1. F-points well vars, C-points cell-centered pressure
 * 2. F-points smoother: jacobi (soon, direct solver)
 * 3. C-points coarse-grid/Schur complement solver: BoomerAMG
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
    // Level 0: eliminate the well variables, and just keep the cell-centered pressures
    m_labels[0].push_back( 0 );

    setupLabels();

    // Level 0
    m_levelFRelaxMethod[0]     = MGRFRelaxationMethod::singleLevel;
    m_levelFRelaxType[0]       = MGRFRelaxationType::gsElimWInverse; //default, i.e. Jacobi
    m_levelInterpType[0]       = MGRInterpolationType::blockJacobi;
    m_levelRestrictType[0]     = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0] = MGRCoarseGridMethod::galerkin;

    m_levelSmoothIters[0] = 0;
  }

  /**
   * @brief Setup the MGR strategy.
   * @param mgrParams parameters for the configuration of the MGR recipe
   * @param precond preconditioner wrapper
   * @param mgrData auxiliary MGR data
   */
  void setup( LinearSolverParameters::MGR const & mgrParams,
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
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelSmoothType( precond.ptr, m_levelSmoothType ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelSmoothIters( precond.ptr, m_levelSmoothIters ) );

    // if the wells are shut, using Gaussian elimination as F-relaxation for the well block is an overkill
    // in that case, we just use Jacobi
    if( mgrParams.areWellsShut )
    {
      m_levelFRelaxType[0] = MGRFRelaxationType::jacobi;
    }

    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelFRelaxType( precond.ptr, toUnderlyingPtr( m_levelFRelaxType ) ));
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetNumRelaxSweeps( precond.ptr, 1 ));

#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetRelaxType( precond.ptr, getAMGRelaxationType( LinearSolverParameters::AMG::SmootherType::l1jacobi ) ) ); // l1-Jacobi
#endif

    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &mgrData.coarseSolver.ptr ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( mgrData.coarseSolver.ptr, 0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( mgrData.coarseSolver.ptr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( mgrData.coarseSolver.ptr, 0.0 ) );
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCoarsenType( mgrData.coarseSolver.ptr, toUnderlying( AMGCoarseningType::PMIS ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( mgrData.coarseSolver.ptr, getAMGRelaxationType( LinearSolverParameters::AMG::SmootherType::l1jacobi ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( mgrData.coarseSolver.ptr, 2 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxRowSum( mgrData.coarseSolver.ptr, 1.0 ) );
#else
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxOrder( mgrData.coarseSolver.ptr, 1 ) );
#endif

    mgrData.coarseSolver.setup = HYPRE_BoomerAMGSetup;
    mgrData.coarseSolver.solve = HYPRE_BoomerAMGSolve;
    mgrData.coarseSolver.destroy = HYPRE_BoomerAMGDestroy;
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASERESERVOIRFVM_HPP_*/
