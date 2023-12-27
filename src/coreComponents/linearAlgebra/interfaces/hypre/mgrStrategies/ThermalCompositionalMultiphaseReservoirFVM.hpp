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

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRTHERMALCOMPOSITIONALMULTIPHASERESERVOIRFVM_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRTHERMALCOMPOSITIONALMULTIPHASERESERVOIRFVM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

class ThermalCompositionalMultiphaseReservoirFVM : public MGRStrategyBase< 3 >
{
public:

  explicit ThermalCompositionalMultiphaseReservoirFVM( arrayView1d< int const > const & numComponentsPerField )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] + numComponentsPerField[1] ) )
  {
    HYPRE_Int const numResLabels = LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] );

    // Level 0: eliminate the well block
    m_labels[0].resize( numResLabels );
    std::iota( m_labels[0].begin(), m_labels[0].end(), 0 );
    // Level 1: eliminate the last density of the reservoir block
    m_labels[1].resize( numResLabels - 2 );
    std::iota( m_labels[1].begin(), m_labels[1].end(), 0 );
    m_labels[1].push_back( numResLabels - 1 );
    // Level 2: eliminate the rest of the densities
    m_labels[2].push_back( 0 );
    m_labels[2].push_back( numResLabels - 1 );

    setupLabels();

    m_levelFRelaxType[0] = MGRFRelaxationType::gsElimWInverse; // gaussian elimination for the well block
    m_levelInterpType[0] = MGRInterpolationType::blockJacobi;
    m_levelRestrictType[0] = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0] = MGRCoarseGridMethod::galerkin;

    m_levelInterpType[1] = MGRInterpolationType::jacobi; // Diagonal scaling (Jacobi)
    m_levelRestrictType[1] = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1] = MGRCoarseGridMethod::galerkin; // Standard Galerkin

    m_levelInterpType[2] = MGRInterpolationType::injection; // Injection
    m_levelRestrictType[2] = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[2] = MGRCoarseGridMethod::cprLikeBlockDiag; // Non-Galerkin Quasi-IMPES CPR

    // Block Jacobi for the system made of pressure, temperature, and densities
    m_levelGlobalSmootherType[1] = MGRGlobalSmootherType::blockGaussSeidel;
    m_levelGlobalSmootherIters[1] = 1;
#if GEOS_USE_HYPRE_DEVICE != GEOS_USE_HYPRE_HIP
    // ILU smoothing for the system made of pressure, temperature, and densities (except the last one)
    m_levelGlobalSmootherType[2] = MGRGlobalSmootherType::ilu0;
    m_levelGlobalSmootherIters[2] = 1;
#endif
  }

  void setup( LinearSolverParameters::MGR const & mgrParams,
              HyprePrecWrapper & precond,
              HypreMGRData & mgrData )
  {

    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetCpointsByPointMarkerArray( precond.ptr,
                                                                 m_numBlocks, numLevels,
                                                                 m_numLabels, m_ptrLabels,
                                                                 mgrData.pointMarkers.data() ) );

    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( precond.ptr, toUnderlyingPtr( m_levelInterpType ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelRestrictType( precond.ptr, toUnderlyingPtr( m_levelRestrictType ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseGridMethod( precond.ptr, toUnderlyingPtr( m_levelCoarseGridMethod ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( precond.ptr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelSmoothType( precond.ptr, toUnderlyingPtr( m_levelGlobalSmootherType ) ));
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelSmoothIters( precond.ptr, m_levelGlobalSmootherIters ));
    // if the wells are shut, using Gaussian elimination as F-relaxation for the well block is an overkill
    // in that case, we just use Jacobi
    if( mgrParams.areWellsShut )
    {
      m_levelFRelaxType[0] = MGRFRelaxationType::jacobi;
    }

    // Note: uncomment HYPRE_MGRSetLevelFRelaxMethod and comment HYPRE_MGRSetRelaxType breaks the recipe, this requires further
    // investigation
    //GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelFRelaxMethod( precond.ptr, toUnderlyingPtr( m_levelFRelaxMethod ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelFRelaxType( precond.ptr, toUnderlyingPtr( m_levelFRelaxType ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetNumRelaxSweeps( precond.ptr, 1 ) );

#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetRelaxType( precond.ptr, getAMGRelaxationType( LinearSolverParameters::AMG::SmootherType::l1jacobi ) ) ); // l1-Jacobi
#endif

    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &mgrData.coarseSolver.ptr ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( mgrData.coarseSolver.ptr, 0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( mgrData.coarseSolver.ptr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetAggNumLevels( mgrData.coarseSolver.ptr, 1 ) );//??
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( mgrData.coarseSolver.ptr, 0.0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumFunctions( mgrData.coarseSolver.ptr, 2 ) ); // pressure and temperature (CPTR)
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCoarsenType( mgrData.coarseSolver.ptr, toUnderlying( AMGCoarseningType::PMIS ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( mgrData.coarseSolver.ptr, getAMGRelaxationType( LinearSolverParameters::AMG::SmootherType::l1jacobi ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( mgrData.coarseSolver.ptr, 2 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxRowSum( mgrData.coarseSolver.ptr, 1.0 ) );
#else
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxOrder( mgrData.coarseSolver.ptr, 1 ));
#endif

    mgrData.coarseSolver.setup = HYPRE_BoomerAMGSetup;
    mgrData.coarseSolver.solve = HYPRE_BoomerAMGSolve;
    mgrData.coarseSolver.destroy = HYPRE_BoomerAMGDestroy;

  }

};

}// namespace mgr
}// namespace hypre
}// namespace geos

#endif
