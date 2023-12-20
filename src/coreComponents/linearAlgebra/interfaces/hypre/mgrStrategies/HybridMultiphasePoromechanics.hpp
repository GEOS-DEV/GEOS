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
 * @file HybridMultiphasePoromechanics.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRHYBRIDMULTIPHASEPOROMECHANICS_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRHYBRIDMULTIPHASEPOROMECHANICS_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief HybridMultiphasePoromechanics strategy.
 *
 * Labels description stored in point_marker_array
 *   - dofLabel: 0             = nodal displacement, x-component
 *   - dofLabel: 1             = nodal displacement, y-component
 *   - dofLabel: 2             = nodal displacement, z-component
 *   - dofLabel: 3             = pressure
 *   - dofLabel: 4             = density
 *             ...             = densities
 *   - dofLabel: numLabels - 2 = density
 *   - dofLabel: numLabels - 1 = face pressure
 *
 * 4-level MGR reduction strategy based on CompositionalMultiphaseFVM
 *   - 1st level: eliminate displacements (0,1,2)
 *   - 2nd level: eliminate the reservoir density associated with the volume constraint (numLabels - 1)
 *   - 3rd level: eliminate the other reservoir densities
 *   - 4th level: eliminate the cell-centered pressure
 *   - The coarse grid is the interface pressure system and is solved with BoomerAMG
 *
 */
class HybridMultiphasePoromechanics : public MGRStrategyBase< 4 >
{
public:
  /**
   * @brief Constructor.
   * @param numComponentsPerField array with number of components for each field
   */
  explicit HybridMultiphasePoromechanics( arrayView1d< int const > const & numComponentsPerField )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] + numComponentsPerField[1] + numComponentsPerField[2] ) )
  {
    // Level 0: eliminate displacement degrees of freedom
    m_labels[0].resize( m_numBlocks - 3 );
    std::iota( m_labels[0].begin(), m_labels[0].end(), 3 );
    // Level 1: eliminate last density which corresponds to the volume constraint equation
    m_labels[1].resize( m_numBlocks - 5 );
    std::iota( m_labels[1].begin(), m_labels[1].end(), 3 );
    m_labels[1].push_back( m_numBlocks - 1 );
    // Level 2: eliminate the remaining reservoir densities
    m_labels[2].resize( 2 );
    m_labels[2][0] = 3;
    m_labels[2][1] = m_numBlocks - 1;
    // Level 3: eliminate reservoir cell-centered pressure
    m_labels[3].resize( 1 );
    m_labels[3][0] = m_numBlocks - 1;


    setupLabels();

    // Level 0
    m_levelFRelaxMethod[0]     = MGRFRelaxationMethod::amgVCycle;
    m_levelInterpType[0]       = MGRInterpolationType::jacobi;
    m_levelRestrictType[0]     = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0] = MGRCoarseGridMethod::nonGalerkin;

    // Level 1
    m_levelFRelaxMethod[1]     = MGRFRelaxationMethod::singleLevel; //default, i.e. Jacobi
    m_levelInterpType[1]       = MGRInterpolationType::jacobi;
    m_levelRestrictType[1]     = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1] = MGRCoarseGridMethod::galerkin;

    // Level 2
    m_levelFRelaxMethod[2]     = MGRFRelaxationMethod::singleLevel; //default, i.e. Jacobi
    m_levelInterpType[2]       = MGRInterpolationType::injection;
    m_levelRestrictType[2]     = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[2] = MGRCoarseGridMethod::cprLikeBlockDiag;

    // Level 3
    m_levelFRelaxMethod[3]     = MGRFRelaxationMethod::singleLevel; //default, i.e. Jacobi
    m_levelInterpType[3]       = MGRInterpolationType::jacobi;
    m_levelRestrictType[3]     = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[3] = MGRCoarseGridMethod::galerkin;


    // ILU smoothing for the system made of pressure and densities (except the last one)
    m_levelSmoothType[3]  = 1;
    m_levelSmoothIters[3] = 1;
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
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetCpointsByPointMarkerArray( precond.ptr,
                                                                 m_numBlocks, numLevels,
                                                                 m_numLabels, m_ptrLabels,
                                                                 mgrData.pointMarkers.data() ) );

    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelFRelaxMethod( precond.ptr, toUnderlyingPtr( m_levelFRelaxMethod ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( precond.ptr, toUnderlyingPtr( m_levelInterpType ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelRestrictType( precond.ptr, toUnderlyingPtr( m_levelRestrictType ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseGridMethod( precond.ptr, toUnderlyingPtr( m_levelCoarseGridMethod ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( precond.ptr, 1 ));
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetPMaxElmts( precond.ptr, 0 ));

    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelSmoothType( precond.ptr, m_levelSmoothType ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelSmoothIters( precond.ptr, m_levelSmoothIters ) );

    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &mgrData.coarseSolver.ptr ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( mgrData.coarseSolver.ptr, 0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( mgrData.coarseSolver.ptr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( mgrData.coarseSolver.ptr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( mgrData.coarseSolver.ptr, hypre::getAMGRelaxationType( LinearSolverParameters::AMG::SmootherType::l1jacobi ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( mgrData.coarseSolver.ptr, 0.0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxOrder( mgrData.coarseSolver.ptr, 1 ) );

    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetSmoothType( mgrData.coarseSolver.ptr, 5 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetSmoothNumLevels( mgrData.coarseSolver.ptr, 8 ) );

    mgrData.coarseSolver.setup = HYPRE_BoomerAMGSetup;
    mgrData.coarseSolver.solve = HYPRE_BoomerAMGSolve;
    mgrData.coarseSolver.destroy = HYPRE_BoomerAMGDestroy;

    // Configure the BoomerAMG solver used as F-relaxation for the first level
    setMechanicsFSolver( precond, mgrData );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRMULTIPHASEPOROMECHANICS_HPP_*/
