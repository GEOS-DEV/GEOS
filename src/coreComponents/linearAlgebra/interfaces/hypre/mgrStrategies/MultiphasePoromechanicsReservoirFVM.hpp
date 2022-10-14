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
 * @file MultiphasePoromechanics.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRMULTIPHASEPOROMECHANICSRESERVOIRFVM_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRMULTIPHASEPOROMECHANICSRESERVOIRFVM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geosx
{

namespace hypre
{

namespace mgr
{

/**
 * @brief MultiphasePoromechanicsReservoirFVM strategy.
 *
 * Labels description stored in point_marker_array
 *   - dofLabel: 0                = nodal displacement, x-component
 *   - dofLabel: 1                = nodal displacement, y-component
 *   - dofLabel: 2                = nodal displacement, z-component
 *   - dofLabel: 3                = pressure
 *   - dofLabel: 4                = density
 *             ...                = densities
 *   - dofLabel: numResLabels + 2 = density
 *   - dofLabel: numResLabels + 3 = well pressure
 *   - dofLabel: numResLabels + 4 = well density
 *             ...                = well densities
 *   - dofLabel: numResLabels + numWellLabels + 1 = well density
 *   - dofLabel: numResLabels + numWellLabels + 2 = well rate
 *
 * 4-level MGR reduction strategy
 *   - 1st level: eliminate displacements (0,1,2)
 *   - 2nd level: eliminate the well block
 *   - 3nd level: eliminate the reservoir density associated with the volume constraint
 *   - 4nd level: eliminate the other reservoir densities
 *   - The coarse grid is solved with BoomerAMG.
 *
 */
class MultiphasePoromechanicsReservoirFVM : public MGRStrategyBase< 4 >
{
public:
  /**
   * @brief Constructor.
   * @param numComponentsPerField array with number of components for each field
   */
  explicit MultiphasePoromechanicsReservoirFVM( arrayView1d< int const > const & numComponentsPerField )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] + numComponentsPerField[1] + numComponentsPerField[2] ) )
  {
    // Level 0: eliminate displacement degrees of freedom
    m_labels[0].resize( m_numBlocks - 3 );
    std::iota( m_labels[0].begin(), m_labels[0].end(), 3 );

    HYPRE_Int const numResLabels = LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[1] );

    // Level 1: eliminate the well block
    m_labels[1].resize( numResLabels );
    std::iota( m_labels[1].begin(), m_labels[1].end(), 3 );
    // Level 2: eliminate reservoir last density which corresponds to the volume constraint equation
    m_labels[2].resize( numResLabels - 1 );
    std::iota( m_labels[2].begin(), m_labels[2].end(), 3 );
    // Level 3: eliminate the remaining reservoir densities
    m_labels[3].push_back( 3 );

    setupLabels();

    // Level 0
    m_levelFRelaxMethod[0]     = MGRFRelaxationMethod::amgVCycle;
    m_levelInterpType[0]       = MGRInterpolationType::jacobi;
    m_levelRestrictType[0]     = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0] = MGRCoarseGridMethod::nonGalerkin;

    // Level 1
    m_levelFRelaxMethod[1]     = MGRFRelaxationMethod::singleLevel;
    m_levelInterpType[1]       = MGRInterpolationType::blockJacobi;
    m_levelRestrictType[1]     = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1] = MGRCoarseGridMethod::galerkin;

    // Level 2
    m_levelFRelaxMethod[2]     = MGRFRelaxationMethod::singleLevel; //default, i.e. Jacobi
    m_levelInterpType[2]       = MGRInterpolationType::jacobi;
    m_levelRestrictType[2]     = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[2] = MGRCoarseGridMethod::galerkin;

    // Level 3
    m_levelFRelaxMethod[3]     = MGRFRelaxationMethod::singleLevel; //default, i.e. Jacobi
    m_levelInterpType[3]       = MGRInterpolationType::injection;
    m_levelRestrictType[3]     = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[3] = MGRCoarseGridMethod::cprLikeBlockDiag;

    // ILU smoothing for the system made of pressure and densities (except the last one)
    m_levelSmoothType[3] = 16;
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
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCpointsByPointMarkerArray( precond.ptr,
                                                                  m_numBlocks, numLevels,
                                                                  m_numLabels, m_ptrLabels,
                                                                  mgrData.pointMarkers.data() ) );

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelFRelaxMethod( precond.ptr, toUnderlyingPtr( m_levelFRelaxMethod ) ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( precond.ptr, toUnderlyingPtr( m_levelInterpType ) ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelRestrictType( precond.ptr, toUnderlyingPtr( m_levelRestrictType ) ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseGridMethod( precond.ptr, toUnderlyingPtr( m_levelCoarseGridMethod ) ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( precond.ptr, 1 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetPMaxElmts( precond.ptr, 0 ));

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelSmoothType( precond.ptr, m_levelSmoothType ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelSmoothIters( precond.ptr, m_levelSmoothIters ) );

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

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRMULTIPHASEPOROMECHANICSRESERVOIRFVM_HPP_*/
