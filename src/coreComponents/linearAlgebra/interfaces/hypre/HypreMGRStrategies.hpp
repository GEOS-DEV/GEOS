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
 * @file HypreMGRStrategies.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRSTRATEGIES_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRSTRATEGIES_HPP_

#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"

#include <_hypre_utilities.h>

namespace geosx
{

/**
 * @brief Container for hypre preconditioner auxiliary data for MGR.
 */
struct HypreMGRData
{
  array1d< HYPRE_Int > pointMarkers;  ///< array1d of unique tags for local degrees of freedom
  HyprePrecWrapper coarseSolver;      ///< MGR coarse solver pointer and functions
  HyprePrecWrapper mechSolver;        ///< MGR mechanics fine solver pointer and functions
};

namespace hypre
{

namespace mgr
{

/**
 * @brief Helper to simplify MGR setup
 * @tparam STRATEGY strategy class (one of structs defined below)
 * @param params MGR parameters
 * @param numComponentsPerField array containing number of components in each field
 * @param precond the preconditioner wrapper
 * @param mgrData additional MGR data struct with marker array populated
 */
template< typename STRATEGY >
void setStrategy( LinearSolverParameters::MGR const & params,
                  arrayView1d< localIndex const > const & numComponentsPerField,
                  HyprePrecWrapper & precond,
                  HypreMGRData & mgrData )
{
  STRATEGY strategy( numComponentsPerField );
  strategy.setup( params, precond, mgrData );
}

/**
 * @brief Helper struct for strategies that provides some basic parameter arrays needed by MGR.
 * @tparam NLEVEL number of reduction levels (not including the coarsest level)
 */
template< int NLEVEL >
class MGRStrategyBase
{
public:

  static constexpr HYPRE_Int numLevels = NLEVEL;       ///< Number of levels

protected:

  HYPRE_Int m_numBlocks{ 0 };                          ///< Number of different matrix blocks treated separately

  std::vector< HYPRE_Int > m_labels[numLevels]{};      ///< Dof labels kept at each level
  HYPRE_Int m_numLabels[numLevels]{ -1 };              ///< Number of dof labels kept
  HYPRE_Int * m_ptrLabels[numLevels]{ nullptr };       ///< Pointers to each level's labels, as consumed by MGR

  HYPRE_Int m_levelFRelaxMethod[numLevels]{ -1 };      ///< F-relaxation method for each level
  HYPRE_Int m_levelInterpType[numLevels]{ -1 };        ///< Interpolation type for each level
  HYPRE_Int m_levelRestrictType[numLevels]{ -1 };      ///< Restriction type for each level
  HYPRE_Int m_levelCoarseGridMethod[numLevels]{ -1 };  ///< Coarse grid method for each level

  HYPRE_Int m_numRestrictSweeps{ -1 };                 ///< Number of restrict sweeps
  HYPRE_Int m_numInterpSweeps{ -1 };                   ///< Number of interpolation sweeps

  HYPRE_Int m_fRelaxMethod{ -1 };                      ///< F-relaxation method
  HYPRE_Int m_relaxType{ -1 };                         ///< F-relaxation type
  HYPRE_Int m_numRelaxSweeps{ -1 };                    ///< F-relaxation number of sweeps

  HYPRE_Int m_globalSmoothType{ -1 };                  ///< Global smoothing type
  HYPRE_Int m_numGlobalSmoothSweeps{ -1 };             ///< Global smoothing number of iterations

  /**
   * @brief Constructor.
   * @param numBlocks number of blocks
   */
  explicit MGRStrategyBase( HYPRE_Int const numBlocks )
    : m_numBlocks( numBlocks )
  {}

  /**
   * @brief Call this after populating lv_cindexes.
   */
  void setupLabels()
  {
    for( HYPRE_Int i = 0; i < numLevels; ++i )
    {
      m_numLabels[i] = m_labels[i].size();
      m_ptrLabels[i] = m_labels[i].data();
    }
  }
};

/**
 * @brief SinglePhasePoromechanics strategy.
 *
 * dofLabel: 0 = displacement, x-component
 * dofLabel: 1 = displacement, y-component
 * dofLabel: 2 = displacement, z-component
 * dofLabel: 3 = pressure
 *
 * Ingredients:
 * 1. F-points displacement (0,1,2), C-points pressure (3)
 * 2. F-points smoother: AMG, single V-cycle, separate displacement components
 * 3. C-points coarse-grid/Schur complement solver: boomer AMG
 * 4. Global smoother: none
 */
class SinglePhasePoromechanics : public MGRStrategyBase< 1 >
{
public:

  /**
   * @brief Constructor.
   */
  explicit SinglePhasePoromechanics( arrayView1d< localIndex const > const & )
    : MGRStrategyBase( 4 )
  {
    m_labels[0].push_back( 3 );
    setupLabels();

    m_levelInterpType[0] = 2; // diagonal scaling (Jacobi)
    m_levelCoarseGridMethod[0] = 1; // diagonal sparsification

    m_fRelaxMethod = 2; // AMG V-cycle
    m_numGlobalSmoothSweeps = 0;
  }

  /**
   * @brief Setup the MGR strategy.
   * @param params MGR parameters
   * @param precond preconditioner wrapper
   * @param mgrData auxiliary MGR data
   */
  void setup( LinearSolverParameters::MGR const & params,
              HyprePrecWrapper & precond,
              HypreMGRData & mgrData )
  {
    if( params.strategy == LinearSolverParameters::MGR::StrategyType::hydrofracture )
    {
      m_levelCoarseGridMethod[0] = 0; // Galerkin coarse grid computation using RAP
    }

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCpointsByPointMarkerArray( precond.ptr,
                                                                  m_numBlocks, numLevels,
                                                                  m_numLabels, m_ptrLabels,
                                                                  mgrData.pointMarkers.data() ) );

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetFRelaxMethod( precond.ptr, m_fRelaxMethod ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( precond.ptr, 1 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetPMaxElmts( precond.ptr, 0 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( precond.ptr, m_levelInterpType ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseGridMethod( precond.ptr, m_levelCoarseGridMethod ) );
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

/**
 * @brief HybridSinglePhasePoromechanics strategy.
 *
 * Labels description stored in point_marker_array
 *   - dofLabel: 0 = nodal displacement, x-component
 *   - dofLabel: 1 = nodal displacement, y-component
 *   - dofLabel: 2 = nodal displacement, z-component
 *   - dofLabel: 3 = cell pressure
 *   - dofLabel: 4 = face pressure
 *
 * 2-level MGR reduction:
 *   - 1st level: eliminate displacements (0,1,2)
 *   - 2nd level: eliminate cell pressure (3)
 *   - The coarse grid solved with boomer AMG
 *
 */
class HybridSinglePhasePoromechanics : public MGRStrategyBase< 2 >
{
public:

  /**
   * @brief Constructor.
   */
  explicit HybridSinglePhasePoromechanics( arrayView1d< localIndex const > const & )
    : MGRStrategyBase( 5 )
  {
    // Level 0: eliminate displacement degrees of freedom
    m_labels[0].push_back( 3 );
    m_labels[0].push_back( 4 );
    // Level 1: eliminate cell pressure degrees of freedom
    m_labels[1].push_back( 4 );

    setupLabels();

    m_levelFRelaxMethod[0] = 2;  // AMG V-cycle
    m_levelFRelaxMethod[1] = 18; // l1-Jacobi

    m_levelInterpType[0] = 2;       // diagonal scaling (Jacobi)
    m_levelCoarseGridMethod[0] = 1; // diagonal sparsification
    m_levelInterpType[1] = 2;       // diagonal scaling (Jacobi)
    m_levelCoarseGridMethod[1] = 0; // Galerkin coarse grid computation using RAP
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
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetPMaxElmts( precond.ptr, 0 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( precond.ptr, m_levelInterpType ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseGridMethod( precond.ptr, m_levelCoarseGridMethod ) );

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

/**
 * @brief CompositionalMultiphaseReservoir strategy.
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
class CompositionalMultiphaseReservoir : public MGRStrategyBase< 3 >
{
public:

  /**
   * @brief Constructor.
   * @param numComponentsPerField array with number of components for each field
   */
  explicit CompositionalMultiphaseReservoir( arrayView1d< localIndex const > const & numComponentsPerField )
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

    m_levelFRelaxMethod[0] = 0; // Jacobi
    m_levelFRelaxMethod[1] = 0; // Jacobi
    m_levelFRelaxMethod[2] = 2; // AMG V-cycle

    m_levelCoarseGridMethod[0] = 1;
    m_levelCoarseGridMethod[1] = 1;
    m_levelCoarseGridMethod[2] = 0;

    m_levelInterpType[0] = 2;
    m_levelInterpType[1] = 2;
    m_levelInterpType[2] = 2;

    m_globalSmoothType = 16; // ILU(0)
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
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetPMaxElmts( precond.ptr, 15 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( precond.ptr, m_levelInterpType ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseGridMethod( precond.ptr, m_levelCoarseGridMethod ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetGlobalsmoothType( precond.ptr, m_globalSmoothType ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetMaxGlobalsmoothIters( precond.ptr, m_numGlobalSmoothSweeps ) );

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRDirectSolverCreate( &mgrData.coarseSolver.ptr ) );

    mgrData.coarseSolver.setup = HYPRE_MGRDirectSolverSetup;
    mgrData.coarseSolver.solve = HYPRE_MGRDirectSolverSolve;
    mgrData.coarseSolver.destroy = HYPRE_MGRDirectSolverDestroy;
  }
};

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
  explicit CompositionalMultiphaseHybridFVM( arrayView1d< localIndex const > const & numComponentsPerField )
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

    m_levelFRelaxMethod[0] = 0; // Jacobi
    m_levelFRelaxMethod[1] = 0; // Jacobi
    m_levelFRelaxMethod[2] = 0; // Jacobi

    m_levelInterpType[0] = 2;
    m_levelInterpType[1] = 2;
    m_levelInterpType[2] = 2;

    m_globalSmoothType = 16; // ILU(0)
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
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( precond.ptr, m_levelInterpType ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetGlobalsmoothType( precond.ptr, m_globalSmoothType ) );
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

/**
 * @brief LagrangianContactMechanics strategy
 *
 * Contact mechanics with face-centered lagrangian multipliers
 *
 * dofLabel: 0 = displacement, x-component
 * dofLabel: 1 = displacement, y-component
 * dofLabel: 2 = displacement, z-component
 * dofLabel: 3 = pressure
 *
 * Ingredients:
 * 1. F-points displacement (0,1,2), C-points pressure (3)
 * 2. F-points smoother: AMG, single V-cycle, separate displacemente components
 * 3. C-points coarse-grid/Schur complement solver: boomer AMG
 * 4. Global smoother: none
 *
 * @todo: (Sergey) I changed this from contiguous block setup to an equivalent marker array setup
 *        to make it a little easier (avoid passing DofManager, etc.). Will need to come up with a
 *        way to support both ways equally well. Maybe pass block offsets as constructor input.
 */
class LagrangianContactMechanics : public MGRStrategyBase< 1 >
{
public:

  /**
   * @brief Constructor.
   */
  explicit LagrangianContactMechanics( arrayView1d< localIndex const > const & )
    : MGRStrategyBase( 2 )
  {
    // Level 0: all three displacements kept
    m_labels[0] = { 0, 1, 2 };
    setupLabels();

    m_levelFRelaxMethod[0] = 0; // Jacobi
    m_levelInterpType[0] = 2; // diagonal scaling (Jacobi)
    m_levelRestrictType[0] = 0; // injection
    m_levelCoarseGridMethod[0] = 0; // all

    m_numRestrictSweeps = 0;
    m_numInterpSweeps = 0;
    m_relaxType = 0;
    m_numRelaxSweeps = 0; // skip F-relaxation
    m_globalSmoothType = 0; // block Jacobi
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

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelRestrictType( precond.ptr, m_levelRestrictType ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNumRestrictSweeps( precond.ptr, m_numRestrictSweeps ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNumInterpSweeps( precond.ptr, m_numInterpSweeps ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNumRelaxSweeps( precond.ptr, m_numRelaxSweeps ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( precond.ptr, 1 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetRelaxType( precond.ptr, m_relaxType ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelFRelaxMethod( precond.ptr, m_levelFRelaxMethod ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( precond.ptr, m_levelInterpType ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseGridMethod( precond.ptr, m_levelCoarseGridMethod ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetGlobalsmoothType( precond.ptr, m_globalSmoothType ) );
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

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREMGRSTRATEGIES_HPP_*/
