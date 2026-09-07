/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HypreMGR.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSTRATEGIES_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSTRATEGIES_HPP_

#include "linearAlgebra/common/common.hpp"

#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#include <_hypre_utilities.h>

namespace geos
{

/**
 * @brief Container for hypre preconditioner auxiliary data for MGR.
 */
struct HypreMGRData
{
  array1d< HYPRE_Int > pointMarkers;  ///< array1d of unique tags for local degrees of freedom
  HyprePrecWrapper coarseSolver;      ///< MGR coarse solver pointer and functions
  HyprePrecWrapper mechSolver;        ///< MGR mechanics fine solver pointer and functions
  HyprePrecWrapper nestedSolver;      ///< Optional nested MGR F-relaxation wrapper
};

namespace hypre
{

namespace mgr
{

/**
 * @brief MGR settings shared by the legacy HYPRE setup and generated YAML.
 *
 * The values in this descriptor are the settings used when MGR is configured
 * as a one-step preconditioner.  Strategy-specific reduction metadata remains
 * in MGRStrategyBase, while these common settings are serialized by the
 * hypredrive adapter and applied to HYPRE handles from the same object.
 */
struct MGRParameters
{
  HYPRE_Real tolerance{ 0.0 };
  HYPRE_Int maxIterations{ 1 };
  HYPRE_Int printLevel{ 0 };
  HYPRE_Int cycleType{ 1 };
  HYPRE_Int fRelaxCycle{ 1 };
  HYPRE_Int globalSmoothCycle{ 1 };
  HYPRE_Int nonCpointsToFpoints{ 1 };
  HYPRE_Int nonGalerkinMaxElmts{ 1 };
  HYPRE_Int pMaxElmts{ 0 };
};

inline MGRParameters defaultMGRParameters()
{
  return {};
}

HYPRE_Int constexpr hydrofractureMinCoarseSize = 1000;

inline void setMGRCycleSettings( HYPRE_Solver const solver,
                                 MGRParameters const & params = defaultMGRParameters() )
{
  GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetCycleType( solver, params.cycleType ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetFRelaxCycle( solver, params.fRelaxCycle ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetGlobalSmoothCycle( solver, params.globalSmoothCycle ) );
}

/**
 * @brief The BoomerAMG options used by the MGR strategies.
 *
 * This is the single description of the nested AMG instances. The legacy
 * path applies it through the HYPRE API and the hypredrive path serializes
 * the same values to YAML.
 *
 * A negative value means that an option is not part of a particular flavor.
 */
struct BoomerAMGParameters
{
  HYPRE_Real tolerance{ 0.0 };
  HYPRE_Int maxIterations{ 1 };
  HYPRE_Int printLevel{ 0 };
  HYPRE_Int minCoarseSize{ -1 };
  HYPRE_Int maxCoarseSize{ 9 };
  HYPRE_Int smoothType{ 6 };
  HYPRE_Int smoothNumLevels{ 0 };
  HYPRE_Int smoothNumSweeps{ 1 };
  HYPRE_Int smoothMaxRowNnz{ 20 };
  HYPRE_Int iluLocalReordering{ 0 };

  HYPRE_Real maxRowSum{ -1.0 };
  HYPRE_Real strongThreshold{ -1.0 };
  HYPRE_Int numFunctions{ -1 };
  HYPRE_Int filterFunctions{ -1 };
  HYPRE_Int pMaxElmts{ -1 };

  HYPRE_Int aggressiveNumLevels{ -1 };
  HYPRE_Int aggressiveInterpType{ -1 };
  HYPRE_Int aggressivePMaxElmts{ -1 };
  HYPRE_Int coarseningType{ -1 };

  HYPRE_Int relaxType{ -1 };
  HYPRE_Int downRelaxType{ -1 };
  HYPRE_Int upRelaxType{ -1 };
  HYPRE_Int coarseRelaxType{ -1 };
  HYPRE_Int numSweeps{ -1 };
  HYPRE_Int relaxOrder{ -1 };
};

inline BoomerAMGParameters displacementAMGParameters( integer const separateComponents,
                                                       bool const filterFunctions,
                                                       bool const useALMSmoother = false )
{
  BoomerAMGParameters result;
  result.maxRowSum = 1.0;
  result.strongThreshold = useALMSmoother ? 0.8 : 0.6;
  result.numFunctions = 3;
  result.filterFunctions = filterFunctions ? separateComponents : 0;

  if( useALMSmoother )
  {
    result.pMaxElmts = 20;
    result.aggressiveNumLevels = 1;
#if GEOS_USE_HYPRE_DEVICE != GEOS_USE_HYPRE_CUDA && GEOS_USE_HYPRE_DEVICE != GEOS_USE_HYPRE_HIP
    result.coarseningType = hypre::getAMGCoarseningType( LinearSolverParameters::AMG::CoarseningType::Falgout );
#endif
  }

#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
  result.coarseningType = hypre::getAMGCoarseningType( LinearSolverParameters::AMG::CoarseningType::PMIS );
  result.relaxType = hypre::getAMGRelaxationType( LinearSolverParameters::AMG::SmootherType::chebyshev );
  result.numSweeps = useALMSmoother ? 2 : 1;
#else
  if( useALMSmoother )
  {
    result.downRelaxType = 89;
    result.upRelaxType = 89;
    result.coarseRelaxType = 9;
    result.numSweeps = 2;
    result.relaxOrder = 0;
  }
  else
  {
    result.relaxOrder = 1;
  }
#endif

  return result;
}

inline BoomerAMGParameters almBubbleAMGParameters()
{
  BoomerAMGParameters result = displacementAMGParameters( 0, false, true );
  result.strongThreshold = 0.75;
  result.filterFunctions = 0;
  result.pMaxElmts = 10;
  result.aggressiveNumLevels = -1;
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
  result.coarseningType = hypre::getAMGCoarseningType( LinearSolverParameters::AMG::CoarseningType::PMIS );
#else
  result.coarseningType = -1;
#endif
  result.numSweeps = 1;
  return result;
}

inline BoomerAMGParameters almReservoirDisplacementAMGParameters( integer const separateComponents )
{
  BoomerAMGParameters result = displacementAMGParameters( separateComponents, true );
#if GEOS_USE_HYPRE_DEVICE != GEOS_USE_HYPRE_CUDA && GEOS_USE_HYPRE_DEVICE != GEOS_USE_HYPRE_HIP
  result.downRelaxType = 89;
  result.upRelaxType = 89;
  result.coarseRelaxType = 9;
  result.numSweeps = 1;
  result.relaxOrder = 0;
#endif
  return result;
}

inline BoomerAMGParameters pressureAMGParameters( HYPRE_Int const minCoarseSize = -1 )
{
  BoomerAMGParameters result;
  result.minCoarseSize = minCoarseSize;
  result.aggressiveNumLevels = 1;
  result.aggressivePMaxElmts = 20;
  result.aggressiveInterpType = hypre::getAMGAggressiveInterpolationType( LinearSolverParameters::AMG::AggInterpType::multipass );

#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
  result.aggressiveInterpType = hypre::getAMGAggressiveInterpolationType( LinearSolverParameters::AMG::AggInterpType::modifiedExtendedE );
  result.coarseningType = hypre::getAMGCoarseningType( LinearSolverParameters::AMG::CoarseningType::PMIS );
  result.maxRowSum = 1.0;
  result.relaxType = hypre::getAMGRelaxationType( LinearSolverParameters::AMG::SmootherType::l1jacobi );
  result.numSweeps = 2;
#else
  result.relaxOrder = 1;
#endif

  return result;
}

inline BoomerAMGParameters pressureTemperatureAMGParameters()
{
  BoomerAMGParameters result;
  result.aggressiveNumLevels = 1;
  result.aggressivePMaxElmts = 16;
  result.numFunctions = 2;

#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
  result.aggressiveNumLevels = 0;
  result.coarseningType = hypre::getAMGCoarseningType( LinearSolverParameters::AMG::CoarseningType::PMIS );
  result.maxRowSum = 1.0;
  result.relaxType = hypre::getAMGRelaxationType( LinearSolverParameters::AMG::SmootherType::l1jacobi );
  result.numSweeps = 2;
#else
  result.relaxOrder = 1;
#endif

  return result;
}

inline void configureBoomerAMG( HyprePrecWrapper & solver,
                                BoomerAMGParameters const & params )
{
  GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &solver.ptr ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetILULocalReordering( solver.ptr, params.iluLocalReordering ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetILUMaxRowNnz( solver.ptr, params.smoothMaxRowNnz ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( solver.ptr, params.tolerance ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( solver.ptr, params.maxIterations ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( solver.ptr, params.printLevel ) );
  if( params.minCoarseSize >= 0 )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMinCoarseSize( solver.ptr, params.minCoarseSize ) );
  }
  GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxCoarseSize( solver.ptr, params.maxCoarseSize ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetSmoothType( solver.ptr, params.smoothType ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetSmoothNumLevels( solver.ptr, params.smoothNumLevels ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetSmoothNumSweeps( solver.ptr, params.smoothNumSweeps ) );

  if( params.maxRowSum >= 0.0 )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxRowSum( solver.ptr, params.maxRowSum ) );
  }
  if( params.strongThreshold >= 0.0 )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetStrongThreshold( solver.ptr, params.strongThreshold ) );
  }
  if( params.numFunctions >= 0 )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumFunctions( solver.ptr, params.numFunctions ) );
  }
  if( params.filterFunctions >= 0 )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetFilterFunctions( solver.ptr, params.filterFunctions ) );
  }
  if( params.pMaxElmts >= 0 )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPMaxElmts( solver.ptr, params.pMaxElmts ) );
  }
  if( params.aggressiveNumLevels >= 0 )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetAggNumLevels( solver.ptr, params.aggressiveNumLevels ) );
  }
  if( params.aggressiveInterpType >= 0 )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetAggInterpType( solver.ptr, params.aggressiveInterpType ) );
  }
  if( params.aggressivePMaxElmts >= 0 )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetAggPMaxElmts( solver.ptr, params.aggressivePMaxElmts ) );
  }
  if( params.coarseningType >= 0 )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCoarsenType( solver.ptr, params.coarseningType ) );
  }
  if( params.relaxType >= 0 )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( solver.ptr, params.relaxType ) );
  }
  if( params.numSweeps >= 0 )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( solver.ptr, params.numSweeps ) );
  }
  if( params.relaxOrder >= 0 )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxOrder( solver.ptr, params.relaxOrder ) );
  }
  if( params.downRelaxType >= 0 )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleRelaxType( solver.ptr, params.downRelaxType, 1 ) );
  }
  if( params.upRelaxType >= 0 )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleRelaxType( solver.ptr, params.upRelaxType, 2 ) );
  }
  if( params.coarseRelaxType >= 0 )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleRelaxType( solver.ptr, params.coarseRelaxType, 3 ) );
  }
}

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
                  arrayView1d< int const > const & numComponentsPerField,
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

  static constexpr HYPRE_Int numLevels = NLEVEL;              ///< Number of levels

protected:

  HYPRE_Int m_numBlocks{ 0 };                                 ///< Number of different matrix blocks treated separately

  stdVector< HYPRE_Int > m_labels[numLevels]{};             ///< Dof labels kept at each level
  HYPRE_Int m_numLabels[numLevels]{ -1 };                     ///< Number of dof labels kept
  HYPRE_Int * m_ptrLabels[numLevels]{ nullptr };              ///< Pointers to each level's labels, as consumed by MGR

  MGRFRelaxationType m_levelFRelaxType[numLevels];            ///< F-relaxation type for each level
  HYPRE_Int m_levelFRelaxIters[numLevels]{ -1 };              ///< Number of F-relaxation iterations for each level
  MGRInterpolationType m_levelInterpType[numLevels];          ///< Interpolation type for each level
  MGRRestrictionType m_levelRestrictType[numLevels];          ///< Restriction type for each level
  MGRCoarseGridMethod m_levelCoarseGridMethod[numLevels];     ///< Coarse grid method for each level
  MGRGlobalSmootherType m_levelGlobalSmootherType[numLevels]; ///< Global smoother type for each level
  HYPRE_Int m_levelGlobalSmootherIters[numLevels]{ -1 };      ///< Number of global smoother iterations for each level
  HYPRE_Real m_coarseGridThreshold{ 1.0e-20 };                ///< Coarse grid truncation threshold

  // TODO: the following options are currently commented out in MGR's code.
  //       Let's consider their use when re-enable in hypre
  // HYPRE_Int m_numRestrictSweeps{ -1 }; ///< Number of restrict sweeps
  // HYPRE_Int m_numInterpSweeps{ -1 };   ///< Number of interpolation sweeps

  /**
   * @brief Total number of dof labels, i.e. the sum of all fields' components.
   * @param numComponentsPerField number of components of each field
   * @return the total number of blocks
   */
  static HYPRE_Int totalNumBlocks( arrayView1d< int const > const & numComponentsPerField )
  {
    HYPRE_Int result = 0;
    for( localIndex i = 0; i < numComponentsPerField.size(); ++i )
    {
      result += LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[i] );
    }
    return result;
  }

  /**
   * @brief Constructor.
   * @param numBlocks number of blocks
   */
  explicit MGRStrategyBase( HYPRE_Int const numBlocks )
    : m_numBlocks( numBlocks )
  {
    for( HYPRE_Int i = 0; i < numLevels; ++i )
    {
      m_levelFRelaxType[i]         = MGRFRelaxationType::jacobi;
      m_levelInterpType[i]         = MGRInterpolationType::jacobi;
      m_levelRestrictType[i]       = MGRRestrictionType::injection;
      m_levelCoarseGridMethod[i]   = MGRCoarseGridMethod::galerkin;
      m_levelGlobalSmootherType[i] = MGRGlobalSmootherType::none;
    }
  }

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

  /**
   * @brief Helper function that sets the reduction features common to all mgr strategies
   * @param precond the preconditioner wrapper
   * @param mgrData auxiliary MGR data
   */
  void setReduction( HyprePrecWrapper & precond,
                     HypreMGRData & mgrData )

  {
    normalizeReductionParameters();
    MGRParameters const mgrParameters = defaultMGRParameters();

    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetCpointsByPointMarkerArray( precond.ptr,
                                                                 m_numBlocks, numLevels,
                                                                 m_numLabels, m_ptrLabels,
                                                                 mgrData.pointMarkers.data() ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelFRelaxType( precond.ptr, toUnderlyingPtr( m_levelFRelaxType ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelNumRelaxSweeps( precond.ptr, m_levelFRelaxIters ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( precond.ptr, toUnderlyingPtr( m_levelInterpType ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelRestrictType( precond.ptr, toUnderlyingPtr( m_levelRestrictType ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseGridMethod( precond.ptr, toUnderlyingPtr( m_levelCoarseGridMethod ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelSmoothType( precond.ptr, toUnderlyingPtr( m_levelGlobalSmootherType ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelSmoothIters( precond.ptr, m_levelGlobalSmootherIters ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetTruncateCoarseGridThreshold( precond.ptr, m_coarseGridThreshold ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( precond.ptr, mgrParameters.nonCpointsToFpoints ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetNonGalerkinMaxElmts( precond.ptr, mgrParameters.nonGalerkinMaxElmts ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetPMaxElmts( precond.ptr, mgrParameters.pMaxElmts ) );
  }

public:
  /**
   * @brief Normalize MGR iteration counts without touching a HYPRE handle.
   *
   * Generated hypredrive YAML needs the same normalization as the legacy
   * setup, but does not need to construct a temporary MGR object.
   */
  void normalizeReductionParameters()
  {
    // Ensure that if no F-relaxation or global smoothing is chosen the corresponding number
    // of iteration is set to 0
    for( HYPRE_Int i = 0; i < numLevels; ++i )
    {
      if( m_levelFRelaxType[i] == MGRFRelaxationType::none )
      {
        m_levelFRelaxIters[i] = 0;
      }
      if( m_levelGlobalSmootherType[i] == MGRGlobalSmootherType::none )
      {
        m_levelGlobalSmootherIters[i] = 0;
      }
    }
  }

protected:

  /**
   * @brief Set up BoomerAMG to perform the solve for the displacement system
   * @param solver the solver wrapper
   * @param separateComponents flag controlling the use of the separate displacement component (SDC) approximation
   */
  void setDisplacementAMG( HyprePrecWrapper & solver,
                           integer const & separateComponents )
  {
    configureBoomerAMG( solver, displacementAMGParameters( separateComponents, true ) );

    solver.setup = HYPRE_BoomerAMGSetup;
    solver.solve = HYPRE_BoomerAMGSolve;
    solver.destroy = HYPRE_BoomerAMGDestroy;
  }

  /**
   * @brief Set up one of the two BoomerAMG instances used by the fully
   *        coupled single-phase ALM hierarchy.
   * @param solver solver wrapper to initialize
   * @param separateComponents whether displacement components are filtered
   * @param bubbleCoarse true for the inner (bubble-displacement) coarse solve
   *
   * The CPU values mirror the reference nested-MGR YAML.  The device branch
   * retains the device-safe smoother choices used by the existing strategies.
   */
  void setALMDisplacementAMG( HyprePrecWrapper & solver,
                              integer const separateComponents,
                              bool const bubbleCoarse )
  {
    configureBoomerAMG( solver, bubbleCoarse
                        ? almBubbleAMGParameters()
                        : displacementAMGParameters( separateComponents, true, true ) );

    solver.setup = HYPRE_BoomerAMGSetup;
    solver.solve = HYPRE_BoomerAMGSolve;
    solver.destroy = HYPRE_BoomerAMGDestroy;
  }

  /**
   * @brief Set up the displacement F-solver used by the reservoir ALM strategy.
   * @param solver solver wrapper to initialize
   * @param separateComponents whether displacement components are filtered
   */
  void setALMReservoirDisplacementAMG( HyprePrecWrapper & solver,
                                       integer const separateComponents )
  {
    configureBoomerAMG( solver, almReservoirDisplacementAMGParameters( separateComponents ) );

    solver.setup = HYPRE_BoomerAMGSetup;
    solver.solve = HYPRE_BoomerAMGSolve;
    solver.destroy = HYPRE_BoomerAMGDestroy;
  }

  /**
   * @brief Set up BoomerAMG to perform the solve for the pressure system
   * @param solver the solver wrapper
   */
  void setPressureAMG( HyprePrecWrapper & solver,
                       HYPRE_Int const minCoarseSize = -1 )
  {
    configureBoomerAMG( solver, pressureAMGParameters( minCoarseSize ) );

    solver.setup = HYPRE_BoomerAMGSetup;
    solver.solve = HYPRE_BoomerAMGSolve;
    solver.destroy = HYPRE_BoomerAMGDestroy;
  }

  /**
   * @brief Set up BoomerAMG to perform the solve for the pressure/temperature system
   * @param solver the solver wrapper
   */
  void setPressureTemperatureAMG( HyprePrecWrapper & solver )
  {
    configureBoomerAMG( solver, pressureTemperatureAMGParameters() );

    solver.setup = HYPRE_BoomerAMGSetup;
    solver.solve = HYPRE_BoomerAMGSolve;
    solver.destroy = HYPRE_BoomerAMGDestroy;
  }

  /**
   * @brief Set up BoomerAMG to perform the mechanics F-solve for the first F-relaxation
   * @param precond the preconditioner wrapper
   * @param mgrData auxiliary MGR data
   * @param separateComponents flag controlling the use of the separate displacement component (SDC) approximation
   *
   * @note This function should be rethought once MGR allows for customizing boomerAMG (or
   *       any other solver) for F-relaxation at any level
   */
  void setMechanicsFSolver( HyprePrecWrapper & precond,
                            HypreMGRData & mgrData,
                            integer const & separateComponents )
  {
    setDisplacementAMG( mgrData.mechSolver, separateComponents );
    HYPRE_MGRSetFSolver( precond.ptr, mgrData.mechSolver.solve, mgrData.mechSolver.setup, mgrData.mechSolver.ptr );
  }

  /**
   * @brief Configure the displacement F-solver attached to a specific MGR level.
   */
  void setMechanicsFSolverAtLevel( HyprePrecWrapper & precond,
                                    HypreMGRData & mgrData,
                                    integer const & separateComponents,
                                    HYPRE_Int const level )
  {
    setDisplacementAMG( mgrData.mechSolver, separateComponents );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetFSolverAtLevel( precond.ptr, mgrData.mechSolver.ptr, level ) );
  }
  /**
   * @brief Set up an explicitly configured ILU(0) F-solver for a given MGR level
   * @param level the MGR level whose F system is solved with ILU
   * @param precond the preconditioner wrapper
   * @param mgrData auxiliary MGR data
   *
   * @note hypre's internal ILU F-relaxation (MGRFRelaxationType::ilu without an attached
   *       F-solver) is configured with different defaults and has been observed to stall
   *       on hybrid FVM cell-block eliminations; an explicitly configured ILU F-solver
   *       matches the configuration hypredrive uses and is robust.
   */
  void setILUFSolverAtLevel( HYPRE_Int const level,
                             HyprePrecWrapper & precond,
                             HypreMGRData & mgrData )
  {
    HYPRE_ILUCreate( &mgrData.mechSolver.ptr );
    HYPRE_ILUSetType( mgrData.mechSolver.ptr, 0 );
    HYPRE_ILUSetLevelOfFill( mgrData.mechSolver.ptr, 0 );
    HYPRE_ILUSetMaxIter( mgrData.mechSolver.ptr, 1 );
    HYPRE_ILUSetTol( mgrData.mechSolver.ptr, 0.0 );
    HYPRE_ILUSetLocalReordering( mgrData.mechSolver.ptr, 0 );
    HYPRE_ILUSetPrintLevel( mgrData.mechSolver.ptr, 0 );

    mgrData.mechSolver.setup = HYPRE_ILUSetup;
    mgrData.mechSolver.solve = HYPRE_ILUSolve;
    mgrData.mechSolver.destroy = HYPRE_ILUDestroy;

    HYPRE_MGRSetFSolverAtLevel( precond.ptr, mgrData.mechSolver.ptr, level );
  }

  /**
   * @brief
   *
   * @param solver
   */
  void setILUCoarseSolver( HyprePrecWrapper & solver )
  {
    /* (Required) Create ILU solver */
    HYPRE_ILUCreate( &solver.ptr );

    /* (Recommended) General solver options */
    int const ilu_type = 0; /* 0, 1, 10, 11, 20, 21, 30, 31, 40, 41, 50 */
    int const max_iter = 1;
    double const tol = 0.0;
    int const reordering = 0; /* 0: none, 1: RCM */
    int const print_level = 0;
    HYPRE_ILUSetType( solver.ptr, ilu_type );
    HYPRE_ILUSetMaxIter( solver.ptr, max_iter );
    HYPRE_ILUSetTol( solver.ptr, tol );
    HYPRE_ILUSetLocalReordering( solver.ptr, reordering );
    HYPRE_ILUSetPrintLevel( solver.ptr, print_level );

    solver.setup = HYPRE_ILUSetup;
    solver.solve = HYPRE_ILUSolve;
    solver.destroy = HYPRE_ILUDestroy;
  }

};

/**
 * @brief Create the MGR preconditioner object.
 * @param params preconditioner parameters
 * @param dofManager pointer to DofManager for the linear system
 * @param precond the preconditioner
 * @param mgrData auxiliary data for MGR
 */
void createMGR( LinearSolverParameters const & params,
                DofManager const * const dofManager,
                HyprePrecWrapper & precond,
                HypreMGRData & mgrData );

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSTRATEGIES_HPP_*/
