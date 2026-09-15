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
 * @file SinglePhasePoromechanicsConformingFracturesALM.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASEPOROMECHANICSCONFORMINGFRACTURESALM_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASEPOROMECHANICSCONFORMINGFRACTURESALM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

#include <cstddef>

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief Adapter that lets an MGR instance be used as an MGR F-relaxation
 * solver.
 *
 * HYPRE's public MGR API accepts callbacks for the first F-relaxation, while
 * subsequent solves dispatch through the base hypre_Solver layout. Keep the
 * callback members first so that the adapter is valid in both paths. The
 * nested MGR and its two user-owned AMG solvers are released by destroy().
 */
struct SinglePhaseALMNestedMGR
{
  HyprePrecWrapper::SetupFunc setup{};
  HyprePrecWrapper::SolveFunc solve{};
  HyprePrecWrapper::DestroyFunc destroy{};
  HYPRE_Int isSetup{ 0 };
  HYPRE_Solver mgr{};
  HyprePrecWrapper displacementSolver{};
  HyprePrecWrapper bubbleSolver{};
  stdVector< HYPRE_Int > pointMarkers;
};

inline SinglePhaseALMNestedMGR * singlePhaseALMNestedMGR( HYPRE_Solver solver )
{
  return reinterpret_cast< SinglePhaseALMNestedMGR * >( solver );
}

inline HYPRE_Int singlePhaseALMNestedMGRSetup( HYPRE_Solver solver,
                                               HYPRE_ParCSRMatrix A,
                                               HYPRE_ParVector b,
                                               HYPRE_ParVector x )
{
  SinglePhaseALMNestedMGR * const nested = singlePhaseALMNestedMGR( solver );
  if( nested == nullptr || nested->mgr == nullptr )
  {
    return 1;
  }

  HYPRE_Int const ierr = HYPRE_MGRSetup( nested->mgr, A, b, x );
  if( ierr == 0 )
  {
    nested->isSetup = 1;
  }
  return ierr;
}

inline HYPRE_Int singlePhaseALMNestedMGRSolve( HYPRE_Solver solver,
                                               HYPRE_ParCSRMatrix A,
                                               HYPRE_ParVector b,
                                               HYPRE_ParVector x )
{
  SinglePhaseALMNestedMGR * const nested = singlePhaseALMNestedMGR( solver );
  if( nested == nullptr || nested->mgr == nullptr )
  {
    return 1;
  }
  return HYPRE_MGRSolve( nested->mgr, A, b, x );
}

inline HYPRE_Int singlePhaseALMNestedMGRDestroy( HYPRE_Solver solver )
{
  SinglePhaseALMNestedMGR * const nested = singlePhaseALMNestedMGR( solver );
  if( nested == nullptr )
  {
    return 0;
  }

  HYPRE_Int ierr = 0;
  if( nested->mgr != nullptr )
  {
    ierr |= HYPRE_MGRDestroy( nested->mgr );
    nested->mgr = nullptr;
  }
  if( nested->displacementSolver.ptr != nullptr && nested->displacementSolver.destroy != nullptr )
  {
    ierr |= nested->displacementSolver.destroy( nested->displacementSolver.ptr );
    nested->displacementSolver.ptr = nullptr;
  }
  if( nested->bubbleSolver.ptr != nullptr && nested->bubbleSolver.destroy != nullptr )
  {
    ierr |= nested->bubbleSolver.destroy( nested->bubbleSolver.ptr );
    nested->bubbleSolver.ptr = nullptr;
  }
  delete nested;
  return ierr;
}

static_assert( offsetof( SinglePhaseALMNestedMGR, isSetup ) ==
               3 * sizeof( HYPRE_PtrToSolverFcn ),
               "The nested MGR adapter must start with the hypre_Solver callback layout" );

/**
 * @brief SinglePhasePoromechanicsConformingFracturesALM strategy.
 *
 * dofLabel: 0 = displacement, x-component
 * dofLabel: 1 = displacement, y-component
 * dofLabel: 2 = displacement, z-component
 * dofLabel: 3 = displacement bubble function, x-component
 * dofLabel: 4 = displacement bubble function, y-component
 * dofLabel: 5 = displacement bubble function, z-component
 * dofLabel: 6 = pressure (cell elem + fracture elems)
 *
 * Well unknowns are not handled here: see
 * SinglePhasePoromechanicsConformingFracturesALMReservoirFVM, which adds a
 * third level to eliminate the well block before the coarse solve.
 *
 * Ingredients:
 * 1. The outer MGR F-block is the fully coupled displacement/bubble block.
 * 2. Its F-relaxation is a nested MGR that eliminates nodal displacement
 *    (0,1,2) and leaves bubble displacement (3,4,5) for the inner coarse AMG.
 * 3. Both displacement AMG solves use three functions and the ALM smoother
 *    settings used by the reference strategy.
 * 4. The pressure Schur complement is solved with BoomerAMG.
 * 5. Both MGR V-cycles use pre-relaxation only and no global smoother.
 */
class SinglePhasePoromechanicsConformingFracturesALM : public MGRStrategyBase< 1 >
{
public:

  /**
   * @brief Constructor.
   */
  explicit SinglePhasePoromechanicsConformingFracturesALM( arrayView1d< int const > const & numComponentsPerField )
    : MGRStrategyBase( totalNumBlocks( numComponentsPerField ) )
  {
    GEOS_ERROR_IF_NE_MSG( numComponentsPerField.size(), 3,
                          "singlePhasePoromechanicsConformingFracturesALM requires exactly the displacement, "
                          "bubble-displacement, and pressure fields. Use "
                          "singlePhasePoromechanicsConformingFracturesALMReservoirFVM when wells are present." );
    GEOS_ERROR_IF_NE_MSG( numComponentsPerField[0], 3,
                          "singlePhasePoromechanicsConformingFracturesALM requires three displacement components" );
    GEOS_ERROR_IF_NE_MSG( numComponentsPerField[1], 3,
                          "singlePhasePoromechanicsConformingFracturesALM requires three bubble-displacement components" );

    // The outer MGR eliminates the fully coupled displacement/bubble block.
    // Its F-relaxation is a nested MGR that eliminates nodal displacement and
    // leaves the bubble displacement block for the inner coarse AMG solve.
    HYPRE_Int const numDisplacementLabels = LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] );
    HYPRE_Int const pressureLabel = numDisplacementLabels +
                                    LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[1] );
    for( HYPRE_Int label = pressureLabel; label < m_numBlocks; ++label )
    {
      m_labels[0].push_back( label );
    }

    setupLabels();

    // Level 0: nested MGR F-relaxation over all six displacement unknowns.
    // HYPRE requires the AMG F-relaxation enum for the callback-based solver.
    m_levelFRelaxType[0]          = MGRFRelaxationType::amgVCycle;
    m_levelFRelaxIters[0]         = 1;
    m_levelGlobalSmootherType[0]  = MGRGlobalSmootherType::none;
    m_levelGlobalSmootherIters[0] = 0;
    m_levelInterpType[0]          = MGRInterpolationType::jacobi;
    m_levelRestrictType[0]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]    = MGRCoarseGridMethod::galerkin;

  }

  /**
   * @brief Setup the MGR strategy.
   * @param precond preconditioner wrapper
   * @param mgrData auxiliary MGR data
   */
  void setup( LinearSolverParameters::MGR const & mgrParams,
              HyprePrecWrapper & precond,
              HypreMGRData & mgrData )
  {
    GEOS_UNUSED_VAR( mgrParams );
    setReduction( precond, mgrData );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetCycleType( precond.ptr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetFRelaxCycle( precond.ptr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetGlobalSmoothCycle( precond.ptr, 1 ) );

    auto * const nested = new SinglePhaseALMNestedMGR;
    nested->setup = singlePhaseALMNestedMGRSetup;
    nested->solve = singlePhaseALMNestedMGRSolve;
    nested->destroy = singlePhaseALMNestedMGRDestroy;
    mgrData.nestedSolver.ptr = reinterpret_cast< HYPRE_Solver >( nested );
    mgrData.nestedSolver.setup = nested->setup;
    mgrData.nestedSolver.solve = nested->solve;
    mgrData.nestedSolver.destroy = nested->destroy;

    // The outer F block contains labels 0..5. The marker order must match
    // HYPRE's projected A_FF row order, hence the filtering of the original
    // local marker array rather than constructing a synthetic marker list.
    nested->pointMarkers.reserve( mgrData.pointMarkers.size() );
    for( HYPRE_Int const marker : mgrData.pointMarkers )
    {
      if( marker < 6 )
      {
        nested->pointMarkers.push_back( marker );
      }
    }

    GEOS_LAI_CHECK_ERROR( HYPRE_MGRCreate( &nested->mgr ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetTol( nested->mgr, 0.0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetMaxIter( nested->mgr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetPrintLevel( nested->mgr, 0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetCycleType( nested->mgr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetFRelaxCycle( nested->mgr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetGlobalSmoothCycle( nested->mgr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( nested->mgr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetNonGalerkinMaxElmts( nested->mgr, 1 ) );

    HYPRE_Int innerCoarseLabels[3] = { 3, 4, 5 };
    HYPRE_Int * innerCoarseLabelsPtr[1] = { innerCoarseLabels };
    HYPRE_Int numInnerCoarseLabels[1] = { 3 };
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetCpointsByPointMarkerArray( nested->mgr,
                                                                 6, 1,
                                                                 numInnerCoarseLabels,
                                                                 innerCoarseLabelsPtr,
                                                                 nested->pointMarkers.data() ) );

    HYPRE_Int fRelaxType[1] = { static_cast< HYPRE_Int >( MGRFRelaxationType::amgVCycle ) };
    HYPRE_Int fRelaxIters[1] = { 1 };
    HYPRE_Int interpolationType[1] = { static_cast< HYPRE_Int >( MGRInterpolationType::injection ) };
    HYPRE_Int restrictionType[1] = { static_cast< HYPRE_Int >( MGRRestrictionType::injection ) };
    HYPRE_Int coarseGridMethod[1] = { static_cast< HYPRE_Int >( MGRCoarseGridMethod::galerkin ) };
    HYPRE_Int globalSmoothType[1] = { static_cast< HYPRE_Int >( MGRGlobalSmootherType::none ) };
    HYPRE_Int globalSmoothIters[1] = { 0 };
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelFRelaxType( nested->mgr, fRelaxType ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelNumRelaxSweeps( nested->mgr, fRelaxIters ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( nested->mgr, interpolationType ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelRestrictType( nested->mgr, restrictionType ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseGridMethod( nested->mgr, coarseGridMethod ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelSmoothType( nested->mgr, globalSmoothType ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetLevelSmoothIters( nested->mgr, globalSmoothIters ) );

    // Inner level 0: AMG for nodal displacement; inner coarsest level:
    // AMG for the bubble-displacement block.
    setALMDisplacementAMG( nested->displacementSolver, 1, false );
    setALMDisplacementAMG( nested->bubbleSolver, 0, true );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetFSolverAtLevel( nested->mgr, nested->displacementSolver.ptr, 0 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseSolver( nested->mgr,
                                                    nested->bubbleSolver.solve,
                                                    nested->bubbleSolver.setup,
                                                    nested->bubbleSolver.ptr ) );

    // HYPRE invokes the callback during setup and subsequently dispatches via
    // the adapter's hypre_Solver-compatible first four fields.
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetFSolver( precond.ptr,
                                               nested->solve,
                                               nested->setup,
                                               mgrData.nestedSolver.ptr ) );

    // Configure the BoomerAMG solver used as the outer coarse solver for the
    // pressure reduced system.
    setPressureAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASEPOROMECHANICSCONFORMINGFRACTURESALM_HPP_*/
