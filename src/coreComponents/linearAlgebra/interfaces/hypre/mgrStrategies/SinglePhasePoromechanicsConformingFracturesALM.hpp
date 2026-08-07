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

namespace geos
{

namespace hypre
{

namespace mgr
{

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
 * Ingredients:
 * 1. Level 0 eliminates bubble displacement (3,4,5) with L1-Jacobi.
 * 2. Level 1 eliminates nodal displacement (0,1,2) with one BoomerAMG V-cycle.
 * 3. The displacement AMG uses three functions, separate-component filtering,
 *    and symmetric L1 hybrid Gauss-Seidel relaxation on CPU.
 * 4. The pressure Schur complement is solved with BoomerAMG.
 * 5. The MGR V-cycle uses pre-relaxation only and no global smoother.
 */
class SinglePhasePoromechanicsConformingFracturesALM : public MGRStrategyBase< 2 >
{
public:

  /**
   * @brief Constructor.
   */
  explicit SinglePhasePoromechanicsConformingFracturesALM( arrayView1d< int const > const & )
    : MGRStrategyBase( 7 )
  {

    // Eliminate bubble displacement first, then nodal displacement.
    m_labels[0] = { 0, 1, 2, 6 };
    m_labels[1] = { 6 };

    setupLabels();

    // Level 0
    m_levelFRelaxType[0]          = MGRFRelaxationType::l1jacobi;
    m_levelFRelaxIters[0]         = 1;
    m_levelGlobalSmootherType[0]  = MGRGlobalSmootherType::none;
    m_levelGlobalSmootherIters[0] = 0;
    m_levelInterpType[0]          = MGRInterpolationType::blockJacobi;
    m_levelRestrictType[0]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]    = MGRCoarseGridMethod::galerkin;

    // Level 1
    m_levelFRelaxType[1]         = MGRFRelaxationType::amgVCycle;
    m_levelFRelaxIters[1]        = 1;
    m_levelGlobalSmootherType[1] = MGRGlobalSmootherType::none;
    m_levelInterpType[1]         = MGRInterpolationType::jacobi;
    m_levelRestrictType[1]       = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1]   = MGRCoarseGridMethod::nonGalerkin;

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
    setReduction( precond, mgrData );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetCycleType( precond.ptr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetFRelaxCycle( precond.ptr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetGlobalSmoothCycle( precond.ptr, 1 ) );

    // Configure the BoomerAMG solver used as level-1 F-relaxation.
    setDisplacementAMG( mgrData.mechSolver, mgrParams.separateComponents );

#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCoarsenType( mgrData.mechSolver.ptr, hypre::getAMGCoarseningType( LinearSolverParameters::AMG::CoarseningType::PMIS ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( mgrData.mechSolver.ptr, hypre::getAMGRelaxationType( LinearSolverParameters::AMG::SmootherType::chebyshev ) ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( mgrData.mechSolver.ptr, 1 ) );
#else
    HYPRE_Int constexpr l1SymmetricHybridGaussSeidel = 89;
    HYPRE_Int constexpr gaussianElimination = 9;
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleRelaxType( mgrData.mechSolver.ptr, l1SymmetricHybridGaussSeidel, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleRelaxType( mgrData.mechSolver.ptr, l1SymmetricHybridGaussSeidel, 2 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleRelaxType( mgrData.mechSolver.ptr, gaussianElimination, 3 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( mgrData.mechSolver.ptr, 1 ) );
    GEOS_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxOrder( mgrData.mechSolver.ptr, 0 ) );
#endif
    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetFSolverAtLevel( precond.ptr, mgrData.mechSolver.ptr, 1 ) );

    // Configure the BoomerAMG solver used as mgr coarse solver for the pressure reduced system
    setPressureAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASEPOROMECHANICSCONFORMINGFRACTURESALM_HPP_*/
