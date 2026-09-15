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
 * @file SinglePhasePoromechanicsConformingFracturesALMReservoirFVM.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASEPOROMECHANICSCONFORMINGFRACTURESALMRESERVOIRFVM_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASEPOROMECHANICSCONFORMINGFRACTURESALMRESERVOIRFVM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief SinglePhasePoromechanicsConformingFracturesALMReservoirFVM strategy.
 *
 * dofLabel: 0 = displacement, x-component
 * dofLabel: 1 = displacement, y-component
 * dofLabel: 2 = displacement, z-component
 * dofLabel: 3 = displacement bubble function, x-component
 * dofLabel: 4 = displacement bubble function, y-component
 * dofLabel: 5 = displacement bubble function, z-component
 * dofLabel: 6 = pressure (cell elem + fracture elems)
 * dofLabel: 7 = well pressure
 * dofLabel: 8 = well rate
 *
 * Ingredients:
 * 1. Level 0 eliminates bubble displacement (3,4,5) with L1-Jacobi.
 * 2. Level 1 eliminates nodal displacement (0,1,2) with one BoomerAMG V-cycle.
 * 3. Level 2 eliminates the well block (7,8), as in
 *    SinglePhasePoromechanicsReservoirFVM. Leaving the well unknowns in the
 *    coarse grid would hand BoomerAMG the well rate and BHP constraint rows,
 *    which have no elliptic structure.
 * 4. The displacement AMG uses three functions, separate-component filtering,
 *    and symmetric L1 hybrid Gauss-Seidel relaxation on CPU.
 * 5. The reservoir pressure Schur complement is solved with BoomerAMG.
 * 6. The MGR V-cycle uses pre-relaxation only and no global smoother.
 */
class SinglePhasePoromechanicsConformingFracturesALMReservoirFVM : public MGRStrategyBase< 3 >
{
public:

  /**
   * @brief Constructor.
   * @param numComponentsPerField array with number of components for each field
   */
  explicit SinglePhasePoromechanicsConformingFracturesALMReservoirFVM( arrayView1d< int const > const & numComponentsPerField )
    : MGRStrategyBase( totalNumBlocks( numComponentsPerField ) )
  {
    GEOS_ERROR_IF_NE_MSG( numComponentsPerField.size(), 4,
                          "singlePhasePoromechanicsConformingFracturesALMReservoirFVM requires exactly the "
                          "displacement, bubble-displacement, pressure, and well fields. Any further field would "
                          "be swept into the well block eliminated on level 2." );

    // Eliminate bubble displacement first, then nodal displacement, then the
    // well block. Only the reservoir pressure survives on the coarsest grid.
    HYPRE_Int const numDisplacementLabels = LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] );
    HYPRE_Int const pressureLabel = numDisplacementLabels +
                                    LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[1] );
    HYPRE_Int const wellLabel = pressureLabel +
                                LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[2] );
    for( HYPRE_Int label = 0; label < numDisplacementLabels; ++label )
    {
      m_labels[0].push_back( label );
    }
    for( HYPRE_Int label = pressureLabel; label < m_numBlocks; ++label )
    {
      m_labels[0].push_back( label );
      m_labels[1].push_back( label );
    }
    for( HYPRE_Int label = pressureLabel; label < wellLabel; ++label )
    {
      m_labels[2].push_back( label );
    }

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

    // Level 2
    m_levelFRelaxType[2]         = MGRFRelaxationType::gsElimWInverse;
    m_levelFRelaxIters[2]        = 1;
    m_levelGlobalSmootherType[2] = MGRGlobalSmootherType::none;
    m_levelInterpType[2]         = MGRInterpolationType::blockJacobi;
    m_levelRestrictType[2]       = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[2]   = MGRCoarseGridMethod::galerkin;
  }

  /**
   * @brief Setup the MGR strategy.
   * @param mgrParams MGR configuration parameters
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

    // Configure the BoomerAMG solver used as mgr coarse solver for the reservoir pressure reduced system
    setPressureAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRSINGLEPHASEPOROMECHANICSCONFORMINGFRACTURESALMRESERVOIRFVM_HPP_*/
