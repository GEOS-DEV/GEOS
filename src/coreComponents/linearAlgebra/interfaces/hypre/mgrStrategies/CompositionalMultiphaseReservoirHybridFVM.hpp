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

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASERESERVOIRHYBRIDFVM_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASERESERVOIRHYBRIDFVM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief CompositionalMultiphaseReservoirHybridFVM strategy.
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
 *   - 1st level: eliminate the well block
 *   - 2nd level: eliminate the reservoir density associated with the volume constraint
 *   - 3rd level: eliminate the remaining the reservoir densities
 *   - 4th level: eliminate the cell-centered pressure
 *   - The coarse grid is the face pressure system and is solved with BoomerAMG
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

    // Level 0: eliminate the well block
    m_labels[0].resize( numResLabels );
    std::iota( m_labels[0].begin(), m_labels[0].end(), 0 );
    // Level 1: eliminate the last density of the reservoir block
    m_labels[1].resize( numResLabels - 1 );
    std::iota( m_labels[1].begin(), m_labels[1].begin() + numResCellCenteredLabels - 1, 0 );
    m_labels[1][numResCellCenteredLabels-1] = numResCellCenteredLabels;
    // Level 2: eliminate remaining densities of the reservoir block
    m_labels[2].resize( 2 );
    m_labels[2][0] = 0;
    m_labels[2][1] = numResCellCenteredLabels;
    // Level 3: eliminate reservoir cell centered pressure
    m_labels[3].resize( 1 );
    m_labels[3][0] = numResCellCenteredLabels;

    setupLabels();

    // Level 0
#if GEOS_USE_HYPRE_DEVICE != GEOS_USE_HYPRE_CPU
    m_levelFRelaxType[0]          = MGRFRelaxationType::gsElimWInverse;
#else
    GEOS_ERROR( "General support for Gaussian elimination on GPU not available yet" );
#endif
    m_levelFRelaxIters[0]         = 1;
    m_levelInterpType[0]          = MGRInterpolationType::blockJacobi;
    m_levelRestrictType[0]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[0]  = MGRGlobalSmootherType::none;

    // Level 1
    m_levelFRelaxType[1]          = MGRFRelaxationType::jacobi;
    m_levelFRelaxIters[1]         = 1;
    m_levelInterpType[1]          = MGRInterpolationType::jacobi;
    m_levelRestrictType[1]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[1]  = MGRGlobalSmootherType::none;

    // Level 2
    m_levelFRelaxType[2]          = MGRFRelaxationType::jacobi;
    m_levelFRelaxIters[2]         = 1;
    m_levelInterpType[2]          = MGRInterpolationType::injection;
    m_levelRestrictType[2]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[2]    = MGRCoarseGridMethod::cprLikeBlockDiag;
    m_levelGlobalSmootherType[2]  = MGRGlobalSmootherType::none;

    // Level 3
    m_levelFRelaxType[3]          = MGRFRelaxationType::jacobi;
    m_levelFRelaxIters[3]         = 1;
    m_levelInterpType[3]          = MGRInterpolationType::jacobi;
    m_levelRestrictType[3]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[3]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[3]  = MGRGlobalSmootherType::blockGaussSeidel;
    m_levelGlobalSmootherIters[3] = 1;
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
    // if the wells are shut, using Gaussian elimination as F-relaxation for the well block is an overkill
    // in that case, we just use Jacobi
    if( mgrParams.areWellsShut )
    {
      m_levelFRelaxType[0] = MGRFRelaxationType::jacobi;
    }

    setReduction( precond, mgrData );

    // Configure the BoomerAMG solver used as mgr coarse solver for the pressure reduced system
    setPressureAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASERESERVOIRHYBRIDFVM_HPP_*/
