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
 * 2-level MGR reduction strategy
 *   - 1st level: eliminate the well block
 *   - 2nd level: eliminate all cell-centered reservoir dofs (pressure and densities);
 *     the cell-centered block is block-diagonal (cell dofs couple only within a cell
 *     and to the faces), so block-Jacobi interpolation with a Galerkin coarse grid
 *     yields the exact face-pressure Schur complement, and an ILU F-relaxation
 *     solves the block-diagonal F system at linear cost.
 *   - The coarse grid is the face pressure system and is solved with BoomerAMG
 *
 * Note: a multi-level variant that eliminates the densities first cannot use the
 * True-IMPES (blockColLumped) restriction: it requires every C point to pair
 * positionally with one F block of uniform size, while the intermediate C spaces
 * contain both cell and face pressures (hypre's block column-lumped helper fails
 * on the resulting row/column block-count mismatch).
 */
class CompositionalMultiphaseReservoirHybridFVM : public MGRStrategyBase< 2 >
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
    // Level 1: eliminate all cell-centered reservoir dofs, keep the face-centered pressure
    m_labels[1].resize( 1 );
    m_labels[1][0] = numResCellCenteredLabels;

    setupLabels();

    // Level 0
    m_levelFRelaxType[0]          = MGRFRelaxationType::gsElimWInverse;
    m_levelFRelaxIters[0]         = 1;
    m_levelInterpType[0]          = MGRInterpolationType::blockJacobi;
    m_levelRestrictType[0]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[0]  = MGRGlobalSmootherType::none;

    // Level 1
    m_levelFRelaxType[1]          = MGRFRelaxationType::ilu;
    m_levelFRelaxIters[1]         = 1;
    m_levelInterpType[1]          = MGRInterpolationType::blockJacobi;
    m_levelRestrictType[1]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[1]  = MGRGlobalSmootherType::none;
  }

  /**
   * @brief Setup the MGR strategy.
   * @param mgrParams parameters for the configuration of the MGR recipe
   * @param precond preconditioner wrapper
   * @param mgrData auxiliary MGR data
   */
  void setup( LinearSolverParameters::MGR const & GEOS_UNUSED_PARAM( mgrParams ),
              HyprePrecWrapper & precond,
              HypreMGRData & mgrData )
  {
    // Note: unlike other reservoir strategies, the well-block F-relaxation is not downgraded
    // to Jacobi when the wells are shut: a single Jacobi sweep on the shut well block leads
    // to Krylov stagnation for this hybrid FVM reduction, while exact elimination of the
    // small well block is cheap.
    setReduction( precond, mgrData );

    // Attach an explicitly configured ILU F-solver for the cell-centered block elimination
    setILUFSolverAtLevel( 1, precond, mgrData );

    // Configure the BoomerAMG solver used as mgr coarse solver for the pressure reduced system
    setPressureAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASERESERVOIRHYBRIDFVM_HPP_*/
