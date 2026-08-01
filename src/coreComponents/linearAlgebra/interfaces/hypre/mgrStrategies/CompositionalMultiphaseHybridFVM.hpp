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
 * @file CompositionalMultiphaseHybridFVM.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASEHYBRIDFVM_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASEHYBRIDFVM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

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
 * Single-level MGR reduction strategy, mirroring the static condensation of the
 * hybrid FVM discretization:
 *   - 1st level: eliminate all cell-centered dofs (pressure and densities); the
 *     cell-centered block is block-diagonal (cell dofs couple only within a cell
 *     and to the faces), so block-Jacobi interpolation with a Galerkin coarse grid
 *     yields the exact face-pressure Schur complement, and an ILU F-relaxation
 *     solves the block-diagonal F system at linear cost.
 *   - The coarse grid is the interface pressure system and is solved with BoomerAMG
 *
 * Note: a multi-level variant that eliminates the densities first cannot use the
 * True-IMPES (blockColLumped) restriction: it requires every C point to pair
 * positionally with one F block of uniform size, while the intermediate C spaces
 * contain both cell and face pressures (hypre's block column-lumped helper fails
 * on the resulting row/column block-count mismatch).
 */
class CompositionalMultiphaseHybridFVM : public MGRStrategyBase< 1 >
{
public:
  /**
   * @brief Constructor.
   * @param numComponentsPerField array with number of components for each field
   */
  explicit CompositionalMultiphaseHybridFVM( arrayView1d< int const > const & numComponentsPerField )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] + numComponentsPerField[1] ) )
  {
    // Level 0: eliminate all cell-centered dofs, keep the face-centered pressure
    m_labels[0].resize( 1 );
    m_labels[0][0] = m_numBlocks - 1;

    setupLabels();

    // Level 0
    m_levelFRelaxType[0]          = MGRFRelaxationType::ilu;
    m_levelFRelaxIters[0]         = 1;
    m_levelInterpType[0]          = MGRInterpolationType::blockJacobi;
    m_levelRestrictType[0]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[0]  = MGRGlobalSmootherType::none;
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
    setReduction( precond, mgrData );

    // Attach an explicitly configured ILU F-solver for the cell-centered block elimination
    setILUFSolverAtLevel( 0, precond, mgrData );

    // Configure the BoomerAMG solver used as mgr coarse solver for the pressure reduced system
    setPressureAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRCOMPOSITIONALMULTIPHASEHYBRIDFVM_HPP_*/
