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
 * @file ImmiscibleMultiphaseMFD.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRIMMISCIBLEMULTIPHASEMFD_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRIMMISCIBLEMULTIPHASEMFD_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief ImmiscibleMultiphaseMFD strategy.
 *
 * Labels description stored in point_marker_array
 *   dofLabel: 0 = cell-centered pressure
 *   dofLabel: 1 = phase volume fraction (first phase, second is eliminated by volume constraint)
 *   dofLabel: 2 = face-centered pressure (Lagrange multiplier)
 *
 * 2-level MGR reduction strategy inspired from CompositionalMultiphaseHybridFVM:
 *   - 1st level: eliminate the phase volume fractions (transport problem)
 *   - 2nd level: eliminate the cell-centered pressure
 *   - The coarse grid is the face-centered pressure (Lagrange multiplier) system
 *   - The coarse grid solved with BoomerAMG
 *   - Saturations are eliminated first using True-IMPES approach, then pressure elimination
 */
class ImmiscibleMultiphaseMFD : public MGRStrategyBase< 2 >
{
public:
  /**
   * @brief Constructor.
   * @param numComponentsPerField array with number of components for each field
   */
  explicit ImmiscibleMultiphaseMFD( arrayView1d< int const > const & numComponentsPerField )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] + numComponentsPerField[1] ) )
  {
    // Level 0: eliminate phase volume fractions (transport problem)
    m_labels[0].resize( 2 );
    m_labels[0][0] = 0;                  // cell-centered pressure
    m_labels[0][1] = m_numBlocks - 1;    // face-centered pressure (last block)

    // Level 1: eliminate cell-centered pressure (similar to SinglePhaseReservoirHybridFVM)
    m_labels[1].resize( 1 );
    m_labels[1][0] = m_numBlocks - 1;    // face-centered pressure

    setupLabels();

    // Level 0: Transport problem - eliminate phase volume fractions
    m_levelFRelaxType[0]          = MGRFRelaxationType::none;
    m_levelFRelaxIters[0]         = 1;
    m_levelInterpType[0]          = MGRInterpolationType::injection;
    m_levelRestrictType[0]        = MGRRestrictionType::blockColLumped; // True-IMPES
    m_levelCoarseGridMethod[0]    = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[0]  = MGRGlobalSmootherType::ilu0;
    m_levelGlobalSmootherIters[0] = 1;

    // Level 1: Pressure problem - eliminate cell-centered pressure (similar to SinglePhaseReservoirHybridFVM)
    m_levelFRelaxType[1]          = MGRFRelaxationType::l1jacobi;
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
  void setup( LinearSolverParameters::MGR const & mgrParams,
              HyprePrecWrapper & precond,
              HypreMGRData & mgrData )
  {
    // Similar to SinglePhaseReservoirHybridFVM: adjust F-relaxation based on well status
    // For the pressure elimination level (Level 1)
    if( mgrParams.areWellsShut )
    {
      m_levelFRelaxType[1] = MGRFRelaxationType::jacobi;
    }

    setReduction( precond, mgrData );

    // Configure the BoomerAMG solver used as mgr coarse solver for the pressure reduced system
    setPressureAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRIMMISCIBLEMULTIPHASEMFD_HPP_*/
