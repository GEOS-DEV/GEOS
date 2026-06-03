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
 * @file ThermalSinglePhasePoromechanicsReservoirFVM.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRTHERMALSINGLEPHASEPOROMECHANICSRESERVOIRFVM_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRTHERMALSINGLEPHASEPOROMECHANICSRESERVOIRFVM_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief ThermalSinglePhasePoromechanicsReservoirFVM strategy.
 *
 * dofLabel: 0 = displacement, x-component
 * dofLabel: 1 = displacement, y-component
 * dofLabel: 2 = displacement, z-component
 * dofLabel: 3 = pressure
 * dofLabel: 4 = temperature
 * dofLabel: 5 = well pressure
 * dofLabel: 6 = well rate
 * dofLabel: 7 = well temperature
 *
 * 2-level MGR reduction strategy
 *   - 1st level: eliminate displacements (0,1,2)
 *   - 2nd level: eliminate the thermal well block
 *   - The coarse grid (pressure and temperature) is solved with BoomerAMG.
 *
 */
class ThermalSinglePhasePoromechanicsReservoirFVM : public MGRStrategyBase< 2 >
{
public:

  /**
   * @brief Constructor.
   */
  explicit ThermalSinglePhasePoromechanicsReservoirFVM( arrayView1d< int const > const & )
    : MGRStrategyBase( 8 )
  {
    // Level 0: eliminate displacements
    m_labels[0].push_back( 3 );
    m_labels[0].push_back( 4 );
    m_labels[0].push_back( 5 );
    m_labels[0].push_back( 6 );
    m_labels[0].push_back( 7 );

    // Level 1: eliminate the well block
    m_labels[1].push_back( 3 );
    m_labels[1].push_back( 4 );

    setupLabels();

    // Level 0
    m_levelFRelaxType[0]         = MGRFRelaxationType::amgVCycle;
    m_levelFRelaxIters[0]        = 1;
    m_levelInterpType[0]         = MGRInterpolationType::jacobi;
    m_levelRestrictType[0]       = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]   = MGRCoarseGridMethod::nonGalerkin;
    m_levelGlobalSmootherType[0] = MGRGlobalSmootherType::none;

    // Level 1
    m_levelFRelaxType[1]         = MGRFRelaxationType::gsElimWInverse;
    m_levelFRelaxIters[1]        = 1;
    m_levelInterpType[1]         = MGRInterpolationType::blockJacobi;
    m_levelRestrictType[1]       = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1]   = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[1] = MGRGlobalSmootherType::none;
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
    // If the wells are shut, Gaussian elimination on the well block is unnecessary.
    if( mgrParams.areWellsShut )
    {
      m_levelFRelaxType[1] = MGRFRelaxationType::jacobi;
    }

    setReduction( precond, mgrData );

    // Configure the BoomerAMG solver used as F-relaxation for the first level.
    setMechanicsFSolver( precond, mgrData, mgrParams.separateComponents );

    // Configure the BoomerAMG solver used as coarse solver for the pressure/temperature system.
    setPressureTemperatureAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRTHERMALSINGLEPHASEPOROMECHANICSRESERVOIRFVM_HPP_*/
