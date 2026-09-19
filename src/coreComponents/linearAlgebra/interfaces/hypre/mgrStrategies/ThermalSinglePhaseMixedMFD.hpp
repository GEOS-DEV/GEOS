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
 * @file ThermalSinglePhaseMixedMFD.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRTHERMALSINGLEPHASEMIXEDMFD_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRTHERMALSINGLEPHASEMIXEDMFD_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief MGR strategy of the thermal mixed mimetic single-phase solver.
 *
 * Unknowns: face mass and heat fluxes, cell pressure, temperature and enthalpy. Level 0 eliminates
 * every face flux as in the isothermal strategy (symmetric Gauss-Seidel on the flux blocks, Jacobi
 * interpolation from their diagonals); level 1 eliminates the enthalpy, whose closure rows
 * h - h_eos(p, T) = 0 have a unit diagonal, exactly; the coarse solve is one BoomerAMG V-cycle on
 * the pressure-temperature system.
 *
 * Point markers provided by the solver (LinearSolverParameters::MGR::customPointMarkers):
 *  0 = condensed face mass flux, 1 = saddle-point face mass flux,
 *  2 = cell pressure, 3 = cell temperature, 4 = cell enthalpy,
 *  5 = condensed face heat flux, 6 = saddle-point face heat flux
 */
class ThermalSinglePhaseMixedMFD : public MGRStrategyBase< 2 >
{
public:
  /**
   * @brief Constructor.
   */
  explicit ThermalSinglePhaseMixedMFD( arrayView1d< int const > const & )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( 7 ) )
  {
    // Level 0: eliminate every face flux, keep the cell unknowns
    m_labels[0].push_back( 2 );
    m_labels[0].push_back( 3 );
    m_labels[0].push_back( 4 );
    // Level 1: eliminate the enthalpy, keep pressure and temperature
    m_labels[1].push_back( 2 );
    m_labels[1].push_back( 3 );
    setupLabels();

    m_levelFRelaxType[0]         = MGRFRelaxationType::hybridSymmetricGaussSeidel;
    m_levelFRelaxIters[0]        = 1;
    m_levelInterpType[0]         = MGRInterpolationType::jacobi;
    m_levelRestrictType[0]       = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]   = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[0] = MGRGlobalSmootherType::none;

    m_levelFRelaxType[1]         = MGRFRelaxationType::jacobi;
    m_levelFRelaxIters[1]        = 1;
    m_levelInterpType[1]         = MGRInterpolationType::jacobi;
    m_levelRestrictType[1]       = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1]   = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[1] = MGRGlobalSmootherType::none;
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

    // one BoomerAMG V-cycle on the pressure-temperature Schur complement, without aggressive coarsening.
    // C/F Jacobi relaxation: independent of the partition, unlike the hybrid Gauss-Seidel default
    BoomerAMGParameters amgParameters = pressureTemperatureAMGParameters();
    amgParameters.aggressiveNumLevels = 0;
    amgParameters.relaxType = hypre::getAMGRelaxationType( LinearSolverParameters::AMG::SmootherType::jacobi );
    amgParameters.numSweeps = 2;
    configureBoomerAMG( mgrData.coarseSolver, amgParameters );
    mgrData.coarseSolver.setup = HYPRE_BoomerAMGSetup;
    mgrData.coarseSolver.solve = HYPRE_BoomerAMGSolve;
    mgrData.coarseSolver.destroy = HYPRE_BoomerAMGDestroy;
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRTHERMALSINGLEPHASEMIXEDMFD_HPP_*/
