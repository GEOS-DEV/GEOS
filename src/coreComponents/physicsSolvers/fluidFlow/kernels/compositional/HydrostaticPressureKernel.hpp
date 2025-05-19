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
 * @file HydrostaticPressureKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_KERNELS_COMPOSITIONAL_HYDROSTATICPRESSUREKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_KERNELS_COMPOSITIONAL_HYDROSTATICPRESSUREKERNEL_HPP

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "functions/TableFunction.hpp"

namespace geos
{
namespace constitutive
{
class MultiFluidBase;
}
namespace isothermalCompositionalMultiphaseBaseKernels
{

/******************************** HydrostaticPressureKernel ********************************/

struct HydrostaticPressureStackVariables
{
  static integer constexpr maxNumPhase = constitutive::MultiFluidConstants::MAX_NUM_PHASES;
  static integer constexpr maxNumComp = constitutive::MultiFluidConstants::MAX_NUM_COMPONENTS;

  HydrostaticPressureStackVariables( integer const numPhases, integer const numComps ):
    phaseFraction( 1, 1, numPhases ),
    phaseDensity( 1, 1, numPhases ),
    phaseMassDensity( 1, 1, numPhases ),
    phaseViscocity( 1, 1, numPhases ),
    phaseEnthalpy( 1, 1, numPhases ),
    phaseInternalEnergy( 1, 1, numPhases ),
    phaseComposition( 1, 1, numPhases, numComps ),
    totalDensity( 1 )
  {}

  constitutive::MultiFluidBase::PhaseProp::StackValueType< maxNumPhase > phaseFraction;
  constitutive::MultiFluidBase::PhaseProp::StackValueType< maxNumPhase > phaseDensity;
  constitutive::MultiFluidBase::PhaseProp::StackValueType< maxNumPhase > phaseMassDensity;
  constitutive::MultiFluidBase::PhaseProp::StackValueType< maxNumPhase > phaseViscocity;
  constitutive::MultiFluidBase::PhaseProp::StackValueType< maxNumPhase > phaseEnthalpy;
  constitutive::MultiFluidBase::PhaseProp::StackValueType< maxNumPhase > phaseInternalEnergy;
  constitutive::MultiFluidBase::PhaseComp::StackValueType< maxNumPhase * maxNumComp > phaseComposition;
  constitutive::MultiFluidBase::FluidProp::StackValueType< 1 > totalDensity;
};

/**
 * @brief Wrapper for the fluid query
 * @tparam FLUID_WRAPPER The type of fluid wrapper
 */
template< typename FLUID_WRAPPER >
struct HydrostaticPressureFluidUpdate
{
  /**
   * @brief Calculate the fluid mass density at specified conditions
   * @param[in] numPhases -the number of phases
   * @param[in] fluidWrapper - a wrapper to the fluid
   * @param[in] pressure - the pressure
   * @param[in] temperature - the temperature
   * @param[in] composition - the total composition
   * @param[out] massDensity - the mass density
   * @param[out] stackVars - workspace for fluid calculation
   */
  static void getMassDensity( integer const numPhases,
                              FLUID_WRAPPER fluidWrapper,
                              real64 const pressure,
                              real64 const temperature,
                              arraySlice1d< real64 const, compflow::USD_COMP - 1 > const composition,
                              arraySlice1d< real64 > const massDensity,
                              HydrostaticPressureStackVariables & stackVars );
};

struct HydrostaticPressureKernel
{

  /**
   * @brief A return type for the hydrostatic pressure computation
   */
  enum class ReturnType : integer
  {
    FAILED_TO_CONVERGE = 0,
    DETECTED_MULTIPHASE_FLOW = 1,
    SUCCESS = 2
  };

  /**
   * @brief Perform hydrostatic pressure initialisation
   * @details This kernel is purely serial
   */
  static ReturnType
  launch( localIndex const size,
          integer const numComps,
          integer const numPhases,
          integer const ipInit,
          integer const maxNumEquilIterations,
          real64 const equilTolerance,
          real64 const (&gravVector)[ 3 ],
          real64 const & minElevation,
          real64 const & elevationIncrement,
          real64 const & datumElevation,
          real64 const & datumPres,
          constitutive::MultiFluidBase & fluid,
          arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappers,
          TableFunction::KernelWrapper tempTableWrapper,
          arrayView1d< arrayView1d< real64 > const > elevationValues,
          arrayView1d< real64 > pressureValues );
};

} // namespace isothermalCompositionalMultiphaseBaseKernels
} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_KERNELS_COMPOSITIONAL_HYDROSTATICPRESSUREKERNEL_HPP
