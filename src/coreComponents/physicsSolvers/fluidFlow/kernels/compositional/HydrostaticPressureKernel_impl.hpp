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
 * @file HydrostaticPressureKernel_impl.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_KERNELS_COMPOSITIONAL_HYDROSTATICPRESSUREKERNEL_IMPL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_KERNELS_COMPOSITIONAL_HYDROSTATICPRESSUREKERNEL_IMPL_HPP

#include "HydrostaticPressureKernel.hpp"

namespace geos
{
namespace isothermalCompositionalMultiphaseBaseKernels
{

template< typename FLUID_WRAPPER >
void HydrostaticPressureFluidUpdate< FLUID_WRAPPER >::getMassDensity( integer const numPhases,
                                                                      FLUID_WRAPPER fluidWrapper,
                                                                      real64 const pressure,
                                                                      real64 const temperature,
                                                                      arraySlice1d< real64 const, compflow::USD_COMP - 1 > const composition,
                                                                      arraySlice1d< real64 > const massDensity,
                                                                      HydrostaticPressureStackVariables & stackVars )
{
  constitutive::MultiFluidBase::KernelWrapper::computeValues( fluidWrapper,
                                                              pressure,
                                                              temperature,
                                                              composition,
                                                              stackVars.phaseFraction[0][0],
                                                              stackVars.phaseDensity[0][0],
                                                              stackVars.phaseMassDensity[0][0],
                                                              stackVars.phaseViscocity[0][0],
                                                              stackVars.phaseEnthalpy[0][0],
                                                              stackVars.phaseInternalEnergy[0][0],
                                                              stackVars.phaseComposition[0][0],
                                                              stackVars.totalDensity[0][0] );

  for( integer ip = 0; ip < numPhases; ++ip )
  {
    massDensity[ip] = stackVars.phaseMassDensity[0][0][ip];
  }
}

} // namespace isothermalCompositionalMultiphaseBaseKernels
} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_KERNELS_COMPOSITIONAL_HYDROSTATICPRESSUREKERNEL_IMPL_HPP
