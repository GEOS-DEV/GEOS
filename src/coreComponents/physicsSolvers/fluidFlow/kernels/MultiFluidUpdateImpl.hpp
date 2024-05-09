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
 * @file MultiFluidUpdateImpl.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_KERNELS_MULTIFLUIDUPDATEIMPL_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_KERNELS_MULTIFLUIDUPDATEIMPL_HPP_

#include "MultiFluidUpdate.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "physicsSolvers/fluidFlow/ThermalCompositionalMultiphaseBaseKernels.hpp"

namespace geos
{

template< typename FLUID_WRAPPER >
GEOS_HOST_DEVICE
void MultiFluidUpdate::KernelWrapper< FLUID_WRAPPER >::update( FLUID_WRAPPER const & fluidWrapper,
                                                               real64 const & pressure,
                                                               real64 const & temperature,
                                                               arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                                                               arraySlice1d< real64, constitutive::multifluid::USD_PHASE - 2 > const & phaseFraction,
                                                               arraySlice1d< real64, constitutive::multifluid::USD_PHASE - 2 > const & phaseDensity,
                                                               arraySlice1d< real64, constitutive::multifluid::USD_PHASE - 2 > const & phaseMassDensity,
                                                               arraySlice1d< real64, constitutive::multifluid::USD_PHASE - 2 > const & phaseViscosity,
                                                               arraySlice1d< real64, constitutive::multifluid::USD_PHASE - 2 > const & phaseEnthalpy,
                                                               arraySlice1d< real64, constitutive::multifluid::USD_PHASE - 2 > const & phaseInternalEnergy,
                                                               arraySlice2d< real64, constitutive::multifluid::USD_PHASE_COMP-2 > const & phaseCompFraction,
                                                               real64 & totalDensity )
{
  constitutive::MultiFluidBase::KernelWrapper::computeValues( fluidWrapper,
                                                              pressure,
                                                              temperature,
                                                              composition,
                                                              phaseFraction,
                                                              phaseDensity,
                                                              phaseMassDensity,
                                                              phaseViscosity,
                                                              phaseEnthalpy,
                                                              phaseInternalEnergy,
                                                              phaseCompFraction,
                                                              totalDensity );
}

template< typename FLUID_TYPE >
GEOS_HOST_DEVICE
void MultiFluidUpdate::Updater< FLUID_TYPE >::update( typename FLUID_TYPE::KernelWrapper const & fluidWrapper,
                                                      localIndex const index,
                                                      integer const node,
                                                      real64 const & pressure,
                                                      real64 const & temperature,
                                                      arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition )
{
  fluidWrapper.update( index, node, pressure, temperature, composition );
}

template< typename FLUID_TYPE >
void MultiFluidUpdate::Updater< FLUID_TYPE >::update( typename FLUID_TYPE::KernelWrapper const & fluidWrapper,
                                                      localIndex const size,
                                                      arrayView1d< real64 const > const & pressure,
                                                      arrayView1d< real64 const > const & temperature,
                                                      arrayView2d< real64 const, compflow::USD_COMP > const & composition )
{
  thermalCompositionalMultiphaseBaseKernels::FluidUpdateKernel::
    launch< typename FLUID_TYPE::exec_policy >( size,
                                                fluidWrapper,
                                                pressure,
                                                temperature,
                                                composition );
}

template< typename FLUID_TYPE >
void MultiFluidUpdate::Updater< FLUID_TYPE >::update( typename FLUID_TYPE::KernelWrapper const & fluidWrapper,
                                                      SortedArrayView< localIndex const > const & targetSet,
                                                      arrayView1d< real64 const > const & pressure,
                                                      arrayView1d< real64 const > const & temperature,
                                                      arrayView2d< real64 const, compflow::USD_COMP > const & composition )
{
  thermalCompositionalMultiphaseBaseKernels::FluidUpdateKernel::
    launch< typename FLUID_TYPE::exec_policy >( targetSet,
                                                fluidWrapper,
                                                pressure,
                                                temperature,
                                                composition );
}

} //namespace geos

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_KERNELS_MULTIFLUIDUPDATEIMPL_HPP_
