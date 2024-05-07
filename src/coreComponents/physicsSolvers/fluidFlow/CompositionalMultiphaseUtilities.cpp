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
 * @file CompositionalMultiphaseUtilities.cpp
 */

#include "CompositionalMultiphaseUtilities.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidSelector.hpp"
#include "physicsSolvers/fluidFlow/ThermalCompositionalMultiphaseBaseKernels.hpp"

namespace geos
{

namespace compositionalMultiphaseUtilities
{

void updateFluidModel( constitutive::MultiFluidBase & fluid,
                       localIndex const size,
                       arrayView1d< real64 const > const & pressure,
                       arrayView1d< real64 const > const & temperature,
                       arrayView2d< real64 const, compflow::USD_COMP > const & composition )
{
  constitutiveUpdatePassThru( fluid, [&] ( auto & castFluid )
  {
    using FluidType = TYPEOFREF( castFluid );
    using ExecPolicy = typename FluidType::exec_policy;
    typename FluidType::KernelWrapper fluidWrapper = castFluid.createKernelWrapper();

    thermalCompositionalMultiphaseBaseKernels::
      FluidUpdateKernel::
      launch< ExecPolicy >( size,
                            fluidWrapper,
                            pressure,
                            temperature,
                            composition );
  } );
}

void updateFluidModel( constitutive::MultiFluidBase & fluid,
                       SortedArrayView< localIndex const > const & targetSet,
                       arrayView1d< real64 const > const & pressure,
                       arrayView1d< real64 const > const & temperature,
                       arrayView2d< real64 const, compflow::USD_COMP > const & composition )
{
  constitutiveUpdatePassThru( fluid, [&] ( auto & castFluid )
  {
    using FluidType = TYPEOFREF( castFluid );
    using ExecPolicy = typename FluidType::exec_policy;
    typename FluidType::KernelWrapper fluidWrapper = castFluid.createKernelWrapper();

    thermalCompositionalMultiphaseBaseKernels::
      FluidUpdateKernel::
      launch< ExecPolicy >( targetSet,
                            fluidWrapper,
                            pressure,
                            temperature,
                            composition );
  } );
}

} // namespace compositionalMultiphaseUtilities

} // namespace geos
