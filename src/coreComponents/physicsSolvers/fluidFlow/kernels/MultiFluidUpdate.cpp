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
 * @file MultiFluidUpdate.cpp
 */

#include "MultiFluidUpdate.hpp"
#include "constitutive/fluid/multifluid/MultiFluidSelector.hpp"

namespace geos
{

void MultiFluidUpdate::update( constitutive::MultiFluidBase & fluid,
                               localIndex const size,
                               arrayView1d< real64 const > const & pressure,
                               arrayView1d< real64 const > const & temperature,
                               arrayView2d< real64 const, compflow::USD_COMP > const & composition )
{
  constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    using FluidType = TYPEOFREF( castedFluid );
    typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

    fluidUpdate< FluidType >( fluidWrapper,
                              size,
                              pressure,
                              temperature,
                              composition );
  } );
}

void MultiFluidUpdate::update( constitutive::MultiFluidBase & fluid,
                               SortedArrayView< localIndex const > const & targetSet,
                               arrayView1d< real64 const > const & pressure,
                               arrayView1d< real64 const > const & temperature,
                               arrayView2d< real64 const, compflow::USD_COMP > const & composition )
{
  constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    using FluidType = TYPEOFREF( castedFluid );
    typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

    fluidUpdate< FluidType >( fluidWrapper,
                              targetSet,
                              pressure,
                              temperature,
                              composition );
  } );
}

} // namespace geos
