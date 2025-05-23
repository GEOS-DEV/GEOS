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
 * @file FluidUpdateKernel.cpp
 */

#include "FluidUpdateKernel.hpp"
#include "constitutive/fluid/multifluid/MultiFluidSelector.hpp"

namespace geos
{
namespace thermalCompositionalMultiphaseBaseKernels
{

void FluidUpdate::update( localIndex const size,
                          constitutive::MultiFluidBase & fluid,
                          arrayView1d< real64 const > const & pres,
                          arrayView1d< real64 const > const & temp,
                          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac )
{
  constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castFluid )
  {
    using FluidType = TYPEOFREF( castFluid );
    using FluidUpdateType = typename FluidType::KernelWrapper;
    FluidUpdateType fluidWrapper = castFluid.createKernelWrapper();
    FluidUpdateKernel< typename FluidType::exec_policy, FluidUpdateType >::launch( size,
                                                                                   fluidWrapper,
                                                                                   pres,
                                                                                   temp,
                                                                                   compFrac );

  } );
}

template< typename POLICY >
void FluidUpdate::localUpdate( localIndex const size,
                               constitutive::MultiFluidBase & fluid,
                               arrayView1d< real64 const > const & pres,
                               arrayView1d< real64 const > const & temp,
                               arrayView2d< real64 const, compflow::USD_COMP > const & compFrac )
{
  constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castFluid )
  {
    using FluidType = TYPEOFREF( castFluid );
    using FluidUpdateType = typename FluidType::KernelWrapper;
    FluidUpdateType fluidWrapper = castFluid.createKernelWrapper();
    FluidUpdateKernel< POLICY, FluidUpdateType >::launch( size,
                                                          fluidWrapper,
                                                          pres,
                                                          temp,
                                                          compFrac );

  } );
}

void FluidUpdate::update( SortedArrayView< localIndex const > const & targetSet,
                          constitutive::MultiFluidBase & fluid,
                          arrayView1d< real64 const > const & pres,
                          arrayView1d< real64 const > const & temp,
                          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac )
{
  constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castFluid )
  {
    using FluidType = TYPEOFREF( castFluid );
    using FluidUpdateType = typename FluidType::KernelWrapper;
    FluidUpdateType fluidWrapper = castFluid.createKernelWrapper();

    FluidUpdateKernel< typename FluidType::exec_policy, FluidUpdateType >::launch( targetSet,
                                                                                   fluidWrapper,
                                                                                   pres,
                                                                                   temp,
                                                                                   compFrac );
  } );
}

template
void FluidUpdate::localUpdate< serialPolicy >( localIndex const,
                                               constitutive::MultiFluidBase &,
                                               arrayView1d< real64 const > const &,
                                               arrayView1d< real64 const > const &,
                                               arrayView2d< real64 const, compflow::USD_COMP > const & );

} // namespace thermalCompositionalMultiphaseBaseKernels
} // namespace geos
