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

template< typename FLUID_TYPE >
void MultiFluidUpdate::fluidUpdate( typename FLUID_TYPE::KernelWrapper const & fluidWrapper,
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
void MultiFluidUpdate::fluidUpdate( typename FLUID_TYPE::KernelWrapper const & fluidWrapper,
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
