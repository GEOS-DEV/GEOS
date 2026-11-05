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
 * @file FluidUpdateKernel_impl.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_FLUIDUPDATEKERNEL_IMPL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_FLUIDUPDATEKERNEL_IMPL_HPP

#include "FluidUpdateKernel.hpp"

namespace geos
{

namespace thermalCompositionalMultiphaseBaseKernels
{

template< typename POLICY, typename FLUID >
void
FluidUpdateKernel< POLICY, FLUID >::launch( localIndex const size,
                                            typename FLUID::KernelWrapper const & fluidWrapper,
                                            arrayView1d< real64 const > const & pres,
                                            arrayView1d< real64 const > const & temp,
                                            arrayView2d< real64 const, compflow::USD_COMP > const & compFrac )
{
  forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
    {
      fluidWrapper.update( k, q, pres[k], temp[k], compFrac[k] );
    }
  } );
}

template< typename POLICY, typename FLUID >
void
FluidUpdateKernel< POLICY, FLUID >::launch( SortedArrayView< localIndex const > const & targetSet,
                                            typename FLUID::KernelWrapper const & fluidWrapper,
                                            arrayView1d< real64 const > const & pres,
                                            arrayView1d< real64 const > const & temp,
                                            arrayView2d< real64 const, compflow::USD_COMP > const & compFrac )
{
  forAll< POLICY >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
  {
    localIndex const k = targetSet[a];
    for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
    {
      fluidWrapper.update( k, q, pres[k], temp[k], compFrac[k] );
    }
  } );
}

} // namespace thermalCompositionalMultiphaseBaseKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_FLUIDUPDATEKERNEL_IMPL_HPP
