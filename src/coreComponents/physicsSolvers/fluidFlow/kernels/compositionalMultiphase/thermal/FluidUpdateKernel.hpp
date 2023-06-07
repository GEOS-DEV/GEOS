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
 * @file FluidUpdateKernel.hpp
 */

#ifndef GEOSX_FLUIDUPDATEKERNEL_HPP
#define GEOSX_FLUIDUPDATEKERNEL_HPP


namespace geos
{

namespace thermalCompositionalMultiphaseBaseKernels
{

/******************************** FluidUpdateKernel ********************************/

struct FluidUpdateKernel
{
  template< typename POLICY, typename FLUID_WRAPPER >
  static void
  launch( geos::localIndex const size,
          FLUID_WRAPPER const & fluidWrapper,
          geos::arrayView1d< geos::real64 const > const & pres,
          geos::arrayView1d< geos::real64 const > const & temp,
          arrayView2d< geos::real64 const, geos::compflow::USD_COMP > const & compFrac )
  {
    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.update( k, q, pres[k], temp[k], compFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename FLUID_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          FLUID_WRAPPER const & fluidWrapper,
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
};

}

}

#endif //GEOSX_FLUIDUPDATEKERNEL_HPP
