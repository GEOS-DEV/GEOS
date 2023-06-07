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

namespace thermalSinglePhaseBaseKernels
{

/******************************** FluidUpdateKernel ********************************/

struct FluidUpdateKernel
{
  template< typename FLUID_WRAPPER >
  static void launch( FLUID_WRAPPER const & fluidWrapper,
                      geos::arrayView1d< geos::real64 const > const & pres,
                      geos::arrayView1d< geos::real64 const > const & temp )
  {
    geos::forAll< geos::parallelDevicePolicy<> >( fluidWrapper.numElems(), [=] GEOS_HOST_DEVICE ( geos::localIndex const k )
    {
      for( geos::localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.update( k, q, pres[k], temp[k] );
      }
    } );
  }
};

}

}

#endif //GEOSX_FLUIDUPDATEKERNEL_HPP
