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
 * @file SolidInternalEnergyUpdateKernel.hpp
 */

#ifndef GEOSX_SOLIDINTERNALENERGYUPDATEKERNEL_HPP
#define GEOSX_SOLIDINTERNALENERGYUPDATEKERNEL_HPP


namespace geos
{

namespace thermalSinglePhaseBaseKernels
{

/******************************** SolidInternalEnergyUpdateKernel ********************************/

struct SolidInternalEnergyUpdateKernel
{

  template< typename POLICY, typename SOLID_INTERNAL_ENERGY_WRAPPER >
  static void
  launch( geos::localIndex const size,
          SOLID_INTERNAL_ENERGY_WRAPPER const & solidInternalEnergyWrapper,
          geos::arrayView1d< geos::real64 const > const & temp )
  {
    geos::forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( geos::localIndex const k )
    {
      solidInternalEnergyWrapper.update( k, temp[k] );
    } );
  }
};

}

}

#endif //GEOSX_SOLIDINTERNALENERGYUPDATEKERNEL_HPP
