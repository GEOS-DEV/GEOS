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
 * @file RelativePermeabilityUpdateKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_RELATIVEPERMEABILITYUPDATEKERNEL_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_RELATIVEPERMEABILITYUPDATEKERNEL_HPP


namespace geos
{

namespace isothermalCompositionalMultiphaseBaseKernels
{

/******************************** RelativePermeabilityUpdateKernel ********************************/

struct RelativePermeabilityUpdateKernel
{
  template< typename POLICY, typename RELPERM_WRAPPER >
  static void
  launch( geos::localIndex const size,
          RELPERM_WRAPPER const & relPermWrapper,
          arrayView2d< geos::real64 const, geos::compflow::USD_PHASE > const & phaseVolFrac )
  {
    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < relPermWrapper.numGauss(); ++q )
      {
        relPermWrapper.update( k, q, phaseVolFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename RELPERM_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          RELPERM_WRAPPER const & relPermWrapper,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < relPermWrapper.numGauss(); ++q )
      {
        relPermWrapper.update( k, q, phaseVolFrac[k] );
      }
    } );
  }
};

}

}

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_RELATIVEPERMEABILITYUPDATEKERNEL_HPP
