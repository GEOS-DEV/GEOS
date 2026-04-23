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

#ifndef GEOS_FRICTIONDRIVERRUNTEST_HPP_
#define GEOS_FRICTIONDRIVERRUNTEST_HPP_

#include "constitutiveDrivers/contact/FrictionDriver.hpp"
#include "physicsSolvers/solidMechanics/contact/FractureState.hpp"
#include "constitutive/solid/SolidFields.hpp"


namespace geos
{


template< typename FRICTION_TYPE >
void
FrictionDriver::runTest( FRICTION_TYPE & friction,
                         const arrayView2d< real64 > & table )
{

  array2d< real64 > jumps, tractions;
  jumps.resize( table.size( 0 ), 3 );
  tractions.resize( table.size( 0 ), 3 );

  for( integer n = 0; n < table.size( 0 ); ++n )
  {
    jumps[n][0] = table( n, NJUMP );
    jumps[n][1] = table( n, SLIP0 );
    jumps[n][2] = table( n, SLIP1 );

    tractions[n][0] = table( n, NTRAC );
    tractions[n][1] = table( n, STRAC0 );
    tractions[n][2] = table( n, STRAC1 );

  }


  // create kernel wrapper
  typename FRICTION_TYPE::KernelWrapper const kernelWrapper = friction.createKernelUpdates();

  forAll< parallelDevicePolicy<> >( 1,
                                    [ kernelWrapper, table, jumps, tractions ]
                                    GEOS_HOST_DEVICE ( integer const GEOS_UNUSED_PARAM( ei ) )
  {

    for( integer i = 1; i < table.size( 0 ); ++i )
    {
      integer fs = fields::contact::FractureState::Stick;
      kernelWrapper.updateFractureState( jumps[i],
                                         tractions[i],
                                         fs );

      table( i, FS ) = fs;
    }
  } );

}

}


#endif //GEOS_FRICTIONDRIVERRUNTEST_HPP_
