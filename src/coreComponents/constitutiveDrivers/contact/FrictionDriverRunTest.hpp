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

#include "constitutive/contact/FrictionDriver.hpp"
#include "physicsSolvers/solidMechanics/contact/FractureState.hpp"
#include "constitutive/solid/SolidFields.hpp"


namespace geos
{


template< typename FRICTION_TYPE >
void
FrictionDriver::runTest( FRICTION_TYPE & relperm,
                        const arrayView2d< real64 > & table )
{
  // get number of phases and components

  integer const numPhases = relperm.numFluidPhases();

  // create kernel wrapper

  typename RELPERM_TYPE::KernelWrapper const kernelWrapper = relperm.createKernelWrapper();

  // set saturation to user specified feed
  // it is more convenient to provide input in molar, so perform molar to mass conversion here

  array2d< real64, compflow::LAYOUT_PHASE > saturationValues;
  if( numPhases > 2 )
  {
    saturationValues.resize(( m_numSteps + 1 ) * ( m_numSteps + 1 ), numPhases );
  }
  else
  {
    saturationValues.resize( m_numSteps + 1, numPhases );
  }
  integer const ipWater = relperm.getPhaseOrder()[PT::WATER];
  integer const ipOil = relperm.getPhaseOrder()[PT::OIL];
  integer const ipGas = relperm.getPhaseOrder()[PT::GAS];
  const localIndex offset = std::max( std::max( ipOil, ipWater ), std::max( ipOil, ipGas ) ) + 1;

  for( integer n = 0; n < table.size( 0 ); ++n )
  {


    if( m_numPhases > 2 )
    {
      saturationValues[n][ipWater] = table( n, ipWater + 1 );
      saturationValues[n][ipOil] = table( n, ipOil + 1 );
      saturationValues[n][ipGas] = table( n, ipGas + 1 );
    }
    else//two-phase
    {
      if( ipWater < 0 )
      {
        saturationValues[n][ipOil] = table( n, ipOil + 1 );
        saturationValues[n][ipGas] = table( n, ipGas + 1 );
      }
      else if( ipGas < 0 )
      {
        saturationValues[n][ipWater] = table( n, ipWater + 1 );
        saturationValues[n][ipOil] = table( n, ipOil + 1 );
      }
    }

  }

  arrayView2d< real64 const, compflow::USD_PHASE > const saturation = saturationValues.toViewConst();

  // perform relperm update using table (Swet,Snonwet) and save resulting total density, etc.
  // note: column indexing should be kept consistent with output file header below.

  forAll< parallelDevicePolicy<> >( saturation.size( 0 ),
                                    [numPhases, kernelWrapper, saturation, table,
                                     offset] GEOS_HOST_DEVICE ( integer const n )
  {
    kernelWrapper.update( 0, 0, saturation[n] );
    for( integer p = 0; p < numPhases; ++p )
    {
      table( n, offset + 1 + p ) = kernelWrapper.relperm()( 0, 0, p );
    }
  } );

}


}


#endif //GEOS_FRICTIONDRIVERRUNTEST_HPP_
