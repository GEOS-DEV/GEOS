//
// Created by root on 10/24/22.
//

#ifndef GEOSX_RELPERMDRIVERRUNTEST_HPP_
#define GEOSX_RELPERMDRIVERRUNTEST_HPP_

#include "RelpermDriver.hpp"
#include "layouts.hpp"


namespace geosx
{

using namespace constitutive;

//specific to Hysteresis
template< typename RELPERM_TYPE >
std::enable_if_t< std::is_same< TableRelativePermeabilityHysteresis, RELPERM_TYPE >::value, void >
RelpermDriver::runTest( RELPERM_TYPE & relperm,
                        arrayView3d< real64 > const & table )
{
  // get number of phases and components
  integer const numPhases = relperm.numFluidPhases();

  // create kernel wrapper

  typename TableRelativePermeabilityHysteresis::KernelWrapper const kernelWrapper = relperm.createKernelWrapper();

  // set saturation to user specified feed
  // it is more convenient to provide input in molar, so perform molar to mass conversion here

  array3d< real64, relperm::LAYOUT_RELPERM > saturationValues;
  if( numPhases > 2 )
  {
    saturationValues.resize( 1, ( m_numSteps + 1 ) * ( m_numSteps + 1 ), numPhases );
  }
  else
  {
    saturationValues.resize( 1, m_numSteps + 1, numPhases );
  }
  using PT = typename RELPERM_TYPE::PhaseType;
  integer const ipWater = relperm.getPhaseOrder()[PT::WATER];
  integer const ipOil = relperm.getPhaseOrder()[PT::OIL];
  integer const ipGas = relperm.getPhaseOrder()[PT::GAS];
  localIndex offset = std::max( std::max( ipOil, ipWater ), std::max( ipOil, ipGas ) ) + 1;

  for( integer n = 0; n < table.size( 0 ); ++n )
  {


    if( m_numPhases > 2 )
    {
      saturationValues[0][n][ipWater] = table( n, 0, ipWater + 1 );
      saturationValues[0][n][ipOil] = table( n, 0, ipOil + 1 );
      saturationValues[0][n][ipGas] = table( n, 0, ipGas + 1 );
    }
    else//two-phase
    {
      if( ipWater < 0 )
      {
        saturationValues[0][n][ipOil] = table( n, 0, ipOil + 1 );
        saturationValues[0][n][ipGas] = table( n, 0, ipGas + 1 );
      }
      else if( ipGas < 0 )
      {
        saturationValues[0][n][ipWater] = table( n, 0, ipWater + 1 );
        saturationValues[0][n][ipOil] = table( n, 0, ipOil + 1 );
      }
      else if( ipOil < 0 )
      {
        saturationValues[0][n][ipWater] = table( n, 0, ipWater + 1 );
        saturationValues[0][n][ipGas] = table( n, 0, ipGas + 1 );
      }
    }
  }


  arrayView3d< real64 const, relperm::USD_RELPERM > const saturation = saturationValues.toViewConst();

  // perform relperm update using table (Swet,Snonwet) and save resulting relative permeabilities
  // note: column indexing should be kept consistent with output file header below.

  relperm.setMinMaxToDrainage( this );
  forAll< parallelDevicePolicy<> >( saturation.size( 0 ),
                                    [numPhases, kernelWrapper, saturation, table,
                                     offset] GEOSX_HOST_DEVICE ( localIndex const i )
  {
    for( integer n = 0; n < saturation.size( 1 ); ++n )
    {

      // nw phase set max to snw_max to get the imbibition bounding curve
      kernelWrapper.update( i, 0, saturation[0][n] );
      for( integer p = 0; p < numPhases; ++p )
      {
        table( n, 0, offset + 1 + p ) = kernelWrapper.relperm()( i, 0, p );
      }
    }
  } );

  //loop in charge of hysteresis values
  offset += numPhases;
  relperm.setMinMaxToImbibition( this );
  forAll< parallelDevicePolicy<> >( saturation.size( 0 ),
                                    [numPhases, kernelWrapper, saturation, table,
                                     offset] GEOSX_HOST_DEVICE ( localIndex const i )
  {
    for( integer n = 0; n < saturation.size( 1 ); ++n )
    {

      // nw phase set max to snw_max to get the imbibition bounding curve

      kernelWrapper.update( i, 0, saturation[0][n] );
      for( integer p = 0; p < numPhases; ++p )
      {
        table( n, 0, offset + 1 + p ) = kernelWrapper.relperm()( i, 0, p );
      }
    }
  } );


}

template< typename RELPERM_TYPE >
std::enable_if_t< !std::is_same< TableRelativePermeabilityHysteresis, RELPERM_TYPE >::value, void >
RelpermDriver::runTest( RELPERM_TYPE & relperm,
                        arrayView3d< real64 > const & table )
{
  // get number of phases and components

  integer const numPhases = relperm.numFluidPhases();

  // create kernel wrapper

  typename RELPERM_TYPE::KernelWrapper const kernelWrapper = relperm.createKernelWrapper();

  // set saturation to user specified feed
  // it is more convenient to provide input in molar, so perform molar to mass conversion here

  array3d< real64, relperm::LAYOUT_RELPERM > saturationValues;
  if( numPhases > 2 )
  {
    saturationValues.resize( 1, ( m_numSteps + 1 ) * ( m_numSteps + 1 ), numPhases );
  }
  else
  {
    saturationValues.resize( 1, m_numSteps + 1, numPhases );
  }
  using PT = typename RELPERM_TYPE::PhaseType;
  integer const ipWater = relperm.getPhaseOrder()[PT::WATER];
  integer const ipOil = relperm.getPhaseOrder()[PT::OIL];
  integer const ipGas = relperm.getPhaseOrder()[PT::GAS];
  const localIndex offset = std::max( std::max( ipOil, ipWater ), std::max( ipOil, ipGas ) ) + 1;

  for( integer n = 0; n < table.size( 0 ); ++n )
  {


    if( m_numPhases > 2 )
    {
      saturationValues[0][n][ipWater] = table( n, 0, ipWater + 1 );
      saturationValues[0][n][ipOil] = table( n, 0, ipOil + 1 );
      saturationValues[0][n][ipGas] = table( n, 0, ipGas + 1 );
    }
    else//two-phase
    {
      if( ipWater < 0 )
      {
        saturationValues[0][n][ipOil] = table( n, 0, ipOil + 1 );
        saturationValues[0][n][ipGas] = table( n, 0, ipGas + 1 );
      }
      else if( ipGas < 0 )
      {
        saturationValues[0][n][ipWater] = table( n, 0, ipWater + 1 );
        saturationValues[0][n][ipOil] = table( n, 0, ipOil + 1 );
      }
    }

  }

  arrayView3d< real64 const, relperm::USD_RELPERM > const saturation = saturationValues.toViewConst();

  // perform relperm update using table (Swet,Snonwet) and save resulting total density, etc.
  // note: column indexing should be kept consistent with output file header below.

  forAll< parallelDevicePolicy<> >( saturation.size( 0 ),
                                    [numPhases, kernelWrapper, saturation, table,
                                     offset] GEOSX_HOST_DEVICE ( localIndex const i )
  {
    for( integer n = 0; n < saturation.size( 1 ); ++n )
    {

      kernelWrapper.update( i, 0, saturation[0][n] );
      for( integer p = 0; p < numPhases; ++p )
      {
        table( n, 0, offset + 1 + p ) = kernelWrapper.relperm()( i, 0, p );
      }
    }
  } );

}


}


#endif //GEOSX_RELPERMDRIVERRUNTEST_HPP_
