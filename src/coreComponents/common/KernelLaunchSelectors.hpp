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
 * @file KernelLaunchSelectors.hpp
 */

#ifndef GEOS_COMMON_KERNELLAUNCHSELECTORS_HPP
#define GEOS_COMMON_KERNELLAUNCHSELECTORS_HPP

#include "common/Logger.hpp"

namespace geos
{
namespace internal {

template< typename T, typename LAMBDA >
void kernelLaunchSelectorThermalSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "kernelLaunchSelectorThermalSwitch: type should be integral" );

  switch( value )
  {
    case 0:
    { lambda( std::integral_constant< T, 0 >() ); return; }
    case 1:
    { lambda( std::integral_constant< T, 1 >() ); return; }

    default:
    { GEOS_ERROR( "Unsupported thermal state: " << value ); }
  }
}

template< typename T, typename LAMBDA >
void kernelLaunchSelectorCompThermSwitch( T value, bool const isThermal, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "kernelLaunchSelectorCompSwitch: value type should be integral" );

  //constexpr T a = isThermal ? std::integral_constant< T, 1 >() : std::integral_constant< T, 0 >();
  if( isThermal )
  {
    switch( value )
    {
      case 1:
      {
        lambda( std::integral_constant< T, 1 >(), std::integral_constant< T, 1 >() ); return;
      }
      case 2:
      { lambda( std::integral_constant< T, 2 >(), std::integral_constant< T, 1 >() ); return; }
      case 3:
      { lambda( std::integral_constant< T, 3 >(), std::integral_constant< T, 1 >() ); return; }
      case 4:
      { lambda( std::integral_constant< T, 4 >(), std::integral_constant< T, 1 >() ); return; }
      case 5:
      { lambda( std::integral_constant< T, 5 >(), std::integral_constant< T, 1 >()); return; }
      default:
      { GEOS_ERROR( "Unsupported number of components: " << value ); }
    }
  }
  else
  {
    switch( value )
    {
      case 1:
      {
        lambda( std::integral_constant< T, 1 >(), std::integral_constant< T, 0 >() ); return;
      }
      case 2:
      { lambda( std::integral_constant< T, 2 >(), std::integral_constant< T, 0 >() ); return; }
      case 3:
      { lambda( std::integral_constant< T, 3 >(), std::integral_constant< T, 0 >() ); return; }
      case 4:
      { lambda( std::integral_constant< T, 4 >(), std::integral_constant< T, 0 >() ); return; }
      case 5:
      { lambda( std::integral_constant< T, 5 >(), std::integral_constant< T, 0 >() ); return; }
      default:
      { GEOS_ERROR( "Unsupported number of components: " << value ); }
    }
  }
}


template< typename T, typename LAMBDA >
void kernelLaunchSelectorCompPhaseSwitch( T value, T n_phase, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "kernelLaunchSelectorCompSwitch: value type should be integral" );

  //constexpr T a = isThermal ? std::integral_constant< T, 1 >() : std::integral_constant< T, 0 >();
  if ( n_phase == 1 )
  {
      switch( value )
      {
        case 1:
        {
          lambda( std::integral_constant< T, 1 >(), std::integral_constant< T, 1 >() ); return;
        }
        case 2:
        { lambda( std::integral_constant< T, 2 >(), std::integral_constant< T, 1 >() ); return; }
        case 3:
        { lambda( std::integral_constant< T, 3 >(), std::integral_constant< T, 1 >() ); return; }
        case 4:
        { lambda( std::integral_constant< T, 4 >(), std::integral_constant< T, 1 >() ); return; }
        case 5:
        { lambda( std::integral_constant< T, 5 >(), std::integral_constant< T, 1 >() ); return; }
        default:
        { GEOS_ERROR( "Unsupported number of components: " << value ); }
      }
  }
  else if ( n_phase == 2 )
  {
      switch( value )
      {
        case 1:
        {
          lambda( std::integral_constant< T, 1 >(), std::integral_constant< T, 2 >() ); return;
        }
        case 2:
        { lambda( std::integral_constant< T, 2 >(), std::integral_constant< T, 2 >() ); return; }
        case 3:
        { lambda( std::integral_constant< T, 3 >(), std::integral_constant< T, 2 >() ); return; }
        case 4:
        { lambda( std::integral_constant< T, 4 >(), std::integral_constant< T, 2 >() ); return; }
        case 5:
        { lambda( std::integral_constant< T, 5 >(), std::integral_constant< T, 2 >() ); return; }
        default:
        { GEOS_ERROR( "Unsupported number of components: " << value ); }
      }

  }
  else if ( n_phase == 3 )
  {
      switch( value )
      {
        case 1:
        {
          lambda( std::integral_constant< T, 1 >(), std::integral_constant< T, 3 >()  ); return;
        }
        case 2:
        { lambda( std::integral_constant< T, 2 >(), std::integral_constant< T, 3 >()  ); return; }
        case 3:
        { lambda( std::integral_constant< T, 3 >(), std::integral_constant< T, 3 >()  ); return; }
        case 4:
        { lambda( std::integral_constant< T, 4 >(), std::integral_constant< T, 3 >() ); return; }
        case 5:
        { lambda( std::integral_constant< T, 5 >(), std::integral_constant< T, 3 >() ); return; }
        default:
        { GEOS_ERROR( "Unsupported number of components: " << value ); }
      }
  }
  else
  {
    { GEOS_ERROR( "Unsupported number of phases: " << n_phase ); }
  }
}
} // end namspace internal
} // end namespace geos

 
#endif // GEOS_COMMON_KERNELLAUNCHSELECTORS_HPP
