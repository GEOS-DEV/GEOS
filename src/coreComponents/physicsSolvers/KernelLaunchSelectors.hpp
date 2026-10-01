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
 * @file KernelLaunchSelectors.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_KERNELLAUNCHSELECTORS_HPP
#define GEOS_PHYSICSSOLVERS_KERNELLAUNCHSELECTORS_HPP

#include "common/GeosxConfig.hpp"
#include "common/logger/Logger.hpp"

#include <type_traits>
#include <utility>

namespace geos
{
namespace internal
{

template< typename S, typename T, typename LAMBDA >
void invokePhaseDispatchLambda ( S val, T numPhases, LAMBDA && lambda )
{
  if( numPhases == 1 )
  {
    lambda( val, std::integral_constant< T, 1 >());
    return;
  }
  else if( numPhases == 2 )
  {
    lambda( val, std::integral_constant< T, 2 >());
    return;
  }
  else if( numPhases == 3 )
  {
    lambda( val, std::integral_constant< T, 3 >());
    return;
  }
  else
  {
    GEOS_ERROR( GEOS_FMT( "Unsupported state: {}", numPhases ) );
  }
}

template< typename S, typename T, typename LAMBDA >
void invokeThermalDispatchLambda ( S val, T isThermal, LAMBDA && lambda )
{
  if( isThermal == 1 )
  {
    lambda( val, std::integral_constant< T, 1 >());
    return;
  }
  else if( isThermal == 0 )
  {
    lambda( val, std::integral_constant< T, 0 >());
    return;
  }
  else
  {
    GEOS_ERROR( GEOS_FMT( "Unsupported state: {}", isThermal ) );
  }
}

template< typename T, typename LAMBDA >
void kernelLaunchSelectorThermalSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "kernelLaunchSelectorThermalSwitch: type should be integral" );

  switch( value )
  {
    case 0:
    {
      lambda( std::integral_constant< T, 0 >() );
      return;
    }
    case 1:
    {
      lambda( std::integral_constant< T, 1 >() );
      return;
    }
    default:
    {
      GEOS_ERROR( GEOS_FMT( "Unsupported thermal state: {}", value ) );
    }
  }
}

/**
 * @brief Dispatch a runtime component count to a compile-time constant.
 * @note Counts above GEOS_MAX_COMPONENTS are not instantiated.
 */
template< typename T, typename LAMBDA >
void kernelLaunchSelectorCompSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "kernelLaunchSelectorCompSwitch: type should be integral" );
  switch( value )
  {
    #define GEOS_DISPATCH_COMPONENT( NC ) \
      case NC: \
      { lambda( std::integral_constant< T, NC >() ); return; \
      }
    GEOS_FOR_EACH_COMPONENT( GEOS_DISPATCH_COMPONENT )
#undef GEOS_DISPATCH_COMPONENT
    default:
    {
      GEOS_ERROR( GEOS_FMT( "Unsupported number of components: {}. This build instantiates 1 through {} (GEOS_MAX_COMPONENTS).",
                            value, GEOS_MAX_COMPONENTS ) );
    }
  }
}

template< typename T, typename LAMBDA >
void kernelLaunchSelectorCompThermSwitch( T value, bool const isThermal, LAMBDA && lambda )
{
  kernelLaunchSelectorCompSwitch( value, [&] ( auto NC )
  {
    invokeThermalDispatchLambda( NC, isThermal, lambda );
  } );
}

template< typename T, typename LAMBDA >
void kernelLaunchSelectorCompPhaseSwitch( T value, T n_phase, LAMBDA && lambda )
{
  kernelLaunchSelectorCompSwitch( value, [&] ( auto NC )
  {
    invokePhaseDispatchLambda( NC, n_phase, lambda );
  } );
}


} // end namspace internal
} // end namespace geos


#endif // GEOS_PHYSICSSOLVERS_KERNELLAUNCHSELECTORS_HPP
