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
 * @file CompositionalMultiphaseHybridFVMKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMKERNELS_HPP


namespace geos
{

namespace compositionalMultiphaseHybridFVMKernels
{

using namespace constitutive;

// struct to specify local and neighbor derivatives
struct Pos
{
  static constexpr integer LOCAL = 0;
  static constexpr integer NEIGHBOR = 1;
};

/******************************** Kernel switches ********************************/

namespace internal
{

template< typename T, typename LAMBDA >
void kernelLaunchSelectorFaceSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "KernelLaunchSelectorFaceSwitch: type should be integral" );

  switch( value )
  {
    case 4:
    { lambda( std::integral_constant< T, 4 >() ); return;}
    case 5:
    { lambda( std::integral_constant< T, 5 >() ); return;}
    case 6:
    { lambda( std::integral_constant< T, 6 >() ); return;}
    default: GEOS_ERROR( "Unknown numFacesInElem value: " << value );
  }
}

} // namespace internal

template< typename KERNELWRAPPER, typename INNER_PRODUCT, typename ... ARGS >
void simpleKernelLaunchSelector( localIndex numFacesInElem, ARGS && ... args )
{
  internal::kernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NUM_FACES )
  {
    KERNELWRAPPER::template launch< INNER_PRODUCT, NUM_FACES() >( std::forward< ARGS >( args )... );
  } );
}


template< typename KERNELWRAPPER, typename IP_TYPE, typename ... ARGS >
void kernelLaunchSelector( integer numFacesInElem, integer numComps, integer numPhases, ARGS && ... args )
{
  // Ideally this would be inside the dispatch, but it breaks on Summit with GCC 9.1.0 and CUDA 11.0.3.
  if( numPhases == 2 )
  {
    if( numComps == 2 )
    {
      internal::kernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NF )
      { KERNELWRAPPER::template launch< NF(), 2, 2, IP_TYPE >( std::forward< ARGS >( args )... ); } );
    }
    else if( numComps == 3 )
    {
      internal::kernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NF )
      { KERNELWRAPPER::template launch< NF(), 3, 2, IP_TYPE >( std::forward< ARGS >( args )... ); } );
    }
    else if( numComps == 4 )
    {
      internal::kernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NF )
      { KERNELWRAPPER::template launch< NF(), 4, 2, IP_TYPE >( std::forward< ARGS >( args )... ); } );
    }
    else if( numComps == 5 )
    {
      internal::kernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NF )
      { KERNELWRAPPER::template launch< NF(), 5, 2, IP_TYPE >( std::forward< ARGS >( args )... ); } );
    }
    else
    {
      GEOS_ERROR( "Unsupported number of components: " << numComps );
    }
  }
  else if( numPhases == 3 )
  {
    if( numComps == 2 )
    {
      internal::kernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NF )
      { KERNELWRAPPER::template launch< NF(), 2, 3, IP_TYPE >( std::forward< ARGS >( args )... ); } );
    }
    else if( numComps == 3 )
    {
      internal::kernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NF )
      { KERNELWRAPPER::template launch< NF(), 3, 3, IP_TYPE >( std::forward< ARGS >( args )... ); } );
    }
    else if( numComps == 4 )
    {
      internal::kernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NF )
      { KERNELWRAPPER::template launch< NF(), 4, 3, IP_TYPE >( std::forward< ARGS >( args )... ); } );
    }
    else if( numComps == 5 )
    {
      internal::kernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NF )
      { KERNELWRAPPER::template launch< NF(), 5, 3, IP_TYPE >( std::forward< ARGS >( args )... ); } );
    }
    else
    {
      GEOS_ERROR( "Unsupported number of components: " << numComps );
    }
  }
  else
  {
    GEOS_ERROR( "Unsupported number of phases: " << numPhases );
  }
}

} // namespace compositionalMultiphaseHybridFVMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMKERNELS_HPP
