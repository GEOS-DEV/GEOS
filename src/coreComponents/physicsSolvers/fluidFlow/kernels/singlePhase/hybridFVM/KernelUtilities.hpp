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
 * @file SinglePhaseHybridFVMKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVMKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVMKERNELS_HPP


namespace geos
{

namespace singlePhaseHybridFVMKernels
{

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
    case 7:
    { lambda( std::integral_constant< T, 7 >() ); return;}
    case 8:
    { lambda( std::integral_constant< T, 8 >() ); return;}
    case 9:
    { lambda( std::integral_constant< T, 9 >() ); return;}
    case 10:
    { lambda( std::integral_constant< T, 10 >() ); return;}
    case 11:
    { lambda( std::integral_constant< T, 11 >() ); return;}
    case 12:
    { lambda( std::integral_constant< T, 12 >() ); return;}
    case 13:
    { lambda( std::integral_constant< T, 13 >() ); return;}
    default: GEOS_ERROR( "Unknown numFacesInElem value: " << value );
  }
}

} // namespace internal

} // namespace singlePhaseHybridFVMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVMKERNELS_HPP
