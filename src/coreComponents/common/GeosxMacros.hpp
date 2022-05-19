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
 * @file GeosxMacros.hpp
 *
 * This file contains various macro definitions.
 */

#include "common/GeosxConfig.hpp"
#include "LvArray/src/Macros.hpp"

#ifndef GEOSX_COMMON_GEOSXMACROS_HPP_
#define GEOSX_COMMON_GEOSXMACROS_HPP_

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
  #define GEOSX_DEVICE_COMPILE
#endif


/**
 * @name Host-device markers
 *
 * These macros are used to denote host/device/inline functions in a compiler-specific way.
 * They must be prepended to a function or lambda declaration/definition.
 * They will be defined differently when compiled by e.g. a CUDA compiler.
 */
///@{

#if defined(GEOSX_USE_DEVICE)
#define GEOSX_HOST __host__
#define GEOSX_DEVICE __device__
#define GEOSX_HOST_DEVICE __host__ __device__
#define GEOSX_FORCE_INLINE __forceinline__
#define PRAGMA_UNROLL _Pragma("unroll")
#else
/// Marks a host-only function.
#define GEOSX_HOST
/// Marks a device-only function.
#define GEOSX_DEVICE
/// Marks a host-device function.
#define GEOSX_HOST_DEVICE
/// Marks a function or lambda for inlining
#define GEOSX_FORCE_INLINE inline
/// Compiler directive specifying to unroll the loop.
#define PRAGMA_UNROLL
#endif

///@}

/**
 * @name Unused variable markers.
 *
 * These macros are used to explicitly mark a variable/argument as unused
 * and thus silence compiler warnings.
 */
///@{

/// Mark an unused argument and silence compiler warnings.
#define GEOSX_UNUSED_PARAM( X )

/**
 * @brief Used to silence unused variable warnings, cuda doesn't respect casting to void.
 * @tparam ARGS argument types
 * @param ...
 */
template< typename ... ARGS >
GEOSX_HOST_DEVICE inline constexpr
void i_g_n_o_r_e( ARGS const & ... ) {}

/// Mark an unused variable and silence compiler warnings.
#define GEOSX_UNUSED_VAR( ... ) i_g_n_o_r_e( __VA_ARGS__ )

/// Mark a debug variable and silence compiler warnings.
#define GEOSX_DEBUG_VAR( ... ) GEOSX_UNUSED_VAR( __VA_ARGS__ )

///@}

#if defined(GEOSX_USE_OPENMP)
/// Wrap a pragma clause in the _Pragma statement. We seek to make this include the omp portion of the clause.
#define PRAGMA_OMP( clause ) _Pragma( clause )
//  #define PRAGMA_OMP( clause ) _Pragma( STRINGIZE( omp clause ) )
#else
/// No-op version of PRAGMA_OMP
#define PRAGMA_OMP( clause )
#endif

/// preprocessor variable for the C99 restrict keyword for use with pointers
#define GEOSX_RESTRICT LVARRAY_RESTRICT

/// preprocessor variable for the C99 restrict keyword for use with the "this" pointer
#define GEOSX_RESTRICT_THIS LVARRAY_RESTRICT_THIS

/// Doxygen can't parse a `decltype( auto )` return type, using this gets around that.
#define GEOSX_DECLTYPE_AUTO_RETURN decltype( auto )

/// Macro to concatenate two tokens (low level)
#define GEOSX_CONCAT_IMPL( A, B ) A ## B

/// Macro to concatenate two tokens (user level)
#define GEOSX_CONCAT( A, B ) GEOSX_CONCAT_IMPL( A, B )

#endif // GEOSX_COMMON_GEOSXMACROS_HPP_
