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

#ifndef GEOS_COMMON_GEOSXMACROS_HPP_
#define GEOS_COMMON_GEOSXMACROS_HPP_

/**
 * @name Host-device markers
 *
 * These macros are used to denote host/device/inline functions in a compiler-specific way.
 * They must be prepended to a function or lambda declaration/definition.
 * They will be defined differently when compiled by e.g. a CUDA compiler.
 */
///@{

#if defined(GEOS_USE_DEVICE)
#define GEOS_HOST __host__
#define GEOS_DEVICE __device__
#define GEOS_HOST_DEVICE __host__ __device__
#define GEOS_FORCE_INLINE __forceinline__
#else
/// Marks a host-only function.
#define GEOS_HOST
/// Marks a device-only function.
#define GEOS_DEVICE
/// Marks a host-device function.
#define GEOS_HOST_DEVICE
/// Marks a function or lambda for inlining
#define GEOS_FORCE_INLINE inline
#endif

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
/// Macro defined when currently compiling on device (only defined in the device context).
#define GEOS_DEVICE_COMPILE
#endif

///@}


#define GEOS_PRAGMA( PRAGMA_STRING ) _Pragma( #PRAGMA_STRING )

//#if defined(__GNUC__) && !defined( __clang__) &&  __GNUC__>=11
//#define GEOS_PRAGMA_UNROLL_DEPTH( NUMITER ) GEOS_PRAGMA( GCC unroll NUMITER )
// #else
 #define GEOS_PRAGMA_UNROLL_DEPTH( NUMITER ) GEOS_PRAGMA( unroll NUMITER )
// #endif

#define GEOS_PRAGMA_UNROLL GEOS_PRAGMA_UNROLL_DEPTH( 16 )

// #if defined( GEOS_DEVICE_COMPILE )
// #define GEOS_PRAGMA_UNROLL #pragma unroll
// #else
// #define GEOS_PRAGMA_UNROLL _Pragma( "unroll" )
// #endif

#include <utility>
#include <iostream>

template<typename Lambda, std::size_t... Idx>
constexpr void static_device_for_impl(Lambda&& lambda, std::index_sequence<Idx...>)
{
    [[clang::always_inline]] (lambda(std::integral_constant<std::size_t, Idx>{}), ...);
}

template<std::size_t N, typename Lambda>
constexpr void static_device_for(Lambda&& lambda)
{
#ifdef GEOS_DEVICE_COMPILE
    static_device_for_impl(std::forward<Lambda>(lambda), std::make_index_sequence<N>{});
#else
    for( std::size_t ii = 0; ii < N; ++ii )
    {
        [[clang::always_inline]] lambda( ii );
    }
#endif
}

/**
 * @name Unused variable markers.
 *
 * These macros are used to explicitly mark a variable/argument as unused
 * and thus silence compiler warnings.
 */
///@{

/// Mark an unused argument and silence compiler warnings.
#define GEOS_UNUSED_PARAM( X )

/**
 * @brief Used to silence unused variable warnings, cuda doesn't respect casting to void.
 * @tparam ARGS argument types
 * @param ...
 */
template< typename ... ARGS >
GEOS_HOST_DEVICE inline constexpr
void i_g_n_o_r_e( ARGS const & ... ) {}

/// Mark an unused variable and silence compiler warnings.
#define GEOS_UNUSED_VAR( ... ) i_g_n_o_r_e( __VA_ARGS__ )

/// Mark a debug variable and silence compiler warnings.
#define GEOS_DEBUG_VAR( ... ) GEOS_UNUSED_VAR( __VA_ARGS__ )

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
#define GEOS_RESTRICT LVARRAY_RESTRICT

/// preprocessor variable for the C99 restrict keyword for use with the "this" pointer
#define GEOS_RESTRICT_THIS LVARRAY_RESTRICT_THIS

/// Doxygen can't parse a `decltype( auto )` return type, using this gets around that.
#define GEOS_DECLTYPE_AUTO_RETURN decltype( auto )

/// Macro to concatenate two tokens (low level)
#define GEOS_CONCAT_IMPL( A, B ) A ## B

/// Macro to concatenate two tokens (user level)
#define GEOS_CONCAT( A, B ) GEOS_CONCAT_IMPL( A, B )

/**
 * @brief [[maybe_unused]] when >= C++17, or compiler-specific implementations
 *        when < C++17
 */
#if __cplusplus >= 201703L
#define GEOS_MAYBE_UNUSED [[maybe_unused]]
#else
// If not C++17 or later, check the compiler.
    #ifdef _MSC_VER
// Microsoft Visual Studio
#define GEOS_MAYBE_UNUSED __pragma(warning(suppress: 4100))
    #elif defined(__GNUC__) || defined(__clang__)
// GCC or Clang
#define GEOS_MAYBE_UNUSED __attribute__((unused))
    #else
// If the compiler is unknown, we can't suppress the warning,
// so we define GEOS_MAYBE_UNUSED as an empty macro.
#define GEOS_MAYBE_UNUSED
    #endif
#endif

#endif // GEOS_COMMON_GEOSXMACROS_HPP_
