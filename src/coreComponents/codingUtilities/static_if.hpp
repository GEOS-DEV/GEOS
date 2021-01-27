/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file static_if.hpp
 */

#ifndef GEOSX_CODINGUTILITIES_STATIC_IF_HPP_
#define GEOSX_CODINGUTILITIES_STATIC_IF_HPP_

// UNCRUSTIFY-OFF

namespace geosx
{

template<bool CONDITION>
struct static_if_wrapper
{
  template<typename LAMBDA_BODY>
  constexpr inline static void ifFunction(LAMBDA_BODY&&) { }
};

template<>
struct static_if_wrapper<true>
{
  template<typename LAMBDA_BODY>
  constexpr inline static void ifFunction(LAMBDA_BODY&& lambda)
  {
    lambda();
  }
};

/// static_if is sometimes used to contain code that is only valid on host
/// host_device_static_if can be used when the contained code is valid both on
///                       host and device
template <bool COND>
struct host_device_static_if
{
  template <typename LAMBDA_BODY>
  constexpr inline GEOSX_HOST_DEVICE static void ifFunction(LAMBDA_BODY &&) { }
};

template <>
struct host_device_static_if<true>
{
  template <typename LAMBDA_BODY>
  constexpr inline GEOSX_HOST_DEVICE static void ifFunction(LAMBDA_BODY && lambda)
  {
    lambda();
  }
};

}

#define static_if( condition ) \
  geosx::static_if_wrapper<condition>::ifFunction( [&] () -> void

#define static_if_host_device( condition ) \
  geosx::host_device_static_if<condition>::ifFunction( [&] () -> void

#define end_static_if );

// UNCRUSTIFY-ON

#endif /* GEOSX_CODINGUTILITIES_STATIC_IF_HPP_ */
