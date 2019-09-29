/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file static_if.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CODINGUTILITIES_STATIC_IF_HPP_
#define SRC_COMPONENTS_CORE_SRC_CODINGUTILITIES_STATIC_IF_HPP_


namespace geosx
{

template<bool CONDITION>
struct static_if_wrapper
{
  template<typename LAMBDA_BODY>
  constexpr inline static void if_function(LAMBDA_BODY&&) { }
};

template<>
struct static_if_wrapper<true>
{
  template<typename LAMBDA_BODY>
  constexpr inline static void if_function(LAMBDA_BODY&& lambda)
  {
    lambda();
  }
};

}

#define static_if( condition ) \
  geosx::static_if_wrapper<condition>::if_function( [&] () -> void

#define end_static_if );

#endif /* SRC_COMPONENTS_CORE_SRC_CODINGUTILITIES_STATIC_IF_HPP_ */
