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
 * @file MathOperators.hpp
 */

#ifndef GEOS_MATH_OPERATIONS_HPP
#define GEOS_MATH_OPERATIONS_HPP

#include <variant>
#include "common/logger/Logger.hpp"
#include "common/TypesHelpers.hpp"

namespace geos
{

namespace math
{
struct Sum
{
  template< typename T >
  T operator()( T const & a, T const & b ) const
  {
    return a + b;
  }
};

struct Product
{
  template< typename T >
  T operator()( T const & a, T const & b ) const
  {
    return a * b;
  }
};

struct Difference
{
  template< typename T >
  T operator()( T const & a, T const & b ) const
  {
    return a - b;
  }
};

struct Division
{
  template< typename T >
  T operator()( T const & a, T const & b ) const
  {
    return a / b;
  }
};

using Operations = std::variant< Sum, Product, Difference, Division >;

namespace internal
{

template< typename T, typename Operation >
std::enable_if_t< std::is_arithmetic_v< T >, T >
apply_impl( Operation const & op, T const & first, T const & second )
{
  return op( first, second );
}

template< typename T, typename Operation, typename ... Args >
std::enable_if_t< std::conjunction_v< std::is_arithmetic< T >, std::is_arithmetic< Args >... >, T >
apply_impl( Operation const & op, T const & first, Args const & ... args )
{
  return op( first, apply_impl( op, args ... ));
}

template< typename Operation, typename CONTAINER_TYPE >
get_value_type_t< CONTAINER_TYPE >
apply_impl( Operation const & op, CONTAINER_TYPE const & container )
{
  auto result = container[0];
  for( std::size_t i = 1; i < container.size(); ++i )
  {
    result = op( result, container[i] );
  }
  return result;
}

} /* namespace internal */

template< typename T, typename ... Args >
std::enable_if_t< std::conjunction_v< std::is_arithmetic< T >, std::is_arithmetic< Args >... >, T >
apply( Operations const & op, T const & first, Args const & ... args )
{
  return std::visit( [&]( const auto & operation ) {
    return internal::apply_impl( operation, first, args ... );
  }, op );
}

template< typename CONTAINER_TYPE >
get_value_type_t< CONTAINER_TYPE >
apply( Operations const & op, CONTAINER_TYPE const & vector )
{
  static_assert( std::is_arithmetic_v< typename get_value_type< CONTAINER_TYPE >::type >, "The type in the container must be an arithmetic type" );
  return std::visit( [&]( const auto & operation )
  {
    return internal::apply_impl( operation, vector );
  }, op );
}

} /* namespace math */

} /* namespace geos */

#endif /* GEOS_MATH_OPERATIONS_HPP */
