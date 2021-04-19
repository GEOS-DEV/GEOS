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

#ifndef GEOSX_COMMON_TENSOR_HPP_
#define GEOSX_COMMON_TENSOR_HPP_

#include "GeosxMacros.hpp"

namespace geosx
{

/**
 * @brief Lightweight wrapper around a c-array
 * @tparam T data type
 * @tparam SIZE_TPARAM number of values
 */
template< typename T, int SIZE_TPARAM >
class Tensor
{
public:

  static_assert( SIZE_TPARAM > 0, "Tensor size must be a positive value" );

  /// Alias for type template parameter
  using value_type = T;

  /// Alias for size template parameter
  static constexpr int SIZE = SIZE_TPARAM;

  /**
   * @brief Const element access.
   * @param i element index
   * @return const reference to the i-th element
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  T const & operator[]( std::ptrdiff_t const i ) const
  {
    return data[i];
  }

  /**
   * @brief Non-const element access.
   * @param i element index
   * @return reference to the i-th element
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  T & operator[]( std::ptrdiff_t const i )
  {
    return data[i];
  }

  /**
   * @brief Equality comparison operator.
   * @tparam U dummy parameter to enable SFINAE, do not provide explicitly
   * @param rhs the right-hand side value to compare to
   * @return true iff @code data[i] == rhs.data[i] @endcode
   * @note Version for floating point types that avoids direct equality comparison.
   */
  template< typename U = T >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  std::enable_if_t< std::is_floating_point< U >::value, bool >
  operator==( Tensor< U, SIZE > const & rhs ) const
  {
    for( int i = 0; i < SIZE; ++i )
    {
      if( (data[i] > rhs.data[i]) || (data[i] < rhs.data[i]) )
      {
        return false;
      }
    }
    return true;
  }

  /**
   * @brief Equality comparison operator.
   * @tparam U dummy parameter to enable SFINAE, do not provide explicitly
   * @param rhs the right-hand side value to compare to
   * @return true iff @code data[i] == rhs.data[i] @endcode
   * @note Version for all types except floating point.
   */
  template< typename U = T >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  std::enable_if_t< !std::is_floating_point< U >::value, bool >
  operator==( Tensor< U, SIZE > const & rhs ) const
  {
    for( int i = 0; i < SIZE; ++i )
    {
      if( data[i] != rhs.data[i] )
      {
        return false;
      }
    }
    return true;
  }

  /**
   * @brief Returns the size of the tensor
   * @param junk Unused
   * @return The value of the template parameter SIZE.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr int size( int junk ) const
  {
    GEOSX_UNUSED_VAR( junk )
    return SIZE;
  }

  /// Underlying array
  T data[SIZE] = {};

private:

  /**
   * @brief Stream insertion operator for Tensor.
   * @param os the output stream
   * @param t the tensor value
   * @return reference to @p os
   */
  friend inline std::ostream & operator<<( std::ostream & os, Tensor< T, SIZE > const & t )
  {
    os << t.data[0];
    for( int i = 1; i < SIZE; ++i )
    {
      os << ',' << t.data[i];
    }
    return os;
  }
};

}

#endif //GEOSX_COMMON_TENSOR_HPP_
