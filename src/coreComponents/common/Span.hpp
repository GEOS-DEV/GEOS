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
 * @file Span.hpp
 */

#ifndef GEOSX_CODINGUTILITIES_SPAN_HPP
#define GEOSX_CODINGUTILITIES_SPAN_HPP

#include "codingUtilities/traits.hpp"
#include "common/Logger.hpp"

#include <memory>
#include <iterator>

namespace geosx
{

/**
 * @brief Lightweight non-owning wrapper over a contiguous range of elements.
 * @tparam T type of range element
 *
 * This is a simplified version of C++17's std::span<> that doesn't support
 * compile-time extents.
 */
template< typename T >
class Span
{
public:

  /// Type of range element
  using value_type = T;

  /// Type used for indexing the range
  using size_type = std::ptrdiff_t;

  /**
   * @brief Construct an empty span.
   */
  constexpr Span() noexcept = default;

  /**
   * @brief Construct a span from pointer and size.
   * @param ptr pointer to beginning of contiguous range
   * @param size number of elements
   */
  Span( T * const ptr, size_type const size ) noexcept
    : m_data( ptr ),
    m_size( size )
  {
    GEOSX_ASSERT_GE( size, 0 );
  }

  /**
   * @brief Construct a span from pair of iterators.
   * @tparam ITER type of iterators
   * @param begin iterator to start of the range
   * @param end iterator to end of the range
   */
  template< typename ITER >
  constexpr Span( ITER const begin, ITER const end ) noexcept
    : Span( std::addressof( *begin ), std::distance( begin, end ) )
  {}

  /**
   * @brief Construct a span from a c-array.
   * @tparam N compile-time size of array
   * @param arr reference to array
   */
  template< int N >
  constexpr Span( T (& arr)[N] ) noexcept
    : Span( std::addressof( arr[0] ), N )
  {}

  /**
   * @brief Construct a span from a range-like object (anything that has begin() and end()).
   * @tparam R the range type
   * @param range reference to the range
   * @note The user must guarantee that the range is contiguous in memory.
   */
  template< typename R, typename std::enable_if_t< traits::is_range_like< R > > * = nullptr >
  constexpr Span( R const & range )
    : Span( range.begin(), range.end() )
  {}

  /**
   * @brief @return size of the range
   */
  constexpr size_type size() const noexcept
  {
    return m_size;
  }

  /**
   * @brief @return @p true iff range is empty
   */
  constexpr bool empty() const noexcept
  {
    return m_size == 0;
  }

  /**
   * @brief @return raw pointer to the data
   */
  constexpr T * data() const noexcept
  {
    return m_data;
  }

  /**
   * @brief @return iterator to start of the range
   */
  constexpr T * begin() const noexcept
  {
    return m_data;
  }

  /**
   * @brief @return iterator past-the-end of the range
   */
  constexpr T * end() const noexcept
  {
    return m_data + m_size;
  }

  /**
   * @brief @return reverse iterator to start or the range
   */
  constexpr std::reverse_iterator< T * > rbegin() const noexcept
  {
    return std::reverse_iterator< T * >( m_data + m_size - 1 );
  }

  /**
   * @brief @return reverse iterator to end of the range
   */
  constexpr std::reverse_iterator< T * > rend() const noexcept
  {
    return std::reverse_iterator< T * >( m_data - 1 );
  }

  /**
   * @brief @return reference to the first element
   * @pre !empty()
   */
  T & front() const
  {
    GEOSX_ASSERT_GT( m_size, 0 );
    return m_data[0];
  }

  /**
   * @brief @return reference to the last element
   * @pre !empty()
   */
  T & back() const
  {
    GEOSX_ASSERT_GT( m_size, 0 );
    return m_data[m_size-1];
  }

  /**
   * @brief @return reference to @p i-th element
   * @param i element index
   * @pre 0 <= i < size()
   */
  T & operator[]( size_type const i ) const
  {
    GEOSX_ASSERT_GE( i, 0 );
    GEOSX_ASSERT_GT( m_size, i );
    return m_data[i];
  }

private:

  /// Pointer to contiguous range of elements
  T * m_data{};

  /// Number of elements represented by this span
  size_type m_size{};

};

}

#endif //GEOSX_CODINGUTILITIES_SPAN_HPP
