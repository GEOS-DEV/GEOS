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
#ifndef FIXEDSIZEDEQUE_HPP
#define FIXEDSIZEDEQUE_HPP

#if defined( GEOSX_USE_CUDA )
#  include <cuda.h>
#  include <cuda_runtime.h>
#endif

#include <chai/ArrayManager.hpp>
#include <chai/ChaiMacros.hpp>
#include <string>
#include <functional>
#include "LvArray/src/Array.hpp"
#include "LvArray/src/ArrayView.hpp"
#include "LvArray/src/memcpy.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "common/Logger.hpp"

/// Get the positive value of a module b
#define POSITIVE_MODULO( a, b ) ( ( ( a ) % ( b ) ) + b ) % ( b )
namespace geosx
{
template< typename T, typename INDEX_TYPE >
class fixedSizeDeque
{
  /// The integer type used for indexing.
  using IndexType = INDEX_TYPE;
  /// 1D array slice
  using ArraySlice1D = LvArray::ArraySlice< T const, 1, 0, INDEX_TYPE >;
  /// 2D array type. See LvArray:Array for details.
  using Array2D = LvArray::Array< T, 2, camp::make_idx_seq_t< 2 >, IndexType, LvArray::ChaiBuffer >;
public:
  fixedSizeDeque( IndexType maxEntries, IndexType valuesPerEntry, LvArray::MemorySpace space )
  {
    GEOSX_THROW_IF( maxEntries < 0 , "Fixed sized queue size must be positive", std::runtime_error );
    GEOSX_THROW_IF( valuesPerEntry < 0 , "Fixed sized queue array size must be positive", std::runtime_error );
    m_storage.resizeWithoutInitializationOrDestruction( space, maxEntries, valuesPerEntry );
  }
  
  bool empty() const {
    return m_begin > m_end;
  }

  bool full() const
  {
    return size() == m_storage.size( 0 );
  }

  size_t size() const {
    return  (size_t)( m_end - m_begin + 1 );
  }

  size_t capacity() const {
    return  m_storage.size( 0 );
  }

  ArraySlice1D front() const {
    GEOSX_THROW_IF( empty(), "Can't get front from empty queue", std::runtime_error );
    return m_storage[ POSITIVE_MODULO( m_begin, m_storage.size( 0 ) ) ];
  }

  ArraySlice1D back() const {
    GEOSX_THROW_IF( empty(), "Can't get back from empty queue", std::runtime_error );
    return m_storage[ POSITIVE_MODULO( m_end, m_storage.size( 0 ) ) ];
  }

  void pop_back() {
    GEOSX_THROW_IF( empty(), "Can't pop back from empty queue", std::runtime_error );
    m_end--;
  }

  void inc_back() {
    GEOSX_THROW_IF( full(), "Can't increase in a full queue", std::runtime_error );
    m_end++;
  }

  void inc_front() {
    GEOSX_THROW_IF( full(), "Can't increase in a full queue", std::runtime_error );
    m_begin--;
  }

  void pop_front() {
    GEOSX_THROW_IF( empty(), "Can't pop front from empty queue", std::runtime_error );
    m_begin++;
  }

  void emplace_back( const ArraySlice1D & src ) {
    GEOSX_THROW_IF( full(), "Can't emplace in a full queue", std::runtime_error );
    ++m_end;
    LvArray::memcpy( m_storage[ POSITIVE_MODULO( m_end, m_storage.size( 0 ) ) ], src);

  }

  void emplace_front( const ArraySlice1D & src ) {
    GEOSX_THROW_IF( full(), "Can't emplace in a full  queue", std::runtime_error );
    --m_begin;
    LvArray::memcpy( m_storage[ POSITIVE_MODULO( m_begin, m_storage.size( 0 ) ) ], src );
  }
  
private:
  Array2D m_storage;
  IndexType m_begin = 0;
  IndexType m_end = -1;
};
}
#endif // FIXEDSIZEDEQUE_HPP
