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
    m_storage.resizeWithoutInitializationOrDestruction( space, maxEntries, valuesPerEntry );
  }
  
  bool full() const
  {
    return m_end - m_begin == m_storage.size( 0 );
  }

  ArraySlice1D front() const {
    return m_storage[m_begin];
    //return ArraySlice1D<T const>(&m_storage[m_begin][0], m_storage.size(1), 1);
  }

  ArraySlice1D back() const {
    return m_storage[m_end];
    //return ArraySlice1D<T const>(&m_storage[m_end][0], m_storage.size(1), 1);
  }

  void pop_back() {
    assert( m_end > m_begin );
    m_end--;
    m_end = m_end % m_storage.size(0);
  }

  void pop_front() {
    assert( m_end > m_begin );
    m_begin++;
    m_begin = m_begin % m_storage.size(0);
  }

  void emplace_back( const ArraySlice1D & src ) {
    LvArray::memcpy( m_storage[++m_end], src);
    m_end = m_end % m_storage.size(0);
  }

  void emplace_front( const ArraySlice1D & src ) {
    if ( m_begin == 0 ) m_begin = m_storage.size(0);
    LvArray::memcpy( m_storage[--m_begin], src );
  }
  
private:
  Array2D m_storage;
  IndexType m_begin = 0;
  IndexType m_end = 0;
};
}
#endif // FIXEDSIZEDEQUE_HPP
