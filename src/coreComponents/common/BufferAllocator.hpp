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
 * @file BufferAllocator.hpp
 */

#ifndef GEOSX_BUFFER_ALLOCATOR_HPP
#define GEOSX_BUFFER_ALLOCATOR_HPP

#ifdef USE_CHAI
#include <umpire/ResourceManager.hpp>
#include <umpire/TypedAllocator.hpp>

/// Wrapper class for umpire allocator, only used to determine which umpire allocator to use based on availability.
template< typename T >
class buffer_allocator
{
public:
  // The type used to instantiate the class, and the underlying umpire allocator.
  typedef T value_type;
private:
  // An umpire allocator allocating the type for which this class is instantiated.
  umpire::TypedAllocator< value_type > alloc;
public:
  /**
   * @brief Default behavior is to allocate host memory, if there is a pinned memory allocator
   *        provided by umpire for the target platform, use that instead.
   */
  buffer_allocator()
    : alloc( umpire::TypedAllocator< value_type >( umpire::ResourceManager::getInstance().getAllocator( umpire::resource::Host )))
  {
    auto & rm = umpire::ResourceManager::getInstance();
    if( rm.isAllocator( "PINNED" ) )
      alloc = umpire::TypedAllocator< value_type >( rm.getAllocator( umpire::resource::Pinned ));
  }
  /**
   * @brief Allocate a buffer.
   * @param sz The number of elements of type value_type to allocate a buffer for.
   * @return A pointer to the allocated buffer.
   */
  value_type * allocate( size_t sz ) { return alloc.allocate( sz ); }
  /**
   * @brief Deallocate a buffer.
   * @param buffer A pointer to the buffer to deallocate
   * @param sz The size of the buffer to deallocate.
   */
  void deallocate( value_type * buffer, size_t sz )
  {
    if( buffer != nullptr )
      alloc.deallocate( buffer, sz );
  }
  /**
   * @brief Inequality operator.
   * @param The other buffer_allocator to test against this buffer allocator for inequality.
   * @return Always false. Since the actual umpire allocator is a singleton, so any properly-typed
   *         buffer can be deallocated from any properly-typed buffer_allocator.
   */
  bool operator!=( const buffer_allocator & ) { return false; }
  /**
   * @brief Equality operator.
   * @param other The other buffer_allocator to test against this buffer allocator for equality.
   * return Always true. Since the actual umpire allocator is a singleton, so any properly-typed
   *        buffer can be deallocated from any properly-typed buffer_allocator.
   */
  bool operator==( const buffer_allocator & other ) { return !operator!=( other ); }
};

#endif

#endif
