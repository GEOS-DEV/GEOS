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


/**
 * @brief Wrapper class for umpire allocator, only used to determine which umpire allocator to use based on
 * availability.
 * @note It would be better in some ways to have this be a singleton simply because the underlying umpire allocators
 *       are also singletons. It currently isn't mostly because we use the default constructor in the buffer_type
 *       instantiations in the neighbor comm (which can be refactored to use a singleton instance-getter fairly easily).
 *       In general though it would mean any trivial usage of the buffer_type this is used in more cumbersome.
 *       At the moment our usage looks like:
 *         buffer_type(size)
 *       but if this was a singleton it would necessitate that we create buffers as: (outside of neighbor comm)
 *         buffer_type(size,buffer_allocator<buffer_type>::getInstance())
 *       which is messy. Changing the buffer_type can fix the issue but would introduce additional refactoring so for
 *       the moment this implementation suffices.
 */
template< typename T >
class buffer_allocator
{
public:
  // The type used to instantiate the class, and the underlying umpire allocator.
  typedef T value_type;
private:
  // An umpire allocator allocating the type for which this class is instantiated.
  static bool prefer_pinned;
  umpire::TypedAllocator< value_type > alloc;
  bool prefer_pinned_l;
public:
  /**
   * @brief Changed the default behavior of the class to either prefer using pinned memory
   *        ( if a pinned memory allocator is available ), or to prefer host memory.
   *        This option can be changed at any time and will not effect the behavior of existing
   *        buffer_allocator objects.
   * @note It would be preferable in some ways to have a
   *       factory returning pointers to the requester allocator but the buffer_allocator
   *       needs to have the default constructor usable to allow the vector and vector of vectors
   *       containing the buffers to be built without needing additional reworking.
   *       If we decide we want a richer set of runtime decisions other than just pinned/host
   *       that refactoring would likely be worthwhile.
   */
  static void preferPinned( bool p ) { prefer_pinned = p; }
  /**
   * @brief Default behavior is to allocate host memory, if there is a pinned memory allocator
   *        provided by umpire for the target platform, use that instead.
   */
  buffer_allocator()
    : alloc( umpire::TypedAllocator< value_type >( umpire::ResourceManager::getInstance().getAllocator( umpire::resource::Host )))
    , prefer_pinned_l( prefer_pinned )
  {
    auto & rm = umpire::ResourceManager::getInstance();

    if( rm.isAllocator( "PINNED" ) && prefer_pinned_l )
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

template< typename T >
bool buffer_allocator< T >::prefer_pinned = false;

#endif

#endif
