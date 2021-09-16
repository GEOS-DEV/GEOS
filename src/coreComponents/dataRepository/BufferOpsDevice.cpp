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

#include "BufferOpsDevice.hpp"

namespace geosx
{

namespace bufferOps
{

template< bool DO_PACKING, typename T >
GEOSX_HOST_DEVICE
localIndex
PackPointerDataDevice( buffer_unit_type * & buffer,
                       T const * const GEOSX_RESTRICT var,
                       localIndex const length )
{
  localIndex const sizeOfPackedChars = length * sizeof( T );
  if( DO_PACKING )
  {
    memcpy( buffer, var, length * sizeof( T ) );
    buffer += length * sizeof( T );
  }
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename T >
GEOSX_HOST_DEVICE
localIndex
PackPointerDevice( buffer_unit_type * & buffer,
                   T const * const GEOSX_RESTRICT var,
                   localIndex const length )
{
  localIndex sizeOfPackedChars = sizeof( localIndex );
  if( DO_PACKING )
  {
    memcpy( buffer, &length, sizeof( localIndex ) );
    buffer += sizeof( localIndex );
  }
  sizeOfPackedChars += PackPointerDataDevice< DO_PACKING >( buffer, var, length );
  return sizeOfPackedChars;
}

template< typename T >
GEOSX_HOST_DEVICE
localIndex
UnpackPointerDataDevice( buffer_unit_type const * & buffer,
                         T * const GEOSX_RESTRICT var,
                         localIndex const expectedLength )
{
  localIndex sizeOfUnpackedChars = expectedLength * sizeof(T);
  memcpy( var, buffer, expectedLength * sizeof(T) );
  buffer += expectedLength * sizeof(T);
  return sizeOfUnpackedChars;
}

template< typename T >
GEOSX_HOST_DEVICE
localIndex
UnpackPointerDevice( buffer_unit_type const * & buffer,
                     T * const GEOSX_RESTRICT var,
                     localIndex const expectedLength )
{
  localIndex length = 0;
  localIndex sizeOfUnpackedChars = sizeof( localIndex );
  memcpy( &length, buffer, sizeof( localIndex ) );
  buffer += sizeof( localIndex );
  GEOSX_ASSERT_EQ( length, expectedLength );
  GEOSX_DEBUG_VAR( expectedLength );
  sizeOfUnpackedChars += UnpackPointerDataDevice( buffer, var, length );
  return sizeOfUnpackedChars;
}

template< bool DO_PACKING, typename T, int NDIM, int USD >
typename std::enable_if< can_memcpy< T >, localIndex >::type
PackDataDevice( buffer_unit_type * & buffer,
                ArrayView< T const, NDIM, USD > const & var,
                parallelDeviceEvents & events )
{
  if( DO_PACKING )
  {
    parallelDeviceStream stream;
    events.emplace_back( forAll< parallelDeviceAsyncPolicy<> >( stream, var.size(), [=] GEOSX_DEVICE ( localIndex ii )
    {
      reinterpret_cast< std::remove_const_t< T > * >( buffer )[ ii ] = var.data()[ ii ];
    } ) );
  }
  localIndex packedSize = var.size() * sizeof(T);
  if( DO_PACKING )
  {
    buffer += var.size() * sizeof(T);
  }
  return packedSize;
}

template< bool DO_PACKING, typename T, int NDIM, int USD >
typename std::enable_if< can_memcpy< T >, localIndex >::type
PackDevice( buffer_unit_type * & buffer,
            ArrayView< T const, NDIM, USD > const & var,
            parallelDeviceEvents & events )
{
  localIndex packedSize = PackPointerDevice< DO_PACKING >( buffer, var.dims(), NDIM );
  packedSize += PackPointerDevice< DO_PACKING >( buffer, var.strides(), NDIM );
  packedSize += PackDataDevice< DO_PACKING >( buffer, var, events );
  return packedSize;
}

template< typename T, int NDIM, int USD >
typename std::enable_if< can_memcpy< T >, localIndex >::type
UnpackDataDevice( buffer_unit_type const * & buffer,
                  ArrayView< T, NDIM, USD > const & var,
                  parallelDeviceEvents & events )
{
  parallelDeviceStream stream;
  events.emplace_back( forAll< parallelDeviceAsyncPolicy<> >( stream, var.size(), [=] GEOSX_DEVICE ( localIndex ii )
  {
    var.data()[ ii ] = reinterpret_cast< const T * >( buffer )[ ii ];
  } ) );
  localIndex packedSize = var.size() * sizeof(T);
  buffer += var.size() * sizeof(T);
  return packedSize;
}

template< typename T, int NDIM, int USD >
typename std::enable_if< can_memcpy< T >, localIndex >::type
UnpackDevice( buffer_unit_type const * & buffer,
              ArrayView< T, NDIM, USD > const & var,
              parallelDeviceEvents & events )
{
  localIndex dims[NDIM];
  localIndex packedSize = UnpackPointerDevice( buffer, dims, NDIM );
  localIndex strides[NDIM];
  packedSize += UnpackPointerDevice( buffer, strides, NDIM );
  for( int dd = 0; dd < NDIM; ++dd )
  {
    GEOSX_ERROR_IF_NE( strides[dd], var.strides()[dd] );
  }
  packedSize += UnpackDataDevice( buffer, var, events );
  return packedSize;
}

template< bool DO_PACKING, typename T, int NDIM, int USD, typename T_INDICES >
typename std::enable_if< can_memcpy< T >, localIndex >::type
PackDataByIndexDevice ( buffer_unit_type * & buffer,
                        ArrayView< T const, NDIM, USD > const & var,
                        const T_INDICES & indices,
                        parallelDeviceEvents & events )
{
  localIndex const numIndices = indices.size();
  localIndex const sliceSize = var.size() / var.size( 0 );
  localIndex packedSize = numIndices * sliceSize * sizeof( T );
  if( DO_PACKING )
  {
    T * devBuffer = reinterpret_cast< T * >( buffer );
    parallelDeviceStream stream;
    events.emplace_back( forAll< parallelDevicePolicy< > >( stream, numIndices, [=] GEOSX_DEVICE ( localIndex const ii )
    {
      T * threadBuffer = &devBuffer[ ii * sliceSize ];
      LvArray::forValuesInSlice( var[ indices[ ii ] ], [&] GEOSX_DEVICE ( T const & value )
      {
        *threadBuffer = value;
        ++threadBuffer;
      } );
    } ) );

    buffer += packedSize;
  }

  return packedSize;
}

template< bool DO_PACKING, typename T, int NDIM, int USD, typename T_INDICES >
typename std::enable_if< can_memcpy< T >, localIndex >::type
PackByIndexDevice ( buffer_unit_type * & buffer,
                    ArrayView< T const, NDIM, USD > const & var,
                    const T_INDICES & indices,
                    parallelDeviceEvents & events )
{
  localIndex packedSize = PackPointerDevice< DO_PACKING >( buffer, var.strides(), NDIM );
  size_t typeSize = sizeof( T );
  uintptr_t misalignment = 0;
  if( DO_PACKING )
  {
    uintptr_t address = reinterpret_cast< uintptr_t >( buffer );
    misalignment = address % typeSize;
    if( misalignment != 0 )
    {
      buffer += typeSize - misalignment;
    }
  }
  packedSize += typeSize - 1;
  packedSize += PackDataByIndexDevice< DO_PACKING >( buffer, var, indices, events );
  return packedSize;
}

template< typename T, int NDIM, int USD, typename T_INDICES >
typename std::enable_if< can_memcpy< T >, localIndex >::type
UnpackDataByIndexDevice ( buffer_unit_type const * & buffer,
                          ArrayView< T, NDIM, USD > const & var,
                          T_INDICES const & indices,
                          parallelDeviceEvents & events )
{
  localIndex numIndices = indices.size();
  localIndex sliceSize = var.size() / var.size( 0 );
  localIndex unpackSize = numIndices * sliceSize * sizeof( T );
  T const * devBuffer = reinterpret_cast< T const * >( buffer );
  parallelDeviceStream stream;
  events.emplace_back( forAll< parallelDeviceAsyncPolicy<> >( stream, numIndices, [=] GEOSX_DEVICE ( localIndex const ii )
  {
    T const * threadBuffer = &devBuffer[ ii * sliceSize ];
    LvArray::forValuesInSlice( var[ indices[ ii ] ], [&threadBuffer] GEOSX_DEVICE ( T & value )
    {
      value = *threadBuffer;
      ++threadBuffer;
    } );
  } ) );

  buffer += unpackSize;
  return unpackSize;
}

template< typename T, int NDIM, int USD, typename T_INDICES >
typename std::enable_if< can_memcpy< T >, localIndex >::type
UnpackByIndexDevice ( buffer_unit_type const * & buffer,
                      ArrayView< T, NDIM, USD > const & var,
                      T_INDICES const & indices,
                      parallelDeviceEvents & events )
{
  size_t typeSize = sizeof( T );
  localIndex strides[NDIM];
  localIndex unpackSize = UnpackPointerDevice( buffer, strides, NDIM );
  uintptr_t address = reinterpret_cast< uintptr_t >( buffer );
  uintptr_t misalignment = address % typeSize;
  if( misalignment != 0 )
  {
    buffer += typeSize - misalignment;
  }
  unpackSize += typeSize - 1;
  unpackSize += UnpackDataByIndexDevice( buffer, var, indices, events );
  return unpackSize;
}

#define DECLARE_PACK_UNPACK( TYPE, NDIM, USD ) \
  template localIndex PackDevice< true, TYPE, NDIM, USD > \
    ( buffer_unit_type * &buffer, \
    ArrayView< TYPE const, NDIM, USD > const & var, \
    parallelDeviceEvents & events ); \
  template localIndex PackDevice< false, TYPE, NDIM, USD > \
    ( buffer_unit_type * &buffer, \
    ArrayView< TYPE const, NDIM, USD > const & var, \
    parallelDeviceEvents & events ); \
  template localIndex UnpackDevice< TYPE, NDIM, USD > \
    ( buffer_unit_type const * & buffer, \
    ArrayView< TYPE, NDIM, USD > const & var, \
    parallelDeviceEvents & events ); \
  template localIndex PackByIndexDevice< true, TYPE, NDIM, USD > \
    ( buffer_unit_type * &buffer, \
    ArrayView< TYPE const, NDIM, USD > const & var, \
    arrayView1d< const localIndex > const & indices, \
    parallelDeviceEvents & events ); \
  template localIndex PackByIndexDevice< false, TYPE, NDIM, USD > \
    ( buffer_unit_type * &buffer, \
    ArrayView< TYPE const, NDIM, USD > const & var, \
    arrayView1d< const localIndex > const & indices, \
    parallelDeviceEvents & events ); \
  template localIndex UnpackByIndexDevice< TYPE, NDIM, USD > \
    ( buffer_unit_type const * & buffer, \
    ArrayView< TYPE, NDIM, USD > const & var, \
    arrayView1d< const localIndex > const & indices, \
    parallelDeviceEvents & events ); \
    \
  template localIndex PackDataDevice< true, TYPE, NDIM, USD > \
    ( buffer_unit_type * &buffer, \
    ArrayView< TYPE const, NDIM, USD > const & var, \
    parallelDeviceEvents & events ); \
  template localIndex PackDataDevice< false, TYPE, NDIM, USD > \
    ( buffer_unit_type * &buffer, \
    ArrayView< TYPE const, NDIM, USD > const & var, \
    parallelDeviceEvents & events ); \
  template localIndex UnpackDataDevice< TYPE, NDIM, USD > \
    ( buffer_unit_type const * & buffer, \
    ArrayView< TYPE, NDIM, USD > const & var, \
    parallelDeviceEvents & events ); \
  template localIndex PackDataByIndexDevice< true, TYPE, NDIM, USD > \
    ( buffer_unit_type * &buffer, \
    ArrayView< TYPE const, NDIM, USD > const & var, \
    arrayView1d< const localIndex > const & indices, \
    parallelDeviceEvents & events ); \
  template localIndex PackDataByIndexDevice< false, TYPE, NDIM, USD > \
    ( buffer_unit_type * &buffer, \
    ArrayView< TYPE const, NDIM, USD > const & var, \
    arrayView1d< const localIndex > const & indices, \
    parallelDeviceEvents & events ); \
  template localIndex UnpackDataByIndexDevice< TYPE, NDIM, USD > \
    ( buffer_unit_type const * & buffer, \
    ArrayView< TYPE, NDIM, USD > const & var, \
    arrayView1d< const localIndex > const & indices, \
    parallelDeviceEvents & events )

#define DECLARE_PACK_UNPACK_UP_TO_2D( TYPE ) \
  DECLARE_PACK_UNPACK( TYPE, 1, 0 ); \
  DECLARE_PACK_UNPACK( TYPE, 2, 0 ); \
  DECLARE_PACK_UNPACK( TYPE, 2, 1 )

#define DECLARE_PACK_UNPACK_UP_TO_3D( TYPE ) \
  DECLARE_PACK_UNPACK_UP_TO_2D( TYPE ); \
  DECLARE_PACK_UNPACK( TYPE, 3, 0 ); \
  DECLARE_PACK_UNPACK( TYPE, 3, 1 ); \
  DECLARE_PACK_UNPACK( TYPE, 3, 2 )

#define DECLARE_PACK_UNPACK_UP_TO_4D( TYPE ) \
  DECLARE_PACK_UNPACK_UP_TO_3D( TYPE ); \
  DECLARE_PACK_UNPACK( TYPE, 4, 0 ); \
  DECLARE_PACK_UNPACK( TYPE, 4, 1 ); \
  DECLARE_PACK_UNPACK( TYPE, 4, 2 ); \
  DECLARE_PACK_UNPACK( TYPE, 4, 3 )

#define DECLARE_PACK_UNPACK_UP_TO_5D( TYPE ) \
  DECLARE_PACK_UNPACK_UP_TO_4D( TYPE ); \
  DECLARE_PACK_UNPACK( TYPE, 5, 0 ); \
  DECLARE_PACK_UNPACK( TYPE, 5, 1 ); \
  DECLARE_PACK_UNPACK( TYPE, 5, 2 ); \
  DECLARE_PACK_UNPACK( TYPE, 5, 3 ); \
  DECLARE_PACK_UNPACK( TYPE, 5, 4 )

DECLARE_PACK_UNPACK_UP_TO_3D( int );
DECLARE_PACK_UNPACK_UP_TO_3D( long int );
DECLARE_PACK_UNPACK_UP_TO_3D( long long int );
DECLARE_PACK_UNPACK_UP_TO_3D( real32 );
DECLARE_PACK_UNPACK_UP_TO_5D( real64 );
DECLARE_PACK_UNPACK_UP_TO_3D( R1Tensor );

} // namespace bufferOps

} // namespace geosx
