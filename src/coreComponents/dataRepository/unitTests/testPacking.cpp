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

#include <gtest/gtest.h>

#include "LvArray/src/Array.hpp"
#include "dataRepository/BufferOps.hpp"
#include "dataRepository/BufferOpsDevice.hpp"
#include "dataRepository/wrapperHelpers.hpp"

#include <ctime>
#include <cstdlib>
#include <unistd.h>

using namespace geos;

real64 drand( real64 min = 0.0, real64 max = 1.0 )
{
  real64 f = (real64)rand() / RAND_MAX;
  return min + f * (max - min);
}

TEST( testPacking, testPacking )
{
  std::srand( std::time( nullptr ));
  constexpr localIndex size = 10000;
  array1d< R1Tensor > veloc( size );
  array1d< R1Tensor > unpacked( size );

  for( localIndex ii = 0; ii < size; ++ii )
    for( localIndex jj = 0; jj < 3; ++jj )
      veloc[ii][jj] = drand();

  buffer_unit_type * null_buf = NULL;
  localIndex calc_size = bufferOps::Pack< false >( null_buf, veloc );
  buffer_type buf( calc_size );
  buffer_unit_type * buffer = &buf[0];
  bufferOps::Pack< true >( buffer, veloc );
  buffer_unit_type const * cbuffer = &buf[0];
  bufferOps::Unpack( cbuffer, unpacked );
  for( localIndex ii = 0; ii < size; ++ii )
    EXPECT_TRUE( veloc[ii] == unpacked[ii] );
}

TEST( testPacking, testPackByIndex )
{

  std::srand( std::time( nullptr ));
  constexpr localIndex size = 10000;
  localIndex pack_count = std::rand() % size;
  array1d< R1Tensor > veloc( size );
  array1d< localIndex > indices( pack_count );
  array1d< R1Tensor > unpacked( size );

  for( localIndex ii = 0; ii < size; ++ii )
    for( localIndex jj = 0; jj < 3; ++jj )
      veloc[ii][jj] = drand();

  for( localIndex ii = 0; ii < pack_count; ++ii )
    indices[ii] = std::rand() % size;

  // std::sort(indices.begin(),indices.end())

  buffer_unit_type * null_buf = NULL;
  // [ ndim, din_i_stride, .. , num_tens_in_arr_i, tens_len_i, data_i, ... ]
  localIndex calc_size = bufferOps::PackByIndex< false >( null_buf, veloc, indices );
  buffer_type buf( calc_size );
  buffer_unit_type * buffer = &buf[0];
  bufferOps::PackByIndex< true >( buffer, veloc, indices );
  buffer_unit_type const * cbuffer = &buf[0];
  bufferOps::UnpackByIndex( cbuffer, unpacked, indices );
  for( localIndex ii = 0; ii < size; ++ii )
  {
    if( std::find( indices.begin(), indices.end(), ii ) != indices.end() )
    {
      EXPECT_EQ( veloc[ii], unpacked[ii] );
    }
  }
}

TEST( testPacking, testTensorPacking )
{
  std::srand( std::time( nullptr ));
  array1d< R1Tensor > tns( 1 );
  for( localIndex ii = 0; ii < 3; ++ii )
    tns[0][ii] = drand();

  buffer_unit_type * null_buf = nullptr;
  parallelDeviceEvents packEvents;
  localIndex calc_size = bufferOps::PackDevice< false >( null_buf, tns.toViewConst(), packEvents );
  buffer_type buf( calc_size );
  buffer_unit_type * b = &buf[0];
  bufferOps::PackDevice< true >( b, tns.toViewConst(), packEvents );
  waitAllDeviceEvents( packEvents );

  array1d< R1Tensor > unp( 1 );
  buffer_unit_type const * bc = &buf[0];
  parallelDeviceEvents unpackEvents;
  bufferOps::UnpackDevice( bc, unp.toView(), unpackEvents );
  waitAllDeviceEvents( unpackEvents );
  unp.move( hostMemorySpace );
  for( localIndex ii = 0; ii < 3; ++ii )
    EXPECT_TRUE( tns[0][ii] = unp[0][ii] );
}

// void printArray( arrayView1d< R1Tensor const > const & arr,
//                  arrayView1d< R1Tensor const > const & unpackArray )
// {
//   printf( "arr.size() = %d\n", arr.size() );
//   printf( "unpackArray.size() = %d\n", unpackArray.size() );

//   for( localIndex ii = 0; ii < unpackArray.size(); ++ii )
//   {
//     if( !( arr[ii] == unpackArray[ii] ) )
//     {
//       printf( "arr[%d]         = ( %f, %f, %f )\n", ii, arr[ii][0], arr[ii][1], arr[ii][2] );
//       printf( "unPackarray[%d] = ( %f, %f, %f ) : ", ii, unpackArray[ii][0], unpackArray[ii][1], unpackArray[ii][2] );

//       forAll< geos::parallelDevicePolicy<1> >( 1, [=] GEOS_DEVICE ( localIndex )
//       {
//         printf( "( %f, %f, %f ) : ", unpackArray[ii][0], unpackArray[ii][1], unpackArray[ii][2] );
//       } );

//       forAll< serialPolicy >( 1, [=]( localIndex )
//       {
//         printf( "( %f, %f, %f )\n", unpackArray[ii][0], unpackArray[ii][1], unpackArray[ii][2] );
//       } );
//       }
//   }
// }

TEST( testPacking, testPackingDevice )
{
  std::srand( std::time( nullptr ));
  constexpr localIndex size = 10000;
  array1d< R1Tensor > veloc( size );
  array1d< R1Tensor > unpacked( size );

  for( localIndex ii = 0; ii < size; ++ii )
    for( localIndex jj = 0; jj < 3; ++jj )
      veloc[ii][jj] = drand();

  buffer_unit_type * null_buf = NULL;
  parallelDeviceEvents packEvents;
  localIndex calc_size = bufferOps::PackDevice< false >( null_buf, veloc.toViewConst(), packEvents );

  buffer_type buf( calc_size );
  buffer_unit_type * buffer = &buf[0];
  bufferOps::PackDevice< true >( buffer, veloc.toViewConst(), packEvents );
  // waitAllDeviceEvents( packEvents );

  // R1Tensor const * const castedBuffer = reinterpret_cast< R1Tensor const * >( &buf[16] );
  // for( localIndex ii = 0; ii < size; ++ii )
  // {
  //   printf( " %d = ( %f, %f, %f ) != ( %f, %f, %f )\n", ii, veloc[ii][0], veloc[ii][1], veloc[ii][2], castedBuffer[ii][0],
  // castedBuffer[ii][1], castedBuffer[ii][2] );
  // }

  buffer_unit_type const * cbuffer = &buf[0];
  parallelDeviceEvents unpackEvents;
  bufferOps::UnpackDevice( cbuffer, unpacked.toView(), unpackEvents );
  waitAllDeviceEvents( unpackEvents );
  unpacked.move( hostMemorySpace );

//  printArray( veloc, unpacked.toViewConst() );

  for( localIndex ii = 0; ii < size; ++ii )
  {
    EXPECT_EQ( veloc[ii], unpacked[ii] );
//    if( !(veloc[ii] == unpacked[ii]) )
//    {
//      printf( " veloc[%d]    = ( %f, %f, %f )\n", ii, veloc[ii][0], veloc[ii][1], veloc[ii][2] );
//      printf( " unpacked[%d] = ( %f, %f, %f )\n", ii, unpacked[ii][0], unpacked[ii][1], unpacked[ii][2] );
//    }
  }
}

TEST( testPacking, testPackingDeviceHelper )
{
  std::srand( std::time( nullptr ));
  constexpr localIndex size = 10000;
  array1d< double > veloc( size );
  array1d< double > unpacked( size );

  for( localIndex ii = 0; ii < size; ++ii )
    //for( localIndex jj = 0; jj < 3; ++jj )
    veloc[ii] = drand();

  buffer_unit_type * null_buf = NULL;
  parallelDeviceEvents packEvents;
  localIndex calc_size = bufferOps::PackDevice< false >( null_buf, veloc.toViewConst(), packEvents );

  buffer_type buf( calc_size );
  buffer_unit_type * buffer = &buf[0];
  dataRepository::wrapperHelpers::PackDevice< true >( buffer, veloc.toViewConst(), packEvents );
  waitAllDeviceEvents( packEvents );

  parallelDeviceEvents unpackEvents;
  buffer_unit_type const * cbuffer = &buf[0];
  dataRepository::wrapperHelpers::UnpackDevice( cbuffer, unpacked.toView(), unpackEvents );
  waitAllDeviceEvents( unpackEvents );
  unpacked.move( hostMemorySpace );
  for( localIndex ii = 0; ii < size; ++ii )
    EXPECT_EQ( veloc[ii], unpacked[ii] );
}

TEST( testPacking, testPackByIndexDevice )
{
  std::srand( std::time( nullptr ));
  constexpr localIndex size = 10000;
  localIndex pack_count = 5000;
  array1d< R1Tensor > veloc( size );
  array1d< localIndex > indices( pack_count );
  array1d< R1Tensor > unpacked( size );

  for( localIndex ii = 0; ii < size; ++ii )
    for( localIndex jj = 0; jj < 3; ++jj )
      veloc[ii][jj] = drand();

  for( localIndex ii = 0; ii < pack_count; ++ii )
    indices[ii] = std::rand() % size;

  std::sort( indices.begin(), indices.end());

  buffer_unit_type * null_buf = NULL;
  // [ num_dim, stride_i.. , tensor_0, tensor_1, ..., tensor_n ]
  arrayView1d< R1Tensor const > const & veloc_view = veloc.toViewConst();
  parallelDeviceEvents packEvents;
  localIndex calc_size = bufferOps::PackByIndexDevice< false >( null_buf, veloc_view, indices.toViewConst(), packEvents );
  buffer_type buf( calc_size );
  buffer_unit_type * buffer = &buf[0];
  localIndex packed_size = bufferOps::PackByIndexDevice< true >( buffer, veloc_view, indices.toViewConst(), packEvents );
  EXPECT_EQ ( calc_size, packed_size );
  waitAllDeviceEvents( packEvents );
  buffer_unit_type const * cbuffer = &buf[0];
  parallelDeviceEvents unpackEvents;
  localIndex unpacked_size = bufferOps::UnpackByIndexDevice( cbuffer, unpacked.toView(), indices.toViewConst(), unpackEvents );
  EXPECT_EQ ( unpacked_size, packed_size );
  waitAllDeviceEvents( unpackEvents );
  unpacked.move( hostMemorySpace );
  for( localIndex ii = 0; ii < size; ++ii )
  {
    if( std::find( indices.begin(), indices.end(), ii ) != indices.end() )
    {
      EXPECT_EQ( veloc[ii], unpacked[ii] );
    }
  }
}

#define CHECK_EXPLICIT_DEVICE_BUFFER_OPS( TYPE, NDIM, PERMUTATION ) \
  { \
    using ArrayType = Array< TYPE, NDIM, PERMUTATION >; \
    localIndex dimensions[NDIM]; \
    for( localIndex dimension = 0; dimension < NDIM; ++dimension ) \
    { \
      dimensions[dimension] = 2; \
    } \
    ArrayType source; \
    source.resize( NDIM, dimensions ); \
    for( localIndex valueIndex = 0; valueIndex < source.size(); ++valueIndex ) \
    { \
      source.data()[valueIndex] = static_cast< double >( valueIndex + 1 ); \
    } \
    \
    buffer_unit_type * nullBuffer = nullptr; \
    parallelDeviceEvents sizeEvents; \
    localIndex const fullCalculatedSize = \
      bufferOps::PackDevice< false >( nullBuffer, source.toViewConst(), sizeEvents ); \
    EXPECT_EQ( nullBuffer, nullptr ); \
    EXPECT_GT( fullCalculatedSize, 0 ); \
    \
    buffer_type fullBuffer( fullCalculatedSize ); \
    buffer_unit_type * fullWriteBuffer = fullBuffer.data(); \
    parallelDeviceEvents fullPackEvents; \
    localIndex const fullPackedSize = \
      bufferOps::PackDevice< true >( fullWriteBuffer, source.toViewConst(), fullPackEvents ); \
    EXPECT_EQ( fullPackedSize, fullCalculatedSize ); \
    EXPECT_EQ( fullWriteBuffer, fullBuffer.data() + fullPackedSize ); \
    waitAllDeviceEvents( fullPackEvents ); \
    \
    ArrayType fullDestination; \
    fullDestination.resize( NDIM, dimensions ); \
    buffer_unit_type const * fullReadBuffer = fullBuffer.data(); \
    parallelDeviceEvents fullUnpackEvents; \
    localIndex const fullUnpackedSize = \
      bufferOps::UnpackDevice( fullReadBuffer, fullDestination.toView(), fullUnpackEvents ); \
    EXPECT_EQ( fullUnpackedSize, fullPackedSize ); \
    EXPECT_EQ( fullReadBuffer, fullBuffer.data() + fullUnpackedSize ); \
    waitAllDeviceEvents( fullUnpackEvents ); \
    source.move( hostMemorySpace ); \
    fullDestination.move( hostMemorySpace ); \
    for( localIndex valueIndex = 0; valueIndex < source.size(); ++valueIndex ) \
    { \
      EXPECT_EQ( fullDestination.data()[valueIndex], source.data()[valueIndex] ); \
    } \
    \
    array1d< localIndex > indices( dimensions[0] ); \
    for( localIndex index = 0; index < indices.size(); ++index ) \
    { \
      indices[index] = index; \
    } \
    nullBuffer = nullptr; \
    parallelDeviceEvents indexedSizeEvents; \
    localIndex const indexedCalculatedSize = bufferOps::PackByIndexDevice< false >( \
      nullBuffer, source.toViewConst(), indices.toViewConst(), indexedSizeEvents ); \
    EXPECT_EQ( nullBuffer, nullptr ); \
    EXPECT_GT( indexedCalculatedSize, 0 ); \
    \
    buffer_type indexedBuffer( indexedCalculatedSize ); \
    buffer_unit_type * indexedWriteBuffer = indexedBuffer.data(); \
    parallelDeviceEvents indexedPackEvents; \
    localIndex const indexedPackedSize = bufferOps::PackByIndexDevice< true >( \
      indexedWriteBuffer, source.toViewConst(), indices.toViewConst(), indexedPackEvents ); \
    EXPECT_EQ( indexedPackedSize, indexedCalculatedSize ); \
    localIndex const indexedWrittenSize = indexedWriteBuffer - indexedBuffer.data(); \
    EXPECT_LE( indexedWrittenSize, indexedPackedSize ); \
    EXPECT_LT( indexedPackedSize - indexedWrittenSize, static_cast< localIndex >( sizeof( TYPE ) ) ); \
    waitAllDeviceEvents( indexedPackEvents ); \
    \
    ArrayType replaceDestination; \
    replaceDestination.resize( NDIM, dimensions ); \
    for( localIndex valueIndex = 0; valueIndex < replaceDestination.size(); ++valueIndex ) \
    { \
      replaceDestination.data()[valueIndex] = -1.0; \
    } \
    buffer_unit_type const * replaceReadBuffer = indexedBuffer.data(); \
    parallelDeviceEvents replaceEvents; \
    localIndex const replaceSize = bufferOps::UnpackByIndexDevice( \
      replaceReadBuffer, replaceDestination.toView(), indices.toViewConst(), replaceEvents, MPI_REPLACE ); \
    EXPECT_EQ( replaceSize, indexedPackedSize ); \
    localIndex const replaceReadSize = replaceReadBuffer - indexedBuffer.data(); \
    EXPECT_LE( replaceReadSize, replaceSize ); \
    EXPECT_LT( replaceSize - replaceReadSize, static_cast< localIndex >( sizeof( TYPE ) ) ); \
    waitAllDeviceEvents( replaceEvents ); \
    replaceDestination.move( hostMemorySpace ); \
    source.move( hostMemorySpace ); \
    for( localIndex valueIndex = 0; valueIndex < source.size(); ++valueIndex ) \
    { \
      EXPECT_EQ( replaceDestination.data()[valueIndex], source.data()[valueIndex] ); \
    } \
    \
    ArrayType sumDestination; \
    sumDestination.resize( NDIM, dimensions ); \
    for( localIndex valueIndex = 0; valueIndex < sumDestination.size(); ++valueIndex ) \
    { \
      sumDestination.data()[valueIndex] = 2.0; \
    } \
    buffer_unit_type const * sumReadBuffer = indexedBuffer.data(); \
    parallelDeviceEvents sumEvents; \
    localIndex const sumSize = bufferOps::UnpackByIndexDevice( \
      sumReadBuffer, sumDestination.toView(), indices.toViewConst(), sumEvents, MPI_SUM ); \
    EXPECT_EQ( sumSize, indexedPackedSize ); \
    localIndex const sumReadSize = sumReadBuffer - indexedBuffer.data(); \
    EXPECT_LE( sumReadSize, sumSize ); \
    EXPECT_LT( sumSize - sumReadSize, static_cast< localIndex >( sizeof( TYPE ) ) ); \
    waitAllDeviceEvents( sumEvents ); \
    sumDestination.move( hostMemorySpace ); \
    source.move( hostMemorySpace ); \
    for( localIndex valueIndex = 0; valueIndex < source.size(); ++valueIndex ) \
    { \
      TYPE expected = source.data()[valueIndex]; \
      expected += 2.0; \
      EXPECT_EQ( sumDestination.data()[valueIndex], expected ); \
    } \
    \
    ArrayType maxIncomingDestination; \
    maxIncomingDestination.resize( NDIM, dimensions ); \
    for( localIndex valueIndex = 0; valueIndex < maxIncomingDestination.size(); ++valueIndex ) \
    { \
      maxIncomingDestination.data()[valueIndex] = 0.0; \
    } \
    buffer_unit_type const * maxIncomingReadBuffer = indexedBuffer.data(); \
    parallelDeviceEvents maxIncomingEvents; \
    localIndex const maxIncomingSize = bufferOps::UnpackByIndexDevice( \
      maxIncomingReadBuffer, maxIncomingDestination.toView(), indices.toViewConst(), maxIncomingEvents, MPI_MAX ); \
    EXPECT_EQ( maxIncomingSize, indexedPackedSize ); \
    localIndex const maxIncomingReadSize = maxIncomingReadBuffer - indexedBuffer.data(); \
    EXPECT_LE( maxIncomingReadSize, maxIncomingSize ); \
    EXPECT_LT( maxIncomingSize - maxIncomingReadSize, static_cast< localIndex >( sizeof( TYPE ) ) ); \
    waitAllDeviceEvents( maxIncomingEvents ); \
    maxIncomingDestination.move( hostMemorySpace ); \
    source.move( hostMemorySpace ); \
    for( localIndex valueIndex = 0; valueIndex < source.size(); ++valueIndex ) \
    { \
      EXPECT_EQ( maxIncomingDestination.data()[valueIndex], source.data()[valueIndex] ); \
    } \
    \
    ArrayType maxExistingDestination; \
    maxExistingDestination.resize( NDIM, dimensions ); \
    for( localIndex valueIndex = 0; valueIndex < maxExistingDestination.size(); ++valueIndex ) \
    { \
      maxExistingDestination.data()[valueIndex] = 1000.0; \
    } \
    buffer_unit_type const * maxExistingReadBuffer = indexedBuffer.data(); \
    parallelDeviceEvents maxExistingEvents; \
    localIndex const maxExistingSize = bufferOps::UnpackByIndexDevice( \
      maxExistingReadBuffer, maxExistingDestination.toView(), indices.toViewConst(), maxExistingEvents, MPI_MAX ); \
    EXPECT_EQ( maxExistingSize, indexedPackedSize ); \
    localIndex const maxExistingReadSize = maxExistingReadBuffer - indexedBuffer.data(); \
    EXPECT_LE( maxExistingReadSize, maxExistingSize ); \
    EXPECT_LT( maxExistingSize - maxExistingReadSize, static_cast< localIndex >( sizeof( TYPE ) ) ); \
    waitAllDeviceEvents( maxExistingEvents ); \
    maxExistingDestination.move( hostMemorySpace ); \
    for( localIndex valueIndex = 0; valueIndex < maxExistingDestination.size(); ++valueIndex ) \
    { \
      TYPE expected; \
      expected = 1000.0; \
      EXPECT_EQ( maxExistingDestination.data()[valueIndex], expected ); \
    } \
  }

#define CHECK_EXPLICIT_DEVICE_BUFFER_OPS_UP_TO_2D( TYPE ) \
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS( TYPE, 1, RAJA::PERM_I ); \
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS( TYPE, 2, RAJA::PERM_JI ); \
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS( TYPE, 2, RAJA::PERM_IJ )

#define CHECK_EXPLICIT_DEVICE_BUFFER_OPS_UP_TO_3D( TYPE ) \
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS_UP_TO_2D( TYPE ); \
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS( TYPE, 3, RAJA::PERM_JKI ); \
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS( TYPE, 3, RAJA::PERM_IKJ ); \
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS( TYPE, 3, RAJA::PERM_IJK )

#define CHECK_EXPLICIT_DEVICE_BUFFER_OPS_UP_TO_4D( TYPE ) \
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS_UP_TO_3D( TYPE ); \
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS( TYPE, 4, RAJA::PERM_JKLI ); \
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS( TYPE, 4, RAJA::PERM_IKLJ ); \
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS( TYPE, 4, RAJA::PERM_IJLK ); \
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS( TYPE, 4, RAJA::PERM_IJKL )

#define CHECK_EXPLICIT_DEVICE_BUFFER_OPS_UP_TO_5D( TYPE ) \
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS_UP_TO_4D( TYPE ); \
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS( TYPE, 5, RAJA::PERM_JKLMI ); \
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS( TYPE, 5, RAJA::PERM_IKLMJ ); \
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS( TYPE, 5, RAJA::PERM_IJLMK ); \
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS( TYPE, 5, RAJA::PERM_IJKML ); \
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS( TYPE, 5, RAJA::PERM_IJKLM )

TEST( testPacking, testExplicitDeviceInstantiationMatrix )
{
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS_UP_TO_3D( int );
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS_UP_TO_3D( long int );
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS_UP_TO_3D( long long int );
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS_UP_TO_5D( real32 );
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS_UP_TO_5D( real64 );
  CHECK_EXPLICIT_DEVICE_BUFFER_OPS_UP_TO_3D( R1Tensor );
}

#undef CHECK_EXPLICIT_DEVICE_BUFFER_OPS_UP_TO_5D
#undef CHECK_EXPLICIT_DEVICE_BUFFER_OPS_UP_TO_4D
#undef CHECK_EXPLICIT_DEVICE_BUFFER_OPS_UP_TO_3D
#undef CHECK_EXPLICIT_DEVICE_BUFFER_OPS_UP_TO_2D
#undef CHECK_EXPLICIT_DEVICE_BUFFER_OPS


int main( int ac, char * av[] )
{
  ::testing::InitGoogleTest( &ac, av );
  int const result = RUN_ALL_TESTS();
  return result;
}
