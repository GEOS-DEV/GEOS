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

#include <gtest/gtest.h>

#include "LvArray/src/Array.hpp"
#include "dataRepository/BufferOps.hpp"
#include "dataRepository/BufferOpsDevice.hpp"
#include "mainInterface/initialization.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

#include "dataRepository/wrapperHelpers.hpp"

#include <ctime>
#include <cstdlib>

using namespace geosx;

#ifndef GTEST_SKIP
#define GTEST_SKIP() return
#endif

#define SKIP_TEST_IF( COND, REASON ) \
  do \
  { \
    if( COND ) \
    { \
      GEOSX_WARNING( "This test is currently known to fail when " #COND " because:\n" REASON "\n" \
                                                                                             "Therefore, we skip it entirely for this run (may show as PASSED or SKIPPED)" ); \
      GTEST_SKIP(); \
    } \
  } while( 0 )

#define SKIP_TEST_IN_SERIAL( REASON ) \
  do \
  { \
    int const mpiSize = MpiWrapper::commSize( MPI_COMM_GEOSX ); \
    SKIP_TEST_IF( mpiSize == 1, REASON ); \
  } while( 0 )

#define SKIP_TEST_IN_PARALLEL( REASON ) \
  do \
  { \
    int const mpiSize = MpiWrapper::commSize( MPI_COMM_GEOSX ); \
    SKIP_TEST_IF( mpiSize != 1, REASON ); \
  } while( 0 )

real64 drand( real64 min = 0.0, real64 max = 1.0 )
{
  real64 f = (real64)rand() / RAND_MAX;
  return min + f * (max - min);
}

// all the types we instantiate packing functions for, execpting r1tensor which is being removed
using PackingTemplateTestTypes = ::testing::Types< int, long int, long long int, real32, real64 >;
using PackingTemplateDims = ::testing::Types< 1, 2, 3 >;

template< typename TYPE, int DIM >
class PackingTemplateTest : public ::testing::Test
{
  public:
    void pack( )
    {
      std::srand( std::time( nullptr ));
      constexpr localIndex size = 10000;
      array< TYPE, DIM > veloc( size, 3 );
      array< TYPE, DIM > unpacked( size, 3 );

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

        for( localIndex jj = 0; jj < 3; ++jj )
          EXPECT_EQ( veloc[ii][jj], unpacked[ii][jj] );
    }
    void packByIndex( )
    {

    }
    void packDevice( )
    {

    }
    void packByIndexDevice( )
    {

    }
}

TYPED_TEST_SUITE( PackingTemplateTest, PackingTemplateTestTypes, );

TYPED_TEST( testPacking, testPacking )
{

}

TEST( testPacking, testPackByIndex )
{

  std::srand( std::time( nullptr ));
  constexpr localIndex size = 10000;
  localIndex pack_count = std::rand() % size;
  array< TYPE, DIM > veloc( size, 3 );
  array1d< localIndex > indices( pack_count );
  array< TYPE, DIM > unpacked( size, 3 );

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
      for( localIndex jj = 0; jj < 3; ++jj )
      {
        EXPECT_EQ( veloc[ii][jj], unpacked[ii][jj] );
      }
    }
  }
}

TEST( testPacking, testPackingDevice )
{
  std::srand( std::time( nullptr ));
  constexpr localIndex size = 10000;
  array< TYPE, DIM > veloc( size, 3 );
  array< TYPE, DIM > unpacked( size, 3 );

  for( localIndex ii = 0; ii < size; ++ii )
    for( localIndex jj = 0; jj < 3; ++jj )
      veloc[ii][jj] = drand();

  buffer_unit_type * null_buf = NULL;
  parallelDeviceEvents packEvents;
  localIndex calc_size = bufferOps::PackDevice< false >( null_buf, veloc.toViewConst(), packEvents );

  buffer_type buf( calc_size );
  buffer_unit_type * buffer = &buf[0];
  bufferOps::PackDevice< true >( buffer, veloc.toViewConst(), packEvents );
  waitAllDeviceEvents( packEvents );

  buffer_unit_type const * cbuffer = &buf[0];
  parallelDeviceEvents unpackEvents;
  bufferOps::UnpackDevice( cbuffer, unpacked.toView(), unpackEvents );
  waitAllDeviceEvents( unpackEvents );
  unpacked.move( LvArray::MemorySpace::host );
  for( localIndex ii = 0; ii < size; ++ii )
    for( localIndex jj = 0; jj < 3; ++jj )
      EXPECT_EQ( veloc[ii][jj], unpacked[ii][jj] );
}

TEST( testPacking, testPackingDeviceHelper )
{
  std::srand( std::time( nullptr ));
  constexpr localIndex size = 10000;
  array1d< double > veloc( size );
  array1d< double > unpacked( size );

  for( localIndex ii = 0; ii < size; ++ii )
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
  unpacked.move( LvArray::MemorySpace::host );
  for( localIndex ii = 0; ii < size; ++ii )
    EXPECT_EQ( veloc[ii], unpacked[ii] );
}

TEST( testPacking, testPackByIndexDevice )
{
  std::srand( std::time( nullptr ));
  constexpr localIndex size = 10000;
  localIndex pack_count = 5000;
  array< TYPE, DIM > veloc( size, 3 );
  array1d< localIndex > indices( pack_count );
  array< TYPE, DIM > unpacked( size, 3 );

  for( localIndex ii = 0; ii < size; ++ii )
    for( localIndex jj = 0; jj < 3; ++jj )
      veloc[ii][jj] = drand();

  for( localIndex ii = 0; ii < pack_count; ++ii )
    indices[ii] = std::rand() % size;

  std::sort( indices.begin(), indices.end());

  buffer_unit_type * null_buf = NULL;
  // [ num_dim, stride_i.. , tensor_0, tensor_1, ..., tensor_n ]
  arrayView2d< const double > const & veloc_view = veloc.toViewConst();
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
  unpacked.move( LvArray::MemorySpace::host );
  for( localIndex ii = 0; ii < size; ++ii )
  {
    if( std::find( indices.begin(), indices.end(), ii ) != indices.end() )
    {
      for( localIndex jj = 0; jj < 3; ++jj )
      {
        EXPECT_EQ( veloc[ii][jj], unpacked[ii][jj] );
      }
    }
  }
}


int main( int ac, char * av[] )
{
  ::testing::InitGoogleTest( &ac, av );
  geosx::basicSetup( ac, av );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
