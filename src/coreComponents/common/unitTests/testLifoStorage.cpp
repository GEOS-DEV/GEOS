#include "mainInterface/initialization.hpp"
#include "common/lifoStorage.hpp"
#include "LvArray/src/Array.hpp"
#include "LvArray/src/MallocBuffer.hpp"
#if defined(LVARRAY_USE_CHAI)
  #include "LvArray/src/ChaiBuffer.hpp"
#endif

#include <gtest/gtest.h>

#ifndef __CUDA_ARCH__
#define PORTABLE_EXPECT_EQ( L, R ) EXPECT_EQ( L, R )
#define PORTABLE_EXPECT_NEAR( L, R, EPSILON ) EXPECT_LE( math::abs( ( L ) -( R ) ), EPSILON ) << \
    STRINGIZE( L ) " = " << ( L ) << "\n" << STRINGIZE( R ) " = " << ( R );
#else
#define PORTABLE_EXPECT_EQ( L, R ) LVARRAY_ERROR_IF_NE( L, R )
#define PORTABLE_EXPECT_NEAR( L, R, EPSILON ) LVARRAY_ERROR_IF_GE_MSG( math::abs( ( L ) -( R ) ), EPSILON, \
                                                                       STRINGIZE( L ) " = " << ( L ) << "\n" << STRINGIZE( R ) " = " << ( R ) );
#endif


template< typename >
struct RAJAHelper
{};

using serialPolicy = RAJA::loop_exec;

template<>
struct RAJAHelper< serialPolicy >
{
  using ReducePolicy = RAJA::seq_reduce;
  using AtomicPolicy = RAJA::seq_atomic;
  static constexpr LvArray::MemorySpace space = LvArray::MemorySpace::host;
};

#if defined(RAJA_ENABLE_OPENMP)

using parallelHostPolicy = RAJA::omp_parallel_for_exec;

template<>
struct RAJAHelper< parallelHostPolicy >
{
  using ReducePolicy = RAJA::omp_reduce;
  using AtomicPolicy = RAJA::omp_atomic;
  static constexpr LvArray::MemorySpace space = LvArray::MemorySpace::host;
};

#endif

#if defined(LVARRAY_USE_CUDA)

template< unsigned long THREADS_PER_BLOCK >
using devicePolicy = RAJA::cuda_exec< THREADS_PER_BLOCK >;

template< unsigned long N >
struct RAJAHelper< RAJA::cuda_exec< N > >
{
  using ReducePolicy = RAJA::cuda_reduce;
  using AtomicPolicy = RAJA::cuda_atomic;
  static constexpr LvArray::MemorySpace space = LvArray::MemorySpace::cuda;
};

#endif

using namespace geosx;


template< typename POLICY >
void testLifoStorage( )
{
  int elemCnt = 10;
  int numberOfElementsOnDevice = 2;
  int numberOfElementsOnHost = 3;
  int totalNumberOfBuffers = 10;

  array1d< float > array( elemCnt );
  array.move( RAJAHelper< POLICY >::space );
  lifoStorage< float > lifo( "lifo", array, numberOfElementsOnDevice, numberOfElementsOnHost, totalNumberOfBuffers );

  for( int j = 0; j < totalNumberOfBuffers; j++ )
  {

    float * dataPointer = array.data();
    forAll< POLICY >( elemCnt, [dataPointer, j, elemCnt] GEOSX_HOST_DEVICE ( int i ) { dataPointer[ i ] = j*elemCnt+i; } );
    lifo.push( array );
  }

  for( int j = 0; j < totalNumberOfBuffers; j++ )
  {
    lifo.pop( array );
    float * dataPointer = array.data();
    forAll< POLICY >( elemCnt, [dataPointer, totalNumberOfBuffers, j, elemCnt] GEOSX_HOST_DEVICE ( int i )
    {
      GEOSX_ERROR_IF_NE( dataPointer[ i ], (totalNumberOfBuffers-j-1)*elemCnt+i );
    } );
  }
}

template< typename POLICY >
void testLifoStorageAsync( )
{
  int elemCnt = 10;
  int numberOfElementsOnDevice = 2;
  int numberOfElementsOnHost = 3;
  int totalNumberOfBuffers = 10;

  array1d< float > array( elemCnt );
  array.move( RAJAHelper< POLICY >::space );
  lifoStorage< float > lifo( "lifo", array, numberOfElementsOnDevice, numberOfElementsOnHost, totalNumberOfBuffers );

  for( int j = 0; j < totalNumberOfBuffers; j++ )
  {

    float * dataPointer = array.data();
    lifo.pushWait( );
    forAll< POLICY >( elemCnt, [dataPointer, j, elemCnt] GEOSX_HOST_DEVICE ( int i ) { dataPointer[ i ] = j*elemCnt+i; } );
    lifo.pushAsync( array );
  }

  for( int j = 0; j < totalNumberOfBuffers; j++ )
  {
    lifo.popAsync( array );
    lifo.popWait( );
    float * dataPointer = array.data();
    forAll< POLICY >( elemCnt, [dataPointer, totalNumberOfBuffers, j, elemCnt] GEOSX_HOST_DEVICE ( int i )
    {
      GEOSX_ERROR_IF_NE( dataPointer[ i ], (totalNumberOfBuffers-j-1)*elemCnt+i );
    } );
  }
}


TEST( LifoStorageTest, LifoStorageHost )
{
  testLifoStorage< serialPolicy >( );
}


TEST( LifoStorageTest, LifoStorageCUDA )
{
  testLifoStorage< parallelDevicePolicy< > >( );
}

TEST( LifoStorageTest, LifoStorageAsyncHost )
{
  testLifoStorageAsync< serialPolicy >( );
}


TEST( LifoStorageTest, LifoStorageAsyncCUDA )
{
  testLifoStorageAsync< parallelDevicePolicy< > >( );
}

int main( int ac, char * av[] )
{
  ::testing::InitGoogleTest( &ac, av );
  int const result = RUN_ALL_TESTS();
  return result;
}
