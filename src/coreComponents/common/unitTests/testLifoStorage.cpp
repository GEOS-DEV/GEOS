#include "mainInterface/initialization.hpp"
#include "common/lifoStorage.hpp"
#include "LvArray/src/Array.hpp"
#include "LvArray/src/MallocBuffer.hpp"
#if defined(LVARRAY_USE_CHAI)
  #include "LvArray/src/ChaiBuffer.hpp"
#endif

#include <gtest/gtest.h>

using namespace geosx;


TEST( LifoStorageTest, LifoStorage )
{
  int elemCnt = 10;
  int numberOfElementsOnDevice = 2;
  int numberOfElementsOnHost = 3;
  int totalNumberOfBuffers = 10;

  array1d< float > array( elemCnt );
  lifoStorage< float > lifo( "lifo", array, numberOfElementsOnDevice, numberOfElementsOnHost, totalNumberOfBuffers );

  for( int j = 0; j < totalNumberOfBuffers; j++ )
  {
    for( int i = 0; i < elemCnt; i++ )
      array[i] = j*elemCnt+i;
    lifo.push( array );
  }

  for( int j = 0; j < totalNumberOfBuffers; j++ )
  {
    lifo.pop( array );
    array.move( LvArray::MemorySpace::host, false );
    for( int i = 0; i < elemCnt; i++ )
      EXPECT_EQ( array[i], (totalNumberOfBuffers-j-1)*elemCnt+i );
  }


}

int main( int ac, char * av[] )
{
  ::testing::InitGoogleTest( &ac, av );
  int const result = RUN_ALL_TESTS();
  return result;
}
