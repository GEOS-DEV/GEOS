#include "managers/Outputs/TimeHistoryOutput.hpp"
#include "fileIO/hdf/HDFFile.hpp"
#include "dataRepository/BufferOpsDevice.hpp"
#include "managers/initialization.hpp"

#include <gtest/gtest.h>

using namespace geosx;

TEST( testHDFIO, SingleValueHistory )
{
  string filename("single_value");
  HistoryMetadata spec("Time History", 1, std::type_index(typeid(real64)));

  real64 time = 0.0;
  HDFHistIO io( filename, spec );
  io.Init( true, false );
  for( localIndex tidx = 0; tidx < 100; ++tidx )
  {
    time += 0.333;
    buffer_unit_type * buffer = io.GetBufferHead( );
    memcpy(buffer,&time,sizeof(real64));
  }

  io.Write( );

  //read and check the data using hdf api
  // remove( filename.c_str() );
}

TEST( testHDFIO, ArrayHistory )
{
  srand(time(NULL));

  {
    string filename("array1d_history");
    Array<real64, 1> arr(4096);
    real64 count = 0.0;
    forValuesInSlice(arr.toSlice(),[&count](real64 & value)
      {
        value = count++;
      });

    HistoryMetadata spec = getHistoryMetadata("Array1d History",arr);
    HDFHistIO io( filename, spec );
    io.Init( true, false );

    buffer_unit_type * buffer = io.GetBufferHead( );
    bufferOps::PackDevice<true>(buffer,arr.toViewConst( ));

    io.Write( );

    //read and check the data using hdf api
    // remove( filename.c_str() );
  }
  {
    string filename("array2d_history");
    Array<real64, 2> arr(1024,4);
    real64 count = 0.0;
    forValuesInSlice(arr.toSlice(),[&count](real64 & value)
      {
        value = count++;
      });

    HistoryMetadata spec = getHistoryMetadata("Array2d History",arr);
    HDFHistIO io(filename,spec);
    io.Init( true, false );

    buffer_unit_type * buffer = io.GetBufferHead( );
    bufferOps::PackDevice<true>(buffer,arr.toViewConst( ));

    io.Write( );

    //read and check the data using hdf api
    // remove( filename.c_str() );
  }
}

TEST( testHDFIO, IdxArrayHistory )
{
  srand(time(NULL));
  {
    string filename("array1d_idx_history");
    Array<localIndex, 1> idx(256);
    Array<real64, 2> arr(1024,4);
    real64 count = 0.0;
    forValuesInSlice(arr.toSlice(),[&count](real64 & value)
      {
        value = count++;
      });
    forValuesInSlice(idx.toSlice(),[](localIndex & value)
      {
        value = rand() % 1024;
      });

    HistoryMetadata spec = getHistoryMetadata("Array1d Idx History",arr,idx.size( ));
    HDFHistIO io( filename, spec );
    io.Init( true, false );

    buffer_unit_type * buffer = io.GetBufferHead( );
    bufferOps::PackByIndexDevice<true>(buffer,arr.toViewConst( ),idx.toViewConst( ) );

    io.Write( );
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
