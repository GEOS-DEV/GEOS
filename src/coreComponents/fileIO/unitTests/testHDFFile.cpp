#include "managers/Outputs/TimeHistoryOutput.hpp"
#include "fileIO/hdf/HDFFile.hpp"
#include "dataRepository/BufferOps.hpp"
#include "dataRepository/BufferOpsDevice.hpp"
#include <gtest/gtest.h>


// todo: template
#define DIM 1
#define ARR_TYPE real64
#define IND_TYPE localIndex

using namespace geosx;

TEST( testHDFIO, SingleValueTable )
{
  DataSpec spec;
  spec.SetTitleID("Time History", "time_hist");
  spec.Append<real64>(1,1,"Time");
  spec.Finalize();

  real64 time = 0.0;
  HDFTableIO io;
  for( localIndex tidx = 0; tidx < 100; ++tidx )
  {
    time += 0.333;
    buffer_unit_type * buffer = io.GetBufferHead( &spec );
    memcpy(buffer,&time,sizeof(real64));
  }

  string filename( "single_value" );
  io.Init( filename, &spec, false );
  io.Write( filename, &spec );

  //read and check the data using hdf api
  remove( filename.c_str() );
}

TEST( testHDFIO, ArrayTableIO )
{
  srand(time(NULL));

  {
    Array<real64, 1> arr(4096);
    real64 count = 0.0;
    forValuesInSlice(arr.toSlice(),[&count](real64 & value)
      {
        value = count++;
      });

    DataSpec spec;
    spec.SetTitleID("Array1d History", "arr1_hist");
    AppendArraySpec(spec,arr,"Array1d Data");
    spec.Finalize();

    HDFTableIO io;

    buffer_unit_type * buffer = io.GetBufferHead( &spec );
    bufferOps::PackDevice<true>(buffer,arr);

    string filename( "array1d_value" );
    io.Init( filename, &spec, false );
    io.Write( filename, &spec );

    //read and check the data using hdf api
    remove( filename.c_str() );
  }
  {
    Array<real64, 2> arr(1024,4);
    real64 count = 0.0;
    forValuesInSlice(arr.toSlice(),[&count](real64 & value)
      {
        value = count++;
      });

    DataSpec spec;
    spec.SetTitleID("Array2d History", "arr2_hist");
    AppendArraySpec(spec,arr,"Array2d Data");
    spec.Finalize();

    HDFTableIO io;

    buffer_unit_type * buffer = io.GetBufferHead( &spec );
    bufferOps::PackDevice<true>(buffer,arr);

    string filename( "array2d_value" );
    io.Init( filename, &spec, false );
    io.Write( filename, &spec );

    //read and check the data using hdf api
    remove( filename.c_str() );
  }
}

TEST( testHDFIO, IdxArrayTable )
{
  srand(time(NULL));
  {
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

    DataSpec spec;
    spec.SetTitleID("Array1d History", "arr1_hist");
    AppendArrayIndicesSpec(spec,arr,"Array1d Data",idx.size());
    spec.Finalize();

    HDFTableIO io;

    buffer_unit_type * buffer = io.GetBufferHead( &spec );
    bufferOps::PackByIndexDevice<true>(buffer,arr,idx);

    string filename( "array1d_value" );
    io.Init( filename, &spec, false );
    io.Write( filename, &spec );
  }
}