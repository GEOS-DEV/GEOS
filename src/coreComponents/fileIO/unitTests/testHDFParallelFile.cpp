#include "managers/Outputs/TimeHistoryOutput.hpp"
#include "fileIO/hdf/HDFFile.hpp"
#include "dataRepository/BufferOpsDevice.hpp"
#include "managers/initialization.hpp"

#include <gtest/gtest.h>

using namespace geosx;

TEST( testHDFIO_parallel, SingleValueHistory )
{
  string filename("single_value_parallel");
  HistoryMetadata spec("Time History", 1, std::type_index(typeid(real64)));

  int rank = MpiWrapper::Comm_rank( );

  HDFHistIO io( filename, spec );
  io.Init( true, false );
  real64 val = 0.0;
  for( localIndex tidx = 0; tidx < 100; ++tidx )
  {
    val += 0.5 * (rank+1);
    buffer_unit_type * buffer = io.GetBufferHead( );
    memcpy(buffer,&val,sizeof(real64));
  }

  io.Write( );

  //read and check the data using hdf api
  // remove( filename.c_str() );
}

int main( int ac, char * av[] )
{
  ::testing::InitGoogleTest( &ac, av );
  geosx::basicSetup( ac, av );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
