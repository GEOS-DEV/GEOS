//#include "fileIO/Outputs/TimeHistoryOutput.hpp"
#include "fileIO/timeHistory/TimeHistHDF.hpp"
//#include "dataRepository/BufferOpsDevice.hpp"
#include "mainInterface/initialization.hpp"

#include <gtest/gtest.h>

using namespace geosx;

TEST( testHDFIO_parallel, SingleValueHistory )
{
  GEOSX_MARK_FUNCTION;
  string filename( "single_value_parallel" );
  HistoryMetadata spec( "Time History", 1, std::type_index( typeid(real64)));

  int rank = MpiWrapper::commRank( );

  HDFHistIO io( filename, spec );
  io.init( true );
  real64 val = 0.0;
  for( localIndex tidx = 0; tidx < 100; ++tidx )
  {
    val += 0.5 * (rank+1);
    buffer_unit_type * buffer = io.getBufferHead( );
    memcpy( buffer, &val, sizeof(real64));
  }

  io.write( );
}

int main( int ac, char * av[] )
{
  ::testing::InitGoogleTest( &ac, av );
  int const result = RUN_ALL_TESTS();
  return result;
}
