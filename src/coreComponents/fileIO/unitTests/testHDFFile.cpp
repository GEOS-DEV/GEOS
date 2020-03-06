#include "fileIO/hdf/HDFFile.hpp"
#include "dataRepository/BufferOps.hpp"
#include "dataRepository/BufferOpsDevice.hpp"
#include <gtest/gtest.h>


// todo: template
#define DIM 1
#define ARR_TYPE real64
#define IND_TYPE localIndex

using namespace geosx;

TEST( testHDFTraits, can_hdf_io )
{
  static_assert( can_hdf_io< real32 >, "Should be true.");
  static_assert( can_hdf_io< real64 >, "Should be true.");
  static_assert( can_hdf_io< integer >, "Should be true.");
  static_assert( can_hdf_io< localIndex >, "Should be true.");
  static_assert( can_hdf_io< globalIndex >, "Should be true.");

  static_assert( can_hdf_io< const real32 >, "Should be true.");
  static_assert( can_hdf_io< const real64 >, "Should be true.");
  static_assert( can_hdf_io< const integer >, "Should be true.");
  static_assert( can_hdf_io< const localIndex >, "Should be true.");
  static_assert( can_hdf_io< const globalIndex >, "Should be true.");
}

// need versions that use raw arrays as well

TEST( testHDFIO, WholeTabularIO )
{
  srand(time(NULL));
  localIndex dims[DIM] = {0};
  for(integer dd = 0; dd < DIM; ++dd )
  {
    dims[dd] = rand();
  }
  Array<ARR_TYPE, DIM> arr(rand());

  {
    HDFFile file("whole_array");

    // specify the tables (done during init)
    HDFTable time_table("Array1 Time", "arr1_t");
    time_table.AddCols(1,1,sizeof(real64),typeid(real64),[](localIndex idx){ return std::string("nd_") + std::to_string(idx); });
    time_table.Finalize( );
    time_table.CreateInTarget( file );

    HDFTable table("Array1","arr1");
    SpecFromArray(table,arr);
    table.Finalize( );
    table.CreateInTarget( file );

    localIndex indices[] = { 0, 5, 9 };
    HDFTable idx_table("Array1 Indexed","arr1_idx");
    SpecFromArrayIndices(idx_table,arr,3,&indices[0]);
    idx_table.Finalize( );
    idx_table.CreateInTarget( file );

    HDFTableIO table_io( table );
    HDFTableIO time_table_io( time_table );

    // collect time history data into buffers
    buffer_unit_type * buf_head = NULL;
    size_t bufferSize = bufferOps::Pack<false>(buf_head,arr.toView());
    std::vector<buffer_unit_type> buf(bufferSize);
    buf_head = &buf[0];
    bufferOps::Pack<true>(buf_head,arr.toView());
    // we pack arr.strides() before the actual array
    size_t pack_meta_size = DIM * sizeof(localIndex);

    real64 time = 0.4;
    std::vector<buffer_unit_type> tbuf(sizeof(decltype(time)));
    buf_head = &tbuf[0];
    bufferOps::Pack<true>(buf_head,time);

    // done on collection
    // buffer io on the table
    table_io.BufferRow( &buf[pack_meta_size] );
    time_table_io.BufferRow( &tbuf[0] );

    // done on write
    // open the table to write
    table_io.Open( file );
    table_io.WriteBuffered( );
    table_io.Close();

    time_table_io.Open( file );
    time_table_io.WriteBuffered( );
    time_table_io.Close( );

    // clear the table in the target from row 0
    table_io.Open( file );
    table_io.ClearAfterWriteHead( );
    table_io.Close();
  }
}

// TEST( testHDFIO, IndexedTabularIO )
// {
//   srand(time(NULL));
//   localIndex dims[DIM] = {0};
//   for(integer dd = 0; dd < DIM; ++dd )
//   {
//     dims[dd] = rand();
//   }
//   Array<ARR_TYPE, DIM> arr(rand());
//   // rand primary dim
//   Array<IND_TYPE, 1> ind_arr(rand() % dims[0]);

//   {
//     HDFFile file("indexed_array");
//     HDFTableIO<decltype(arr)> table_out("Array1","arr1",ind_arr,arr.size( 0 ),0,"nd_");
//   }
// }

// //
// TEST( testHDFIO, PartialTabularIO )
// {
//   srand(time(NULL));
//   localIndex dims[DIM] = {0};
//   for(integer dd = 0; dd < DIM; ++dd )
//   {
//     dims[dd] = rand();
//   }
//   Array<ARR_TYPE, DIM> arr(rand());
//   array1d<IND_TYPE> ind_arr(rand() % dims[0]);
//   array2d<IND_TYPE> cmp_arr(3,DIM-1);

//   {
//     HDFFile file("arr_output");
//     HDFTableIO<decltype(arr)> table_out("Stress","s",0,ind_arr,cmp_arr,"nd");
//     table_out.OpenTable( file );
//     table_out.AppendRow( arr );
//     table_out.AppendRow( arr );

//     table_out.ClearFrom( 1 );
//   }

//   // check that there is only one row and that it has the correct content
//   {

//   }
// }

// //TEST( testHDFIO, IndexedPartialTabularIO )

// TEST( testHDFIO, WholeArrayTimeHistory )
// {
//   srand(time(NULL));
//   localIndex dims[DIM] = {0};
//   for(integer dd = 0; dd < DIM; ++dd )
//   {
//     dims[dd] = rand();
//   }
//   Array<ARR_TYPE, DIM> arr(rand());

//   {
//     HDFFile file("arr_time_hist");
//     HDFTimeHistoryTable<decltype(arr)> time_hist("Scalar","scl",arr.size(),arr.size(0),0,"nd_");
//     time_hist.OpenTable( file );
//     time_hist.AppendRow( 0.3, arr );
//     time_hist.AppendRow( 0.5, arr );
//     time_hist.ClearFromTime( 0.4 );
//   }

// }

//TEST( testHDFIO, IndexedTabularTimeHistory )
//TEST( testHDFIO, PartialTabularTimeHistory )
//TEST( testHDFIO, IndexedPartialTimeHistory )
//