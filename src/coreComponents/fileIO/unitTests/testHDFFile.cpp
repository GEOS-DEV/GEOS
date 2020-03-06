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

TEST( testHDFIO, SingleValueTable )
{
  {
    HDFTable spec("Array1","arr1");
    spec.AddCols(1,1,sizeof(real64),typeid(real64));
    spec.Finalize();

    real64 time = 1.0;
    HDFTableIO out( spec );
    out.BufferRow( reinterpret_cast<buffer_unit_type const * >(&time) );

    HDFFile file( "single_value" );
    out.CreateInTarget( file );
    out.WriteBuffered( file );
  }
}

TEST( testHDFIO, WholeArrayTableIO )
{
  srand(time(NULL));
  localIndex dims[DIM] = {0};
  for(integer dd = 0; dd < DIM; ++dd )
  {
    dims[dd] = rand();
  }
  Array<ARR_TYPE, DIM> arr(rand());

  {
    HDFTable spec("Array1","arr1");
    SpecFromArray(spec,arr);
    spec.Finalize();

    HDFTableIO out( spec );
    out.BufferRow( reinterpret_cast<buffer_unit_type const * >(arr.data()) );

    HDFFile file( "whole_array" );
    out.CreateInTarget( file );
    out.WriteBuffered( file );
  }

  {

    // TimeHistoryCollector * arr_collector = new ArrayTimeHistoryCollector<decltype(arr)>( arr );

    // TimeHistoryUpdate update_event("TimeHistoryUpdate",NULL);
    // TimeHistory & arr_hist = update_event.getTimeHistory( );
    // arr_hist.AddHistory("arr1","Array1",arr_collector);

    // update_event.Execute(1.0,0.5,1,4,0.6,NULL);

    // TimeHistoryOutput output_event("whole_array_hist","TimeHistoryOutput",NULL);
    // output_event.InitHistoryFile( );

    // output_event.Execute(1.0,0.5,1,4,0.6,NULL);

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