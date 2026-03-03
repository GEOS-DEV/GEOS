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

// Source includes

// #include "common/GeosxConfig.hpp"

// // #include "BufferAllocator.hpp" 
// #include "common/DataLayouts.hpp"
// #include "common/GeosxMacros.hpp"
// #include "common/Path.hpp"
// #include "common/StdContainerWrappers.hpp"
// #include "common/Tensor.hpp"

// #include "LvArray/src/Macros.hpp"
// #include "LvArray/src/Array.hpp"
// #include "LvArray/src/ArrayOfArrays.hpp"
// #include "LvArray/src/ArrayOfSets.hpp"
// #include "LvArray/src/SparsityPattern.hpp"
// #include "LvArray/src/CRSMatrix.hpp"
// #include "LvArray/src/SortedArray.hpp"
// #include "LvArray/src/StackBuffer.hpp"
// #include "LvArray/src/ChaiBuffer.hpp"


// // TPL includes
// #include <camp/camp.hpp>

// // System includes
// #ifdef GEOS_USE_MPI
//   #include <mpi.h>
// #endif

// #include <cassert>
// //#include <cmath>
// #include <cstdint>
// #include <iostream>
// #include <optional>
// #include <set>
// #include <string>
// #include <string_view>
// #include <typeindex>
// #include <typeinfo>

// // TPL includes
#include <gtest/gtest.h>
#include <gtest/gtest-spi.h>

// using namespace geos;

// TEST( testDataTypes, testBoundChecking )
// {
//   internal::StdVectorWrapper< std::string,
//                               internal::DefaultAllocator< std::string >,
//                               true > vectorBoundsChecking = {"test"};
//   EXPECT_THROW( {
//     try
//     {
//       for( integer i = 0; i <= 1; i++ )
//       {
//         std::string crash = vectorBoundsChecking[i];
//         GEOS_UNUSED_VAR( crash );
//       }
//     }
//     catch( const std::out_of_range & e )
//     {
//       throw;
//     }
//   }, std::out_of_range );

//   internal::StdMapWrapper< map< integer, integer >, true > mapBoundsChecking{{0, 1}};
//   EXPECT_THROW( {
//     try
//     {
//       for( integer i = 0; i <= 1; i++ )
//       {
//         integer crash = mapBoundsChecking[i];
//         GEOS_UNUSED_VAR( crash );
//       }
//     }
//     catch( const std::out_of_range & e )
//     {
//       throw;
//     }
//   }, std::out_of_range );


//   internal::StdMapWrapper< std::unordered_map< integer, integer >, true > unorderedMapBoundsChecking{{0, 1}};
//   EXPECT_THROW( {
//     try
//     {
//       for( integer i = 0; i <= 1; i++ )
//       {
//         integer crash = unorderedMapBoundsChecking[i];
//         GEOS_UNUSED_VAR( crash );
//       }
//     }
//     catch( const std::out_of_range & e )
//     {
//       throw;
//     }
//   }, std::out_of_range );



// }

int main( int argc, char * * argv )
{
  testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
