/*	
 * ------------------------------------------------------------------------------------------------------------	
 * SPDX-License-Identifier: LGPL-2.1-only	
 *	
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC	
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University	
 * Copyright (c) 2018-2020 Total, S.A	
 * Copyright (c) 2020-     GEOSX Contributors	
 * All right reserved	
 *	
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.	
 * ------------------------------------------------------------------------------------------------------------	
 */

#include "../VTKPolyDataWriterInterface.hpp"

#include "mesh/ElementRegionManager.hpp"
#include "mesh/NodeManager.hpp"

#include "managers/initialization.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace geosx;

// TODO implement the ElementRegionManagerABC, NodeManagerABC
// TODO VTKPolyDataWriterInterface must compile with these...
class MockElementRegionManager: public ElementRegionManager
{
public:
  MOCK_METHOD(
    (std::list< std::reference_wrapper< const CellElementRegionABC > >),
    getCellElementRegions,
    (),
    (const, override)
  );
};

class MockNodeManager: public NodeManager
{
public:

};

TEST( testVTK, dieIfFileAlreadyExists )
{
//  geosx::vtk::VTKPolyDataWriterInterface writer( "/I/cannot/be/created/simply" );

//  auto instanciate = []() {
//    geosx::vtk::VTKPolyDataWriterInterface writer( "/I/cannot/be/created/simply" );
//  };

  testing::FLAGS_gtest_death_test_style="fast";
//  testing::FLAGS_gtest_death_test_style="threadsafe";

  auto instanciate2 = []() {
//    std::cerr << "/I/cannot/be/created/simply" << std::endl ;
//    exit(1);
    geosx::vtk::VTKPolyDataWriterInterface writer( "/I/cannot/be/created/simply" );
  };

//  const real64 time = 0;
//  const integer cycle = 0;
  const MockElementRegionManager elementRegionManager;
  const MockNodeManager nodeManager;
  using ::testing::HasSubstr;
//  writer.Write( time, cycle, elementRegionManager, nodeManager );
  // I cannot use the msg match since our error logs sends on stdout and not stderr.
  ASSERT_DEATH( instanciate2(), HasSubstr("/I/cannot/be/created/simply"));
//  ASSERT_DEATH( writer.Write( time, cycle, elementRegionManager, nodeManager ), "HasSubstr(already)");
}

//TEST( testVTK, doIt )
//{
//  EXPECT_TRUE(true);
//
//  geosx::vtk::VTKPolyDataWriterInterface writer("/tmp/test");
//
//
//  const real64 time = 0;
//  const integer cycle = 0;
//
//  const MockElementRegionManager elementRegionManager;
//  const MockNodeManager nodeManager;
//
//  writer.Write(time, cycle, elementRegionManager, nodeManager);
//}

int main( int ac, char * av[] )
{
  ::testing::InitGoogleTest( &ac, av );
  geosx::basicSetup( ac, av );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
