/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// TPL includes
#include <gtest/gtest.h>

// Source includes
#include "codingUtilities/UnitTestUtilities.hpp"
#include "dataRepository/xmlWrapper.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/generators/CellBlockManagerABC.hpp"
#include "mesh/generators/CellBlockABC.hpp"

// special CMake-generated include
#include "tests/meshDirName.hpp"

using namespace geos;
using namespace geos::testing;
using namespace geos::dataRepository;


CommandLineOptions g_commandLineOptions;

struct TestCase
{
  string name;
  bool isExpectedToPass = true;
  std::vector< string > stringsToMention;
  string_view xmlRegions;
};

class ElementRegionTestFixture : public ::testing::TestWithParam< TestCase >
{
public:
  ElementRegionTestFixture():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

  virtual ~ElementRegionTestFixture() = default;
private:
  GeosxState state;
};

TEST_P( ElementRegionTestFixture, testVTKImportRegionSyntaxes )
{
  TestCase const testCase = GetParam();

  string const pattern =
    R"xml(
      <Problem>
        <Mesh>
          <VTKMesh name="mesh" 
                   file="{}" />
        </Mesh>
        <ElementRegions>
          {}
        </ElementRegions>
      </Problem>
    )xml";
  string const xmlInput = GEOS_FMT( pattern,
                                    testMeshDir + "/box_hybrid_mesh.vtu",
                                    testCase.xmlRegions );

  ProblemManager & problem = getGlobalState().getProblemManager();
  problem.parseInputString( xmlInput );

  if( testCase.isExpectedToPass )
  {
    try
    {
      EXPECT_NO_FATAL_FAILURE( problem.problemSetup() ) << GEOS_FMT( "Test case '{}' did throw an error.",
                                                                     testCase.name );
    }
    catch( std::exception const & e )
    {
      GTEST_FAIL() << GEOS_FMT( "Test case '{}' did throw an exception ({}) but was not expected to :\n{}",
                                testCase.name, LvArray::system::demangle( typeid( e ).name() ), e.what() );
    }
  }
  else
  {
    try
    {
      EXPECT_NO_FATAL_FAILURE( problem.problemSetup() ) << GEOS_FMT( "Test case '{}' did throw an error.",
                                                                     testCase.name );
      GTEST_FAIL() << GEOS_FMT( "Test case '{}' did not thrown any exception but was expected to.",
                                testCase.name );
    }
    catch( InputError const & e )
    {
      string const expStr = e.what();
      for( auto const & str : testCase.stringsToMention )
      {
        bool isExceptionContainingStr = expStr.find( str ) != string::npos;
        EXPECT_TRUE( isExceptionContainingStr ) << GEOS_FMT( "Test case '{}' exception did not mention the string '{}'. Exception string:\n{}",
                                                             testCase.name, str, e.what() );
      }
    }
    catch( std::exception const & e )
    {
      GTEST_FAIL() << GEOS_FMT( "Test case '{}' did throw an exception ({}) but was not expected to :\n{}",
                                testCase.name, LvArray::system::demangle( typeid( e ).name() ), e.what() );
    }
    catch( ... )
    {
      GTEST_FAIL() << GEOS_FMT( "Test case '{}': Unexpected exception.",
                                testCase.name );
    }
  }
}

TestCase const vtkImportRegionSyntaxCases[] = {
  { // should not crash
    "regular subRegion list", true, { },
    R"xml(
      <CellElementRegion name="overburden" materialList="{ }"
                         cellBlocks="{ 3_hexahedra, 3_tetrahedra, 3_pyramids, 5_hexahedra, 5_tetrahedra, 5_pyramids }" />
      <CellElementRegion name="reservoir" materialList="{ }"
                         cellBlocks="{ 1_hexahedra, 1_tetrahedra, 1_pyramids, 6_hexahedra, 6_tetrahedra, 6_pyramids }" />
      <CellElementRegion name="underburden" materialList="{ }"
                         cellBlocks="{ 2_hexahedra, 2_tetrahedra, 2_pyramids, 4_hexahedra, 4_tetrahedra, 4_pyramids }" />
    )xml"
  },
  { // crash because of not existing primitive (3_pendecagonalPrism)
    "non existing 3_pendecagonalPrism", false, { "overburden", "3_pendecagonalPrism" },
    R"xml(
      <CellElementRegion name="overburden" materialList="{ }" 
                         cellBlocks="{ 3_pendecagonalPrism, 3_hexahedra, 3_tetrahedra, 3_pyramids, 5_hexahedra, 5_tetrahedra, 5_pyramids }" />
      <CellElementRegion name="reservoir" materialList="{ }"
                         cellBlocks="{ 1_hexahedra, 1_tetrahedra, 1_pyramids, 6_hexahedra, 6_tetrahedra, 6_pyramids }" />
      <CellElementRegion name="underburden" materialList="{ }"
                         cellBlocks="{ 2_hexahedra, 2_tetrahedra, 2_pyramids, 4_hexahedra, 4_tetrahedra, 4_pyramids }" />
    )xml"
  },
  { // crash because of not unknown primitive name (helloWorld)
    "non existing helloWorld", false, { "overburden", "helloWorld" },
    R"xml(
      <CellElementRegion name="overburden" materialList="{ }" 
                          cellBlocks="{ helloWorld, 3_hexahedra, 3_tetrahedra, 3_pyramids, 5_hexahedra, 5_tetrahedra, 5_pyramids }" />
      <CellElementRegion name="reservoir" materialList="{ }"
                          cellBlocks="{ 1_hexahedra, 1_tetrahedra, 1_pyramids, 6_hexahedra, 6_tetrahedra, 6_pyramids }" />
      <CellElementRegion name="underburden" materialList="{ }"
                          cellBlocks="{ 2_hexahedra, 2_tetrahedra, 2_pyramids, 4_hexahedra, 4_tetrahedra, 4_pyramids }" />
    )xml"
  },
  { // crash because of lacking one primitive (3_hexahedra)
    "lacking 3_hexahedra", false, { "3_hexahedra" },
    R"xml(
      <CellElementRegion name="overburden" materialList="{ }"
                          cellBlocks="{ 3_tetrahedra, 3_pyramids, 5_hexahedra, 5_tetrahedra, 5_pyramids }" />
      <CellElementRegion name="reservoir" materialList="{ }"
                          cellBlocks="{ 1_hexahedra, 1_tetrahedra, 1_pyramids, 6_hexahedra, 6_tetrahedra, 6_pyramids }" />
      <CellElementRegion name="underburden" materialList="{ }"
                          cellBlocks="{ 2_hexahedra, 2_tetrahedra, 2_pyramids, 4_hexahedra, 4_tetrahedra, 4_pyramids }" />
    )xml"
  },
  { // mentioning the same sub-regions in multiple cellBlocks (3_hexahedra)
    "multiple 3_hexahedra", false, { "overburden", "reservoir", "3_hexahedra" },
    R"xml(
      <CellElementRegion name="overburden" materialList="{ }"
                          cellBlocks="{ 3_hexahedra, 3_tetrahedra, 3_pyramids, 5_hexahedra, 5_tetrahedra, 5_pyramids }" />
      <CellElementRegion name="reservoir" materialList="{ }"
                          cellBlocks="{ 3_hexahedra, 1_hexahedra, 1_tetrahedra, 1_pyramids, 6_hexahedra, 6_tetrahedra, 6_pyramids }" />
      <CellElementRegion name="underburden" materialList="{ }"
                          cellBlocks="{ 2_hexahedra, 2_tetrahedra, 2_pyramids, 4_hexahedra, 4_tetrahedra, 4_pyramids }" />
    )xml"
  },
  { // should not crash
    "regular region list", true, { },
    R"xml(
      <CellElementRegion name="overburden" materialList="{ }" cellBlockAttributeValues="{ 3, 5 }" />
      <CellElementRegion name="reservoir" materialList="{ }" cellBlockAttributeValues="{ 1, 6 }" />
      <CellElementRegion name="underburden" materialList="{ }" cellBlockAttributeValues="{ 2, 4 }" />
    )xml"
  },
  { // mentioning the same region in multiple cellBlocks (6 in overburden & reservoir)
    "multiple region 1", false, { "6", "overburden", "reservoir" },
    R"xml(
      <CellElementRegion name="overburden" materialList="{ }" cellBlockAttributeValues="{ 3, 5, 6 }" />
      <CellElementRegion name="reservoir" materialList="{ }" cellBlockAttributeValues="{ 1, 6 }" />
      <CellElementRegion name="underburden" materialList="{ }" cellBlockAttributeValues="{ 2, 4 }" />
    )xml"
  },
  { // should not crash
    "regular * wildcard", true, { },
    R"xml(
      <CellElementRegion name="everything" materialList="{ }" cellBlocksMatch="{ * }" />
    )xml"
  },
  { // mentioning the same regions in multiple cellBlocks (because of "*")
    "* wildcard + region list", false, { "everything" },
    R"xml(
      <CellElementRegion name="everything" materialList="{ }" cellBlocksMatch="{ * }" />
      <CellElementRegion name="reservoir" materialList="{ }" cellBlockAttributeValues="{ 1, 6 }" />
    )xml"
  },
  { // using multiple selection methods on the same region
    "multiple selection methods", false, { "everything", "cellBlocksMatch", "cellBlockAttributeValues" },
    R"xml(
      <CellElementRegion name="everything" materialList="{ }" cellBlocksMatch="{ * }" cellBlockAttributeValues="{ 1, 6 }" />
    )xml"
  },
  { // should not crash
    "mixing all selection methods", true, { },
    R"xml(
      <CellElementRegion name="overburden" materialList="{ }" cellBlockAttributeValues="{ 3, 5 }" />
      <CellElementRegion name="reservoir" materialList="{ }" cellBlocksMatch="{ 1_*, 6_* }" />
      <CellElementRegion name="underburden" materialList="{ }"
                         cellBlocks="{ 2_hexahedra, 2_tetrahedra, 2_pyramids, 4_hexahedra, 4_tetrahedra, 4_pyramids }" />
    )xml"
  }
};
INSTANTIATE_TEST_SUITE_P( testElementRegions, ElementRegionTestFixture,
                          ::testing::ValuesIn( vtkImportRegionSyntaxCases ) );


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
