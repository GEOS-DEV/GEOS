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
#include "dataRepository/DeprecatedGroup.hpp"
#include "dataRepository/Group.hpp"
#include "common/logger/Logger.hpp"
#include "dataRepository/DataContext.hpp"
#include "mainInterface/initialization.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

using namespace geos;
using namespace dataRepository;

// Test fixture for DeprecatedGroup
class DeprecatedGroupTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    // Redirect cout to capture log output
    m_oldCoutBuf = std::cout.rdbuf();
    std::cout.rdbuf( m_buffer.rdbuf() );
  }

  void TearDown() override
  {
    // Restore cout
    std::cout.rdbuf( m_oldCoutBuf );
  }

  std::string getCapturedOutput()
  {
    return m_buffer.str();
  }

  void clearCapturedOutput()
  {
    m_buffer.str( "" );
    m_buffer.clear();
  }

private:
  std::stringstream m_buffer;
  std::streambuf * m_oldCoutBuf;
};

// Test class that uses DeprecatedGroup with custom message
// DeprecatedGroup<Group> injects itself into the inheritance hierarchy
class TestDeprecatedGroupWithMessage : public DeprecatedGroup< Group >
{
public:
  TestDeprecatedGroupWithMessage( string const & name, Group * const parent )
    : DeprecatedGroup< Group >( name, parent )
  {}

  TestDeprecatedGroupWithMessage( string const & name, conduit::Node & rootNode )
    : DeprecatedGroup< Group >( name, rootNode )
  {}

  virtual string getDeprecationMessage() const override
  {
    return "TestDeprecatedGroupWithMessage is deprecated. Use NewTestGroup instead.";
  }

  static string catalogName() { return "TestDeprecatedGroupWithMessage"; }
};

// Register the deprecated group in the catalog to test catalog compatibility
REGISTER_CATALOG_ENTRY( Group, TestDeprecatedGroupWithMessage, string const &, Group * const );

// Test class that uses DeprecatedGroup without custom message (uses default)
class TestDeprecatedGroupDefaultMessage : public DeprecatedGroup< Group >
{
public:
  TestDeprecatedGroupDefaultMessage( string const & name, Group * const parent )
    : DeprecatedGroup< Group >( name, parent )
  {}

  TestDeprecatedGroupDefaultMessage( string const & name, conduit::Node & rootNode )
    : DeprecatedGroup< Group >( name, rootNode )
  {}
};

// Test that DeprecatedGroup logs a warning on construction
TEST_F( DeprecatedGroupTest, DeprecatedGroupLogsWarning )
{
  conduit::Node node;
  clearCapturedOutput();

  // Create a deprecated group - should log a warning
  TestDeprecatedGroupWithMessage deprecatedGroup( "testDeprecated", node );

  std::string output = getCapturedOutput();

  // Check that output contains deprecation warning markers
  EXPECT_NE( output.find( "DEPRECATION WARNING" ), std::string::npos );
  EXPECT_NE( output.find( "TestDeprecatedGroupWithMessage" ), std::string::npos );
  EXPECT_NE( output.find( "Use NewTestGroup instead" ), std::string::npos );
}

// Test that DeprecatedGroup with default message works
TEST_F( DeprecatedGroupTest, DeprecatedGroupDefaultMessage )
{
  conduit::Node node;
  clearCapturedOutput();

  // Create a deprecated group with default message
  TestDeprecatedGroupDefaultMessage deprecatedGroup( "testDeprecatedDefault", node );

  std::string output = getCapturedOutput();

  // Check that output contains deprecation warning with default message
  EXPECT_NE( output.find( "DEPRECATION WARNING" ), std::string::npos );
  EXPECT_NE( output.find( "TestDeprecatedGroupDefaultMessage" ), std::string::npos );
  EXPECT_NE( output.find( "deprecated and may be removed in a future version" ), std::string::npos );
}

// Test that DeprecatedGroup maintains Group functionality
TEST_F( DeprecatedGroupTest, DeprecatedGroupFunctionality )
{
  conduit::Node node;
  clearCapturedOutput(); // Clear warning output

  TestDeprecatedGroupWithMessage deprecatedGroup( "testDeprecated", node );

  // Test that basic Group operations still work
  EXPECT_EQ( deprecatedGroup.getName(), "testDeprecated" );

  // Test registering wrappers
  deprecatedGroup.registerWrapper< int >( "testInt" ).setDefaultValue( 42 );
  EXPECT_TRUE( deprecatedGroup.hasWrapper( "testInt" ) );
  EXPECT_EQ( deprecatedGroup.getWrapper< int >( "testInt" ).reference(), 42 );

  // Test registering subgroups
  auto & subGroup = deprecatedGroup.registerGroup< Group >( "subGroup" );
  EXPECT_TRUE( deprecatedGroup.hasGroup( "subGroup" ) );
  EXPECT_EQ( &subGroup, deprecatedGroup.getGroupPointer( "subGroup" ) );
}

// Test that DeprecatedGroup can be used with parent constructor
TEST_F( DeprecatedGroupTest, DeprecatedGroupWithParentConstructor )
{
  conduit::Node rootNode;
  Group parentGroup( "parent", rootNode );

  clearCapturedOutput();

  // Create deprecated group as child of another group
  TestDeprecatedGroupWithMessage deprecatedChild( "deprecatedChild", &parentGroup );

  std::string output = getCapturedOutput();

  // Should log warning
  EXPECT_NE( output.find( "DEPRECATION WARNING" ), std::string::npos );

  // Should be properly registered
  EXPECT_EQ( deprecatedChild.getName(), "deprecatedChild" );
  EXPECT_EQ( &deprecatedChild.getParent(), &parentGroup );
}

// Test multiple instantiations of deprecated group
TEST_F( DeprecatedGroupTest, MultipleDeprecatedGroupInstances )
{
  conduit::Node node1;
  conduit::Node node2;

  clearCapturedOutput();

  // Create multiple instances - each should log a warning
  TestDeprecatedGroupWithMessage group1( "deprecated1", node1 );
  std::string output1 = getCapturedOutput();
  EXPECT_NE( output1.find( "DEPRECATION WARNING" ), std::string::npos );

  clearCapturedOutput();

  TestDeprecatedGroupWithMessage group2( "deprecated2", node2 );
  std::string output2 = getCapturedOutput();
  EXPECT_NE( output2.find( "DEPRECATION WARNING" ), std::string::npos );
}

// Test that DeprecatedGroup types are registered in the catalog and can be looked up
TEST_F( DeprecatedGroupTest, DeprecatedGroupCatalogRegistration )
{
  // Check that the deprecated type is registered in the catalog
  EXPECT_TRUE( Group::CatalogInterface::hasKeyName( "TestDeprecatedGroupWithMessage" ) );

  // Verify we can get the list of catalog keys and our deprecated type is present
  auto keys = Group::CatalogInterface::getKeys();
  bool found = false;
  for( auto const & key : keys )
  {
    if( key == "TestDeprecatedGroupWithMessage" )
    {
      found = true;
      break;
    }
  }
  EXPECT_TRUE( found );

  // Test that we can construct a deprecated group via the catalog factory
  conduit::Node node;
  Group parentGroup( "parent", node );

  clearCapturedOutput();

  // Use the catalog factory to create the deprecated group
  DataFileContext const context = DataFileContext( "Test Deprecated Catalog", __FILE__, __LINE__ );
  std::unique_ptr< Group > deprecatedGroup = Group::CatalogInterface::factory(
    "TestDeprecatedGroupWithMessage",
    context,
    "catalogCreatedDeprecated",
    &parentGroup
  );

  // Verify the group was created
  ASSERT_NE( deprecatedGroup, nullptr );
  EXPECT_EQ( deprecatedGroup->getName(), "catalogCreatedDeprecated" );

  // Verify warning was logged during construction
  std::string output = getCapturedOutput();
  EXPECT_NE( output.find( "DEPRECATION WARNING" ), std::string::npos );
  EXPECT_NE( output.find( "TestDeprecatedGroupWithMessage" ), std::string::npos );

  clearCapturedOutput();

  // Verify that the deprecated group still functions normally
  deprecatedGroup->registerWrapper< int >( "testInt" ).setDefaultValue( 99 );
  EXPECT_TRUE( deprecatedGroup->hasWrapper( "testInt" ) );
  EXPECT_EQ( deprecatedGroup->getWrapper< int >( "testInt" ).reference(), 99 );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  geos::basicSetup( argc, argv, false );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
