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
#include "dataRepository/ObsoleteGroup.hpp"
#include "dataRepository/Group.hpp"
#include "common/logger/Logger.hpp"
#include "dataRepository/DataContext.hpp"
#include "mainInterface/initialization.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

using namespace geos;
using namespace dataRepository;

// Test fixture for ObsoleteGroup and NoopWrapper
class ObsoleteGroupTest : public ::testing::Test
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

// Test class that uses ObsoleteGroup with custom message
// ObsoleteGroup<Group> injects itself into the inheritance hierarchy
class TestObsoleteGroupWithMessage : public ObsoleteGroup< Group >
{
public:
  TestObsoleteGroupWithMessage( string const & name, Group * const parent )
    : ObsoleteGroup< Group >( name, parent )
  {}

  TestObsoleteGroupWithMessage( string const & name, conduit::Node & rootNode )
    : ObsoleteGroup< Group >( name, rootNode )
  {}

  virtual string getObsoleteMessage() const override
  {
    return "TestObsoleteGroupWithMessage has been removed. Use ReplacementGroup instead.";
  }

  static string catalogName() { return "TestObsoleteGroupWithMessage"; }
};

// Register the obsolete group in the catalog to test catalog compatibility
REGISTER_CATALOG_ENTRY( Group, TestObsoleteGroupWithMessage, string const &, Group * const );

// Test that ObsoleteGroup allows construction but logs warning
TEST_F( ObsoleteGroupTest, ObsoleteGroupAllowsConstructionWithWarning )
{
  conduit::Node node;
  clearCapturedOutput();

  // Creating an obsolete group should succeed (for parsing) but log warning
  EXPECT_NO_THROW(
    {
      TestObsoleteGroupWithMessage obsoleteGroup( "testObsolete", node );

      std::string output = getCapturedOutput();

      // Check that warning was logged
      EXPECT_NE( output.find( "OBSOLETE GROUP WARNING" ), std::string::npos );
      EXPECT_NE( output.find( "TestObsoleteGroupWithMessage" ), std::string::npos );
      EXPECT_NE( output.find( "Use ReplacementGroup instead" ), std::string::npos );
      EXPECT_NE( output.find( "input parsing compatibility only" ), std::string::npos );
    }
  );
}

// Test that ObsoleteGroup throws error on initialization
TEST_F( ObsoleteGroupTest, ObsoleteGroupThrowsErrorOnInitialization )
{
  conduit::Node node;
  clearCapturedOutput();

  TestObsoleteGroupWithMessage obsoleteGroup( "testObsolete", node );

  // Clear the construction warning
  clearCapturedOutput();

  // Trying to initialize should throw an error (initialize() calls initializePreSubGroups())
  EXPECT_THROW(
    {
      obsoleteGroup.initialize();
    },
    std::exception
  );
}

// Test that ObsoleteGroup uses NoopWrappers for registration
TEST_F( ObsoleteGroupTest, ObsoleteGroupUsesNoopWrappers )
{
  conduit::Node node;
  clearCapturedOutput();

  TestObsoleteGroupWithMessage obsoleteGroup( "testObsolete", node );
  clearCapturedOutput(); // Clear construction warning

  // Register a wrapper - should be NoopWrapper
  auto & wrapper = obsoleteGroup.registerWrapper< int >( "testInt" );

  // The wrapper should exist but have no memory
  EXPECT_TRUE( obsoleteGroup.hasWrapper( "testInt" ) );
  WrapperBase const & wrapperBase = obsoleteGroup.getWrapperBase( "testInt" );
  EXPECT_EQ( wrapperBase.size(), 0 );
  EXPECT_EQ( wrapperBase.bytesAllocated(), 0 );
}

// Test NoopWrapper size and memory methods
TEST_F( ObsoleteGroupTest, NoopWrapperHasNoMemoryFootprint )
{
  conduit::Node node;
  Group testGroup( "test", node );

  NoopWrapper< int > noopWrapper( "noopInt", testGroup );

  // Check that NoopWrapper reports no memory usage
  EXPECT_EQ( noopWrapper.size(), 0 );
  EXPECT_EQ( noopWrapper.bytesAllocated(), 0 );
  EXPECT_EQ( noopWrapper.capacity(), 0 );
  EXPECT_EQ( noopWrapper.voidPointer(), nullptr );
  EXPECT_EQ( noopWrapper.elementByteSize(), sizeof( int ) );
}

// Test NoopWrapper operations are no-ops
TEST_F( ObsoleteGroupTest, NoopWrapperOperationsAreNoOps )
{
  conduit::Node node;
  Group testGroup( "test", node );

  NoopWrapper< int > noopWrapper( "noopInt", testGroup );

  // These should not throw and should do nothing
  EXPECT_NO_THROW( noopWrapper.resize( 100 ) );
  EXPECT_NO_THROW( noopWrapper.reserve( 200 ) );
  EXPECT_NO_THROW( noopWrapper.copy( 0, 1 ) );
  EXPECT_NO_THROW( noopWrapper.erase( std::set< localIndex >{ 0, 1, 2 } ) );

  // Size should still be 0
  EXPECT_EQ( noopWrapper.size(), 0 );
  EXPECT_EQ( noopWrapper.capacity(), 0 );
}

// Test NoopWrapper packing operations
TEST_F( ObsoleteGroupTest, NoopWrapperPackingReturnsZero )
{
  conduit::Node node;
  Group testGroup( "test", node );

  NoopWrapper< int > noopWrapper( "noopInt", testGroup );

  // Check that packing operations return 0
  EXPECT_FALSE( noopWrapper.isPackable( false ) );
  EXPECT_FALSE( noopWrapper.isPackable( true ) );

  parallelDeviceEvents events;
  buffer_unit_type * buffer = nullptr;
  buffer_unit_type const * constBuffer = nullptr;

  EXPECT_EQ( noopWrapper.pack< false >( buffer, false, false, events ), 0 );
  EXPECT_EQ( noopWrapper.pack< true >( buffer, false, false, events ), 0 );
  EXPECT_EQ( noopWrapper.unpack( constBuffer, false, false, events ), 0 );
}

// Test NoopWrapper type information
TEST_F( ObsoleteGroupTest, NoopWrapperTypeInformation )
{
  conduit::Node node;
  Group testGroup( "test", node );

  NoopWrapper< int > noopWrapper( "noopInt", testGroup );

  // Check type information
  EXPECT_EQ( noopWrapper.getTypeId(), typeid( int ) );
  EXPECT_EQ( noopWrapper.numArrayDims(), 0 );
  EXPECT_EQ( noopWrapper.numArrayComp(), 0 );
  EXPECT_FALSE( noopWrapper.hasDefaultValue() );
  EXPECT_EQ( noopWrapper.getDefaultValueString(), "" );
}

// Test NoopWrapper cloning
TEST_F( ObsoleteGroupTest, NoopWrapperCloning )
{
  conduit::Node node;
  Group testGroup( "test", node );

  NoopWrapper< int > noopWrapper( "noopInt", testGroup );

  // Clone the wrapper
  auto cloned = noopWrapper.clone( "clonedNoop", testGroup );

  EXPECT_NE( cloned.get(), nullptr );
  EXPECT_EQ( cloned->getName(), "clonedNoop" );
  EXPECT_EQ( cloned->size(), 0 );
}

// Test that ObsoleteGroup types are registered in the catalog and can be looked up
TEST_F( ObsoleteGroupTest, ObsoleteGroupCatalogRegistration )
{
  // Check that the obsolete type is registered in the catalog
  EXPECT_TRUE( Group::CatalogInterface::hasKeyName( "TestObsoleteGroupWithMessage" ) );

  // Verify we can get the list of catalog keys and our obsolete type is present
  auto keys = Group::CatalogInterface::getKeys();
  bool found = false;
  for( auto const & key : keys )
  {
    if( key == "TestObsoleteGroupWithMessage" )
    {
      found = true;
      break;
    }
  }
  EXPECT_TRUE( found );

  // Test that we can construct an obsolete group via the catalog factory
  conduit::Node node;
  Group parentGroup( "parent", node );

  clearCapturedOutput();

  // Use the catalog factory to create the obsolete group
  DataFileContext const context = DataFileContext( "Test Obsolete Catalog", __FILE__, __LINE__ );
  std::unique_ptr< Group > obsoleteGroup = Group::CatalogInterface::factory(
    "TestObsoleteGroupWithMessage",
    context,
    "catalogCreatedObsolete",
    &parentGroup
  );

  // Verify the group was created
  ASSERT_NE( obsoleteGroup, nullptr );
  EXPECT_EQ( obsoleteGroup->getName(), "catalogCreatedObsolete" );

  // Verify warning was logged during construction
  std::string output = getCapturedOutput();
  EXPECT_NE( output.find( "OBSOLETE GROUP WARNING" ), std::string::npos );
  EXPECT_NE( output.find( "TestObsoleteGroupWithMessage" ), std::string::npos );

  clearCapturedOutput();

  // Verify that trying to initialize still throws error (initialize() calls initializePreSubGroups())
  EXPECT_THROW(
    {
      obsoleteGroup->initialize();
    },
    std::exception
  );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  basicSetup( argc, argv, false );
  int const result = RUN_ALL_TESTS();
  basicCleanup();
  return result;
}
