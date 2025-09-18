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

#include "fileIO/timeHistory/ScatterDataHistoryCollection.hpp"
#include "fileIO/Outputs/TimeHistoryOutput.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/initialization.hpp"
#include <gtest/gtest.h>
#include <conduit.hpp>
#include <hdf5.h>
#include <filesystem>

using namespace geos;

// Minimal mock ScatterDataProvider
class MockScatterDataProvider : public ScatterDataProvider, public dataRepository::Group
{
public:
  MockScatterDataProvider( string const & name, dataRepository::Group *parent ): dataRepository::Group( name, parent )
  {
    m_scatterData.resize( 3 );
    m_scatterData[0] = 10; m_scatterData[1] = 20; m_scatterData[2] = 30;
    m_coordinates.resize( 3, 3 );
    m_coordinates( 0, 0 ) = 1000.0; m_coordinates( 0, 1 ) = 2000.0; m_coordinates( 0, 2 ) = 100.0;
    m_coordinates( 1, 0 ) = 1500.0; m_coordinates( 1, 1 ) = 2500.0; m_coordinates( 1, 2 ) = 200.0;
    m_coordinates( 2, 0 ) = 2000.0; m_coordinates( 2, 1 ) = 3000.0; m_coordinates( 2, 2 ) = 300.0;
    m_metadata.resize( 3 );
    m_metadata[0] = "Station_A"; m_metadata[1] = "Station_B"; m_metadata[2] = "Station_C";
  }
  localIndex getNumScatterPoints() const override { return 3; }
  array1d< real64 > const & getScatterData() const override { return m_scatterData; }
  array2d< real64 > const & getScatterCoordinates() const override { return m_coordinates; }
  string_array const & getScatterMetadata() const override { return m_metadata; }
  static string catalogName() { return "MockScatterDataProvider"; }
private:
  array1d< real64 > m_scatterData;
  array2d< real64 > m_coordinates;
  string_array m_metadata;
};

REGISTER_CATALOG_ENTRY( dataRepository::Group, MockScatterDataProvider, string const &, dataRepository::Group * const )


class ScatterDataHistoryIntegrationTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    m_conduitNode = std::make_unique< conduit::Node >();
    m_rootGroup = std::make_unique< dataRepository::Group >( "Problem", *m_conduitNode );
    auto & solvers = m_rootGroup->registerGroup< PhysicsSolverManager >( "Solvers" );
    m_domain = &m_rootGroup->registerGroup< DomainPartition >( "domain" );
    solvers.registerGroup< MockScatterDataProvider >( "gravitySolver" );
    m_collection = &m_rootGroup->registerGroup< ScatterDataHistoryCollection >( "gravityHistory" );
    m_collection->getReference< string >( "solverName" ) = "gravitySolver";
    m_collection->getReference< integer >( "includeCoordinates" ) = 1;
    m_collection->getReference< integer >( "includeMetadata" ) = 0;
  }
  std::unique_ptr< conduit::Node > m_conduitNode;
  std::unique_ptr< dataRepository::Group > m_rootGroup;
  DomainPartition *m_domain;
  ScatterDataHistoryCollection *m_collection;
};

TEST_F( ScatterDataHistoryIntegrationTest, CompleteWorkflow )
{
  // Test initialization
  m_collection->initializePostInitialConditionsPreSubGroups();
  EXPECT_EQ( m_collection->numCollectors(), 4 ); // data + 3 coordinates
  EXPECT_EQ( m_collection->getTargetName(), "gravitySolver" );
  EXPECT_EQ( m_collection->catalogName(), "ScatterDataHistoryCollection" );
  EXPECT_EQ( m_collection->numMetaDataCollectors(), 0 );
  EXPECT_THROW( m_collection->getMetaDataCollector( 0 ), std::runtime_error );

  // Test metadata for all collectors
  for( localIndex idx = 0; idx < 4; ++idx )
  {
    auto metadata = m_collection->getMetaData( *m_domain, idx );
    EXPECT_EQ( metadata.size(), 3 );
    EXPECT_EQ( metadata.getType(), std::type_index( typeid(real64)));
    if( idx == 0 )
      EXPECT_EQ( metadata.getName(), "scatterData" );
    else
      EXPECT_TRUE( metadata.getName().find( "coordinate" ) != string::npos );
  }
}

TEST_F( ScatterDataHistoryIntegrationTest, ErrorHandling )
{
  m_collection->getReference< string >( "solverName" ) = "invalidSolver";
  EXPECT_THROW( m_collection->initializePostInitialConditionsPreSubGroups(), InputError );
}

// HDF5 test
class ScatterDataHistoryHDF5Test : public ::testing::Test
{
protected:
  void SetUp() override
  {
    m_conduitNode = std::make_unique< conduit::Node >();
    m_rootGroup = std::make_unique< dataRepository::Group >( "Problem", *m_conduitNode );
    auto & solvers = m_rootGroup->registerGroup< PhysicsSolverManager >( "Solvers" );
    m_provider = &solvers.registerGroup< MockScatterDataProvider >( "gravitySolver" );
    m_domain = &m_rootGroup->registerGroup< DomainPartition >( "domain" );
    auto & outputs = m_rootGroup->registerGroup< dataRepository::Group >( "Outputs" );
    m_timeHistoryOutput = &outputs.registerGroup< TimeHistoryOutput >( "timeHistory" );
    m_timeHistoryOutput->getReference< string >( "filename" ) = "testOutput";  // No path, just filename
    m_timeHistoryOutput->getReference< string >( "format" ) = "hdf5";
    m_collection = &m_timeHistoryOutput->registerGroup< ScatterDataHistoryCollection >( "gravityHistory" );
    m_collection->getReference< string >( "solverName" ) = "gravitySolver";
    m_collection->getReference< integer >( "includeCoordinates" ) = 1;
    // Register the collection as a source for TimeHistoryOutput using string_array
    m_timeHistoryOutput->getReference< string_array >( "sources" ).resize( 1 );
    m_timeHistoryOutput->getReference< string_array >( "sources" )[0] = "gravityHistory";
  }
  void TearDown() override
  {
    std::filesystem::remove( "testOutput.hdf5" );
  }
  std::unique_ptr< conduit::Node > m_conduitNode;
  std::unique_ptr< dataRepository::Group > m_rootGroup;
  DomainPartition *m_domain;
  TimeHistoryOutput *m_timeHistoryOutput;
  ScatterDataHistoryCollection *m_collection;
  MockScatterDataProvider *m_provider;
};

TEST_F( ScatterDataHistoryHDF5Test, HDF5FileOutput )
{
  // Set up output directory for TimeHistoryOutput
  OutputBase::setOutputDirectory( "." );

  // Initialize collection first
  m_collection->initializePostInitialConditionsPreSubGroups();

  // Debug: Check the sources array is properly set
  auto & sources = m_timeHistoryOutput->getReference< string_array >( "sources" );
  EXPECT_EQ( sources.size(), 1 );
  EXPECT_EQ( sources[0], "gravityHistory" );

  // Initialize TimeHistoryOutput - this should now work without directory errors
  EXPECT_NO_THROW( m_timeHistoryOutput->initializePostInitialConditionsPostSubGroups() );

  // Verify the collection is configured correctly
  EXPECT_EQ( m_collection->numCollectors(), 4 ); // data + 3 coordinates
  EXPECT_EQ( m_collection->getTargetName(), "gravitySolver" );

  // Verify that we can get metadata for all collectors
  for( localIndex idx = 0; idx < 4; ++idx )
  {
    auto metadata = m_collection->getMetaData( *m_domain, idx );
    EXPECT_EQ( metadata.size(), 3 );  // 3 scatter points
    EXPECT_EQ( metadata.getType(), std::type_index( typeid(real64) ) );
    if( idx == 0 )
      EXPECT_EQ( metadata.getName(), "scatterData" );
    else
      EXPECT_TRUE( metadata.getName().find( "coordinate" ) != string::npos );
  }

  // Trigger buffer population before writing
  m_collection->collect( *m_domain );

  // Now test the actual I/O execution!
  EXPECT_NO_THROW( m_timeHistoryOutput->execute( 0.0, 1.0, 0, 1, 1.0, *m_domain ) );

  // Verify HDF5 file was created
  EXPECT_TRUE( std::filesystem::exists( "testOutput.hdf5" ) );

  // Read back the HDF5 file and verify the contents
  hid_t file_id = H5Fopen( "testOutput.hdf5", H5F_ACC_RDONLY, H5P_DEFAULT );
  EXPECT_GE( file_id, 0 );

  // Check scatter data values
  {
    const char * dataset_names[] = { "/scatterData", "/coordinateX", "/coordinateY", "/coordinateZ" };
    double values[3];
    auto & scatterData = m_provider->getScatterData();
    auto & coordinates = m_provider->getScatterCoordinates();
    for( int i = 0; i < 4; ++i )
    {
      hid_t dataset_id = H5Dopen2( file_id, dataset_names[i], H5P_DEFAULT );
      EXPECT_GE( dataset_id, 0 ) << "Failed to open dataset " << dataset_names[i];
      herr_t status = H5Dread( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values );
      EXPECT_GE( status, 0 ) << "Failed to read dataset " << dataset_names[i];
      for( int j = 0; j < 3; ++j )
      {
        double expected = 0.0;
        if( i == 0 )
          expected = scatterData[j];
        else
          expected = coordinates( j, i-1 );
        EXPECT_DOUBLE_EQ( values[j], expected ) << "Mismatch in " << dataset_names[i] << " at index " << j;
      }
      H5Dclose( dataset_id );
    }
  }
  H5Fclose( file_id );
}

int main( int argc, char * *argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  basicSetup( argc, argv );
  int result = RUN_ALL_TESTS();
  basicCleanup();
  return result;
}
