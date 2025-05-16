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

#include "common/DataTypes.hpp"
#include "common/initializeEnvironment.hpp"
#include "codingUtilities/UnitTestUtilities.hpp"
#include "gtest/gtest.h"
#include "functions/HDF5Utilities.hpp"

using namespace geos;
using namespace geos::hdf5Utils;

#include <hdf5.h>

namespace geos
{

namespace testing
{

class HDF5UtilsTestScope
{
public:

  HDF5UtilsTestScope( int argc, char * * argv )
  {
    ::testing::InitGoogleTest( &argc, argv );
    geos::setupEnvironment( argc, argv );
  }

  ~HDF5UtilsTestScope()
  {
    geos::cleanupEnvironment();
  }
};

}
}

// Test Fixture for SerialHDF5File
class SerialHDF5FileTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    // Setup code before each test
    testFilename = "test_file.h5";

    // Create a dummy HDF5 file for testing
    hid_t fileId = H5Fcreate( testFilename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
    ASSERT_GE( fileId, 0 );
    H5Fclose( fileId );
  }

  void TearDown() override
  {
    // Cleanup code after each test
    ASSERT_EQ( std::remove( testFilename.c_str()), 0 ) << "Failed to delete file " << testFilename;
  }

  std::string testFilename;
};

// Test for opening and closing an HDF5 file
TEST_F( SerialHDF5FileTest, OpenAndCloseFile )
{
  hid_t file_id = -1;
  {
    SerialHDF5File file( testFilename );

    // Check that the file ID is valid
    file_id = file.getFileId();
    ASSERT_GE( file_id, 0 );
    EXPECT_TRUE( H5Iis_valid( file_id ));
  }

  // After file goes out of scope, file is closed
  EXPECT_FALSE( H5Iis_valid( file_id ));
}

// Test for move constructor
TEST_F( SerialHDF5FileTest, MoveConstructor )
{
  SerialHDF5File file1( testFilename );
  SerialHDF5File file2( std::move( file1 ));

  // Check that file1 is invalidated
  EXPECT_EQ( file1.getFileId(), -1 );

  // Check that file2 has a valid file ID
  EXPECT_GE( file2.getFileId(), 0 );
}

// Test for move assignment operator
TEST_F( SerialHDF5FileTest, MoveAssignment )
{
  SerialHDF5File file1( testFilename );
  SerialHDF5File file2 = std::move( file1 );

  // Check that file1 is invalidated
  EXPECT_EQ( file1.getFileId(), -1 );

  // Check that file2 has a valid file ID
  EXPECT_GE( file2.getFileId(), 0 );
}

// Test Fixture for SerialHDF5Reader
class SerialHDF5ReaderTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    // Setup code before each test
    testFilename = "test_file.h5";

    // Create a dummy HDF5 file with a dataset for testing
    hid_t fileId = H5Fcreate( testFilename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
    ASSERT_GE( fileId, 0 );

    hsize_t dims[1] = {10};
    hid_t dataspaceId = H5Screate_simple( 1, dims, nullptr );
    ASSERT_GE( dataspaceId, 0 );

    hid_t datasetId = H5Dcreate2( fileId, "test_dataset", H5T_NATIVE_INT, dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    ASSERT_GE( datasetId, 0 );

    H5Dclose( datasetId );
    H5Sclose( dataspaceId );
    H5Fclose( fileId );
  }

  void TearDown() override
  {
    // Cleanup code after each test
    std::remove( testFilename.c_str());
  }

  std::string testFilename;
};

// Test for reading a non-existent dataset
TEST_F( SerialHDF5ReaderTest, ReadNonExistentDataset )
{
  SerialHDF5Reader reader( testFilename );

  // Attempt to read a non-existent dataset and expect an exception
  EXPECT_THROW( reader.read1DAs< int >( "non_existent_dataset", 1 ), std::runtime_error );
}



// === Test Configuration Struct ===
template< typename DatasetType, typename ArrayType >
struct ReadTestConfig
{
  using dataset_type = DatasetType;
  using array_type = ArrayType;

  static std::string label()
  {
    std::string in, out;
    if constexpr (std::is_same_v< DatasetType, uint8_t >) in = "uint8_t";
    else if constexpr (std::is_same_v< DatasetType, localIndex >) in = "localIndex";
    else if constexpr (std::is_same_v< DatasetType, real32 >) in = "real32";
    else if constexpr (std::is_same_v< DatasetType, real64 >) in = "real64";

    if constexpr (std::is_same_v< ArrayType, uint8_t >) out = "uint8_t";
    else if constexpr (std::is_same_v< ArrayType, localIndex >) out = "localIndex";
    else if constexpr (std::is_same_v< ArrayType, real32 >) out = "real32";
    else if constexpr (std::is_same_v< ArrayType, real64 >) out = "real64";

    return in + "_to_" + out;
  }
};

// === HDF5 Type Mapping ===
template< typename T > hid_t getH5Type();
template<> hid_t getH5Type< uint8_t >()  { return H5T_NATIVE_UINT8; }
template<> hid_t getH5Type< localIndex >()  { return H5T_NATIVE_INT; }
template<> hid_t getH5Type< real32 >()    { return H5T_NATIVE_FLOAT; }
template<> hid_t getH5Type< real64 >()   { return H5T_NATIVE_DOUBLE; }

void createTestDataset( hid_t file_id, const std::string & name, hid_t hdf5_type,
                        const std::vector< hsize_t > & dims, const void * data )
{
  hid_t space_id = H5Screate_simple( static_cast< int >(dims.size()), dims.data(), nullptr );
  hid_t dset_id = H5Dcreate( file_id, name.c_str(), hdf5_type, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
  H5Dwrite( dset_id, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );
  H5Dclose( dset_id );
  H5Sclose( space_id );
}

template< typename T >
void createTestFileForAllRanks( const std::string & filename )
{
  hid_t file_id = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
  T data[16]{};
  for( int i = 0; i < 16; ++i )
    data[i] = static_cast< T >(i);

  std::string prefix;
  if constexpr (std::is_same_v< T, uint8_t >) prefix = "uint8_t";
  else if constexpr (std::is_same_v< T, localIndex >) prefix = "localIndex";
  else if constexpr (std::is_same_v< T, real32 >) prefix = "real32";
  else if constexpr (std::is_same_v< T, real64 >) prefix = "real64";

  createTestDataset( file_id, prefix + "_1D", getH5Type< T >(), {16}, data );
  createTestDataset( file_id, prefix + "_2D", getH5Type< T >(), {4, 4}, data );
  createTestDataset( file_id, prefix + "_3D", getH5Type< T >(), {2, 2, 4}, data );
  createTestDataset( file_id, prefix + "_4D", getH5Type< T >(), {2, 2, 2, 2}, data );
  H5Fclose( file_id );
}

// === Valid Conversion Test Fixture ===
template< typename Config >
class SerialHDF5ReaderTypedTest : public ::testing::Test
{
protected:
  using DatasetType = typename Config::dataset_type;
  using ArrayType = typename Config::array_type;
  const std::string filename = "typed_test.h5";

  void SetUp() override
  {
    createTestFileForAllRanks< DatasetType >( filename );
  }

  void TearDown() override
  {
    std::remove( filename.c_str());
  }

  template< typename T >
  void runReadTest( const std::string & datasetName, int expectedDims )
  {
    SerialHDF5Reader reader( filename );
    array1d< T > result = reader.read1DAs< T >( datasetName, expectedDims );
    ASSERT_EQ( result.size(), 16 );
    for( int i = 0; i < result.size(); ++i )
    {
      EXPECT_EQ( result[i], static_cast< T >(i) );
    }
  }
};

// === Supported Type Pairs ===
using MyTestTypes = ::testing::Types<
  ReadTestConfig< uint8_t, uint8_t >,
  ReadTestConfig< uint8_t, localIndex >,
  ReadTestConfig< uint8_t, real32 >,
  ReadTestConfig< uint8_t, real64 >,
  ReadTestConfig< localIndex, localIndex >,
  ReadTestConfig< localIndex, real32 >,
  ReadTestConfig< localIndex, real64 >,
  ReadTestConfig< real32, real32 >,
  ReadTestConfig< real32, real64 >,
  ReadTestConfig< real64, real64 >
  >;

TYPED_TEST_SUITE( SerialHDF5ReaderTypedTest, MyTestTypes );

TYPED_TEST( SerialHDF5ReaderTypedTest, CanRead1DTo4DFlattenedCorrectly )
{
  using OutType = typename TypeParam::array_type;
  std::string label = TypeParam::label();
  std::string marker = "to_";
  std::size_t pos = label.find( "to_" );
  std::string prefix;
  if( pos != std::string::npos )
  {
    prefix = label.substr( 0, pos-1 );
  }

  for( int rank = 1; rank <= 4; ++rank )
  {
    std::string name = prefix + "_" + std::to_string( rank ) + "D";

    this->template runReadTest< OutType >( name, rank );

  }
}

// === Invalid Dimension Tests ===
template< typename Config >
class SerialHDF5ReaderInvalidDimsTest : public ::testing::Test
{
protected:
  using DatasetType = typename Config::dataset_type;
  using ArrayType = typename Config::array_type;
  const std::string filename = "invalid_dims_test.h5";

  void SetUp() override
  {
    createTestFileForAllRanks< DatasetType >( filename );
  }

  void TearDown() override
  {
    std::remove( filename.c_str());
  }

  void expectInvalidDims( const std::string & datasetName, int wrongDims )
  {
    SerialHDF5Reader reader( filename );
    EXPECT_THROW( reader.read1DAs< ArrayType >( datasetName, wrongDims ), std::runtime_error );
  }
};

TYPED_TEST_SUITE( SerialHDF5ReaderInvalidDimsTest, MyTestTypes );

TYPED_TEST( SerialHDF5ReaderInvalidDimsTest, ThrowsOnWrongExpectedDims )
{
  std::string label = TypeParam::label();
  std::string marker = "to_";
  std::size_t pos = label.find( "to_" );
  std::string prefix;
  if( pos != std::string::npos )
  {
    prefix = label.substr( 0, pos-1 );
  }

  for( int actualDim = 1; actualDim <= 4; ++actualDim )
  {
    std::string dataset = prefix + "_" + std::to_string( actualDim ) + "D";
    for( int wrongDim = 1; wrongDim <= 4; ++wrongDim )
    {
      if( wrongDim == actualDim )
        continue;
      this->expectInvalidDims( dataset, wrongDim );
    }
  }
}

// === Unsupported Conversion Tests ===
using InvalidConversionTypes = ::testing::Types<
  ReadTestConfig< localIndex, uint8_t >,
  ReadTestConfig< real32, uint8_t >,
  ReadTestConfig< real32, localIndex >,
  ReadTestConfig< real64, uint8_t >,
  ReadTestConfig< real64, localIndex >,
  ReadTestConfig< real64, real32 >
  >;

template< typename Config >
class SerialHDF5ReaderUnsupportedConversionTest : public ::testing::Test
{
protected:
  using DatasetType = typename Config::dataset_type;
  using ArrayType = typename Config::array_type;
  const std::string filename = "unsupported_conversion.h5";

  void SetUp() override
  {
    createTestFileForAllRanks< DatasetType >( filename );
  }

  void TearDown() override
  {
    std::remove( filename.c_str());
  }

  void expectUnsupportedRead( const std::string & datasetName, int dims )
  {
    SerialHDF5Reader reader( filename );
    EXPECT_THROW( reader.read1DAs< ArrayType >( datasetName, dims ), std::runtime_error );
  }
};

TYPED_TEST_SUITE( SerialHDF5ReaderUnsupportedConversionTest, InvalidConversionTypes );

TYPED_TEST( SerialHDF5ReaderUnsupportedConversionTest, ThrowsOnUnsupportedConversion )
{
  std::string label = TypeParam::label();
  std::string marker = "to_";
  std::size_t pos = label.find( "to_" );
  std::string prefix;
  if( pos != std::string::npos )
  {
    prefix = label.substr( 0, pos-1 );
  }
  std::string dataset = prefix + "_4D";
  this->expectUnsupportedRead( dataset, 4 );
}

// === Non-Existent Dataset Test ===
template< typename Config >
class SerialHDF5ReaderMissingDatasetTest : public ::testing::Test
{
protected:
  using DatasetType = typename Config::dataset_type;
  using ArrayType = typename Config::array_type;
  const std::string filename = "missing_dataset_test.h5";

  void SetUp() override
  {
    createTestFileForAllRanks< DatasetType >( filename );
  }

  void TearDown() override
  {
    std::remove( filename.c_str());
  }

  void expectMissingDatasetRead( const std::string & datasetName, int dims )
  {
    SerialHDF5Reader reader( filename );
    EXPECT_THROW( reader.read1DAs< ArrayType >( datasetName, dims ), std::runtime_error );
  }
};

TYPED_TEST_SUITE( SerialHDF5ReaderMissingDatasetTest, MyTestTypes );

TYPED_TEST( SerialHDF5ReaderMissingDatasetTest, ThrowsOnMissingDataset )
{
  std::string missingDataset = "non_existent_dataset";
  this->expectMissingDatasetRead( missingDataset, 1 );
}

int main( int argc, char * * argv )
{
  geos::testing::HDF5UtilsTestScope scope( argc, argv );
  return RUN_ALL_TESTS();
}
