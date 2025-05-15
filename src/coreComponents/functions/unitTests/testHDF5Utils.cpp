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
    hid_t fileId = H5Fcreate(testFilename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    ASSERT_GE(fileId, 0);
    H5Fclose(fileId);
  }

  void TearDown() override
  {
    // Cleanup code after each test
    std::remove(testFilename.c_str());
  }

  std::string testFilename;
};

// Test for opening and closing an HDF5 file
TEST_F(SerialHDF5FileTest, OpenAndCloseFile)
{
  hid_t file_id = -1;
  {
    SerialHDF5File file(testFilename);

    // Check that the file ID is valid
    file_id = file.getFileId();
    ASSERT_GE( file_id, 0);
    EXPECT_TRUE(H5Iis_valid(file_id));
  }

  // After file goes out of scope, file is closed
  EXPECT_FALSE(H5Iis_valid(file_id));
}

// Test for move constructor
TEST_F(SerialHDF5FileTest, MoveConstructor)
{
  SerialHDF5File file1(testFilename);
  SerialHDF5File file2(std::move(file1));

  // Check that file1 is invalidated
  EXPECT_EQ(file1.getFileId(), -1);

  // Check that file2 has a valid file ID
  EXPECT_GE(file2.getFileId(), 0);
}

// Test for move assignment operator
TEST_F(SerialHDF5FileTest, MoveAssignment)
{
  SerialHDF5File file1(testFilename);
  SerialHDF5File file2 = std::move(file1);

  // Check that file1 is invalidated
  EXPECT_EQ(file1.getFileId(), -1);

  // Check that file2 has a valid file ID
  EXPECT_GE(file2.getFileId(), 0);
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
    hid_t fileId = H5Fcreate(testFilename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    ASSERT_GE(fileId, 0);

    hsize_t dims[1] = {10};
    hid_t dataspaceId = H5Screate_simple(1, dims, nullptr);
    ASSERT_GE(dataspaceId, 0);

    hid_t datasetId = H5Dcreate2(fileId, "test_dataset", H5T_NATIVE_INT, dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    ASSERT_GE(datasetId, 0);

    H5Dclose(datasetId);
    H5Sclose(dataspaceId);
    H5Fclose(fileId);
  }

  void TearDown() override
  {
    // Cleanup code after each test
    std::remove(testFilename.c_str());
  }

  std::string testFilename;
};

// Test for reading a 1D dataset
TEST_F(SerialHDF5ReaderTest, Read1DAsInt)
{
  SerialHDF5Reader reader(testFilename);

  // Read the dataset as an array of integers
  array1d<int> data = reader.read1DAs<int>("test_dataset", 1);

  // Check that the data size matches the expected size
  EXPECT_EQ(data.size(), 10);
}

// Test for reading a non-existent dataset
TEST_F(SerialHDF5ReaderTest, ReadNonExistentDataset)
{
  SerialHDF5Reader reader(testFilename);

  // Attempt to read a non-existent dataset and expect an exception
  EXPECT_THROW(reader.read1DAs<int>("non_existent_dataset", 1), std::runtime_error);
}


int main( int argc, char * * argv )
{
  geos::testing::HDF5UtilsTestScope scope( argc, argv );
  return RUN_ALL_TESTS();
}