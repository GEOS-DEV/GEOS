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

#include "fileIO/Outputs/ArchiveInputDeck.hpp"

#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>

using namespace geos;

namespace fs = std::filesystem;

class ArchiveInputDeckTest : public ::testing::Test
{
protected:
  fs::path m_inputDir;
  fs::path m_outputDir;
  string_array m_tagOrder;

  void SetUp() override
  {
    m_inputDir  = fs::temp_directory_path() / "geos_archiveInputDeck_test_input";
    m_outputDir = fs::temp_directory_path() / "geos_archiveInputDeck_test_output";
    fs::create_directories( m_inputDir );
    fs::create_directories( m_outputDir );
    m_tagOrder = string_array{ "Mesh", "ElementRegions", "FieldSpecifications", "Outputs" };
  }

  void TearDown() override
  {
    fs::remove_all( m_inputDir );
    fs::remove_all( m_outputDir );
  }

  string inputPath( string const & filename )
  {
    return ( m_inputDir / filename ).string();
  }

  void writeXML( string const & filename, string const & content )
  {
    fs::path const fullPath = m_inputDir / filename;
    fs::create_directories( fullPath.parent_path() );
    std::ofstream f( fullPath );
    f << content;
  }

  /// Return the single timestamped subdir in archive_inputFiles/ or an empty path
  fs::path findArchiveDir()
  {
    fs::path const inputFilesDir = m_outputDir / "archive_inputFiles";
    if( !fs::exists( inputFilesDir ) )
    {
      return {};
    }
    for( auto const & entry : fs::directory_iterator( inputFilesDir ) )
    {
      if( entry.is_directory() )
      {
        return entry.path();
      }
    }
    return {};
  }

  std::set< string > collectArchiveFiles( fs::path const & archiveDir )
  {
    std::set< string > files;

    for( auto const & entry : fs::recursive_directory_iterator( archiveDir ) )
    {
      if( entry.is_regular_file() )
      {
        files.insert( fs::relative( entry.path(), archiveDir ).string() );
      }
    }

    return files;
  }
};

TEST_F( ArchiveInputDeckTest, levelZeroDoesNotArchive )
{
  writeXML( "base.xml", "<Problem>"
                        "</Problem>" );
  string_array const inputFiles{ inputPath( "base.xml" ) };

  archiveInputDeck::archiveInputDeck( inputFiles, m_outputDir.string(), m_tagOrder, 0 );

  EXPECT_TRUE( findArchiveDir().empty() );
}

TEST_F( ArchiveInputDeckTest, emptyInputFileDoesNotArchive )
{
  string_array const inputFiles;

  archiveInputDeck::archiveInputDeck( inputFiles, m_outputDir.string(), m_tagOrder, 1 );

  EXPECT_TRUE( findArchiveDir().empty() );
}

TEST_F( ArchiveInputDeckTest, levelOnecreatesFlatXmlWithoutSchema )
{
  writeXML( "base.xml", "<Problem>"
                        "  <Mesh/>"
                        "</Problem>" );
  string_array inputFiles { inputPath( "base.xml" ) };

  archiveInputDeck::archiveInputDeck( inputFiles, m_outputDir.string(), m_tagOrder, 1 );

  fs::path const archiveDir = findArchiveDir();
  ASSERT_FALSE( archiveDir.empty() ) << "Archive directory was not created";
  EXPECT_TRUE( fs::exists( archiveDir / "input.xml" ) );
  EXPECT_FALSE( fs::exists( archiveDir / "schema.xsd" ) ) << "Schema should not be archived at level 1";
}

TEST_F( ArchiveInputDeckTest, levelTwoCopiesSchema )
{
  writeXML( "base.xml", "<Problem>"
                        "</Problem>" );
  string_array inputFiles { inputPath( "base.xml" ) };

  archiveInputDeck::archiveInputDeck( inputFiles, m_outputDir.string(), m_tagOrder, 2 );

  fs::path const archiveDir = findArchiveDir();
  ASSERT_FALSE( archiveDir.empty() ) << "Archive directory was not created";

  EXPECT_TRUE( fs::exists( archiveDir / "input.xml" ) );
  EXPECT_TRUE( fs::exists( archiveDir / "schema.xsd" ) );
}

TEST_F( ArchiveInputDeckTest, archiveContainsContentOfIncludes )
{
  writeXML( "base.xml", "<Problem>"
                        "  <Included>"
                        "    <File name=\"subdir/include.xml\"/>"
                        "  </Included>"
                        "</Problem>" );
  writeXML( "subdir/include.xml", "<Problem>"
                                  "  <FieldSpecifications>"
                                  "    <FieldSpecification"
                                  "      name=\"test_name\""
                                  "    />"
                                  "  </FieldSpecifications>"
                                  "</Problem>" );
  string_array inputFiles { inputPath( "base.xml" ) };

  archiveInputDeck::archiveInputDeck( inputFiles, m_outputDir.string(), m_tagOrder, 1 );

  fs::path const archiveDir = findArchiveDir();
  ASSERT_FALSE( archiveDir.empty() ) << "Archive directory was not created";

  ASSERT_TRUE( fs::exists( archiveDir / "input.xml" ) );

  std::ifstream ifs( archiveDir / "input.xml" );
  std::string content( (std::istreambuf_iterator< char >( ifs ) ),
                       (std::istreambuf_iterator< char >()    ) );

  EXPECT_NE( content.find( "test_name" ), std::string::npos )
    << "input.xml does not contain content from include.xml";
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  int const result = RUN_ALL_TESTS();

  return result;
}
