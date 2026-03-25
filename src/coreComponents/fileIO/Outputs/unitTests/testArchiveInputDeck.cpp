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

  void SetUp() override
  {
    m_inputDir  = fs::temp_directory_path() / "geos_archiveInputDeck_test_input";
    m_outputDir = fs::temp_directory_path() / "geos_archiveInputDeck_test_output";
    fs::create_directories( m_inputDir );
    fs::create_directories( m_outputDir );
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

	/// Return the single timestamped subdir in inputFiles/ 
	/// or an empty path
	fs::path findArchiveDir()
	{
		fs::path const inputFilesDir = m_outputDir / "archive_inputFiles";
		if ( !fs::exists( inputFilesDir ) )
		{
			return {};
		}
		for ( auto const & entry : fs::directory_iterator( inputFilesDir ) )
		{
			if ( entry.is_directory() )
			{
				return entry.path();
			}
		}
		return {};
	}

	std::set< string > collectArchiveFiles( fs::path const & archiveDir )
	{
		std::set< string > files;

		for ( auto const & entry : fs::recursive_directory_iterator( archiveDir ) )
		{
			if ( entry.is_regular_file() )
			{
				files.insert( fs::relative( entry.path(), archiveDir ).string() );
			}
		}

		return files;
	}
};


TEST_F( ArchiveInputDeckTest, singleFile )
{
	writeXML( "base.xml", "<Problem>"
												"</Problem>" );
	string_array inputFiles { inputPath( "base.xml" ) };

	archiveInputDeck::archiveInputDeck( inputFiles, m_outputDir.string() );

	fs::path const archiveDir = findArchiveDir();
	EXPECT_FALSE( archiveDir.empty() ) << "Archive directory was not created";

	auto const files = collectArchiveFiles( archiveDir );
	EXPECT_NE( files.find( "base.xml" ), files.end() );
}

TEST_F( ArchiveInputDeckTest, fileWithSubdirInclude )
{
	writeXML( "subdir/child.xml", "<Problem>"
							 				      		"</Problem>" );
  writeXML( "base.xml", "<Problem>"
                        " <Included>"
                        "   <File name=\"subdir/child.xml\"/>"
                        " </Included>"
                        "</Problem>" );
	string_array inputFiles { inputPath( "base.xml" ) };

	archiveInputDeck::archiveInputDeck( inputFiles, m_outputDir.string() );

	fs::path const archiveDir = findArchiveDir();
	EXPECT_FALSE( archiveDir.empty() ) << "Archive directory was not created";

	auto const files = collectArchiveFiles( archiveDir );
	EXPECT_NE( files.find( "subdir/child.xml" ), files.end() );
}

TEST_F( ArchiveInputDeckTest, fileOneDirectoryAbove )
{
	writeXML( "other.xml", "<Problem>"
							 		   		 "</Problem>" );
  writeXML( "subdir/base.xml", "<Problem>"
                               " <Included>"
                               "   <File name=\"../other.xml\"/>"
                               " </Included>"
                               "</Problem>" );
	string_array inputFiles { inputPath( "subdir/base.xml" ) };

	archiveInputDeck::archiveInputDeck( inputFiles, m_outputDir.string() );

	fs::path const archiveDir = findArchiveDir();
	EXPECT_FALSE( archiveDir.empty() ) << "Archive directory was not created";

	auto const files = collectArchiveFiles( archiveDir );
	EXPECT_NE( files.find( "__other.xml" ), files.end() );
}

TEST_F( ArchiveInputDeckTest, fileTwoDirectoriesAbove )
{
	writeXML( "other.xml", "<Problem>"
							 		   		 "</Problem>" );
  writeXML( "subdir/subdir/base.xml", "<Problem>"
                                      " <Included>"
                                      "   <File name=\"../../other.xml\"/>"
                                      " </Included>"
                                      "</Problem>" );
	string_array inputFiles { inputPath( "subdir/subdir/base.xml" ) };

	archiveInputDeck::archiveInputDeck( inputFiles, m_outputDir.string() );

	fs::path const archiveDir = findArchiveDir();
	EXPECT_FALSE( archiveDir.empty() ) << "Archive directory was not created";

	auto const files = collectArchiveFiles( archiveDir );
	EXPECT_NE( files.find( "____other.xml" ), files.end() );
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  int const result = RUN_ALL_TESTS();

  return result;
}
