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

/**
 * @file ArchiveInputDeck.cpp
 */

#include "ArchiveInputDeck.hpp"

#include "common/Path.hpp"
#include "dataRepository/xmlWrapper.hpp"

#include <chrono>
#include <filesystem>


namespace geos
{

using namespace dataRepository;

namespace archiveInputDeck
{

namespace
{

string makeTimestamp()
{
  auto const now = std::chrono::system_clock::now();
  auto const time_t_now = std::chrono::system_clock::to_time_t( now );
  std::ostringstream timestampStream;
  timestampStream << std::put_time( std::localtime( &time_t_now ), "%Y%m%d_%H%M%S" );
  return timestampStream.str();
}

std::set< string > collectAbsFilePaths( string_array const & fileNames )
{
  std::set< string > collection;
  for ( string const & fileName : fileNames )
  {
    xmlWrapper::collectIncludedRecursive( fileName, collection );
  }
  return collection;
}

/// @brief Prefixes a file path string if it is located "behind" the 
///        specified directory
/// @param absFilePath The absolute path to the file 
/// @param absDirPath The absolute path to the directory
/// @return A relative path of the file prefixed with "__" for every "../" 
///         from the directory location
///
/// Example:
/// @code
///   std::string foo = prefixBackwardPath( "/usr/foo/file.txt", "/usr/bar/buzz" )
///   assert( foo == "____foo/file.txt" )
/// @endcode
string prefixBackwardPath( string const & absFilePath, string const & absDirPath )
{
    string relPath = std::filesystem::relative( std::filesystem::path( absFilePath ), 
                                                std::filesystem::path( absDirPath ) );

    string prefix;
    while( relPath.size() >= 3 && relPath.substr( 0, 3 ) == "../" )
    {
      prefix += "__";
      relPath = relPath.substr( 3 );
    }

    return prefix + relPath;
}

}

void archiveInputDeck( string_array const & inputFileNames,
                       string const & outputDirectory )
{
  if ( inputFileNames.empty() )
  {
    return;
  }

  string const timestamp = makeTimestamp();
  string const archiveDir = joinPath( outputDirectory, "inputFiles", timestamp );
  makeDirsForPath( archiveDir + "/" );

  string const baseDir = splitPath( getAbsolutePath(inputFileNames[0]) ).first;
  std::set< string > absFilePaths = collectAbsFilePaths( inputFileNames );

  for ( string const & absFilePath : absFilePaths )
  {
    string const destPath = joinPath( archiveDir, prefixBackwardPath( absFilePath, baseDir ) );
    makeDirsForPath( splitPath( destPath ).first + "/" );

    std::error_code ec;
    bool copied = std::filesystem::copy_file( absFilePath, 
                                              destPath, 
                                              std::filesystem::copy_options::overwrite_existing,
                                              ec );
    GEOS_LOG_IF( !copied,
                 GEOS_FMT( "Failed to copy archive file '{}' into '{}': {}",
                           absFilePath, destPath, ec.message() ) );
  }

}

} /* namespace archiveInputDeck */

} /* namespace geos */
