/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "Path.hpp"
#include "Logger.hpp"

#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <vector>

namespace geosx
{

std::string getAbsolutePath( std::string const & path )
{
  char absFilePath[ PATH_MAX + 1 ];
  if( realpath( path.data(), absFilePath ) )
  {
    return absFilePath;
  }
  else
  {
    char const * ret = getcwd( absFilePath, PATH_MAX + 1 );
    if( ret != nullptr )
    {
      GEOSX_THROW( "Could not get the absolute path for " << path << " from " << absFilePath, InputError );
    }
    GEOSX_THROW( "Could not get the absolute path for " << path, InputError );
  }
}

std::istream & operator>>( std::istream & is, Path & p )
{
  std::string & s = static_cast< std::string & >( p );
  is >> s;
  if( !isAbsolutePath( s ) && !Path::pathPrefix().empty() )
  {
    s = getAbsolutePath( Path::pathPrefix() + '/' + p );
  }
  return is;
}

std::pair< std::string, std::string > splitPath( std::string const & path )
{
  std::pair< std::string, std::string > parts;
  size_t pos = path.find_last_of( '/' );

  if( pos == std::string::npos )
  {
    parts.first = ".";
    parts.second = path;
  }
  else if( pos == 0 )
  {
    parts.first = "/";
    parts.second = path.substr( 1 );
  }
  else if( pos == path.size() - 1 )
  {
    // trim trailing slash
    parts = splitPath( path.substr( 0, pos ) );
  }
  else
  {
    parts.first = path.substr( 0, pos );
    parts.second = path.substr( pos + 1 );
  }
  return parts;
}

std::vector< std::string > readDirectory( std::string const & path )
{
  // Taken from http://www.martinbroadhurst.com/list-the-files-in-a-directory-in-c.html
  std::vector< std::string > files;
  DIR * dirp = opendir( path.c_str() );
  struct dirent * dp;
  while( (dp = readdir( dirp )) != nullptr )
  {
    files.emplace_back( dp->d_name );
  }
  closedir( dirp );
  return files;
}

void makeDirectory( std::string const & path )
{
  constexpr mode_t mode = 0770; // user and group rwx permissions
  int const err = mkdir( path.c_str(), mode );
  GEOSX_THROW_IF( err && ( errno != EEXIST ), "Failed to create directory: " << path, std::runtime_error );
}

void makeDirsForPath( std::string const & path )
{
  std::string::size_type pos = 0;
  do
  {
    pos = path.find( '/', pos + 1 );
    makeDirectory( path.substr( 0, pos ) );
  }
  while( pos != std::string::npos );
}

} /* end namespace geosx */
