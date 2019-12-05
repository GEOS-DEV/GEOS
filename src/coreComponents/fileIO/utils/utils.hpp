/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_FILEIO_UTILS_UTILS_HPP
#define GEOSX_FILEIO_UTILS_UTILS_HPP

#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <vector>
#include <regex>

namespace geosx
{

/* Taken from http://www.martinbroadhurst.com/list-the-files-in-a-directory-in-c.html */
inline void readDirectory( std::string const & name, std::vector< std::string >& v)
{
  DIR * dirp = opendir( name.c_str() );
  struct dirent * dp;
  while( (dp = readdir( dirp )) != nullptr )
  {
    v.push_back( dp->d_name );
  }
  closedir( dirp );
}


inline void getAbsolutePath( std::string const & path, std::string & absolute_path )
{
  char abs_file_path[ PATH_MAX + 1 ];
  if( realpath( path.data(), abs_file_path ) )
  {
    absolute_path = abs_file_path;
  }
  else
  {
    char const * ret = getcwd( abs_file_path, PATH_MAX + 1 );
    if ( ret != nullptr )
      GEOS_ERROR( "Could not get the absolute path for " << path << " from " << abs_file_path );
    else
      GEOS_ERROR( "Could not get the absolute path for " << path );
  }
}


inline void splitPath( std::string const & path, std::string & dirname, std::string & basename )
{
  size_t pos = path.find_last_of( '/' );
  
  if( pos == string::npos )
  {
    dirname = ".";
    basename = path;
  }
  else if( pos == 0 )
  {
    dirname = "/";
    basename = path.substr( 1 );
    return;
  }
  else
  {
    dirname = path.substr( 0, pos );
    basename = path.substr( pos + 1 );
  }
}

inline bool isAbsolutePath( const std::string & path )
{
  return path[ 0 ] == '/';
}

template< typename REGEX >
inline bool regexMatch( const std::string & str, REGEX regex )
{
  return std::regex_match(str, regex);
}

} /* end namespace geosx */

#endif /* GEOSX_FILEIO_UTILS_UTILS_HPP */
