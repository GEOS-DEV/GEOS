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

#ifndef GEOSX_COMMON_PATH_HPP
#define GEOSX_COMMON_PATH_HPP


#include <unistd.h>
#include "common/Logger.hpp"
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <vector>

namespace geosx
{
/*!
 * @brief Class describing a file Path
 * @details Purpose of this class is to
 * be used as a type to specify file path within the XML input files,
 * through the operator >>
 */
class Path : public std::string
{
public:

  ~Path()
  {}

  using std::string::string;

  /*!
   * @brief Get the path prefix of the file
   * @details The path prefix is usually a folder path
   * in which the XML file is located
   * @return the path prefix
   */
  static std::string & pathPrefix()
  {
    static std::string m_pathPrefix;
    return m_pathPrefix;
  }
};

/*!
 * @brief Gets the absolute path of a file
 * @param[in] path the relative path to the file
 * @param[out] absolutePath the absolute path to the file
 */
inline void getAbsolutePath( std::string const & path, std::string & absolutePath )
{
  char absFilePath[ PATH_MAX + 1 ];
  if( realpath( path.data(), absFilePath ) )
  {
    absolutePath = absFilePath;
  }
  else
  {
    char const * ret = getcwd( absFilePath, PATH_MAX + 1 );
    if( ret != nullptr )
      GEOSX_ERROR( "Could not get the absolute path for " << path << " from " << absFilePath );
    else
      GEOSX_ERROR( "Could not get the absolute path for " << path );
  }
}

/*!
 * @brief Tells wether the path is absolute of not
 * @param[in] path the input path
 * @retval true if the path is absolute
 * @retval false if the path is relative
 */
inline bool isAbsolutePath( const std::string & path )
{
  return path[ 0 ] == '/';
}

/*!
 * @brief Operator use with the class Path while parsing the XML file
 * @param[in,out] is the input stream
 * @param[in,out] p the path that will be set to an absolute path relative to the xml file
 * @return the input stream
 */
inline std::istream & operator>>( std::istream & is, Path & p )
{
  is >> static_cast< std::string & >(p);
  if( !isAbsolutePath( p ) && !p.pathPrefix().empty())
  {
    getAbsolutePath( p.pathPrefix() + '/' + p, p );
  }
  return is;
}

/*!
 * @brief Split the path in two parts : directory name and file name
 * @param[in] path the input path
 * @param[out] dirName name of the directory
 * @param[out] baseName the name of the file
 */
inline void splitPath( std::string const & path, std::string & dirName, std::string & baseName )
{
  size_t pos = path.find_last_of( '/' );

  if( pos == std::string::npos )
  {
    dirName = ".";
    baseName = path;
  }
  else if( pos == 0 )
  {
    dirName = "/";
    baseName = path.substr( 1 );
    return;
  }
  else
  {
    dirName = path.substr( 0, pos );
    baseName = path.substr( pos + 1 );
  }
}

/*!
 * @brief List all the files of one directory
 * @details Taken from http://www.martinbroadhurst.com/list-the-files-in-a-directory-in-c.html
 * @param[in] path path to the directory
 * @param[out] files vector containing allt the file path
 */
inline void readDirectory( std::string const & path, std::vector< std::string > & files )
{
  DIR * dirp = opendir( path.c_str() );
  struct dirent * dp;
  while( (dp = readdir( dirp )) != nullptr )
  {
    files.push_back( dp->d_name );
  }
  closedir( dirp );
}

} /* end namespace geosx */


#endif
