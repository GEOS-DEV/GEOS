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

#ifndef GEOSX_COMMON_PATH_HPP
#define GEOSX_COMMON_PATH_HPP

// System includes
#include <string>
#include <sstream>
#include <vector>

namespace geosx
{
/*!
 * @brief Class describing a file Path
 *
 * Purpose of this class is to be used as a type to specify file path
 * within the XML input files, through the operator>>.
 */
class Path : public std::string
{
public:

  /// Default constructor.
  Path():
    std::string()
  {}

  /**
   * @brief Copy constructor, creates a copy of @p src.
   * @param src the Path to copy.
   */
  Path( Path const & src ):
    std::string( src )
  {}

  /// Destructor.
  ~Path()
  {}

  /**
   * @brief Copy Constructor
   * @param rhs Reference to the Path that will be copied.
   * @return *this
   */
  Path & operator=( Path const & rhs )
  {
    std::string::operator=( rhs );
    return *this;
  }

  using std::string::string;

  /*!
   * @brief Get the path prefix of the file
   * @details The path prefix is usually a folder path in which the XML file is located
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
 * @return the absolute path to the file resolved from current working directory
 * @exception InputError if absolute path could not be resolved
 */
std::string getAbsolutePath( std::string const & path );

/*!
 * @brief Tells whether the path is absolute of not
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
std::istream & operator>>( std::istream & is, Path & p );

/**
 * @brief Remove the trailing slash is if present
 * @param path the input path
 * @return the trimmed path
 */
inline std::string trimPath( std::string const & path )
{
  return path.back() == '/' ? path.substr( 0, path.size() - 1 ) : path;
}

/*!
 * @brief Split the path in two parts: directory name and file name
 * @param[in] path the input path
 * @return a pair of strings, containing names of the directory and the file
 */
std::pair< std::string, std::string > splitPath( std::string const & path );

/*!
 * @brief Join two parts of a path
 * @param baseName base name (e.g. directory)
 * @param name name to add (e.g. relative path to file or directory)
 * @return the combined path
 */
inline std::string joinPath( std::string const & baseName, std::string const & name )
{
  return trimPath( baseName ) + '/' + trimPath( name );
}

/*!
 * @brief Join parts of a path
 * @tparam ARGS types of arguments
 * @param args parts of the path
 * @return the combined path
 */
template< typename ... ARGS >
inline std::string joinPath( ARGS const & ... args )
{
  size_t constexpr numParts = sizeof...(args);
  static_assert( numParts > 0, "Must provide arguments" );
  std::string parts[numParts] { trimPath( args ) ... };
  std::ostringstream oss;
  oss << parts[0];
  for( size_t i = 1; i < numParts; ++i )
  {
    oss << '/' << parts[i];
  }
  return oss.str();
}

/*!
 * @brief List all the files of one directory
 * @param[in] path path to the directory
 * @return vector containing all the file paths
 */
std::vector< std::string > readDirectory( std::string const & path );

/*!
 * @brief Create a directory @p path, where parent directories must already exist.
 * @param path The path to create.
 */
void makeDirectory( std::string const & path );

/*!
 * @brief Make directories for @p path.
 * @param path The path to create.
 *
 * This function operates similarly to 'mkdir -p'.
 * Everything in @p path is intended to be a directory.
 * If a directory in the path already exists nothing is done.
 * If a directory doesn't exist it is created.
 */
void makeDirsForPath( std::string const & path );

} /* end namespace geosx */


#endif
