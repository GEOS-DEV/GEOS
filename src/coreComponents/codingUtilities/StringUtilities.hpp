/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file StringUtilities.hpp
 */

#ifndef GEOSX_CODINGUTILITIES_STRINGUTILITIES_HPP_
#define GEOSX_CODINGUTILITIES_STRINGUTILITIES_HPP_

#include <iomanip>
#include <sstream>

#include "common/DataTypes.hpp"

namespace geosx
{
namespace stringutilities
{

/**
 * @brief Return a copy of the string in lower case.
 * @param input The input string which is not modified.
 * @return A new string instance.
 */
string toLower( string const & input );

/**
 * @brief Join strings or other printable objects with a delimiter.
 * @tparam IT   type of iterator into the range of objects to join
 * @tparam S    type of delimiter, usually char, char const * or string
 * @param first iterator to start of the range
 * @param last  iterator past-the-end of the range
 * @param delim delimiter used to glue together strings
 * @return a string containing input values concatenated with a delimiter
 */
template< typename IT, typename S = char >
string join( IT first, IT last, S const & delim = S() )
{
  if( first == last )
  {
    return {};
  }
  std::ostringstream oss;
  oss << *first;
  while( ++first != last )
  {
    oss << delim << *first;
  }
  return oss.str();
}

/**
 * @brief Join strings or other printable objects with a delimiter.
 * @tparam CONTAINER type of container to join
 * @tparam S    type of delimiter, usually char, char const * or string
 * @param container container to join
 * @param delim delimiter used to glue together strings
 * @return a string containing input values concatenated with a delimiter
 */
template< typename CONTAINER, typename S = char >
string join( CONTAINER const & cont, S const & delim = S() )
{
  return join( std::begin( cont ), std::end( cont ), delim );
}

/**
 * @brief Concatenate variadic arguments into a string with a delimiter.
 * @tparam S type of delimiter (printable to std::ostringstream)
 * @tparam T type of first argument (printable to std::ostringstream)
 * @tparam Ts types of remaining arguments (printable to std::ostringstream)
 * @param delim delimiter
 * @param v first value
 * @param vs remaining values
 * @return string containing concatenated printed arguments
 */
template< typename S = char, typename T, typename ... Ts >
string concat( S const & delim, T const & v, Ts const & ... vs )
{
  std::ostringstream oss;
  oss << v;
  // Use array initializer and comma trick to get "fold expression" pre C++-17
  using expander = int[];
  (void) expander{ 0, ( void ( oss << delim << vs ), 0) ... };
  return oss.str();
}

/// Subdivide string by delimiters
template< typename RETURN_TYPE = string_array >
RETURN_TYPE tokenize( string const & str,
                      string const & delimiters,
                      bool const treatConsecutiveDelimAsOne = true );

/**
 * @brief Trim the string
 * @param[in] str the string to trim
 * @param[in] charsToRemove the list of characters to remove
 * @return the trimmed string
 */
string trim( string const & str,
             string const & charsToRemove );

/**
 * @brief Search for a string in the line, and return the line truncated before the string
 * @param[in] str the line to truncate
 * @param[in] strToRemove the string to search for in the line
 * @return the new (truncated) string
 */
string removeStringAndFollowingContent( string const & str,
                                        string const & strToRemove );

/**
 * @brief Take a string, and return a array1d with the cast values
 * @tparam T the type to which the string will be cast
 * @param[in] str the string to turn into an array1d
 * @return the array1d that stores the cast values
 */
template< typename T >
array1d< T > fromStringToArray( string const & str )
{
  array1d< T > v;
  T sub;

  std::istringstream iss( str );
  while( iss >> sub )
  {
    v.emplace_back( sub );
  }
  return v;
}

} // namespace stringutilities
} // namespace geosx

#endif /* GEOSX_CODINGUTILITIES_STRINGUTILITIES_HPP_ */
