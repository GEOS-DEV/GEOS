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

/**
 * @file StringUtilities.hpp
 */

#ifndef GEOSX_CODINGUTILITIES_STRINGUTILITIES_OLD_HPP_
#define GEOSX_CODINGUTILITIES_STRINGUTILITIES_OLD_HPP_

#include <algorithm>
#include <cstring>
#include <cxxabi.h>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>


#include "common/DataTypes.hpp"
#include "LvArray/src/limits.hpp"

namespace geosx
{
namespace stringutilities
{

/// Overloaded function to check equality between strings and char arrays
/// Mainly used to avoid char*==char* mistakes
inline bool streq( string const & strA, string const & strB )
{ return strA == strB; }

inline bool streq( string const & strA, char const * const strB )
{ return strA == strB; }

inline bool streq( char const * const strA, string const & strB )
{ return strA == strB; }

inline bool streq( char const * const strA, char const * const strB )
{ return !strcmp( strA, strB ); }

/**
 * @brief Return a copy of the string in lower case.
 * @param input The input string which is not modified.
 * @return A new string instance.
 */
string toLower( string const & input );

/**
 * @brief Join strings or other printable objects with a delimiter.
 * @tparam S    type of delimiter, usually char, char const * or string
 * @tparam IT   type of iterator into the range of objects to join
 * @param delim delimiter used to glue together strings
 * @param first iterator to start of the range
 * @param last  iterator past-the-end of the range
 * @return a new string containing input strings concatenated with a delimiter
 */
template< typename IT, typename S = char >
string strjoin( IT first, IT last, S const & delim = S() )
{
  std::ostringstream oss;
  if( first != last )
  {
    oss << *first;
  }
  while( ++first != last )
  {
    oss << delim << *first;
  }
  return oss.str();
}

/// Subdivide string by delimiters
string_array Tokenize( string const & str, string const & delimiters );

/**
 * @brief Retuns a string containing a padded value
 * @param[in] value to be padded
 * @param[in] size size of the padding
 */
template< typename T >
string PadValue( T value, int size )
{
  std::stringstream paddedStringStream;
  paddedStringStream << std::setfill( '0' ) << std::setw( size ) << value;
  return paddedStringStream.str();
}
}
}

#endif /* GEOSX_CODINGUTILITIES_STRINGUTILITIES_OLD_HPP_ */
