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
 * @file StringUtilities.cpp
 */

#include "StringUtilities.hpp"

#include <algorithm>

namespace geosx
{
namespace stringutilities
{

string toLower( string const & input )
{
  string output;
  output.resize( input.size() );
  auto const toLowerCase = []( unsigned char c )
  { return std::tolower( c ); };
  std::transform( input.cbegin(), input.cend(), output.begin(), toLowerCase );
  return output;
}

/**
 * String tokenizing function
 **/
template< typename RETURN_TYPE >
RETURN_TYPE tokenize( string const & str,
                      string const & delimiters )
{
  if( str.empty() )
  {
    return {};
  }

  auto const isNonSpace = []( char const c ){ return !isspace( c ); };
  bool const usesNonWhitespaceDelimiters = std::any_of( delimiters.begin(), delimiters.end(), isNonSpace );

  // When only whitespace delimiters, skip multiple adjacent delimiters; otherwise don't and keep empty tokens
  RETURN_TYPE tokens;
//  size_t lastPos = usesNonWhitespaceDelimiters ? 0 : str.find_first_not_of( delimiters, 0 );
  size_t lastPos = str.find_first_not_of( delimiters, 0 );
  size_t newPos;
  while( ( newPos = str.find_first_of( delimiters, lastPos ) ) != string::npos )
  {
    tokens.emplace_back( str.substr( lastPos, newPos - lastPos ) );
    lastPos = usesNonWhitespaceDelimiters ? newPos + 1 : str.find_first_not_of( delimiters, newPos );
  }
  if( lastPos != string::npos )
  {
    tokens.emplace_back( str.substr( lastPos ) );
  }

  return tokens;
}

template string_array tokenize< string_array >( string const & str,
                                                string const & delimiters );

template std::vector< string > tokenize< std::vector< string > >( string const & str,
                                                                  string const & delimiters );


bool areAnyCharsInString( string const & str, string const & chars )
{
  bool rval = false;

  for( size_t a=0; a<chars.size(); ++a )
  {
    if( str.find( chars[a] )!=string::npos )
    {
      rval = true;
    }
  }
  return rval;
}



string trim( string const & str,
             string const & charsToRemove )
{
  std::size_t const first = str.find_first_not_of( charsToRemove );
  if( first != string::npos )
  {
    std::size_t const last = str.find_last_not_of( charsToRemove );
    return str.substr( first, ( last - first + 1 ) );
  }
  return {};
}


string removeStringAndFollowingContent( string const & str,
                                        string const & strToRemove )
{
  string newStr = str;

  // check if the line contains the string to remove
  std::size_t const pos = newStr.find( strToRemove );

  if( pos != string::npos )
  {
    // remove the character and everything afterwards
    newStr = newStr.substr( 0, pos );
  }
  return newStr;
}

}
}
