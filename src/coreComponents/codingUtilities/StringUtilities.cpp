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
                      string const & delimiters,
                      bool const treatConsecutiveDelimAsOne )
{
  if( str.empty() )
  {
    return {};
  }

  RETURN_TYPE tokens;
  size_t lastPos = 0;
  size_t newPos;
  while( ( newPos = str.find_first_of( delimiters, lastPos ) ) != string::npos )
  {
    tokens.emplace_back( str.substr( lastPos, newPos - lastPos ) );
    lastPos = !treatConsecutiveDelimAsOne ? newPos + 1 : str.find_first_not_of( delimiters, newPos );
  }
  if( lastPos != string::npos )
  {
    tokens.emplace_back( str.substr( lastPos ) );
  }

  return tokens;
}
/**
 * String tokenizing by whitespace function
 **/
template< typename RETURN_TYPE >
RETURN_TYPE tokenizeBySpaces( string const & str )
{
  if( str.empty() )
  {
    return {};
  } else {
    const size_t strEnd=str.size();
    RETURN_TYPE tokens;
    size_t pos=0;

    while (pos<strEnd) {
      //search for a token begining
      while ( pos<strEnd && ( std::isspace( str[pos] ) ) )
      {
        ++pos;
      }

      if (pos<strEnd)
      {
        size_t lastPos=pos;

        //search for the token ending
        do
        {
          ++pos;
        } while ( pos<strEnd && !std::isspace( str[pos] ) );

        tokens.emplace_back( str.substr( lastPos, pos - lastPos ) );
      }
    }
    
    return tokens;
  }
}

template string_array tokenize< string_array >( string const & str,
                                                string const & delimiters,
                                                bool const treatConsecutiveDelimAsOne );

template std::vector< string > tokenize< std::vector< string > >( string const & str,
                                                                  string const & delimiters,
                                                                  bool const treatConsecutiveDelimAsOne );


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
