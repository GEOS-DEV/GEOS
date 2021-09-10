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

namespace geosx
{
namespace stringutilities
{

string toLower( string const & input )
{
  string output;
  output.resize( input.size() );
  auto const toLowerCase = []( unsigned char c ) { return std::tolower( c ); };
  std::transform( input.cbegin(), input.cend(), output.begin(), toLowerCase );
  return output;
}

/**
 * String tokenizing function
 **/
string_array tokenize( const string & str, const string & delimiters )
{
  string_array tokens;

  if( str.length() == 0 )
  {
    tokens.emplace_back( str );
  }
  else
  {

    bool usesNonWhitespaceDelimiters = false;
    string::size_type i =0;
    while( delimiters[i] && !usesNonWhitespaceDelimiters )
    {
      usesNonWhitespaceDelimiters |= !isspace( int(delimiters[i]) );
      ++i;
    }

    if( usesNonWhitespaceDelimiters )
    {
      // do not skip multiple adjacent delimiters - indicates empty strings
      size_t lastPos = 0;

      size_t newPos = lastPos;
      while( (newPos=str.find_first_of( delimiters, lastPos )) != string::npos )
      {
        tokens.emplace_back( str.substr( lastPos, newPos-lastPos ));
        lastPos = newPos + 1;
      }
      tokens.emplace_back( str.substr( lastPos, str.length()-lastPos ));
    }
    else
    {
      // whitespace delimiters
      // skip multiple adjacent delimiters
      size_t lastPos = str.find_first_not_of( delimiters, 0 );
      lastPos = (lastPos == string::npos) ? 0 : lastPos;

      size_t newPos = lastPos;
      while( (newPos=str.find_first_of( delimiters, lastPos )) != string::npos )
      {
        tokens.emplace_back( str.substr( lastPos, newPos-lastPos ));
        lastPos = str.find_first_not_of( delimiters, newPos );
      }
      if( lastPos!= string::npos )
        tokens.emplace_back( str.substr( lastPos, str.length()-lastPos ));

    }
  }
  return tokens;
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
