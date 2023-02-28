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
template< template< class ... > class LIST_T = std::vector >
LIST_T< string > tokenize( string const & str,
                           string const & delimiters,
                           bool const treatConsecutiveDelimAsOne,
                           bool const preTrimStr )
{
  if( str.empty())
  {
    return {};
  }

  LIST_T< string > tokens;
  size_t tokenBegin, tokenEnd, strEnd;

  if( preTrimStr )
  {
    tokenBegin = str.find_first_not_of( delimiters );
    strEnd = str.find_last_not_of( delimiters ) + 1;
  }
  else
  {
    tokenBegin = 0;
    strEnd = str.size();
  }

  while( ( ( tokenEnd = str.find_first_of( delimiters, tokenBegin ) ) < strEnd ) && tokenBegin < strEnd )
  {
    tokens.emplace_back( str.substr( tokenBegin, tokenEnd - tokenBegin ) );
    tokenBegin = !treatConsecutiveDelimAsOne ? tokenEnd + 1 : str.find_first_not_of( delimiters, tokenEnd );
  }

  if( tokenBegin < strEnd )
  {
    tokens.emplace_back( str.substr( tokenBegin, strEnd-tokenBegin ));
  }
  else if( !preTrimStr )
  {
    if( str.find_first_of( delimiters, strEnd - 1 ) != string::npos )
    {
      tokens.emplace_back( "" );
    }
  }

  return tokens;
}
/**
 * String tokenizing by whitespace function
 **/
template< template< class ... > class LIST_T = std::vector >
LIST_T< string > tokenizeBySpaces( string const & str )
{
  return tokenize< LIST_T >( str, " \f\n\r\t\v", true, true );
}

//template specialization for tokenize and tokenizeBySpaces
template std::vector< string > tokenize< std::vector >( string const & str,
                                                        string const & delimiters,
                                                        bool const treatConsecutiveDelimAsOne,
                                                        bool const preTrimStr );
template std::list< string > tokenize< std::list >( string const & str,
                                                    string const & delimiters,
                                                    bool const treatConsecutiveDelimAsOne,
                                                    bool const preTrimStr );
template array1d< string > tokenize< array1d >( string const & str,
                                                string const & delimiters,
                                                bool const treatConsecutiveDelimAsOne,
                                                bool const preTrimStr );
template std::vector< string > tokenizeBySpaces< std::vector >( string const & str );
template std::list< string > tokenizeBySpaces< std::list >( string const & str );
template array1d< string > tokenizeBySpaces< array1d >( string const & str );

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
