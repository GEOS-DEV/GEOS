/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file StringUtilities.cpp
 */

#include "StringUtilities.hpp"
#include "common/logger/Logger.hpp"
#include "limits.h"

#include <algorithm>
#include <cstdlib>
#include <cmath>

namespace geos
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

string_view trim( string_view str,
                  string_view charsToRemove )
{
  std::size_t const first = str.find_first_not_of( charsToRemove );
  if( first != string::npos )
  {
    std::size_t const last = str.find_last_not_of( charsToRemove );
    return str.substr( first, ( last - first + 1 ) );
  }
  return {};
}
string_view trimSpaces( string_view str )
{
  return trim( str, " \f\n\r\t\v" );
}


std::string_view ltrimSpaces( std::string_view s )
{
  std::size_t const first = s.find_first_not_of( " \f\n\r\t\v" );
  if( first != std::string::npos )
  {
    return s.substr( first, ( s.size() - first ) );
  }
  return {};
}


string removeStringAndFollowingContent( string_view const str,
                                        string_view const strToRemove )
{
  string_view newStr = str;

  // check if the line contains the string to remove
  std::size_t const pos = newStr.find( strToRemove );

  if( pos != string::npos )
  {
    // remove the character and everything afterwards
    newStr = newStr.substr( 0, pos );
  }
  return string( newStr );
}

// Add comma separators for thousands
template< typename T >
string addCommaSeparators( T const & num )
{
  static_assert( std::is_integral< T >::value, "addCommaSeparators only supports integral types" );

  string const numStr = std::to_string( num );
  string result;

  for( std::size_t i = 0; i < numStr.size(); ++i )
  {
    result += numStr[i];
    if((numStr.size() - i - 1) % 3 == 0 && i != numStr.size() - 1 )
    {
      result += ",";
    }
  }
  return result;
}
template string addCommaSeparators( int const & num );
template string addCommaSeparators( long int const & num );
template string addCommaSeparators( long long int const & num );

// put definition here so we can control the allowable values of T and
// modication of this function triggers a whole code recompile...which
// should be avoided.
template< typename T >
string toMetricPrefixString( T const & value )
{
  if( std::fpclassify( value ) == FP_ZERO )
  {
    return " 0.0  ";
  }

  // These are the metric prefixes corrosponding to kilo, mega, giga...etc.
  char const prefixes[12] = { 'f', 'p', 'n', 'u', 'm', ' ', 'K', 'M', 'G', 'T', 'P', 'E'};
  string rval;

  int const power = floor( log10( std::abs( (double)value ) ) );
  int const a = floor( power / 3.0 );

  real64 const scaledValue = value * pow( 10.0, -a * 3 );

  // format the output of the value to 3 significant digits and append the
  // metric prefix.
  int const p = 2-std::abs( power - a * 3 );
  char temp[10];
  snprintf( temp, 8, "%5.*f %c", p, scaledValue, prefixes[a+5] );
  rval = temp;

  GEOS_ERROR_IF( rval.empty(),
                 GEOS_FMT( "The value of {} was not able to be converted with a metric prefix", value ) );


  return rval;
}
template string toMetricPrefixString( int const & );
template string toMetricPrefixString( long int const & );
template string toMetricPrefixString( long long int const & );
template string toMetricPrefixString( unsigned long int const & );
template string toMetricPrefixString( unsigned long long int const & );
template string toMetricPrefixString( float const & );
template string toMetricPrefixString( double const & );

template< typename STRING_T >
stdVector< STRING_T > divideLines( size_t & linesWidth, string_view value )
{
  size_t current = 0;
  size_t end = value.find( '\n' );

  stdVector< STRING_T > lines;
  linesWidth = 0;

  // Process each line until no more newlines are found
  while( end != STRING_T::npos )
  {
    lines.push_back( STRING_T( value.substr( current, end - current ) ) );
    current = end + 1;
    end = value.find( '\n', current );
    linesWidth = std::max( linesWidth, lines.back().size() );
  }
  // Add the last part
  if( current <= value.size())
  {
    lines.push_back( STRING_T( value.substr( current )  ) );
    linesWidth = std::max( linesWidth, lines.back().size() );
  }

  return lines;
}

template< typename STRING_T >
std::vector< STRING_T > divideLines( string_view value )
{
  size_t current = 0;
  size_t end = value.find( '\n' );

  std::vector< STRING_T > lines;

  // Process each line until no more newlines are found
  while( end != STRING_T::npos )
  {
    lines.push_back( STRING_T( value.substr( current, end - current ) ) );
    current = end + 1;
    end = value.find( '\n', current );
  }
  // Add the last part
  if( current <= value.size())
    lines.push_back( STRING_T( value.substr( current )  ) );

  return lines;
}
template std::vector< string > divideLines( size_t &, string_view );
template std::vector< string_view > divideLines( size_t &, string_view );


template< typename STRING_T >
stdVector< STRING_T > wrapTextToMaxLength( stdVector< STRING_T > const & lines,
                                           size_t & maxLineLength )
{
  if( lines.empty())
    return lines;

  size_t effectiveMaxLineLength = 0;

  stdVector< STRING_T > formattedLines;
  formattedLines.reserve( lines.size() );
  for( const auto & line : lines )
  {
    size_t startPos = 0;

    while( startPos < line.size())
    {
      // if the remaining part is shorter than maxLineLength
      if( startPos + maxLineLength >= line.size())
      {
        formattedLines.push_back( STRING_T( ltrimSpaces( line.substr( startPos ))));
        effectiveMaxLineLength = std::max( effectiveMaxLineLength, formattedLines.back().size() );
        break;
      }

      // find last space occurence before maxLineLength
      size_t const endPos = startPos + maxLineLength;
      size_t const spacePos = line.rfind( ' ', endPos );
      if( spacePos != STRING_T::npos && spacePos > startPos )
      {
        // cut and push at the last space found
        formattedLines.push_back( STRING_T( ltrimSpaces( line.substr( startPos, spacePos - startPos ))));
        startPos = spacePos + 1;
      }
      else
      {
        // no space found, cut in the middle of the word with maxLineLength
        formattedLines.push_back( STRING_T( ltrimSpaces( line.substr( startPos, maxLineLength ))));
        startPos += maxLineLength;
      }
      effectiveMaxLineLength = std::max( effectiveMaxLineLength, formattedLines.back().size() );
    }
  }

  maxLineLength = effectiveMaxLineLength;
  return formattedLines;
}
template stdVector< string > wrapTextToMaxLength( stdVector< string > const &, size_t & );
template stdVector< string_view > wrapTextToMaxLength( stdVector< string_view > const &, size_t & );

}
}
