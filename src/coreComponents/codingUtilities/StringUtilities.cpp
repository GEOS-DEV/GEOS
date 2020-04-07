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

/**
 * @file StringUtilities.cpp
 */

#include "math/TensorT/TensorT.h"
#include "codingUtilities/StringUtilities.hpp"
//#include "codingUtilities/UnitManager.h"

#include <stdarg.h>

namespace geosx
{
namespace stringutilities
{

/**
 * String tokenizing function
 **/
string_array Tokenize( const std::string & str, const std::string & delimiters )
{
  string_array tokens;

  if( str.length() == 0 )
  {
    tokens.push_back( str );
  }
  else
  {

    bool usesNonWhitespaceDelimiters = false;
    std::string::size_type i =0;
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
      while( (newPos=str.find_first_of( delimiters, lastPos )) != std::string::npos )
      {
        tokens.push_back( str.substr( lastPos, newPos-lastPos ));
        lastPos = newPos + 1;
      }
      tokens.push_back( str.substr( lastPos, str.length()-lastPos ));
    }
    else
    {
      // whitespace delimiters
      // skip multiple adjacent delimiters
      size_t lastPos = str.find_first_not_of( delimiters, 0 );
      lastPos = (lastPos == std::string::npos) ? 0 : lastPos;

      size_t newPos = lastPos;
      while( (newPos=str.find_first_of( delimiters, lastPos )) != std::string::npos )
      {
        tokens.push_back( str.substr( lastPos, newPos-lastPos ));
        lastPos = str.find_first_not_of( delimiters, newPos );
      }
      if( lastPos!= std::string::npos )
        tokens.push_back( str.substr( lastPos, str.length()-lastPos ));

    }
  }
  return tokens;
}

integer FindBase64StringLength( integer dataSize )
{
  integer base64StringLength = (dataSize * 8) / 6;
  while( base64StringLength % 4 )
  {
    base64StringLength++;
  }
  return base64StringLength;
}

static const std::string base64Chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                       "abcdefghijklmnopqrstuvwxyz"
                                       "0123456789+/";

string EncodeBase64( unsigned char const * const bytes,
                     integer dataSize )
{
  string output;
  output.reserve( FindBase64StringLength( dataSize ) );
  integer val = 0;
  integer valB = -6;
  integer size = 0;

  for( integer i = 0; i < dataSize; i++ )
  {
    val = ( val << 8 ) + bytes[i];
    valB += 8;
    while( valB >= 0 )
    {
      output.push_back( base64Chars[ ( val>>valB ) &0x3F ] );  //0x3f is the Hexadecimal for 63
      ++size;
      valB -= 6;
    }
  }
  if( valB > -6 )
  {
    output.push_back( base64Chars[ ( ( val << 8 ) >> ( valB + 8 ) ) &0x3F ] );
    ++size;
  }
  while( size % 4 )
  {
    output.push_back( '=' );
    ++size;
  }
  return output;
}

}
}
