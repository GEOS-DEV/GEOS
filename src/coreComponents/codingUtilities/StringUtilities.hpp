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
 * @file StringUtilities.hpp
 */

#ifndef GEOSX_CODINGUTILITIES_STRINGUTILITIES_HPP_
#define GEOSX_CODINGUTILITIES_STRINGUTILITIES_HPP_

#include <cxxabi.h>
#include <cstring>
#include <memory>
#include <sstream>
#include <algorithm>
#include <map>


#include "common/DataTypes.hpp"
#include "cxx-utilities/src/IntegerConversion.hpp"

namespace geosx
{
namespace stringutilities
{
template <class Type>
std::string toString(Type theVar);

template <class Type>
Type fromString(std::string theString);

//////////////////////////////////////////////////

//
// toString
//
/// convert a variable to a string
template <class Type>
std::string toString(Type theVar){
  std::ostringstream oss;
  oss << theVar;
  return oss.str();
}
// override the template for string->string
template <>
inline std::string toString<std::string>(std::string theVar)
{
  return theVar;
}

//
// fromString
//
/// returns a variable from a string
template <class Type>
Type fromString(std::string theString){
  std::istringstream iss(theString);
  Type theVar;
  iss >> theVar;
  return theVar;
}
// override the template for string->string
template <>
inline std::string fromString<std::string>(std::string theVar)
{
  return theVar;
}

// override the template for string->real64
// Allows unit manager to convert units
template <>
real64 fromString<real64>(std::string theVar);

// override the template for FieldType
//template <>
//inline FieldType fromString<FieldType>(std::string theString){
//  using namespace FieldInfo;
//  if(theString==FieldInfo::IntegerStr){
//    return integerField;
//  } else if(theString==FieldInfo::GlobalIndexStr){
//  return globalIndexField;
//  } else if(theString==FieldInfo::LocalIndexStr){
//  return localIndexField;
//  } else if(theString==FieldInfo::RealStr){
//    return realField;
//  } else if(theString==FieldInfo::R1TensorStr){
//    return R1TensorField;
//  } else if(theString==FieldInfo::R2TensorStr){
//    return R2TensorField;
//  } else if(theString==FieldInfo::R2SymTensorStr){
//    return R2SymTensorField;
//  }else {
//    throw GPException("Error fromString: unrecognized field type: " +
// theString +".");
//  }
//}

/// Convert a string to lowercase
void toLower(std::string& theString);
std::string lowercase(std::string theString);

/// Convert a string to uppercase
void toUpper(std::string& theString);
std::string uppercase(std::string theString);

/// Check for case insensitive equality between strings
bool ieq(std::string strA,std::string strB);

/// Overloaded function to check equality between strings and char arrays
/// Mainly used to avoid char*==char* mistakes
inline bool streq(const std::string& strA, const std::string& strB){
  return strA == strB;
}
inline bool streq(const std::string& strA, const char * strB){
  return strA == strB;
}
inline bool streq(const char * strA, const std::string& strB){
  return strA == strB;
}
inline bool streq(const char * strA, const char * strB){
  return !strcmp(strA, strB);
}

/// string is integer
inline bool strIsInt(std::string theString){
  std::istringstream iss(theString);
  int dummy;
  iss >> dummy;
  return !iss.fail();
}

/// Subdivide string by delimiters
string_array Tokenize(const std::string& str, const std::string& delimiters);

/// Subdivide string delimited by sequence of characters
string_array TokenizeSeq(const std::string& str, const std::string& seq);

/// Split string at first token
string_array Split(const std::string& str, const std::string& delimiters);


/// Expand string vector based on multiple tokens eg [a, b**3, c] => [a,b,b,b,c]
inline void ExpandMultipleTokens(string_array& sVector, const std::string& multipleToken="**"){
  localIndex n= integer_conversion<localIndex>(sVector.size());
  string_array newVec;
  for( localIndex i =0 ; i < n ; ++i)
  {
    string_array keyMult = TokenizeSeq(sVector[i], multipleToken);
    if( (keyMult.size() == 2) && strIsInt(keyMult[1]) )
    {
      int numMult = fromString<int>(keyMult[1]);
      for(int j=0 ; j < numMult ; ++j )
      {
        newVec.push_back(keyMult[0]);
      }
    }
    else
    {
      newVec.push_back(sVector[i]);
    }

  }
  sVector = newVec;
}


static const std::string base64Chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                       "abcdefghijklmnopqrstuvwxyz"
                                       "0123456789+/";

inline string & EncodeBase64( unsigned char const * const bytes,
                              string & outputString,
                              int len)
{
  char * out = &outputString[0];
  integer val = 0;
  integer valB = -6;
  integer size = 0;

  for( integer i = 0 ; i < len ; i++ )
  {
    val = ( val << 8 ) + bytes[i];
    valB += 8;
    while ( valB >= 0 )
    {
      *out = base64Chars[ ( val>>valB ) &0x3F ] ; //0x3f is the Hexadecimal for 63
      ++out;
      ++size;
      valB -= 6;
    }
  }
  if( valB > -6 )
  {
    *out = base64Chars[ ( ( val << 8 ) >> ( valB + 8 ) ) &0x3F ];
    ++out;
    ++size;
  }
  while( size % 4 )
  {
    *out = '=';
    ++out;
    ++size;
  }
  return outputString;
}

}
}

#endif /* GEOSX_CODINGUTILITIES_STRINGUTILITIES_HPP_ */
