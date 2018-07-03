/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef STRINGUTILITIES_H_
#define STRINGUTILITIES_H_

/**
 * @file StringUtilities.h
 * @author walsh24
 */

#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

#include "legacy/Common/Common.h"

#ifdef USE_ATK
#include <slic/slic.hpp>
#endif


/////////////////////////////////////////////////
// Forward declaration of templated functions

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
{ return theVar; }

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
{ return theVar; }

// override the template for string->realT
// Allows unit manager to convert units
template <>
realT fromString<realT>(std::string theVar);

// override the template for FieldType
template <>
inline FieldType fromString<FieldType>(std::string theString){
  using namespace FieldInfo;
  if(theString==FieldInfo::IntegerStr)
  {
    return integerField;
  }
  else if(theString==FieldInfo::GlobalIndexStr)
  {
    return globalIndexField;
  }
  else if(theString==FieldInfo::LocalIndexStr)
  {
    return localIndexField;
  }
  else if(theString==FieldInfo::RealStr)
  {
    return realField;
  }
  else if(theString==FieldInfo::R1TensorStr)
  {
    return R1TensorField;
  }
  else if(theString==FieldInfo::R2TensorStr)
  {
    return R2TensorField;
  }
  else if(theString==FieldInfo::R2SymTensorStr)
  {
    return R2SymTensorField;
  }
  else
  {
#ifdef USE_ATK
    SLIC_ERROR("Error fromString: unrecognized field type: " + theString +".");
#endif
  }
}

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
inline bool streq(const std::string& strA, const std::string& strB){ return strA == strB;}
inline bool streq(const std::string& strA, const char * strB){ return strA == strB;}
inline bool streq(const char * strA, const std::string& strB){ return strA == strB;}
inline bool streq(const char * strA, const char * strB){ return !strcmp(strA, strB);}

/// string is integer
inline bool strIsInt(std::string theString){
  std::istringstream iss(theString);
  int dummy;
  iss >> dummy;
  return !iss.fail();
}

/// Subdivide string by delimiters
array<string> Tokenize(const std::string& str, const std::string& delimiters);

/// Subdivide string delimited by sequence of characters
array<string> TokenizeSeq(const std::string& str, const std::string& seq);

/// Split string at first token
array<string> Split(const std::string& str, const std::string& delimiters);

/// Remove comments from end of string
void RemoveComments(std::string& str, char d='%');

/// Remove all spaces ' ' from a string
inline void RemoveSpaces(std::string& aString){
  aString.erase(std::remove(aString.begin(), aString.end(), ' '), aString.end());
}

/// Expand string vector based on multiple tokens eg [a, b**3, c] => [a,b,b,b,c]
inline void ExpandMultipleTokens(array<string>& sVector, const std::string& multipleToken="**"){
  int n= sVector.size();
  array<string> newVec;
  for(int i =0 ; i < n ; ++i)
  {
    array<string> keyMult = TokenizeSeq(sVector[i], multipleToken);
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
  sVector.swap(newVec);
}

/// Trim whitespace from string
void TrimLeft(std::string& str, const std::string& d=" \t\n\r");
void TrimRight(std::string& str, const std::string& d=" \t\n\r");
void Trim(std::string& str, const std::string& d=" \t\n\r");

inline void Trim(array<string>& strVect, const std::string& d=" \t\n\r"){
  for(unsigned i =0 ; i < strVect.size() ; ++i)
    Trim(strVect[i],d);
}

/// Replace parameters of form "$:NAME" in a string, returns true if a parameter
// is detected
/// @param lineStr - the string containing the parameters
/// @param parameterMap - map of parameter names and replacement strings,
/// @param prefix - Character sequence used to signal start of parameter ("$:"
// by default),
bool ReplaceParameters(std::string& lineStr, const std::map<std::string,std::string>& parameterMap,const std::string& prefix= "$:");


#endif /*STRINGUTILITIES_H_*/
