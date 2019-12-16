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

void toLower(std::string& theString){
  std::transform(theString.begin(), theString.end(), theString.begin(), ::tolower);
}
std::string lowercase(std::string theString){
  std::transform(theString.begin(), theString.end(), theString.begin(), ::tolower);
  return theString;
}

void toUpper(std::string& theString){
  std::transform(theString.begin(), theString.end(), theString.begin(), ::toupper);
}
std::string uppercase(std::string theString){
  std::transform(theString.begin(), theString.end(), theString.begin(), ::toupper);
  return theString;
}

bool ieq(std::string strA,std::string strB){
  toLower(strA); toLower(strB); return strA == strB;
}



/**
 * String tokenizing function
 **/
string_array Tokenize(const std::string& str, const std::string& delimiters)
{
  string_array tokens;

  if ( str.length() == 0 )
  {
    tokens.push_back(str);
  }
  else
  {

    bool usesNonWhitespaceDelimiters = false;
    std::string::size_type i =0;
    while(delimiters[i] && !usesNonWhitespaceDelimiters)
    {
      usesNonWhitespaceDelimiters |= !isspace( int(delimiters[i]) );
      ++i;
    }

    if(usesNonWhitespaceDelimiters)
    {
      // do not skip multiple adjacent delimiters - indicates empty strings
      size_t lastPos = 0;

      size_t newPos = lastPos;
      while ( (newPos=str.find_first_of(delimiters, lastPos)) != std::string::npos )
      {
        tokens.push_back(str.substr(lastPos, newPos-lastPos));
        lastPos = newPos + 1;
      }
      tokens.push_back(str.substr(lastPos, str.length()-lastPos));
    }
    else
    {
      // whitespace delimiters
      // skip multiple adjacent delimiters
      size_t lastPos = str.find_first_not_of(delimiters,0);
      lastPos = (lastPos == std::string::npos) ? 0 : lastPos;

      size_t newPos = lastPos;
      while ( (newPos=str.find_first_of(delimiters, lastPos)) != std::string::npos )
      {
        tokens.push_back(str.substr(lastPos, newPos-lastPos));
        lastPos = str.find_first_not_of(delimiters,newPos);
      }
      if(lastPos!= std::string::npos)
        tokens.push_back(str.substr(lastPos, str.length()-lastPos));

    }
  }
  return tokens;
}

/**
 * String tokenizing function using a character sequence,
 * ie. Tokenize "1.0000 HPO4-- +1.0000 Cu++" with " +" gives {"1.0000
 * HPO4--","1.0000 Cu++"}
 **/
string_array TokenizeSeq(const std::string& str, const std::string& seq)
{
  string_array tokens;

  if ( str.length() == 0 )
  {
    tokens.push_back(str);
  }
  else
  {

    size_t lastPos = 0;

    size_t newPos = lastPos;
    while ( (newPos=str.find(seq, lastPos)) != std::string::npos )
    {
      tokens.push_back(str.substr(lastPos, newPos-lastPos));
      lastPos = newPos + seq.size();
    }
    tokens.push_back(str.substr(lastPos, str.length()-lastPos));
  }
  return tokens;
}

/**
 * Split string at deliminator
 **/
string_array Split(const std::string& str, const std::string& delimiters)
{
  string_array tokens;

  if ( str.length() == 0 )
  {
    tokens.push_back(str);
  }
  else
  {

    size_t pos = str.find_first_of(delimiters,0);

    tokens.push_back(str.substr(0, pos));
    pos++;
    if(pos < str.size())
      tokens.push_back(str.substr(pos));
  }
  return tokens;
}

/**
   Remove everything in a string after a comment character (eg '%')
 **/
void RemoveComments(std::string& str, char d)
{
  size_t indx = str.find(d);
  str = str.substr(0,indx);
}

/**
   Remove white space from left of string
 **/
void TrimLeft(std::string& str, const std::string& trimChars)
{
  str.erase(0,str.find_first_not_of(trimChars));
}

/**
   Remove white space from right of string
 **/
void TrimRight(std::string& str, const std::string& trimChars) {
  str.erase(str.find_last_not_of(trimChars)+1);
}

/**
   Remove white space around string
 **/
void Trim(std::string& str, const std::string& trimChars) {
  TrimLeft(str,trimChars);  TrimRight(str,trimChars);
}

/// Replace parameters of form "$:NAME" in a string, returns true if a parameter
// is detected
/// @param lineStr The string containing the parameters.
/// @param parameterMap Map of parameter names and replacement strings.
/// @param prefix Character sequence used to signal start of parameter ("$:" by
// default).
bool ReplaceParameters(std::string& lineStr, const std::map<std::string,std::string>& parameterMap,const std::string& prefix){

  bool rv = false;
  std::map<std::string,std::string>::const_iterator endMap = parameterMap.end();

  const std::string validParamChars = std::string("abcdefghijklmnopqrstuvwxyz") +
                                      std::string("ABCDEFGHIJKLMNOPQRSTUVWXYZ") +
                                      std::string("0123456789_");

  size_t pSize = prefix.size();
  size_t startIndx = lineStr.find(prefix);
//  size_t endIndx =0;

  while(startIndx < lineStr.size())
  {
    rv = true;
    size_t startIndxB = startIndx+pSize;
    size_t endIndx = lineStr.find_first_not_of(validParamChars,startIndxB);

    std::string paramName = lineStr.substr(startIndxB, endIndx-startIndxB);

    std::map<std::string,std::string>::const_iterator itr = parameterMap.find(paramName);
    if(itr == endMap)
    {
      std::map<std::string,std::string>::const_iterator itrB  = parameterMap.begin();
      while(itrB != parameterMap.end())
      {
        ++itrB;
      }

      GEOSX_ERROR("Error: Undefined model parameter: " << paramName << ".");
    }

    const std::string& replaceStr = itr->second;
    lineStr.replace(startIndx,endIndx-startIndx,replaceStr);

    startIndx = lineStr.find(prefix);
  }

  return rv;
}

}
}
