//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file StringUtilities.cpp
 * @author walsh24
 */

#include "math/TensorT/TensorT.h"
#include "codingUtilities/StringUtilities.hpp"
#include "codingUtilities/UnitManager.h"

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

bool ieq(std::string strA,std::string strB){ toLower(strA); toLower(strB); return strA == strB; }

// override the template for string->real64
// Allows for unit conversion
template<>
real64 fromString<real64>(std::string theVar)
{  
  	Trim(theVar);
	UnitManager& um = UnitManager::Instance();
	//std::cout  << theVar << " " << um.Convert(theVar) <<std::endl; 
	return um.Convert(theVar);
}

template<>
R2SymTensor fromString<R2SymTensor>(std::string theVar)
{
  Trim(theVar);
  UnitManager& um = UnitManager::Instance();
  R2SymTensor var;
  //std::cout  << theVar << " " << um.Convert(theVar) <<std::endl;

  string_array svalues = Tokenize(theVar,",");

  for( int i=0 ; i<R2SymTensor::Length() ; ++i )
  {
    var.Data()[i] = um.Convert(svalues[i]);
  }
  return var;
}

/**
 * String tokenizing function
**/
string_array Tokenize(const std::string& str, const std::string& delimiters)
{
  string_array tokens;

  if ( str.length() == 0 ) {
    tokens.push_back(str);
  } else {
  	
    bool usesNonWhitespaceDelimiters = false;
    int i =0;
    while(delimiters[i] && !usesNonWhitespaceDelimiters)
    {
      usesNonWhitespaceDelimiters |= !isspace( int(delimiters[i]) );
      ++i;
    }
    
    if(usesNonWhitespaceDelimiters){
      // do not skip multiple adjacent delimiters - indicates empty strings
      size_t lastPos = 0;

      size_t newPos = lastPos;
      while ( (newPos=str.find_first_of(delimiters, lastPos)) != std::string::npos ) {
        tokens.push_back(str.substr(lastPos, newPos-lastPos));
        lastPos = newPos + 1;
      }
      tokens.push_back(str.substr(lastPos, str.length()-lastPos));
    } else {
      // whitespace delimiters
      // skip multiple adjacent delimiters
      size_t lastPos = str.find_first_not_of(delimiters,0);
      lastPos = (lastPos == std::string::npos)? 0:lastPos;

      size_t newPos = lastPos;
      while ( (newPos=str.find_first_of(delimiters, lastPos)) != std::string::npos ) {
        tokens.push_back(str.substr(lastPos, newPos-lastPos));
        lastPos = str.find_first_not_of(delimiters,newPos);
      }
      if(lastPos!= std::string::npos) tokens.push_back(str.substr(lastPos, str.length()-lastPos));
      
    }
  }
  return tokens;
}

/**
 * String tokenizing function using a character sequence,
 * ie. Tokenize "1.0000 HPO4-- +1.0000 Cu++" with " +" gives {"1.0000 HPO4--","1.0000 Cu++"}
**/
string_array TokenizeSeq(const std::string& str, const std::string& seq)
{
  string_array tokens;
  
  if ( str.length() == 0 ) {
    tokens.push_back(str);
  } else {

    size_t lastPos = 0;

    size_t newPos = lastPos;
    while ( (newPos=str.find(seq, lastPos)) != std::string::npos ) {
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

  if ( str.length() == 0 ) {
    tokens.push_back(str);
  } else {

    size_t pos = str.find_first_of(delimiters,0);
    
    tokens.push_back(str.substr(0, pos));
    pos++;
    if(pos < str.size()) tokens.push_back(str.substr(pos));
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
void Trim(std::string& str, const std::string& trimChars) { TrimLeft(str,trimChars);  TrimRight(str,trimChars); }

/// Replace parameters of form "$:NAME" in a string, returns true if a parameter is detected
/// @param lineStr The string containing the parameters.
/// @param parameterMap Map of parameter names and replacement strings.
/// @param prefix Character sequence used to signal start of parameter ("$:" by default).
bool ReplaceParameters(std::string& lineStr, const std::map<std::string,std::string>& parameterMap,const std::string& prefix){

  bool rv = false;
  std::map<std::string,std::string>::const_iterator endMap = parameterMap.end();

  const std::string validParamChars = std::string("abcdefghijklmnopqrstuvwxyz") +
                                      std::string("ABCDEFGHIJKLMNOPQRSTUVWXYZ") +
                                      std::string("0123456789_");

  size_t pSize = prefix.size();
  size_t startIndx = lineStr.find(prefix);
//  size_t endIndx =0;

  while(startIndx < lineStr.size()){
  	rv = true;
    size_t startIndxB = startIndx+pSize;
    size_t endIndx = lineStr.find_first_not_of(validParamChars,startIndxB);

    std::string paramName = lineStr.substr(startIndxB, endIndx-startIndxB);

    std::map<std::string,std::string>::const_iterator itr = parameterMap.find(paramName);
    if(itr == endMap){
      // throw GPException("Error: Undefined model parameter: " + paramName + ".");
      std::map<std::string,std::string>::const_iterator itrB  = parameterMap.begin();
      while(itrB != parameterMap.end()){
      	std::cout  << itrB->first << std::endl;
      	++itrB;
      }
      std::cout << "Error: Undefined model parameter: " + paramName + "." << std::endl; exit(1);
    }
    
    const std::string& replaceStr = itr->second;
    lineStr.replace(startIndx,endIndx-startIndx,replaceStr);

    startIndx = lineStr.find(prefix);
  }
  
  return rv;
}

}
}


