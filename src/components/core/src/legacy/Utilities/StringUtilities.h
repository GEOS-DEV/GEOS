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
#include <slic/slic.hpp>


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
  if(theString==FieldInfo::IntegerStr){
  	return integerField;
  } else if(theString==FieldInfo::GlobalIndexStr){
	return globalIndexField;
  } else if(theString==FieldInfo::LocalIndexStr){
	return localIndexField;
  } else if(theString==FieldInfo::RealStr){
  	return realField;
  } else if(theString==FieldInfo::R1TensorStr){
  	return R1TensorField;
  } else if(theString==FieldInfo::R2TensorStr){
  	return R2TensorField;
  } else if(theString==FieldInfo::R2SymTensorStr){
  	return R2SymTensorField;
  }else {
    SLIC_ERROR("Error fromString: unrecognized field type: " + theString +".");
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
sArray1d Tokenize(const std::string& str, const std::string& delimiters);

/// Subdivide string delimited by sequence of characters
sArray1d TokenizeSeq(const std::string& str, const std::string& seq);

/// Split string at first token
sArray1d Split(const std::string& str, const std::string& delimiters);

/// Remove comments from end of string
void RemoveComments(std::string& str, char d='%');

/// Remove all spaces ' ' from a string
inline void RemoveSpaces(std::string& aString){
  aString.erase(std::remove(aString.begin(), aString.end(), ' '), aString.end());
}

/// Expand string vector based on multiple tokens eg [a, b**3, c] => [a,b,b,b,c]
inline void ExpandMultipleTokens(sArray1d& sVector , const std::string& multipleToken="**"){
	int n= sVector.size();
	sArray1d newVec;
	for(int i =0; i < n; ++i){
		sArray1d keyMult = TokenizeSeq(sVector[i], multipleToken);
		if( (keyMult.size() == 2) && strIsInt(keyMult[1]) ){
		  int numMult = fromString<int>(keyMult[1]);
		  for(int j=0; j < numMult; ++j ){
		      newVec.push_back(keyMult[0]);
		  }
	    } else {
	      newVec.push_back(sVector[i]);
		}

	}
	sVector.swap(newVec);
}

/// Trim whitespace from string
void TrimLeft(std::string& str, const std::string& d=" \t\n\r");
void TrimRight(std::string& str, const std::string& d=" \t\n\r");
void Trim(std::string& str, const std::string& d=" \t\n\r");

inline void Trim(sArray1d& strVect, const std::string& d=" \t\n\r"){
  for(unsigned i =0; i < strVect.size(); ++i) Trim(strVect[i],d);	
}

/// Replace parameters of form "$:NAME" in a string, returns true if a parameter is detected
/// @param lineStr - the string containing the parameters 
/// @param parameterMap - map of parameter names and replacement strings, 
/// @param prefix - Character sequence used to signal start of parameter ("$:" by default),
bool ReplaceParameters(std::string& lineStr, const std::map<std::string,std::string>& parameterMap,const std::string& prefix= "$:");


#endif /*STRINGUTILITIES_H_*/
