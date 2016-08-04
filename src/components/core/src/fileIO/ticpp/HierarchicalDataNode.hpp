#ifndef __HIERARCHICALDATAPARSER__
#define __HIERARCHICALDATAPARSER__

//#include "Common/typedefs.h"
#include "codingUtilities/StringUtilities.hpp"
#include <map>
#include <stdarg.h>
#include <limits>
#include <string>
#include <sstream>
#include <iostream>
#include <cctype>
#include <algorithm>

#include "math/TensorT/TensorT.h"
#include "slic/slic.hpp"

#include "common/DataTypes.hpp"
//SLIC_ERROR( msg )

const int nsdof = 3;
const int nsdof2 = nsdof * nsdof;

/*
   The MIT License

   Copyright (c) 2011 Grainflow Dynamics, Inc. (www.grainflow.com)

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
 */
namespace geosx
{
namespace TICPP
{
/**
 * @author Scott Johnson
 * This holds key-child and key-value pairs to form a
 * hierarchical arrangement of data
 */
class HierarchicalDataNode
{
public:
  HierarchicalDataNode() :
    m_level(-1), m_heading(""), m_it(), m_attributes(), m_itc(), m_children()
  {}

  HierarchicalDataNode(const int level, const char* name = 0) :
    m_level(-1), m_heading(""), m_it(), m_attributes(), m_itc(), m_children()
  {
    if( name )
      SetHeading(name);
    this->m_level = level;
  }

  ~HierarchicalDataNode()
  {}

  inline int Level() const
  {
    return this->m_level;
  }

  void SetLevel(int level)
  {
    m_level = level;
    std::vector<HierarchicalDataNode>::iterator itr = m_children.begin();
    std::vector<HierarchicalDataNode>::iterator iend = m_children.end();
    for( ; itr != iend ; ++itr)
    {
      itr->SetLevel(level+1);
    }
  }

  inline const char*
  Heading() const
  {
    return m_heading.c_str();
  }

  inline void SetHeading(const char* heading)
  {
    this->m_heading = std::string(heading);
//		StrToLower(this->heading);
  }

  inline static void StrToLower(std::string& str)
  {
    //std::transform(str.begin(), str.end(), str.begin(), tolower);
    for( uint64 i = 0 ; i < str.length() ; i++ )
      if( str[i] >= 'A' && str[i] <= 'Z' )
        str[i] = str[i] + ('a' - 'A');
  }

  inline static bool StrEq(const std::string& str, const std::string& str2)
  {
    return !strcmp(str.c_str(), str2.c_str());
  }


  template<class T>
  inline static bool StrVal(const std::string& str,
                            const std::map<std::string, std::string>::const_iterator& it, T& val)

  {
    if( StrEq(str, it->first) )
    {
      val = StrVal<T>(it->second);
      return true;
    }
    return false;
  }

  template<class T>
  inline static T StrVal(const std::string& str)
  {
    T val = (T)0;    // explicit casting needed for enums
    std::istringstream iss(str, std::istringstream::in);

    SLIC_ERROR_IF( !(iss >> val), str );

    return val;
  }

  template<class T>
  inline static void StrVal(const std::string& str, T* vval, const uint64 length)

  {
    std::istringstream iss(str, std::istringstream::in);
    for( uint64 i = 0 ; i < length ; i++ )
      if( !(iss >> vval[i]) )
        SLIC_ERROR_IF( !(iss >> vval[i]), str );

  }

  inline static void StrVal(const std::string& str, real64_array& vval)
  {
    real64 temp;
    vval.resize(0);
    std::istringstream iss(str, std::istringstream::in);
    while( (iss >> temp) )
    {
      vval.push_back(temp);
    }
  }

  inline static bool StrVal(const std::string& str,
                            const std::map<std::string, std::string>::const_iterator& it
                            , real64* val,
                            uint64 vectSize = nsdof)
  {
    if( StrEq(str, it->first) )
    {
      std::istringstream iss(it->second, std::istringstream::in);
      for( uint64 i = 0 ; i < vectSize ; i++ )
        SLIC_ERROR_IF( !(iss >> val[i]), it->first );
      return static_cast<bool>(iss);
    }
    return false;
  }

  inline void ResetCursors()
  {
    m_it = m_attributes.begin();
  }

  inline void ResetChildCursor()
  {
    m_itc = m_children.begin();
  }

  inline bool Next(std::string& key, std::string& value, const bool reset = false)
  {
    if( m_attributes.size() == 0 )
      return false;
    if( reset )
      ResetCursors();
    if( m_it == m_attributes.end() )
      return false;
    key = m_it->first;
    value = m_it->second;
    ++m_it;
    return true;
  }

  /**
   * Returns the string associated with attribute "key".
   * Returns empty string if key is not found.
   */
  inline std::string GetAttributeString(const std::string& key) const
  {
    std::map<std::string, std::string>::const_iterator it = m_attributes.find(key);
    if( it == m_attributes.end() )
    {
      return "";
    }
    return it->second;
  }

  inline bool HasAttribute(const std::string& key) const {
    std::map<std::string, std::string>::const_iterator it = m_attributes.find(key);
    return it != m_attributes.end();
  }

  /**
   * Returns the value associated with attribute "key".
   * Should return an error if key is not found.
   */
  template<class T>
  T GetAttributeValue(const std::string& key) const
  {
    std::string str = GetAttributeString(key);
    SLIC_ERROR_IF( str.empty(), "Error: Required value was not given for attribute "+ key );

    return StrVal<T>(str);
  }


  template<class T>
  std::vector<T> GetAttributeVector(const std::string& key, std::string token = " \t\n") const;


  string_array GetStringVector(const std::string& key, std::string token = " \t\n") const
  {
    string_array vec;
    std::string varsStr = GetAttributeString(key);
    if( !varsStr.empty() )
    {
      vec = stringutilities::Tokenize(varsStr, token);
    }
    return vec;
  }

  template<class T>
  std::vector<T> GetAttributeVector(const char * key, std::string token = " \t\n") const
  {
    return GetAttributeVector<T>(std::string(key), token);
  }


  template<class T>
  std::vector<T> GetAttributeVectorOrDefault(const std::string& keyStr, std::vector<T> def, std::string token = " \t\n") const
  {
    if(this->GetAttributeString(keyStr).length() == 0)
    {
      return def;
    }
    else
      return GetAttributeVector<T>(keyStr, token);
  }


  template<class T>
  std::vector<T> GetAttributeVectorOrDefault(const char * key, const std::vector<T>& def, std::string token = " \t\n") const
  {
    return GetAttributeVectorOrDefault<T>(std::string(key), def, token);
  }

  R1Tensor GetAttributeTensor(const char * key) const
  {
    std::vector<real64> tmp = GetAttributeVector<real64>(std::string(key), " ");
    SLIC_ERROR_IF( tmp.size() != nsdof, "HierarchicalDataNode::GetAttributeTensor - Wrong size for a R1Tensor");

    R1Tensor ret(tmp.data());
    return ret;
  }

  R1Tensor GetAttributeTensorOrDefault(const char * key, const R1Tensor& def) const
  {
    std::string keyStr(key);
    if(this->GetAttributeString(keyStr).length() == 0)
    {
      R1Tensor ret(def);
      return ret;
    }
    else
      return GetAttributeTensor(key);
  }

  R2Tensor GetAttributeR2Tensor(const char * key) const
  {
    std::vector<real64> tmp = GetAttributeVector<real64>(std::string(key), " ");
    SLIC_ERROR_IF( tmp.size() != nsdof*nsdof,
                   "HierarchicalDataNode::GetAttributeR2Tensor: Size " + stringutilities::toString( tmp.size() ) +
                   std::string(" - Wrong size for a R2Tensor") );
    R2Tensor ret(tmp.data());
    return ret;
  }

  R2Tensor GetAttributeR2TensorOrDefault(const char * key, const R2Tensor& def) const
  {
    std::string keyStr(key);
    if(this->GetAttributeString(keyStr).length() == 0)
    {
      R2Tensor ret(def);
      return ret;
    }
    else
      return GetAttributeR2Tensor(key);
  }

  /**
   * @author walsh24
   *
   * Returns a real value associated with attribute "key".
   * If key is not found returns value interpreted from the default string.
   *
   * This function can be used to return default values defined in terms of the problem units:
   * eg. real64 pressure = GetAttributeOrDefault("pressure", "10 MPa");
   *
   */
  real64 GetAttributeOrDefault(const std::string& key, const std::string& def);
  real64 GetAttributeOrDefault(const std::string& key, const char * def)
  {
    return GetAttributeOrDefault(key, std::string(def));
  }

  /**
   * @author walsh24
   *
   * Returns a value associated with attribute "key".
   * If key is not found returns a default value.
   *
   */
  template<class T>
  T GetAttributeOrDefault(const std::string& key, T def)
  {
    std::string val = GetAttributeString(key);
    if( val == "" )
    {

      switch(m_defaultReportLevel)
      {
      case disableDefaults:
      {
        std::string exceptionStr = std::string("Error: Attempt to set a default value when default values have been disabled.\n") +
                                   std::string("Parameter: ") + key + "\n" +
                                   std::string("Default Value: ") + geosx::stringutilities::toString<T>(def) + "\n";
        SLIC_ERROR(exceptionStr);
        break;
      }
      case reportDefaults:
        std::cout << "Setting " + key + " to default value: " + stringutilities::toString<T>(def)  << std::endl;
      // fall through to next case
      case recordDefaults:
        AddAttributePair(key, stringutilities::toString<T>(def) );    // add attribute to hdn - default will be included if xml written to file
      // fall through to next case
      case silent:
        break;
      }
      return def;
      /*
         #ifdef DISABLE_DEFAULT_VALUES
         std::string exceptionStr = "Error: Attempt to set a default value when default values have been prohibited.\n" +
         "Parameter: " + key + "\n" +
         "Default Value: " + toString<T>(def) + "\n";
         throw GPException(exceptionStr);
         #else
         #ifdef REPORT_DEFAULT_VALUES
         std::cout << "Setting " + key + " to default value: " + toString<T>(def)  << std::endl;
         #endif
         return def;
       #endif*/
    }
    else
    {
      return StrVal<T>(val);
    }
  }

  std::string GetAttributeStringOrDefault(const std::string& key, const std::string& def)
  {
    std::string val = GetAttributeString(key);
    if( val == "" )
    {
      switch(m_defaultReportLevel)
      {
      case disableDefaults:
      {
        std::string exceptionStr = std::string("Error: Attempt to set a default value when default values have been prohibited.\n") +
                                   "Parameter: " + key + "\n" +
                                   "Default Value: " + def + "\n";
        SLIC_ERROR(exceptionStr);
        break;
      }
      case reportDefaults:
        std::cout << "Setting " + key + " to default value: " + def  << std::endl;
      // fall through to next case
      case recordDefaults:
        AddAttributePair(key,def);    // add attribute to hdn - default will be included if xml written to file
      // fall through to next case
      case silent:
        break;
      default:
        break;
      }
      return def;

    }
    else
    {
      return val;
    }
  }

  inline void AddAttributePair(const std::string key, const std::string value)
  {
    /* std::string keyLC(key);StrToLower(keyLC);attributes[keyLC] = value;*/
    m_attributes[key] = value;
  }

  inline const std::map<std::string, std::string>&
  GetAttributes() const
  {
    return m_attributes;
  }

  inline uint64 Count() const
  {
    return m_attributes.size();
  }

  inline HierarchicalDataNode*
  NewChild(const char* header = 0)
  {
    SLIC_ERROR_IF( !header, "HierarchicalDataNode::NewChild Cannot create a child without a valid header" );

    std::string heading(header);
    int nextLevel = this->Level();
    ++nextLevel;
    m_children.push_back(HierarchicalDataNode(nextLevel, heading.c_str()));
    if( m_children.size() == 1 )
      ResetChildCursor();
    return &(*m_children.rbegin());
  }

  /*
     inline HierarchicalDataNode*
     NewFirstChild(const char* header = 0)
     {
     if( !header )
      throw GPException("HierarchicalDataNode::NewChild Cannot create a child without a valid header");
     std::string heading(header);
     int nextLevel = this->Level();
     ++nextLevel;
     m_children.push_back(HierarchicalDataNode(nextLevel, heading.c_str()));
     if( m_children.size() == 1 )
      ResetChildCursor();
     return &(*m_children.rbegin());
     }
   */

  /// Insert nodes in the range [first,last) after the child node
  /// NB the child node pointer is invalidate by the insertion and the child cursor is reset.
  /// returns the index of the child node
  inline int InsertAfter(HierarchicalDataNode* childNode,
                         std::vector<HierarchicalDataNode>::iterator first,
                         std::vector<HierarchicalDataNode>::iterator last)
  {
    std::vector<HierarchicalDataNode>::iterator itr = m_children.begin();
    int indx = 0;
    while( &(*itr) != childNode && itr != m_children.end() )
    {
      ++itr;
      ++indx;
    }
    if( itr == m_children.end() )
    {
      indx = -1;
      SLIC_ERROR("InsertAfter: Child node not found.");
    }

    ++itr;
    ++indx;
    m_children.insert(itr, first, last);
    ResetChildCursor();
    return indx;
  }

  inline HierarchicalDataNode*
  Next(const bool reset = false)
  {
    if( m_children.size() == 0 )
      return 0;
    if( reset )
      ResetChildCursor();
    if( m_itc == m_children.end() )
      return 0;
    HierarchicalDataNode* ret = &(*m_itc);
    ++m_itc;
    return ret;
  }

  inline HierarchicalDataNode*
  GetChild(const char* childName)
  {
    HierarchicalDataNode* hdn = 0;
    for( hdn = Next(true) ; hdn ; hdn = Next() )
      if( !strcmp(hdn->Heading(), childName) )
        return hdn;
    return 0;
  }

  inline uint64 CountChildren() const
  {
    return m_children.size();
  }

  inline std::vector<HierarchicalDataNode>
  CopyChildren() const
  {
    return m_children;
  }

  /// Removes a child from the parent
  /// NB the child node pointer is invalidated by the removal and the child cursor is reset.
  /// Returns the index where the child was removed (-1 if not found).
  int RemoveChild(HierarchicalDataNode* childNode)
  {
    std::vector<HierarchicalDataNode>::iterator itr = m_children.begin();
    std::vector<HierarchicalDataNode>::iterator itre;
    int indx = 0;
    while( &(*itr) != childNode && itr != m_children.end() )
    {
      ++itr;
      ++indx;
    }
    if(itr != m_children.end())
    {
      itre = itr;
      ++itre;
      m_children.erase(itr,itre);
      ResetChildCursor();
    }
    else
    {
      indx = -1; // child not found
    }
    return indx;
  }


  enum DefaultReportLevel
  {
    silent = 0,
    recordDefaults = 1,
    reportDefaults = 2,
    disableDefaults = 3
  };

  static void SetDefaultReportLevel(DefaultReportLevel rl){
    m_defaultReportLevel = rl;
  }

  // convert to XML
  std::string ToXML(bool doIndent = true){
    std::string rv;

    if(doIndent) rv += indentString(m_level);

    rv += "<" + m_heading;

    // loop over attributes
    if(m_attributes.size() > 0)
    {
      std::map<std::string, std::string>::iterator itr = m_attributes.begin();
      std::map<std::string, std::string>::iterator iend = m_attributes.end();
      for( ; itr != iend ; ++itr)
      {
        rv += "\n";
        if(doIndent) rv += indentString( std::max(m_level+1,1) );
        rv += itr->first + " = \"" + itr->second + "\"";
      }
    }


    // loop over children
    if(m_children.size() > 0)
    {
      rv += ">\n";    // end of heading

      if(m_level < 1 && doIndent) rv+="\n";     // space out lower level elements


      for(unsigned i = 0 ; i < m_children.size() ; ++i)
      {
        rv +=   m_children[i].ToXML(doIndent);

        if(m_level < 1 && doIndent) rv+="\n";        // space out lower level elements

      }

      if(doIndent) rv += indentString(m_level);
      rv += "</" + m_heading + ">";

    }
    else
    {
      rv += "/>";    // end of heading
    }

    rv += "\n";

    return rv;
  }



private:
  int m_level;
  std::string m_heading;

  //attribute data
  std::map<std::string, std::string>::iterator m_it;
  std::map<std::string, std::string> m_attributes;

  //children container data
  std::vector<HierarchicalDataNode>::iterator m_itc;
  std::vector<HierarchicalDataNode> m_children;

  static DefaultReportLevel m_defaultReportLevel;

  std::string indentString(int level){
    std::string rv;
    for(int i =0 ; i < level ; ++i)
    {
      rv += "    ";
    }
    return rv;
  }
};

/**
 *  @author walsh24
 *  @brief Overloaded function for converting strings to real64 - enables unit conversion
 *
 **/
template<>
inline real64 HierarchicalDataNode::StrVal(const std::string& str)
{
  return stringutilities::fromString<real64>(str);
}


inline real64 HierarchicalDataNode::GetAttributeOrDefault(const std::string& key,
                                                          const std::string& def)
{
  std::string val = GetAttributeString(key);
  if( val == "" )
  {
    switch(m_defaultReportLevel)
    {
    case disableDefaults:
    {
      std::string exceptionStr = std::string("Error: Attempt to set a default value when default values have been prohibited.\n") +
                                 "Parameter: " + key + "\n" +
                                 "Default Value: " + def + "\n";
      SLIC_ERROR(exceptionStr);
    }
    break;
    case reportDefaults:
      std::cout << "Setting " + key + " to default value: " + def  << std::endl;
    // fall through to next case
    case recordDefaults:
      AddAttributePair(key, def );  // add attribute to hdn - default will be included if xml written to file
    // fall through to next case
    case silent:
      /*empty*/
      break;
    default:
      break;

    }
    return StrVal<real64>(def);
  }
  else
  {
    return StrVal<real64>(val);
  }
}


template<class T>
std::vector<T> HierarchicalDataNode::GetAttributeVector(const std::string& key, std::string token) const
{
  std::vector<T> vec;
  std::string varsStr = GetAttributeString(key);
  if( !varsStr.empty() )
  {
    string_array vecStrs = stringutilities::Tokenize(varsStr, token);
    for( size_t i = 0 ; i < vecStrs.size() ; ++i )
    {
      vec.push_back(StrVal<T>(vecStrs[i]));
    }
  }
  return vec;
}

}
}
#endif
