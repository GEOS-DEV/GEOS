/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file xmlWrapper.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_FILEIO_XMLWRAPPER_HPP_
#define SRC_COMPONENTS_CORE_SRC_FILEIO_XMLWRAPPER_HPP_

#include "ArrayUtilities.hpp"
#include "codingUtilities/StringUtilities.hpp"
#include "common/DataTypes.hpp"
#include "dataRepository/DefaultValue.hpp"
#include "pugixml.hpp"
#include <sstream>

namespace geosx
{

namespace dataRepository
{
class ManagedGroup;
}


class xmlWrapper
{
public:
  using xmlDocument = pugi::xml_document;
  using xmlResult = pugi::xml_parse_result;
  using xmlNode = pugi::xml_node;
  using xmlAttribute = pugi::xml_attribute;
  using xmlTypes = pugi::xml_node_type;

  static constexpr auto filePathString = "filePath";

  xmlWrapper();
  virtual ~xmlWrapper();

  /**
   * Function to add xml nodes from included files.
   * @param targetNode the node for which to look for included children specifications
   *
   * This function looks for a "Included" node under the targetNode, loops over all subnodes under the "Included"
   * node, and then parses the file specified in those subnodes taking all the nodes in the file and adding them to
   * the targetNode.
   */
  static void addIncludedXML( xmlNode & targetNode );

  template< typename T >
  static void StringToInputVariable( T& target, string value );


  static void StringToInputVariable( R1Tensor & target, string value );

  template< typename T, int NDIM >
  static void StringToInputVariable(  LvArray::Array<T,NDIM,localIndex> & array, string const & str );


  template< typename T >
  static void ReadAttributeAsType( T & rval,
                                   string const & name,
                                   xmlNode const & targetNode,
                                   bool const required );

  template< typename T, typename T_DEF = T >
  static void ReadAttributeAsType( T & rval,
                                   string const & name,
                                   xmlNode const & targetNode ,
                                   T_DEF const & defVal );

  template< typename T >
  static typename std::enable_if_t<!dataRepository::DefaultValue<T>::has_default_value>
  ReadAttributeAsType( T & rval,
                       string const & name,
                       xmlNode const & targetNode ,
                       dataRepository::DefaultValue<T> const & defVal )
  {
    ReadAttributeAsType(rval, name, targetNode, false );
  }

  template< typename T >
  static typename std::enable_if_t<dataRepository::DefaultValue<T>::has_default_value>
  ReadAttributeAsType( T & rval,
                       string const & name,
                       xmlNode const & targetNode ,
                       dataRepository::DefaultValue<T> const & defVal )
  {
    ReadAttributeAsType(rval, name, targetNode, defVal.value );
  }
};


template< typename T >
void xmlWrapper::StringToInputVariable( T & target, string inputValue )
{
  std::istringstream ss( inputValue );
  ss>>target;
}


template< typename T >
void xmlWrapper::ReadAttributeAsType( T & rval,
                                      string const & name,
                                      xmlNode const & targetNode,
                                      bool const required  )
{
  pugi::xml_attribute xmlatt = targetNode.attribute( name.c_str() );

  GEOS_ERROR_IF( xmlatt.empty() && required , "Input variable " + name + " is required in " + targetNode.path() );

  StringToInputVariable( rval, xmlatt.value() );
}


template< typename T, typename T_DEF >
void xmlWrapper::ReadAttributeAsType( T & rval,
                                      string const & name,
                                      xmlNode const & targetNode,
                                      T_DEF const & defVal )
{
  pugi::xml_attribute xmlatt = targetNode.attribute( name.c_str() );
  if( !xmlatt.empty() )
  {
    StringToInputVariable( rval, xmlatt.value() );
  }
  else
  {
    rval = defVal;
  }
}

template< typename T, int NDIM >
void xmlWrapper::StringToInputVariable(  LvArray::Array<T,NDIM,localIndex> & array, string const & str )
{
  string collapsedStr( str );
  collapsedStr.erase(std::remove(collapsedStr.begin(), collapsedStr.end(), ' '), collapsedStr.end());

  size_t openPos = 0;
  size_t closePos = 0;

  GEOS_ERROR_IF( collapsedStr[0]!='{',
                 "First non-space character of input string for an array must be {" );

  GEOS_ERROR_IF( collapsedStr.find("{}")!=string::npos,
                 "Cannot have an empty dimension of an array, i.e. {}" );

  size_t const numOpen = std::count( collapsedStr.begin(), collapsedStr.end(), '{' );
  size_t const numClose = std::count( collapsedStr.begin(), collapsedStr.end(), '}' );

  GEOS_ERROR_IF( numOpen != numClose,
                 "Number of opening { not equal to number of } in processing of string for filling"
                 " an Array. Given string is: \n"<<collapsedStr);

  // get the number of dimensions from the number of { characters that begin the input string
  int const ndims = integer_conversion<int>(collapsedStr.find_first_not_of('{'));
  GEOS_ERROR_IF( ndims!=NDIM,
                 "number of dimensions in string ("<<ndims<<
                 ") does not match dimensions of array("<<NDIM<<
                 "). Original string is:/n"<<str );

  int dimLevel = -1;
  localIndex dims[NDIM] = {0};
  localIndex otherDims[NDIM] = {0};
  for( int i=0 ; i<NDIM ; ++i )
  {
    dims[i]=1;
    otherDims[i] = 1;
  }
  bool dimSet[NDIM] = {false};

  char lastChar = 0;

  for( size_t charCount = 0; charCount<collapsedStr.size() ; ++charCount )
  {
    char const c = collapsedStr[charCount];
    if( c=='{')
    {
      ++dimLevel;
    }
    else if( c=='}')
    {
      dimSet[dimLevel] = true;
      GEOS_ERROR_IF( dims[dimLevel]!=otherDims[dimLevel],
                     "Dimension "<<dimLevel<<" is inconsistent across the expression. "
                     "The first set value of the dimension is "<<dims[dimLevel]<<
                     " while the current value of the dimension is"<<otherDims[dimLevel]<<
                     ". The values that have been parsed prior to the error are:\n"<<
                     collapsedStr.substr(0,charCount+1) );
      otherDims[dimLevel] = 1;
      --dimLevel;
      GEOS_ERROR_IF( dimLevel<0 && charCount<(collapsedStr.size()-1),
                     "In parsing the input string, the current dimension of the array has dropped "
                     "below 0. This means that there are more '}' than '{' at some point in the"
                     " parsing. The values that have been parsed prior to the error are:\n"<<
                     collapsedStr.substr(0,charCount+1) );

    }
    else if( c==',' )
    {
      GEOS_ERROR_IF( lastChar=='{' || lastChar==',',
                     "character of ',' follows '"<<lastChar<<"'. Comma must follow an array value.");
      if( dimSet[dimLevel]==false )
      {
        ++(dims[dimLevel]);
      }
      ++(otherDims[dimLevel]);

    }

    lastChar = c;
  }
  GEOS_ERROR_IF( dimLevel!=-1,
                 "Expression fails to close all '{' with a corresponding '}'. Check your input:"<<
                 collapsedStr );


  array.resize( NDIM, dims );

  T * arrayData = array.data();

  // In order to use the stringstream facility to read in values of a Array<string>,
  // we need to replace all {}, with spaces.
  std::replace( collapsedStr.begin(), collapsedStr.end(), '{', ' ' );
  std::replace( collapsedStr.begin(), collapsedStr.end(), '}', ' ' );
  std::replace( collapsedStr.begin(), collapsedStr.end(), ',', ' ' );
  std::istringstream strstream(collapsedStr);

  // iterate through the stream and insert values into array in a linear fashion. This will be
  // incorrect if we ever have Array with a permuted index capability.
  while( strstream )
  {
    int c = strstream.peek();

    if( c== ' ' )
    {
      strstream.ignore();
    }
    else
    {
      strstream>>*(arrayData++);
    }
  }
}


} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_FILEIO_XMLWRAPPER_HPP_ */
