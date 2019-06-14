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
  string str_nospace( str );
  str_nospace.erase(std::remove(str_nospace.begin(), str_nospace.end(), ' '), str_nospace.end());


  string const openDelim("{");
  string const closeDelim("}");
  size_t openPos = 0;
  size_t closePos = 0;


  size_t const numOpen = std::count( str_nospace.begin(), str_nospace.end(), '{' );
  size_t const numClose = std::count( str_nospace.begin(), str_nospace.end(), '}' );

  GEOS_ERROR_IF( numOpen != numClose,
                 "Number of opening { not equal to number of } in processing of string for filling"
                 " an Array. Given string is: \n"<<str_nospace);


  int ndims = 0;
  int dimLevel = -1;
  localIndex dims[NDIM] = {0};
  for( int i=0 ; i<NDIM ; ++i ) dims[i]=1;
  bool dimSet[NDIM] = {false};
  bool reachedClose = false;

  for( char const & c : str_nospace )
  {
    if( c=='{')
    {
      if( reachedClose==false )
      {
        ++ndims;
      }
      ++dimLevel;
    }
    else if( c=='}')
    {
      dimSet[dimLevel] = true;
      --dimLevel;
      reachedClose = true;
    }
    else if( c==',' )
    {
      if( dimSet[dimLevel]==false )
      {
        ++dims[dimLevel];
      }
    }
  }

  GEOS_ERROR_IF( ndims!=NDIM,
                 "number of dimensions in string ("<<ndims<<
                 ") does not match dimensions of array("<<NDIM<<
                 "). Original string is:/n"<<str );
  array.resize( NDIM, dims );

  T * arrayData = array.data();

  // we need this because array<string> will capture the } in the value.
  std::replace( str_nospace.begin(), str_nospace.end(), '}', ' ');
  std::replace( str_nospace.begin(), str_nospace.end(), ',', ' ');
  std::istringstream strstream(str_nospace);

  ndims = 0;
  dimLevel = -1;
  while( strstream )
  {
    int c = strstream.peek();

    if( c=='{' || c== ' ' )
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
