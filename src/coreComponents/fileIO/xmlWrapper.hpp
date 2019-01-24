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

/*
 * xmlWrapper.hpp
 *
 *  Created on: Jun 3, 2017
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_FILEIO_XMLWRAPPER_HPP_
#define SRC_COMPONENTS_CORE_SRC_FILEIO_XMLWRAPPER_HPP_

#include "ArrayUtilities.hpp"
#include "common/DataTypes.hpp"
#include "dataRepository/DefaultValue.hpp"
#include "pugixml.hpp"
#include <sstream>

namespace cxx_utilities
{
class DocumentationNode;
}
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

  xmlWrapper();
  virtual ~xmlWrapper();

  template< typename T >
  static void StringToInputVariable( T& target, string value );

  template< typename T >
  static void StringToInputVariable( array1d<T> & target, string value );

  template< typename T >
  static void StringToInputVariable( array2d<T> & target, string value );

  static void StringToInputVariable( R1Tensor & target, string value );

//  static void StringToInputVariable( R2Tensor & target, string value );
//
//  static void StringToInputVariable( R2SymTensor & target, string value );

//  template< typename T >
//  static T StringToInputVariable( xmlNode const & node, string const name, T defValue
// );

//  static R1Tensor StringToInputVariable( xmlNode const & node, string const name, R1Tensor defValue );

//  static void ReadAttributeAsType( dataRepository::ManagedGroup & group,
//                                   cxx_utilities::DocumentationNode const & subDocNode,
//                                   xmlNode const & targetNode );

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
void xmlWrapper::StringToInputVariable( array1d<T> & target, string inputValue )
{
  string csvstr = inputValue;
  std::istringstream ss( csvstr );

  T value;

  while( ss.peek() == ',' || ss.peek() == ' ' )
  {
    ss.ignore();
  }
  while( !((ss>>value).fail()) )
  {
    target.push_back( value );
    while( ss.peek() == ',' || ss.peek() == ' ' )
    {
      ss.ignore();
    }
  }
}

template< typename T >
void xmlWrapper::StringToInputVariable( array2d<T> & target, string inputValue )
{
  array1d<T> temp;
  StringToInputVariable( temp, inputValue );

  target.resize(1,temp.size());
  for( localIndex i=0 ; i<temp.size() ; ++i )
  {
    target[0][i] = temp[i];
  }
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




} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_FILEIO_XMLWRAPPER_HPP_ */
