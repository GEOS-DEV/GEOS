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
 * @file SchemaUtilities.hpp
 */

#ifndef SCHEMAUTILITIES_HPP_
#define SCHEMAUTILITIES_HPP_

#include "dataRepository/xmlWrapper.hpp"
#include "common/DataTypes.hpp"
#include "dataRepository/DefaultValue.hpp"
#include <iostream>
#include <sstream>
#include <string>

namespace geosx
{

class SchemaUtilities
{
public:

SchemaUtilities();
virtual ~SchemaUtilities();

static void ConvertDocumentationToSchema(std::string const & fname, dataRepository::Group * const group, integer documentationType);
static void BuildSimpleSchemaTypes(xmlWrapper::xmlNode schemaRoot);
static void SchemaConstruction(dataRepository::Group * const group, xmlWrapper::xmlNode schemaRoot, xmlWrapper::xmlNode schemaParent, integer documentationType);



/// These functions are used to convert default values into strings for the schema
template< typename T >
static string DefaultValueToString( T& target )
{
  std::stringstream ss;
  ss << target;
  return ss.str();
}


template< typename T >
static string DefaultValueToString( array1d<T> & target )
{
  string csvstr = DefaultValueToString(target[0]);
  for (integer ii=1; ii<target.size(); ++ii)
  {
    csvstr += ", " + DefaultValueToString(target[ii]);
  }

  return csvstr;
}


template< typename T >
static string DefaultValueToString( array2d<T> & target )
{
  string csvstr = "{{" + DefaultValueToString< array1d<T> >( target[0] ) + "}";

  for( integer ii=1 ; ii<target.size(); ++ii )
  {
    csvstr += ", {" + DefaultValueToString< array1d<T> >( target[ii] ) + "}";
  }
  csvstr += "}";

  return csvstr;
}



template< typename T >
static typename std::enable_if_t<dataRepository::DefaultValue<T>::has_default_value>
SetDefaultValueString( dataRepository::DefaultValue<T> const & defVal, xmlWrapper::xmlNode attributeNode )
{
  attributeNode.append_attribute("default") = DefaultValueToString(defVal.value).c_str();
}


};


}

#endif /* SCHEMAUTILITIES_HPP_ */
