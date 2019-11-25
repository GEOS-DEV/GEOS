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
 * @file SchemaUtilities.hpp
 */

#ifndef GEOSX_FILEIO_SCHEMA_SCHEMAUTILITIES_HPP_
#define GEOSX_FILEIO_SCHEMA_SCHEMAUTILITIES_HPP_

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

#endif /* GEOSX_FILEIO_SCHEMA_SCHEMAUTILITIES_HPP_ */
