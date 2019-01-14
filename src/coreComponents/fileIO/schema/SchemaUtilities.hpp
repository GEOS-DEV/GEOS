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
 * SchemaUtilities.hpp
 *
 *  Created on: Sep 15, 2016
 *      Author: sherman
 */

#ifndef SCHEMAUTILITIES_HPP_
#define SCHEMAUTILITIES_HPP_

#include "fileIO/xmlWrapper.hpp"
#include "common/DataTypes.hpp"
#include <iostream>
#include <string>

namespace geosx
{

void ConvertDocumentationToSchema(std::string const & fname, dataRepository::ManagedGroup * const group);
void BuildSimpleSchemaTypes(xmlWrapper::xmlNode schemaRoot);
void SchemaConstruction(dataRepository::ManagedGroup * const group, xmlWrapper::xmlNode schemaRoot, xmlWrapper::xmlNode schemaParent);


/*
/// These functions are used to convert default values into strings for the schema
template< typename T >
static string DefaultValueToString( T& target );

template< typename T >
static string DefaultValueToString( array1d<T> & target );

template< typename T >
static string DefaultValueToString( array2d<T> & target );



template< typename T >
string DefaultValueToString( T& target )
{
  return std::to_string(target);
}


template< typename T >
string DefaultValueToString( array1d<T> & target )
{
  string csvstr = std::to_string(target[0]);
  for (integer ii=1; ii<target.size(); ++ii)
  {
    csvstr += ", " + std::to_string(target[ii]);
  }

  return csvstr;
}


template< typename T >
string DefaultValueToString( array2d<T> & target )
{
  string csvstr = "{{" + DefaultValueToString< array1d<T> >( target[0] ) + "}";

  for( integer ii=1 ; ii<target.size(); ++ii )
  {
    csvstr += ", {" + DefaultValueToString< array1d<T> >( target[ii] ) + "}";
  }
  csvstr += "}";

  return csvstr;
}
*/


}

#endif /* SCHEMAUTILITIES_HPP_ */
