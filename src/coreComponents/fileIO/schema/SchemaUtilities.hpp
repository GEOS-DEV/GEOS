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

};


}

#endif /* GEOSX_FILEIO_SCHEMA_SCHEMAUTILITIES_HPP_ */
