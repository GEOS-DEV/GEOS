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

namespace geosx
{

void ConvertDocumentationToSchema(std::string const & fname, dataRepository::ManagedGroup * const group);
void BuildSimpleSchemaTypes(xmlWrapper::xmlNode schemaRoot);
void SchemaConstruction(dataRepository::ManagedGroup * const group, xmlWrapper::xmlNode schemaRoot, xmlWrapper::xmlNode schemaParent);

}

#endif /* SCHEMAUTILITIES_HPP_ */
