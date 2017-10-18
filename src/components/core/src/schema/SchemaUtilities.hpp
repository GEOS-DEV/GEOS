/*
 * SchemaUtilities.hpp
 *
 *  Created on: Sep 15, 2016
 *      Author: sherman
 */

#ifndef SCHEMAUTILITIES_HPP_
#define SCHEMAUTILITIES_HPP_

#include "DocumentationNode.hpp"
#include "common/DataTypes.hpp"
#include "fileIO/xmlWrapper.hpp"
namespace geosx
{


void ConvertDocumentationToSchema(std::string const & fname, cxx_utilities::DocumentationNode const & inputDocumentationHead, integer verbosityLevel);
void BuildSimpleSchemaTypes(xmlWrapper::xmlNode schemaRoot);
void SchemaConstruction(cxx_utilities::DocumentationNode const & docNode, xmlWrapper::xmlNode schemaNode, xmlWrapper::xmlNode schemaRoot, integer verbosityLevel);

}

#endif /* SCHEMAUTILITIES_HPP_ */
