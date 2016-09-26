/*
 * SchemaUtilities.hpp
 *
 *  Created on: Sep 15, 2016
 *      Author: sherman
 */

#ifndef SCHEMAUTILITIES_HPP_
#define SCHEMAUTILITIES_HPP_

#include "DocumentationNode.hpp"
#include "pugixml/src/pugixml.hpp"

namespace geosx
{

void ConvertDocumentationToSchema(std::string const & fname, cxx_utilities::DocumentationNode const & inputDocumentationHead);
void SchemaConstruction(cxx_utilities::DocumentationNode const & docNode, pugi::xml_node schemaNode, pugi::xml_node schemaRoot);

}

#endif /* SCHEMAUTILITIES_HPP_ */
