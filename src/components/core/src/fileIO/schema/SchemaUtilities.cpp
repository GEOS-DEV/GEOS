/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
 * @file SchemaUtilities.cpp
 * @author sherman
 */

#include "SchemaUtilities.hpp"

namespace geosx
{

void ConvertDocumentationToSchema(std::string const & fname, cxx_utilities::DocumentationNode const & inputDocumentationHead, integer verbosityLevel)
{
  std::string schemaBase=
    "<?xml version=\"1.1\" encoding=\"ISO-8859-1\" ?>\
  <xsd:schema xmlns:xsd=\"http://www.w3.org/2001/XMLSchema\">\
  <xsd:annotation>\
  <xsd:documentation xml:lang=\"en\">New schema for GEOS</xsd:documentation>\
  </xsd:annotation>\
  </xsd:schema>";

  xmlWrapper::xmlDocument schemaTree;
  schemaTree.load_string(schemaBase.c_str());
  xmlWrapper::xmlNode schemaRoot = schemaTree.child("xsd:schema");

  // Build the simple schema types
  BuildSimpleSchemaTypes(schemaRoot);

  // Recursively build the schema from the documentation string
  SchemaConstruction(inputDocumentationHead, schemaRoot, schemaRoot, verbosityLevel);

  // Write the schema to file
  schemaTree.save_file(fname.c_str());
}


void BuildSimpleSchemaTypes(xmlWrapper::xmlNode schemaRoot)
{
  rtTypes::typeRegex typeRegex;

  for( auto regex=typeRegex.begin() ; regex!=typeRegex.end() ; ++regex )
  {
    xmlWrapper::xmlNode newNode = schemaRoot.append_child("xsd:simpleType");
    newNode.append_attribute("name") = regex->first.c_str();
    xmlWrapper::xmlNode restrictionNode = newNode.append_child("xsd:restriction");
    restrictionNode.append_attribute("base") = "xsd:string";
    xmlWrapper::xmlNode patternNode = restrictionNode.append_child("xsd:pattern");

    // Default regex to string
    if( regex->second.empty())
    {
      std::cout << "Warning: schema regex not defined for " << regex->first << "...  Defaulting to limited string" << std::endl;
      patternNode.append_attribute("value") = "[a-zA-Z0-9_,\\(\\)+-/\\*]*";
    }
    else
    {
      patternNode.append_attribute("value") = regex->second.c_str();
    }
  }
}


void SchemaConstruction(cxx_utilities::DocumentationNode const & docNode, xmlWrapper::xmlNode schemaNode, xmlWrapper::xmlNode schemaRoot,
                        integer verbosityLevel)
{
  if( docNode.getVerbosity() > verbosityLevel )
  {
    // Do nothing
  }
  else if( docNode.getSchemaType().find("Node") != std::string::npos )
  {
    // Special case for root nodes:
    xmlWrapper::xmlNode targetNode = schemaRoot;
    if( docNode.getSchemaType().find("Root") == std::string::npos )
    {
      targetNode = schemaNode.child("xsd:choice");
      if( targetNode.empty() )
      {
        targetNode = schemaNode.prepend_child("xsd:choice");
        targetNode.append_attribute("maxOccurs") = "unbounded";
      }
    }

    // Insert the schema node if not present, then iterate over children
    if( targetNode.find_child_by_attribute("xsd:element", "name", docNode.m_name.c_str()).empty())
    {
      // Add the entries to the current and root nodes
      xmlWrapper::xmlNode newNode = targetNode.append_child("xsd:element");
      newNode.append_attribute("name") = docNode.m_name.c_str();
      newNode.append_attribute("type") = (docNode.m_name+"Type").c_str();

      // Set the occurance limits
      if( docNode.getSchemaType().find("Required") != std::string::npos )
      {
        newNode.append_attribute("minOccurs") = "1";
      }
      if( docNode.getSchemaType().find("Unique") != std::string::npos )
      {
        newNode.append_attribute("maxOccurs") = "1";
      }

      // Insert a new type into the root node if not present
      newNode = schemaRoot.find_child_by_attribute("xsd:complexType", "name", (docNode.m_name+"Type").c_str());
      if( newNode.empty())
      {
        newNode = schemaRoot.append_child("xsd:complexType");
        newNode.append_attribute("name") = (docNode.m_name+"Type").c_str();
      }

      // Set the type characteristics
      for( auto const & subDocNode : docNode.m_child )
      {
        SchemaConstruction(subDocNode.second, newNode, schemaRoot, verbosityLevel);
      }
    }
  }
  else if( docNode.getSchemaType().empty() == 0 )
  {
    // Insert a new attribute if not present
    if( schemaNode.find_child_by_attribute("xsd:attribute", "name", docNode.m_name.c_str()).empty())
    {
      xmlWrapper::xmlNode newNode = schemaNode.append_child("xsd:attribute");
      newNode.append_attribute("name") = docNode.m_name.c_str();
      newNode.append_attribute("type") = (docNode.getSchemaType()).c_str();

      if( !strcmp(docNode.getDefault().c_str(), "REQUIRED"))
      {
        newNode.append_attribute("use") = "required";
      }
    }
  }
}

}
