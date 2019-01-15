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
 * @file SchemaUtilities.cpp
 * @author sherman
 */

#include "SchemaUtilities.hpp"
#include "common/DataTypes.hpp"
#include "dataRepository/ViewWrapper.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "dataRepository/SchemaFlags.hpp"
#include "dataRepository/InputFlags.hpp"

namespace geosx
{

using namespace dataRepository;


SchemaUtilities::SchemaUtilities()
{
  // TODO Auto-generated constructor stub

}

SchemaUtilities::~SchemaUtilities()
{
  // TODO Auto-generated destructor stub
}


void SchemaUtilities::ConvertDocumentationToSchema(std::string const & fname, ManagedGroup * const group)
{
  std::cout << "\nGenerating XML Schema..." << std::endl;

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
  std::cout << "  Basic datatypes" << std::endl;
  BuildSimpleSchemaTypes(schemaRoot);

  // Recursively build the schema from the data structure skeleton
  std::cout << "  Data structure layout" << std::endl;
  SchemaConstruction(group, schemaRoot, schemaRoot);

  // Write the schema to file
  std::cout << "  Saving file"<< std::endl;
  schemaTree.save_file(fname.c_str());

  std::cout << "  Done!" << std::endl;
}


void SchemaUtilities::BuildSimpleSchemaTypes(xmlWrapper::xmlNode schemaRoot)
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
      GEOS_WARNING("schema regex not defined for " << regex->first << "...  Defaulting to limited string");
      patternNode.append_attribute("value") = "[a-zA-Z0-9_,\\(\\)+-/\\*]*";
    }
    else
    {
      patternNode.append_attribute("value") = regex->second.c_str();
    }
  }
}


void SchemaUtilities::SchemaConstruction(ManagedGroup * const group, xmlWrapper::xmlNode schemaRoot, xmlWrapper::xmlNode schemaParent)
{
  // Get schema details
  SchemaFlags schemaType = group->getSchemaFlags();
  
  if (schemaType != SchemaFlags::IGNORE)
  {
    string targetName = group->getName();
    string typeName = targetName + "Type";

    // Insert the schema node if not present, then iterate over children
    if( schemaParent.find_child_by_attribute("xsd:element", "name", targetName.c_str()).empty())
    {
      // Add the entries to the current and root nodes
      xmlWrapper::xmlNode targetIncludeNode = schemaParent.append_child("xsd:element");
      targetIncludeNode.append_attribute("name") = targetName.c_str();
      targetIncludeNode.append_attribute("type") = typeName.c_str();

      // Add occurence conditions
      if((schemaType == SchemaFlags::REQUIRED_NODE) || (schemaType == SchemaFlags::REQUIRED_UNIQUE_NODE))
      {
        targetIncludeNode.append_attribute("minOccurs") = "1";
      }
      if((schemaType == SchemaFlags::UNIQUE_NODE) || (schemaType == SchemaFlags::REQUIRED_UNIQUE_NODE))
      {
        targetIncludeNode.append_attribute("maxOccurs") = "1";
      }

      // Insert a new type into the root node if not present
      xmlWrapper::xmlNode targetTypeDefNode = schemaRoot.find_child_by_attribute("xsd:complexType", "name", typeName.c_str());
      if( targetTypeDefNode.empty())
      {
        targetTypeDefNode = schemaRoot.append_child("xsd:complexType");
        targetTypeDefNode.append_attribute("name") = typeName.c_str();

        // Add subgroups
        if (group->numSubGroups() > 0)
        {
          xmlWrapper::xmlNode targetChoiceNode = targetTypeDefNode.child("xsd:choice");
          if( targetChoiceNode.empty() )
          {
            targetChoiceNode = targetTypeDefNode.prepend_child("xsd:choice");

            // Add children of the group
            group->forSubGroups<ManagedGroup>([&]( ManagedGroup * subGroup ) -> void
            {
              SchemaConstruction(subGroup, schemaRoot, targetChoiceNode);
            });

            // Add schema deviations
            group->SetSchemaDeviations(schemaRoot, targetChoiceNode);
          }
        }

        // Add attributes
        for ( auto wrapperPair : group->wrappers() )
        {
          ViewWrapperBase * const wrapper = wrapperPair.second;
          InputFlags flag = wrapper->getInputFlag();
          
          if (( flag == InputFlags::OPTIONAL ) || ( flag == InputFlags::REQUIRED ))
          {
            // Basic attributes
            xmlWrapper::xmlNode attributeNode = targetTypeDefNode.append_child("xsd:attribute");
            attributeNode.append_attribute("name") = wrapper->getName().c_str();
            attributeNode.append_attribute("type") = (rtTypes::typeNames(wrapper->get_typeid()).c_str());

            // (Optional) Description
            string description = wrapper->getDescription();
            if (!description.empty())
            {
              attributeNode.append_attribute("description") = description.c_str();
            }

            // (Optional) Default Value
            if ( flag == InputFlags::OPTIONAL )
            {
              rtTypes::TypeIDs const wrapperTypeID = rtTypes::typeID(wrapper->get_typeid());
              rtTypes::ApplyIntrinsicTypeLambda2( wrapperTypeID,
                                                  [&]( auto a, auto b ) -> void
              {
                using COMPOSITE_TYPE = decltype(a);
                ViewWrapper<COMPOSITE_TYPE>& typedWrapper = ViewWrapper<COMPOSITE_TYPE>::cast( *wrapper );
                dataRepository::DefaultValue<COMPOSITE_TYPE> tmp = typedWrapper.getDefaultValueStruct();

                if( typedWrapper.getDefaultValueStruct().has_default_value )
                {
                  // string stringVal = std::to_string(typedWrapper.getDefaultValue());
                  attributeNode.append_attribute("default") = "default";

                  GetDefaultValueString( typedWrapper.getDefaultValueStruct() );
                }
              });
            }
          }
        }
      }
    }
  }
}

}
